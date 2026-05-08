"""
PDH プロセスの 1 パス計算 (リサイクル収束の 1 反復に相当)。

入力:
  - tear streams (Mem retentate / Dist3 bottom それぞれの A=C3H8, B=C3H6)
  - Fresh LPG 流量 (C3H8, C4H10)
  - 設計変数バンドル (FlowsheetDesignVars)
  - 運転条件 (OperatingConfig)

戻り値: dict (各ユニットの結果と tear stream 更新値・温度)

外側のソルバ (flowsheet/solver.py) がこの関数を反復呼び出しして
リサイクルを収束させる。
"""

import warnings

from stream.stream import ProcessStream
from units.utils.mixer import mix_streams
from units.utils.cooler import simulate_cooler
from units.utils.compressor import simulate_compressor
from units.separators.column1.fake_column1 import simulate_column1
from units.separators.column2.fake_column2 import simulate_column2
from units.separators.column3.fake_column3 import simulate_column3
from units.reactors.swing import (
    FeedStream as SwingFeed, FixedParams as SwingFixed,
    simulate_swing_reactor_system,
)
from units.separators.psa.psa_system import (
    PSAFeedStream, PSAFixedParams, simulate_psa_system,
)
from units.separators.membrane.membrane_system import (
    MemFeedStream, MemFixedParams, simulate_membrane_system,
)

from flowsheet.design import FlowsheetDesignVars
from config.load import OperatingConfig

_ZERO = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'D': 0.0, 'E': 0.0, 'F': 0.0, 'Z': 0.0}


def run_one_pass(
    tear_dist3:   dict,
    tear_mem:     dict,
    T_d3:         float,
    T_mem:        float,
    F_C3H8_feed:  float,
    F_C4H10_feed: float,
    design:       FlowsheetDesignVars,
    config:       OperatingConfig,
) -> dict:
    """1 パスのプロセスシミュレーションを実行。

    Parameters
    ----------
    tear_dist3 : dict {'A','B'} [kmol/h]
        Dist3 塔底由来のリサイクル (前反復値)
    tear_mem : dict {'A','B'} [kmol/h]
        Mem 保留側由来のリサイクル (前反復値)
    T_d3, T_mem : float [K]
        各リサイクルの温度 (前反復値)
    F_C3H8_feed, F_C4H10_feed : float [kmol/h]
        Fresh LPG 中 C3H8, C4H10 流量 (外側で調整)
    design : FlowsheetDesignVars
        設計変数 (SWING, PSA, MEM)
    config : OperatingConfig
        運転条件 (圧力, 温度, 製品仕様, 原料状態)

    Returns
    -------
    dict
        全ユニット結果 + tear_*_new + T_*_new
    """
    P_atm = config.pressure.reactor_inlet_Pa  # 反応器入口・膨張弁後の圧力

    # ---- Fresh LPG (常時固定: 25°C, 1 atm) ----
    fresh = ProcessStream(
        F_in={'A': F_C3H8_feed, 'Z': F_C4H10_feed,
              'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=config.feed.T_K, P_in=config.feed.P_Pa,
    )

    # ---- Step 1: Comp1 → Dist1 (Fresh のみ) ----
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp1 = simulate_compressor(fresh, P_out_target=config.pressure.comp1_out_Pa)
        r1 = simulate_column1(comp1.outlet)

    # 塔頂を 1 atm に膨張 (C4 除去済みの C3 主成分)
    dist1_top_1atm = ProcessStream(
        F_in=dict(r1.top.F_in), T_in=r1.top.T_in, P_in=P_atm,
    )

    # ---- リサイクルストリーム (1 atm 膨張弁経由、コストなし) ----
    recycle_dist3 = ProcessStream(
        F_in={**_ZERO, 'A': tear_dist3['A'], 'B': tear_dist3['B']},
        T_in=T_d3, P_in=P_atm,
    )
    recycle_mem = ProcessStream(
        F_in={**_ZERO, 'A': tear_mem['A'], 'B': tear_mem['B']},
        T_in=T_mem, P_in=P_atm,
    )

    # ---- Reactor 入口で合流 ----
    reactor_inlet = mix_streams([dist1_top_1atm, recycle_dist3, recycle_mem])

    # ---- Step 2: Swing Reactor ----
    swing_feed = SwingFeed(
        F_in=reactor_inlet.F_in,
        T_feed=reactor_inlet.T_in,
        P_in=reactor_inlet.P_in,
    )
    r_rx = simulate_swing_reactor_system(design.swing, swing_feed, SwingFixed())

    rx_out = ProcessStream(
        F_in=r_rx.effluent.F_out_avg,
        T_in=r_rx.effluent.T_out_avg,
        P_in=r_rx.effluent.P_out,
    )

    # ---- Step 3: Cooler → Comp2 → Dist2 ----
    cooled = simulate_cooler(rx_out, T_out_target=config.temperature.cooler_after_reactor_K)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp2 = simulate_compressor(cooled.outlet, P_out_target=config.pressure.comp2_out_Pa)
        r2 = simulate_column2(comp2.outlet)

    # ---- Step 4: PSA ----
    psa_feed = PSAFeedStream(
        F_in=r2.top.F_in, T_in=r2.top.T_in, P_in=r2.top.P_in,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_psa = simulate_psa_system(design.psa, psa_feed, PSAFixedParams())

    # ---- Step 5: Membrane ----
    # Dist2 を 8.5 bar 運転にした影響で塔底液の泡点が ~20°C まで下がり、冷却水
    # では液状態を保てない。膜の P_H <= 9.5 bar (Hua et al. 2024) を守るため、
    # 塔底液を mem_feed_K まで気化・過熱してガスフィードで膜へ送る。
    mem_precool = simulate_cooler(r2.bottom, T_out_target=config.temperature.mem_feed_K)
    mem_feed = MemFeedStream(
        F_C3H6=mem_precool.outlet.F_in.get('B', 0.),
        F_C3H8=mem_precool.outlet.F_in.get('A', 0.),
        T_in=mem_precool.outlet.T_in,
        P_in=mem_precool.outlet.P_in,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_mem = simulate_membrane_system(design.mem, mem_feed, MemFixedParams(vapor_feed=True))

    # ---- Step 6: Dist3 ----
    mem_to_dist3 = ProcessStream(
        F_in={'A': r_mem.product.F_C3H8, 'B': r_mem.product.F_C3H6,
              'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=r_mem.product.T_out, P_in=r_mem.product.P_out,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r3 = simulate_column3(mem_to_dist3)

    # ---- tear stream の更新値 ----
    tear_dist3_new = {
        'A': r3.bottom.F_in.get('A', 0.0),
        'B': r3.bottom.F_in.get('B', 0.0),
    }
    tear_mem_new = {
        'A': r_mem.retentate.F_C3H8,
        'B': r_mem.retentate.F_C3H6,
    }
    T_d3_new  = r3.bottom.T_in
    T_mem_new = r_mem.retentate.T_out

    return dict(
        comp1=comp1, r1=r1, dist1_top_1atm=dist1_top_1atm,
        reactor_inlet=reactor_inlet,
        r_rx=r_rx, rx_out=rx_out,
        cooled=cooled, comp2=comp2, r2=r2,
        r_psa=r_psa, mem_precool=mem_precool, r_mem=r_mem, r3=r3,
        tear_dist3_new=tear_dist3_new, tear_mem_new=tear_mem_new,
        T_d3_new=T_d3_new, T_mem_new=T_mem_new,
    )
