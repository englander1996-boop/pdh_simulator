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

import math
import warnings

from stream.stream import ProcessStream
from units.utils.mixer import mix_streams
from units.utils.cooler import simulate_cooler
from flowsheet.heat_integration import StreamPhase
from units.utils.compressor import simulate_compressor
from units.utils.pump import simulate_pump
from units.utils.expansion_valve import simulate_jt_expansion
from units.separators.column1.column1 import simulate_column1
from units.separators.column2.column2 import simulate_column2
from units.separators.column3.column3 import simulate_column3
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
    P_rx = config.pressure.reactor_inlet_Pa  # 反応器入口・膨張弁後の圧力 (contest 規定 0.5 bar)

    # ---- Fresh LPG 原料 (30°C 飽和液、C3H8:C4H10 = 9:1, ~9.97 bar) ----
    fresh = ProcessStream(
        F_in={'A': F_C3H8_feed, 'Z': F_C4H10_feed,
              'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=config.feed.T_K, P_in=config.feed.P_Pa,
    )

    # ---- Step 1: Pump1 → Dist1 (Fresh のみ) ----
    # 原料は液 (30°C 飽和液) のため、Dist1 (17 bar) への昇圧は液送ポンプで行う。
    # contest §3-3-3 「加圧すべき箇所には、ポンプ(液) を入れること」に従う。
    # 設計判断 (2026-05-09): 出口圧力は design.dist1.P_col に同期する
    # (operating.toml の pump1_out_Pa は backward compat のため残置だが未使用)。
    pump1 = simulate_pump(fresh, P_out_target=design.dist1.P_col)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r1 = simulate_column1(pump1.outlet, tunables=design.dist1)

    # 塔頂を反応器圧力 (0.5 bar) に膨張 (C4 除去済みの C3 主成分)
    # 設計判断 (2026-05-08): 旧版は P を書き換えるだけで T を維持していたが
    # 物理的には等エンタルピー (JT) 膨張で温度低下する。反応器入口プレヒート
    # Q_preheat の見積もり精度を上げるため simulate_jt_expansion を経由する。
    # 膨張弁本体はコストフリー (装置 CAPEX/OPEX に計上しない、配管中の絞り弁)。
    dist1_top_rx = simulate_jt_expansion(
        ProcessStream(F_in=dict(r1.top.F_in), T_in=r1.top.T_in, P_in=r1.top.P_in),
        P_out=P_rx,
    )

    # ---- リサイクルストリーム (膨張弁経由、コストなし) ----
    # tear stream の元圧力:
    #   recycle_dist3 : Dist3 塔底 = design.dist3.P_col
    #   recycle_mem   : Mem 保留   = design.mem.P_H
    recycle_dist3 = simulate_jt_expansion(
        ProcessStream(
            F_in={**_ZERO, 'A': tear_dist3['A'], 'B': tear_dist3['B']},
            T_in=T_d3, P_in=design.dist3.P_col,
        ),
        P_out=P_rx,
    )
    recycle_mem = simulate_jt_expansion(
        ProcessStream(
            F_in={**_ZERO, 'A': tear_mem['A'], 'B': tear_mem['B']},
            T_in=T_mem, P_in=design.mem.P_H,
        ),
        P_out=P_rx,
    )

    # ---- Reactor 入口で合流 ----
    reactor_inlet = mix_streams([dist1_top_rx, recycle_dist3, recycle_mem])

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

    # ---- Step 3: Cooler → Comp2 (2 段+段間冷却) → Dist2 ----
    # 反応器圧力 0.5 bar → Dist2 圧力 8.5 bar = 圧縮比 17:1。
    # 単段では断熱温度上昇が過大になる (~T_in×4) ため等圧縮比 √17≈4.12 の 2 段
    # に分割し、段間で T_intercool まで冷却して動力と機械的負荷を低減する。
    # 設計判断 (2026-05-09): Comp2 最終出口圧力は design.dist2.P_col に同期する
    # (operating.toml の comp2_out_Pa は backward compat のため残置だが未使用)。
    cooled = simulate_cooler(
        rx_out,
        T_out_target=config.temperature.cooler_after_reactor_K,
        process_phase=StreamPhase.GAS,
    )
    P_in_comp2  = cooled.outlet.P_in
    P_out_final = design.dist2.P_col
    P_mid       = math.sqrt(P_in_comp2 * P_out_final)
    T_intercool = config.temperature.cooler_after_reactor_K
    # 設計判断 (2026-05-09): Comp2b 出口 (~151°C) は Dist2 dew point (~50°C @ 8.5 bar)
    # を遥かに上回る超加熱蒸気。Dist2 partial condenser に直接入れると顕熱を冷凍冷媒で
    # 処理することになるため、工業実機の desuperheater (冷却水 HE) で dew 直上まで冷却。
    # Q ~12 MW を冷却水 (60 円/GJ) で除去 → 冷凍冷媒 (~1820 円/GJ) より遥かに安い。
    T_dist2_feed_K = 323.15   # 50°C: 8.5bar dew (~40-50°C) より少し上、5K margin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp2a    = simulate_compressor(cooled.outlet, P_out_target=P_mid)
        intercool = simulate_cooler(comp2a.outlet, T_out_target=T_intercool,
                                    process_phase=StreamPhase.GAS)
        comp2b    = simulate_compressor(intercool.outlet, P_out_target=P_out_final)
        desuper   = simulate_cooler(comp2b.outlet, T_out_target=T_dist2_feed_K,
                                    process_phase=StreamPhase.GAS)
        r2        = simulate_column2(desuper.outlet, tunables=design.dist2)

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
    # 設計判断 (2026-05-08): 旧版は感熱のみで潜熱無視 → Mem 気化器 OPEX が 0
    # になっていた既知バグ。phase_change=True で潜熱を加算する。
    mem_precool = simulate_cooler(
        r2.bottom,
        T_out_target=config.temperature.mem_feed_K,
        phase_change=True,
        process_phase=StreamPhase.LIQUID,    # 顕熱区間: 液相加熱、潜熱区間は EVAPORATING に自動切替
    )
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
        r3 = simulate_column3(mem_to_dist3, tunables=design.dist3)

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
        pump1=pump1, r1=r1, dist1_top_rx=dist1_top_rx,
        reactor_inlet=reactor_inlet,
        r_rx=r_rx, rx_out=rx_out,
        cooled=cooled,
        comp2a=comp2a, intercool=intercool, comp2b=comp2b,
        desuper=desuper,
        r2=r2,
        r_psa=r_psa, mem_precool=mem_precool, r_mem=r_mem, r3=r3,
        tear_dist3_new=tear_dist3_new, tear_mem_new=tear_mem_new,
        T_d3_new=T_d3_new, T_mem_new=T_mem_new,
    )
