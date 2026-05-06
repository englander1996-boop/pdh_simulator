"""
exp2.py — リサイクルあり PDH プロセス全体フロー シミュレーション

リサイクル構成:
  - Membrane 保留側 (C3H8 富化, 残留 C3H6 含)        ─┐
  - Dist3 塔底       (未透過 C3H8, 残留 C3H6 含)     ─┴→ Reactor 直前で合流
  Dist1 (脱ブタン) には戻さない
    ↳ 軽質ガス・C4 は既に系外へ抜けているためエネルギーロスを避ける

全体フロー:
  Fresh LPG (1atm) → Comp1 → Dist1 (脱ブタン)
       ↓ 塔頂(C3) を 1atm 膨張
   ┌→ Mixer ←┬─ Membrane 保留側 (1atm 膨張)
   │         └─ Dist3 塔底     (1atm 膨張)
   │
   └→ Reactor (PDH 600°C, 1atm)
       → Cooler → Comp2 → Dist2 (脱エタン)
            塔頂 → PSA (H2 製品 + オフガス)
            塔底 → Mem前冷却 → Membrane
                     透過側  → Dist3 → 塔頂 C3H6 製品
                                       塔底 → Recycle ↑
                     保留側 → Recycle ↑

収束方式:
  逐次置換 + アンダーリラックス (alpha=0.7)
  tear stream: Mem retentate (A,B) と Dist3 bottom (A,B) の 2 本（合計 4 変数）
  収束基準  : ‖tear_new − tear‖_∞ < TOL_F
  ペナルティガード: PSA/Mem の CAPEX_total ≥ 1e8 で打ち切り
  暴走ガード      : tear_mem.A が Fresh × RECYCLE_GUARD_RATIO 超で打ち切り

差し替えポイント:
  1. 蒸留塔モデル (fake_column1/2/3 → 正式 VLE モデル):
       同一インターフェイス (ProcessStream → DistResult) なら import 行のみ変更
  2. ユーティリティ単価:
       src/cost_parameters.py の以下 7 定数を実値に置換するだけ
         ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
         COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
         CATALYST_PTSN_JPY_PER_KG, CATALYST_PTSN_LIFE_YEARS,
         OPERATING_HOURS_PER_YEAR
"""

import os
import sys
import warnings

# Windows コンソール (cp932) でも記号を表示
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from units.utils.process_stream import ProcessStream
from units.utils.mixer import mix_streams
from units.utils.cooler import simulate_cooler
from units.utils.compressor import simulate_compressor
from units.separators.column1.fake_column1 import simulate_column1
from units.separators.column2.fake_column2 import simulate_column2
from units.separators.column3.fake_column3 import simulate_column3
from units.reactors.swing import (
    DesignVars as SwingDesign, FeedStream as SwingFeed,
    FixedParams as SwingFixed, simulate_swing_reactor_system,
)
from units.separators.psa.psa_system import (
    PSADesignVars, PSAFeedStream, PSAFixedParams, simulate_psa_system,
)
from units.separators.membrane.membrane_system import (
    MemDesignVars, MemFeedStream, MemFixedParams, simulate_membrane_system,
)
from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
    CATALYST_PTSN_JPY_PER_KG, CATALYST_PTSN_LIFE_YEARS,
    OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS,
)

# ===========================================================================
# 設計変数（最適化対象）
# ===========================================================================

# 原料 [kmol/h]
F_C3H8_FEED  = 100.0
F_C4H10_FEED =  10.0

# スイング反応器
# T_in=900K でシングルパス転化率 ~22%。リサイクル系で C3H6 が蓄積するのを
# 「反応器で吐き出す」だけの転化率が無いと膜分離を通り抜けて閉ループになる。
SWING = SwingDesign(T_in=900.0, z_cat=5.0, t_cyc=15.0, D=2.0)

# PSA
PSA = PSADesignVars(D_col=0.7, L_bed=6.0, desorption_target=0.35)

# 膜分離 — リサイクル系では膜フィード流量が単パスの 10 倍程度に膨らむため
# A_mem=1000 では stage_cut が 5% 以下に低下し C3H6 が retentate 経由で
# 反応器に戻り蓄積する。A_mem=10000 で stage_cut ~30% 以上を確保する。
MEM = MemDesignVars(P_H=25.0e5, P_L=1.0e5, A_mem=10000.0, P_dist=20.0e5)

# 膜入口冷却目標温度（Dist2 塔底をここまで冷却して液相にする）
T_MEM_FEED = 323.15  # [K] = 50°C

# ===========================================================================
# リサイクル収束パラメータ
# ===========================================================================

MAX_ITER            = 200     # 最大反復
TOL_F               = 0.10    # 収束判定 [kmol/h]
RELAX_ALPHA         = 0.30    # アンダーリラックス係数 (0<α≤1, 小さいほど安定)
RECYCLE_GUARD_RATIO = 15.0    # tear_mem.A が Fresh × この値超で暴走打ち切り

ZERO = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'D': 0.0, 'E': 0.0, 'F': 0.0, 'Z': 0.0}

# ===========================================================================
# ユーティリティ
# ===========================================================================

def hdr(title):
    print(f"\n{'='*64}")
    print(f"  {title}")
    print('='*64)

def show_stream(label, stream):
    comp_names = {'A':'C3H8','B':'C3H6','C':'H2','D':'C2H4','E':'CH4','F':'C2H6','Z':'C4H10'}
    parts = [f"{comp_names.get(k,k)}:{v:.1f}" for k, v in sorted(stream.F_in.items()) if v > 0.01]
    print(f"  {label}: {', '.join(parts)}")
    print(f"  {' '*len(label)}  T={stream.T_in-273.15:.0f}°C  P={stream.P_in/1e5:.1f}bar  F={stream.total_flow():.1f}kmol/h")

# ===========================================================================
# 1 周回シミュレーション
# ===========================================================================

def run_one_pass(tear_dist3, tear_mem, T_d3, T_mem):
    """
    Parameters
    ----------
    tear_dist3 : dict {'A','B'} [kmol/h]  Dist3 塔底由来のリサイクル
    tear_mem   : dict {'A','B'} [kmol/h]  Mem 保留側由来のリサイクル
    T_d3, T_mem: float [K]                各リサイクルの温度（前ループの値）

    Returns
    -------
    dict  全ユニット結果＋ tear_*_new ＋ T_*_new
    """
    # ---- Fresh LPG (1 atm, 25°C) ----
    fresh = ProcessStream(
        F_in={'A': F_C3H8_FEED, 'Z': F_C4H10_FEED,
              'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=298.15, P_in=101325.,
    )

    # ---- Step 1: Comp1 → Dist1 (Fresh のみ) ----
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp1 = simulate_compressor(fresh, P_out_target=17.0e5)
        r1 = simulate_column1(comp1.outlet)

    # 塔頂を 1 atm に膨張（C4 除去済みの C3 主成分）
    dist1_top_1atm = ProcessStream(
        F_in=dict(r1.top.F_in), T_in=r1.top.T_in, P_in=101325.,
    )

    # ---- リサイクルストリーム（1 atm 膨張弁経由、コストなし）----
    recycle_dist3 = ProcessStream(
        F_in={**ZERO, 'A': tear_dist3['A'], 'B': tear_dist3['B']},
        T_in=T_d3, P_in=101325.,
    )
    recycle_mem = ProcessStream(
        F_in={**ZERO, 'A': tear_mem['A'], 'B': tear_mem['B']},
        T_in=T_mem, P_in=101325.,
    )

    # ---- Reactor 入口で合流 ----
    reactor_inlet = mix_streams([dist1_top_1atm, recycle_dist3, recycle_mem])

    # ---- Step 2: Swing Reactor ----
    swing_feed = SwingFeed(
        F_in=reactor_inlet.F_in,
        T_feed=reactor_inlet.T_in,
        P_in=reactor_inlet.P_in,
    )
    r_rx = simulate_swing_reactor_system(SWING, swing_feed, SwingFixed())

    rx_out = ProcessStream(
        F_in=r_rx.effluent.F_out_avg,
        T_in=r_rx.effluent.T_out_avg,
        P_in=r_rx.effluent.P_out,
    )

    # ---- Step 3: Cooler → Comp2 → Dist2 ----
    cooled = simulate_cooler(rx_out, T_out_target=320.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp2 = simulate_compressor(cooled.outlet, P_out_target=20.0e5)
        r2 = simulate_column2(comp2.outlet)

    # ---- Step 4: PSA ----
    psa_feed = PSAFeedStream(
        F_in=r2.top.F_in, T_in=r2.top.T_in, P_in=r2.top.P_in,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_psa = simulate_psa_system(PSA, psa_feed, PSAFixedParams())

    # ---- Step 5: Membrane ----
    mem_precool = simulate_cooler(r2.bottom, T_out_target=T_MEM_FEED)
    mem_feed = MemFeedStream(
        F_C3H6=mem_precool.outlet.F_in.get('B', 0.),
        F_C3H8=mem_precool.outlet.F_in.get('A', 0.),
        T_in=mem_precool.outlet.T_in,
        P_in=mem_precool.outlet.P_in,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_mem = simulate_membrane_system(MEM, mem_feed, MemFixedParams())

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


# ===========================================================================
# 収束ループ
# ===========================================================================

hdr("リサイクル収束ループ (逐次置換 + アンダーリラックス)")
print(f"  TOL_F = {TOL_F} kmol/h    RELAX_ALPHA = {RELAX_ALPHA}    MAX_ITER = {MAX_ITER}")
print()
print(f"  iter | tear_d3 (A,B) | tear_mem (A,B) | Δ_max  | Reactor転化率")
print(f"  -----+---------------+----------------+--------+--------------")

tear_dist3 = {'A': 0.0, 'B': 0.0}
tear_mem   = {'A': 0.0, 'B': 0.0}
T_d3       = 333.15
T_mem      = 323.15

results       = None
converged     = False
penalty_hit   = False
guard_hit     = False

for it in range(1, MAX_ITER + 1):
    results = run_one_pass(tear_dist3, tear_mem, T_d3, T_mem)

    # ペナルティガード（PSA/Mem が無効解を返した時）
    if (results['r_psa'].equipment.CAPEX_total >= 1e8 or
        results['r_mem'].equipment.CAPEX_total >= 1e8):
        print(f"  {it:4d} | --- PSA/Mem ペナルティ発火 → 設計変数の見直しが必要 ---")
        penalty_hit = True
        break

    diff_d3  = max(abs(results['tear_dist3_new'][k] - tear_dist3[k]) for k in tear_dist3)
    diff_mem = max(abs(results['tear_mem_new'][k]   - tear_mem[k])   for k in tear_mem)
    diff = max(diff_d3, diff_mem)

    conv = results['r_rx'].performance.Conversion
    print(f"  {it:4d} | {results['tear_dist3_new']['A']:5.2f},{results['tear_dist3_new']['B']:5.2f}"
          f"   | {results['tear_mem_new']['A']:6.2f},{results['tear_mem_new']['B']:5.2f}"
          f"   | {diff:6.3f} | {conv:5.1f}%")

    # 暴走ガード
    if results['tear_mem_new']['A'] > F_C3H8_FEED * RECYCLE_GUARD_RATIO:
        print(f"  → リサイクル暴走ガード発火 (tear_mem.A > {F_C3H8_FEED * RECYCLE_GUARD_RATIO} kmol/h)")
        guard_hit = True
        break

    if diff < TOL_F:
        converged = True
        print(f"  → 収束 (Δ={diff:.4f} < TOL_F={TOL_F})")
        break

    # アンダーリラックス更新
    tear_dist3 = {k: RELAX_ALPHA * results['tear_dist3_new'][k]
                     + (1 - RELAX_ALPHA) * tear_dist3[k] for k in tear_dist3}
    tear_mem   = {k: RELAX_ALPHA * results['tear_mem_new'][k]
                     + (1 - RELAX_ALPHA) * tear_mem[k] for k in tear_mem}
    T_d3  = RELAX_ALPHA * results['T_d3_new']  + (1 - RELAX_ALPHA) * T_d3
    T_mem = RELAX_ALPHA * results['T_mem_new'] + (1 - RELAX_ALPHA) * T_mem

if not (converged or penalty_hit or guard_hit):
    print(f"  → 未収束 ({MAX_ITER} 回打ち切り、最終状態で集計)")

if penalty_hit:
    print("\n設計NG。設計変数を見直してから再実行してください。")
    sys.exit(1)


# ===========================================================================
# CAPEX/OPEX/TAC 集計
# ===========================================================================

def collect_capex_opex(R):
    capex = {}
    capex['Comp1']       = R['comp1'].equipment.CAPEX
    capex['Dist1']       = R['r1'].equipment.CAPEX
    capex['Reactor']     = R['r_rx'].equipment.Reactor_CAPEX
    capex['Cooler']      = R['cooled'].equipment.CAPEX
    capex['Comp2']       = R['comp2'].equipment.CAPEX
    capex['Dist2']       = R['r2'].equipment.CAPEX
    capex['PSA容器']     = R['r_psa'].equipment.CAPEX_vessels
    capex['PSA活性炭']   = R['r_psa'].equipment.CAPEX_adsorbent
    capex['MemPrecool']  = R['mem_precool'].equipment.CAPEX
    capex['Mem気化器']   = R['r_mem'].equipment.CAPEX_vap
    capex['Mem F圧縮機'] = R['r_mem'].equipment.CAPEX_comp_feed
    capex['Mem P圧縮機'] = R['r_mem'].equipment.CAPEX_comp_prod
    capex['Mem冷却器']   = R['r_mem'].equipment.CAPEX_cond
    capex['Mem膜本体']   = R['r_mem'].equipment.CAPEX_mem
    capex['Dist3']       = R['r3'].equipment.CAPEX

    # OPEX 換算ヘルパー
    def ele(W_kW):
        return W_kW * ELECTRICITY_JPY_PER_KWH * OPERATING_HOURS_PER_YEAR / 1.0e8
    def heat(Q_kW, jpy_per_GJ):
        # 1 kW × 1 h = 3.6 MJ = 3.6e-3 GJ
        return Q_kW * 3.6e-3 * OPERATING_HOURS_PER_YEAR * jpy_per_GJ / 1.0e8

    opex = {}
    opex['Comp1電力']        = ele(R['comp1'].equipment.W_kW)
    opex['Comp2電力']        = ele(R['comp2'].equipment.W_kW)
    opex['MemF圧縮機電力']   = ele(R['r_mem'].equipment.W_feed_kW)
    opex['MemP圧縮機電力']   = ele(R['r_mem'].equipment.W_prod_kW)
    opex['Dist1リボイラ蒸気']= heat(R['r1'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ)
    opex['Dist2リボイラ蒸気']= heat(R['r2'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ)
    opex['Dist3リボイラ蒸気']= heat(R['r3'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ)
    opex['Mem気化器蒸気']    = heat(R['r_mem'].equipment.Q_vap_kW, LP_STEAM_JPY_PER_GJ)

    # 反応器プリヒーター（GJ/h → kW 換算してから heat()）
    Q_preheat_kW = R['r_rx'].effluent.Q_preheat * 1.0e9 / 3600.0 / 1000.0
    opex['Reactor予熱燃料']  = heat(Q_preheat_kW, FUEL_JPY_PER_GJ)

    opex['Cooler冷水']       = heat(abs(R['cooled'].equipment.Q_duty_kW), COOLING_WATER_JPY_PER_GJ)
    opex['Dist1コンデンサ冷水']= heat(R['r1'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['Dist2コンデンサ冷水']= heat(R['r2'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['Dist3コンデンサ冷水']= heat(R['r3'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['MemPrecool冷水']   = heat(abs(R['mem_precool'].equipment.Q_duty_kW), COOLING_WATER_JPY_PER_GJ)
    opex['Mem冷却器冷水']    = heat(R['r_mem'].equipment.Q_cond_kW, COOLING_WATER_JPY_PER_GJ)

    opex['Reactor触媒交換']  = (R['r_rx'].equipment.Catalyst_Weight_Total
                                * CATALYST_PTSN_JPY_PER_KG
                                / CATALYST_PTSN_LIFE_YEARS / 1.0e8)
    opex['PSA活性炭交換']    = R['r_psa'].equipment.OPEX_adsorbent_okuyen_per_year

    return capex, opex


capex, opex = collect_capex_opex(results)


# ===========================================================================
# 結果サマリ表示
# ===========================================================================

hdr("収束時のフロー一覧")
fresh = ProcessStream(
    F_in={'A': F_C3H8_FEED, 'Z': F_C4H10_FEED},
    T_in=298.15, P_in=101325.,
)
show_stream("Fresh LPG", fresh)
show_stream("Dist1 塔頂 (1atm 膨張後)", results['dist1_top_1atm'])
recycle_total_F = (results['tear_dist3_new']['A'] + results['tear_dist3_new']['B']
                  + results['tear_mem_new']['A']  + results['tear_mem_new']['B'])
print(f"  Recycle 合計 : Mem(A={results['tear_mem_new']['A']:.2f},B={results['tear_mem_new']['B']:.2f})"
      f" + Dist3(A={results['tear_dist3_new']['A']:.2f},B={results['tear_dist3_new']['B']:.2f})"
      f"  = {recycle_total_F:.2f} kmol/h")
show_stream("Reactor 入口 (Fresh + Recycle)", results['reactor_inlet'])
show_stream("Reactor 出口", results['rx_out'])
show_stream("Dist2 塔頂 (→ PSA)",      results['r2'].top)
show_stream("Dist2 塔底 (→ Mem)",      results['r2'].bottom)
show_stream("Membrane 透過 (→ Dist3)", ProcessStream(
    F_in={'A': results['r_mem'].product.F_C3H8, 'B': results['r_mem'].product.F_C3H6,
          'C':0.,'D':0.,'E':0.,'F':0.},
    T_in=results['r_mem'].product.T_out, P_in=results['r_mem'].product.P_out))
show_stream("Dist3 塔頂 (C3H6 製品)",   results['r3'].top)


hdr("生産・収率（収束時）")
C3H6_product = results['r3'].top.F_in.get('B', 0.0)
H2_product   = results['r_psa'].product.get('C', 0.0)
yield_pct    = C3H6_product / F_C3H8_FEED * 100.0
recycle_C3H8 = results['tear_dist3_new']['A'] + results['tear_mem_new']['A']
recycle_C3H6 = results['tear_dist3_new']['B'] + results['tear_mem_new']['B']

print(f"  Fresh C3H8     : {F_C3H8_FEED:7.2f} kmol/h")
print(f"  Fresh C4H10    : {F_C4H10_FEED:7.2f} kmol/h")
print(f"  C3H6 製品      : {C3H6_product:7.2f} kmol/h   (収率 {yield_pct:5.1f}%)")
print(f"  H2 副産物      : {H2_product:7.2f} kmol/h")
print(f"  Recycle C3H8   : {recycle_C3H8:7.2f} kmol/h   (Fresh比 {recycle_C3H8/F_C3H8_FEED*100:.0f}%)")
print(f"  Recycle C3H6   : {recycle_C3H6:7.2f} kmol/h")
print(f"  Reactor 入口   : {sum(results['reactor_inlet'].F_in.values()):7.2f} kmol/h"
      f"  (Fresh の {sum(results['reactor_inlet'].F_in.values())/F_C3H8_FEED:.2f} 倍)")
print(f"  Reactor 転化率 : {results['r_rx'].performance.Conversion:5.1f}%")
print(f"  Reactor 選択率 : {results['r_rx'].performance.Selectivity:5.1f}%")


hdr("CAPEX 内訳 [億円]")
for n, v in capex.items():
    if v < 1e6:
        print(f"  {n:<14}: {v:8.4f}")
    else:
        print(f"  {n:<14}:   ペナルティ")
total_capex = sum(v for v in capex.values() if v < 1e6)
print(f"  {'-'*26}")
print(f"  {'合計':<14}: {total_capex:8.4f}")


hdr("OPEX 内訳 [億円/年]")
for n, v in opex.items():
    print(f"  {n:<24}: {v:8.4f}")
total_opex = sum(opex.values())
print(f"  {'-'*36}")
print(f"  {'合計':<24}: {total_opex:8.4f}")


hdr("TAC（年間総費用）")
TAC = total_capex / DEPRECIATION_YEARS + total_opex
print(f"  CAPEX/{DEPRECIATION_YEARS}年(償却) : {total_capex/DEPRECIATION_YEARS:8.4f}  億円/年")
print(f"  OPEX 合計          : {total_opex:8.4f}  億円/年")
print(f"  ──────────────────────────────")
print(f"  TAC                : {TAC:8.4f}  億円/年")
print()
if C3H6_product > 0:
    # C3H6 分子量 42.08 g/mol = 0.04208 kg/mol = 42.08 kg/kmol
    annual_kg = C3H6_product * OPERATING_HOURS_PER_YEAR * 42.08
    unit_jpy_per_t = TAC * 1.0e8 / (annual_kg / 1000.0)
    print(f"  C3H6 年間生産量    : {annual_kg/1000.0:.0f}  ton/年")
    print(f"  C3H6 製品単価      : {unit_jpy_per_t:.0f}  円/ton  ({unit_jpy_per_t/1000:.1f} 円/kg)")
print()
print(f"  仮ユーティリティ単価:")
print(f"    電力 {ELECTRICITY_JPY_PER_KWH}円/kWh  /  LP蒸気 {LP_STEAM_JPY_PER_GJ}円/GJ")
print(f"    冷水 {COOLING_WATER_JPY_PER_GJ}円/GJ  /  燃料 {FUEL_JPY_PER_GJ}円/GJ")
print(f"    PtSn触媒 {CATALYST_PTSN_JPY_PER_KG:.0f}円/kg / 寿命{CATALYST_PTSN_LIFE_YEARS:.0f}年")
print(f"    稼働 {OPERATING_HOURS_PER_YEAR:.0f}h/年  /  償却 {DEPRECIATION_YEARS}年")
print()
print(f"  ★ 仮値は src/cost_parameters.py に集約。コンテスト課題 Ver.2.0 のサイト仕様に置換のこと。")
print(f"  ★ 蒸留塔は fake_columnX (split_fracs ベースのダミー)。正式 VLE モデル実装後に置換のこと。")
