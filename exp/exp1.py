"""
exp1.py — PDH プロセス全体フロー シングルパス確認実験

目的:
    実装済みの全ユニットを一本のパイプラインとして接続し、
    物質・エネルギー・コストの流れが端から端まで通ることを確認する。

全体フロー (contest §3-3-2/§3-3-3 準拠):
    Fresh LPG (303K, ~9.97bar 飽和液、C3H8:C4H10 = 9:1)
     → [Pump1]   17bar へ液送 (contest §3-3-3)
     → [Dist1]   脱ブタン塔 (C4H10 を除去)
     → [膨張弁]  0.5bar へ減圧 (contest §3-3-2 反応器圧力)
     → [Reactor] PDH スイング反応 (600°C, 0.5bar)
     → [Cooler]  320K へ冷却
     → [Comp2a]  圧縮比 √17 で 1段目
     → [Intercool] 段間冷却 320K
     → [Comp2b]  8.5bar へ 2段目
     → [Dist2]   脱エタン (H2/CH4/C2H4 を頂部に分離)
          塔頂 → [PSA]      H2 製品 + オフガス(燃料)
          塔底 → [Membrane] C3H6/C3H8 分離
                   透過側  → [Dist3]  C3H6 製品精製
                   保留側  → リサイクル (今回は廃棄扱い)

注意:
    リサイクルなしのシングルパスのため、実際の C3H6 収率はこれより高い。
    リサイクル収束ループは別の実験（exp2.py 予定）で実施する。
"""

import math
import os
import sys
import warnings

# Windows コンソール (cp932) でも日本語・記号を表示するため stdout を UTF-8 に
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from stream.stream import ProcessStream
from units.utils.cooler import simulate_cooler
from units.utils.compressor import simulate_compressor
from units.utils.pump import simulate_pump
from units.separators.column1.column1 import simulate_column1
from units.separators.column2.column2 import simulate_column2
from units.separators.column3.column3 import simulate_column3
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

# ===========================================================================
# 設計変数（ここを変えて実験する）
# ===========================================================================

# LPG 原料 (30°C 飽和液、C3H8:C4H10 = 9:1)
F_C3H8_FEED  = 100.0   # プロパン [kmol/h]
F_C4H10_FEED =  10.0   # ブタン (LPG 不純物) [kmol/h]

# スイング反応器
SWING = SwingDesign(T_in=873.15, z_cat=5.0, t_cyc=15.0, D=2.0)

# PSA — 非C3 流量が小さい (~14 kmol/h) ため塔径も小さく取る。D=1.5 だと u_0 が
# 過小で破過に達せずペナルティ条件に該当する。D=0.7, L=6.0 で H2 回収率 ~84%。
PSA = PSADesignVars(D_col=0.7, L_bed=6.0, desorption_target=0.35)

# Membrane — P_dist は Dist3 P_col (=20 bar) に揃える (Mem→Dist3 直結のため)。
# C3H6 ~99wt% の透過側を冷却水で凝縮するには P_dist >= ~17 bar が必要 (low-P 化すると
# 泡点が冷却水出口温度 40°C 以下になり温度クロス → ペナルティ)。
# P_H は Hua et al. (2024) で検証された圧力範囲 ≤9.5 bar に合わせる。
# Dist2 を 8.5 bar 運転にしたため P_H > feed.P_in 制約を 1 bar マージンで満たす 9.5 bar。
# A_mem は 100 kmol/h スケールでフラックスを取るため 1000 m²。
MEM = MemDesignVars(P_H=9.5e5, P_L=1.0e5, A_mem=1000.0, P_dist=20.0e5)

# Membrane 入口冷却 — Dist2 塔底の T が新運転圧力 8.5 bar の泡点を下回る場合があるため
# mem_feed_K (50°C) で気化フィード化する (run_one_pass.py と同じ phase_change=True 想定)。
T_MEM_FEED = 323.15  # [K] = 50°C

# ===========================================================================
# ユーティリティ
# ===========================================================================

def hdr(title):
    print(f"\n{'='*58}")
    print(f"  {title}")
    print('='*58)

def show_stream(label, stream):
    comp_names = {'A':'C3H8','B':'C3H6','C':'H2','D':'C2H4','E':'CH4','F':'C2H6','Z':'C4H10'}
    parts = [f"{comp_names.get(k,k)}:{v:.1f}" for k, v in sorted(stream.F_in.items()) if v > 0.01]
    print(f"  {label}: {', '.join(parts)}")
    print(f"  {'':>{len(label)}}  T={stream.T_in-273.15:.0f}°C  P={stream.P_in/1e5:.1f}bar  F={stream.total_flow():.1f}kmol/h")

capex = {}   # 各ユニットの CAPEX [億円] を集計する辞書

# ===========================================================================
# Step 0: LPG 原料 (30°C 飽和液、C3H8:C4H10 = 9:1)
# ===========================================================================
hdr("Step 0: LPG 原料")

# LPG (9:1) の 30°C bubble point は Raoult 近似で ≈ 9.97 bar
fresh = ProcessStream(
    F_in={'A': F_C3H8_FEED, 'Z': F_C4H10_FEED,
          'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0.},
    T_in=303.15, P_in=9.97e5,
)
show_stream("原料", fresh)

# ===========================================================================
# Step 1: Pump1 + Dist1 (脱ブタン塔)
# ===========================================================================
hdr("Step 1: Pump1 + Dist1 (脱ブタン塔)")

# contest §3-3-3: 加圧すべき箇所には、ポンプ(液) を入れること
# 30°C 飽和液 → 17 bar (Dist1 圧力) は液送ポンプで行う
pump1 = simulate_pump(fresh, P_out_target=17.0e5)

capex['Pump1'] = pump1.equipment.CAPEX
print(f"  Pump1: {pump1.equipment.W_kW:.2f} kW  ρ={pump1.equipment.rho_liq:.0f}kg/m³"
      f"  V={pump1.equipment.V_dot*3600:.2f}m³/h  CAPEX={pump1.equipment.CAPEX:.4f}億円")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r1 = simulate_column1(pump1.outlet)

capex['Dist1'] = r1.equipment.CAPEX
print(f"  Dist1: D={r1.equipment.D_col:.2f}m  H={r1.equipment.H_col:.1f}m  CAPEX={r1.equipment.CAPEX:.3f}億円")
print(f"         Q_cond={r1.equipment.Q_cond:.0f}kW  Q_reb={r1.equipment.Q_reb:.0f}kW")
show_stream("塔頂 (→反応器)", r1.top)
show_stream("塔底 (C4製品)", r1.bottom)

# 塔頂を 0.5 bar (反応器圧力) に膨張 (膨張弁、機器コストなし)
P_RX = 50000.0   # 反応器圧力 [Pa] = 0.5 bar (contest §3-3-2)
dist1_top = ProcessStream(F_in=r1.top.F_in, T_in=r1.top.T_in, P_in=P_RX)

# ===========================================================================
# Step 2: Swing Reactor
# ===========================================================================
hdr("Step 2: Swing Reactor (PDH 600°C, 0.5bar)")

swing_feed = SwingFeed(F_in=dist1_top.F_in, T_feed=dist1_top.T_in, P_in=dist1_top.P_in)
r_rx = simulate_swing_reactor_system(SWING, swing_feed, SwingFixed())

capex['Reactor'] = r_rx.equipment.Reactor_CAPEX
print(f"  転化率: {r_rx.performance.Conversion:.1f}%   選択率: {r_rx.performance.Selectivity:.1f}%")
print(f"  T_out:  {r_rx.effluent.T_out_avg-273.15:.0f}°C  Q_preheat: {r_rx.effluent.Q_preheat:.2f}GJ/h")
print(f"  CAPEX:  {r_rx.equipment.Reactor_CAPEX:.3f}億円  ({r_rx.equipment.N_reactors_total}基)")

rx_out = ProcessStream(F_in=r_rx.effluent.F_out_avg, T_in=r_rx.effluent.T_out_avg, P_in=r_rx.effluent.P_out)
show_stream("Reactor出口", rx_out)

# ===========================================================================
# Step 3: Cooler + Comp2 + Dist2 (脱エタン塔)
# ===========================================================================
hdr("Step 3: Cooler + Comp2 (2段) + Dist2 (脱エタン塔)")

# 反応器出口 (~700°C) を 320K まで冷却してから昇圧
cooled = simulate_cooler(rx_out, T_out_target=320.0)
capex['Cooler'] = cooled.equipment.CAPEX
print(f"  Cooler: Q={cooled.equipment.Q_duty_kW:.0f}kW  A={cooled.equipment.A_est_m2:.0f}m2  CAPEX={cooled.equipment.CAPEX:.3f}億円")

# Comp2 多段化: 0.5bar → 8.5bar = 圧縮比 17:1。等圧縮比 √17≈4.12 で 2段+段間冷却。
P_OUT_FINAL = 8.5e5    # Dist2 圧力
P_MID = math.sqrt(cooled.outlet.P_in * P_OUT_FINAL)
T_INTERCOOL = 320.0    # 段間冷却ターゲット

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    comp2a    = simulate_compressor(cooled.outlet, P_out_target=P_MID)
    intercool = simulate_cooler(comp2a.outlet, T_out_target=T_INTERCOOL)
    comp2b    = simulate_compressor(intercool.outlet, P_out_target=P_OUT_FINAL)

capex['Comp2a']    = comp2a.equipment.CAPEX
capex['Intercool'] = intercool.equipment.CAPEX
capex['Comp2b']    = comp2b.equipment.CAPEX
print(f"  Comp2a:    {comp2a.equipment.W_kW:.0f}kW  T_out={comp2a.equipment.T_out-273.15:.0f}°C"
      f"  P_mid={P_MID/1e5:.2f}bar  CAPEX={comp2a.equipment.CAPEX:.3f}億円")
print(f"  Intercool: Q={intercool.equipment.Q_duty_kW:.0f}kW  A={intercool.equipment.A_est_m2:.0f}m2"
      f"  CAPEX={intercool.equipment.CAPEX:.3f}億円")
print(f"  Comp2b:    {comp2b.equipment.W_kW:.0f}kW  T_out={comp2b.equipment.T_out-273.15:.0f}°C"
      f"  CAPEX={comp2b.equipment.CAPEX:.3f}億円")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r2 = simulate_column2(comp2b.outlet)

capex['Dist2'] = r2.equipment.CAPEX
print(f"  Dist2:  D={r2.equipment.D_col:.2f}m  H={r2.equipment.H_col:.1f}m  CAPEX={r2.equipment.CAPEX:.3f}億円")
print(f"          Q_cond={r2.equipment.Q_cond:.0f}kW  Q_reb={r2.equipment.Q_reb:.0f}kW")
show_stream("塔頂 (→PSA)", r2.top)
show_stream("塔底 (→Membrane)", r2.bottom)

# ===========================================================================
# Step 4: PSA (H2 分離)
# ===========================================================================
hdr("Step 4: PSA (H2 分離)")

psa_feed = PSAFeedStream(F_in=r2.top.F_in, T_in=r2.top.T_in, P_in=r2.top.P_in)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r_psa = simulate_psa_system(PSA, psa_feed, PSAFixedParams())

capex['PSA'] = r_psa.equipment.CAPEX_total
if r_psa.equipment.CAPEX_total < 1e8:
    print(f"  H2 回収率: {r_psa.H2_recovery*100:.1f}%")
    print(f"  H2 製品:   {r_psa.product.get('C', 0):.2f} kmol/h")
    print(f"  オフガス(燃料): {sum(r_psa.offgas.values()):.2f} kmol/h")
    print(f"  CAPEX: {r_psa.equipment.CAPEX_total:.3f}億円")
else:
    print("  PSA: ペナルティ条件 → 設計変数を調整してください")

# ===========================================================================
# Step 5: Membrane (C3H6/C3H8 分離)
# ===========================================================================
hdr("Step 5: Membrane (C3H6/C3H8 分離)")

# Dist2 塔底を mem_feed_K まで気化・過熱してガスフィードで膜へ送る (run_one_pass と同じ)。
# 8.5bar 運転下では塔底液の泡点が低く、潜熱を含めた相変化計算が必要。
mem_precool = simulate_cooler(r2.bottom, T_out_target=T_MEM_FEED, phase_change=True)
capex['MemPrecool'] = mem_precool.equipment.CAPEX
print(f"  Mem前冷却器: Q={mem_precool.equipment.Q_duty_kW:.0f}kW  A={mem_precool.equipment.A_est_m2:.0f}m2  CAPEX={mem_precool.equipment.CAPEX:.3f}億円")

# Dist2 塔底の C3H6(B) と C3H8(A) のみ渡す（C2H6等は別系統、今回は無視）
mem_feed = MemFeedStream(
    F_C3H6 = mem_precool.outlet.F_in.get('B', 0.),
    F_C3H8 = mem_precool.outlet.F_in.get('A', 0.),
    T_in   = mem_precool.outlet.T_in,
    P_in   = mem_precool.outlet.P_in,
)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r_mem = simulate_membrane_system(MEM, mem_feed, MemFixedParams(vapor_feed=True))

capex['Membrane'] = r_mem.equipment.CAPEX_total
print(f"  透過率(stage cut): {r_mem.stage_cut*100:.1f}%")
print(f"  透過純度(C3H6):    {r_mem.perm_purity*100:.1f}%")
print(f"  透過側 (→Dist3):   C3H6={r_mem.product.F_C3H6:.2f}  C3H8={r_mem.product.F_C3H8:.2f} kmol/h")
print(f"  保留側 (→Recycle): C3H6={r_mem.retentate.F_C3H6:.2f}  C3H8={r_mem.retentate.F_C3H8:.2f} kmol/h")
print(f"  CAPEX: {r_mem.equipment.CAPEX_total:.3f}億円")

# ===========================================================================
# Step 6: Dist3 (C3 スプリッター, C3H6 製品精製)
# ===========================================================================
hdr("Step 6: Dist3 (C3 スプリッター)")

# Membrane 透過側を Dist3 フィードに変換
mem_to_dist3 = ProcessStream(
    F_in={'A': r_mem.product.F_C3H8, 'B': r_mem.product.F_C3H6,
          'C': 0., 'D': 0., 'E': 0., 'F': 0.},
    T_in=r_mem.product.T_out,
    P_in=r_mem.product.P_out,
)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r3 = simulate_column3(mem_to_dist3)

capex['Dist3'] = r3.equipment.CAPEX
print(f"  Dist3:  D={r3.equipment.D_col:.2f}m  H={r3.equipment.H_col:.1f}m  CAPEX={r3.equipment.CAPEX:.3f}億円")
print(f"          Q_cond={r3.equipment.Q_cond:.0f}kW  Q_reb={r3.equipment.Q_reb:.0f}kW")
show_stream("塔頂 → C3H6 製品", r3.top)
show_stream("塔底 → C3H8 リサイクル(今回は廃棄扱い)", r3.bottom)

# ===========================================================================
# 総括
# ===========================================================================
hdr("総括")

C3H6_product = r3.top.F_in.get('B', 0.)
H2_product   = r_psa.product.get('C', 0.) if r_psa.equipment.CAPEX_total < 1e8 else 0.
yield_pct    = C3H6_product / F_C3H8_FEED * 100.

print(f"  フィード  C3H8  : {F_C3H8_FEED:.1f} kmol/h")
print(f"  製品      C3H6  : {C3H6_product:.2f} kmol/h  (全体収率 {yield_pct:.1f}%)")
print(f"  副産物    H2    : {H2_product:.2f} kmol/h")
print()
print(f"  CAPEX 内訳 [億円]:")
for name, val in capex.items():
    try:
        print(f"    {name:<12}: {val:.3f}")
    except (ValueError, TypeError):
        print(f"    {name:<12}: ---")
total = sum(v for v in capex.values() if isinstance(v, float) and v < 1e8)
print(f"    {'合計':<12}: {total:.3f}")
print()
print("  ※ リサイクルなしのシングルパス。")
print("    Dist3塔底(C3H8)・Membrane保留側(C3H8)をリサイクルすると収率が向上する。")


# ===========================================================================
# 監査セクション: 各ユニットが返す CAPEX/OPEX 計算用変数の一覧
# ===========================================================================
hdr("CAPEX 計算入力一覧")
# Bare Module Cost 法のサイジングパラメータ:
#   - 圧縮機:      W_kW (流体動力)
#   - 熱交換器:    A_m2 (伝熱面積)
#   - 容器/塔/反応器: V_m3, P_abs_pa, D_m, N
#   - 膜:         A_mem, MEM_UNIT_PRICE_USD_PER_M2 (★仮)
#   - PSA活性炭:   W_adsorbent_kg, ACTIVATED_CARBON_PRICE_USD_PER_KG (★仮)

print(f"  {'装置':<14} {'タイプ':<8} {'サイジング変数':<32} {'値':<14} {'CAPEX[億円]':>10}")
print(f"  {'-'*14} {'-'*8} {'-'*32} {'-'*14} {'-'*10}")
import math as _math
def _row(name, kind, var, val, capex_v):
    print(f"  {name:<14} {kind:<8} {var:<32} {val:<14} {capex_v:>10.3f}")

# Pump1 (液送、contest §3-3-3)
_row('Pump1', 'ポンプ', f'W_kW={pump1.equipment.W_kW:.2f}, P=17bar', '', pump1.equipment.CAPEX)
# Dist1
_dist1_V = _math.pi/4 * r1.equipment.D_col**2 * r1.equipment.H_col
_row('Dist1', '塔', f'V≈{_dist1_V:.1f}m³, D={r1.equipment.D_col:.2f}m, P=17bar', '', r1.equipment.CAPEX)
# Reactor
_row('Reactor', '塔/反応器',
     f'V_vessel={r_rx.equipment.V_vessel_actual:.1f}m³ × {r_rx.equipment.N_reactors_total}基',
     '', r_rx.equipment.Reactor_CAPEX)
# Cooler
_row('Cooler', '熱交', 'A_m2', f'{cooled.equipment.A_est_m2:.1f}', cooled.equipment.CAPEX)
# Comp2 (2段+段間冷却)
_row('Comp2a', '圧縮機', 'W_kW [kW]', f'{comp2a.equipment.W_kW:.1f}', comp2a.equipment.CAPEX)
_row('Intercool', '熱交', 'A_m2', f'{intercool.equipment.A_est_m2:.1f}', intercool.equipment.CAPEX)
_row('Comp2b', '圧縮機', 'W_kW [kW]', f'{comp2b.equipment.W_kW:.1f}', comp2b.equipment.CAPEX)
# Dist2
_dist2_V = _math.pi/4 * r2.equipment.D_col**2 * r2.equipment.H_col
_row('Dist2', '塔', f'V≈{_dist2_V:.2f}m³, D={r2.equipment.D_col:.2f}m, P=20bar', '', r2.equipment.CAPEX)
# PSA
_row('PSA容器', '塔', f'V/塔={_math.pi/4*PSA.D_col**2*PSA.L_bed:.2f}m³ × {r_psa.equipment.N_total_columns}基',
     '', r_psa.equipment.CAPEX_vessels)
_row('PSA活性炭', '材料', f'W={r_psa.equipment.W_adsorbent_kg:.0f}kg × 5USD/kg(★仮)', '',
     r_psa.equipment.CAPEX_adsorbent)
# MemPrecool
_row('MemPrecool', '熱交', 'A_m2', f'{mem_precool.equipment.A_est_m2:.1f}',
     mem_precool.equipment.CAPEX)
# Membrane
_row('Mem気化器', '熱交', 'A_m2', f'{r_mem.equipment.A_vap:.1f}', r_mem.equipment.CAPEX_vap)
_row('Mem F圧縮機', '圧縮機', 'W_kW [kW]', f'{r_mem.equipment.W_feed_kW:.1f}',
     r_mem.equipment.CAPEX_comp_feed)
_row('Mem P圧縮機', '圧縮機', 'W_kW [kW]', f'{r_mem.equipment.W_prod_kW:.1f}',
     r_mem.equipment.CAPEX_comp_prod)
_row('Mem冷却器', '熱交', 'A_m2', f'{r_mem.equipment.A_cond:.1f}', r_mem.equipment.CAPEX_cond)
_row('Mem膜本体', '膜', f'A_mem={MEM.A_mem:.0f}m² × {r_mem.equipment.n_modules}本(★仮単価)', '',
     r_mem.equipment.CAPEX_mem)
# Dist3
_dist3_V = _math.pi/4 * r3.equipment.D_col**2 * r3.equipment.H_col
_row('Dist3', '塔', f'V≈{_dist3_V:.1f}m³, D={r3.equipment.D_col:.2f}m, P=20bar', '', r3.equipment.CAPEX)


hdr("OPEX 計算入力一覧（消費熱量・電力・触媒）")
# 仮ユーティリティ単価 ★ コンテスト課題 Ver.2.0 の値を別途確認後に置き換え
ELECTRICITY_JPY_PER_KWH = 15.0   # ★ 仮 (購入電力)
LP_STEAM_JPY_PER_GJ     = 1800.0 # ★ 仮 (LP Steam 160°C ≈ 3円/kg, λ=2200kJ/kg)
COOLING_WATER_JPY_PER_GJ =  60.0 # ★ 仮 (冷却水 ~0.05円/MJ)
FUEL_JPY_PER_GJ         = 1500.0 # ★ 仮 (LNG 燃料、反応器プリヒーター)
CATALYST_PTSN_JPY_PER_KG = 50000.0  # ★ 仮 (PtSn 触媒)
CATALYST_PTSN_LIFE_YEAR  = 3.0      # ★ 仮 (触媒交換間隔)
OPERATING_HOURS_PER_YEAR = 8000.0   # 稼働時間 [h/年] ≒ 333日

ELE = ELECTRICITY_JPY_PER_KWH * OPERATING_HOURS_PER_YEAR / 1.0e8  # [億円/年/kW]
def _ele(W_kW): return W_kW * ELE
def _heat(Q_kW, jpy_per_GJ):
    # 1 kW × 1 h = 3.6 MJ = 3.6e-3 GJ
    GJ_per_year = Q_kW * 3.6e-3 * OPERATING_HOURS_PER_YEAR
    return GJ_per_year * jpy_per_GJ / 1.0e8

opex = {}  # [億円/年]
print(f"  {'装置':<14} {'費目':<10} {'消費量':<24} {'年間OPEX[億円]':>14}")
print(f"  {'-'*14} {'-'*10} {'-'*24} {'-'*14}")
def _orow(name, kind, qty, val):
    print(f"  {name:<14} {kind:<10} {qty:<24} {val:>14.4f}")
    opex[f'{name}_{kind}'] = val

# 電力
_orow('Pump1', '電力', f'{pump1.equipment.W_kW:.2f}kW', _ele(pump1.equipment.W_kW))
_orow('Comp2a', '電力', f'{comp2a.equipment.W_kW:.1f}kW', _ele(comp2a.equipment.W_kW))
_orow('Comp2b', '電力', f'{comp2b.equipment.W_kW:.1f}kW', _ele(comp2b.equipment.W_kW))
_orow('Mem F圧縮機', '電力', f'{r_mem.equipment.W_feed_kW:.1f}kW', _ele(r_mem.equipment.W_feed_kW))
_orow('Mem P圧縮機', '電力', f'{r_mem.equipment.W_prod_kW:.1f}kW', _ele(r_mem.equipment.W_prod_kW))

# スチーム (リボイラー、Mem気化器)
_orow('Dist1リボイラ', 'LP蒸気', f'{r1.equipment.Q_reb:.0f}kW',
      _heat(r1.equipment.Q_reb, LP_STEAM_JPY_PER_GJ))
_orow('Dist2リボイラ', 'LP蒸気', f'{r2.equipment.Q_reb:.0f}kW',
      _heat(r2.equipment.Q_reb, LP_STEAM_JPY_PER_GJ))
_orow('Dist3リボイラ', 'LP蒸気', f'{r3.equipment.Q_reb:.0f}kW',
      _heat(r3.equipment.Q_reb, LP_STEAM_JPY_PER_GJ))
_orow('Mem気化器', 'LP蒸気', f'{r_mem.equipment.Q_vap_kW:.0f}kW',
      _heat(r_mem.equipment.Q_vap_kW, LP_STEAM_JPY_PER_GJ))

# 燃料 (反応器プリヒーター: 600°C は LP 蒸気では不可、燃焼炉)
_Q_preheat_kW = r_rx.effluent.Q_preheat * 1e9 / 3600.0 / 1000.0  # GJ/h → kW
_orow('Reactor予熱', '燃料', f'{_Q_preheat_kW:.0f}kW',
      _heat(_Q_preheat_kW, FUEL_JPY_PER_GJ))

# 冷却水
_orow('Cooler', '冷水', f'{abs(cooled.equipment.Q_duty_kW):.0f}kW',
      _heat(abs(cooled.equipment.Q_duty_kW), COOLING_WATER_JPY_PER_GJ))
_orow('Intercool', '冷水', f'{abs(intercool.equipment.Q_duty_kW):.0f}kW',
      _heat(abs(intercool.equipment.Q_duty_kW), COOLING_WATER_JPY_PER_GJ))
_orow('Dist1コンデンサ', '冷水', f'{r1.equipment.Q_cond:.0f}kW',
      _heat(r1.equipment.Q_cond, COOLING_WATER_JPY_PER_GJ))
_orow('Dist2コンデンサ', '冷水', f'{r2.equipment.Q_cond:.0f}kW',
      _heat(r2.equipment.Q_cond, COOLING_WATER_JPY_PER_GJ))
_orow('Dist3コンデンサ', '冷水', f'{r3.equipment.Q_cond:.0f}kW',
      _heat(r3.equipment.Q_cond, COOLING_WATER_JPY_PER_GJ))
_orow('MemPrecool', '冷水', f'{abs(mem_precool.equipment.Q_duty_kW):.0f}kW',
      _heat(abs(mem_precool.equipment.Q_duty_kW), COOLING_WATER_JPY_PER_GJ))
_orow('Mem冷却器', '冷水', f'{r_mem.equipment.Q_cond_kW:.0f}kW',
      _heat(r_mem.equipment.Q_cond_kW, COOLING_WATER_JPY_PER_GJ))

# 触媒交換
_cat_kg = r_rx.equipment.Catalyst_Weight_Total
_cat_opex = _cat_kg * CATALYST_PTSN_JPY_PER_KG / CATALYST_PTSN_LIFE_YEAR / 1.0e8
_orow('Reactor触媒', '交換',
      f'{_cat_kg:.0f}kg/{CATALYST_PTSN_LIFE_YEAR:.0f}年(★仮単価)', _cat_opex)
# PSA 活性炭交換 (psa_system が直接出力)
_orow('PSA活性炭', '交換',
      f'{r_psa.equipment.W_adsorbent_kg:.0f}kg/{4.0:.0f}年(★仮)',
      r_psa.equipment.OPEX_adsorbent_okuyen_per_year)


hdr("TAC（年間総費用）")
total_capex = sum(v for v in capex.values() if isinstance(v, float) and v < 1e8)
total_opex  = sum(opex.values())
DEPR = 8  # 償却年数
TAC = total_capex / DEPR + total_opex

print(f"  仮ユーティリティ単価:")
print(f"    電力 {ELECTRICITY_JPY_PER_KWH}円/kWh, LP蒸気 {LP_STEAM_JPY_PER_GJ}円/GJ,")
print(f"    冷水 {COOLING_WATER_JPY_PER_GJ}円/GJ, 燃料 {FUEL_JPY_PER_GJ}円/GJ")
print(f"    稼働時間 {OPERATING_HOURS_PER_YEAR:.0f}h/年, 償却 {DEPR}年")
print()
print(f"  CAPEX 合計        : {total_capex:8.3f}  億円")
print(f"  CAPEX/{DEPR}年(償却) : {total_capex/DEPR:8.3f}  億円/年")
print(f"  OPEX 合計          : {total_opex:8.3f}  億円/年")
print(f"  ──────────────────────────────────")
print(f"  TAC                : {TAC:8.3f}  億円/年")
print()
print(f"  C3H6 製品単価換算 : {TAC*1e8 / max(C3H6_product*OPERATING_HOURS_PER_YEAR*42.08, 1e-6):.0f} 円/kg")
print(f"    (C3H6 分子量 42.08 g/mol, {OPERATING_HOURS_PER_YEAR:.0f}h/年稼働基準)")
print()
print("  ★ 単価は全て仮置き — コンテスト課題 Ver.2.0 のサイト仕様に置き換えること。")
