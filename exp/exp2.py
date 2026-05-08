"""
exp2.py — リサイクルあり PDH プロセス全体フロー シミュレーション

このスクリプトは flowsheet.evaluate() を呼ぶ薄いランナー。
構成:
  - 運転条件 (固定値)  → config/operating.toml
  - 設計変数 (最適化対象) → 下記 SWING/PSA/MEM
  - 物理計算・収束・経済計算 → flowsheet/* に移譲

最適化器を組むときは flowsheet.evaluate(design, config) を直接呼ぶ。

リサイクル構成:
  - Membrane 保留側 (C3H8 富化, 残留 C3H6 含)        ─┐
  - Dist3 塔底       (未透過 C3H8, 残留 C3H6 含)     ─┴→ Reactor 直前で合流
  Dist1 (脱ブタン) には戻さない (軽質ガス・C4 は既に系外へ抜けているため)

差し替えポイント:
  1. 蒸留塔モデル (fake_column1/2/3 → 正式 VLE モデル):
       同一インターフェイス (ProcessStream → DistResult) なら import 行のみ変更
  2. ユーティリティ単価:
       src/cost_parameters.py の定数を実値に置換
"""

import os
import sys

# Windows コンソール (cp932) でも記号を表示
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from stream.stream import ProcessStream
from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate

from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
    CATALYST_PTSN_JPY_PER_KG, CATALYST_PTSN_LIFE_YEARS,
    OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS,
)

# ===========================================================================
# 設計変数 (最適化対象)
# ===========================================================================

# スイング反応器
# 物質収支の定常解 Fresh ≈ S × (0.02 + 0.98X) より、リサイクル系では反応器入口
# C3H8 流量 S が Fresh の 3〜4 倍に達する (X=22% で S≈6200, X=15% で S≈8700)。
# 設計判断 (2026-05-08): contest §3-3-2 準拠で反応器圧力を 0.5 bar に変更したため、
# 速度式の P_C3H8 が 1/2 に落ち X が大幅低下 (1atm: ~14% → 0.5bar: ~5%)、リサイクル
# 暴走に至る。同等 X を維持するため触媒量と入口温度を引き上げ:
#   - T_in: 900 → 950 K (k_1 を増やして速度補償、副反応 cracking も増えるが T 上限 700°C 内)
#   - z_cat: 15 → 30 m³ (触媒体積倍増で接触時間倍増)
#   - D: 5 → 7 m (断面積拡大で u_0 を抑え圧損を回避、S~10000 kmol/h を捌く)
SWING = SwingDesign(T_in=950.0, z_cat=30.0, t_cyc=15.0, D=7.0)

# PSA — exp2 はスケールが exp1 の 15 倍 + リサイクル系のため副生 H2 も増え、
# PSA フィード非C3 流量が exp1 から 50 倍以上に膨らむ。t_abs_css ≥ 60s を満たす
# には容積を相応に拡大する必要がある (容積 ∝ t_abs)。
PSA = PSADesignVars(D_col=3.0, L_bed=20.0, desorption_target=0.35)

# 膜分離 — A_mem=10000 では C3H6 透過量が ~120 kmol/h で頭打ちになり、リサイクル
# 系の物質収支 (反応器で生成する C3H6 を毎時系外に抜く必要量 ≈ 1450 kmol/h) を
# 満たせず C3H8 が保留側に蓄積して発散する。
# P_H: Hua et al. (2024) で検証された圧力範囲上限 9.5 bar に合わせる。
# 但し simulate_membrane_system は P_H > feed.P_in を要求するため、
# Dist2 P_col (= mem feed P) を 8.5 bar に下げて 1 bar の差を確保。
MEM = MemDesignVars(P_H=9.5e5, P_L=1.0e5, A_mem=100000.0, P_dist=20.0e5)


# ===========================================================================
# 表示ヘルパ
# ===========================================================================

def hdr(title):
    print(f"\n{'='*64}")
    print(f"  {title}")
    print('='*64)


def show_stream(label, stream):
    comp_names = {'A':'C3H8','B':'C3H6','C':'H2','D':'C2H4','E':'CH4','F':'C2H6','Z':'C4H10'}
    parts = [f"{comp_names.get(k,k)}:{v:.1f}" for k, v in sorted(stream.F_in.items()) if v > 0.01]
    print(f"  {label}: {', '.join(parts)}")
    print(f"  {' '*len(label)}  T={stream.T_in-273.15:.0f}°C  P={stream.P_in/1e5:.1f}bar"
          f"  F={stream.total_flow():.1f}kmol/h")


# ===========================================================================
# 実行
# ===========================================================================

config = load_operating_config()
design = FlowsheetDesignVars(swing=SWING, psa=PSA, mem=MEM)

hdr("外側ループ: 製品流量厳密化 (Fresh を調整)")
result = evaluate(design, config, verbose=True)

# solver-level 失敗 (economics=None) は結果が信頼できないため早期終了。
# spec 違反 (economics 計算済み) は最後まで走らせて TAC + ペナルティ内訳を表示する。
if result.economics is None:
    print(f"\n設計NG: {result.failure_reason}")
    print(f"  effective_TAC = {result.effective_TAC:.0f} 億円/年 (固定打ち切り値)")
    print("  設計変数を見直してから再実行してください。")
    sys.exit(1)


# ===========================================================================
# 結果サマリ表示
# ===========================================================================

R            = result.solver.one_pass
F_C3H8_feed  = result.solver.fresh_C3H8
F_C4H10_feed = result.solver.fresh_C4H10

hdr("収束時のフロー一覧")
fresh = ProcessStream(
    F_in={'A': F_C3H8_feed, 'Z': F_C4H10_feed},
    T_in=config.feed.T_K, P_in=config.feed.P_Pa,
)
show_stream("Fresh LPG", fresh)
show_stream("Dist1 塔頂 (0.5bar 膨張後)", R['dist1_top_rx'])

tear_d3  = R['tear_dist3_new']
tear_mem = R['tear_mem_new']
recycle_total_F = tear_d3['A'] + tear_d3['B'] + tear_mem['A'] + tear_mem['B']
print(f"  Recycle 合計 : Mem(A={tear_mem['A']:.2f},B={tear_mem['B']:.2f})"
      f" + Dist3(A={tear_d3['A']:.2f},B={tear_d3['B']:.2f})"
      f"  = {recycle_total_F:.2f} kmol/h")

show_stream("Reactor 入口 (Fresh + Recycle)", R['reactor_inlet'])
show_stream("Reactor 出口", R['rx_out'])
show_stream("Dist2 塔頂 (→ PSA)",      R['r2'].top)
show_stream("Dist2 塔底 (→ Mem)",      R['r2'].bottom)
show_stream("Membrane 透過 (→ Dist3)", ProcessStream(
    F_in={'A': R['r_mem'].product.F_C3H8, 'B': R['r_mem'].product.F_C3H6,
          'C':0.,'D':0.,'E':0.,'F':0.},
    T_in=R['r_mem'].product.T_out, P_in=R['r_mem'].product.P_out))
show_stream("Dist3 塔頂 (C3H6 製品)", R['r3'].top)


hdr("生産・収率(収束時)")
target_kmol_h = (config.product.target_mta * 1000.0
                 / config.product.mw_kg_per_kmol / OPERATING_HOURS_PER_YEAR)
C3H6_product  = R['r3'].top.F_in.get('B', 0.0)
H2_product    = R['r_psa'].product.get('C', 0.0)
yield_pct     = C3H6_product / F_C3H8_feed * 100.0
recycle_C3H8  = tear_d3['A'] + tear_mem['A']
recycle_C3H6  = tear_d3['B'] + tear_mem['B']

print(f"  C3H6 目標      : {target_kmol_h:7.2f} kmol/h"
      f"  ({config.product.target_mta:.0f} t/年 @ {OPERATING_HOURS_PER_YEAR:.0f} h/年)")
print(f"  Fresh C3H8     : {F_C3H8_feed:7.2f} kmol/h   (外側ループで決定)")
print(f"  Fresh C4H10    : {F_C4H10_feed:7.2f} kmol/h"
      f"   (LPG 組成 C3H8={config.feed.lpg_c3h8_mol_fraction} より)")
print(f"  C3H6 製品      : {C3H6_product:7.2f} kmol/h"
      f"   (実収率 {yield_pct:5.2f}%, 仮定 {config.feed.yield_assumed*100:.1f}%)")
print(f"  目標との差     : {C3H6_product - target_kmol_h:+7.3f} kmol/h")
print(f"  H2 副産物      : {H2_product:7.2f} kmol/h")
print(f"  Recycle C3H8   : {recycle_C3H8:7.2f} kmol/h"
      f"   (Fresh比 {recycle_C3H8/F_C3H8_feed*100:.0f}%)")
print(f"  Recycle C3H6   : {recycle_C3H6:7.2f} kmol/h")
print(f"  Reactor 入口   : {sum(R['reactor_inlet'].F_in.values()):7.2f} kmol/h"
      f"  (Fresh の {sum(R['reactor_inlet'].F_in.values())/F_C3H8_feed:.2f} 倍)")
print(f"  Reactor 転化率 : {R['r_rx'].performance.Conversion:5.1f}%")
print(f"  Reactor 選択率 : {R['r_rx'].performance.Selectivity:5.1f}%")


hdr("CAPEX 内訳 [億円]")
for n, v in result.economics.capex.items():
    if v < 1e6:
        print(f"  {n:<14}: {v:8.4f}")
    else:
        print(f"  {n:<14}:   ペナルティ")
print(f"  {'-'*26}")
print(f"  {'合計':<14}: {result.economics.total_capex:8.4f}")


hdr("OPEX 内訳 [億円/年]")
for n, v in result.economics.opex.items():
    print(f"  {n:<24}: {v:8.4f}")
print(f"  {'-'*36}")
print(f"  {'合計':<24}: {result.economics.total_opex:8.4f}")


hdr("製品仕様 compliance")
specs = result.specs
print(f"  C3H6 純度 (Dist3 塔頂) : {specs.c3h6_purity_wtfrac*100:6.3f} wt%"
      f"   {'✓' if specs.c3h6_pass else '✗'}"
      f" (spec ≥ {99.5:.1f} wt%)")
print(f"  H2 純度 (PSA 製品)     : {specs.h2_purity_molfrac*100:6.3f} mol%"
      f"   {'✓' if specs.h2_pass else '✗'}"
      f" (spec ≥ {99.9:.1f} mol%)")
print(f"  生産量                 : {specs.production_kmol_h:7.2f} kmol/h"
      f" {'✓' if specs.production_pass else '✗'}"
      f" (片側 spec ≥ {specs.target_kmol_h * 0.99:.2f})")
if not specs.all_pass:
    print(f"  違反内訳 : {result.failure_reason}")


hdr("TAC(年間総費用)")
TAC = result.economics.TAC
print(f"  CAPEX/{DEPRECIATION_YEARS}年(償却) : {result.economics.total_capex/DEPRECIATION_YEARS:8.4f}  億円/年")
print(f"  OPEX 合計          : {result.economics.total_opex:8.4f}  億円/年")
print(f"  ──────────────────────────────")
print(f"  TAC (実コスト)     : {TAC:8.4f}  億円/年")
# 設計判断: 最適化器は effective_TAC を最小化対象にする。
# spec 違反時はソフトペナルティが上乗せされ、違反 0 のときは TAC と一致する。
penalty_amount = result.effective_TAC - TAC
if penalty_amount > 0:
    print(f"  + spec違反ペナルティ : {penalty_amount:8.4f}  億円/年")
print(f"  ──────────────────────────────")
print(f"  effective_TAC      : {result.effective_TAC:8.4f}  億円/年  (← 最適化器の目的関数)")
print()
if C3H6_product > 0:
    print(f"  C3H6 年間生産量    : {result.economics.annual_kg_C3H6/1000.0:.0f}  ton/年")
    print(f"  C3H6 製品単価      : {result.economics.unit_jpy_per_t:.0f}  円/ton"
          f"  ({result.economics.unit_jpy_per_t/1000:.1f} 円/kg)")
print()
print(f"  仮ユーティリティ単価:")
print(f"    電力 {ELECTRICITY_JPY_PER_KWH}円/kWh  /  LP蒸気 {LP_STEAM_JPY_PER_GJ}円/GJ")
print(f"    冷水 {COOLING_WATER_JPY_PER_GJ}円/GJ  /  燃料 {FUEL_JPY_PER_GJ}円/GJ")
print(f"    PtSn触媒 {CATALYST_PTSN_JPY_PER_KG:.0f}円/kg / 寿命{CATALYST_PTSN_LIFE_YEARS:.0f}年")
print(f"    稼働 {OPERATING_HOURS_PER_YEAR:.0f}h/年  /  償却 {DEPRECIATION_YEARS}年")
print()
print(f"  ★ 仮値は src/cost_parameters.py に集約。コンテスト課題 Ver.2.0 のサイト仕様に置換のこと。")
print(f"  ★ 蒸留塔は fake_columnX (split_fracs ベースのダミー)。正式 VLE モデル実装後に置換のこと。")
