# -*- coding: utf-8 -*-
r"""_run_lpg_breakeven_catofin.py — lpg_price_breakeven.ipynb の catofin 版を実走して
損益分岐数値を再計算し、レポート用の図 (PNG/PDF) を再生成する。

ノート本体 (lpg_price_breakeven.ipynb) のセル1-5/11/13 と同一ロジックを catofin best
(trial #227, outputs/main_20260604_014318/best.json) で実行する。図は monitor/ に保存。
"""
import os, sys, json, dataclasses as dc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

for _f in ['Yu Gothic', 'Meiryo', 'MS Gothic']:
    if any(_f.lower() == e.name.lower() for e in matplotlib.font_manager.fontManager.ttflist):
        matplotlib.rcParams['font.family'] = _f; break
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['axes.grid'] = False
matplotlib.rcParams['xtick.direction'] = 'in'; matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True; matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['font.size'] = 11

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
import flowsheet.economics as econ_mod
from flowsheet.economics import calculate_economics
from src.distillation_core import ColumnTunables
from units.reactors.catofin import CatofinDesignVars
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from src.component_data import MW
from src import cost_parameters as cp

P_L_Pa = 1.0e5
APPLY_HI, APPLY_STAGE2, HI_DT_MIN_K = True, False, 10.0
_cfg = load_operating_config()
config = dc.replace(_cfg, spec=dc.replace(_cfg.spec, c3h6_min_wtfrac=0.9945))
BASE_PRICE = cp.LPG_C3H8_JPY_PER_KG

_FEED_STAGE_ABS = {"col2": (2, 9999), "col3": (70, 180)}
def _feed_stage_from_ratio(ratio, n, lo, hi):
    fs = int(round(ratio * n)); hi_eff = min(hi, n - 2); lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))

def build_design(p):
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages']); fs2 = _feed_stage_from_ratio(p['col2_feed_ratio'], n2, *_FEED_STAGE_ABS['col2'])
    n3 = int(p['col3_n_stages']); fs3 = _feed_stage_from_ratio(p['col3_feed_ratio'], n3, *_FEED_STAGE_ABS['col3'])
    p3 = float(p['col3_p_kpa'])
    reactor = CatofinDesignVars(T_in=p['T_in_K'], t_cyc=p['t_cyc_min'], D=p['D_reactor_m'],
                                L_bed=p['L_bed_m'], N_online=int(p['N_online']),
                                d_p=float(p['d_p_mm']) / 1000.0)
    return FlowsheetDesignVars(
        swing=reactor,
        psa=PSADesignVars(D_col=p['D_psa_col_m'], L_bed=p['L_psa_bed_m'], desorption_target=p['desorption_target']),
        mem=MemDesignVars(P_H=p['P_H_Pa'], P_L=P_L_Pa, A_mem=p['A_mem_m2'], P_dist=p3 * 1000.0),
        dist1=ColumnTunables(P_col=float(p['col1_p_kpa']) * 1000.0, N_stages=n1, N_feed=1, reflux_ratio=2.0,
                             solver_method='sm', hysys_spec_value=float(p['col1_comp_frac_2']), hysys_feed_stage=fs1),
        dist2=ColumnTunables(P_col=float(p['col2_p_kpa']) * 1000.0, N_stages=n2, N_feed=1,
                             reflux_ratio=float(p['col2_reflux_ratio']), solver_method='hysys',
                             hysys_spec_value=float(p['col2_reflux_ratio']), hysys_feed_stage=fs2),
        dist3=ColumnTunables(P_col=p3 * 1000.0, N_stages=n3, N_feed=1, reflux_ratio=12.0,
                             solver_method='sm', hysys_spec_value=0.99, hysys_feed_stage=fs3),
    )

BEST = os.path.join(ROOT, 'outputs', 'main_20260604_014318', 'best.json')
best = json.load(open(BEST, encoding='utf-8'))
params = best['params']; design = build_design(params)
F_fresh = float(params['F_C3H8_fresh_kmol_h'])
print(f"最適点 trial #{best['number']}  BO TAC={best['effective_TAC']:.2f}  F_fresh={F_fresh:.1f}  反応器=catofin")

res = evaluate(design, config, verbose=False, apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
               apply_stage2=APPLY_STAGE2, F_C3H8_override=F_fresh)
assert res.economics_hi is not None, f"infeasible: {res.failure_reason}"
econ_hi = res.economics_hi; one_pass = res.solver.one_pass
prod_kmolh = one_pass['r3'].top.F_in.get('B', 0.0)
print(f"  TAC={econ_hi.TAC:.2f}  Rev={econ_hi.total_revenue:.2f}  Profit={econ_hi.profit:.2f}  "
      f"原単価={econ_hi.unit_jpy_per_t/1000:.2f}円/kg  prod={prod_kmolh:.1f}({prod_kmolh/F_fresh*100:.1f}%)")

MW_C3H6 = MW['B']
def econ_preHI_at(price):
    o3, o4 = econ_mod.LPG_C3H8_JPY_PER_KG, econ_mod.LPG_C4H10_JPY_PER_KG
    try:
        econ_mod.LPG_C3H8_JPY_PER_KG = price; econ_mod.LPG_C4H10_JPY_PER_KG = price
        return calculate_economics(one_pass, MW_C3H6)
    finally:
        econ_mod.LPG_C3H8_JPY_PER_KG, econ_mod.LPG_C4H10_JPY_PER_KG = o3, o4

econ_pre_base = econ_preHI_at(BASE_PRICE)
econ_pre_15 = econ_preHI_at(BASE_PRICE * 1.5)
changed = [(k, econ_pre_base.opex[k], econ_pre_15.opex.get(k)) for k in econ_pre_base.opex
           if abs(econ_pre_base.opex[k] - econ_pre_15.opex.get(k, 0.0)) > 1e-9]
raw_total_base = sum(a for _, a, _ in changed)

def metrics_at_price(price):
    dTAC = econ_preHI_at(price).TAC - econ_pre_base.TAC
    tac = econ_hi.TAC + dTAC
    return econ_hi.total_revenue - tac, tac, tac * 1e8 / (econ_hi.annual_kg_C3H6 / 1000.0)

prices = np.linspace(40, 140, 201)
profits = np.array([metrics_at_price(p)[0] for p in prices])
units = np.array([metrics_at_price(p)[2] for p in prices]) / 1000.0
price_be = BASE_PRICE * (1.0 + econ_hi.profit / raw_total_base)
price_be_interp = float(np.interp(0.0, profits[::-1], prices[::-1]))
sens = raw_total_base / BASE_PRICE

print("\n=== 損益分岐サマリ (catofin) ===")
print(f"  基準LPG価格        : {BASE_PRICE:.1f} 円/kg")
print(f"  現行Profit         : {econ_hi.profit:.1f} 億円/年")
print(f"  損益分岐(解析)     : {price_be:.1f} 円/kg")
print(f"  損益分岐(数値)     : {price_be_interp:.1f} 円/kg")
print(f"  下落幅             : {(1-price_be/BASE_PRICE)*100:.1f} %")
print(f"  価格感度           : {sens:.2f} 億円/年 per 円/kg")
print(f"  原料費/TAC         : {raw_total_base/econ_hi.TAC*100:.1f} %")
print(f"  現行製造原単価     : {econ_hi.unit_jpy_per_t/1000:.1f} 円/kg")

# ---- 図1: 価格 vs 利益 ----
fig1, ax = plt.subplots(figsize=(7.0, 4.6))
ax.axhline(0, color='black', lw=0.8)
ax.fill_between(prices, profits, 0, where=(profits > 0), facecolor='none', hatch='////', edgecolor='0.45', linewidth=0.0)
ax.fill_between(prices, profits, 0, where=(profits < 0), facecolor='0.88', edgecolor='none')
ax.plot(prices, profits, lw=2.2, color='black')
ax.axvline(price_be, color='black', ls='--', lw=1.4); ax.axvline(BASE_PRICE, color='black', ls=':', lw=1.6)
ax.plot(BASE_PRICE, econ_hi.profit, 'o', color='black', ms=7, zorder=5)
ax.annotate(f'損益分岐 {price_be:.1f} 円/kg', xy=(price_be, 0), xytext=(price_be + 11, max(profits) * 0.62),
            ha='left', color='black', fontsize=11, arrowprops=dict(arrowstyle='->', color='black'))
ax.annotate(f'現行 {BASE_PRICE:.0f} 円/kg\n{econ_hi.profit:.0f} 億円/年', xy=(BASE_PRICE, econ_hi.profit),
            xytext=(BASE_PRICE + 7, econ_hi.profit - 110), ha='left', color='black', fontsize=10,
            arrowprops=dict(arrowstyle='->', color='black'))
ax.text(46, max(profits) * 0.26, '黒字', color='black', fontsize=13, fontweight='bold')
ax.text(123, min(profits) * 0.72, '赤字', color='black', fontsize=13, fontweight='bold')
ax.set_xlabel('LPG 原料単価 [円/kg]', fontsize=12); ax.set_ylabel('年間利益 Profit [億円/年]', fontsize=12)
ax.set_title('LPG 原料価格と年間利益(最適設計を固定)', fontsize=12); ax.set_xlim(prices[0], prices[-1])
fig1.tight_layout()
fig1.savefig(os.path.join(HERE, 'lpg_breakeven_profit.png'), dpi=200, bbox_inches='tight')
fig1.savefig(os.path.join(HERE, 'lpg_breakeven_profit.pdf'), bbox_inches='tight')

# ---- 図2: 価格 vs 製造原単価 ----
fig2, ax = plt.subplots(figsize=(7.0, 4.6))
ax.plot(prices, units, lw=2.2, color='black', label='C3H6 製造原単価(TAC 基準)')
ax.axhline(cp.C3H6_PRODUCT_JPY_PER_KG, color='black', ls='--', lw=1.4,
           label=f'C3H6 出荷単価 {cp.C3H6_PRODUCT_JPY_PER_KG:.0f} 円/kg')
ax.axvline(price_be, color='0.5', ls='--', lw=1.2); ax.axvline(BASE_PRICE, color='0.5', ls=':', lw=1.4)
u_now = econ_hi.unit_jpy_per_t / 1000
ax.plot(BASE_PRICE, u_now, 'o', color='black', ms=6, zorder=5)
ax.annotate(f'現行 {u_now:.0f} 円/kg', xy=(BASE_PRICE, u_now), xytext=(BASE_PRICE + 5, u_now + 18),
            color='black', fontsize=10, arrowprops=dict(arrowstyle='->', color='black'))
ax.set_xlabel('LPG 原料単価 [円/kg]', fontsize=12); ax.set_ylabel('C3H6 製造原単価 [円/kg]', fontsize=12)
ax.set_title('LPG 原料価格と C3H6 製造原単価(最適設計を固定)', fontsize=12); ax.set_xlim(prices[0], prices[-1])
ax.legend(fontsize=10, loc='upper left', framealpha=1.0, edgecolor='0.3')
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, 'lpg_breakeven_unitcost.png'), dpi=200, bbox_inches='tight')
fig2.savefig(os.path.join(HERE, 'lpg_breakeven_unitcost.pdf'), bbox_inches='tight')
print("\n図を再生成: monitor/lpg_breakeven_profit.{png,pdf} / lpg_breakeven_unitcost.{png,pdf}")
print("=== DONE ===")
