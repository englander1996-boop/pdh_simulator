"""
Part 2: thick (#115) vs thin (#246) deep dive + searching dist2 sweep
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import pandas as pd
import numpy as np
from pathlib import Path

csv_path = Path(r"z:\pdh_simulator\outputs\main_20260520_124508\trials.csv")
df = pd.read_csv(csv_path)
df['N_dist2'] = pd.to_numeric(df['param_N_dist2'], errors='coerce')
df['R_dist2'] = pd.to_numeric(df['param_reflux_dist2'], errors='coerce')
df['rec_LK']  = pd.to_numeric(df['param_rec_LK_top_dist2'], errors='coerce')
df['rec_HK']  = pd.to_numeric(df['param_rec_HK_bot_dist2'], errors='coerce')
df['TAC_bo']  = pd.to_numeric(df['value'], errors='coerce')
df['feas']    = df['attr_is_feasible'] == True

# --- Q1: c3h6_purity_wtfrac は実は (BO 段階で) Dist2 厚薄でほぼ同じ?
print("§ feasible 群: c3h6_purity_wtfrac の分布")
fa = df[df['feas']].copy()
print(fa['attr_c3h6_purity_wtfrac'].describe())
print()

# --- Q2: rev_okuyen と TAC_okuyen (Hasebe + Hi 後) の関係
print("§ feasible 群: revenue, production, profit_raw 分布")
for col in ['attr_revenue_okuyen','attr_revenue_hi_okuyen','attr_TAC_okuyen','attr_TAC_hi_okuyen','attr_profit_raw_okuyen','attr_profit_hi_okuyen','attr_production_kmol_h']:
    fa[col] = pd.to_numeric(fa[col], errors='coerce')
print(fa[['attr_revenue_okuyen','attr_TAC_okuyen','attr_profit_hi_okuyen','attr_production_kmol_h']].describe().round(2))

# --- Q3: 厚い trial #115 と 薄い trial #246 を直接比較
print()
print("§ trial #246 (薄 best) vs #115 (厚 唯一 feasible) 詳細")
for tn in [246, 115]:
    r = df[df['trial_number']==tn].iloc[0]
    print(f"\n--- trial #{tn} ---")
    for k in ['N_dist2','R_dist2','rec_LK','rec_HK','TAC_bo',
             'attr_TAC_okuyen','attr_TAC_hi_okuyen',
             'attr_revenue_okuyen','attr_revenue_hi_okuyen',
             'attr_profit_raw_okuyen','attr_profit_hi_okuyen',
             'attr_c3h6_purity_wtfrac','attr_production_kmol_h',
             'attr_proxy_penalty_total_okuyen',
             'attr_proxy_penalty_r1_okuyen','attr_proxy_penalty_r2_okuyen','attr_proxy_penalty_r3_okuyen']:
        v = r.get(k, np.nan)
        try:
            v = float(v)
            print(f"  {k:42s} = {v:12.4f}")
        except Exception:
            print(f"  {k:42s} = {v}")

# --- Q4: BO objective formula (TAC_bo = ?)
print()
print("§ TAC_bo = TAC_hi - revenue_hi (検算 feasible)")
fa['TAC_check'] = fa['attr_TAC_hi_okuyen'] - fa['attr_revenue_hi_okuyen']
print(fa[['trial_number','TAC_bo','attr_TAC_hi_okuyen','attr_revenue_hi_okuyen','TAC_check']].head(5).to_string(index=False))

# --- Q5: 厚 vs 薄 で TAC_hi (CAPEX+OPEX) はどう違う? revenue はどう違う?
print()
print("§ feasible 群を N_dist2 で分けて TAC_hi / revenue_hi / profit_hi の差")
for label, mask in [('N<=24',fa['N_dist2']<=24),('N=25-31',(fa['N_dist2']>24)&(fa['N_dist2']<32)),('N>=32',fa['N_dist2']>=32)]:
    g = fa[mask]
    if len(g)==0: continue
    print(f"\n  {label} (n={len(g)})")
    print(f"    TAC_hi    med = {g['attr_TAC_hi_okuyen'].median():8.2f}")
    print(f"    TAC       med = {g['attr_TAC_okuyen'].median():8.2f}")
    print(f"    revenue_hi med = {g['attr_revenue_hi_okuyen'].median():8.2f}")
    print(f"    profit_hi (= -TAC_bo) med = {g['attr_profit_hi_okuyen'].median():8.2f}")

# --- Q6: BO が「薄い」を選ぶ trial を「厚い同等条件」と差分:
# trial #115 と他厚 trial があれば
thick_all = df[(df['N_dist2']>=33) & (df['state']=='COMPLETE')]
print()
print(f"§ N>=33 の COMPLETE trial 数 = {len(thick_all)}")
print(f"   うち feasible = {(thick_all['attr_is_feasible']==True).sum()}")
print(f"   TAC_bo 統計 (全 COMPLETE):")
print(thick_all['TAC_bo'].describe().round(2))

# --- Q7: 厚 trial で feasible にならなかった主因
print()
print("§ N>=33 で infeasible になった trial の failure_reason (頭出し)")
infeas_thick = thick_all[thick_all['attr_is_feasible']!=True]
reasons = infeas_thick['attr_failure_reason'].fillna('').str[:80].value_counts().head(10)
print(reasons)
