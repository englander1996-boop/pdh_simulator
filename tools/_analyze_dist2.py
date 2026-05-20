"""
Analysis: naze BO ha Dist2 atsui sekkei wo erabanaika
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import pandas as pd
import numpy as np
from pathlib import Path

csv_path = Path(r"z:\pdh_simulator\outputs\main_20260520_124508\trials.csv")
df = pd.read_csv(csv_path)

print("=" * 80)
print(f"Total trials: {len(df)}  COMPLETE: {(df['state']=='COMPLETE').sum()}")
print(f"  feasible (BO): {(df['attr_is_feasible']==True).sum()}")
print()

# value 列 (BO TAC): solve 失敗時 10000
df['TAC_bo'] = pd.to_numeric(df['value'], errors='coerce')
df['N_dist2'] = pd.to_numeric(df['param_N_dist2'], errors='coerce')
df['R_dist2'] = pd.to_numeric(df['param_reflux_dist2'], errors='coerce')
df['rec_LK'] = pd.to_numeric(df['param_rec_LK_top_dist2'], errors='coerce')
df['rec_HK'] = pd.to_numeric(df['param_rec_HK_bot_dist2'], errors='coerce')
df['feas']   = df['attr_is_feasible'] == True

# ---------------------------------------------------------------- §1 分布
print("=" * 80)
print("§1. N_dist2 の探索分布 (5 bins)")
print("=" * 80)
bins = [19.5, 24.5, 28.5, 32.5, 36.5, 40.5]
labels = ['20-24', '25-28', '29-32', '33-36', '37-40']
df['N_bin'] = pd.cut(df['N_dist2'], bins=bins, labels=labels)

# 全体
all_hist = df['N_bin'].value_counts().reindex(labels)
# startup (0..49) / main 前半 (50..199) / main 後半 (200..299)
startup = df[df['trial_number'] < 50]['N_bin'].value_counts().reindex(labels)
mid     = df[(df['trial_number'] >= 50) & (df['trial_number'] < 200)]['N_bin'].value_counts().reindex(labels)
late    = df[df['trial_number'] >= 200]['N_bin'].value_counts().reindex(labels)

hist = pd.DataFrame({
    'all(300)':   all_hist,
    'startup(50)': startup,
    'mid(150)':    mid,
    'late(100)':   late,
})
hist['late_ratio_%'] = (late / 100 * 100).round(1)
print(hist.fillna(0).astype({'all(300)':int,'startup(50)':int,'mid(150)':int,'late(100)':int}))
print()
print(f"  N_dist2 全期間 mean = {df['N_dist2'].mean():.2f}, median = {df['N_dist2'].median():.1f}")
print(f"  late 100 trial   mean = {df[df['trial_number']>=200]['N_dist2'].mean():.2f}, median = {df[df['trial_number']>=200]['N_dist2'].median():.1f}")

# ---------------------------------------------------------------- §2 経済性
print()
print("=" * 80)
print("§2. Dist2 厚い (N≥32) vs 薄い (N≤24) の TAC_bo 比較 (feasible のみ)")
print("=" * 80)

feas = df[df['feas']].copy()
print(f"  feasible total: {len(feas)}")

thick = feas[feas['N_dist2'] >= 32]
thin  = feas[feas['N_dist2'] <= 24]
mid_g = feas[(feas['N_dist2'] > 24) & (feas['N_dist2'] < 32)]

def stats(g, name):
    if len(g) == 0:
        return {'group': name, 'n': 0}
    return {
        'group': name,
        'n': len(g),
        'TAC_med': round(g['TAC_bo'].median(), 2),
        'TAC_min': round(g['TAC_bo'].min(), 2),
        'TAC_max': round(g['TAC_bo'].max(), 2),
        'TAC_p25': round(g['TAC_bo'].quantile(0.25), 2),
        'rec_LK_med': round(g['rec_LK'].median(), 4),
        'rec_HK_med': round(g['rec_HK'].median(), 5),
        'R_med': round(g['R_dist2'].median(), 2),
        'N_med': round(g['N_dist2'].median(), 1),
    }

comp = pd.DataFrame([
    stats(thin,  'thin (N<=24)'),
    stats(mid_g, 'mid  (25-31)'),
    stats(thick, 'thick(N>=32)'),
])
print(comp.to_string(index=False))

# ---------------------------------------------------------------- §3 価値逆算
print()
print("=" * 80)
print("§3. C3H6 漏れ 1pp の経済的価値 (近い trial ペアでの差)")
print("=" * 80)

# まず top-feasible (TAC<210) の中で N=22 (薄) と N≥32 (厚) を比較
top_thin  = feas[(feas['N_dist2']==22) & (feas['TAC_bo']<220)]
top_thick = feas[(feas['N_dist2']>=32) & (feas['TAC_bo']<300)]
print(f"\n  N=22, TAC<220 の trial:  {len(top_thin)} 件")
print(f"  N>=32, TAC<300 の trial: {len(top_thick)} 件")

# TAC vs N_dist2 (feasible) を 5 ビン化
print("\n  feasible TAC_bo bin 別 N_dist2 平均:")
feas_lowtac = feas[feas['TAC_bo'] < 500].copy()
feas_lowtac['TAC_bin'] = pd.cut(feas_lowtac['TAC_bo'], bins=[0,200,250,300,400,500], labels=['<200','200-250','250-300','300-400','400-500'])
agg = feas_lowtac.groupby('TAC_bin', observed=True).agg(
    n=('TAC_bo','size'),
    N_med=('N_dist2','median'),
    N_mean=('N_dist2','mean'),
    R_med=('R_dist2','median'),
    rec_LK_med=('rec_LK','median'),
    rec_HK_med=('rec_HK','median'),
).round(3)
print(agg)

# 同じような他パラの trial で N_dist2 だけ違う peer を探す
# best trial #246 と近い: N_dist1=26, P_dist2~540kPa, F~1420 を絞り込み
ref = df[df['trial_number'] == 246].iloc[0]
near = feas[
    (feas['param_N_dist1'] == ref['param_N_dist1']) &
    (feas['param_F_C3H8_fresh_kmol_h'].between(ref['param_F_C3H8_fresh_kmol_h']-30, ref['param_F_C3H8_fresh_kmol_h']+30)) &
    (feas['param_P_dist2_Pa'].between(ref['param_P_dist2_Pa']-60000, ref['param_P_dist2_Pa']+120000))
].copy()
print(f"\n  best trial 246 と近傍 (N_dist1, F, P_dist2 近似) feasible: {len(near)} 件")
near_show = near[['trial_number','N_dist2','R_dist2','rec_LK','rec_HK','TAC_bo','attr_c3h6_purity_wtfrac','attr_production_kmol_h']].sort_values('N_dist2')
print(near_show.to_string(index=False))

# ---------------------------------------------------------------- §4 構造的相関
print()
print("=" * 80)
print("§4. SEARCH_SPACE 構造: 厚い N と タイトな rec の組合せ実現率")
print("=" * 80)

# all complete (non-failed)
comp_df = df[df['state']=='COMPLETE'].copy()
combos = pd.crosstab(
    pd.cut(comp_df['N_dist2'], bins=bins, labels=labels),
    pd.cut(comp_df['rec_LK'], bins=[0.949,0.97,0.985,0.995,1.001], labels=['0.95-0.97','0.97-0.985','0.985-0.995','0.995+']),
)
print("  全 COMPLETE trial: N_dist2 × rec_LK_top_dist2 件数")
print(combos)

combos_f = pd.crosstab(
    pd.cut(feas['N_dist2'], bins=bins, labels=labels),
    pd.cut(feas['rec_LK'], bins=[0.949,0.97,0.985,0.995,1.001], labels=['0.95-0.97','0.97-0.985','0.985-0.995','0.995+']),
)
print("\n  feasible trial: N_dist2 × rec_LK_top_dist2 件数")
print(combos_f)

# N>=33 で rec_LK >= 0.985 (厚×タイト) の経済性
strict = feas[(feas['N_dist2'] >= 33) & (feas['rec_LK'] >= 0.985)]
print(f"\n  N>=33 & rec_LK>=0.985 の feasible: {len(strict)} 件")
if len(strict) > 0:
    print(strict[['trial_number','N_dist2','R_dist2','rec_LK','rec_HK','TAC_bo','attr_c3h6_purity_wtfrac']].head(10).to_string(index=False))

# ---------------------------------------------------------------- §5 revenue
print()
print("=" * 80)
print("§5. PSA offgas 燃料 revenue の支配率 (rec_LK 緩 vs タイト)")
print("=" * 80)

rev_loose  = feas[feas['rec_LK'] < 0.97]
rev_tight  = feas[feas['rec_LK'] >= 0.985]
for g, name in [(rev_loose, 'rec_LK<0.97 (緩・C3H6リサイクル少)'),
                (rev_tight, 'rec_LK>=0.985 (タイト・厳格分離)')]:
    if len(g)>0:
        print(f"\n  {name} (n={len(g)})")
        print(f"    TAC_bo med    = {g['TAC_bo'].median():.2f}")
        print(f"    revenue med   = {pd.to_numeric(g['attr_revenue_okuyen'],errors='coerce').median():.2f}")
        print(f"    revenue_HI med= {pd.to_numeric(g['attr_revenue_hi_okuyen'],errors='coerce').median():.2f}")
        print(f"    profit_raw med= {pd.to_numeric(g['attr_profit_raw_okuyen'],errors='coerce').median():.2f}")
        print(f"    production med= {pd.to_numeric(g['attr_production_kmol_h'],errors='coerce').median():.2f}")
