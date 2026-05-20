"""Analyze F_fresh vs yield vs TAC from BO trials."""
import json
import pandas as pd
import numpy as np

CSV = "z:/pdh_simulator/outputs/main_20260520_003551/trials.csv"
TOPK = "z:/pdh_simulator/outputs/main_20260520_003551/topk.txt"
BEST = "z:/pdh_simulator/outputs/main_20260520_003551/best.json"

df = pd.read_csv(CSV)
print(f"Total trials: {len(df)}")
print(f"State value_counts: {df['state'].value_counts().to_dict()}")

# Definitions
df["F_fresh"] = df["param_F_C3H8_fresh_kmol_h"]
df["production"] = df["attr_production_kmol_h"]
df["TAC"] = df["value"]
df["TAC_re"] = df["attr_TAC_okuyen"]
df["is_feas"] = df["attr_is_feasible"]
df["yield"] = df["production"] / df["F_fresh"]

# Feasible: value(TAC_bo) < 10000 sentinel and not NaN
fea = df[(df["TAC"] < 10000) & df["TAC"].notna()].copy()
print(f"\nFeasible (TAC<10000) trials: {len(fea)} / {len(df)}")
print(f"is_feasible==True trials: {df['is_feas'].sum()}")
print(f"production not-null trials: {df['production'].notna().sum()}")

# 2. Yield distribution among feasible
y = fea["yield"]
print("\n=== [Feasible] yield (production / F_fresh) distribution ===")
print(f"  n      : {len(y)}")
print(f"  min    : {y.min():.4f}")
print(f"  25%    : {y.quantile(0.25):.4f}")
print(f"  median : {y.median():.4f}")
print(f"  75%    : {y.quantile(0.75):.4f}")
print(f"  max    : {y.max():.4f}")
print(f"  mean   : {y.mean():.4f}")

# Yield vs TAC: top50 (lowest TAC) vs bot50 (highest TAC)
fea_sorted = fea.sort_values("TAC")
top50 = fea_sorted.head(50)
bot50 = fea_sorted.tail(50)
print("\n=== Yield: top-50 (lowest TAC) vs bot-50 (highest TAC) among feasible ===")
print(f"  top50  yield median = {top50['yield'].median():.4f}  (TAC median={top50['TAC'].median():.2f})")
print(f"  bot50  yield median = {bot50['yield'].median():.4f}  (TAC median={bot50['TAC'].median():.2f})")
# Spearman
rho = fea[["yield", "TAC"]].corr(method="spearman").iloc[0, 1]
pearson = fea[["yield", "TAC"]].corr(method="pearson").iloc[0, 1]
print(f"  Spearman(yield,TAC) = {rho:.3f}  Pearson = {pearson:.3f}  (negative = higher yield -> lower TAC)")

# 3. Bucket F_fresh
bins = [1200, 1300, 1400, 1500, 1600, 1700]
labels = ["1200-1300", "1300-1400", "1400-1500", "1500-1600", "1600-1700"]
df["bucket"] = pd.cut(df["F_fresh"], bins=bins, labels=labels, include_lowest=True, right=False)
fea["bucket"] = pd.cut(fea["F_fresh"], bins=bins, labels=labels, include_lowest=True, right=False)

print("\n=== F_fresh bucket summary ===")
print(f"{'bucket':<12} {'n_all':>6} {'n_fea':>6} {'feas%':>7} {'TAC_med':>10} {'TAC_min':>10} {'yld_med':>9} {'yld_min':>9} {'yld_max':>9}")
for lab in labels:
    all_b = df[df["bucket"] == lab]
    fea_b = fea[fea["bucket"] == lab]
    n_all = len(all_b)
    n_fea = len(fea_b)
    feas_pct = (n_fea / n_all * 100) if n_all else 0
    if n_fea > 0:
        tac_med = fea_b["TAC"].median()
        tac_min = fea_b["TAC"].min()
        y_med = fea_b["yield"].median()
        y_min = fea_b["yield"].min()
        y_max = fea_b["yield"].max()
        print(f"{lab:<12} {n_all:>6} {n_fea:>6} {feas_pct:>6.1f}% {tac_med:>10.2f} {tac_min:>10.2f} {y_med:>9.4f} {y_min:>9.4f} {y_max:>9.4f}")
    else:
        print(f"{lab:<12} {n_all:>6} {n_fea:>6} {feas_pct:>6.1f}%        n/a        n/a       n/a       n/a       n/a")

# 4. topk: best.json contains topk array
with open(BEST, "r", encoding="utf-8") as f:
    best = json.load(f)
print("\n=== best.json top-level keys ===")
print(list(best.keys()) if isinstance(best, dict) else type(best))

# Try to read topk.txt
with open(TOPK, "r", encoding="utf-8") as f:
    topk_txt = f.read()
print("\n=== topk.txt head ===")
print(topk_txt[:2000])

# Also derive top 10 from CSV (lowest TAC among feasible)
top10 = fea.sort_values("TAC").head(10)[
    ["trial_number", "F_fresh", "production", "yield", "TAC", "TAC_re", "is_feas",
     "attr_c3h6_purity_wtfrac", "attr_h2_purity_molfrac"]
]
print("\n=== Top-10 feasible trials (by TAC_bo) from CSV ===")
print(top10.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

# 5. Where does TAC-min trial come from? Check yield < 0.80 region
print("\n=== TAC-min trials low-yield check ===")
fea["yield_band"] = pd.cut(fea["yield"], bins=[0, 0.70, 0.75, 0.80, 0.85, 0.90, 1.0],
                            labels=["<.70", ".70-.75", ".75-.80", ".80-.85", ".85-.90", ".90+"])
yb = fea.groupby("yield_band", observed=True).agg(
    n=("TAC", "size"),
    TAC_med=("TAC", "median"),
    TAC_min=("TAC", "min"),
    F_med=("F_fresh", "median"),
).reset_index()
print(yb.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

# Top-10 TAC trials yield band membership
top10_yld = fea.sort_values("TAC").head(10)
print("\nYield of top-10 feasible (lowest TAC):")
for _, r in top10_yld.iterrows():
    print(f"  trial={int(r['trial_number']):4d}  F_fresh={r['F_fresh']:.1f}  yield={r['yield']:.4f}  TAC={r['TAC']:.2f}")

# Distribution: production vs target
print("\n=== production vs target=1188.21 ===")
print(f"feasible production min  = {fea['production'].min():.2f}")
print(f"feasible production med  = {fea['production'].median():.2f}")
print(f"feasible production max  = {fea['production'].max():.2f}")
short = fea[fea["production"] < 1188.21]
print(f"feasible w/ production < 1188.21: {len(short)} / {len(fea)} (should be 0 for true feas)")

# F_fresh value range observed for feasible
print("\n=== F_fresh range among feasible ===")
print(f"  min={fea['F_fresh'].min():.2f}  med={fea['F_fresh'].median():.2f}  max={fea['F_fresh'].max():.2f}")
print(f"  count with F_fresh<1300: {(fea['F_fresh']<1300).sum()}")
print(f"  count with F_fresh<1400: {(fea['F_fresh']<1400).sum()}")
print(f"  count with F_fresh<1500: {(fea['F_fresh']<1500).sum()}")
print(f"  count with F_fresh<1600: {(fea['F_fresh']<1600).sum()}")
