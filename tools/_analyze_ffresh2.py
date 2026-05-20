"""Stricter analysis: distinguish 'BO-feasible (TAC<10000)' from 'true feasible (production>=target & is_feas)'."""
import pandas as pd

CSV = "z:/pdh_simulator/outputs/main_20260520_003551/trials.csv"
TARGET = 1188.21

df = pd.read_csv(CSV)
df["F_fresh"] = df["param_F_C3H8_fresh_kmol_h"]
df["production"] = df["attr_production_kmol_h"]
df["TAC"] = df["value"]
df["TAC_re"] = df["attr_TAC_okuyen"]
df["is_feas"] = df["attr_is_feasible"]
df["yield"] = df["production"] / df["F_fresh"]

fea_bo = df[(df["TAC"] < 10000) & df["TAC"].notna()].copy()
fea_true = fea_bo[fea_bo["is_feas"] == True].copy()
prod_ok = fea_bo[fea_bo["production"] >= TARGET].copy()

print(f"BO-feasible (TAC<10000)            : n={len(fea_bo)}")
print(f"True-feasible (attr_is_feasible==T) : n={len(fea_true)}")
print(f"production >= target (1188.21)     : n={len(prod_ok)}")

print("\n=== True-feasible (n=3) detail ===")
cols = ["trial_number", "F_fresh", "production", "yield", "TAC", "TAC_re"]
print(fea_true[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

print("\n=== production>=target (overshoot OK) detail ===")
print(prod_ok.sort_values("TAC")[cols].head(15).to_string(index=False, float_format=lambda x: f"{x:.4f}"))

print("\n=== F_fresh of production>=target trials ===")
print(f"  min={prod_ok['F_fresh'].min():.2f}  med={prod_ok['F_fresh'].median():.2f}  max={prod_ok['F_fresh'].max():.2f}")
print(f"  yield: min={prod_ok['yield'].min():.4f}  med={prod_ok['yield'].median():.4f}  max={prod_ok['yield'].max():.4f}")

# Need F_fresh >= ? to have any chance of hitting target
# yield max ~ 0.89; so F_fresh_min ~ target / 0.89 ~ 1335 to be feasible at all.
print(f"\nMin F_fresh needed assuming yield<=0.89: {TARGET/0.89:.1f}")
print(f"Min F_fresh needed assuming yield<=0.88: {TARGET/0.88:.1f}")
print(f"Min F_fresh needed assuming yield<=0.85: {TARGET/0.85:.1f}")

# Bucket by F_fresh: fraction that hit target
bins = [1200, 1300, 1400, 1500, 1600, 1700]
labels = ["1200-1300", "1300-1400", "1400-1500", "1500-1600", "1600-1700"]
df["bucket"] = pd.cut(df["F_fresh"], bins=bins, labels=labels, include_lowest=True, right=False)

print("\n=== Bucket: hit-target rate (production>=1188.21) and is_feas rate ===")
print(f"{'bucket':<12} {'n_all':>6} {'prod>=T':>9} {'%hit':>7} {'is_feas':>8} {'%feas':>7} {'min_TAC_prod>=T':>17}")
for lab in labels:
    sub = df[df["bucket"] == lab]
    n = len(sub)
    hit = (sub["production"] >= TARGET).sum()
    isf = (sub["is_feas"] == True).sum()
    sub_hit = sub[sub["production"] >= TARGET]
    min_tac = sub_hit["TAC"].min() if len(sub_hit) else float("nan")
    print(f"{lab:<12} {n:>6} {hit:>9} {hit/n*100:>6.1f}% {isf:>8} {isf/n*100:>6.1f}% {min_tac:>17.2f}")

# Among feasible-by-production trials, distribution of F_fresh
print("\n=== F_fresh histogram of trials with production>=target ===")
if len(prod_ok):
    pd_hist = prod_ok["F_fresh"].describe()
    print(pd_hist)

# Best TAC by region of F_fresh upper cap
print("\n=== Best TAC achievable if we cap F_fresh upper at X ===")
for cap in [1300, 1400, 1500, 1600, 1700]:
    sub = prod_ok[prod_ok["F_fresh"] <= cap]
    if len(sub):
        print(f"  cap={cap}: n={len(sub)}  best TAC={sub['TAC'].min():.2f}  (trial #{int(sub.loc[sub['TAC'].idxmin(),'trial_number'])}, F_fresh={sub.loc[sub['TAC'].idxmin(),'F_fresh']:.1f})")
    else:
        print(f"  cap={cap}: NO feasible trial below cap")

# Yield in true feasibility group
print("\n=== Yield distribution among production>=target ===")
y = prod_ok["yield"]
if len(y):
    print(f"  n={len(y)} min={y.min():.4f} 25%={y.quantile(0.25):.4f} med={y.median():.4f} 75%={y.quantile(0.75):.4f} max={y.max():.4f}")
