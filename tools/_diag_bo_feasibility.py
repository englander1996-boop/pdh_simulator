"""Diagnose BO feasibility drop: 124508 (fug, 14%) vs 182716 (rigorous, 1%)."""
import pandas as pd
import numpy as np

NEW = r"outputs/main_20260520_182716/trials.csv"
OLD = r"outputs/main_20260520_124508/trials.csv"

dfn = pd.read_csv(NEW)
dfo = pd.read_csv(OLD)

print("=== Shape & state ===")
print("new:", dfn.shape, dfn['state'].value_counts().to_dict())
print("old:", dfo.shape, dfo['state'].value_counts().to_dict())

# is_feasible
fn = dfn[dfn['attr_is_feasible'] == True]
io = dfo[dfo['attr_is_feasible'] == True]
print(f"\nfeasible: new={len(fn)} / {len(dfn)}  old={len(io)} / {len(dfo)}")

infn = dfn[dfn['attr_is_feasible'] != True]
print(f"infeasible: new={len(infn)}")

# === 1. failure_reason distribution (new run) ===
print("\n=== 1. failure_reason distribution (new run infeasible 297) ===")
reasons = infn['attr_failure_reason'].fillna('(none)')
def bucket(r):
    s = str(r)
    if '生産量' in s or 'production' in s.lower() or '未達' in s: return '生産量未達'
    if 'proxy' in s.lower():
        if 'r1' in s.lower() or 'dist1' in s.lower(): return 'proxy(Dist1)'
        if 'r2' in s.lower() or 'dist2' in s.lower(): return 'proxy(Dist2)'
        if 'r3' in s.lower() or 'dist3' in s.lower(): return 'proxy(Dist3)'
        return 'proxy(other)'
    if 'recovery' in s.lower() or 'strict' in s.lower(): return 'strict recovery'
    if 'psa' in s.lower() or 'mem' in s.lower() or 'capex' in s.lower(): return 'PSA/Mem CAPEX'
    if 'reactor' in s.lower() or 'sv' in s.lower(): return 'Reactor SV'
    if 'solver' in s.lower() or 'fail' in s.lower() or 'converge' in s.lower(): return 'solver fail'
    if s == '(none)' or s == 'nan': return '(none/null)'
    return f'OTHER: {s[:60]}'

infn = infn.copy()
infn['_bucket'] = reasons.map(bucket)
vc = infn['_bucket'].value_counts()
for k, v in vc.items():
    print(f"  {k:40s} {v:4d}  ({100*v/len(infn):.1f}%)")

# raw top-10 unique strings (truncated)
print("\n  -- raw top-15 failure_reason strings --")
raw_vc = reasons.value_counts().head(15)
for s, c in raw_vc.items():
    print(f"   [{c:4d}] {str(s)[:120]}")

# === 2. feasible 3 vs infeasible 297 (new) ===
params = [c for c in dfn.columns if c.startswith('param_')]
print("\n=== 2. new run: feasible-3 vs infeasible-297 design vars ===")
print(f"{'var':30s} {'feas_med':>12s} {'inf_med':>12s} {'feas_min':>12s} {'feas_max':>12s}")
for p in params:
    fm = fn[p].median() if len(fn) else float('nan')
    im = infn[p].median()
    fmin = fn[p].min() if len(fn) else float('nan')
    fmax = fn[p].max() if len(fn) else float('nan')
    print(f"{p:30s} {fm:>12.4g} {im:>12.4g} {fmin:>12.4g} {fmax:>12.4g}")

# show the 3 feasible trial dicts
print("\n  -- 3 feasible trials (new run) full params --")
print(fn[['trial_number','attr_TAC_okuyen','attr_production_kmol_h','attr_failure_reason'] + params].to_string(index=False))

# === 3. old run: feasible 41 distribution ===
print(f"\n=== 3. OLD run: feasible {len(io)} trials  param distribution (q25 / med / q75) ===")
print(f"{'var':30s} {'q25':>12s} {'med':>12s} {'q75':>12s} {'min':>12s} {'max':>12s}")
old_stats = {}
for p in params:
    if p not in io.columns: continue
    q25 = io[p].quantile(0.25)
    med = io[p].median()
    q75 = io[p].quantile(0.75)
    mn  = io[p].min()
    mx  = io[p].max()
    old_stats[p] = (q25, med, q75, mn, mx)
    print(f"{p:30s} {q25:>12.4g} {med:>12.4g} {q75:>12.4g} {mn:>12.4g} {mx:>12.4g}")

# === 4. Suggested narrowing ===
SEARCH = {
    'param_T_in_K':              (900.0,  970.0),
    'param_z_cat_m':              (15.0,   40.0),
    'param_t_cyc_min':            (12.0,   25.0),
    'param_D_reactor_m':          (7.0,    10.0),
    'param_D_psa_col_m':          (2.5,    5.0),
    'param_L_psa_bed_m':          (15.0,   30.0),
    'param_desorption_target':    (0.15,   0.40),
    'param_P_H_Pa':               (7.5e5,  9.5e5),
    'param_A_mem_m2':             (5.0e4,  3.0e5),  # 2026-05-22 E-plan で下限引上げ
    'param_P_dist1_Pa':           (12e5,   25e5),
    'param_N_dist1':              (16,     30),
    'param_reflux_dist1':         (1.3,    3.0),
    'param_P_dist2_Pa':           (5e5,    7e5),
    'param_N_dist2':              (20,     50),  # 2026-05-22 改良 3 で上限拡張
    'param_reflux_dist2':         (5.0,    10.0),
    'param_P_dist3_Pa':           (16e5,   25e5),  # 2026-05-22 E-plan で下限引上げ
    'param_N_dist3':              (80,     200),
    'param_reflux_dist3':         (11.0,   20.0),
    'param_F_C3H8_fresh_kmol_h':  (1200.0, 1700.0),
    'param_rec_LK_top_dist2':     (0.95,   0.999),
    'param_rec_HK_bot_dist2':     (0.9995, 0.9999),
}
print("\n=== 4. proposed SEARCH_SPACE narrowing (based on OLD feasible 25-75% +/- margin) ===")
print(f"{'var':30s} {'cur_low':>12s} {'cur_high':>12s} {'old_q25':>10s} {'old_q75':>10s} {'sug_low':>12s} {'sug_high':>12s}")
for p, (lo, hi) in SEARCH.items():
    if p not in old_stats: continue
    q25, med, q75, mn, mx = old_stats[p]
    # margin: 10% of (q75 - q25), or 5% of (hi - lo) -- whichever larger
    span = q75 - q25
    marg = max(0.10 * span, 0.05 * (hi - lo))
    sug_lo = max(lo, q25 - marg)
    sug_hi = min(hi, q75 + marg)
    # ints
    if 'N_dist' in p:
        sug_lo = int(np.floor(sug_lo))
        sug_hi = int(np.ceil(sug_hi))
    print(f"{p:30s} {lo:>12.4g} {hi:>12.4g} {q25:>10.4g} {q75:>10.4g} {sug_lo:>12.4g} {sug_hi:>12.4g}")

# === 5. cross-check: are old-feasible param ranges populated in new run? ===
# i.e. did the new run actually sample near old-feasible region?
print("\n=== 5. new run sampling coverage of OLD feasible q25-q75 ===")
print(f"{'var':30s} {'old_q25':>10s} {'old_q75':>10s} {'new_in_range':>14s} {'new_total':>10s}")
for p in params:
    if p not in old_stats: continue
    q25, med, q75, mn, mx = old_stats[p]
    inrng = ((dfn[p] >= q25) & (dfn[p] <= q75)).sum()
    print(f"{p:30s} {q25:>10.4g} {q75:>10.4g} {inrng:>14d} {len(dfn):>10d}")
