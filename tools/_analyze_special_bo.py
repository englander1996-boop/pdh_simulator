# -*- coding: utf-8 -*-
"""special trials.csv を解析し、BO 精度に効く構造を定量化する一時診断スクリプト。

目的:
  1. Stage1(HI) と Stage2(HEN) の乖離分布 = 目的関数ノイズの定量化
  2. Stage1基準 vs Stage2基準でランキングがどれだけ入れ替わるか (ノイズ→誤選択)
  3. feasible 分布・TAC を動かす変数の確認
使い方: .venv\\Scripts\\python.exe tools\\_analyze_special_bo.py outputs\\special_<ts>_trials.csv
"""
import sys, csv, statistics as st

path = sys.argv[1] if len(sys.argv) > 1 else "outputs/special_20260526_172608_trials.csv"

rows = []
with open(path, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append(r)

def fnum(r, k):
    v = r.get(k, "")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

comp = [r for r in rows if r.get("state") == "COMPLETE"]
feas = [r for r in comp if r.get("attr.is_feasible") == "True"]

print(f"=== {path} ===")
print(f"全trial={len(rows)}  COMPLETE={len(comp)}  feasible={len(feas)} ({100*len(feas)/max(len(comp),1):.0f}%)")

# --- Stage1 vs Stage2 乖離 (feasible のみ economics 有効) ---
gaps = []
for r in feas:
    s1 = fnum(r, "attr.TAC_hi_okuyen")
    s2 = fnum(r, "attr.TAC_stage2_okuyen")
    if s1 is not None and s2 is not None:
        gaps.append((r, s1, s2, s2 - s1))

print("\n=== Stage1(HI) vs Stage2(HEN) 乖離 [億円]  (feasible n={}) ===".format(len(gaps)))
if gaps:
    d = [g[3] for g in gaps]
    print(f"  gap = Stage2 - Stage1:  min={min(d):.1f}  中央={st.median(d):.1f}  "
          f"平均={st.mean(d):.1f}  max={max(d):.1f}  stdev={st.pstdev(d):.1f}")
    # gap が大きい = Stage2 が回収を取りこぼした (greedy flip 候補)
    big = sorted(gaps, key=lambda g: -g[3])[:8]
    print("  gap 上位 (Stage2 が Stage1 比で悪化=HEN取りこぼし):")
    for r, s1, s2, g in big:
        print(f"    #{r['number']:>3}  Stage1={s1:7.1f}  Stage2={s2:7.1f}  gap={g:6.1f}")

# --- ランキング安定性: Stage1基準 top10 vs Stage2基準 top10 ---
print("\n=== ランキング安定性 (feasible) ===")
by_s2 = sorted(gaps, key=lambda g: g[2])  # Stage2(目的関数)
by_s1 = sorted(gaps, key=lambda g: g[1])  # Stage1(安定)
top_s2 = [g[0]["number"] for g in by_s2[:10]]
top_s1 = [g[0]["number"] for g in by_s1[:10]]
print(f"  Stage2基準 top10 (現目的関数): {top_s2}")
print(f"  Stage1基準 top10 (安定指標)  : {top_s1}")
overlap = set(top_s2) & set(top_s1)
print(f"  共通: {len(overlap)}/10  → {sorted(overlap)}")
if by_s2 and by_s1:
    b2 = by_s2[0]; b1 = by_s1[0]
    print(f"  Stage2-best = #{b2[0]['number']} (S2={b2[2]:.1f}, S1={b2[1]:.1f})")
    print(f"  Stage1-best = #{b1[0]['number']} (S1={b1[1]:.1f}, S2={b1[2]:.1f})")
    # Stage1-best を Stage2 で見たら何位? 逆も
    s2_rank_of_s1best = [g[0]["number"] for g in by_s2].index(b1[0]["number"]) + 1
    s1_rank_of_s2best = [g[0]["number"] for g in by_s1].index(b2[0]["number"]) + 1
    print(f"  Stage1-best(#{b1[0]['number']}) は Stage2基準で {s2_rank_of_s1best} 位")
    print(f"  Stage2-best(#{b2[0]['number']}) は Stage1基準で {s1_rank_of_s2best} 位")

# --- failure_unit tally ---
print("\n=== 失敗ユニット集計 (COMPLETE) ===")
from collections import Counter
fu = Counter(r.get("attr.failure_unit", "") for r in comp)
for k, v in fu.most_common():
    print(f"  {k or '(空)':<12} {v}")

# --- production 方向 (feasible 化を阻む要因) ---
print("\n=== production_direction 集計 (COMPLETE) ===")
pd = Counter(r.get("attr.production_direction", "") for r in comp)
for k, v in pd.most_common():
    print(f"  {k or '(空)':<12} {v}")
