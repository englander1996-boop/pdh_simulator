r"""
comparing/case_p02_pinch/main.py — P02「熱統合の後置 (post-hoc pinch)」の忠実再現。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p02_pinch\main.py`

再現する問題のあるやり方 (最適化手法分析 C04。多数のレポートがこの手順):
  1. まず熱統合を考えずにプロセスを設計・最適化する (apply_hi=False)。
  2. 設計を確定した後で、後付けでピンチ解析/熱交換器網を当てる (apply_hi=True)。
     その際 ΔTmin は慣用値 10K に固定する。

なぜ問題か (定量化のねらい):
  - P02: 熱統合を最適化ループの外 (後段) に置くと、プロセス温度・流量が熱回収を考えずに
    決まってしまい、回収しきれない熱が残る。熱統合をループ内で同時に扱う BO との TAC 差が
    「後置熱統合の損失」。さらに後付け ΔTmin の選び方 (10K 固定 vs 掃引最良) でも差が出る。

実装 (2 フェーズ):
  Phase1: 熱統合なし (apply_hi=False, apply_stage2=False) でコストドライバを座標降下最適化 → 設計 D*。
  Phase2: D* を固定し ΔTmin を掃引 (apply_hi=True, apply_stage2=True)。慣用 10K と掃引最良を比較。

出力 (comparing/results/pinch_posthoc_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式、Phase2=後置HI 結果)
  dtmin_sweep.csv  : D* を ΔTmin 各点で評価した (hi_dT_min_K, effective_TAC, ...) — 10K と最良の差
  phase1_no_hi.csv : Phase1 (熱統合なし) の設計最適化の掃引軌跡

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
BO との比較はユーザ側: best.json と outputs/special_*/best.json を突合し ΔTAC を算出。
"""

import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

from comparing.shared import space, simulator, reporting, harness


# ===========================================================================
# 設定 (編集可)
# ===========================================================================
# Phase1 (熱統合なし) で座標降下するコストドライバ。
PHASE1_VARS = [
    'col2_reflux_ratio',
    'col3_n_stages',
    'F_C3H8_fresh_kmol_h',
    'T_in_K',
]
K_POINTS = 3                               # Phase1 の 1 変数あたり掃引点数
DTMIN_SWEEP_K = [5.0, 8.0, 10.0, 15.0, 20.0, 30.0]  # Phase2 で当てる ΔTmin [K] (10K=慣用)
POSTHOC_DTMIN_K = 10.0                     # 後置慣用値 (= 学生がそのまま採る値)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P02 熱統合の後置 の忠実再現。Phase1 で熱統合なし (apply_hi=False) にプロセスを座標降下最適化して
設計 D* を確定し、Phase2 で D* を固定したまま後付けで ΔTmin を当てる (apply_hi=True)。
熱統合をループ外に置くことで取り逃す熱回収を、BO (ループ内同時熱統合) best との ΔTAC で定量化する。
dtmin_sweep.csv に慣用 10K と掃引最良の差も記録。
"""


def run(*, phase1_vars=None, k=K_POINTS, dtmin_sweep=None, start=None,
        backend=None, seed=SEED, top_n=TOP_N):
    phase1_vars = phase1_vars or PHASE1_VARS
    dtmin_sweep = dtmin_sweep or DTMIN_SWEEP_K
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    # ---- Phase 1: 熱統合なしで設計を座標降下最適化 → D* ----
    def _phase1_opts(_p):
        return dict(apply_hi=False, hi_dT_min_K=POSTHOC_DTMIN_K, apply_stage2=False)

    n_total_p1 = sum(len(space.grid_points(v, k)) for v in phase1_vars)
    print("=" * 72, flush=True)
    print(f"  case_p02_pinch pinch_posthoc (P02) Phase1: 熱統合なしで設計最適化 "
          f"({len(phase1_vars)} 変数 × ~{k} 点 = ~{n_total_p1} 評価)", flush=True)
    print(f"  順序: {' → '.join(phase1_vars)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    obj1 = harness.make_objective(backend=backend, eval_opts=_phase1_opts)
    study1 = harness.new_study(seed)
    cb1 = reporting.make_callback(n_total_p1)

    phase1_rows = []
    for step, var in enumerate(phase1_vars, 1):
        pts = space.grid_points(var, k)
        points = [dict(current, **{var: val}) for val in pts]
        print(f"\n--- Phase1 step {step}/{len(phase1_vars)}: '{var}' を {len(pts)} 点掃引 ---", flush=True)
        trials = harness.run_batch(study1, obj1, points, cb1)
        for t, val in zip(trials, pts):
            phase1_rows.append({
                'step': step, 'var': var, 'value': val,
                'effective_TAC': t.value,
                'feasible': t.user_attrs.get('is_feasible'),
                'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                'production_kmol_h': t.user_attrs.get('production_kmol_h'),
                'trial_number': t.number,
            })
        b = harness.best_of(trials)
        if b is not None and var in b.params:
            current[var] = space.clamp(var, b.params[var])
            print(f"    → '{var}' = {current[var]} を採用 (TAC_noHI={b.value:.2f})", flush=True)

    d_star = {name: current[name] for name in space.PARAM_NAMES}  # 確定設計

    # ---- Phase 2: D* を固定し ΔTmin を掃引 (後置熱統合) ----
    def _phase2_opts(p):
        return dict(apply_hi=True, hi_dT_min_K=float(p['hi_dT_min_K']), apply_stage2=True)

    lo, hi = min(dtmin_sweep), max(dtmin_sweep)
    obj2 = harness.make_objective(backend=backend,
                                  control_space={'hi_dT_min_K': (lo, hi)},
                                  eval_opts=_phase2_opts)
    study2 = harness.new_study(seed)
    cb2 = reporting.make_callback(len(dtmin_sweep))
    print(f"\n{'='*72}\n  Phase2: 確定設計 D* を固定し ΔTmin を {len(dtmin_sweep)} 点掃引 "
          f"(後置熱統合)\n{'='*72}", flush=True)

    points2 = [dict(d_star, hi_dT_min_K=dt) for dt in dtmin_sweep]
    trials2 = harness.run_batch(study2, obj2, points2, cb2)

    dtmin_rows = []
    for t, dt in zip(trials2, dtmin_sweep):
        dtmin_rows.append({
            'hi_dT_min_K': dt,
            'effective_TAC': t.value,
            'feasible': t.user_attrs.get('is_feasible'),
            'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
            'production_kmol_h': t.user_attrs.get('production_kmol_h'),
            'is_posthoc_default': (dt == POSTHOC_DTMIN_K),
            'trial_number': t.number,
        })

    # 慣用 10K と掃引最良の差 (後置 ΔTmin 選択の損失)。
    posthoc = next((r for r in dtmin_rows if r['is_posthoc_default']), None)
    best_dt = min((r for r in dtmin_rows if r['effective_TAC'] is not None),
                  key=lambda r: r['effective_TAC'], default=None)
    extra = [
        "## 手法固有の出力",
        "- `dtmin_sweep.csv` … 確定設計 D* を後付け ΔTmin 各点で評価。",
        "- `phase1_no_hi.csv` … Phase1 (熱統合なし) の設計最適化の掃引軌跡。",
        "",
        "## 後置熱統合の損失",
        "P02 の本質的損失は **ループ内同時熱統合 (BO) との ΔTAC** (best.json と "
        "outputs/special_* を突合)。",
    ]
    if posthoc and best_dt and posthoc['effective_TAC'] is not None and best_dt['effective_TAC'] is not None:
        gap = posthoc['effective_TAC'] - best_dt['effective_TAC']
        extra.append(
            f"参考: 後置慣用 ΔTmin={POSTHOC_DTMIN_K:.0f}K の TAC={posthoc['effective_TAC']:.2f} に対し、"
            f"掃引最良 ΔTmin={best_dt['hi_dT_min_K']:.0f}K で TAC={best_dt['effective_TAC']:.2f} "
            f"(ΔTmin 固定の損失 ≈ {gap:.2f} 億円/年)。"
        )

    settings = {
        'phase1_vars': ' → '.join(phase1_vars),
        'k_points': k,
        'dtmin_sweep_K': dtmin_sweep,
        'posthoc_dtmin_K': POSTHOC_DTMIN_K,
        'n_total_evals': n_total_p1 + len(dtmin_sweep),
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'note': 'Phase1=apply_hi:False で設計確定 → Phase2=後置でΔTmin掃引',
    }
    # 成果物の best.json/top-N は Phase2 (後置熱統合した実プロセス) を対象にする。
    out_dir, best = harness.finalize(
        study2, method='pinch_posthoc',
        p_codes='P02 熱統合の後置',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        dtmin_rows, os.path.join(out_dir, 'dtmin_sweep.csv'),
        ['hi_dT_min_K', 'effective_TAC', 'feasible', 'purity_wt',
         'production_kmol_h', 'is_posthoc_default', 'trial_number'],
    )
    harness.save_table_csv(
        phase1_rows, os.path.join(out_dir, 'phase1_no_hi.csv'),
        ['step', 'var', 'value', 'effective_TAC', 'feasible',
         'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'dtmin_sweep.csv')} / phase1_no_hi.csv", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
