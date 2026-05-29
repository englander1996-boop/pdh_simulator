r"""
comparing/case_combo_typical/main.py — 「典型的な学部レポート」の忠実再現 (C×P の意味ある組み合わせ)。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_combo_typical\main.py`

再現する documented practice (最適化手法分析 C01+C04。コーパスの 96.7% がこの 2 手法):
  実在レポートの大多数は、次の手順を 1 本のワークフローで実行している。
    Phase1: C01 1次元逐次最適化 — 熱統合は考えず、コストドライバを 1 変数ずつ粗い格子で
            掃引→最小採用 (座標降下) を、単一始点から 1 巡だけ回して打ち切る。
    Phase2: C04 ピンチ解析 — 確定した設計に後付けで熱統合を当てる (ΔTmin=10K 慣用固定)。

この 1 本が同時に体現する欠陥 (= C×P の意味ある組み合わせ):
  - P04 1次元逐次   : 相互作用する変数を 1 次元ずつ最適化 (Phase1)
  - P05 整数の連続扱い: 段数 N を粗い格子で掃引するだけ (Phase1)
  - P06 大域性未保証 : 単一始点の座標降下 (Phase1)
  - P12 検証なし     : 1 巡で打ち切り、2 巡目/始点感度を確認しない (Phase1)
  - P02 熱統合の後置 : プロセス確定後に後付けピンチ、ΔTmin も慣用値固定 (Phase2)

なぜ問題か (定量化のねらい):
  この「典型レポート」設計の最終 TAC を、BO (全変数 + 熱統合を同時最適化) の best と突合した
  ΔTAC が、「学部標準ワークフロー全体の損失」。個別欠陥 (case_p01_subsystem/2/4/5/6/12) の損失の積み上げが
  ここに集約される。

出力 (comparing/results/typical_workflow_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式、Phase2=後置HI 結果)
  phase1_cost_curves.csv : Phase1 (C01 逐次, HIなし) の 1 次元掃引軌跡

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
# Phase1 (C01 1次元逐次, 熱統合なし) で掃引する変数と順序 (上流→下流)。
VARS_ORDER = [
    'col2_p_kpa',
    'col2_n_stages',
    'col2_reflux_ratio',
    'col3_n_stages',
    'col3_p_kpa',
    'F_C3H8_fresh_kmol_h',
    'T_in_K',
]
K_POINTS = 5               # Phase1 の 1 変数あたり掃引点数 (粗い = 過去手法に忠実)
POSTHOC_DTMIN_K = 10.0     # Phase2 後置ピンチの慣用 ΔTmin (学生がそのまま採る値)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
典型的な学部レポート (C01 1次元逐次 → C04 後置ピンチ) の忠実再現。Phase1 で熱統合なしに
コストドライバを単一始点・1 巡の座標降下で最適化し (P04/P05/P06/P12)、Phase2 で確定設計に
後付けピンチを ΔTmin=10K 固定で当てる (P02)。この学部標準ワークフロー全体の TAC を、
BO (同時最適化) best との ΔTAC で定量化する。
"""


def run(*, vars_order=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    # ---- Phase 1: C01 1次元逐次 (熱統合なし、単一始点 1 巡) → D* ----
    def _phase1_opts(_p):
        return dict(apply_hi=False, hi_dT_min_K=POSTHOC_DTMIN_K, apply_stage2=False)

    n_p1 = sum(len(space.grid_points(v, k)) for v in vars_order)
    print("=" * 72, flush=True)
    print(f"  case_combo_typical (C01+C04): 典型レポート再現", flush=True)
    print(f"  Phase1 = C01 1次元逐次 (HIなし, 単一始点1巡): {len(vars_order)} 変数 × ~{k} 点 = ~{n_p1} 評価", flush=True)
    print(f"  順序: {' → '.join(vars_order)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    obj1 = harness.make_objective(backend=backend, eval_opts=_phase1_opts)
    study1 = harness.new_study(seed)
    cb1 = reporting.make_callback(n_p1)

    cost_curves = []
    for step, var in enumerate(vars_order, 1):
        pts = space.grid_points(var, k)
        points = [dict(current, **{var: val}) for val in pts]
        print(f"\n--- Phase1 step {step}/{len(vars_order)}: '{var}' を {len(pts)} 点掃引 (他固定) ---", flush=True)
        trials = harness.run_batch(study1, obj1, points, cb1)
        for t, val in zip(trials, pts):
            cost_curves.append({
                'step': step, 'var': var, 'value': val,
                'effective_TAC_noHI': t.value,
                'feasible': t.user_attrs.get('is_feasible'),
                'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                'production_kmol_h': t.user_attrs.get('production_kmol_h'),
                'trial_number': t.number,
            })
        b = harness.best_of(trials)
        if b is not None and var in b.params:
            current[var] = space.clamp(var, b.params[var])
            print(f"    → '{var}' = {current[var]} を採用 (TAC_noHI={b.value:.2f})", flush=True)

    d_star = {name: current[name] for name in space.PARAM_NAMES}
    b1 = harness.best_of(study1.trials)
    tac_no_hi = b1.value if b1 is not None else None

    # ---- Phase 2: C04 後置ピンチ (ΔTmin=10K 固定) を D* に当てる ----
    def _phase2_opts(_p):
        return dict(apply_hi=True, hi_dT_min_K=POSTHOC_DTMIN_K, apply_stage2=True)

    obj2 = harness.make_objective(backend=backend, eval_opts=_phase2_opts)
    study2 = harness.new_study(seed)
    cb2 = reporting.make_callback(1)
    print(f"\n{'='*72}\n  Phase2 = C04 後置ピンチ: 確定設計 D* に ΔTmin={POSTHOC_DTMIN_K:.0f}K を後付け\n{'='*72}", flush=True)
    trials2 = harness.run_batch(study2, obj2, [dict(d_star)], cb2)
    b2 = harness.best_of(trials2)
    tac_hi = b2.value if b2 is not None else None

    extra = [
        "## このワークフローが同時に体現する欠陥",
        "- P04 1次元逐次 / P05 整数の連続扱い / P06 大域性未保証 / P12 検証なし (Phase1)",
        "- P02 熱統合の後置 (Phase2、ΔTmin=10K 固定)",
        "",
        "## 手法固有の出力",
        "- `phase1_cost_curves.csv` … Phase1 (C01 逐次・HIなし) の 1 次元掃引軌跡。",
    ]
    if tac_no_hi is not None and tac_hi is not None:
        extra += [
            "",
            "## 段階別 TAC",
            f"- Phase1 確定 (熱統合なし) effective_TAC = {tac_no_hi:.2f} 億円/年",
            f"- Phase2 後置ピンチ後 effective_TAC = {tac_hi:.2f} 億円/年",
            "- 本質的損失は **BO (全変数+熱統合 同時最適化) best との ΔTAC** (best.json と outputs/special_* を突合)。",
        ]

    settings = {
        'phase1_vars_order': ' → '.join(vars_order),
        'k_points': k,
        'posthoc_dtmin_K': POSTHOC_DTMIN_K,
        'n_total_evals': n_p1 + 1,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'workflow': 'C01 1次元逐次(HIなし,単一始点1巡) → C04 後置ピンチ(ΔTmin=10K)',
    }
    # best.json/top-N は Phase2 (後置ピンチ済みの最終プロセス = レポートの提出設計) を対象にする。
    out_dir, best = harness.finalize(
        study2, method='typical_workflow',
        p_codes='C01+C04 典型ワークフロー (P02/P04/P05/P06/P12)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        cost_curves, os.path.join(out_dir, 'phase1_cost_curves.csv'),
        ['step', 'var', 'value', 'effective_TAC_noHI', 'feasible',
         'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'phase1_cost_curves.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
