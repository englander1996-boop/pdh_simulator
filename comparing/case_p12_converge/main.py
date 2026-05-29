r"""
comparing/case_p12_converge/main.py — P12「収束・頑健性検証なし (1 巡で打ち切り)」の忠実再現。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p12_converge\main.py`

再現する問題のあるやり方:
  case_p04_sequential (1次元逐次) のような素朴手法は、座標降下を **1 巡** 回した時点で打ち切り、
  「2 巡目を回したら結果が変わらないか (= 本当に収束したか)」を検証しない (P12)。
  Aspen の『収束済み』ステータスを最適化の収束と取り違えるのが典型。

このスクリプトが示すこと (定量化のねらい):
  同じ座標降下を **複数巡 (既定 3 巡)** 回し、各巡終了時点の best TAC を記録する。
  - 2 巡目以降も best が更新される = 1 巡で打ち切った素朴手法は未収束の点を「答え」と誤認している。
  - 巡ごとの改善量 (pass1_best − pass2_best …) = 「検証していれば回収できた分」(P12 の損失)。
  (始点感度=大域性 P06 は case_p06_multistart で別途扱う)

出力 (comparing/results/converge_check_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式、全巡の全 trial)
  pass_curves.csv : (pass, step, var, 掃引値, effective_TAC, ...) 各巡の掃引軌跡
  pass_best.csv   : (pass, best_so_far_TAC, improved_from_prev) 巡ごとの best 推移

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
# case_p04_sequential と同じコストドライバを座標降下する (1 巡で打ち切らず複数巡して収束を「検証」する)。
VARS_ORDER = [
    'col2_p_kpa',
    'col2_n_stages',
    'col2_reflux_ratio',
    'col3_n_stages',
    'col3_p_kpa',
    'F_C3H8_fresh_kmol_h',
    'T_in_K',
]
N_PASSES = 3               # 巡回数 (素朴手法は 1 巡で打ち切る。ここでは検証のため複数巡)
K_POINTS = 5               # 1 変数あたり掃引点数
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P12 検証なし の忠実再現。case_p04_sequential と同じ座標降下を複数巡し、各巡終了時点の best TAC を記録する。
2 巡目以降も best が更新されるなら、1 巡で打ち切った素朴手法は未収束の点を「答え」と誤認している。
巡ごとの改善量が「収束を検証していれば回収できた分」= P12 の損失。pass_best.csv に巡ごとの推移を記録。
"""


def run(*, vars_order=None, n_passes=N_PASSES, k=K_POINTS, start=None,
        backend=None, seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    per_pass_evals = sum(len(space.grid_points(v, k)) for v in vars_order)
    n_total = n_passes * per_pass_evals

    print("=" * 72, flush=True)
    print(f"  case_p12_converge converge_check (P12): {n_passes} 巡 × 座標降下({len(vars_order)}変数×~{k}点) "
          f"= ~{n_total} 評価", flush=True)
    print(f"  順序: {' → '.join(vars_order)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    pass_curves = []
    pass_best_rows = []
    best_so_far = float('inf')
    for p in range(1, n_passes + 1):
        print(f"\n==== 巡 {p}/{n_passes} ====", flush=True)
        for step, var in enumerate(vars_order, 1):
            pts = space.grid_points(var, k)
            points = [dict(current, **{var: val}) for val in pts]
            trials = harness.run_batch(study, objective, points, cb)
            for t, val in zip(trials, pts):
                pass_curves.append({
                    'pass': p, 'step': step, 'var': var, 'value': val,
                    'effective_TAC': t.value,
                    'feasible': t.user_attrs.get('is_feasible'),
                    'trial_number': t.number,
                })
            b = harness.best_of(trials)
            if b is not None and var in b.params:
                current[var] = space.clamp(var, b.params[var])

        # この巡終了時点の best (feasible 優先, TAC 最小) を全 trial から取る。
        _comp, _feas, gbest = reporting.summarize(study)
        cur = gbest.value if gbest is not None else float('inf')
        improved = (best_so_far - cur) if (best_so_far < float('inf')) else None
        pass_best_rows.append({
            'pass': p,
            'best_so_far_TAC': (cur if cur < float('inf') else None),
            'improved_from_prev': (improved if improved is not None else ''),
            'best_feasible': (gbest.user_attrs.get('is_feasible') if gbest is not None else None),
        })
        if cur < float('inf'):
            msg = f"    巡{p} 終了時 best TAC = {cur:.2f}"
            if improved is not None:
                msg += f" (前巡から {improved:+.2f} 改善)" if improved > 1e-9 else " (前巡から変化なし=収束)"
            print(msg, flush=True)
        best_so_far = min(best_so_far, cur)

    # P12 の損失: 1 巡で打ち切った場合 (pass1) と複数巡後 (最終) の差。
    p1 = next((r['best_so_far_TAC'] for r in pass_best_rows if r['pass'] == 1), None)
    pN = next((r['best_so_far_TAC'] for r in pass_best_rows if r['pass'] == n_passes), None)
    extra = [
        "## 手法固有の出力",
        "- `pass_best.csv` … 各巡終了時の best TAC 推移 (2 巡目以降の改善が「1 巡打ち切り」の見落とし)。",
        "- `pass_curves.csv` … 各巡・各 step の掃引軌跡。",
        "",
        "## 収束未検証の損失 (P12)",
    ]
    if p1 is not None and pN is not None:
        gap = p1 - pN
        if gap > 1e-9:
            extra.append(
                f"1 巡で打ち切った場合の best={p1:.2f} に対し、{n_passes} 巡後 best={pN:.2f} "
                f"→ 収束を検証していれば回収できた分 ≈ {gap:.2f} 億円/年 (= 1 巡打ち切りの損失)。"
            )
        else:
            extra.append(
                f"1 巡 best={p1:.2f} と {n_passes} 巡後 best={pN:.2f} は一致 (この設定では 1 巡で収束)。"
                "ただし素朴手法はこの確認自体を行わない点が P12 の本質。"
            )

    settings = {
        'vars_order': ' → '.join(vars_order),
        'n_passes': n_passes,
        'k_points': k,
        'n_total_evals': n_total,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'start': 'midpoint' if not start else 'custom',
    }
    out_dir, best = harness.finalize(
        study, method='converge_check',
        p_codes='P12 検証なし (収束未確認)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        pass_curves, os.path.join(out_dir, 'pass_curves.csv'),
        ['pass', 'step', 'var', 'value', 'effective_TAC', 'feasible', 'trial_number'],
    )
    harness.save_table_csv(
        pass_best_rows, os.path.join(out_dir, 'pass_best.csv'),
        ['pass', 'best_so_far_TAC', 'improved_from_prev', 'best_feasible'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'pass_best.csv')} / pass_curves.csv", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
