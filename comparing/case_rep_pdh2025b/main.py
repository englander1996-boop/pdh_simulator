r"""
comparing/case_rep_pdh2025b/main.py — 実在 PDH レポート (2025, C04+C01) の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_pdh2025b\main.py`
  (出典著者名は repo に残さない方針。対応は Claude メモリのみ)

出典 (PDH=プロパン脱水素→プロピレン、2025、成熟度 4/12、検出 P03/P05/P06/P07/P08/P09/P11/P12):
  そのレポートが実際にやった手法 (§6.2.2 脱ブタン塔の最適化、原文より):
  塔建設コスト+用役コストの合計を評価関数とし、**段数だけ** を 1 次元掃引して最小段数を採用
  (フィード段は還流比最小の段)。プロピレン精製は膜分離+蒸留を併用。
  報告された最適: 実段数 17 / フィード段 15 / 還流比 0.30 / 塔高 10.0m / 塔径 4.22m。

手法はテーマ非依存。PDH sim では「1 塔の段数のみを 1 次元掃引」を Dist1 段数に写す:
  col1_n_stages を単独で掃引 (他は固定)。膜分離は A_mem に対応 (本 case では段数掃引に集中)。

体現する欠陥: P05 整数の連続扱い (段数を粗く掃引するだけ) / P06 / P12。
1 変数しか動かさないため P04 (変数間相互作用の取りこぼし) も最も極端な形で現れる。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P03 / P07 / P08 / P09 / P11
→ 低成熟度ほど ◎ を多く束ね BO 比 ΔTAC が大きい。◎ の複合損失のみ定量化する。

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
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
# 設定 — 「1 塔の段数のみ 1 次元掃引」をそのまま写す
# ===========================================================================
VARS_ORDER = ['col1_n_stages']   # 段数だけを掃引 (報告: 最小段数 17 段)
K_POINTS = 7                      # 段数の掃引点数 (1 変数なので少し細かめ)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND

REPORTED = {'n_stages': 17, 'feed_stage': 15, 'reflux': 0.30, 'height_m': 10.0, 'diameter_m': 4.22}

_DESCRIPTION = """\
実在 PDH レポート (2025) の最適化手法の忠実再現。塔の段数だけを 1 次元掃引して最小コスト段数を
採用する (他変数は固定、単一始点1巡)。1 変数しか動かさないため P04/P05/P06/P12 が最も極端に現れる。
BO best との ΔTAC、および報告最適 (17段/RR0.30) との突合で定量化する。
"""


def run(*, vars_order=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    n_total = sum(len(space.grid_points(v, k)) for v in vars_order)
    print("=" * 72, flush=True)
    print(f"  case_rep_pdh2025b (実在レポート再現, PDH): 塔段数のみ 1 次元掃引 ({', '.join(vars_order)})", flush=True)
    print(f"  ~{n_total} 評価 / backend: {backend} / 報告最適(参考): 17段/RR0.30", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    cost_curves = []
    for step, var in enumerate(vars_order, 1):
        pts = space.grid_points(var, k)
        points = [dict(current, **{var: val}) for val in pts]
        print(f"\n--- step {step}/{len(vars_order)}: '{var}' を {len(pts)} 点掃引 ---", flush=True)
        trials = harness.run_batch(study, objective, points, cb)
        for t, val in zip(trials, pts):
            cost_curves.append({
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
            print(f"    → '{var}' = {current[var]} を採用 (TAC={b.value:.2f})", flush=True)

    settings = {
        'source_report': 'PDH レポート 2025 (C04+C01, 成熟度4/12)',
        'method': '塔の段数のみを 1 次元掃引 (単一始点1巡)',
        'vars_order': ' → '.join(vars_order),
        'k_points': k, 'n_total_evals': n_total, 'backend': backend,
        'reported_optimum': REPORTED,
    }
    extra = [
        "## 出典レポートと手法",
        "- PDH レポート (2025): 脱ブタン塔の段数だけを 1 次元掃引して最小段数を採用。",
        "- 報告最適: 実段数17 / フィード段15 / 還流比0.30 / 塔高10.0m / 塔径4.22m。",
        "- 体現する欠陥: P05 整数の連続扱い / P06 / P12 (+ 1 変数のみで P04 が極端)。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC。報告最適値とも sim 上で突合可。",
    ]
    out_dir, best = harness.finalize(
        study, method='pdh2025b',
        p_codes='PDH2025b (1塔段数のみ→P05/P06/P12)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        cost_curves, os.path.join(out_dir, 'cost_curves.csv'),
        ['step', 'var', 'value', 'effective_TAC', 'feasible',
         'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'cost_curves.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
