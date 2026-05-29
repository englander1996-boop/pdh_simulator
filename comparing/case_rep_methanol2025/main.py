r"""
comparing/case_rep_methanol2025/main.py — 実在レポート「CO2水素化によるメタノール製造 2025」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_methanol2025\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: CO2水素化→メタノール、2025、C04+C01、成熟度 5/12、
  検出 P04/P05/P06/P07/P08/P11/P12):
  そのレポートが実際にやった手法 (§3.4 最適化方針、原文より):
  反応器の **入口温度・圧力・体積** を最適化変数に (総括収率 ≥0.98 制約下)、評価関数=
  触媒+装置+加圧コストを 1 次元的に最適化。続いて蒸留塔は **段数 40 固定→フィード圧を最適化**、
  最終 42 段/フィード12段/還流比0.72。気液分離工程も別途最適化。

手法はテーマ非依存。PDH sim へは「反応器(入口温度) → 塔(フィード圧/段数) を 1 つずつ逐次」に写す:
  T_in_K → col2_p_kpa → col2_n_stages の逐次最適化。

体現する欠陥: P04 1次元逐次 / P05 整数の連続扱い / P06 / P12。成熟度 5/12 の中位の対照例。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P07 / P08 / P11
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


# 成熟度 5/12。反応器(入口温度)→塔(フィード圧→段数) の逐次1次元。
VARS_ORDER = ['T_in_K', 'col2_p_kpa', 'col2_n_stages']
K_POINTS = 5
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '5/12'
PROBLEMS = 'P04/P05/P06/P07/P08/P11/P12'

_DESCRIPTION = """\
実在レポート メタノール2025 (成熟度5/12) の手法の忠実再現。反応器(入口温度/圧力/体積)を 1 次元最適化後、
塔はフィード圧→段数を順に最適化 (逐次1次元)。PDH へは T_in_K→col2_p_kpa→col2_n_stages の逐次に写す。
P04/P05/P06/P12 を BO best との ΔTAC で定量化。
"""


def run(*, vars_order=None, k=K_POINTS, start=None, backend=None, seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND
    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    n_total = sum(len(space.grid_points(v, k)) for v in vars_order)
    print("=" * 72, flush=True)
    print(f"  case_rep_methanol2025 (実在レポート再現, メタノール, 成熟度{MATURITY}): {' → '.join(vars_order)} を逐次最適化", flush=True)
    print(f"  ~{n_total} 評価 / backend: {backend}", flush=True)
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
            cost_curves.append({'step': step, 'var': var, 'value': val, 'effective_TAC': t.value,
                                'feasible': t.user_attrs.get('is_feasible'),
                                'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                                'production_kmol_h': t.user_attrs.get('production_kmol_h'),
                                'trial_number': t.number})
        b = harness.best_of(trials)
        if b is not None and var in b.params:
            current[var] = space.clamp(var, b.params[var])
            print(f"    → '{var}' = {current[var]} を採用 (TAC={b.value:.2f})", flush=True)

    settings = {'source_report': 'メタノール レポート 2025 (C04+C01)', 'maturity': MATURITY, 'problems': PROBLEMS,
                'method': '反応器(入口温度)→塔(フィード圧→段数) の逐次1次元',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 出典レポートと手法",
        f"- メタノール 2025 (成熟度 {MATURITY}、検出 {PROBLEMS})。反応器(温度/圧力/体積)→塔(フィード圧/段数) 逐次。",
        "- PDH 写像: T_in_K→col2_p_kpa→col2_n_stages の逐次1次元。",
        "- 体現する欠陥: P04 1次元逐次 / P05 / P06 / P12。",
    ]
    out_dir, best = harness.finalize(study, method='rep_methanol2025_mat5',
        p_codes='メタノール2025 成熟度5/12 (P04/P05/P06/P12)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n)
    harness.save_table_csv(cost_curves, os.path.join(out_dir, 'cost_curves.csv'),
        ['step', 'var', 'value', 'effective_TAC', 'feasible', 'purity_wt', 'production_kmol_h', 'trial_number'])
    print(f"  手法固有: {os.path.join(out_dir, 'cost_curves.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
