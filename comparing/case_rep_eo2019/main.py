r"""
comparing/case_rep_eo2019/main.py — 実在レポート「エチレンオキシド製造 2019」の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_eo2019\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: エチレンオキシド(EO)、2019、C04+C01、**成熟度 5/12**、
  検出 P03/P05/P06/P07/P08/P11/P12):
  そのレポートが実際にやった手法 (§5.2 全体最適化、原文より):
  最適化変数 4 つ (O2/エチレン比・EO/エチレン比・吸収塔液入口温度・**パージ率**) を、
  「同時に動かすと膨大」として **1 つ固定→次を動かす座標降下** で順に最適化。
  蒸留部圧力は 3 通りの組合せを比較 (脱水塔/精留塔とも 1 atm を採用)。

手法はテーマ非依存。ただし O2 比・吸収塔・パージ率は PDH sim に対応変数が無い (パージ率=P03 未露出)。
**写せるのは「複数変数を 1 つずつ座標降下する構造」**なので、PDH のコストドライバ 3 変数に写す:
  col2_reflux_ratio → col3_n_stages → T_in_K を 1 つずつ掃引 (= 4 変数座標降下の構造を再現)。
  (O2比/吸収塔温度/パージ率は PDH に無いので省略、と明記)。

体現する欠陥: P04 1次元逐次(座標降下) / P05 / P06 / P12。成熟度 5/12 の中位の対照例。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P03 / P07 / P08 / P11
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


# 成熟度 5/12。原典は 4 変数座標降下 (O2比/EO比/吸収塔温度/パージ率)。PDH 写像 = コストドライバ 3 変数の座標降下。
VARS_ORDER = ['col2_reflux_ratio', 'col3_n_stages', 'T_in_K']
K_POINTS = 4
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '5/12'
PROBLEMS = 'P03/P05/P06/P07/P08/P11/P12'

_DESCRIPTION = """\
実在レポート EO2019 (成熟度5/12) の手法の忠実再現。原典は 4 変数 (O2比/EO比/吸収塔温度/パージ率) を
「同時は膨大」として 1 つずつ固定→次を動かす座標降下。PDH に無い変数 (O2比/吸収塔/パージ率) は省略し、
コストドライバ 3 変数の座標降下構造に写す。P04/P05/P06/P12 を BO best との ΔTAC で定量化。
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
    print(f"  case_rep_eo2019 (実在レポート再現, EO, 成熟度{MATURITY}): "
          f"4変数座標降下を {' → '.join(vars_order)} に写して逐次最適化", flush=True)
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

    settings = {'source_report': 'EO レポート 2019 (C04+C01)', 'maturity': MATURITY, 'problems': PROBLEMS,
                'method': '4 変数座標降下 (原典: O2比/EO比/吸収塔温度/パージ率) を PDH 3 変数に写す',
                'unmapped': 'O2比/吸収塔温度/パージ率 は PDH に対応変数なし→省略',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 出典レポートと手法",
        f"- EO 2019 (成熟度 {MATURITY}、検出 {PROBLEMS})。4 変数を 1 つずつ座標降下 (同時最適化は膨大として回避)。",
        "- PDH 写像: O2比/吸収塔温度/パージ率は対応変数なし→省略。座標降下の構造をコストドライバ 3 変数で再現。",
        "- 体現する欠陥: P04 1次元逐次 / P05 / P06 / P12。",
    ]
    out_dir, best = harness.finalize(study, method='rep_eo2019_mat5',
        p_codes='EO2019 成熟度5/12 (P04/P05/P06/P12)',
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
