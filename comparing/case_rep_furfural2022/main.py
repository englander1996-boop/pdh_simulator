r"""
comparing/case_rep_furfural2022/main.py — 実在レポート「キシロース脱水によるフルフラール製造 2022」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_furfural2022\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: 酢酸触媒・キシロース脱水→フルフラール、2022、C01、**成熟度 2/12 = 欠陥多**、
  検出 P01/P02/P03/P04/P05/P06/P07/P08/P11/P12 ← 10 個):
  そのレポートが実際にやった手法 (§5.1 反応器最適化、原文より):
  単通反応率・濃縮塔の D/F 比・還流比 (0.10–0.70 等) を 1 つずつ振って利益/コストを調べ最小点を採用。
  反応器と分離を分けて順に決める (部分最適化)。

手法はテーマ非依存。PDH sim へは「D/F比→還流比→反応率 を 1 つずつ逐次」に写す:
  col3_feed_ratio (留出/供給比) → col2_reflux_ratio (還流比) → T_in_K (反応率) の逐次1次元。

## 含まれる欠陥部品 (このレポートが体現する P と、対応する case_p## 部品)
- sim で再現される: **P01** 部分最適化 (反応器/分離を分けて順に→ case_p01_subsystem の性質),
  **P04** 1次元逐次 (→ case_p04_sequential), **P05** 整数の連続扱い (→ case_p05_grid),
  **P06** 大域性未保証 (→ case_p06_multistart), **P12** 検証なし (→ case_p12_converge),
  **P02** 熱統合後置 (→ case_p02_pinch の性質)。
- 検出されたが sim 非対応: P03(リサイクル/パージ), P07(単目的), P08(NPV), P11(動的)。
→ 低成熟度 (2/12) なので欠陥部品を最も多く束ねる = BO との ΔTAC が最大級の対照例。

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


VARS_ORDER = ['col3_feed_ratio', 'col2_reflux_ratio', 'T_in_K']
K_POINTS = 5
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '2/12 (欠陥多)'
PROBLEMS_DETECTED = 'P01/P02/P03/P04/P05/P06/P07/P08/P11/P12'
PROBLEMS_REPRODUCED = 'P01/P02/P04/P05/P06/P12'

_DESCRIPTION = """\
実在レポート フルフラール2022 (成熟度2/12=欠陥多) の手法の忠実再現。濃縮塔 D/F 比・還流比・反応率を
1 つずつ振って最小コスト点を採用 (反応器/分離を分けて順に=部分最適化)。PDH へは
col3_feed_ratio→col2_reflux_ratio→T_in_K の逐次に写す。含む欠陥部品: P01/P02/P04/P05/P06/P12 と多い。
低成熟度なので BO との ΔTAC が最大級になる対照例。
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
    print(f"  case_rep_furfural2022 (実在レポート再現, フルフラール, 成熟度{MATURITY}): {' → '.join(vars_order)} を逐次最適化", flush=True)
    print(f"  含む欠陥部品: {PROBLEMS_REPRODUCED} (検出 {PROBLEMS_DETECTED}) / backend: {backend}", flush=True)
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

    settings = {'source_report': 'フルフラール レポート 2022 (C01)', 'maturity': MATURITY,
                'problems_detected': PROBLEMS_DETECTED, 'problems_reproduced_on_sim': PROBLEMS_REPRODUCED,
                'method': 'D/F比→還流比→反応率 の逐次1次元 (反応器/分離を分けて=部分最適化)',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 含まれる欠陥部品 (このレポートが束ねる P) — 低成熟度で最多",
        f"- sim で再現: **{PROBLEMS_REPRODUCED}** = P01(→case_p01_subsystem)/P02(→case_p02_pinch)/P04(→case_p04_sequential)/P05(→case_p05_grid)/P06(→case_p06_multistart)/P12(→case_p12_converge)。",
        "- 検出だが sim 非対応: P03/P07/P08/P11。",
        f"- 出典検出 P: {PROBLEMS_DETECTED} (10個)、成熟度 {MATURITY}。",
        "## 定量化",
        "本手法 best.json と BO を突合し ΔTAC = 多数欠陥の複合損失 (最大級が期待される)。",
    ]
    out_dir, best = harness.finalize(study, method='rep_furfural2022_mat2',
        p_codes='フルフラール2022 成熟度2/12 (含 P01/P02/P04/P05/P06/P12)',
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
