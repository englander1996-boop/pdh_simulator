r"""
comparing/case_rep_ammonia2025/main.py — 実在レポート「グリーンアンモニア製造 2025」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_ammonia2025\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: グリーンアンモニア (N2+H2→NH3)、2025、C04+C01、成熟度 4/12、
  検出 P03/P04/P05/P06/P07/P08/P11/P12):
  そのレポートが実際にやった手法 (§3.5 最適化 / パージ率の最適化、原文より):
  パージ率を最適化しつつ、多段断熱反応器の入口温度・触媒層体積などを 1 次元的に最適化。

手法はテーマ非依存。ただし **パージ率は PDH sim に対応変数が無い (P03 の温床=未露出)**。
写せるのは反応器側の 1 次元最適化なので、入口温度→反応器径を逐次に写す:
  T_in_K → D_reactor_m の逐次1次元 (パージ率の最適化部分は省略、と明記)。

## 含まれる欠陥部品 (このレポートが体現する P と、対応する case_p## 部品)
- sim で再現される: **P04** 1次元逐次 (→ case_p04_sequential), **P05** 整数の連続扱い (→ case_p05_grid),
  **P06** 大域性未保証 (→ case_p06_multistart), **P12** 検証なし (→ case_p12_converge)。
- 検出されたが sim 非対応 (写せない): **P03**(パージ率=未露出。皮肉にも本レポートはパージ率を
  最適化しているが、それでも recycle 流量自体は固定とされ P03 判定), P07(単目的), P08(NPV), P11(動的)。

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


VARS_ORDER = ['T_in_K', 'D_reactor_m']
K_POINTS = 5
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '4/12'
PROBLEMS_DETECTED = 'P03/P04/P05/P06/P07/P08/P11/P12'
PROBLEMS_REPRODUCED = 'P04/P05/P06/P12'

_DESCRIPTION = """\
実在レポート グリーンアンモニア2025 (成熟度4/12) の手法の忠実再現。パージ率＋反応器を 1 次元最適化。
パージ率は PDH sim 未露出のため反応器側 (入口温度→反応器径) の逐次に写す。含む欠陥部品:
P04(1次元逐次)/P05/P06/P12。検出 P03(パージ)/P07/P08/P11 は sim 非対応。BO との ΔTAC で定量化。
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
    print(f"  case_rep_ammonia2025 (実在レポート再現, アンモニア, 成熟度{MATURITY}): {' → '.join(vars_order)} を逐次最適化", flush=True)
    print(f"  含む欠陥部品: {PROBLEMS_REPRODUCED} (検出 {PROBLEMS_DETECTED}、パージ率は sim 未露出で省略) / backend: {backend}", flush=True)
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

    settings = {'source_report': 'グリーンアンモニア レポート 2025 (C04+C01)', 'maturity': MATURITY,
                'problems_detected': PROBLEMS_DETECTED, 'problems_reproduced_on_sim': PROBLEMS_REPRODUCED,
                'method': 'パージ率+反応器1次元 (パージ率は未露出→反応器 T_in→D の逐次に写す)',
                'unmapped': 'パージ率 (P03) は sim 未露出のため省略',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 含まれる欠陥部品 (このレポートが束ねる P)",
        f"- sim で再現: **{PROBLEMS_REPRODUCED}** = P04(→case_p04_sequential)/P05(→case_p05_grid)/P06(→case_p06_multistart)/P12(→case_p12_converge)。",
        "- 検出だが sim 非対応: P03(パージ率=未露出)/P07(単目的)/P08(NPV)/P11(動的)。",
        f"- 出典検出 P: {PROBLEMS_DETECTED}、成熟度 {MATURITY}。",
        "## 定量化",
        "本手法 best.json と BO を突合し ΔTAC = P04/P05/P06/P12 複合の損失。",
    ]
    out_dir, best = harness.finalize(study, method='rep_ammonia2025_mat4',
        p_codes='アンモニア2025 成熟度4/12 (含 P04/P05/P06/P12)',
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
