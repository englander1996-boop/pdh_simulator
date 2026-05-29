r"""
comparing/case_rep_mek2019/main.py — 実在レポート「2-ブタノールからの MEK 製造 2019」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_mek2019\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: 2-ブタノール脱水素→メチルエチルケトン(MEK)、2019、C01、成熟度 5/12、
  検出 P03/P05/P06/P08/P09/P11/P12):
  そのレポートが実際にやった手法 (§3.8 反応器最適化、原文より):
  反応器を **内径 D (1,2,4 inch) × 熱媒流量 Fh (1000–15000, 100刻み) = 4200 通りのグリッド** で
  総当たり最適化 (反応率・入口温度は全体最適化で別途決定)。報告: 反応率0.98/入口温度700K/D=4inch/Fh=6700。

手法はテーマ非依存。PDH sim へは「反応器の 2 変数を粗グリッド総当たり」に写す:
  反応器径 D_reactor_m × 触媒層厚 z_cat_m の粗グリッド (熱媒流量は PDH に対応変数なし→触媒層厚に置換)。

## 含まれる欠陥部品 (このレポートが体現する P と、対応する case_p## 部品)
- sim で再現される: **P05** 整数/離散の連続扱い・粗グリッド (→ case_p05_grid),
  **P06** 単一グリッドで大域性なし (→ case_p06_multistart), **P12** 検証なし (→ case_p12_converge)。
- 検出されたが sim 非対応 (写せない): P03(パージ/リサイクル不変数化), P08(経済単純化=NPV),
  P09(環境/CO2), P11(定常のみ)。
→ この case 単体では P05/P06/P12 の複合損失を、case_p05_grid 等の単体部品と対比して定量化できる。

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
"""

import os
import sys
import itertools

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

from comparing.shared import space, simulator, reporting, harness


GRID_VARS = ['D_reactor_m', 'z_cat_m']
K_POINTS = 4
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '5/12'
PROBLEMS_DETECTED = 'P03/P05/P06/P08/P09/P11/P12'
PROBLEMS_REPRODUCED = 'P05/P06/P12'  # sim で再現できる部分

_DESCRIPTION = """\
実在レポート MEK2019 (成熟度5/12) の手法の忠実再現。反応器の 内径×熱媒流量 を 4200 点グリッドで総当たり。
PDH へは反応器径×触媒層厚の粗グリッドに写す。含む欠陥部品: P05(粗い離散総当たり)/P06(大域性なし)/P12(検証なし)。
検出 P03/P08/P09/P11 は sim 非対応。BO best との ΔTAC で複合損失を定量化。
"""


def run(*, grid_vars=None, k=K_POINTS, start=None, backend=None, seed=SEED, top_n=TOP_N):
    grid_vars = grid_vars or GRID_VARS
    backend = backend or BACKEND
    base = space.midpoint_params()
    if start:
        base.update(start)
    base = {name: space.clamp(name, base[name]) for name in space.PARAM_NAMES}

    axes = [space.grid_points(v, k) for v in grid_vars]
    combos = list(itertools.product(*axes))
    n_total = len(combos)
    print("=" * 72, flush=True)
    print(f"  case_rep_mek2019 (実在レポート再現, MEK, 成熟度{MATURITY}): "
          f"反応器 {' × '.join(f'{v}[{len(a)}]' for v, a in zip(grid_vars, axes))} = {n_total} 点 粗グリッド", flush=True)
    print(f"  含む欠陥部品: {PROBLEMS_REPRODUCED} (検出 {PROBLEMS_DETECTED}) / backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    points = [dict(base, **{v: val for v, val in zip(grid_vars, combo)}) for combo in combos]
    trials = harness.run_batch(study, objective, points, cb)

    grid_rows = []
    for t, combo in zip(trials, combos):
        row = {f'var_{v}': val for v, val in zip(grid_vars, combo)}
        row.update({'effective_TAC': t.value, 'feasible': t.user_attrs.get('is_feasible'),
                    'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                    'production_kmol_h': t.user_attrs.get('production_kmol_h'), 'trial_number': t.number})
        grid_rows.append(row)

    settings = {'source_report': 'MEK レポート 2019 (C01)', 'maturity': MATURITY,
                'problems_detected': PROBLEMS_DETECTED, 'problems_reproduced_on_sim': PROBLEMS_REPRODUCED,
                'method': '反応器 内径×熱媒流量 の 4200 点グリッド (PDH: 径×触媒層厚に写す)',
                'grid_vars': ' × '.join(grid_vars), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 含まれる欠陥部品 (このレポートが束ねる P)",
        f"- sim で再現: **{PROBLEMS_REPRODUCED}** = P05(→case_p05_grid)/P06(→case_p06_multistart)/P12(→case_p12_converge)。",
        f"- 検出だが sim 非対応: P03(パージ)/P08(NPV)/P09(CO2)/P11(動的)。",
        f"- 出典検出 P (問題点レポート): {PROBLEMS_DETECTED}、成熟度 {MATURITY}。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC = P05/P06/P12 複合の損失。",
    ]
    out_dir, best = harness.finalize(study, method='rep_mek2019_mat5',
        p_codes='MEK2019 成熟度5/12 (含 P05/P06/P12)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n)
    harness.save_table_csv(grid_rows, os.path.join(out_dir, 'grid.csv'),
        [f'var_{v}' for v in grid_vars] + ['effective_TAC', 'feasible', 'purity_wt', 'production_kmol_h', 'trial_number'])
    print(f"  手法固有: {os.path.join(out_dir, 'grid.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
