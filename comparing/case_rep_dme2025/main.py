r"""
comparing/case_rep_dme2025/main.py — 実在レポート「CO2 を原料とした DME 製造 2025」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_dme2025\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: CO2→ジメチルエーテル(DME)、2025、Primary C04 / Secondary C02・C01、成熟度 4/12、
  検出 P02/P03/P05/P06/P07/P08/P11/P12):
  そのレポートが実際にやった手法 (§3.8 / §4 反応器最適化、原文より):
  多段断熱反応器について、**段数 (2–5 段) × 1 段目入口温度 (290–360℃) × 全圧** を変え、
  各段の入口温度は **総当たり (グリッド)** で最適化し、総コスト最小の組合せを採用。

手法はテーマ非依存。PDH sim には多段反応器の「段数」概念が無いので、反応器の転化率を支配する
レバー **入口温度 × 触媒層厚 (滞留時間)** の粗グリッドに写す (= 反応器条件を総当たりする C02 構造)。

体現する欠陥: P05 (条件を粗い格子で総当たり) / P06 (単一グリッドで大域性なし) / P02 / P12。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P02 / P03 / P07 / P08 / P11
→ 低成熟度ほど ◎ を多く束ね BO 比 ΔTAC が大きい。◎ の複合損失のみ定量化する。

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


# 成熟度 4/12。反応器の 入口温度 × 触媒層厚 の粗グリッド (原典: 段数×入口温度×圧力の総当たり)。
GRID_VARS = ['T_in_K', 'z_cat_m']
K_POINTS = 4
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '4/12'
PROBLEMS = 'P02/P03/P05/P06/P07/P08/P11/P12'

_DESCRIPTION = """\
実在レポート DME2025 (成熟度4/12, C02副) の手法の忠実再現。原典は多段反応器の 段数×入口温度×圧力 を
総当たりグリッド。PDH には段数概念が無いので反応器の 入口温度×触媒層厚 の粗グリッドに写す。
粗グリッド総当たりが大域最適を保証しない欠陥 (P05/P06) を BO best との ΔTAC で定量化。
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
    print(f"  case_rep_dme2025 (実在レポート再現, DME, 成熟度{MATURITY}): "
          f"反応器 {' × '.join(f'{v}[{len(a)}]' for v, a in zip(grid_vars, axes))} = {n_total} 点 粗グリッド", flush=True)
    print(f"  backend: {backend}", flush=True)
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

    settings = {'source_report': 'DME レポート 2025 (C04/C02/C01)', 'maturity': MATURITY, 'problems': PROBLEMS,
                'method': '反応器条件 (原典: 段数×入口温度×圧力) の粗グリッド総当たり',
                'grid_vars': ' × '.join(grid_vars), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 出典レポートと手法",
        f"- DME 2025 (成熟度 {MATURITY}、検出 {PROBLEMS})。多段反応器の 段数×入口温度×圧力 を総当たりグリッド。",
        "- PDH 写像: 段数概念なし→反応器の 入口温度×触媒層厚 の粗グリッドに置換。",
        "- 体現する欠陥: P05 粗い離散総当たり / P06 大域性未保証 / P02 / P12。",
    ]
    out_dir, best = harness.finalize(study, method='rep_dme2025_mat4',
        p_codes='DME2025 成熟度4/12 (C02グリッド→P05/P06)',
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
