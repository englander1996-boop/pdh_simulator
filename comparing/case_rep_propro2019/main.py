r"""
comparing/case_rep_propro2019/main.py — 実在レポート「n-プロピルプロピオネート製造 (2019, C02)」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_propro2019\main.py`
  (出典著者名は repo に残さない方針。対応は Claude メモリのみ)

出典 (テーマ: n-プロピルプロピオネート反応蒸留、2019、C02 多次元グリッド、成熟度 3/12、
  検出 P02/P03/P05/P06/P07/P08/P09/P11/P12):
  そのレポートが実際にやった手法 (§3.8 最適化、原文より):
  PA 気泡塔の反応率を 0.7〜0.95 まで 0.05 刻みで 1 次元掃引 → コスト低い 3 点に絞る →
  各々について POH 気泡塔の反応率を同刻みで掃引し、**反応率 3×3 = 9 通りの粗いグリッド** で総コスト最小を採用。
  報告された最適: PA 反応率 0.9 / POH 反応率 0.85。

手法はテーマ非依存。PDH sim には反応が 1 つしかないので、「2 つの反応率を粗グリッドで総当たり」を
**転化率を支配する反応器レバー 2 つ (入口温度 × 触媒層厚) の 3×3 粗グリッド** に写す
(= 「反応条件を粗い格子で総当たり」という C02 の手法構造をそのまま再現)。

体現する欠陥: P05 整数/離散の連続扱い・粗い刻み / P06 大域性未保証 (単一グリッド) /
P02 熱統合後置 / P12 検証なし。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P02 / P03 / P07 / P08 / P09 / P11
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


# ===========================================================================
# 設定 — C02 「反応条件の粗い 3×3 グリッド総当たり」をそのまま写す
# ===========================================================================
# 原典は PA/POH の 2 反応率を 3×3。PDH は 1 反応なので転化率レバー 2 つに写す。
GRID_VARS = [
    'T_in_K',     # 入口温度 (転化率の主レバー、原典の反応率1に対応)
    'z_cat_m',    # 触媒層厚=滞留時間 (転化率の第2レバー、原典の反応率2に対応)
]
K_POINTS = 3                  # 粗い 3 点 (原典の 3×3=9 を再現)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND

_DESCRIPTION = """\
実在レポート ProPro2019 (C02) の手法の忠実再現。原典は PA/POH の 2 反応率を 0.05 刻みの 3×3=9
グリッドで総当たり (報告最適 PA0.9/POH0.85)。PDH は 1 反応なので転化率レバー 2 つ (入口温度×触媒層厚)
の 3×3 粗グリッドに写す。粗い格子総当たりが大域最適を保証しない欠陥 (P05/P06) を BO best との ΔTAC で定量化。
"""


def run(*, grid_vars=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
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
    print(f"  case_rep_propro2019 (実在レポート再現, C02): "
          f"{' × '.join(f'{v}[{len(a)}]' for v, a in zip(grid_vars, axes))} = {n_total} 点 粗グリッド総当たり", flush=True)
    print(f"  backend: {backend} / 原典: 反応率 3×3=9 (報告最適 PA0.9/POH0.85)", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    points = []
    for combo in combos:
        p = dict(base)
        for v, val in zip(grid_vars, combo):
            p[v] = val
        points.append(p)
    trials = harness.run_batch(study, objective, points, cb)

    grid_rows = []
    for t, combo in zip(trials, combos):
        row = {f'var_{v}': val for v, val in zip(grid_vars, combo)}
        row.update({
            'effective_TAC': t.value,
            'feasible': t.user_attrs.get('is_feasible'),
            'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
            'production_kmol_h': t.user_attrs.get('production_kmol_h'),
            'trial_number': t.number,
        })
        grid_rows.append(row)

    settings = {
        'source_report': 'ProPro2019 (n-プロピルプロピオネート, C02, 成熟度3/12)',
        'method': '反応条件の粗い 3×3 グリッド総当たり (原典は PA/POH 反応率)',
        'grid_vars': ' × '.join(grid_vars),
        'k_points': k, 'n_total_evals': n_total, 'backend': backend,
        'reported_optimum': 'PA反応率0.9 / POH反応率0.85',
    }
    extra = [
        "## 出典レポートと手法",
        "- ProPro2019 (C02): PA→top3→POH の反応率 3×3=9 粗グリッド総当たり。報告最適 PA0.9/POH0.85。",
        "- PDH へのマッピング: 1 反応なので転化率レバー 2 つ (入口温度×触媒層厚) の 3×3 に写す。",
        "- 体現する欠陥: P05 粗い離散掃引 / P06 大域性未保証 / P02 / P12。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC。grid.csv で谷の位置・刻みの粗さを確認。",
    ]
    out_dir, best = harness.finalize(
        study, method='propro2019',
        p_codes='ProPro2019 (C02 粗グリッド→P05/P06)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        grid_rows, os.path.join(out_dir, 'grid.csv'),
        [f'var_{v}' for v in grid_vars] +
        ['effective_TAC', 'feasible', 'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'grid.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
