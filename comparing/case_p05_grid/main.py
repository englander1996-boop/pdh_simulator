r"""
comparing/case_p05_grid/main.py — P05「整数の連続扱い/粗い掃引」+ P06「大域性未保証」の忠実再現。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p05_grid\main.py`

再現する問題のあるやり方 (最適化手法分析 C02。例: 2019 のレポートは PA 転化率 × POH 転化率を
3×3 の粗いグリッドで総当たり):
  1. コストドライバを 2〜3 変数選び、各変数を粗い格子に切る (他は妥当な固定値=範囲中央)。
  2. その直積 (全組合せ) を総当たり評価し、最小コスト点を採用。
  3. 1 枚のグリッドで打ち切る (細分化・別領域の確認はしない)。

なぜ問題か (定量化のねらい):
  - P05: 整数の段数 (col3_n_stages) を「粗い格子で掃引して最小」と扱い、MINLP の整数性を解かない。
  - P06: 単一の粗いグリッドは大域最適を保証しない。BO (special.py) が別領域で見つけた点との
    ΔTAC が「粗いグリッド総当たりの損失」。

出力 (comparing/results/grid_search_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式)
  grid.csv : 各格子点の (掃引変数の値, effective_TAC, feasible, purity, production)

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
BO との比較はユーザ側: best.json と outputs/special_*/best.json を突合し ΔTAC を算出。
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
# 設定 (編集可)
# ===========================================================================
# 総当たりするグリッド変数。過去手法に倣い 2〜3 変数に絞る (組合せ爆発 + HYSYS が重いため)。
# 既定: Dist3 段数 (整数=P05・CAPEX ドライバ) × Dist2 還流比 (OPEX ドライバ、N と相互作用)。
GRID_VARS = [
    'col3_n_stages',      # 整数段数 (P05 の核)
    'col2_reflux_ratio',  # 還流比 (OPEX)
]
K_POINTS = 4               # 各変数の格子点数 (粗い)。直積 = K^len(GRID_VARS)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P05 整数の連続扱い + P06 大域性未保証 の忠実再現。コストドライバ 2〜3 変数を粗い格子に切り、
その直積を総当たり評価して最小を採用する (1 枚のグリッドで打ち切り)。整数段数を粗く掃引するだけで
MINLP として解かない欠陥、単一グリッドが大域最適を保証しない欠陥を、BO (special.py) best との
ΔTAC で定量化する。grid.csv に各格子点の effective_TAC を記録 (谷の位置・刻みの粗さを可視化)。
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
    print(f"  case_p05_grid grid_search (P05/P06): {' × '.join(f'{v}[{len(a)}]' for v, a in zip(grid_vars, axes))} "
          f"= {n_total} 点 総当たり", flush=True)
    print(f"  backend: {backend}", flush=True)
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
        'grid_vars': ' × '.join(grid_vars),
        'k_points': k,
        'n_total_evals': n_total,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'start': 'midpoint' if not start else 'custom',
    }
    extra = [
        "## 手法固有の出力",
        "- `grid.csv` … 各格子点 (掃引変数の値 → effective_TAC)。",
        "  P05: 整数段数の格子刻みが粗く、真の整数最適を跨ぎうることを確認できる。",
        "  P06: 単一グリッドの谷が大域最適とは限らない (BO が別領域で更新する)。",
    ]
    out_dir, best = harness.finalize(
        study, method='grid_search',
        p_codes='P05 整数の連続扱い / P06 大域性未保証',
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
