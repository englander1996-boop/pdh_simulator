r"""
comparing/case_p01_subsystem/main.py — P01「部分最適化 (サブシステム別最適化)」の忠実再現。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p01_subsystem\main.py`

再現する問題のあるやり方:
  プロセスをサブシステム (反応器 → 分離(PSA/膜) → 蒸留) に分け、上流から順に
  「そのブロックの変数だけ」を最適化して確定し、次のブロックへ進む。各ブロックは
  下流を考慮せず (=固定したまま) 最適化され、確定後は二度と見直さない。

なぜ問題か (定量化のねらい):
  - P01: ブロック境界をまたぐ相互作用 (例: 反応器の転化率・出口組成が蒸留の段数/還流比に効く)
    を取りこぼす。各ブロックの局所最適を貼り合わせても全体最適にならず、結合後 TAC が
    同時最適 (BO) より悪い。その ΔTAC が「部分最適化の損失」。

実装: ブロックごとに、そのブロックの変数を 1 変数ずつ粗い格子で掃引→最小採用 (ブロック内座標降下)。
  ブロックの best を固定してから次ブロックへ (上流ブロックは下流確定前に凍結される = P01 の核)。

出力 (comparing/results/subsystem_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式)
  block_curves.csv : (ブロック, 変数, 掃引値, effective_TAC, feasible, ...)

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
# サブシステム (上流→下流) と各ブロックで最適化する変数。HYSYS が重いので各ブロックは
# コストドライバに絞る (全 DOF を入れると評価爆発)。ブロック順 = プロセスの流れ順。
BLOCKS = [
    ('反応器',       ['T_in_K', 'D_reactor_m']),
    ('分離(PSA/膜)', ['desorption_target', 'A_mem_m2']),
    ('蒸留',         ['col2_reflux_ratio', 'col2_n_stages', 'col3_n_stages', 'col3_p_kpa']),
]
K_POINTS = 3               # ブロック内 1 変数あたりの掃引点数 (粗い)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P01 部分最適化 の忠実再現。プロセスを 反応器→分離→蒸留 のサブシステムに分け、上流から順に
そのブロックの変数だけをブロック内座標降下で最適化→確定して次へ進む (下流を見ずに凍結、見直しなし)。
ブロック境界の相互作用を取りこぼし、貼り合わせた設計が同時最適 (BO) に劣る欠陥を ΔTAC で定量化する。
block_curves.csv に各ブロックの掃引軌跡を記録。
"""


def run(*, blocks=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
    blocks = blocks or BLOCKS
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    n_total = sum(len(space.grid_points(v, k)) for _, vs in blocks for v in vs)

    print("=" * 72, flush=True)
    print(f"  case_p01_subsystem subsystem (P01 部分最適化): {len(blocks)} ブロックを上流→下流に順次最適化 "
          f"= ~{n_total} 評価", flush=True)
    print(f"  ブロック: {' → '.join(name for name, _ in blocks)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    block_curves = []
    for bi, (bname, bvars) in enumerate(blocks, 1):
        print(f"\n==== ブロック {bi}/{len(blocks)}: {bname} "
              f"(変数 {', '.join(bvars)}) を最適化 (他ブロック固定) ====", flush=True)
        # ブロック内座標降下: そのブロックの変数を 1 つずつ掃引→最小採用。
        for var in bvars:
            pts = space.grid_points(var, k)
            points = []
            for val in pts:
                p = dict(current)
                p[var] = val
                points.append(p)
            print(f"  -- '{var}' を {len(pts)} 点掃引 --", flush=True)
            trials = harness.run_batch(study, objective, points, cb)
            for t, val in zip(trials, pts):
                block_curves.append({
                    'block': bname, 'var': var, 'value': val,
                    'effective_TAC': t.value,
                    'feasible': t.user_attrs.get('is_feasible'),
                    'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                    'production_kmol_h': t.user_attrs.get('production_kmol_h'),
                    'trial_number': t.number,
                })
            b = harness.best_of(trials)
            if b is not None and var in b.params:
                current[var] = space.clamp(var, b.params[var])
                print(f"     → '{var}' = {current[var]} を採用 (TAC={b.value:.2f}, "
                      f"feasible={b.user_attrs.get('is_feasible')})", flush=True)
            else:
                print(f"     → '{var}' は全点 infeasible/評価不能、固定値 {current[var]} のまま", flush=True)

    settings = {
        'blocks': ' → '.join(f"{name}[{','.join(vs)}]" for name, vs in blocks),
        'k_points': k,
        'n_total_evals': n_total,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'order': 'upstream→downstream, no revisit',
        'start': 'midpoint' if not start else 'custom',
    }
    extra = [
        "## 手法固有の出力",
        "- `block_curves.csv` … 各ブロック内の 1 次元掃引 (ブロック×変数×掃引値→effective_TAC)。",
        "  P01: 上流ブロックが下流を見ずに確定されるため、ブロック境界の相互作用が失われる。",
        "  ブロック順を変えると別の答えになる (= 部分最適化の不安定さ)。",
    ]
    out_dir, best = harness.finalize(
        study, method='subsystem',
        p_codes='P01 部分最適化 (サブシステム別)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        block_curves, os.path.join(out_dir, 'block_curves.csv'),
        ['block', 'var', 'value', 'effective_TAC', 'feasible',
         'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'block_curves.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
