r"""
comparing/case_p06_multistart/main.py — P06「大域最適性の未保証 (単一始点)」の忠実再現 (multi-start)。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p06_multistart\main.py`

再現する問題のあるやり方:
  case_p01_subsystem/case_p02_pinch/case_p04_sequential/case_p05_grid の素朴手法は、単一の始点から最適化して打ち切り、
  得られた点が大域最適である保証 (Multi-start / Branch-and-Bound / 凸緩和 等) を一切確認しない。
  プロセス目的関数は非凸 (装置コスト ∝ size^0.6、相変化の不連続) なので、単一始点は局所解に陥る。

このスクリプトが示すこと (定量化のねらい):
  同じ素朴手法 (座標降下) を **複数の始点** から回し、各始点が辿り着く「best」のばらつきを出す。
  - 始点ごとに別の局所解へ収束する事実 = 単一始点の非大域性 (P06)。
  - そのばらつき幅 (range/std) = 単一始点の素朴手法が「答え」と誤認している局所解の散らばり。
  (収束/2 巡目の未検証 = P12 は case_p12_converge で別途扱う)

出力 (comparing/results/multistart_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式、全始点の全 trial)
  multistart.csv : 始点ごとの best (start_id, kind, best_TAC, feasible, 主要 params)

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
BO との比較はユーザ側: best.json と outputs/special_*/best.json を突合し ΔTAC を算出。
"""

import os
import sys
import math
import random
import statistics

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
N_STARTS = 4               # 始点数 (start0 = 範囲中央、以降は乱数始点)
VARS = ['col2_reflux_ratio', 'col3_n_stages', 'T_in_K']  # 各始点で座標降下する変数 (コストドライバ)
K_POINTS = 3               # 1 変数あたり掃引点数
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P06 大域性未保証 の忠実再現。case_p04_sequential と同じ素朴な座標降下を複数の始点から回し、
各始点が辿り着く best のばらつきを出す。始点ごとに別の局所解へ収束する事実が、
単一始点の素朴手法が大域最適を保証していないこと (P06) を定量的に示す。multistart.csv に記録。
"""


def _random_params(rng: random.Random) -> dict:
    """SEARCH_SPACE の範囲内で 1 つの始点 (21 変数) を一様サンプリング (log 変数は対数一様、int は丸め)。"""
    p = {}
    for name, (low, high, scale, typ) in space.SEARCH_SPACE.items():
        if scale == 'log':
            v = math.exp(rng.uniform(math.log(low), math.log(high)))
        else:
            v = rng.uniform(low, high)
        p[name] = int(round(v)) if typ == 'int' else float(v)
    return p


def _coordinate_descent(study, objective, cb, start_params, vars_, k):
    """start_params から vars_ を 1 変数ずつ掃引→最小採用 (1 巡)。この始点で走った trial 列を返す。"""
    current = {name: space.clamp(name, start_params[name]) for name in space.PARAM_NAMES}
    start_trials = []
    for var in vars_:
        pts = space.grid_points(var, k)
        points = [dict(current, **{var: val}) for val in pts]
        trials = harness.run_batch(study, objective, points, cb)
        start_trials.extend(trials)
        b = harness.best_of(trials)
        if b is not None and var in b.params:
            current[var] = space.clamp(var, b.params[var])
    return current, start_trials


def run(*, n_starts=N_STARTS, vars_=None, k=K_POINTS, backend=None,
        seed=SEED, top_n=TOP_N):
    vars_ = vars_ or VARS
    backend = backend or BACKEND
    rng = random.Random(seed)

    per_start_evals = sum(len(space.grid_points(v, k)) for v in vars_)
    n_total = n_starts * per_start_evals

    print("=" * 72, flush=True)
    print(f"  case_p06_multistart multistart (P06): {n_starts} 始点 × 座標降下({len(vars_)}変数×~{k}点) "
          f"= ~{n_total} 評価", flush=True)
    print(f"  座標降下変数: {' → '.join(vars_)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    multistart_rows = []
    for s in range(n_starts):
        if s == 0:
            start_params = space.midpoint_params()
            kind = 'midpoint'
        else:
            start_params = _random_params(rng)
            kind = 'random'
        print(f"\n==== 始点 {s+1}/{n_starts} ({kind}) から座標降下 ====", flush=True)
        _converged, start_trials = _coordinate_descent(study, objective, cb, start_params, vars_, k)
        b = harness.best_of(start_trials)
        row = {'start_id': s, 'kind': kind,
               'best_TAC': (b.value if b is not None else None),
               'feasible': (b.user_attrs.get('is_feasible') if b is not None else None),
               'best_trial_number': (b.number if b is not None else None)}
        for v in vars_:
            row[f'best_{v}'] = (b.params.get(v) if b is not None else None)
        multistart_rows.append(row)
        if b is not None:
            print(f"    → 始点{s+1} の best: TAC={b.value:.2f} "
                  f"(feasible={b.user_attrs.get('is_feasible')})", flush=True)

    # 始点ごとの best のばらつき = 単一始点が見落とす局所性。
    vals = [r['best_TAC'] for r in multistart_rows if r['best_TAC'] is not None]
    spread_lines = ["## 始点間のばらつき (= 単一始点が見落とす局所性 P06)"]
    if len(vals) >= 1:
        vmin, vmax = min(vals), max(vals)
        std = statistics.pstdev(vals) if len(vals) >= 2 else 0.0
        spread_lines += [
            f"- 始点数 = {len(multistart_rows)} (うち best 取得 {len(vals)})",
            f"- best TAC: 最小 {vmin:.2f} / 最大 {vmax:.2f} / 幅 {vmax - vmin:.2f} 億円/年 / std {std:.2f}",
            "- 単一始点の素朴手法は、この幅の中の 1 点を「答え」と誤認している (大域性の保証なし)。",
            "- 始点ごとに別の局所解へ収束する = 非大域性 (P06)。MINLP/B&B/multi-start が必要。",
        ]

    settings = {
        'n_starts': n_starts,
        'coordinate_descent_vars': ' → '.join(vars_),
        'k_points': k,
        'n_total_evals': n_total,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'seed': seed,
    }
    out_dir, best = harness.finalize(
        study, method='multistart',
        p_codes='P06 大域性未保証 (multi-start)',
        description=_DESCRIPTION, settings=settings, extra_lines=spread_lines,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        multistart_rows, os.path.join(out_dir, 'multistart.csv'),
        ['start_id', 'kind', 'best_TAC', 'feasible', 'best_trial_number'] +
        [f'best_{v}' for v in vars_],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'multistart.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
