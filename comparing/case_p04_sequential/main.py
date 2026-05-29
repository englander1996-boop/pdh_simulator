r"""
comparing/case_p04_sequential/main.py — P04「1次元逐次最適化」+ P05「整数の連続扱い」の忠実再現。

★この case フォルダで完結する独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p04_sequential\main.py`
  (case ごとにフォルダを分け、それぞれ独立 main を持つ。中央ディスパッチャは無し)

再現する問題のあるやり方 (最適化手法分析 C01。例: ある PDH レポートは操作圧力→段数→入口温度を
1 変数ずつ掃引（他レポートも同様）):
  1. 初期エンジニアリング推算 (= 全変数を妥当な固定値=範囲中央) から出発。
  2. 変数を決めた順序で 1 つずつ取り上げ、他を固定したまま粗い格子で掃引し、
     最小コストの点を採用して固定 → 次の変数へ (= 座標降下 / Gauss-Seidel)。
  3. 1 巡で打ち切る (2 巡目・始点感度は確認しない = P12 の温床。検証は case_p12_converge(2巡目)・case_p06_multistart(始点感度) で扱う)。

なぜ問題か (定量化のねらい):
  - P04: 段数 N と還流比 R のように相互作用する変数を 1 次元ずつ最適化すると、
    変数間トレードオフを捉えられず局所停留する。同時最適化 (BO) との TAC 差 ΔTAC が
    「1 次元逐次の損失」。順序を変えると別の答えになる (順序依存) ことも cost_curves から読める。
  - P05: 整数の段数を「粗い格子で掃引して最小を採る」だけで、MINLP として整数性を扱わない。

出力 (comparing/results/sequential_1d_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (special と同形式)
  cost_curves.csv : step ごとの (変数, 掃引値, effective_TAC, feasible, purity, production)

蒸留塔バックエンドは special.py と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM (FUG は使わない)。
BO との比較はユーザ側: best.json と outputs/special_*/best.json を突合し ΔTAC を算出。
"""

import os
import sys

# このファイルは comparing/case_p04_sequential/ 直下。repo root = 3 つ上。
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
try:
    sys.stdout.reconfigure(encoding='utf-8')  # ライブログの ★✓✗█░Δ→ 用 (cp932 回避、special.py と同じ)
except Exception:
    pass

from comparing.shared import space, simulator, reporting, harness


# ===========================================================================
# 設定 (編集可)
# ===========================================================================
# 掃引する変数と順序 (上流→下流。documented practice に倣いコストドライバ中心)。
VARS_ORDER = [
    'col2_p_kpa',          # Dist2(HYSYS) 操作圧力 — 深冷コンデンサ費の主レバー
    'col2_n_stages',       # Dist2 段数 (整数 = P05)
    'col2_reflux_ratio',   # Dist2 還流比 — N と相互作用 (P04 の核)
    'col3_n_stages',       # Dist3 段数 (整数 = P05、CAPEX ドライバ)
    'col3_p_kpa',          # Dist3 操作圧力
    'F_C3H8_fresh_kmol_h', # 原料流量 — 生産量帯と収率のトレードオフ
    'T_in_K',              # 反応器入口温度 — 転化率
]
K_POINTS = 5               # 1 変数あたりの掃引点数 (粗い = 過去手法に忠実)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # {dist1:sm, dist2:hysys, dist3:sm} — special.py と同じ

_DESCRIPTION = """\
P04 逐次1次元 + P05 整数の連続扱い の忠実再現。全変数を範囲中央から出発し、
VARS_ORDER の順に 1 変数ずつ他を固定して粗い格子で掃引→最小を採用 (座標降下、1 巡)。
変数間相互作用を捉えられず局所停留する欠陥を、BO (special.py) best との ΔTAC で定量化する。
cost_curves.csv に各 1 次元掃引の軌跡 (極小の拾い方・整数掃引の粗さ) を記録。
"""


def run(*, vars_order=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND

    # 初期点: 範囲中央 (= 妥当な初期推算)。start dict で上書き可。
    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    n_total = sum(len(space.grid_points(v, k)) for v in vars_order)

    print("=" * 72, flush=True)
    print(f"  case_p04_sequential sequential_1d (P04/P05): {len(vars_order)} 変数 × ~{k} 点 = ~{n_total} 評価", flush=True)
    print(f"  順序: {' → '.join(vars_order)}", flush=True)
    print(f"  backend: {backend}", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    cost_curves = []
    for step, var in enumerate(vars_order, 1):
        pts = space.grid_points(var, k)
        points = []
        for val in pts:
            p = dict(current)
            p[var] = val
            points.append(p)
        print(f"\n--- step {step}/{len(vars_order)}: '{var}' を {len(pts)} 点掃引 (他は固定) ---", flush=True)
        trials = harness.run_batch(study, objective, points, cb)
        for t, val in zip(trials, pts):
            cost_curves.append({
                'step': step, 'var': var, 'value': val,
                'effective_TAC': t.value,
                'feasible': t.user_attrs.get('is_feasible'),
                'purity_wt': t.user_attrs.get('c3h6_purity_wtfrac'),
                'production_kmol_h': t.user_attrs.get('production_kmol_h'),
                'trial_number': t.number,
            })
        b = harness.best_of(trials)
        if b is not None and var in b.params:
            current[var] = space.clamp(var, b.params[var])
            print(f"    → '{var}' = {current[var]} を採用 (TAC={b.value:.2f}, "
                  f"feasible={b.user_attrs.get('is_feasible')})", flush=True)
        else:
            print(f"    → '{var}' は全点 infeasible/評価不能、固定値 {current[var]} のまま", flush=True)

    settings = {
        'vars_order': ' → '.join(vars_order),
        'k_points': k,
        'n_total_evals': n_total,
        'backend': backend,
        'penalty_scale': harness.PENALTY_SCALE,
        'passes': 1,
        'start': 'midpoint' if not start else 'custom',
    }
    extra = [
        "## 手法固有の出力",
        "- `cost_curves.csv` … 各 step の 1 次元掃引 (変数×掃引値→effective_TAC)。",
        "  P04: 各掃引が他変数固定下の極小しか見ていないこと、",
        "  P05: 整数段数 N の掃引粒度が粗く真の整数最適を跨ぎうることを確認できる。",
    ]
    out_dir, best = harness.finalize(
        study, method='sequential_1d',
        p_codes='P04 1次元逐次 / P05 整数の連続扱い',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        cost_curves, os.path.join(out_dir, 'cost_curves.csv'),
        ['step', 'var', 'value', 'effective_TAC', 'feasible',
         'purity_wt', 'production_kmol_h', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'cost_curves.csv')}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
