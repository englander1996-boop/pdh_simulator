r"""
comparing/case_rep_styrene2025/main.py — 実在レポート「スチレンモノマー製造 2025」の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_styrene2025\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: エチルベンゼン脱水素→スチレンモノマー、2025、C04+C01、**成熟度 2/12 = 欠陥多**、
  検出 P01/P02/P03/P04/P05/P06/P08/P09/P11/P12 ← 10 個):
  そのレポートが実際にやった手法 (§4.4 反応器の最適化、原文より):
  **EB(エチルベンゼン)反応率** を最適化因子に、反応工程のみのコストの和を反応率ごとに求め、
  反応率 0.45 のとき総コスト最小と決定 (1 次元・反応器ブロック単独・単一始点)。
  分離工程・熱統合は別途・後段で設計 (= 部分最適化 P01 + 後置ピンチ P02)。

手法はテーマ非依存。PDH sim では「反応器ブロックだけを反応率(=入口温度)で 1 次元最適化」に写す:
  T_in_K のみを掃引 (反応器単独最適化、下流は固定)。

**成熟度が低い (2/12) ので欠陥が多い** (P01 部分最適化 + P04 1次元 + P02 後置 + P12 検証なし …)。
反応器を単独・1 変数で決める最も極端な部分最適化で、BO との ΔTAC が最大級になることが期待される対照例。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P01(→case_p01_subsystem) / P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P03 / P08 / P09 / P11
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


# 成熟度 2/12 (欠陥多)。反応器ブロックだけを反応率(=入口温度) 1 次元で最適化 (部分最適化)。
VARS_ORDER = ['T_in_K']
K_POINTS = 7
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '2/12 (欠陥多)'
PROBLEMS = 'P01/P02/P03/P04/P05/P06/P08/P09/P11/P12'

_DESCRIPTION = """\
実在レポート スチレン2025 (成熟度2/12=欠陥多) の手法の忠実再現。反応器ブロックだけを反応率(=入口温度)で
1 次元最適化し下流は固定 (部分最適化 P01 + 1次元 P04 + 後置ピンチ P02 + 検証なし P12)。反応器を単独・
1 変数で決める最も極端な手順で、BO best との ΔTAC が最大級になる対照例。
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
    print(f"  case_rep_styrene2025 (実在レポート再現, スチレン, 成熟度{MATURITY}): "
          f"反応器単独を {' → '.join(vars_order)} で 1 次元最適化", flush=True)
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

    settings = {'source_report': 'スチレン レポート 2025 (C04+C01)', 'maturity': MATURITY, 'problems': PROBLEMS,
                'method': '反応器ブロックのみ反応率(=入口温度) 1 次元最適化、下流固定 (部分最適化)',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 出典レポートと手法",
        f"- スチレン 2025 (成熟度 {MATURITY}、検出 {PROBLEMS} = 10 個)。EB反応率を反応工程のみのコストで 1 次元最適化。",
        "- 反応器を単独・1 変数で決める = P01 部分最適化 + P04 1次元 が最も極端な形。下流/熱統合は後段。",
        "- BO との ΔTAC が最大級になる対照例 (低成熟度=多欠陥)。",
    ]
    out_dir, best = harness.finalize(study, method='rep_styrene2025_mat2',
        p_codes='スチレン2025 成熟度2/12 (P01/P04 ほか10)',
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
