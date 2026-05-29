r"""
comparing/case_rep_eo2025/main.py — 実在レポート「エチレンオキシド製造 2025」の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_eo2025\main.py`
  (著者名は repo に残さない。対応は Claude メモリのみ)

出典 (テーマ: エチレンオキシド(EO)、2025、C04+C01、**成熟度 7/12 = 本コーパス最高**、
  検出 P05/P06/P07/P08/P11 のみ ← P01/P02/P03/P04/P12 は無し = 比較的良い手順):
  そのレポートが実際にやった手法 (§3.4.1 反応器入口温度の最適化、原文より):
  P.I. (=製品利益−原料コスト−塔コスト−用役コスト) を評価関数に、**反応器入口温度** を最適化。
  入口温度↑で単通反応率↑・反応器コスト↓だが、選択率低下で原料コスト↑のトレードオフ。
  併せて EO 吸収塔・蒸留塔段数も最適化 (段数 7/6 段)。後置ピンチではなく熱統合をループ内で扱い、
  検証も行っている (P02/P12 が無い理由)。

手法はテーマ非依存。PDH sim では「反応器入口温度を中心に、塔段数も含めて最適化」に写す:
  T_in_K (入口温度=転化率/選択率トレードオフ) → col3_n_stages (主分離塔段数) を逐次。
  熱統合は apply_hi=True でループ内評価 (後置にしない = 本レポートに P02 が無いことに対応)。

**成熟度が高い (7/12) ので欠陥は少ない**が、残る P05(整数の連続扱い)・P06(大域性)・P07(単目的)・
P08(経済単純化) は BO との ΔTAC で依然定量化できる = 「良い手順でも BO には及ばない」を示す対照例。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P05(→case_p05_grid) / P06(→case_p06_multistart)
- 検出されたが sim 非対応 (△): P07 / P08 / P11
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


# 成熟度 7/12 (最高)。反応器入口温度中心 + 塔段数の逐次最適化。
VARS_ORDER = ['T_in_K', 'col3_n_stages']
K_POINTS = 5
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND
MATURITY = '7/12 (コーパス最高)'
PROBLEMS = 'P05/P06/P07/P08/P11'

_DESCRIPTION = """\
実在レポート EO2025 (成熟度7/12=最高) の手法の忠実再現。反応器入口温度を P.I.(利益)基準で最適化し、
塔段数も含め逐次に決める。熱統合はループ内・検証あり (P02/P12 が無い良い手順)。残る P05/P06/P07/P08 を
BO best との ΔTAC で定量化 = 「良い手順でも BO に及ばない」対照例。
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
    print(f"  case_rep_eo2025 (実在レポート再現, EO, 成熟度{MATURITY}): {' → '.join(vars_order)} を逐次最適化", flush=True)
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

    settings = {'source_report': 'EO レポート 2025 (C04+C01)', 'maturity': MATURITY, 'problems': PROBLEMS,
                'method': '反応器入口温度(利益基準)+塔段数 の逐次最適化、熱統合はループ内',
                'vars_order': ' → '.join(vars_order), 'k_points': k, 'n_total_evals': n_total, 'backend': backend}
    extra = [
        "## 出典レポートと手法",
        f"- EO 2025 (成熟度 {MATURITY}、検出 {PROBLEMS})。反応器入口温度を P.I.(利益)基準で最適化+塔段数。",
        "- P02(後置ピンチ)・P12(検証なし)・P04 が無い = 比較的良い手順 (熱統合ループ内・検証あり)。",
        "- それでも残る P05/P06/P07/P08 を BO との ΔTAC で定量化 (良い手順でも BO 未満の対照例)。",
    ]
    out_dir, best = harness.finalize(study, method='rep_eo2025_mat7',
        p_codes='EO2025 成熟度7/12 (P05/P06/P07/P08)',
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
