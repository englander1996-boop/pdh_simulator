r"""
comparing/case_rep_pdh2024a/main.py — 実在レポート「PDH 2024 (PDH)」の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_pdh2024a\main.py`

出典レポート (最適化手法分析 C04+C01、成熟度 4/12、検出 P03/P04/P05/P06/P07/P08/P09/P12):
  「5 PDH_プロパンの脱水素反応によるプロピレン製造プロセス」2024。

そのレポートが実際にやった手法 (§3.3 最適化、原文より):
  最適化因子に **反応器の単通反応率 XA** を選び、温度範囲 460–540℃・L/D 比 2–4 の制約下で
  **反応器の長さと半径を変えて** XA を 1 次元的に変化させ、総コスト最小点を探した。
  報告された最適: XA=0.261 (総コスト最小)、反応器 高さ9.0m/直径3.0m/L/D=2.4、6 基直列。
  (転化率↑→反応器体積↑だが過剰原料コスト↓ で総コストが下に凸)。

手法はテーマ非依存。PDH sim では反応器の単通転化率を支配する設計変数へ写す:
  入口温度 T_in_K (転化率の主レバー) → 反応器径 D_reactor_m → 触媒層厚 z_cat_m の順に逐次掃引。
  (本レポートは反応器形状で XA を動かしたので、形状変数を含めるのが忠実)。

体現する欠陥: P04 1次元逐次 (反応器を分離と切り離して単独最適) / P05 / P06 / P12。
定量化: 本手法 best と BO best の ΔTAC。XA や反応器寸法も sim 上で報告値と突合できる。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P03 / P07 / P08 / P09
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


# ===========================================================================
# 設定 — PDH2024a の手法をそのまま写す (反応器の単通転化率を1次元最適化)
# ===========================================================================
VARS_ORDER = [
    'T_in_K',        # 入口温度 — 単通転化率 XA の主レバー (本レポート: 460-540℃で XA を変化)
    'D_reactor_m',   # 反応器径 — 本レポートは反応器形状(長さ/半径)で XA を変化
    'z_cat_m',       # 触媒層厚 — 同上 (滞留時間=転化率)
]
K_POINTS = 5
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND

REPORTED = {'XA_single_pass': 0.261, 'reactor_height_m': 9.0, 'reactor_D_m': 3.0,
            'LD_ratio': 2.4, 'n_reactors_series': 6}

_DESCRIPTION = """\
実在レポート PDH2024a (PDH) の最適化手法の忠実再現。反応器の単通転化率 XA を最適化因子とし、
反応器形状/温度を 1 次元的に動かして総コスト最小点を探す (分離工程とは切り離した単独最適、単一始点1巡)。
P04/P05/P06/P12 を体現。BO best との ΔTAC、および報告最適 (XA=0.261) との突合で定量化する。
"""


def run(*, vars_order=None, k=K_POINTS, start=None, backend=None,
        seed=SEED, top_n=TOP_N):
    vars_order = vars_order or VARS_ORDER
    backend = backend or BACKEND

    current = space.midpoint_params()
    if start:
        current.update(start)
    current = {name: space.clamp(name, current[name]) for name in space.PARAM_NAMES}

    n_total = sum(len(space.grid_points(v, k)) for v in vars_order)
    print("=" * 72, flush=True)
    print(f"  case_rep_pdh2024a (実在レポート再現, PDH): 反応器を {' → '.join(vars_order)} の順に逐次最適化", flush=True)
    print(f"  ~{n_total} 評価 / backend: {backend} / 報告最適(参考): 単通転化率 XA=0.261", flush=True)
    print("=" * 72, flush=True)

    objective = harness.make_objective(backend=backend)
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    cost_curves = []
    for step, var in enumerate(vars_order, 1):
        pts = space.grid_points(var, k)
        points = [dict(current, **{var: val}) for val in pts]
        print(f"\n--- step {step}/{len(vars_order)}: '{var}' を {len(pts)} 点掃引 (他固定) ---", flush=True)
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
            print(f"    → '{var}' = {current[var]} を採用 (TAC={b.value:.2f})", flush=True)

    settings = {
        'source_report': 'PDH 2024 (PDH), C04+C01, 成熟度4/12',
        'method': '反応器の単通転化率 XA を 反応器形状/温度の 1 次元掃引で最適化 (分離と切り離し)',
        'vars_order': ' → '.join(vars_order),
        'k_points': k, 'n_total_evals': n_total, 'backend': backend,
        'reported_optimum': REPORTED,
    }
    extra = [
        "## 出典レポートと手法",
        "- PDH2024a (PDH): 反応器の単通転化率 XA を最適化因子に、反応器形状/温度を 1 次元掃引。",
        "- 報告最適: XA=0.261 / 反応器 高さ9.0m・直径3.0m・L/D2.4 / 6 基直列。",
        "- 体現する欠陥: P04 1次元逐次 (反応器を分離と分離) / P05 / P06 / P12。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC。XA・反応器寸法も sim 上で突合可。",
    ]
    out_dir, best = harness.finalize(
        study, method='pdh2024a',
        p_codes='PDH2024a (反応器1次元→P04/P05/P06/P12)',
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
