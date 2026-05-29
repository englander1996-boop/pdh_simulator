r"""
comparing/case_rep_toluene2023/main.py — 実在レポート「トルエン脱アルキル化 (2023, C03 FUG)」の手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_toluene2023\main.py`
  (出典著者名は repo に残さない方針。対応は Claude メモリのみ)

出典 (テーマ: トルエンの脱アルキル化→ベンゼン、2023、C03 短絡計算FUG、成熟度 3/12、
  検出 P01/P02/P03/P05/P06/P07/P08/P11/P12):
  そのレポートが実際にやった手法 (§3.5 反応器の最適化、原文より):
  最適化変数を **パージ率 [0.025,0.05,0.075,0.1] × 単通反応率 [0.7-0.9] × 操作圧力 [20-40bar]** とし、
  各変数を一定刻みで取り、**80 通りの全組合せを「絨毯爆撃的」にフルグリッド** 評価して最小コスト点を採用。
  蒸留塔は Fenske + McCabe-Thiele (FUG 短絡) で N=35 / 還流比 0.60 を算出。

手法はテーマ非依存。PDH sim にはパージ率・反応器操作圧力の設計変数が無いので、「反応器の多変数を
粗い格子でフルグリッド総当たり」という C02/C03 的構造を **反応器の転化率レバー 3 つ
(入口温度 × 反応器径 × 触媒層厚) のフルグリッド** に写す (パージ率は写せないため省略、と明記)。

体現する欠陥: P05 離散の粗い掃引 / P06 大域性未保証 (単一グリッド) / P01 部分最適化 (反応器単独) /
P02 / P12。 (原典は塔に FUG=P10 も使うが、本 sim は FUG 禁止のため塔は通常 backend で評価)。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P01(→case_p01_subsystem) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
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


# ===========================================================================
# 設定 — C03 「反応器の多変数フルグリッド絨毯爆撃」をそのまま写す
# ===========================================================================
# 原典は パージ率×単通転化率×操作圧力 の 80 点。PDH ではパージ率/反応器圧は未露出のため、
# 転化率を支配する反応器レバー 3 つに写す (フルグリッド構造は維持)。
GRID_VARS = [
    'T_in_K',        # 入口温度 (単通転化率レバー)
    'D_reactor_m',   # 反応器径
    'z_cat_m',       # 触媒層厚=滞留時間
]
K_POINTS = 4                  # 各 4 点 → 4^3 = 64 点 (原典 80 点規模の「絨毯爆撃」)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND

_DESCRIPTION = """\
実在レポート トルエン脱アルキル化2023 (C03 FUG) の手法の忠実再現。原典は反応器の パージ率×単通転化率×
操作圧力 を 80 点フルグリッドで絨毯爆撃。PDH ではパージ率/反応器圧が未露出のため、反応器の転化率レバー
3 つ (入口温度×反応器径×触媒層厚) のフルグリッドに写す。多変数フルグリッドが大域最適を保証しない欠陥
(P05/P06)・反応器単独最適 (P01) を BO best との ΔTAC で定量化。
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
    print(f"  case_rep_toluene2023 (実在レポート再現, C03): 反応器 "
          f"{' × '.join(f'{v}[{len(a)}]' for v, a in zip(grid_vars, axes))} = {n_total} 点 フルグリッド絨毯爆撃", flush=True)
    print(f"  backend: {backend} / 原典: パージ率×転化率×圧力 80点 (パージ率は sim 未露出→省略)", flush=True)
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
        'source_report': 'トルエン脱アルキル化2023 (C03 FUG, 成熟度3/12)',
        'method': '反応器の多変数フルグリッド絨毯爆撃 (原典80点、パージ率は省略)',
        'grid_vars': ' × '.join(grid_vars),
        'k_points': k, 'n_total_evals': n_total, 'backend': backend,
        'note_unmapped': 'パージ率・反応器操作圧力は sim 未露出のため転化率レバーに置換',
    }
    extra = [
        "## 出典レポートと手法",
        "- トルエン脱アルキル化2023 (C03): 反応器の パージ率×単通転化率×操作圧力 を 80 点フルグリッド絨毯爆撃。",
        "  塔は Fenske+McCabe-Thiele (FUG) で N=35/RR0.60。",
        "- PDH へのマッピング: パージ率/反応器圧は未露出→反応器の転化率レバー3つ (入口温度×径×触媒層厚) のフルグリッド。",
        "- 体現する欠陥: P05 粗い離散掃引 / P06 大域性未保証 / P01 反応器単独最適 / P02 / P12。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC。grid.csv でグリッドの粗さ・谷を確認。",
    ]
    out_dir, best = harness.finalize(
        study, method='toluene2023',
        p_codes='トルエン2023 (C03フルグリッド→P01/P05/P06)',
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
