r"""
comparing/case_rep_pdh2025a/main.py — 実在レポート「PDH 2025 (PDH)」の最適化手法の忠実再現。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_rep_pdh2025a\main.py`

出典レポート (最適化手法分析 C04+C01、成熟度 4/12、検出 P03/P05/P06/P07/P08/P09/P11/P12):
  「8 PDH_プロパンの脱水素化によるプロピレンの製造プロセス」2025。

そのレポートが実際にやった手法 (§4.3.4 C3 ガス分離塔の最適化、原文より):
  C3 ガス分離塔 (= 本 sim の Dist3) について、最適化因子に **操作圧力・段数・入口温度** を選び、
  **1 つずつ順に** 掃引して最小コスト点を採用 (他は固定、座標降下・1 巡)。
  フィード段は各条件で還流比が最小になる段を選定。報告された最適点:
    操作圧力 17 bar / 理論段数 60 / 入口温度 15℃ / 還流比 1.64 / フィード段 8 / 塔高 52.0m / 塔径 2.7m。

手法はテーマ非依存なので PDH sim にそのまま写せる。sim の Dist3 設計変数へのマッピング:
  操作圧力 → col3_p_kpa, 段数 → col3_n_stages, (入口温度は sim に直接の knob 無し→
  フィード段比 col3_feed_ratio を「還流比最小化のための段選定」として代替掃引)。

この手法が体現する欠陥 (= 報告の検出 P と一致): P04 1次元逐次 / P05 整数の連続扱い /
P06 大域性未保証 (単一始点1巡) / P12 検証なし。

定量化: 本手法 best と BO (special.py、全変数同時最適化) best の ΔTAC。さらに sim 上の結果を
レポート報告値 (17bar/60段/RR1.64) と突合できる (PDH 同テーマなので数値比較が可能)。

## 含まれる欠陥部品 (この再現が束ねる P。単体実演は対応する case_p##)
- sim で再現・定量化 (◎): P04(→case_p04_sequential) / P05(→case_p05_grid) / P06(→case_p06_multistart) / P12(→case_p12_converge)
- 検出されたが sim 非対応 (△): P03 / P07 / P08 / P09 / P11
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
# 設定 — PDH2025a の手法をそのまま写す
# ===========================================================================
# レポートの掃引順: 操作圧力 → 段数 → 入口温度。sim では Dist3 の対応変数へ写す。
VARS_ORDER = [
    'col3_p_kpa',       # 操作圧力 (報告最適 17 bar = 1700 kPa)
    'col3_n_stages',    # 理論段数 (報告最適 60 段)
    'col3_feed_ratio',  # フィード段比 (報告: 還流比最小の段を選定 → feed 8/60)
]
K_POINTS = 5            # 1 変数あたり掃引点数 (レポート同様の粗い 1 次元掃引)
SEED = 42
TOP_N = 3
BACKEND = space.DEFAULT_BACKEND  # Dist1=SM / Dist2=HYSYS / Dist3=SM

# レポート報告値 (sim 結果との突合用、参考)。
REPORTED = {'col3_p_kpa': 1700.0, 'col3_n_stages': 60, 'reflux': 1.64, 'feed_stage': 8,
            'height_m': 52.0, 'diameter_m': 2.7}

_DESCRIPTION = """\
実在レポート PDH2025a (PDH) の最適化手法の忠実再現。C3 ガス分離塔 (sim の Dist3) を
操作圧力→段数→入口温度(=フィード段比) の順に 1 つずつ掃引→最小採用 (座標降下・単一始点・1 巡)。
P04/P05/P06/P12 を体現。BO best との ΔTAC、および報告最適 (17bar/60段/RR1.64) との突合で
1 次元逐次の損失を定量化する。
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
    print(f"  case_rep_pdh2025a (実在レポート再現, PDH): C3分離塔(Dist3) を {' → '.join(vars_order)} の順に逐次最適化", flush=True)
    print(f"  ~{n_total} 評価 / backend: {backend}", flush=True)
    print(f"  報告最適(参考): 17bar/60段/RR1.64/feed8", flush=True)
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
        'source_report': 'PDH 2025 (PDH), C04+C01, 成熟度4/12',
        'method': 'C3分離塔(Dist3) を 操作圧力→段数→フィード段 の逐次1次元 (単一始点1巡)',
        'vars_order': ' → '.join(vars_order),
        'k_points': k, 'n_total_evals': n_total, 'backend': backend,
        'reported_optimum': REPORTED,
    }
    extra = [
        "## 出典レポートと手法",
        "- PDH2025a (PDH): C3 分離塔を 操作圧力→段数→入口温度 の順に 1 次元逐次最適化。",
        "- 報告最適: 操作圧力17bar / 段数60 / 還流比1.64 / フィード段8 / 塔高52.0m / 塔径2.7m。",
        "- 体現する欠陥: P04 1次元逐次 / P05 整数の連続扱い / P06 大域性未保証 / P12 検証なし。",
        "## 定量化",
        "本手法 best.json と BO (outputs/special_*) を突合し ΔTAC。報告最適値とも sim 上で突合可。",
    ]
    out_dir, best = harness.finalize(
        study, method='pdh2025a',
        p_codes='PDH2025a (C01逐次→P04/P05/P06/P12)',
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
