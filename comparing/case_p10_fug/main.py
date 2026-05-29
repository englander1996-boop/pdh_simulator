r"""
comparing/case_p10_fug/main.py — P10「非理想系での短絡計算 (FUG) の精度問題」の忠実再現＋可視化。

★独立スクリプト。実行: `.\.venv\Scripts\python.exe comparing\case_p10_fug\main.py`

再現する問題のあるやり方 (最適化手法分析 C03 / 問題点 P10):
  分離塔の設計・最適化に Fenske-Underwood-Gilliland (FUG) 短絡計算を用い、その結果を
  そのまま信用する。FUG は相対揮発度 α 一定を仮定するため、非理想・narrow-margin 設計では
  HYSYS真値 (MESH/HYSYS) と数十%ずれる。学生は FUG が「feasible・安い」と言えばそれを採用する。

**FUG を使うなら、その精度問題を必ず可視化する** (ユーザ指示 2026-05-29):
  同じ設計を ① FUG backend と ② HYSYS真値 backend (special.py と同じ Dist1=SM/Dist2=HYSYS/Dist3=SM)
  の両方で評価し、両者の乖離を出す。FUG が feasible/安いと判定した設計を HYSYS真値 が
  infeasible/高コストと覆す箇所が「FUG の精度問題」(P10) の定量的証拠。

掃引対象: Dist2 (FUG-HYSYS真値 乖離が最も集中する塔。docs/solver_choice_rationale.md より
  proxy_penalty 発火率 Dist2 8.7% vs Dist1/Dist3 ~0%)。col2_reflux_ratio × col2_n_stages を粗グリッド。

出力 (comparing/results/fug_accuracy_<ts>/):
  trials.csv / best.json / top{1..N}_trial*.txt / README.md  (HYSYS真値=真値 ベースで保存)
  fug_vs_HYSYS真値.csv : 各設計の FUG 予測 vs HYSYS真値 真値 (TAC/純度/feasible) と乖離

注: HYSYS真値 側は Dist2=HYSYS を含むため、実走行は HYSYS のある PC で (このマシンは import 確認まで)。
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
from optimization.objective import _store_diagnostics
from optimization.penalty_scale import set_scale


# ===========================================================================
# 設定
# ===========================================================================
FUG_BACKEND = {'dist1': 'fug', 'dist2': 'fug', 'dist3': 'fug'}   # 学生の短絡計算
RIG_BACKEND = space.DEFAULT_BACKEND                              # HYSYS真値 真値 (SM/HYSYS/SM)
# FUG-HYSYS真値 乖離が集中する Dist2 のコストドライバを粗グリッド。
GRID_VARS = ['col2_reflux_ratio', 'col2_n_stages']
K_POINTS = 3
SEED = 42
TOP_N = 3

_DESCRIPTION = """\
P10 非理想系での FUG 精度問題の忠実再現＋可視化。Dist2 のコストドライバを粗グリッドで掃引し、
各設計を FUG と HYSYS真値(SM/HYSYS/SM) の両方で評価。FUG が feasible/安いと判定した設計を
HYSYS真値 が infeasible/高コストと覆す乖離が P10 の証拠。fug_vs_HYSYS真値.csv に FUG vs 真値を記録。
"""


def _extract(r):
    eff = float(r.effective_TAC) if r.effective_TAC is not None else None
    tac = float(r.economics.TAC) if getattr(r, 'economics', None) is not None else None
    pur = None
    prod = None
    if getattr(r, 'specs', None) is not None:
        pur = (float(r.specs.c3h6_purity_wtfrac)
               if r.specs.c3h6_purity_wtfrac is not None else None)
        prod = (float(r.specs.production_kmol_h)
                if r.specs.production_kmol_h is not None else None)
    return dict(eff_TAC=eff, TAC=tac, purity=pur, prod=prod, feasible=bool(r.is_feasible))


def _make_dual_objective():
    """同じ設計を FUG と HYSYS真値 で評価。trial value = HYSYS真値(真値) effective_TAC。FUG は user_attr。"""
    def objective(trial):
        set_scale(harness.PENALTY_SCALE)
        p = harness.suggest_all(trial)
        design_p = {k: p[k] for k in space.PARAM_NAMES}
        F = float(p['F_C3H8_fresh_kmol_h'])
        opts = dict(simulator.EVAL_KWARGS_DEFAULT)

        r_fug = simulator.raw_evaluate(space.build_design(design_p, FUG_BACKEND), F_fresh=F, **opts)
        e_fug = _extract(r_fug)

        r_rig = simulator.raw_evaluate(space.build_design(design_p, RIG_BACKEND), F_fresh=F, **opts)
        _store_diagnostics(trial, r_rig)  # 標準 attr (is_feasible/純度/生産量) は HYSYS真値=真値 から
        for kk, vv in e_fug.items():
            trial.set_user_attr('fug_' + kk, vv)
        trial.set_user_attr('F_C3H8_fresh_used_kmol_h', F)
        return r_rig.effective_TAC
    return objective


def run(*, grid_vars=None, k=K_POINTS, start=None, seed=SEED, top_n=TOP_N):
    grid_vars = grid_vars or GRID_VARS

    base = space.midpoint_params()
    if start:
        base.update(start)
    base = {name: space.clamp(name, base[name]) for name in space.PARAM_NAMES}

    axes = [space.grid_points(v, k) for v in grid_vars]
    combos = list(itertools.product(*axes))
    n_total = len(combos)

    print("=" * 72, flush=True)
    print(f"  case_p10_fug fug_accuracy (P10): Dist2 {' × '.join(grid_vars)} = {n_total} 設計を "
          f"FUG vs HYSYS真値 で二重評価", flush=True)
    print(f"  FUG={FUG_BACKEND} / HYSYS真値={RIG_BACKEND}", flush=True)
    print("=" * 72, flush=True)

    objective = _make_dual_objective()
    study = harness.new_study(seed)
    cb = reporting.make_callback(n_total)

    points = [dict(base, **{v: val for v, val in zip(grid_vars, combo)}) for combo in combos]
    trials = harness.run_batch(study, objective, points, cb)

    rows = []
    n_fug_says_feasible = 0
    n_fug_misleads = 0     # FUG=feasible だが HYSYS真値=infeasible
    max_tac_gap = 0.0
    for t, combo in zip(trials, combos):
        a = t.user_attrs
        feas_rig = bool(a.get('is_feasible', False))
        feas_fug = bool(a.get('fug_feasible', False))
        tac_fug = a.get('fug_TAC')
        eff_fug = a.get('fug_eff_TAC')
        pur_fug = a.get('fug_purity')
        pur_rig = a.get('c3h6_purity_wtfrac')
        tac_rig = a.get('economics_TAC') if a.get('economics_TAC') is not None else t.value
        if feas_fug:
            n_fug_says_feasible += 1
            if not feas_rig:
                n_fug_misleads += 1
        if tac_fug is not None and t.value is not None:
            max_tac_gap = max(max_tac_gap, abs(float(t.value) - float(tac_fug)))
        row = {f'var_{v}': val for v, val in zip(grid_vars, combo)}
        row.update({
            'TAC_fug': tac_fug, 'eff_TAC_fug': eff_fug, 'purity_fug': pur_fug, 'feasible_fug': feas_fug,
            'eff_TAC_HYSYS真値': t.value, 'purity_HYSYS真値': pur_rig, 'feasible_HYSYS真値': feas_rig,
            'trial_number': t.number,
        })
        rows.append(row)

    extra = [
        "## 手法固有の出力",
        "- `fug_vs_HYSYS真値.csv` … 各設計の FUG 予測 vs HYSYS真値 真値 (TAC/純度/feasible)。",
        "",
        "## FUG の精度問題 (P10) の定量",
        f"- 評価設計数 = {n_total}",
        f"- FUG が feasible と判定 = {n_fug_says_feasible} 件、うち **HYSYS真値 では infeasible = "
        f"{n_fug_misleads} 件** (FUG の誤判定)。",
        f"- |eff_TAC(HYSYS真値) − TAC(FUG)| の最大乖離 ≈ {max_tac_gap:.2f} 億円/年。",
        "- FUG を信用すると、HYSYS真値 で詰む/高コストな設計を「最適」と誤認する = P10。",
    ]
    settings = {
        'fug_backend': FUG_BACKEND, 'HYSYS真値_backend': RIG_BACKEND,
        'grid_vars': ' × '.join(grid_vars), 'k_points': k,
        'n_designs': n_total, 'n_evals': n_total * 2,
        'note': '各設計を FUG と HYSYS真値 で二重評価。trial value = HYSYS真値(真値)。',
    }
    out_dir, best = harness.finalize(
        study, method='fug_accuracy',
        p_codes='P10 非理想でFUG (精度問題の可視化)',
        description=_DESCRIPTION, settings=settings, extra_lines=extra,
        eval_kwargs=dict(simulator.EVAL_KWARGS_DEFAULT), top_n=top_n,
    )
    harness.save_table_csv(
        rows, os.path.join(out_dir, 'fug_vs_HYSYS真値.csv'),
        [f'var_{v}' for v in grid_vars] +
        ['TAC_fug', 'eff_TAC_fug', 'purity_fug', 'feasible_fug',
         'eff_TAC_HYSYS真値', 'purity_HYSYS真値', 'feasible_HYSYS真値', 'trial_number'],
    )
    print(f"  手法固有: {os.path.join(out_dir, 'fug_vs_HYSYS真値.csv')}", flush=True)
    print(f"  FUG 誤判定 (FUG feasible→HYSYS真値 infeasible): {n_fug_misleads}/{n_fug_says_feasible}", flush=True)
    simulator.shutdown()
    return out_dir, best


def main():
    run()


if __name__ == '__main__':
    main()
