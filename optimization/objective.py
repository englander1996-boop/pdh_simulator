"""
optimization/objective.py — Optuna 用 objective 関数ファクトリ

main.py から呼ぶ make_objective() が、設計変数の suggest → FlowsheetDesignVars 構築
→ evaluate() 呼び出し → effective_TAC 返却 までを 1 ステップにラップした関数を返す。

診断情報 (failure_reason, TAC, revenue, spec 値など) は trial.user_attrs に格納し、
後段の reporting で CSV/JSON に出力できるようにする。
"""

from typing import Callable, Dict, Any, Optional

from flowsheet import evaluate, FlowsheetResult
from config.load import OperatingConfig

from optimization.search_space import (
    suggest_params, build_design, extract_F_fresh_override, VarSpec,
)
from optimization.penalty_scale import set_scale, default_schedule


def make_objective(
    search_space:      Dict[str, VarSpec],
    solver_assignment: Dict[str, str],
    config:            OperatingConfig,
    *,
    apply_hi:               bool  = True,
    apply_stage2:           bool  = False,
    hi_dT_min_K:            float = 10.0,
    strict_recovery_check:  bool  = False,
    recovery_tolerance:     float = 0.10,
    baseline:               Dict[str, Any] | None = None,
    n_trials_total:         int = 300,
    penalty_schedule:       Optional[Callable[[int, int], float]] = None,
) -> Callable:
    """Optuna study.optimize() に渡す objective 関数を生成して返す。

    Parameters
    ----------
    search_space : dict
        main.py の SEARCH_SPACE 辞書。
    solver_assignment : dict
        {'dist1': 'fug'|'rigorous'|'sm', ...}。BO ループでは SOLVER_BO を渡す。
    config : OperatingConfig
        load_operating_config() で読み込んだ運転条件。
    apply_hi, apply_stage2, hi_dT_min_K, strict_recovery_check, recovery_tolerance
        evaluate() にそのまま渡す引数。
    baseline : dict | None
        SEARCH_SPACE で suggest されない変数のデフォルト値 (None → DEFAULT_BASELINE)。

    Returns
    -------
    objective : callable
        Optuna study.optimize() に渡す `objective(trial) -> float`。
    """
    eval_kwargs = dict(
        apply_hi=apply_hi,
        apply_stage2=apply_stage2,
        hi_dT_min_K=hi_dT_min_K,
        strict_recovery_check=strict_recovery_check,
        recovery_tolerance=recovery_tolerance,
        verbose=False,
    )
    # adaptive penalty schedule (default は 3 段階 step、custom も渡せる)
    _schedule = penalty_schedule if penalty_schedule is not None else default_schedule

    def objective(trial) -> float:
        # 設計判断 (2026-05-21): trial 開始時に penalty scale を更新。
        # 序盤 (0-30%): scale=0.3 で探索広げる
        # 中盤 (30-70%): scale=1.0 標準
        # 終盤 (70-100%): scale=3.0 で infeas 強制退出
        scale = _schedule(trial.number, n_trials_total)
        set_scale(scale)
        trial.set_user_attr('penalty_scale', scale)

        params = suggest_params(trial, search_space)
        design = build_design(params, solver_assignment, baseline=baseline)
        F_fresh_override = extract_F_fresh_override(params, baseline=baseline)

        # evaluate は内部で solve_flowsheet の例外を catch し、solver_failure_okuyen を
        # 返す実装 (runner.py:118-125 参照)。ここで raise されることは想定しない。
        result: FlowsheetResult = evaluate(
            design, config,
            F_C3H8_override=F_fresh_override,
            **eval_kwargs,
        )

        _store_diagnostics(trial, result)
        if F_fresh_override is not None:
            trial.set_user_attr('F_C3H8_fresh_used_kmol_h', F_fresh_override)
        return result.effective_TAC

    return objective


def _store_diagnostics(trial, result: FlowsheetResult) -> None:
    """Trial.user_attrs に診断情報を格納。Reporting で CSV カラムとして抽出できる。"""
    trial.set_user_attr('failure_reason', result.failure_reason)
    trial.set_user_attr('is_feasible', result.is_feasible)

    if result.economics is not None:
        trial.set_user_attr('TAC_okuyen',     result.economics.TAC)
        trial.set_user_attr('revenue_okuyen', result.economics.total_revenue)
        trial.set_user_attr('profit_raw_okuyen',
                            result.economics.total_revenue - result.economics.TAC)
    if result.economics_hi is not None:
        trial.set_user_attr('TAC_hi_okuyen',     result.economics_hi.TAC)
        trial.set_user_attr('revenue_hi_okuyen', result.economics_hi.total_revenue)
        trial.set_user_attr('profit_hi_okuyen',
                            result.economics_hi.total_revenue - result.economics_hi.TAC)
    if result.economics_synth is not None:
        trial.set_user_attr('TAC_stage2_okuyen',     result.economics_synth.TAC)
        trial.set_user_attr('revenue_stage2_okuyen', result.economics_synth.total_revenue)
        trial.set_user_attr('profit_stage2_okuyen',
                            result.economics_synth.total_revenue - result.economics_synth.TAC)

    if result.specs is not None:
        trial.set_user_attr('c3h6_purity_wtfrac', result.specs.c3h6_purity_wtfrac)
        trial.set_user_attr('h2_purity_molfrac',  result.specs.h2_purity_molfrac)
        trial.set_user_attr('production_kmol_h',  result.specs.production_kmol_h)
        trial.set_user_attr('target_kmol_h',      result.specs.target_kmol_h)

    # rigorous プロキシ罰則の内訳 (Phase A デバッグ用)
    if result.solver is not None and result.solver.one_pass:
        total_proxy = 0.0
        for col_key in ('r1', 'r2', 'r3'):
            col_r = result.solver.one_pass.get(col_key)
            if col_r is None or col_r.equipment is None:
                continue
            p = getattr(col_r.equipment, 'proxy_penalty_okuyen', 0.0)
            if p > 0:
                trial.set_user_attr(f'proxy_penalty_{col_key}_okuyen', p)
                total_proxy += p
        if total_proxy > 0:
            trial.set_user_attr('proxy_penalty_total_okuyen', total_proxy)

        # 設計判断 (2026-05-20): 蒸留塔 infeasibility の連続シグナルを user_attrs に格納。
        # run_one_pass の _build_penalty_after_column が dist{1,2,3}_{N,dT}_shortfall を
        # 一括計算済み。constraints_func (study.py) でこれを TPE に渡す。
        # 正常完走時は 0.0、infeasible 時は連続値 (FUG: N 不足比、Rigorous: log(dT/tol))。
        op = result.solver.one_pass
        for col_idx in ('1', '2', '3'):
            for kind in ('N', 'dT'):
                key = f'dist{col_idx}_{kind}_shortfall'
                v = op.get(key, 0.0) or 0.0
                if v > 0:
                    trial.set_user_attr(key, float(v))
