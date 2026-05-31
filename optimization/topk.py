"""
optimization/topk.py — BO 上位候補の rigorous + Stage 2 再評価

BO ループは高速化のため全塔 FUG で回すが、FUG は narrow-margin 設計で楽観的に
振れることがある (Dist1/Dist3)。BO 終了後、Profit 上位の候補だけを rigorous +
Stage 2 (HEN synthesis) で再評価し、現実的な最終解を提示する。

設計判断:
  - BO study から effective_TAC 昇順 (≒ Profit 降順) で top-k を抽出。
  - 各候補について、SOLVER_TOPK (例: 全塔 rigorous) + apply_stage2=True で再評価。
  - 再評価で feasibility が崩れたり TAC が大きく変動するケースは、FUG bias が
    効いていた候補とみなして警告フラグ。
"""

from dataclasses import dataclass
from typing import Dict, List, Any, Optional
import optuna

from flowsheet import evaluate, FlowsheetResult
from config.load import OperatingConfig

from optimization.search_space import build_design, extract_F_fresh_override, VarSpec


@dataclass
class TopKEntry:
    """1 候補の BO vs 再評価結果。"""
    rank:              int                # 1 始まり
    trial_number:      int                # Optuna trial 番号
    params:            Dict[str, Any]     # 設計変数値
    effective_TAC_bo:  float              # BO ループでの目的関数値 (= effective_TAC)
    effective_TAC_re:  float              # 再評価後の effective_TAC
    is_feasible_bo:    bool
    is_feasible_re:    bool
    failure_reason_re: str
    result:            FlowsheetResult    # 再評価の完全結果 (経済・spec を含む)


def select_topk(study: optuna.Study, k: int) -> List[optuna.trial.FrozenTrial]:
    """完了 trial を effective_TAC 昇順 (= minimize 目的に合致) で k 件取り出す。

    失敗 trial (state != COMPLETE) は除外。infeasible (ペナルティ込み) でも
    完了していれば対象に含まれる (再評価で改善する可能性があるため)。
    """
    completed = [t for t in study.trials
                 if t.state == optuna.trial.TrialState.COMPLETE
                 and t.value is not None]
    completed.sort(key=lambda t: t.value)
    return completed[:k]


def reevaluate_topk(
    study:             optuna.Study,
    k:                 int,
    solver_assignment: Dict[str, str],
    config:            OperatingConfig,
    *,
    apply_hi:               bool  = True,
    apply_stage2:           bool  = False,   # 2026-05-31: HEN(Stage2)不採用、HIのみ
    hi_dT_min_K:            float = 10.0,
    strict_recovery_check:  bool  = True,
    recovery_tolerance:     float = 0.10,
    baseline:               Dict[str, Any] | None = None,
    verbose:                bool  = False,
) -> List[TopKEntry]:
    """BO 上位 k 候補を再評価し、TopKEntry のリストを返す。

    Parameters
    ----------
    study : optuna.Study
        完了済みの BO study。
    k : int
        再評価する候補数。
    solver_assignment : dict
        {'dist1': 'rigorous', 'dist2': 'rigorous', 'dist3': 'rigorous'} など。
    apply_stage2 : bool
        True (デフォルト) で HEN Synthesis (Stage 2) を実行。
        実 HEN 構成・追加 HE CAPEX が反映される。
    strict_recovery_check : bool
        rigorous の non-spec 解を catch。再評価では True を推奨。
    verbose : bool
        evaluate() の verbose に渡す (デバッグ用)。

    Returns
    -------
    list of TopKEntry
        rank=1 が最良。
    """
    candidates = select_topk(study, k)
    entries: List[TopKEntry] = []

    for rank, trial in enumerate(candidates, start=1):
        design = build_design(trial.params, solver_assignment, baseline=baseline)
        F_fresh_override = extract_F_fresh_override(trial.params, baseline=baseline)
        result = evaluate(
            design, config,
            apply_hi=apply_hi,
            apply_stage2=apply_stage2,
            hi_dT_min_K=hi_dT_min_K,
            strict_recovery_check=strict_recovery_check,
            recovery_tolerance=recovery_tolerance,
            F_C3H8_override=F_fresh_override,
            verbose=verbose,
        )
        is_feasible_bo = bool(trial.user_attrs.get('is_feasible', False))
        entries.append(TopKEntry(
            rank              = rank,
            trial_number      = trial.number,
            params            = dict(trial.params),
            effective_TAC_bo  = float(trial.value),
            effective_TAC_re  = float(result.effective_TAC),
            is_feasible_bo    = is_feasible_bo,
            is_feasible_re    = result.is_feasible,
            failure_reason_re = result.failure_reason,
            result            = result,
        ))

    return entries


def best_entry(entries: List[TopKEntry]) -> Optional[TopKEntry]:
    """再評価後の effective_TAC が最良の候補を返す (feasible 優先)。

    feasible エントリがあればその中で最小、無ければ infeasible 中の最小を返す。
    """
    if not entries:
        return None
    feasibles = [e for e in entries if e.is_feasible_re]
    pool = feasibles if feasibles else entries
    return min(pool, key=lambda e: e.effective_TAC_re)
