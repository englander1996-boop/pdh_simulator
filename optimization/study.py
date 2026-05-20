"""
optimization/study.py — Optuna study の生成と最適化ループの薄いラッパ

main.py から `create_study()` で study を作り、`run_optimization()` を回す。
sampler は名前文字列で切替可能 ('tpe' | 'cmaes' | 'random')、SQLite storage 対応で
中断・再開が可能。

採用理由:
  - TPE: 連続/離散混合に強く、本プロセスのような 19 次元問題で実績が豊富。
  - CMA-ES: 連続のみだが収束が速いケースあり (整数比率が低い時に有利)。
  - Random: ベースライン比較・スタートアップ探索用。
"""

import optuna
from typing import Callable, Optional, Sequence


def _default_constraints_func(trial: optuna.trial.FrozenTrial) -> Sequence[float]:
    """Optuna TPE 用 constraints_func (Phase C, 2026-05-19)。

    各要素が**負値なら feasible、正値なら違反**として TPE に解釈される。
    objective.py の _store_diagnostics で trial.user_attrs に設定する値を読む。

    制約:
      [0] proxy_penalty_total_okuyen : rigorous プロキシ罰則合計 [億円/年]
          (>0 でも採用するが TPE が積極的に避けるように働く)
      [1] feasibility flag           : is_feasible=False で 1.0 (非収束等)
      [2] spec_violation_pp_total    : spec 違反の合計 pp (>0 で違反)

    TPE は constraint 違反 trial を「達成不可能と判断するための情報」として使う。
    n_startup_trials 後の TPE モデルに非線形な選好を入れられる。
    """
    proxy = trial.user_attrs.get('proxy_penalty_total_okuyen', 0.0)
    is_feasible = trial.user_attrs.get('is_feasible', True)
    feas_violation = 0.0 if is_feasible else 1.0
    # spec 違反 pp は failure_reason に書かれているが、現状 user_attr 化していない
    # ので簡易的に「proxy + feas」の 2 制約だけにする。
    return [proxy, feas_violation]


def make_sampler(
    name:              str,
    seed:              int,
    n_startup:         int,
    *,
    constraints_func:  Optional[Callable] = None,
) -> optuna.samplers.BaseSampler:
    """Sampler を名前から生成。

    Parameters
    ----------
    name : str
        'tpe' | 'cmaes' | 'random'
    seed : int
        乱数シード。
    n_startup : int
        TPE/CMA-ES の冒頭ランダム探索試行数。
        TPE は通常 n_trials // 6 程度が目安、本プロジェクトでは 50 を既定。
    constraints_func : callable | None
        TPE 用 constraints_func (Phase C, 2026-05-19)。trial → Sequence[float] で
        負値=feasible、正値=violated。TPE 内部で feasible/violated を分けて
        学習させる。None なら _default_constraints_func を使う。CMAES/Random は
        constraints 非対応なので無視される (warning なし)。
    """
    name_lower = name.lower()
    if name_lower == 'tpe':
        cf = constraints_func if constraints_func is not None else _default_constraints_func
        return optuna.samplers.TPESampler(
            seed=seed,
            n_startup_trials=n_startup,
            constraints_func=cf,
        )
    elif name_lower == 'cmaes':
        return optuna.samplers.CmaEsSampler(
            seed=seed,
            n_startup_trials=n_startup,
        )
    elif name_lower == 'random':
        return optuna.samplers.RandomSampler(seed=seed)
    else:
        raise ValueError(
            f"未知の sampler 名: {name!r} (許容: 'tpe' | 'cmaes' | 'random')"
        )


def create_study(
    study_name:    str,
    sampler_name:  str  = 'tpe',
    seed:          int  = 42,
    n_startup:     int  = 50,
    storage_url:   Optional[str] = None,
    load_if_exists: bool = True,
) -> optuna.Study:
    """Optuna Study を生成。

    storage_url を指定すると SQLite (例: `sqlite:///outputs/optuna_xxx.db`) に
    履歴が保存され、同じ study_name で再度呼ぶと続きから再開できる。
    storage_url=None なら in-memory (プロセス終了で履歴消失)。
    """
    sampler = make_sampler(sampler_name, seed, n_startup)
    return optuna.create_study(
        study_name=study_name,
        sampler=sampler,
        direction='minimize',
        storage=storage_url,
        load_if_exists=load_if_exists,
    )


def run_optimization(
    study:               optuna.Study,
    objective:           Callable,
    n_trials:            int,
    show_progress_bar:   bool   = True,
    catch:               tuple  = (Exception,),
    timeout_sec:         Optional[float] = None,
    callbacks:           Optional[list] = None,
) -> optuna.Study:
    """Optuna 最適化ループを実行。

    Parameters
    ----------
    study : optuna.Study
        create_study() で生成した study。
    objective : callable
        objective(trial) -> float (= effective_TAC)。make_objective() で生成。
    n_trials : int
        実行する trial 数。SQLite storage の場合、既存試行数に追加される。
    show_progress_bar : bool
        tqdm ベースの進捗表示。
    catch : tuple
        objective 内で raise された例外を Optuna が「失敗 trial」として
        処理する型のタプル。デフォルト (Exception,) で予期せぬ全例外を吸収し、
        該当 trial を TrialState.FAIL にして次に進む (= study 全体は止まらない)。
        evaluate() は内部で例外を catch する設計だが、import エラーや
        メモリ不足など想定外の例外も拾えるよう Exception 全捕捉が安全。
    timeout_sec : float | None
        study 全体の総時間制限 [秒]。経過すると n_trials 未消化でも停止。
        None なら無制限。BO+top-k で長時間走らせる際の安全弁。
    """
    study.optimize(
        objective,
        n_trials=n_trials,
        show_progress_bar=show_progress_bar,
        catch=catch,
        timeout=timeout_sec,
        callbacks=callbacks,
    )
    return study
