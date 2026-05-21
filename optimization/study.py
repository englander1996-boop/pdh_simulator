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


# ---------------------------------------------------------------------------
# 2-Phase sampler: Sobol で startup → TPE に切替
# ---------------------------------------------------------------------------
# 設計判断 (2026-05-21): Optuna 標準の TPESampler は startup を内部 RandomSampler で
# 生成する仕様で、Sobol/Halton 等の低乖離点列に差し替える正規 API がない。
# 完了 trial 数で sampler を切替える薄いラッパ delegator を用意する。
class _PhaseSwitchSampler(optuna.samplers.BaseSampler):
    """完了 trial 数で sampler を切替える 2-phase wrapper。

    Phase 1 (完了数 < switch_at): phase1_sampler (Sobol QMC) で広域カバレッジ
    Phase 2 (それ以降): phase2_sampler (TPE) で集中探索

    constraint_func は phase2_sampler 側に既に組込まれている前提 (TPE 側で受領済)。
    QMC は constraints 非対応だが、phase1 中も failure_reason は記録され、
    phase2 突入時に TPE が全 trial 履歴 (= phase1 結果含む) から学習する。
    """
    def __init__(self, phase1_sampler, phase2_sampler, switch_at_n_trials: int):
        self._s1 = phase1_sampler
        self._s2 = phase2_sampler
        self._switch_at = switch_at_n_trials

    def _active(self, study):
        n_complete = sum(1 for t in study.trials
                         if t.state == optuna.trial.TrialState.COMPLETE)
        return self._s1 if n_complete < self._switch_at else self._s2

    def infer_relative_search_space(self, study, trial):
        return self._active(study).infer_relative_search_space(study, trial)

    def sample_relative(self, study, trial, search_space):
        return self._active(study).sample_relative(study, trial, search_space)

    def sample_independent(self, study, trial, param_name, param_distribution):
        return self._active(study).sample_independent(study, trial, param_name, param_distribution)

    def reseed_rng(self):
        self._s1.reseed_rng()
        self._s2.reseed_rng()

    def before_trial(self, study, trial):
        # Optuna ≥3.0 の hook
        try:
            self._active(study).before_trial(study, trial)
        except AttributeError:
            pass

    def after_trial(self, study, trial, state, values):
        try:
            self._active(study).after_trial(study, trial, state, values)
        except AttributeError:
            pass


def _default_constraints_func(trial: optuna.trial.FrozenTrial) -> Sequence[float]:
    """Optuna TPE 用 constraints_func (Phase C, 2026-05-19 → 2026-05-20 拡張)。

    各要素が**負値なら feasible、正値なら違反**として TPE に解釈される。
    objective.py の _store_diagnostics で trial.user_attrs に設定する値を読む。

    制約 (拡張):
      [0] proxy_penalty_total_okuyen : rigorous プロキシ罰則合計 [億円/年]
      [1] feasibility flag           : is_feasible=False で 1.0 (非収束等)
      [2] dist1_N_shortfall          : Dist1 FUG Gilliland 不足比 (= max(0,(N_needed-N)/N))
      [3] dist2_N_shortfall          : Dist2 FUG Gilliland 不足比 (rigorous 時は通常 0)
      [4] dist3_N_shortfall          : Dist3 FUG Gilliland 不足比
      [5] dist2_dT_shortfall         : Dist2 Wang-Henke 収束不足 (= log10(dT_max/tol))
      [6] dist1_dT_shortfall         : Dist1 rigorous 不足 (FUG 運用時は 0)
      [7] dist3_dT_shortfall         : Dist3 rigorous 不足 (FUG 運用時は 0)

    設計判断 (2026-05-20): 旧版は [proxy, feas_flag] の 2 制約のみで、ValueError
    (Dist1 FUG 全ゼロ)・Wang-Henke 失敗 (Dist2) 等の異種 infeasibility が同じ
    binary 信号に潰されていた → TPE が「方向」を学習できなかった。本版は塔別に
    N 不足 / dT 不足の **連続値** を追加し、infeasible 領域内でも勾配が立つように
    する (例: Dist1 N=16 で不足比 0.06、N=10 で 0.5 → TPE は 0.06 を相対的に優先)。

    TPE は constraint 違反 trial を「達成不可能と判断するための情報」として使う。
    n_startup_trials 後の TPE モデルに非線形な選好を入れられる。
    """
    proxy = trial.user_attrs.get('proxy_penalty_total_okuyen', 0.0)
    is_feasible = trial.user_attrs.get('is_feasible', True)
    feas_violation = 0.0 if is_feasible else 1.0
    d1_N  = trial.user_attrs.get('dist1_N_shortfall', 0.0)
    d2_N  = trial.user_attrs.get('dist2_N_shortfall', 0.0)
    d3_N  = trial.user_attrs.get('dist3_N_shortfall', 0.0)
    d1_dT = trial.user_attrs.get('dist1_dT_shortfall', 0.0)
    d2_dT = trial.user_attrs.get('dist2_dT_shortfall', 0.0)
    d3_dT = trial.user_attrs.get('dist3_dT_shortfall', 0.0)
    return [proxy, feas_violation, d1_N, d2_N, d3_N, d2_dT, d1_dT, d3_dT]


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
        # 設計判断 (2026-05-21): startup を Sobol QMC で置換。
        # 旧 TPESampler は startup の n_startup_trials を内部 RandomSampler で生成するが、
        # 21 次元空間で pure random は coverage 偏り発生 → 狭い feasible 領域を見逃しがち。
        # Sobol 低乖離点列で序盤 100 trial を網羅的サンプリング、TPE が「とにかく feasible
        # を 1 つでも掴む」確率を底上げする。
        tpe = optuna.samplers.TPESampler(
            seed=seed,
            n_startup_trials=n_startup,
            constraints_func=cf,
        )
        if n_startup > 0:
            try:
                qmc = optuna.samplers.QMCSampler(
                    qmc_type='sobol',
                    seed=seed,
                    warn_independent_sampling=False,
                )
                return _PhaseSwitchSampler(
                    phase1_sampler=qmc,
                    phase2_sampler=tpe,
                    switch_at_n_trials=n_startup,
                )
            except Exception as e:
                # QMCSampler が使えない optuna 旧版なら TPE 単体にフォールバック
                import warnings
                warnings.warn(
                    f"QMCSampler 初期化失敗 ({type(e).__name__}: {e})、TPE 単体にフォールバック。",
                    UserWarning, stacklevel=2,
                )
                return tpe
        return tpe
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
    n_jobs:              int    = 1,
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
    # 設計判断 (2026-05-21): n_jobs で並列化。Optuna は ThreadPoolExecutor で
    # objective を並列実行する。SQLite storage は WAL mode でなければロック競合
    # が発生する可能性あり (n_jobs ≤ 4 推奨)。
    # 注意: penalty_scale は process-local global なので、複数スレッドで scale を
    # 同時更新するとレース発生 (各 trial 開始時に set_scale → 同 thread 内の
    # evaluate が get_scale で読む)。Python の GIL 下では概ね安全だが、完全保証は
    # ないため最終評価では n_jobs=1 推奨。
    study.optimize(
        objective,
        n_trials=n_trials,
        show_progress_bar=show_progress_bar,
        catch=catch,
        timeout=timeout_sec,
        callbacks=callbacks,
        n_jobs=n_jobs,
    )
    return study
