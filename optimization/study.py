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
    def __init__(self, phase1_sampler, phase2_sampler, switch_at_n_trials: int,
                 constraints_func: Optional[Callable] = None,
                 n_constraints: int = 0):
        self._s1 = phase1_sampler
        self._s2 = phase2_sampler
        self._switch_at = switch_at_n_trials
        # 設計判断 (2026-05-21): QMC phase で constraint が記録されないため、TPE が
        # phase 切替時に全 trial を「constraint 欠落」として警告する spam が
        # study.py:51 から大量発生していた。after_trial で phase1 完了 trial に対し
        # constraints_func を実行して system_attrs にダミー constraint を埋め込む。
        self._constraints_func = constraints_func
        self._n_constraints    = n_constraints

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
        # 設計判断 (2026-05-21): QMC phase でも constraints を system_attrs に格納する。
        # TPE は trial.system_attrs[_CONSTRAINTS_KEY] を見て constraint が無い trial を
        # warning + lower priority 扱いするため、QMC trial に対しても同じ constraints_func
        # を実行して保存する必要がある。Optuna は内部的に '_constraints' を使用する。
        # active sampler の after_trial を呼ぶ前に self が constraints を格納すれば、
        # 次回以降の sampling で警告が出なくなる。
        if (self._constraints_func is not None
                and state == optuna.trial.TrialState.COMPLETE):
            try:
                cons = self._constraints_func(trial)
                # Optuna ≥3.0 で TPE が読むキー名 (内部 API)。
                # constraint 値が存在することだけ示せれば良い。
                # 設計判断 (2026-05-21): Optuna 内部の _CONSTRAINTS_KEY = "constraints"
                # (optuna/samplers/_base.py)。これを system_attrs に格納すると TPE の
                # constraint 欠落警告が抑制される。
                study._storage.set_trial_system_attr(
                    trial._trial_id, 'constraints', list(cons),
                )
            except Exception:
                pass
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
      [8] psa_t_abs_shortfall        : PSA 吸着時間 < _T_ABS_MIN の log10 比 (2026-05-21 追加)
      [9] psa_u_0_shortfall          : PSA 空塔速度 > _U0_MAX の log10 比 (2026-05-21 追加)
      [10] psa_feed_shortfall        : PSA feed 異常 (no_non_C3/no_CH4) で 1.0 (2026-05-21 追加)
      [11] reactor_sv_shortfall      : Reactor SV 範囲外への log10 比 (2026-05-21 追加)
      [12] reactor_other_shortfall   : Reactor その他 penalty (sim_failure 等) で 1.0 (2026-05-21 追加)
      [13] production_under_pp       : 生産量下限不足 [%pt] (2026-05-21 追加、F_fresh ↑ シグナル)
      [14] production_over_pp        : 生産量上限超過 [%pt] (2026-05-21 追加、F_fresh ↓ シグナル)
      [15] mem_ph_shortfall          : Mem P_H ≤ feed.P_in の log10 比 (2026-05-22 追加、P_H↑/Dist2P↓)
      [16] mem_bp_shortfall          : Mem 透過 bp ≤ 冷却水出口の log10 比 (2026-05-22 追加、P_dist3↑)
      [17] mem_phase_shortfall       : Mem ガス feed のはずが T_in < 露点 (2026-05-22 追加、Dist2 P↓)
      [18] mem_other_shortfall       : Mem その他 penalty (ODE/comp/cond 失敗等) で 1.0 (2026-05-22 追加)
      [19] trace_bypass_psa_excess   : PSA trace bypass 閾値超過分 [fraction] (2026-05-22 改良2 追加、
                                       TPE が borderline 滞留から離れる方向シグナル: Dist2 厚化/分離強化)
      [20] trace_bypass_mem_excess   : Mem trace bypass 閾値超過分 [fraction] (2026-05-22 改良2 追加、
                                       同上、Mem 経路の非 C3 混入率)
      [21] unknown_failure           : is_feasible=False かつ全 shortfall=0 で 1.0 (silent infeas 保険)
      [22] dist2_cond_shortfall      : HYSYS Dist2 cold-top (凝縮器ΔT不成立) の塔頂温度不足 [K]×0.1 (2026-05-28 追加、col2_p↑/reflux↓ 方向シグナル)

    設計判断 (2026-05-20): 旧版は [proxy, feas_flag] の 2 制約のみで、ValueError
    (Dist1 FUG 全ゼロ)・Wang-Henke 失敗 (Dist2) 等の異種 infeasibility が同じ
    binary 信号に潰されていた → TPE が「方向」を学習できなかった。本版は塔別に
    N 不足 / dT 不足の **連続値** を追加し、infeasible 領域内でも勾配が立つように
    する (例: Dist1 N=16 で不足比 0.06、N=10 で 0.5 → TPE は 0.06 を相対的に優先)。

    設計判断 (2026-05-21): PSA silent _penalty_result() 経路 (t_abs<MIN, u_0>MAX 等)
    を連続値化。main_20260521_131507 で 300/300 trial が PSA silent penalty で stuck し、
    BO が「L/D を上げれば feasible に出る」を学習できなかった問題への対処。

    TPE は constraint 違反 trial を「達成不可能と判断するための情報」として使う。
    n_startup_trials 後の TPE モデルに非線形な選好を入れられる。
    """
    # 設計判断 (2026-05-22 forensic, 施策 1b): constraint 値域を「marginal violation で 1.0」に正規化。
    # 旧版は raw 値をそのまま使っており、scale 差が大きいと TPE 内のランキングが歪む:
    #   - production_pp: 0-15 (大、median violation ~3.6)
    #   - mem_bp_shortfall: 0-0.008 (小、median ~0.007)
    #   - proxy: 0-200 [億円] (極大、median ~30)
    #   - reactor_sv: 0-0.36 (小)
    # raw 値で TPE が constraint 違反量を sum すると proxy + prod_pp が他を圧倒し、
    # 個別装置の方向シグナル (mem_bp 等) が埋もれる。
    # 各 scale を「典型的な marginal 違反」が ~1.0 に揃うよう係数を掛ける。
    # 既に O(1) スケールのもの (d_N, log10 系) はスルー。
    proxy_raw = trial.user_attrs.get('proxy_penalty_total_okuyen', 0.0)
    proxy = proxy_raw * 0.1                                                    # 10 億円 = 1.0
    is_feasible = trial.user_attrs.get('is_feasible', True)
    feas_violation = 0.0 if is_feasible else 1.0
    d1_N  = trial.user_attrs.get('dist1_N_shortfall', 0.0)                     # 0-1, 既に O(1)
    d2_N  = trial.user_attrs.get('dist2_N_shortfall', 0.0)
    d3_N  = trial.user_attrs.get('dist3_N_shortfall', 0.0)
    d1_dT = trial.user_attrs.get('dist1_dT_shortfall', 0.0)                    # log10 ratio, 既に O(1)
    d2_dT = trial.user_attrs.get('dist2_dT_shortfall', 0.0)
    d3_dT = trial.user_attrs.get('dist3_dT_shortfall', 0.0)
    # 設計判断 (2026-05-28): HYSYS Dist2 cold-top (凝縮器ΔT不成立、r2 失敗の86%) の連続シグナル。
    # 塔頂が凝縮可能下限を何K下回るか [K]。10K = 1.0 に正規化 (典型失敗は ~2-15K)。
    # 旧来この失敗は N/dT shortfall が 0 で unknown_failure 扱い=方向勾配なし → feasible 頭打ちの主因。
    d2_cond_raw = trial.user_attrs.get('dist2_cond_shortfall', 0.0)
    d2_cond = d2_cond_raw * 0.1
    psa_t = trial.user_attrs.get('psa_t_abs_shortfall', 0.0)                   # log10, 既に O(1)
    psa_u = trial.user_attrs.get('psa_u_0_shortfall', 0.0)
    psa_f = trial.user_attrs.get('psa_feed_shortfall', 0.0)                    # binary
    rx_sv_raw = trial.user_attrs.get('reactor_sv_shortfall', 0.0)
    # 設計判断 (2026-05-23 forensic, main_20260523_172800): 倍率 5.0 → 20.0 に拡大。
    # 根拠: 172800 run の r_rx fail trial (16件) の reactor_sv_shortfall raw 値は 0.02-0.13、
    #   ×5 後でも 0.1-0.6 で、他の constraint 信号 (prod_pp 1.0, proxy 1.0, d_N O(1)) と
    #   比べて弱信号。これが「D_reactor を bounds で絞ったら r_rx が増えた」(165334→172800
    #   で r_rx 11→16) の構造的説明: 旧 bounds では TPE が「D=9.5-10.0 = SV 安全帯」を
    #   弱信号でも学習できていたが、新 bounds で安全帯を切られると弱信号では D=7.5 への
    #   張付きを振り払えなかった。bounds で切るのではなく constraint で学習させるという
    #   撤退判断 (5/23) と整合させるため、SV 信号自体を強化する。
    # ×20 後の典型値: 0.4-2.6 で他信号と同水準。期待: r_rx 11 件 → 5-7 件 (172800 撤退後ベース)。
    rx_sv = rx_sv_raw * 20.0                                                   # log10 ratio が小さいため強拡大
    rx_ot = trial.user_attrs.get('reactor_other_shortfall', 0.0)               # binary
    prod_under_raw = trial.user_attrs.get('production_under_pp', 0.0)
    prod_under = prod_under_raw * 0.2                                          # 5 pp = 1.0
    prod_over_raw  = trial.user_attrs.get('production_over_pp', 0.0)
    prod_over  = prod_over_raw * 0.2
    mem_ph    = trial.user_attrs.get('mem_ph_shortfall',    0.0)               # log10 ratio, O(1)
    mem_bp_raw = trial.user_attrs.get('mem_bp_shortfall',    0.0)
    mem_bp    = mem_bp_raw * 50.0                                              # 0.02 = 1.0 (T_bp gap ~5K 相当)
    mem_phase = trial.user_attrs.get('mem_phase_shortfall', 0.0)               # log10, O(1)
    mem_other = trial.user_attrs.get('mem_other_shortfall', 0.0)               # binary
    tb_psa_raw = trial.user_attrs.get('trace_bypass_psa_excess', 0.0)
    tb_psa    = tb_psa_raw * 100.0                                             # 1% 超過 = 1.0
    tb_mem_raw = trial.user_attrs.get('trace_bypass_mem_excess', 0.0)
    tb_mem    = tb_mem_raw * 100.0
    # 設計判断 (2026-05-22 forensic, 施策 1a): silent constraint plug。
    # is_feasible=False かつ全 shortfall=0 (= 完全 silent infeas) なる trial は、
    # 旧版だと constraints_func が全要素 0 を返し、TPE は「constraint 上は feas」と
    # 誤判定し得る。is_feasible=False 単独だけが伝わる状態は学習素材として弱い。
    # 214750 run forensic で 50/130=38% が silent → unknown_failure_penalty=1.0 で
    # 「ここは原因不明だが infeas」を TPE に最低限伝える。
    # 観測強化コード (2026-05-22) で大半の silent は failure_unit + 個別 shortfall で
    # 解消するが、未知 / timeout / exception 経路の保険として残す。
    # raw 値で判定 (normalize 後の値だと小さい raw 値が「silent」誤判定するため)
    raw_total = (
        proxy_raw + d1_N + d2_N + d3_N + d1_dT + d2_dT + d3_dT +
        psa_t + psa_u + psa_f + rx_sv_raw + rx_ot +
        prod_under_raw + prod_over_raw +
        mem_ph + mem_bp_raw + mem_phase + mem_other +
        tb_psa_raw + tb_mem_raw + d2_cond_raw
    )
    unknown_failure = 1.0 if (not is_feasible and raw_total == 0.0) else 0.0
    # [22] dist2_cond_shortfall : HYSYS Dist2 塔頂が凝縮可能下限を下回る量 (2026-05-28 追加)
    return [proxy, feas_violation, d1_N, d2_N, d3_N, d2_dT, d1_dT, d3_dT,
            psa_t, psa_u, psa_f, rx_sv, rx_ot, prod_under, prod_over,
            mem_ph, mem_bp, mem_phase, mem_other, tb_psa, tb_mem,
            unknown_failure, d2_cond]


def make_sampler(
    name:              str,
    seed:              int,
    n_startup:         int,
    *,
    constraints_func:  Optional[Callable] = None,
    constant_liar:     bool = False,
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
        # 設計判断 (2026-05-23 forensic, main_20260523_190747): multivariate=True と
        # n_ei_candidates=200 を追加 (施策 P + Q)。
        # 背景: 190747 run の 35 trial 時点で TPE phase が「A_mem=1e5 帯 + N_d3=117 + R_d3=13.85」
        # の局所最適に張り付き、QMC で見つかった best #18 (A_mem=2.84e5, N_d3=107, R_d3=10.94,
        # TAC=1007.35) を捨てて TAC=1068-1100 帯に集中。
        # feasible 9件中 #18 と #3 が「A_mem 大 ↔ N_d3 小 ↔ R_d3 小」の負相関を示しており、
        # 物理的に「Mem 大→C3 損失少→Dist3 薄くて良い」という設計の自由度がある。
        # 旧 TPE は各変数を独立 KDE で扱うため、この変数間相関を学習できなかった。
        #
        # 施策 P (multivariate=True): Parzen estimator を多変量化、A_mem×N_d3×R_d3 の
        # 相関構造を学習させる。学習素材 (feasible) 9 件は最低ライン、TPE が #18 系統を
        # 「高 A_mem パターン」として認識できる土台を作る。
        # 副作用: 計算やや重い (1 trial 中 ~10ms 増)、過学習リスク (feasible 少時)。
        #
        # 施策 Q (n_ei_candidates=200): 各 trial で 200 候補から最良 EI を選ぶ (デフォルト 24)。
        # 24 だと候補が直近成功領域に集中しがち、exploitation 偏重を緩和。
        # 副作用: best 収束ペースが緩む可能性 (exploration 強化のトレードオフ)。
        # 設計判断 (2026-05-26): constant_liar は並列(マルチプロセス)最適化用。
        # 実行中(pending)trial を悲観的な仮値で埋めて扱い、複数 worker が同じ有望領域に
        # 群がる(冗長サンプリング)のを防ぐ。共有 SQLite storage 経由で worker 同士の
        # running trial が見えるので、プロセス並列でも機能する。単一プロセスでは False。
        tpe = optuna.samplers.TPESampler(
            seed=seed,
            n_startup_trials=n_startup,
            constraints_func=cf,
            multivariate=True,                  # 施策 P (変数間相関学習)
            n_ei_candidates=200,                # 施策 Q (exploration 強化、デフォルト 24)
            constant_liar=constant_liar,        # 並列時 True (pending trial を悲観評価)
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
                    # 設計判断 (2026-05-21): QMC trial にも constraints を埋め込み、
                    # TPE 切替後の「Trial X does not have constraint values」spam を抑制。
                    constraints_func=cf,
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
