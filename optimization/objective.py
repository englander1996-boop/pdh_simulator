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
    # 設計判断 (2026-05-22 L1 観測強化): failure_unit (categorical) を保存。
    # callbacks.py が live 表示・tally に使う。CSV groupby 用にも便利。
    # 値の定義は flowsheet/runner.py の FlowsheetResult 内コメント参照。
    fu = getattr(result, 'failure_unit', '') or ''
    if fu:
        trial.set_user_attr('failure_unit', fu)

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
        # 設計判断 (2026-05-21): production violation の方向を TPE に伝える。
        # 旧版は production_violation_pp 単一値で under/over を区別できず、TPE が
        # F_fresh を上下どちらに動かせば良いか学べなかった (main_20260521_160951
        # でソフト fail 6 件中 3 件 over・3 件 under の混在で TPE 混乱)。
        # 別シグナルで提供することで TPE が「under は F_fresh ↑」「over は F_fresh ↓」
        # を独立に学習可能になる。
        trial.set_user_attr('production_direction', result.specs.production_direction)
        if result.specs.production_under_pp > 0:
            trial.set_user_attr('production_under_pp', float(result.specs.production_under_pp))
        if result.specs.production_over_pp > 0:
            trial.set_user_attr('production_over_pp', float(result.specs.production_over_pp))

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
        # 'cond' = HYSYS Dist2 cold-top (凝縮器ΔT不成立) の連続シグナル (2026-05-28 追加)。
        for col_idx in ('1', '2', '3'):
            for kind in ('N', 'dT', 'cond'):
                key = f'dist{col_idx}_{kind}_shortfall'
                v = op.get(key, 0.0) or 0.0
                if v > 0:
                    trial.set_user_attr(key, float(v))

        # 設計判断 (2026-05-21): PSA silent penalty 経路の連続シグナル。
        # psa_system.py で _T_ABS_MIN/_U0_MAX 等を理由ラベル付きで返す改修と対応。
        # run_one_pass._compute_psa_shortfall が以下 3 キーを一括計算済み。
        for psa_key in ('psa_t_abs_shortfall', 'psa_u_0_shortfall', 'psa_feed_shortfall'):
            v = op.get(psa_key, 0.0) or 0.0
            if v > 0:
                trial.set_user_attr(psa_key, float(v))

        # 設計判断 (2026-05-21): Reactor SV silent penalty も同パターンで連続シグナル化。
        # swing.py の SV 範囲外チェックで penalty_reason='sv_out_of_range' を埋め、
        # run_one_pass._compute_reactor_shortfall が log10 比で reactor_sv_shortfall を生成。
        for rx_key in ('reactor_sv_shortfall', 'reactor_other_shortfall'):
            v = op.get(rx_key, 0.0) or 0.0
            if v > 0:
                trial.set_user_attr(rx_key, float(v))

        # 設計判断 (2026-05-22): Mem silent penalty 経路を連続シグナル化。
        # membrane_system.py の 16 個の _penalty_result() に penalty_reason ラベルを付け、
        # run_one_pass._compute_mem_shortfall が ph/bp/phase/other の 4 種に分類。
        # 旧 BO ログ (main_20260522_005631) では 240/300 trial が無方向で死んでおり、
        # うち 176 件は shortfall attr が一切無い完全 silent だった。本シグナルで
        # 「P_H 不足」「Dist3 圧力不足」「Dist2 圧力過大 (露点高)」を独立に学習可能化。
        for mem_key in ('mem_ph_shortfall', 'mem_bp_shortfall',
                        'mem_phase_shortfall', 'mem_other_shortfall'):
            v = op.get(mem_key, 0.0) or 0.0
            if v > 0:
                trial.set_user_attr(mem_key, float(v))

        # 設計判断 (2026-05-22 改良 2): trace_bypass の excess を連続シグナル化。
        # main_20260522_094436 で TPE が trace bypass borderline (TAC=1028 等) に
        # 16 trial 中 8 件はまり込んだ → user_attr に出してなくて binary 「is_feasible=False」
        # としてしか TPE が認識できなかった。connect to constraints_func で TPE が
        # 「あとどれだけ漏れを減らせば良いか」の連続勾配を持てるように。
        # 単位: 閾値 (=1%) 超過分の fraction (例: 0.013 = 1.3pp 超過 = 0.3pp over)
        for tb_key in ('trace_bypass_psa_excess', 'trace_bypass_mem_excess'):
            v = op.get(tb_key, 0.0) or 0.0
            if v > 0:
                trial.set_user_attr(tb_key, float(v))

        # 設計判断 (2026-05-22 L1 観測強化): 装置別の penalty_reason (categorical 文字列) と
        # key actual 値を user_attr に保存。failure_unit (= "どの装置で詰まったか") に対し、
        # こちらは「その装置内のどのラベルか」「実値はいくつか」のサブ情報。
        # 例: failure_unit='r_mem' に対し mem_penalty_reason='bp_le_cold_out',
        # mem_T_bp_perm_actual_K=305.2, mem_T_cold_out_actual_K=313.0。
        # callbacks.py が「Mem.bp_le_cold_out (T_bp=305<313)」のような live 表示に使う。
        # run_one_pass の _extract_unit_diagnostics が one_pass dict に書き込み済み。
        # 0 / '' は default 値なので保存スキップ (CSV を肥大化させない)。
        # str 系
        for str_key in ('first_failed_unit',
                        'reactor_penalty_reason', 'psa_penalty_reason', 'mem_penalty_reason',
                        'r1_penalty_msg', 'r2_penalty_msg', 'r3_penalty_msg'):
            v = op.get(str_key, '') or ''
            if v:
                trial.set_user_attr(str_key, str(v))
        # float 系
        for num_key in (
            'reactor_SV_actual_m_s',
            'psa_t_abs_actual_s', 'psa_u_0_actual_m_s',
            'mem_P_H_actual_Pa', 'mem_P_feed_actual_Pa',
            'mem_T_dew_actual_K', 'mem_T_feed_actual_K',
            'mem_T_bp_perm_actual_K', 'mem_T_cold_out_actual_K',
            'r1_N_needed', 'r1_dT_max_K',
            'r2_N_needed', 'r2_dT_max_K',
            'r3_N_needed', 'r3_dT_max_K',
        ):
            v = op.get(num_key, 0.0) or 0.0
            if v > 0:
                trial.set_user_attr(num_key, float(v))

    # 設計判断 (2026-05-22 L1 観測強化): warning ソース別カウント。
    # run_one_pass が _capture_warnings 経由で全 warning を集約済み (silent fallback
    # 検出用)。1 trial の warning 総数 + source 別カウントを保存し、「Dist2 で
    # 5 warning が出てる trial が増えている」のようなパターンを CSV で追える形に。
    wc = result.warnings_captured or []
    if wc:
        trial.set_user_attr('warnings_count_total', len(wc))
        # source 別カウント (例: "Dist1=2 PSA=1 Mem=3" の形式文字列で保存)
        from collections import Counter
        src_counter: Counter = Counter()
        for w in wc:
            src = getattr(w, 'source', 'unknown') or 'unknown'
            src_counter[src] += 1
        if src_counter:
            trial.set_user_attr(
                'warnings_count_by_source',
                ' '.join(f"{k}={v}" for k, v in src_counter.most_common()),
            )
