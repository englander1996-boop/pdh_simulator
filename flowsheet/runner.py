"""
最適化器から呼ぶ最上位エントリポイント: `evaluate(design, config) → Result`。

最適化器側の典型的な使い方:
    config = load_operating_config()
    design = FlowsheetDesignVars(swing=..., psa=..., mem=...)
    result = evaluate(design, config, verbose=False)
    objective = result.effective_TAC          # 最小化 (= TAC + ペナルティ)
    feasible  = result.is_feasible            # ペナルティ・spec チェックの集約

2階層ペナルティ構造:
  (a) solver-level 失敗 (PSA/Mem CAPEX 1e8, リサイクル暴走, 未収束):
      数値結果が信頼できないので effective_TAC = penalty.solver_failure_okuyen で
      ハード打ち切り。
  (b) spec 違反 (純度・生産量):
      数値結果は信頼できるので、連続的なソフトペナルティで方向感を与える:
        effective_TAC = TAC + spec_base + spec_coef × Σ(violation%pt)
      全 spec 違反は %pt スケールに正規化済み (flowsheet/specs.py)。

採用理由: 化学プロセス最適化は微分不要系 (DE/GA/BO) が実用的で、それらは
feasibility 領域外でもスカラー値を必要とする。陽な制約より連続ソフトペナルティの
方が探索効率が良い。
"""

from dataclasses import dataclass
from typing import Optional

from flowsheet.design import FlowsheetDesignVars
from flowsheet.solver import solve_flowsheet, SolverResult
from flowsheet.economics import calculate_economics, Economics
from flowsheet.specs import check_specs, SpecComplianceResult
from config.load import OperatingConfig


@dataclass
class FlowsheetResult:
    """フローシート評価結果。

    最適化器は effective_TAC を目的関数として最小化する。
    is_feasible / failure_reason / specs は診断用。

    HI (Heat Integration) を apply_hi=True で評価したときのみ、
    economics_hi と hi_result が設定される。BO で apply_hi=False のときは
    両方 None。
    """
    solver:         SolverResult
    economics:      Optional[Economics]            # solver 失敗時は None
    specs:          Optional[SpecComplianceResult] # solver 失敗時は None
    effective_TAC:  float                          # 最適化器が見る目的関数 [億円/年]
    failure_reason: str                            # "" のとき feasible
    # ---- HI (apply_hi=True のときのみ) ----
    economics_hi:   Optional[Economics] = None     # HI 適用後の Economics
    hi_result:      Optional[object]    = None     # pinch_analysis の HIResult
    # ---- Stage 2 = HEN Synthesis (apply_stage2=True のときのみ、top-k 用) ----
    economics_synth: Optional[Economics] = None    # 実 HEN 構成適用後の Economics
    hen_result:      Optional[object]    = None    # synthesize_hen の HENResult
    # ---- 診断用: 1 パス中に捕捉された warning (run_one_pass の集約) ----
    # silent fallback 検出のため warning を catch_warnings(record=True) で捕捉する。
    # 本フィールドで BO ログから fallback 発火を追跡可能。空リストなら全装置が warning なく動作。
    warnings_captured: list = None                 # _CapturedWarning のリスト
    # ---- 観測ラベル ----
    # failure_reason は人間可読の長文字列だが、live ログ・CSV groupby 用に
    # 「どの段階で詰まったか」を categorical 1 列で持つ。値:
    #   'success'                    : feasible 完走
    #   'r1' / 'r_rx' / 'r2'
    #     / 'r_psa' / 'r_mem' / 'r3' : run_one_pass 内の penalty 経路
    #   'timeout'                    : solver per-trial 時間予算超過
    #   'solver_inner_diverge'       : 内側 (recycle) 未収束
    #   'solver_outer_diverge'       : 外側 (生産量調整) 未収束
    #   'solver_guard'               : リサイクル暴走ガード発火
    #   'strict_recovery_<col>_<lk_or_hk>' : strict recovery check 違反
    #   'spec_production_under/over' : 生産量未達/超過
    #   'spec_c3h6_purity'           : C3H6 純度不足
    #   'spec_h2_purity'             : H2 純度不足
    #   'exception:<UnitName>'       : solve_flowsheet 内で未処理例外 (装置名は best-effort)
    failure_unit: str = ""

    @property
    def is_feasible(self) -> bool:
        return self.failure_reason == ""


#  各塔の (LK, HK, recovery_spec) — column1/2/3.py の wrapper と一致。
#  strict_recovery_check で参照する。
_COLUMN_RECOVERY_SPECS = {
    'r1': {'LK': 'A', 'HK': 'Z', 'rec_LK_top': 0.99, 'rec_HK_bot': 0.99,
           'partial_cond': False, 'name': 'Dist1'},
    'r2': {'LK': 'F', 'HK': 'A', 'rec_LK_top': 0.99, 'rec_HK_bot': 0.99,
           'partial_cond': True,  'name': 'Dist2'},
    'r3': {'LK': 'B', 'HK': 'A', 'rec_LK_top': 0.99, 'rec_HK_bot': 0.99,
           'partial_cond': False, 'name': 'Dist3'},
}


# PSA/Mem trace bypass の閾値超過分を BO objective に伝える連続 penalty 係数
# [億円/年・per fraction-超過]。これがないと BO は「PSA に C3H6 を多量に流す」設計を
# 無自覚に選好し、BO ベスト trial が rigorous 再評価で死ぬ原因になる。
#  1000 億円/fraction → 1pp 超過で 10 億、5pp で 50 億。spec 違反 (1pp で 100 億)
#  より弱く、実コスト (PSA 増床 + C3 売り損 = 10-20 億) に近い水準。
_TRACE_BYPASS_PENALTY_COEF_OKUYEN = 1000.0


# solver-level failure の TAC を連続化する。固定値だと silent infeas 同士の TAC が
# 同値になり、TPE は infeas 内序列を constraint_func の violation sum だけに頼る。
# 本連続化で「shortfall 総和が小さい infeas (=feas に近い可能性)」は低い TAC、
# 「大量 shortfall (=明確に違反)」は高い TAC で順位付けし、TPE は constraint_func と
# TAC の両方を見て feas 方向を学習できる。
#
# マップ: silent_base + raw_total × coef、上限 solver_failure_okuyen
# silent_base = 0.5 × solver_failure で「silent infeas は中間 TAC」を表現
# coef は raw_total ≈ 25 で上限 hit するよう調整 (= shortfall 多いほど早く上限へ)
_SILENT_BASE_FRAC          = 0.5    # silent infeas の TAC = solver_failure × この比率
_SHORTFALL_TAC_COEF_OKUYEN = 200.0  # raw_total 1 単位あたりの TAC 加算

# raw_total に集計する shortfall キー (one_pass dict の key)
_SHORTFALL_KEYS_FOR_TAC = (
    'dist1_N_shortfall', 'dist2_N_shortfall', 'dist3_N_shortfall',
    'dist1_dT_shortfall', 'dist2_dT_shortfall', 'dist3_dT_shortfall',
    'psa_t_abs_shortfall', 'psa_u_0_shortfall', 'psa_feed_shortfall', 'psa_dp_shortfall',
    'reactor_sv_shortfall', 'reactor_dp_shortfall', 'reactor_other_shortfall',
    'mem_ph_shortfall', 'mem_bp_shortfall',
    'mem_phase_shortfall', 'mem_other_shortfall',
    'trace_bypass_psa_excess', 'trace_bypass_mem_excess',
)


def _compute_solver_fail_TAC(solver_failure_ceiling: float,
                              one_pass: Optional[dict]) -> float:
    """solver-level failure の連続化 TAC を計算する。

    one_pass が None (= solve_flowsheet 例外経路) なら ceiling 固定。
    """
    if one_pass is None:
        return float(solver_failure_ceiling)
    raw_total = sum(float(one_pass.get(k, 0.0) or 0.0) for k in _SHORTFALL_KEYS_FOR_TAC)
    silent_base = solver_failure_ceiling * _SILENT_BASE_FRAC
    fail_TAC = silent_base + raw_total * _SHORTFALL_TAC_COEF_OKUYEN
    return min(fail_TAC, float(solver_failure_ceiling))


def evaluate(
    design:                 FlowsheetDesignVars,
    config:                 OperatingConfig,
    verbose:                bool  = False,
    apply_hi:               bool  = True,
    hi_dT_min_K:            float = 10.0,
    apply_stage2:           bool  = False,
    strict_recovery_check:  bool  = False,
    recovery_tolerance:     float = 0.10,
    F_C3H8_override:        Optional[float] = None,
) -> FlowsheetResult:
    """設計変数を入力してフローシートを評価し、effective_TAC と状態を返す。

    Parameters
    ----------
    design, config, verbose : 通常評価用
    apply_hi : bool
        True (デフォルト) のとき pinch targeting を実行し、HI 後 OPEX/TAC/Profit
        を計算する。BO ループでも有効にする想定 (pinch は ms オーダーで重くない)。
        False のときは HI なしで raw OPEX を使用 (デバッグ・比較用)。
    hi_dT_min_K : float
        HI の最小接近温度差 [K]。デフォルト 10K (textbook 標準)。BO の設計変数に
        含めず固定が一般的。
    apply_stage2 : bool
        True のとき、Stage 2 (HEN synthesis: 実マッチング・追加 HE CAPEX 計算)
        を実行する。**top-k 候補の re-evaluation 用**で、BO ループには通常含めない
        (greedy アルゴリズムは smooth でなく、計算もやや重め)。
        apply_hi=True のときのみ意味を持つ (apply_hi=False では何もしない)。
    strict_recovery_check : bool
        True のとき、exp1 outer-loop 収束後に各蒸留塔の **実際の recovery が spec
        ±recovery_tolerance 以内か** を検査し、未達なら failure 扱いにする。
        BO で rigorous 使用時に「数値収束したが物理的に non-spec」な解を捕捉する用。
        recycle iter 中の transient state ではなく、全体収束後の最終状態を検査するので
        過敏発動しない。デフォルト False (= 既存挙動維持、product spec check のみ)。
    recovery_tolerance : float
        strict_recovery_check 時の許容偏差 (デフォルト 0.10 = ±10%)。
        partial_condenser の HK_bot は ≥ spec - tolerance を許容
        (ALWAYS_CONDENSABLE 補正で 100% になるため)。
    """
    pen = config.penalty

    # solve_flowsheet が ValueError 等の例外を投げる場合 (例: Dist1 完全 infeasible で
    # 下流の expansion_valve が流量ゼロ受領)、BO 用途では penalty 返却が望ましい。
    # catch して solver_failure として扱う。
    try:
        solver_result = solve_flowsheet(
            design, config, verbose=verbose,
            F_C3H8_override=F_C3H8_override,
        )
    except Exception as e:
        # 観測ラベル: 例外メッセージから装置名を best-effort で抽出。
        # 既知パターン: "[Dist1] ...", "PSA: ...", "membrane ...", "Comp2/Dist2 ..."
        exc_msg = f"{type(e).__name__}: {e}"
        unit_guess = _guess_unit_from_exception(exc_msg)
        return FlowsheetResult(
            solver=None, economics=None, specs=None,
            # 例外経路は one_pass=None なので連続化は不可、ceiling 固定
            effective_TAC=pen.solver_failure_okuyen,
            failure_reason=f"solve_flowsheet で未処理例外: {exc_msg}",
            failure_unit=f"exception:{unit_guess}",
        )

    # ---- (a) solver-level 失敗 → ハード打ち切り (TAC 連続化) ----
    # 結果の数値が信頼できない場合 (発散・暴走・未収束) は失敗扱い。
    # one_pass の shortfall 総和に応じて silent_base 〜 ceiling の幅で TAC を連続化し、
    # TPE が infeas 内序列を学習しやすくする (詳細は _compute_solver_fail_TAC 参照)。
    # penalty_hit の場合は run_one_pass が埋めた first_failed_unit を優先使用
    # (= 'r_psa' / 'r_mem' / 'r_rx' / 'r1' / 'r2' / 'r3' / 'timeout' を区別)。
    s = solver_result.inner_status
    one_pass_dict = solver_result.one_pass or {}
    solver_fail_TAC = _compute_solver_fail_TAC(pen.solver_failure_okuyen, one_pass_dict)
    if s.penalty_hit:
        unit_from_one_pass = one_pass_dict.get('first_failed_unit', '') or 'solver_penalty'
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=solver_fail_TAC,
            failure_reason="solver-level: PSA/Mem CAPEX ペナルティ発火",
            failure_unit=unit_from_one_pass,
        )
    if s.guard_hit:
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=solver_fail_TAC,
            failure_reason="solver-level: リサイクル暴走ガード発火",
            failure_unit="solver_guard",
        )
    if not s.converged or not solver_result.outer_status.converged:
        if not s.converged and not solver_result.outer_status.converged:
            fu = "solver_inner_outer_diverge"
        elif not s.converged:
            fu = "solver_inner_diverge"
        else:
            fu = "solver_outer_diverge"
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=solver_fail_TAC,
            failure_reason=(
                f"solver-level: 内側{'未' if not s.converged else ''}収束 / "
                f"外側{'未' if not solver_result.outer_status.converged else ''}収束"
            ),
            failure_unit=fu,
        )

    # ---- (a') strict recovery check (BO で rigorous 使用時の追加検査) ----
    # outer-loop 収束後の最終状態で各塔 recovery を確認する。recycle iter の
    # transient state では誤発動するので、全体収束後にだけ検査する。
    # non-spec 解 → solver_failure penalty。
    # feed_LK <= 1e-3 kmol/h (≒ 0 流量) の塔は recovery 検査をスキップするが、
    # silent スキップだと「実は分離できていない塔が pass」する可能性があるため
    # warnings で記録する (BO log の追跡用)。
    import warnings as _warnings
    if strict_recovery_check:
        for col_key, spec in _COLUMN_RECOVERY_SPECS.items():
            col_result = solver_result.one_pass.get(col_key)
            if col_result is None:
                continue
            top, bot = col_result.top.F_in, col_result.bottom.F_in
            # フィードを再構築 (top + bot で復元)
            feed_LK = top.get(spec['LK'], 0.0) + bot.get(spec['LK'], 0.0)
            feed_HK = top.get(spec['HK'], 0.0) + bot.get(spec['HK'], 0.0)
            if feed_LK > 1e-3:
                lk_rec = top.get(spec['LK'], 0.0) / feed_LK
                if abs(lk_rec - spec['rec_LK_top']) > recovery_tolerance:
                    return FlowsheetResult(
                        solver=solver_result, economics=None, specs=None,
                        effective_TAC=solver_fail_TAC,
                        failure_reason=(
                            f"strict recovery check: {spec['name']} LK ({spec['LK']}) "
                            f"recovery={lk_rec:.3f} vs spec {spec['rec_LK_top']:.3f} "
                            f"(差 > {recovery_tolerance*100:.0f}%)"
                        ),
                        failure_unit=f"strict_recovery_{col_key}_lk",
                    )
            else:
                _warnings.warn(
                    f"strict_recovery_check: {spec['name']} LK ({spec['LK']}) "
                    f"feed_LK={feed_LK:.3e} kmol/h ≤ 1e-3 のため検査スキップ。"
                    f"塔への流入がほぼゼロ、上流装置の penalty 状態または微小流量設計の可能性。",
                    UserWarning, stacklevel=2,
                )
            if feed_HK > 1e-3 and not spec['partial_cond']:
                hk_rec = bot.get(spec['HK'], 0.0) / feed_HK
                if abs(hk_rec - spec['rec_HK_bot']) > recovery_tolerance:
                    return FlowsheetResult(
                        solver=solver_result, economics=None, specs=None,
                        effective_TAC=solver_fail_TAC,
                        failure_reason=(
                            f"strict recovery check: {spec['name']} HK ({spec['HK']}) "
                            f"bot recovery={hk_rec:.3f} vs spec {spec['rec_HK_bot']:.3f}"
                        ),
                        failure_unit=f"strict_recovery_{col_key}_hk",
                    )
            elif feed_HK <= 1e-3 and not spec['partial_cond']:
                _warnings.warn(
                    f"strict_recovery_check: {spec['name']} HK ({spec['HK']}) "
                    f"feed_HK={feed_HK:.3e} kmol/h ≤ 1e-3 のため検査スキップ。",
                    UserWarning, stacklevel=2,
                )

    # ---- solver 成功: 経済計算 + spec 判定 ----
    economics = calculate_economics(
        solver_result.one_pass,
        mw_C3H6_kg_per_kmol=config.product.mw_kg_per_kmol,
    )
    specs = check_specs(solver_result.one_pass, config)

    # ---- (b) spec 違反 → ソフトペナルティ ----
    # 全 spec 違反量を %pt スケールに揃えてあるため、係数を1つで管理。
    #   penalty = 違反のあった spec ごとに base を加算 + 全違反 %pt の合計に coef を掛ける
    # 違反 0 のときは加算 0 で本来の TAC のみ。違反があるほど線形に増える勾配を持つ。
    failures = []
    total_violation_pp = 0.0
    n_violations = 0

    if not specs.c3h6_pass:
        failures.append(
            f"C3H6 純度 {specs.c3h6_purity_wtfrac*100:.3f}% < spec "
            f"{config.spec.c3h6_min_wtfrac*100:.1f}wt% "
            f"(違反 {specs.c3h6_violation_pp:.3f}pp)"
        )
        total_violation_pp += specs.c3h6_violation_pp
        n_violations += 1
    if not specs.h2_pass:
        failures.append(
            f"H2 純度 {specs.h2_purity_molfrac*100:.3f}% < spec "
            f"{config.spec.h2_min_molfrac*100:.1f}mol% "
            f"(違反 {specs.h2_violation_pp:.3f}pp)"
        )
        total_violation_pp += specs.h2_violation_pp
        n_violations += 1
    if not specs.production_pass:
        # direction で分岐して正確に表示する (上限超過時に「生産量 1252 < 1176」の
        # ような数値矛盾表示にならないように)。
        if specs.production_direction == 'low':
            threshold_low_kmolh = specs.target_kmol_h * (1.0 - config.spec.production_min_relative)
            failures.append(
                f"生産量不足 {specs.production_kmol_h:.2f} < 下限 "
                f"target × {1.0 - config.spec.production_min_relative:.3f} = "
                f"{threshold_low_kmolh:.2f} kmol/h "
                f"(違反 {specs.production_under_pp:.3f}pp)"
            )
        elif specs.production_direction == 'high':
            threshold_high_kmolh = specs.target_kmol_h * (1.0 + config.spec.production_max_relative)
            failures.append(
                f"生産量超過 {specs.production_kmol_h:.2f} > 上限 "
                f"target × {1.0 + config.spec.production_max_relative:.3f} = "
                f"{threshold_high_kmolh:.2f} kmol/h "
                f"(超過 {specs.production_over_pp:.3f}pp)"
            )
        else:  # 'ok' なのに pass=False ということは無いはず (ガード)
            failures.append(
                f"生産量 {specs.production_kmol_h:.2f} (target {specs.target_kmol_h:.2f}) "
                f"direction={specs.production_direction!r}"
            )
        total_violation_pp += specs.production_violation_pp
        n_violations += 1

    # adaptive penalty scale。序盤は係数を弱めて探索を広げ、終盤は強化して
    # infeasible 領域から強制退出させる。詳細は optimization/penalty_scale.py。
    from optimization.penalty_scale import get_scale as _pscale
    _scale = _pscale()
    soft_penalty = (n_violations * pen.spec_base_okuyen * _scale
                    + total_violation_pp * pen.spec_coef_okuyen * _scale)

    # ---- (b') rigorous プロキシ罰則 ----
    # 各蒸留塔の equipment.proxy_penalty_okuyen を合算して soft_penalty にマージ。
    # FUG が「narrow margin / C3 漏れ過大」と判定した設計に対する追加コスト。
    # 詳細は src/distillation_core.py の _compute_proxy_penalty。
    # BO objective が「FUG では feasible でも rigorous で詰む領域」を自動的に避けるように誘導する。
    proxy_penalty_total = 0.0
    proxy_reasons: list[str] = []
    for col_key in ('r1', 'r2', 'r3'):
        col_result = solver_result.one_pass.get(col_key)
        if col_result is None or col_result.equipment is None:
            continue
        p = getattr(col_result.equipment, 'proxy_penalty_okuyen', 0.0)
        if p > 0:
            proxy_penalty_total += p
            reason = getattr(col_result.equipment, 'proxy_penalty_reason', '')
            if reason:
                proxy_reasons.append(f"{col_key}: {reason}")
    if proxy_penalty_total > 0:
        scaled_proxy = proxy_penalty_total * _scale
        soft_penalty += scaled_proxy
        failures.append(
            f"rigorous プロキシ罰則 +{scaled_proxy:.1f} 億円/年 ({' | '.join(proxy_reasons)})"
        )

    # ---- (b'') PSA/Mem trace bypass 連続 penalty ----
    # run_one_pass の _apply_trace_bypass が検出した「閾値超過分」を effective_TAC に
    # 加算。proxy_penalty は rigorous でしか発火しないため BO (FUG) では見えなかった
    # 「Dist2 が C3H6 を PSA に漏らす設計」を BO の探索段階で penalty として伝達する。
    #
    # trace_bypass は failures.append しない (= is_feasible を False にしない)。
    # 経済影響は連続罰金で十分カバー (1pp 超過 ≈ +10 億円 ≈ PSA 増床コスト相当)。
    # 罰金 + binary feasibility 失敗の二重打撃を避け、罰金だけ残し feasibility は許容する。
    one_pass_dict = solver_result.one_pass or {}
    psa_excess = one_pass_dict.get('trace_bypass_psa_excess', 0.0) or 0.0
    mem_excess = one_pass_dict.get('trace_bypass_mem_excess', 0.0) or 0.0
    trace_bypass_excess_total = psa_excess + mem_excess
    if trace_bypass_excess_total > 0:
        # trace_bypass を penalty_scale の緩和ランプ (序盤割引) から切り離す。
        # trace_bypass は「spec を序盤緩める」性質ではなく、PSA design モデルが
        # C3 を 1% 超扱えない = モデル妥当性の壁 (run_one_pass._TRACE_BYPASS_FRAC) なので、
        # 探索序盤でも割引すべきでない。max(_scale, 1.0) で設計意図 (10 億/pp) を
        # 下限保証しつつ終盤の強化は維持する。
        bypass_scale = max(_scale, 1.0)
        bypass_penalty = trace_bypass_excess_total * _TRACE_BYPASS_PENALTY_COEF_OKUYEN * bypass_scale
        soft_penalty += bypass_penalty
        # 注: failures.append は意図的に呼ばない。is_feasible は他の真の spec 違反
        # (純度・生産量) のみで判定する。trace_bypass の情報は user_attr に積んで
        # TPE constraints_func に流す (objective.py)。

    # ---- HI (post-processing) ----
    # apply_hi=True のときのみ pinch targeting を実行し、HI 後の OPEX/TAC/Profit を
    # 別 Economics として保持する。BO ループでは apply_hi=False (高速)、
    # top-k 候補に対してのみ apply_hi=True で再評価する。
    economics_hi = None
    hi_result    = None
    if apply_hi:
        from flowsheet.heat_integration import (
            extract_streams, pinch_analysis,
            get_default_utility_tiers, apply_hi_to_economics,
        )
        from src.cost_parameters import OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS

        streams = extract_streams(solver_result.one_pass, design.swing.T_in)
        heating_tiers, cooling_tiers = get_default_utility_tiers()
        hi_result = pinch_analysis(
            streams, dT_min_K=hi_dT_min_K,
            heating_tiers=heating_tiers, cooling_tiers=cooling_tiers,
        )
        economics_hi = apply_hi_to_economics(
            economics, hi_result, heating_tiers, cooling_tiers,
            operating_hours_per_year=OPERATING_HOURS_PER_YEAR,
            depreciation_years=DEPRECIATION_YEARS,
        )

    # ---- Stage 2 (HEN Synthesis) は不採用 (HI のみ) ----
    # ヒートインテグレーションは Stage 1 (pinch targeting) のみで評価し，
    # 実 HEN ネットワーク合成 (Stage 2) は行わない。apply_stage2 引数は後方互換のため
    # 残すが無効で，economics_synth は常に None (effective_TAC は HI のみの economics_hi を採用)。
    economics_synth = None
    hen_result      = None

    # ---- effective_TAC の選択 ----
    # 優先度: economics_synth (Stage 2) > economics_hi (Stage 1) > economics (raw)
    if economics_synth is not None:
        eff_econ = economics_synth
    elif economics_hi is not None:
        eff_econ = economics_hi
    else:
        eff_econ = economics
    # 純 TAC 最小化 (revenue 項は含めない)。production 制約は spec の下限ペナルティで
    # 担保される (n_violations × spec_base + violation_pp × spec_coef × scale)。
    # revenue を BO 評価から外すことで「装置小型化 → TAC ↓」の方向に純粋にバイアスがかかり、
    # production 下限を満たす最も安い設計に収束する。
    effective_TAC = eff_econ.TAC + soft_penalty

    # run_one_pass で捕捉した warning を取り出す (silent fallback 追跡用)
    one_pass = solver_result.one_pass or {}
    warnings_captured = one_pass.get('warnings_captured', []) or []

    # 観測ラベル: solver は完走したので spec 違反の種別を確定する。
    #   - 全 spec pass → 'success'
    #   - 違反あり → 主要違反 1 件で代表 (優先度: production > c3h6_purity > h2_purity)。
    #     複数違反時はカウントは別途 user_attrs に出すので、ここでは「最初に当てる方向」だけ。
    if not failures:
        failure_unit = 'success'
    elif not specs.production_pass:
        if specs.production_direction == 'low':
            failure_unit = 'spec_production_under'
        elif specs.production_direction == 'high':
            failure_unit = 'spec_production_over'
        else:
            failure_unit = 'spec_production'
    elif not specs.c3h6_pass:
        failure_unit = 'spec_c3h6_purity'
    elif not specs.h2_pass:
        failure_unit = 'spec_h2_purity'
    else:
        # proxy / trace_bypass のみで failures が立つことはあるが、これらは soft で
        # is_feasible に影響しない。ガード: ここに来たら spec 以外の何か。
        failure_unit = 'spec_other'

    return FlowsheetResult(
        solver=solver_result,
        economics=economics,
        specs=specs,
        effective_TAC=effective_TAC,
        failure_reason=" | ".join(failures),
        economics_hi=economics_hi,
        hi_result=hi_result,
        economics_synth=economics_synth,
        hen_result=hen_result,
        warnings_captured=warnings_captured,
        failure_unit=failure_unit,
    )


def _guess_unit_from_exception(exc_msg: str) -> str:
    """例外メッセージから装置名を best-effort で抽出。

    既知パターン:
      - "[Dist1] ...", "[Comp2/Dist2] ..." の角括弧プレフィクス (_capture_warnings の source)
      - "expansion_valve: ...", "simulate_pump: ..." の関数名プレフィクス
      - "PSA", "membrane", "psa_system", "membrane_system" などキーワード
      - "swing", "reactor" などキーワード
    マッチしなければ 'unknown'。
    """
    msg_lower = exc_msg.lower()
    for kw, label in (
        ('dist1', 'Dist1'), ('column1', 'Dist1'),
        ('dist2', 'Dist2'), ('column2', 'Dist2'),
        ('dist3', 'Dist3'), ('column3', 'Dist3'),
        ('psa', 'PSA'),
        ('membrane', 'Mem'), ('mem_system', 'Mem'),
        ('swing', 'Reactor'), ('reactor', 'Reactor'),
        ('expansion_valve', 'ExpansionValve'),
        ('compress', 'Compressor'),
        ('cooler', 'Cooler'),
        ('pump', 'Pump'),
    ):
        if kw in msg_lower:
            return label
    return 'unknown'
