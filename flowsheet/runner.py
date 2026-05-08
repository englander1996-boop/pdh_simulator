"""
最適化器から呼ぶ最上位エントリポイント: `evaluate(design, config) → Result`。

最適化器側の典型的な使い方:
    config = load_operating_config()
    design = FlowsheetDesignVars(swing=..., psa=..., mem=...)
    result = evaluate(design, config, verbose=False)
    objective = result.effective_TAC          # 最小化 (= TAC + ペナルティ)
    feasible  = result.is_feasible            # ペナルティ・spec チェックの集約

設計判断 (2026-05-08, 相談時の合意): 2階層ペナルティ構造
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
    """
    solver:         SolverResult
    economics:      Optional[Economics]            # solver 失敗時は None
    specs:          Optional[SpecComplianceResult] # solver 失敗時は None
    effective_TAC:  float                          # 最適化器が見る目的関数 [億円/年]
    failure_reason: str                            # "" のとき feasible

    @property
    def is_feasible(self) -> bool:
        return self.failure_reason == ""


def evaluate(
    design:  FlowsheetDesignVars,
    config:  OperatingConfig,
    verbose: bool = False,
) -> FlowsheetResult:
    """設計変数を入力してフローシートを評価し、effective_TAC と状態を返す。"""
    solver_result = solve_flowsheet(design, config, verbose=verbose)

    pen = config.penalty

    # ---- (a) solver-level 失敗 → ハード打ち切り ----
    # 設計判断: 結果の数値が信頼できない場合 (発散・暴走・未収束) は固定値で打ち切る。
    # 連続ペナルティを与えても情報量がなく、最適化器を惑わすだけのため。
    s = solver_result.inner_status
    if s.penalty_hit:
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=pen.solver_failure_okuyen,
            failure_reason="solver-level: PSA/Mem CAPEX ペナルティ発火",
        )
    if s.guard_hit:
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=pen.solver_failure_okuyen,
            failure_reason="solver-level: リサイクル暴走ガード発火",
        )
    if not s.converged or not solver_result.outer_status.converged:
        return FlowsheetResult(
            solver=solver_result, economics=None, specs=None,
            effective_TAC=pen.solver_failure_okuyen,
            failure_reason=(
                f"solver-level: 内側{'未' if not s.converged else ''}収束 / "
                f"外側{'未' if not solver_result.outer_status.converged else ''}収束"
            ),
        )

    # ---- solver 成功: 経済計算 + spec 判定 ----
    economics = calculate_economics(
        solver_result.one_pass,
        mw_C3H6_kg_per_kmol=config.product.mw_kg_per_kmol,
    )
    specs = check_specs(solver_result.one_pass, config)

    # ---- (b) spec 違反 → ソフトペナルティ ----
    # 設計判断: 全 spec 違反量を %pt スケールに揃えてあるため、係数を1つで管理。
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
        failures.append(
            f"生産量 {specs.production_kmol_h:.2f} < target × "
            f"{1.0 - config.spec.production_min_relative:.3f} = "
            f"{specs.target_kmol_h * (1.0 - config.spec.production_min_relative):.2f} kmol/h "
            f"(違反 {specs.production_violation_pp:.3f}pp)"
        )
        total_violation_pp += specs.production_violation_pp
        n_violations += 1

    soft_penalty = (n_violations * pen.spec_base_okuyen
                    + total_violation_pp * pen.spec_coef_okuyen)
    effective_TAC = economics.TAC + soft_penalty

    return FlowsheetResult(
        solver=solver_result,
        economics=economics,
        specs=specs,
        effective_TAC=effective_TAC,
        failure_reason=" | ".join(failures),
    )
