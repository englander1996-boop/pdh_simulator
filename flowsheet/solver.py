"""
リサイクル収束ソルバ (内側) + Fresh 調整ループ (外側)。

現状実装: 逐次置換 + アンダーリラックス。
将来的に Wegstein / Anderson 加速や scipy.optimize.root 統合に切り替えやすいよう
solver.py 1 ファイルに閉じ込めている。

設計NG (PSA/Mem CAPEX ペナルティ, リサイクル暴走) は SystemExit せず
SolverResult に状態として返す。最適化器側で penalty として扱えるようにするため。

設計判断 (2026-05-08, 相談時の合意):
  - 内側 TOL: 絶対値 → 相対値。スケール独立にしてノイズ追跡を回避。
    実測で内側 200 反復のうち後半 120 反復が ノイズ追跡だった。相対 0.1% に
    緩めて反復数を 1/3 〜 1/4 に短縮する見込み。
  - 外側 TOL: 両側絶対 → 片側相対。contest 仕様は「target 以上」で overshoot OK。
    片側化で外側反復を 4回 → 2-3 回に短縮見込み。
  - solver は流量の収束だけ判定し、製品品質 (純度) は flowsheet/specs.py で
    独立に評価する。
"""

from dataclasses import dataclass
from typing import Optional

from flowsheet.design import FlowsheetDesignVars
from flowsheet.run_one_pass import run_one_pass
from flowsheet.solvers import make_accelerator
from config.load import OperatingConfig


@dataclass
class TearState:
    """リサイクル収束の暫定状態 (warm-start にも使う)。"""
    tear_dist3: dict
    tear_mem:   dict
    T_d3:       float
    T_mem:      float


@dataclass
class InnerStatus:
    converged:   bool
    penalty_hit: bool
    guard_hit:   bool
    n_iter:      int
    final_diff:  float


@dataclass
class OuterStatus:
    converged:   bool
    n_iter:      int
    final_error: float


@dataclass
class SolverResult:
    one_pass:     dict             # 最終 run_one_pass の戻り値
    fresh_C3H8:   float
    fresh_C4H10:  float
    inner_status: InnerStatus
    outer_status: OuterStatus


def _operating_hours_per_year() -> float:
    """src/cost_parameters.py から取得 (将来 operating.toml へ移すかは要検討)。"""
    from src.cost_parameters import OPERATING_HOURS_PER_YEAR
    return OPERATING_HOURS_PER_YEAR


def _pack_state(tear_d3: dict, tear_mem: dict, T_d3: float, T_mem: float) -> dict:
    """TearState 相当の値を flat dict に変換 (TearAccelerator が扱いやすい形へ)。"""
    return {
        'd3_A':  tear_d3['A'],
        'd3_B':  tear_d3['B'],
        'mem_A': tear_mem['A'],
        'mem_B': tear_mem['B'],
        'T_d3':  T_d3,
        'T_mem': T_mem,
    }


def _unpack_state(flat: dict) -> 'TearState':
    """flat dict を TearState へ戻す。"""
    return TearState(
        tear_dist3={'A': flat['d3_A'],  'B': flat['d3_B']},
        tear_mem  ={'A': flat['mem_A'], 'B': flat['mem_B']},
        T_d3      =flat['T_d3'],
        T_mem     =flat['T_mem'],
    )


def _initial_tear(F_C3H8_feed: float, config: OperatingConfig) -> TearState:
    """Fresh 流量に比例した初期 tear stream を生成。"""
    s = config.solver.init
    scale = F_C3H8_feed / 1500.0
    return TearState(
        tear_dist3={'A': s.tear_d3_A_per_1500_fresh  * scale,
                    'B': s.tear_d3_B_per_1500_fresh  * scale},
        tear_mem  ={'A': s.tear_mem_A_per_1500_fresh * scale,
                    'B': s.tear_mem_B_per_1500_fresh * scale},
        T_d3      =s.T_d3_K,
        T_mem     =s.T_mem_K,
    )


def run_recycle_convergence(
    F_C3H8_feed:  float,
    F_C4H10_feed: float,
    design:       FlowsheetDesignVars,
    config:       OperatingConfig,
    init:         Optional[TearState] = None,
    verbose:      bool = True,
) -> tuple[dict, InnerStatus]:
    """指定 Fresh で内側リサイクルを収束させる。

    Returns
    -------
    (results, status)
        results: 最終 run_one_pass の戻り値 (収束/未収束問わず最後の状態)
        status:  収束フラグ・打ち切り理由・反復数
    """
    s = config.solver.inner
    if init is None:
        state = _initial_tear(F_C3H8_feed, config)
    else:
        state = TearState(
            tear_dist3=dict(init.tear_dist3),
            tear_mem  =dict(init.tear_mem),
            T_d3      =init.T_d3,
            T_mem     =init.T_mem,
        )

    # 設計判断: 加速法は config/operating.toml で選択 (Wegstein 推奨)。
    # 内部状態 (履歴) を持つため、新しい外側 iter ごとに reset() する必要がある。
    accelerator = make_accelerator(s.method, s)
    accelerator.reset()

    if verbose:
        print(f"  Fresh: C3H8={F_C3H8_feed:.2f}, C4H10={F_C4H10_feed:.2f} kmol/h"
              f"   (method={s.method}, TOL_rel={s.tol_relative*100:.3f}%,"
              f" floor={s.tol_floor_kmol_h}, MAX={s.max_iter})")
        print(f"  iter | tear_d3 (A,B) | tear_mem (A,B) | Δ_rel  | Reactor転化率")
        print(f"  -----+---------------+----------------+--------+--------------")

    results     = None
    converged   = False
    penalty_hit = False
    guard_hit   = False
    diff        = 0.0
    it          = 0

    for it in range(1, s.max_iter + 1):
        results = run_one_pass(
            state.tear_dist3, state.tear_mem,
            state.T_d3, state.T_mem,
            F_C3H8_feed, F_C4H10_feed,
            design, config,
        )

        if (results['r_psa'].equipment.CAPEX_total >= 1e8 or
            results['r_mem'].equipment.CAPEX_total >= 1e8):
            if verbose:
                print(f"  {it:4d} | --- PSA/Mem ペナルティ発火 → 設計変数の見直しが必要 ---")
            penalty_hit = True
            break

        # 設計判断: 相対 TOL 計算で小流量に対し過度に厳しくならないよう floor を入れる。
        #   denominator = max(|tear|, tol_floor_kmol_h)
        # tear=7500 → floor 効かず通常の相対精度、tear=30 → floor が効いて緩める。
        floor = s.tol_floor_kmol_h
        diff_d3 = max(
            abs(results['tear_dist3_new'][k] - state.tear_dist3[k])
            / max(abs(state.tear_dist3[k]), floor)
            for k in state.tear_dist3
        )
        diff_mem = max(
            abs(results['tear_mem_new'][k] - state.tear_mem[k])
            / max(abs(state.tear_mem[k]), floor)
            for k in state.tear_mem
        )
        diff = max(diff_d3, diff_mem)   # 相対誤差 [-]

        conv = results['r_rx'].performance.Conversion
        if verbose:
            print(f"  {it:4d} |"
                  f" {results['tear_dist3_new']['A']:5.2f},{results['tear_dist3_new']['B']:5.2f}"
                  f"   | {results['tear_mem_new']['A']:6.2f},{results['tear_mem_new']['B']:5.2f}"
                  f"   | {diff*100:5.3f}% | {conv:5.1f}%")

        if results['tear_mem_new']['A'] > F_C3H8_feed * s.recycle_guard_ratio:
            if verbose:
                print(f"  → リサイクル暴走ガード発火"
                      f" (tear_mem.A > {F_C3H8_feed * s.recycle_guard_ratio:.0f} kmol/h)")
            guard_hit = True
            break

        if diff < s.tol_relative:
            converged = True
            if verbose:
                print(f"  → 内側収束 (Δ_rel={diff*100:.4f}% < TOL_rel={s.tol_relative*100:.3f}%)")
            break

        # 設計判断: TearAccelerator (SS or Wegstein) に次反復の更新を委譲。
        # state を flat dict にパックして渡し、結果を unpack して次の反復へ。
        flat_current  = _pack_state(state.tear_dist3,        state.tear_mem,
                                     state.T_d3,             state.T_mem)
        flat_computed = _pack_state(results['tear_dist3_new'], results['tear_mem_new'],
                                     results['T_d3_new'],   results['T_mem_new'])
        flat_next     = accelerator.step(flat_current, flat_computed)
        state         = _unpack_state(flat_next)

    if verbose and not (converged or penalty_hit or guard_hit):
        print(f"  → 内側未収束 ({s.max_iter} 回打ち切り、最終状態で集計)")

    return results, InnerStatus(
        converged=converged, penalty_hit=penalty_hit,
        guard_hit=guard_hit, n_iter=it, final_diff=diff,
    )


def solve_flowsheet(
    design:  FlowsheetDesignVars,
    config:  OperatingConfig,
    verbose: bool = True,
    F_C3H8_override: Optional[float] = None,
) -> SolverResult:
    """外側ループ (Fresh調整) + 内側ループ (リサイクル収束) を統合して解く。

    Parameters
    ----------
    F_C3H8_override : float | None
        指定した場合、外側ループを skip して内側のみ実行。
        Fresh C3H8 = F_C3H8_override [kmol/h] を強制使用。
        Fresh C4H10 は LPG 組成 (config.feed.lpg_c3h8_mol_fraction) から自動計算。
        BO で F_fresh を直接最適化変数とする場合に使う。
        None なら従来通り外側ループで Fresh を調整。
    """
    so = config.solver.outer
    p  = config.product
    f  = config.feed

    operating_hours = _operating_hours_per_year()
    target_kmol_h   = p.target_mta * 1000.0 / p.mw_kg_per_kmol / operating_hours
    c4_over_c3      = (1.0 - f.lpg_c3h8_mol_fraction) / f.lpg_c3h8_mol_fraction

    if F_C3H8_override is not None:
        F_C3H8_feed = float(F_C3H8_override)
    else:
        F_C3H8_feed = target_kmol_h / f.yield_assumed
    F_C4H10_feed = F_C3H8_feed * c4_over_c3

    if verbose:
        print(f"  C3H6 目標     : {target_kmol_h:.2f} kmol/h"
              f"  ({p.target_mta:.0f} t/年, {operating_hours:.0f} h/年)")
        print(f"  仮定収率      : {f.yield_assumed}"
              f"  → 初期 Fresh C3H8 = {F_C3H8_feed:.2f} kmol/h")
        print(f"  LPG 組成      : C3H8 mol fraction = {f.lpg_c3h8_mol_fraction}"
              f" (Fresh C4H10 / Fresh C3H8 = {c4_over_c3:.4f})")
        print(f"  原料状態      : {f.T_K - 273.15:.1f}°C, {f.P_Pa / 1e5:.3f} bar")
        print(f"  外側収束基準  : actual >= target × {1.0 - so.tol_relative:.3f}"
              f"   (片側 {so.tol_relative*100:.1f}% 不足まで, overshoot 任意)"
              f"  (RELAX={so.relax}, MAX={so.max_iter})")

    inner_init      = None
    results         = None
    inner_status    = None
    outer_converged = False
    error           = 0.0
    it              = 0

    # ---- F_fresh override: 外側ループ skip ----
    # BO で F_fresh を直接最適化変数とした場合。production_min spec の充足は
    # spec check (flowsheet/specs.py) 側で soft penalty として評価される。
    if F_C3H8_override is not None:
        if verbose:
            print(f"\n[F_fresh override] Fresh C3H8={F_C3H8_feed:.2f}, C4H10={F_C4H10_feed:.2f} "
                  f"kmol/h (外側ループ skip)")
        # 設計判断 (2026-05-14): 初期 tear stream は yield_assumed=0.9 ベースの F_fresh
        # (= target/0.9 ≈ 1320) で計算。実 F_fresh をそのまま scale すると、高い F_fresh
        # で初期 tear が過大 → iter 1 で Mem A 過大 → PSA/Mem CAPEX ペナルティ発火。
        # 外側ループ初期化挙動と整合させるため、初期 tear は target ベースで固定 scale。
        F_C3H8_init_for_tear = target_kmol_h / f.yield_assumed
        init_tear = _initial_tear(F_C3H8_init_for_tear, config)
        results, inner_status = run_recycle_convergence(
            F_C3H8_feed, F_C4H10_feed,
            design, config,
            init=init_tear,
            verbose=verbose,
        )
        actual_product = results['r3'].top.F_in.get('B', 0.0) if results is not None else 0.0
        error = actual_product - target_kmol_h
        return SolverResult(
            one_pass     =results,
            fresh_C3H8   =F_C3H8_feed,
            fresh_C4H10  =F_C4H10_feed,
            inner_status =inner_status,
            outer_status =OuterStatus(
                # override 経路は「外側 1 iter」扱いで converged 扱い
                # (実際の spec 充足は spec check 側で評価)
                converged=True, n_iter=1, final_error=error,
            ),
        )

    for it in range(1, so.max_iter + 1):
        if verbose:
            print(f"\n[外側 iter {it}]")

        results, inner_status = run_recycle_convergence(
            F_C3H8_feed, F_C4H10_feed,
            design, config,
            init=inner_init,
            verbose=verbose,
        )

        if inner_status.penalty_hit or inner_status.guard_hit:
            error = float('inf')
            break

        actual_product  = results['r3'].top.F_in.get('B', 0.0)
        error           = actual_product - target_kmol_h
        effective_yield = actual_product / F_C3H8_feed if F_C3H8_feed > 0 else 0.0

        # 設計判断: 片側収束基準。shortfall = max(0, target - actual)。
        # overshoot は許容 (contest 仕様で「target 以上」のため)。
        # Fresh 調整ループは target に張り付く挙動なので、片側でも不必要な反復は
        # 起きず、undershoot 領域のみ厳密に潰すことになる。
        shortfall = max(0.0, target_kmol_h - actual_product)
        threshold = target_kmol_h * so.tol_relative

        if verbose:
            print(f"  → 実生産 {actual_product:7.2f} kmol/h"
                  f" (誤差 {error:+7.2f}, 不足 {shortfall:6.2f},"
                  f" 実収率 {effective_yield*100:5.2f}%)")

        if shortfall < threshold:
            outer_converged = True
            if verbose:
                print(f"  ✓ 外側収束 (不足 {shortfall:.3f} < {threshold:.3f} kmol/h"
                      f" = target × {so.tol_relative*100:.1f}%)")
            break

        if effective_yield <= 0:
            if verbose:
                print("  → 実収率 0 のため Fresh 更新不能。")
            error = float('inf')
            break

        F_C3H8_new   = target_kmol_h / effective_yield
        F_C3H8_feed  = so.relax * F_C3H8_new + (1 - so.relax) * F_C3H8_feed
        F_C4H10_feed = F_C3H8_feed * c4_over_c3

        inner_init = TearState(
            tear_dist3=results['tear_dist3_new'],
            tear_mem  =results['tear_mem_new'],
            T_d3      =results['T_d3_new'],
            T_mem     =results['T_mem_new'],
        )

    if verbose and not outer_converged:
        print(f"\n  → 外側未収束 ({so.max_iter} 回打ち切り、最終状態で集計)")

    return SolverResult(
        one_pass     =results,
        fresh_C3H8   =F_C3H8_feed,
        fresh_C4H10  =F_C4H10_feed,
        inner_status =inner_status,
        outer_status =OuterStatus(
            converged=outer_converged, n_iter=it, final_error=error,
        ),
    )
