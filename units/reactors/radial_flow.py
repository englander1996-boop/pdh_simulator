"""
PDH 径方向流 (radial flow) スイング反応器システム シミュレーター

========================================================================
位置づけ (2026-05-30)
========================================================================
反応器設計レビューで「0.5 bar 低圧固定床では圧力損失が支配的になり、軸流深床
(units/reactors/swing.py) は成立しない」と判明した (詳細は swing.py ヘッダ +
monitor/reactor_pressure_drop_and_geometry.ipynb)。

その物理的な解が **径方向流 (radial flow) ベッド**:
  - ガスを薄い環状床に「半径方向に短く」通す → 流路長 = 床厚 Δr (~0.3-1m) と短い
  - 流路断面積 = 2π r H (塔高 H を稼げば大きい) → 空塔速度が低い
  - 触媒量 (= 塔高 H と床厚 Δr で決まる) と 圧損 (= 床厚 Δr のみで決まる) を **分離**できる
  ⇒ 体積 vs 圧損のトレードオフが解消し、現実的な触媒粒径 (3mm) でも 0.5 bar で ΔP が収まる。
これは実機 UOP Oleflex (径方向流) / Lummus Catofin (減圧・浅床) の設計思想そのもの。

========================================================================
設計方針: swing.py (軸流) は温存し、本ファイルを **新規追加** で並走させる
========================================================================
  - 反応速度・熱力学・Ergun 物性・コスト計算・出力データクラスは全て swing.py から
    import して共有 (二重実装・二重メンテを避ける)。
  - 異なるのは「幾何」だけ: 独立変数を軸 z → 半径 r に、断面積 A=πD²/4(一定) →
    A(r)=2πrH(r 依存) に差し替える。化学(反応・熱・失活)は軸流と同一。
  - したがって同じ V_cat・同じ feed なら 転化率/選択率/温度プロファイルは軸流とほぼ
    同一になり、**差が出るのは ΔP だけ** (= 径方向流が「圧損だけ」を救う、という主張を担保)。
  - 出力は swing と同じ SimulationResult (EffluentStream/EquipmentCost/PerformanceMetrics)
    を返すので、下流 (run_one_pass/economics) は反応器を差し替えるだけで動く。

========================================================================
モデル (centripetal: 外周 r_o から中心 r_i へ流す)
========================================================================
状態 y = [F_A..F_F, T, P]、独立変数 r を r_o → r_i へ積分 (r 減少方向)。
  dF_i/dr = -(1-eps)·2π r H · Σ_j ν_ij r_j
  dT/dr   = +(1-eps)·2π r H · Σ_j ΔH_j r_j / Σ_i F_i Cp_i      (断熱)
  dP/dr   = +[ Ergun(u, ρ, μ) ]                                (流れ方向に圧損)
ここで u = Q_vol / (2π r H · N_parallel) は局所空塔速度 (r 内側ほど面積↓→速度↑)。
符号は「r 減少方向 = 流れ方向」で反応物消費・降温・降圧が起こるよう取ってある。

⚠️ モデル限界 (swing.py と共通、レポート記載): 段間 reheat 未計上、再生動特性なし、
   粒内拡散 (Weisz-Prater) 未考慮、触媒活性は入口温度で代表し空間一定。
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# --- 物理・コスト・出力は swing.py から共有 (二重実装回避) ---
from units.reactors.swing import (
    calc_a, calc_Cp, calc_rate_constants, _reaction_enthalpies,
    _gas_viscosity, _STOICH, _COMPS, _R_GAS, _P_FLOOR_PA,
    FeedStream, FixedParams, EffluentStream, EquipmentCost,
    PerformanceMetrics, SimulationResult, _penalty_result, _PENALTY_CAPEX,
    calc_reactor_capex_okuyen,
)
from src.component_data import MW
from src.config import THERMO_DATA
from src.thermo import PDHThermo

_thermo = PDHThermo()


# ---------------------------------------------------------------------------
# 設計変数 (径方向流 固有)
# ---------------------------------------------------------------------------

@dataclass
class RadialDesignVars:
    """径方向流反応器の設計変数。

    軸流 (swing.DesignVars) の (T_in, z_cat, t_cyc, D) に対し、径方向流は幾何が
    異なるため (T_in, t_cyc, D_inner, bed_thickness, H) を持つ。
    """
    T_in:          float   # 反応器入口温度 [K]
    t_cyc:         float   # 1サイクル反応フェーズ運転時間 [min]
    D_inner:       float   # 環状床の内径 [m] (中心捕集管の外径。r_i = D_inner/2)
    bed_thickness: float   # 環状床の厚さ Δr [m] (r_o = r_i + bed_thickness)
    H:             float   # 環状床の高さ [m]

    @property
    def r_i(self) -> float:
        return self.D_inner / 2.0

    @property
    def r_o(self) -> float:
        return self.r_i + self.bed_thickness


# ---------------------------------------------------------------------------
# 軸方向(半径方向)の常微分方程式
# ---------------------------------------------------------------------------

def _ode_radial(r: float, y: np.ndarray,
                a: float, H: float, eps: float,
                eps_bed: float, d_p: float, sphericity: float, N_parallel: int,
                ) -> np.ndarray:
    """半径方向(r)の常微分方程式 (centripetal: r_o → r_i へ積分)。

    State y = [F_A..F_F, T, P]、単位 F [mol/s], T [K], P [Pa]。
    断面積 A(r) = 2π r H は r 依存。速度式・熱・Ergun は swing._ode_axial と同一物理。
    """
    F_raw = y[:6]
    F = np.maximum(F_raw, 0.0)
    T_local = float(np.clip(y[6], 300.0, 1500.0))
    P_local = max(float(y[7]), _P_FLOOR_PA)

    F_total = float(np.sum(F))
    if F_total <= 0 or r <= 0:
        return np.zeros(8)

    # 局所断面積 (環状床の円筒側面積)
    A_r = 2.0 * math.pi * r * H
    cat_factor = (1.0 - eps) * A_r   # [m²] 触媒断面積相当 (per 単位 r)

    # 分圧 [Pa] (局所全圧基準)
    P = {comp: max(float(F[i]) / F_total * P_local, 0.0)
         for i, comp in enumerate(_COMPS)}

    # --- 反応速度 (swing._ode_axial の計算を踏襲) ---
    rc = calc_rate_constants(T_local)
    k1, k2, k3 = rc['k1'], rc['k2'], rc['k3']
    K_B  = max(rc['K_B'],  1.0)   # [Pa] 低温での吸着項ゼロ除算防止 (swing と同じ数値ガード)
    K_eq = max(rc['K_eq'], 1.0)   # [Pa] 低温での駆動力発散防止
    driving_r1 = P['A'] - P['B'] * P['C'] / K_eq
    r1 = a * k1 * driving_r1 / (1.0 + P['B'] / K_B)
    r2 = k2 * P['A']
    r3 = k3 * P['D'] * P['C']
    rates = np.array([r1, r2, r3])   # [mol/m³_cat/s]

    # 物質収支: dF_i/dr = -(1-eps)·2π r H · Σ ν_ij r_j  (r_o→r_i 積分で反応物が減る符号)
    dFdr = -cat_factor * (_STOICH @ rates)

    # エネルギー収支(断熱): dT/dr = +(1-eps)·2π r H · Σ ΔH_j r_j / Σ F_i Cp_i
    cp_dict = calc_Cp(T_local)
    sum_FCp = sum(max(float(F[i]), 0.0) * cp_dict[comp] for i, comp in enumerate(_COMPS))
    if sum_FCp <= 0:
        dTdr = 0.0
    else:
        dH = _reaction_enthalpies(T_local)
        Q_rxn = float(cat_factor * np.dot(dH, rates))   # [J/(m·s)]
        dTdr = Q_rxn / sum_FCp                           # [K/m] (吸熱で r 減少方向に降温)

    # 圧力損失 (Ergun)。u = Q_vol / (A_r × N_parallel) は 1 基あたり局所空塔速度。
    # r_o→r_i 積分で P が下がるよう dP/dr は正 (流れ方向 = r 減少方向で圧損)。
    if A_r > 0 and N_parallel >= 1:
        Q_vol = F_total * _R_GAS * T_local / P_local
        u = Q_vol / (A_r * N_parallel)
        mass_flow = sum(float(F[i]) * MW[comp] for i, comp in enumerate(_COMPS)) / 1000.0
        rho = mass_flow / Q_vol if Q_vol > 0 else 0.0
        mu = _gas_viscosity(T_local)
        phi_dp = sphericity * d_p
        eb = eps_bed
        visc_term = 150.0 * (1.0 - eb) ** 2 * mu * u / (eb ** 3 * phi_dp ** 2)
        inert_term = 1.75 * (1.0 - eb) * rho * u ** 2 / (eb ** 3 * phi_dp)
        dPdr = (visc_term + inert_term)
        if float(y[7]) <= _P_FLOOR_PA:
            dPdr = 0.0   # 床下限ガード (P が負へ暴走するのを防ぐ。swing と同方針)
    else:
        dPdr = 0.0

    return np.concatenate([dFdr, [dTdr], [dPdr]])


def _simulate_one_time_radial(design: RadialDesignVars, feed: FeedStream,
                              fixed: FixedParams, t_min: float,
                              N_parallel: int) -> tuple:
    """時刻 t_min [min] における半径方向積分。

    Returns (F_out [mol/s], T_out [K], P_out [Pa]) at r_i。失敗時 (None,None,None)。
    """
    a = calc_a(t_min, design.T_in, feed.P_in)
    F0 = np.array([feed.F_in.get(c, 0.0) * 1000.0 / 3600.0 for c in _COMPS])
    y0 = np.concatenate([F0, [design.T_in], [feed.P_in]])

    try:
        # r_o → r_i へ積分 (centripetal、r 減少方向)
        sol = solve_ivp(
            fun=lambda r, y: _ode_radial(
                r, y, a, design.H, fixed.eps,
                fixed.eps_bed, fixed.d_p_m, fixed.sphericity, N_parallel,
            ),
            t_span=(design.r_o, design.r_i),
            y0=y0,
            method='Radau',
            rtol=1e-4,
            atol=1e-7,
        )
    except Exception as e:
        warnings.warn(
            f"radial.solve_ivp 例外: {type(e).__name__}: {e} "
            f"(T_in={design.T_in:.1f}K, r_i={design.r_i:.2f}m, r_o={design.r_o:.2f}m, "
            f"H={design.H:.2f}m, t={t_min:.1f}min)",
            RuntimeWarning, stacklevel=2,
        )
        return None, None, None

    if not sol.success:
        warnings.warn(
            f"radial.solve_ivp 収束失敗 (status={sol.status}, message={sol.message!r}) "
            f"(T_in={design.T_in:.1f}K, r_i={design.r_i:.2f}m, r_o={design.r_o:.2f}m, "
            f"H={design.H:.2f}m, t={t_min:.1f}min)",
            RuntimeWarning, stacklevel=2,
        )
        return None, None, None

    return sol.y[:6, -1], float(sol.y[6, -1]), float(sol.y[7, -1])


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_radial_flow_reactor_system(
    design: RadialDesignVars,
    feed: FeedStream,
    fixed: FixedParams,
    n_time_samples: int = 20,
) -> SimulationResult:
    """径方向流スイング反応器システムをシミュレーションする。

    swing.simulate_swing_reactor_system と同じ SimulationResult を返す
    (下流 run_one_pass/economics はそのまま使える)。差は幾何 (径方向) のみ。

    径方向流の SV チェック: 軸流の [SV_min, SV_max] のうち **SV_max のみ** を、
    速度最大の内側面 (r_i) について課す (高速面でのフルイダイゼーション/同伴の上限)。
    SV_min は課さない: 径方向流は塔高 H を稼ぐと空塔速度が下がるのが本来の利点で、
    低速はディストリビュータ設計で扱う想定 (軸流の channeling 懸念とは事情が異なる)。
    """
    # ---- 入力バリデーション ----
    if (design.t_cyc <= 0 or design.D_inner <= 0 or design.bed_thickness <= 0
            or design.H <= 0):
        return _penalty_result(reason='input_invalid')
    if design.T_in <= 0 or feed.T_feed <= 0 or feed.P_in <= 0:
        return _penalty_result(reason='input_invalid')
    if any(v < 0 for v in feed.F_in.values()):
        return _penalty_result(reason='input_invalid')
    if sum(feed.F_in.values()) <= 0:
        return _penalty_result(reason='input_invalid')

    r_i, r_o, H = design.r_i, design.r_o, design.H

    # ---- ジオメトリ / 並列基数 (時間ループ前に確定、swing と同方針) ----
    # 環状床の触媒体積 V_cat = π(r_o² - r_i²)·H·(1-eps)。
    V_cat_total = math.pi * (r_o ** 2 - r_i ** 2) * H * (1.0 - fixed.eps)
    N_parallel = max(math.ceil(V_cat_total / fixed.V_cat_max_per_vessel), 1)

    # ---- 時間方向サンプリングと半径方向積分 ----
    if n_time_samples < 2:
        warnings.warn(
            f"n_time_samples={n_time_samples} < 2: 台形積分を行わず t=0 の1点を採用します。",
            UserWarning, stacklevel=2,
        )
        F_out, T_out, P_out = _simulate_one_time_radial(design, feed, fixed, 0.0, N_parallel)
        if F_out is None:
            return _penalty_result(reason='sim_failure')
        F_out_avg_kmolh = {comp: float(F_out[i]) * 3600.0 / 1000.0
                           for i, comp in enumerate(_COMPS)}
        T_out_avg = T_out
        P_out_avg = P_out
    else:
        t_samples = np.linspace(0.0, design.t_cyc, n_time_samples)
        F_out_list, T_out_list, P_out_list = [], [], []
        for t in t_samples:
            F_out, T_out, P_out = _simulate_one_time_radial(design, feed, fixed, float(t), N_parallel)
            if F_out is None:
                return _penalty_result(reason='sim_failure')
            F_out_list.append(F_out)
            T_out_list.append(T_out)
            P_out_list.append(P_out)
        F_out_arr = np.array(F_out_list)
        T_out_arr = np.array(T_out_list)
        P_out_arr = np.array(P_out_list)
        _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
        F_out_avg_mol_s = _trapz(F_out_arr, t_samples, axis=0) / design.t_cyc
        T_out_avg = float(_trapz(T_out_arr, t_samples) / design.t_cyc)
        P_out_avg = float(_trapz(P_out_arr, t_samples) / design.t_cyc)
        F_out_avg_kmolh = {comp: float(F_out_avg_mol_s[i]) * 3600.0 / 1000.0
                           for i, comp in enumerate(_COMPS)}

    # ---- 圧力損失 (Ergun) ハード制約 ----
    dP_over_P = (feed.P_in - P_out_avg) / feed.P_in if feed.P_in > 0 else 0.0
    dP_over_P = float(np.clip(dP_over_P, 0.0, 1.0))
    if dP_over_P > fixed.dP_over_P_max:
        warnings.warn(
            f"radial reactor: ΔP/P_in={dP_over_P*100:.1f}% が上限 "
            f"{fixed.dP_over_P_max*100:.0f}% 超 (r_i={r_i:.2f}m, r_o={r_o:.2f}m, "
            f"H={H:.1f}m, d_p={fixed.d_p_m*1e3:.1f}mm) — infeasible 化",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='dP_excess', dP_over_P=dP_over_P)

    # ---- 予熱熱量 Q_preheat [GJ/h] ----
    q_w = 0.0
    for comp in _COMPS:
        F_mol_s = feed.F_in.get(comp, 0.0) * 1000.0 / 3600.0
        q_w += F_mol_s * _thermo.calc_enthalpy_change(comp, feed.T_feed, design.T_in)
    Q_preheat_GJh = q_w * 3600.0 / 1e9

    # ---- 空塔速度チェック (内側面 r_i = 最大速度面に SV_max のみ) ----
    n_inlet_mol_s = sum(feed.F_in.values()) * 1000.0 / 3600.0
    Q_vol_in = n_inlet_mol_s * _R_GAS * design.T_in / feed.P_in if feed.P_in > 0 else 0.0
    A_inner = 2.0 * math.pi * r_i * H        # 最小断面 (最大速度)
    SV_inner = Q_vol_in / (A_inner * N_parallel) if A_inner > 0 else 0.0
    if SV_inner > fixed.SV_max_m_per_s:
        warnings.warn(
            f"radial reactor: 内側面 SV={SV_inner:.2f} m/s が上限 "
            f"{fixed.SV_max_m_per_s} m/s 超 (r_i={r_i:.2f}m, H={H:.1f}m, "
            f"N_parallel={N_parallel}) — infeasible 化 (フルイダイゼーション/同伴懸念)",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='sv_out_of_range', SV_actual=SV_inner)

    N_swing_sets = math.ceil(fixed.t_regen / design.t_cyc) + 1
    N_reactors_total = N_parallel * N_swing_sets

    # 容器体積。CAPEX 計上体積からは中心捕集管 (半径 r_i) の void を除外し、
    # 触媒層を含む環状部のみを計上する (ユーザー決定 2026-05-30)。
    # この定義では V_vessel = V_cat/(1-eps) となり、軸流の容器体積規約と整合する。
    V_vessel_actual = (math.pi * (r_o ** 2 - r_i ** 2) * H) / N_parallel  # [m³] 1基あたり (環状部)
    catalyst_weight_total = V_cat_total * N_swing_sets * fixed.rho_b  # [kg]

    if V_vessel_actual <= 0:
        return _penalty_result(reason='volume_zero')

    # CAPEX: 縦型プロセス容器として推算 (外径 = 2 r_o)
    try:
        reactor_capex = calc_reactor_capex_okuyen(
            V_vessel_m3=V_vessel_actual,
            P_abs_pa=feed.P_in,
            D_m=2.0 * r_o,
            N_reactors_total=N_reactors_total,
        )
    except Exception:
        reactor_capex = _PENALTY_CAPEX

    # ---- パフォーマンス指標 ----
    F_A_in  = feed.F_in.get('A', 0.0)
    F_A_out = F_out_avg_kmolh['A']
    F_B_in  = feed.F_in.get('B', 0.0)
    F_B_out = F_out_avg_kmolh['B']
    conversion  = (F_A_in - F_A_out) / F_A_in * 100.0 if F_A_in > 0 else 0.0
    delta_A     = F_A_in - F_A_out
    selectivity = (F_B_out - F_B_in) / delta_A * 100.0 if delta_A > 0 else 0.0
    conversion  = float(np.clip(conversion,  0.0, 100.0))
    selectivity = float(np.clip(selectivity, 0.0, 100.0))

    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg=F_out_avg_kmolh,
            T_out_avg=T_out_avg,
            Q_preheat=Q_preheat_GJh,
            P_out=P_out_avg,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=V_vessel_actual,
            N_parallel=N_parallel,
            N_swing_sets=N_swing_sets,
            N_reactors_total=N_reactors_total,
            Catalyst_Weight_Total=catalyst_weight_total,
            Reactor_CAPEX=reactor_capex,
            dP_over_P_actual=dP_over_P,
        ),
        performance=PerformanceMetrics(
            Conversion=conversion,
            Selectivity=selectivity,
        ),
    )
