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
ここで u = Q_vol / (2π r H) は局所空塔速度 (r 内側ほど面積↓→速度↑)。
符号は「r 減少方向 = 流れ方向」で反応物消費・降温・降圧が起こるよう取ってある。

【並列分割の規約 (案B, 2026-05-31 確定)】 設計変数 (D_inner, Δr, H) は「作る反応器全体」を
表す。触媒体積 V_cat が V_cat_max(200m³)を超える場合、総断面 2πrH を N_parallel 本のサブ
カラムに機械分割する。総流量 ÷ 総断面 = 各サブカラムの空塔速度は分割で不変なので、u・ΔP・SV
は N_parallel に依存しない (旧実装は u を /N_parallel しており ΔP・SV を過小評価するバグだった)。
N_parallel は容器分割・コスト計上 (V_vessel=V_design/N_parallel, 基数=N_parallel·N_swing) のみに使う。

⚠️ モデル限界 (swing.py と共通、レポート記載): 再生動特性なし、粒内拡散 (Weisz-Prater) 未考慮、
   触媒活性は入口温度で代表し空間一定、床 Ergun 以外の ΔP (分配器/中心管/弁/段間配管) 未計上。
   なお段間 reheat は simulate_radial_multibed_reactor_system で計上済 (本ファイル後半)。
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

    # 圧力損失 (Ergun)。u = Q_vol / A_r は局所空塔速度 (案B: N_parallel で割らない)。
    #   案B (設計環状床 = 反応器全体、200m³超は分割) では、総断面 2πrH を N_parallel 本の
    #   サブカラムに分割しても「総流量 ÷ 総断面」= 各カラムの空塔速度は不変 (Q/N ÷ A/N = Q/A)。
    #   よって ΔP・SV は N_parallel に依存しない (旧実装の /N_parallel は過小評価のバグだった)。
    #   N_parallel は容器分割・コスト計上のみに使う量で、ここ (物理) では使わない。
    # r_o→r_i 積分で P が下がるよう dP/dr は正 (流れ方向 = r 減少方向で圧損)。
    if A_r > 0:
        Q_vol = F_total * _R_GAS * T_local / P_local
        u = Q_vol / A_r
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

    # ---- 圧力損失 (床 Ergun × 内部品マージン) ハード制約 ----
    # 床ΔP に dP_margin_factor を掛けて総反応器ΔP とし、出口圧・feasibility に適用 (化学は真の床P)。
    dP_over_P_bed = (feed.P_in - P_out_avg) / feed.P_in if feed.P_in > 0 else 0.0
    dP_over_P = float(np.clip(dP_over_P_bed * fixed.dP_margin_factor, 0.0, 1.0))
    P_out_avg = feed.P_in * (1.0 - dP_over_P)   # マージン込み出口圧 (→ 下流圧縮機へ伝播)
    if dP_over_P > fixed.dP_over_P_max:
        warnings.warn(
            f"radial reactor: ΔP/P_in={dP_over_P*100:.1f}% (床×{fixed.dP_margin_factor}) が上限 "
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
    # 案B: SV も N_parallel で割らない (総流量 ÷ 総内側断面 = 各サブカラムの空塔速度に等しい)。
    n_inlet_mol_s = sum(feed.F_in.values()) * 1000.0 / 3600.0
    Q_vol_in = n_inlet_mol_s * _R_GAS * design.T_in / feed.P_in if feed.P_in > 0 else 0.0
    A_inner = 2.0 * math.pi * r_i * H        # 最小断面 (最大速度)
    SV_inner = Q_vol_in / A_inner if A_inner > 0 else 0.0
    if SV_inner > fixed.SV_max_m_per_s:
        warnings.warn(
            f"radial reactor: 内側面 SV={SV_inner:.2f} m/s が上限 "
            f"{fixed.SV_max_m_per_s} m/s 超 (r_i={r_i:.2f}m, H={H:.1f}m) "
            f"— infeasible 化 (フルイダイゼーション/同伴懸念。SV は H↑/r_i↑ で下げる)",
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


# ---------------------------------------------------------------------------
# 多段 (Oleflex 型: 径方向流断熱床 N 段直列 + 段間再加熱)
# ---------------------------------------------------------------------------
#
# ■ なぜこの実装に至ったか (2026-05-31, レポート記載用)
#   単段断熱の径方向流は、PDH の強い吸熱で床が降温し、出口温度での平衡に張り付いて
#   per-pass 転化率が ~30% で頭打ちになる (実測: T_in=939K でも 30%、触媒量を増やしても
#   30% で飽和)。この低転化が巨大なプロパン recycle を生み、膜フィードを 18% propylene まで
#   希釈 → 膜回収率 30% → 未回収 propylene が反応器で副生に化け収率 46% + Dist3 SM 分類器の
#   学習域 (In_Flow≥0.36) を割り込んで main BO が feasible 0 になっていた。
#
#   実機 PDH はこれを「単段で解かない」: UOP Oleflex は径方向流断熱反応器を 3〜4 基直列に
#   並べ、各反応器の間に加熱炉 (interstage heater) を置き、吸熱で冷えたガスを反応温度まで
#   戻してから次段へ送る。各段で平衡がリセットされ累積転化が上がる (実測: 3 段で ~59%)。
#   本関数はこの Oleflex 型を、既存の単段ソルバを N 回直列に呼ぶだけで表現する。
#
# ■ 時刻同期で直列計算する (2026-05-31, 反応器設計プロ指摘 #5 で改修)
#   スイング列は時刻 t ごとに Bed1(t)→reheat→Bed2(t)→reheat→Bed3(t) と流れる。各段の活性
#   a(t) は同一 t (同時に新鮮再生) で決まる。反応は非線形なので「各段で時間平均してから次段へ
#   渡す」と f(平均)≠平均(f) のズレ (特に選択率) が出る。よって本実装は **時刻サンプル t ごとに
#   全段を直列積分し、最終段出口を最後に時間平均する**。段間再加熱の熱量も t ごとに評価して平均
#   する (前段の per-t 出口温度 → T_in への昇温)。旧実装 (各段平均出口を次段へ) からの改修。
#
# ■ 段間再加熱の燃料費が乗る仕組み
#   各段の積分は design.T_in から開始する (= reheat 済)。各段 t の昇温熱量「前段 per-t 出口温度
#   → T_in」を合算・時間平均して Q_preheat に計上 → economics の「Reactor予熱炉燃料」に総量が乗る。
#
# ■ モデルの簡略化 (!仮置き的な割り切り、確定時に精緻化)
#   - 各段は同一ジオメトリ (RadialDesignVars の D_inner/bed_thickness/H) とする (Oleflex の
#     ほぼ同型反応器列に対応)。触媒量・CAPEX は N 段ぶん積算する。
#   - 段間加熱炉の CAPEX は単段同様「反応器バンドルに内包」の現行規約を踏襲し、独立 CAPEX 行は
#     立てない (燃料 OPEX は Q_preheat 合算で正しく乗る)。加熱炉 CAPEX の独立計上は今後の精緻化候補。

def simulate_radial_multibed_reactor_system(
    design: RadialDesignVars,
    feed: FeedStream,
    fixed: FixedParams,
    n_beds: int = 3,
    n_time_samples: int = 20,
) -> SimulationResult:
    """Oleflex 型: 同一ジオメトリの径方向流断熱床を n_beds 段直列 + 段間再加熱 (T_in へ)。

    **時刻同期 (2026-05-31)**: 時刻サンプル t ごとに全段を直列積分 (Bed1(t)→reheat→Bed2(t)→…)
    し、最終段出口を最後に時間平均する。反応の非線形性による f(平均)≠平均(f) のズレを避けるため、
    旧実装 (各段を時間平均してから次段へ渡す) から改修。どこか 1 段・1 時刻でも積分が失敗したら
    sim_failure を伝播。n_beds=1 なら単段ソルバに委譲し完全一致。
    """
    if n_beds < 1:
        n_beds = 1
    if n_beds == 1:
        return simulate_radial_flow_reactor_system(design, feed, fixed, n_time_samples)

    # ---- 入力バリデーション (単段ソルバと同条件) ----
    if (design.t_cyc <= 0 or design.D_inner <= 0 or design.bed_thickness <= 0
            or design.H <= 0):
        return _penalty_result(reason='input_invalid')
    if design.T_in <= 0 or feed.T_feed <= 0 or feed.P_in <= 0:
        return _penalty_result(reason='input_invalid')
    if any(v < 0 for v in feed.F_in.values()) or sum(feed.F_in.values()) <= 0:
        return _penalty_result(reason='input_invalid')

    r_i, r_o, H = design.r_i, design.r_o, design.H

    # ---- ジオメトリ / 並列基数 (1 段あたり。全段同一ジオメトリ。時間ループ前に確定) ----
    V_cat_total = math.pi * (r_o ** 2 - r_i ** 2) * H * (1.0 - fixed.eps)   # 1 段あたり
    N_parallel = max(math.ceil(V_cat_total / fixed.V_cat_max_per_vessel), 1)

    # ---- 時間サンプリング × 段直列 (時刻同期) ----
    if n_time_samples < 2:
        t_samples = np.array([0.0])
    else:
        t_samples = np.linspace(0.0, design.t_cyc, n_time_samples)

    F_final_list, T_final_list, P_final_list, Qpre_list = [], [], [], []
    # 各段の入口条件 (時刻ごと) を SV チェック用に保持
    bed_in_n_mol_s = [[] for _ in range(n_beds)]
    bed_in_P = [[] for _ in range(n_beds)]

    for t in t_samples:
        F_cur = dict(feed.F_in)          # [kmol/h]
        T_from = feed.T_feed             # 当段入口 (予熱前) 温度
        P_cur = feed.P_in
        q_w_t = 0.0                      # 当時刻の全段再加熱熱量 [J/s]
        F_out = None
        for k in range(n_beds):
            n_in_mol_s = sum(F_cur.values()) * 1000.0 / 3600.0
            bed_in_n_mol_s[k].append(n_in_mol_s)
            bed_in_P[k].append(P_cur)
            # 当段の予熱/再加熱熱量 (T_from → design.T_in)
            for comp in _COMPS:
                F_mol_s = F_cur.get(comp, 0.0) * 1000.0 / 3600.0
                q_w_t += F_mol_s * _thermo.calc_enthalpy_change(comp, T_from, design.T_in)
            feed_k = FeedStream(F_in=F_cur, T_feed=T_from, P_in=P_cur)
            F_out, T_out, P_out = _simulate_one_time_radial(
                design, feed_k, fixed, float(t), N_parallel)
            if F_out is None:
                return _penalty_result(reason='sim_failure')
            F_cur = {comp: float(F_out[i]) * 3600.0 / 1000.0 for i, comp in enumerate(_COMPS)}
            T_from = T_out               # 次段入口 = 当段 per-t 出口温度 → 次段で T_in へ再加熱
            P_cur = P_out
        F_final_list.append(np.array([F_cur[c] * 1000.0 / 3600.0 for c in _COMPS]))  # [mol/s]
        T_final_list.append(T_out)
        P_final_list.append(P_out)
        Qpre_list.append(q_w_t)

    # ---- 時間平均 (最終段出口・予熱熱量) ----
    if len(t_samples) >= 2:
        _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
        F_avg_mol_s = _trapz(np.array(F_final_list), t_samples, axis=0) / design.t_cyc
        T_out_avg = float(_trapz(np.array(T_final_list), t_samples) / design.t_cyc)
        P_out_avg = float(_trapz(np.array(P_final_list), t_samples) / design.t_cyc)
        q_w_avg = float(_trapz(np.array(Qpre_list), t_samples) / design.t_cyc)
    else:
        F_avg_mol_s = F_final_list[0]
        T_out_avg, P_out_avg, q_w_avg = T_final_list[0], P_final_list[0], Qpre_list[0]
    F_out_final = {comp: float(F_avg_mol_s[i]) * 3600.0 / 1000.0 for i, comp in enumerate(_COMPS)}
    Q_preheat_total = q_w_avg * 3600.0 / 1e9   # [GJ/h] 供給予熱 + 全段間再加熱

    # ---- 累積 ΔP (床 Ergun × 内部品マージン) のハード制約 ----
    # 全段の真の床ΔP 合計 (時間平均) に dP_margin_factor を掛けて総ΔP とする (化学は真の床P)。
    dP_bed_cum = (feed.P_in - P_out_avg) / feed.P_in if feed.P_in > 0 else 0.0
    dP_over_P_cum = float(np.clip(dP_bed_cum * fixed.dP_margin_factor, 0.0, 1.0))
    P_out_avg = feed.P_in * (1.0 - dP_over_P_cum)   # マージン込み最終出口圧 (→ 下流圧縮機)
    if dP_over_P_cum > fixed.dP_over_P_max:
        warnings.warn(
            f"radial multibed: 累積ΔP/P_in={dP_over_P_cum*100:.1f}% (床×{fixed.dP_margin_factor}) が上限 "
            f"{fixed.dP_over_P_max*100:.0f}% 超 (n_beds={n_beds}) — infeasible 化",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='dP_excess', dP_over_P=dP_over_P_cum)

    # ---- 空塔速度チェック (各段の内側面、案B: N_parallel で割らない)。下流段ほど ----
    #   モル数増 + 圧力減で SV が上がりうるので全段の最大を見て SV_max のみ課す。
    A_inner = 2.0 * math.pi * r_i * H
    SV_inner_max = 0.0
    for k in range(n_beds):
        n_avg = float(np.mean(bed_in_n_mol_s[k]))
        P_avg = float(np.mean(bed_in_P[k]))
        Q_vol_in = n_avg * _R_GAS * design.T_in / P_avg if P_avg > 0 else 0.0
        SV_k = Q_vol_in / A_inner if A_inner > 0 else 0.0
        SV_inner_max = max(SV_inner_max, SV_k)
    if SV_inner_max > fixed.SV_max_m_per_s:
        warnings.warn(
            f"radial multibed: 内側面 SV={SV_inner_max:.2f} m/s が上限 "
            f"{fixed.SV_max_m_per_s} m/s 超 (r_i={r_i:.2f}m, H={H:.1f}m, n_beds={n_beds}) "
            f"— infeasible 化 (フルイダイゼーション/同伴懸念)",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='sv_out_of_range', SV_actual=SV_inner_max)

    # ---- 装置 (1 段あたりを算出し全段積算。全段同一ジオメトリ) ----
    N_swing_sets = math.ceil(fixed.t_regen / design.t_cyc) + 1
    N_reactors_per_bed = N_parallel * N_swing_sets
    V_vessel_actual = (math.pi * (r_o ** 2 - r_i ** 2) * H) / N_parallel   # 1 基あたり (環状部)
    if V_vessel_actual <= 0:
        return _penalty_result(reason='volume_zero')
    catalyst_weight_per_bed = V_cat_total * N_swing_sets * fixed.rho_b
    try:
        reactor_capex_per_bed = calc_reactor_capex_okuyen(
            V_vessel_m3=V_vessel_actual, P_abs_pa=feed.P_in,
            D_m=2.0 * r_o, N_reactors_total=N_reactors_per_bed,
        )
    except Exception:
        reactor_capex_per_bed = _PENALTY_CAPEX

    N_reactors_total = N_reactors_per_bed * n_beds
    catalyst_weight_total = catalyst_weight_per_bed * n_beds
    reactor_capex_total = reactor_capex_per_bed * n_beds

    # ---- 累積転化/選択率 (新規 feed → 最終段流出) ----
    F_A_in = feed.F_in.get('A', 0.0)
    F_A_out = F_out_final.get('A', 0.0)
    F_B_in = feed.F_in.get('B', 0.0)
    F_B_out = F_out_final.get('B', 0.0)
    conversion = (F_A_in - F_A_out) / F_A_in * 100.0 if F_A_in > 0 else 0.0
    delta_A = F_A_in - F_A_out
    selectivity = (F_B_out - F_B_in) / delta_A * 100.0 if delta_A > 0 else 0.0
    conversion = float(np.clip(conversion, 0.0, 100.0))
    selectivity = float(np.clip(selectivity, 0.0, 100.0))

    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg=F_out_final,
            T_out_avg=T_out_avg,
            Q_preheat=Q_preheat_total,        # 全段ぶん → economics の加熱炉燃料に総量が乗る
            P_out=P_out_avg,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=V_vessel_actual,  # 1 基あたり (代表)
            N_parallel=N_parallel,
            N_swing_sets=N_swing_sets,
            N_reactors_total=N_reactors_total,                  # 全段合計基数
            Catalyst_Weight_Total=catalyst_weight_total,        # 全段合計
            Reactor_CAPEX=reactor_capex_total,                  # 全段合計
            dP_over_P_actual=dP_over_P_cum,
        ),
        performance=PerformanceMetrics(
            Conversion=conversion,
            Selectivity=selectivity,
        ),
    )
