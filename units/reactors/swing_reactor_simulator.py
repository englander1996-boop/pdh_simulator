"""
PDH スイング反応器システム シミュレーター

触媒失活を伴うプロパン脱水素断熱PFR（スイング操作）を模擬し、
後段分離工程への時間平均ストリームと装置コスト情報を出力する。

使用方法:
    from units.reactors.swing_reactor_simulator import (
        DesignVars, FeedStream, FixedParams, simulate_swing_reactor_system
    )
    result = simulate_swing_reactor_system(design, feed, fixed)
"""

import math
import os
import sys
from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.config import THERMO_DATA
from src.thermo import PDHThermo
from src.catalyst_model import calculate_activity_a

_thermo = PDHThermo()

# 成分順序（状態ベクトルのインデックスと対応）
_COMPS = ['C3H8', 'C3H6', 'H2', 'C2H4', 'CH4', 'C2H6']
_COMP_KEYS = ['A', 'B', 'C', 'D', 'E', 'F']  # THERMO_DATA のキー

# 化学量論行列 stoich[i, j] = 成分i に対する反応j の量論係数
#        r1   r2   r3
_STOICH = np.array([
    [-1,  -1,   0],  # C3H8
    [+1,   0,   0],  # C3H6
    [+1,   0,  -1],  # H2
    [ 0,  +1,  -1],  # C2H4
    [ 0,  +1,   0],  # CH4
    [ 0,   0,  +1],  # C2H6
], dtype=float)

_T_REF = 298.15  # [K] エンタルピー計算基準温度


# ---------------------------------------------------------------------------
# グローバル関数（ODE 内から呼び出す、パラメータ引き回し不要）
# ---------------------------------------------------------------------------

def calc_a(t: float, T: float, P: float) -> float:
    """触媒活性度 a ∈ [0, 1]（t [min], T [K], P [Pa]）"""
    return calculate_activity_a(T - 273.15, t)


def calc_Cp(T: float) -> dict:
    """各成分のモル熱容量 [J/(mol·K)]。keys: 'C3H8','C3H6','H2','C2H4','CH4','C2H6'"""
    return {
        'C3H8': _thermo.calc_cp('A', T),
        'C3H6': _thermo.calc_cp('B', T),
        'H2':   _thermo.calc_cp('C', T),
        'C2H4': _thermo.calc_cp('D', T),
        'CH4':  _thermo.calc_cp('E', T),
        'C2H6': _thermo.calc_cp('F', T),
    }


# ---------------------------------------------------------------------------
# 反応速度定数・平衡定数（ダミー固定値）
# ---------------------------------------------------------------------------
# 外部物性計算関数が用意できたら calc_rate_constants の中身を差し替えること。
# 差し替え手順:
#   1. 外部関数 get_rate_constants(T) を実装
#   2. 下記の固定値辞書を `return get_rate_constants(T)` に置き換える
#   3. ダミー定数 _DUMMY_RC は削除して構わない
_DUMMY_RC: dict = {
    'k1':   9.787e-5,  # [mol m⁻³ s⁻¹ Pa⁻¹]  T=793K 基準値（config.py k01）
    'k2':   8.682e-7,  # [mol m⁻³ s⁻¹ Pa⁻¹]  T=793K 基準値（config.py k02）
    'k3':   4.406e-8,  # [mol m⁻³ s⁻¹ Pa⁻²]  T=793K 基準値（config.py k03）
    'K_B':  3.46e5,    # [Pa]                  T=793K 基準値（config.py K0）
    'K_eq': 2.0e4,     # [Pa]                  概算固定値（要後日実装）
}


def calc_rate_constants(T: float) -> dict:
    """反応速度定数・平衡定数の辞書（ダミー固定値）。差し替え対象。"""
    return _DUMMY_RC


# ---------------------------------------------------------------------------
# データクラス定義
# ---------------------------------------------------------------------------

@dataclass
class DesignVars:
    """最適化アルゴリズムが操作する設計変数"""
    T_in:  float   # 反応器入口温度 [K]
    z_cat: float   # 触媒層長さ [m]
    t_cyc: float   # 1サイクル反応フェーズ運転時間 [min]
    D:     float   # 反応器直径 [m]


@dataclass
class FeedStream:
    """入口流体条件"""
    F_in:   Dict[str, float]  # 各成分入口モル流量 [kmol/h]
                               # keys: 'C3H8','C3H6','H2','C2H4','CH4','C2H6'
    T_feed: float              # 加熱炉入口（予熱前）原料温度 [K]
    P_in:   float              # 反応器入口圧力 [Pa]


@dataclass
class FixedParams:
    """固定定数・制約条件"""
    t_regen:              float = 30.0    # 触媒再生時間 [min]
    V_cat_max_per_vessel: float = 200.0   # 1基最大触媒量 [m³]
    eps:                  float = 0.5     # 空隙率 [-]
    rho_p:                float = 400.0   # 触媒充填密度 [kg/m³]


@dataclass
class EffluentStream:
    """後段分離工程への出口流体情報"""
    F_out_avg: Dict[str, float]  # 各成分出口モル流量の時間平均 [kmol/h]
    T_out_avg: float              # 出口温度の時間平均 [K]
    Q_preheat: float              # T_feed → T_in 予熱熱量 [GJ/h]
    P_out:     float              # 出口圧力 [Pa]


@dataclass
class EquipmentCost:
    """装置・経済性情報"""
    V_vessel_actual:       float  # 1基プロセス容器容積 [m³]
    N_parallel:            int    # 200m³制約による並列基数
    N_swing_sets:          int    # 再生をカバーする切り替えセット数
    N_reactors_total:      int    # 総反応器基数
    Catalyst_Weight_Total: float  # システム全体触媒総量 [kg]
    Reactor_CAPEX:         float  # 全基分建設コスト合計 [億円]（未実装: nan）


@dataclass
class PerformanceMetrics:
    """プロセス指標"""
    Conversion:  float  # プロパン単通反応率（時間平均）[%]
    Selectivity: float  # プロピレン選択率（時間平均）[%]


@dataclass
class SimulationResult:
    """シミュレーション結果（全出力をまとめるルートオブジェクト）"""
    effluent:    EffluentStream
    equipment:   EquipmentCost
    performance: PerformanceMetrics


# ---------------------------------------------------------------------------
# 内部ヘルパー
# ---------------------------------------------------------------------------

def _reaction_enthalpies(T: float) -> np.ndarray:
    """温度T [K] での各反応エンタルピー ΔH_rxn [J/mol]。shape (3,)"""
    H = {}
    for comp, key in zip(_COMPS, _COMP_KEYS):
        H[comp] = (THERMO_DATA[key].dHf_298
                   + _thermo.calc_enthalpy_change(key, _T_REF, T))

    dH1 = H['C3H6'] + H['H2']  - H['C3H8']  # C3H8 → C3H6 + H2
    dH2 = H['C2H4'] + H['CH4'] - H['C3H8']  # C3H8 → C2H4 + CH4
    dH3 = H['C2H6'] - H['C2H4'] - H['H2']   # C2H4 + H2 → C2H6
    return np.array([dH1, dH2, dH3])


def _ode_axial(z: float, y: np.ndarray,
               t_min: float, A_cross: float, eps: float, P_in: float
               ) -> np.ndarray:
    """軸方向(z)の常微分方程式。

    State y = [F_C3H8, F_C3H6, F_H2, F_C2H4, F_CH4, F_C2H6, T]
    単位: F [mol/s], T [K]
    """
    F = y[:6]
    T_local = float(y[6])

    F_total = float(np.sum(F))
    if F_total <= 0:
        return np.zeros(7)

    # 局所分圧 [Pa]（圧力損失なしのため P_total = P_in）
    P_local = P_in
    P = {comp: max(float(F[i]) / F_total * P_local, 0.0)
         for i, comp in enumerate(_COMPS)}

    # 触媒活性（局所温度・局所圧力で毎ステップ評価）
    a = calc_a(t_min, T_local, P_local)

    # 反応速度定数
    rc = calc_rate_constants(T_local)
    k1, k2, k3 = rc['k1'], rc['k2'], rc['k3']
    K_B, K_eq  = rc['K_B'], rc['K_eq']

    # 反応速度 [mol/m³_cat/s]
    driving_r1 = P['C3H8'] - P['C3H6'] * P['H2'] / K_eq
    r1 = a * k1 * driving_r1 / (1.0 + P['C3H6'] / K_B)
    r2 = k2 * P['C3H8']
    r3 = k3 * P['C2H4'] * P['H2']

    rates = np.array([r1, r2, r3])  # [mol/m³_cat/s]

    # 物質収支: dF_i/dz = eps * A * Σ stoich[i,j] * r[j]  [mol/(s·m)]
    dFdz = eps * A_cross * (_STOICH @ rates)

    # エネルギー収支（断熱）: dT/dz = -((1-eps) * A * Σ ΔH_j * r_j) / Σ F_i * Cp_i
    cp_dict = calc_Cp(T_local)
    sum_FCp = sum(float(F[i]) * cp_dict[comp] for i, comp in enumerate(_COMPS))
    if sum_FCp <= 0:
        dTdz = 0.0
    else:
        dH = _reaction_enthalpies(T_local)
        Q_rxn = float((1.0 - eps) * A_cross * np.dot(dH, rates))  # [J/(m·s)]
        dTdz = -Q_rxn / sum_FCp                                     # [K/m]

    return np.concatenate([dFdz, [dTdz]])


def _simulate_one_time(design: DesignVars, feed: FeedStream,
                       fixed: FixedParams, t_min: float) -> tuple:
    """時刻 t_min [min] における空間積分を実施し (F_out [mol/s], T_out [K]) を返す。"""
    A_cross = math.pi / 4.0 * design.D ** 2  # [m²]

    F0 = np.array([feed.F_in.get(c, 0.0) * 1000.0 / 3600.0 for c in _COMPS])
    y0 = np.concatenate([F0, [design.T_in]])

    sol = solve_ivp(
        fun=lambda z, y: _ode_axial(z, y, t_min, A_cross, fixed.eps, feed.P_in),
        t_span=(0.0, design.z_cat),
        y0=y0,
        method='RK45',
        rtol=1e-5,
        atol=1e-8,
    )

    return sol.y[:6, -1], float(sol.y[6, -1])


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_swing_reactor_system(
    design: DesignVars,
    feed: FeedStream,
    fixed: FixedParams,
    n_time_samples: int = 20,
) -> SimulationResult:
    """スイング反応器システムをシミュレーションする。

    Parameters
    ----------
    design         : DesignVars   設計変数
    feed           : FeedStream   入口流体条件
    fixed          : FixedParams  固定パラメータ
    n_time_samples : int          時間軸サンプリング点数（デフォルト20）
    """
    # ---- 時間方向サンプリングと空間積分 ----
    t_samples = np.linspace(0.0, design.t_cyc, n_time_samples)
    F_out_list, T_out_list = [], []

    for t in t_samples:
        F_out, T_out = _simulate_one_time(design, feed, fixed, float(t))
        F_out_list.append(F_out)
        T_out_list.append(T_out)

    F_out_arr = np.array(F_out_list)  # (n_time_samples, 6) [mol/s]
    T_out_arr = np.array(T_out_list)  # (n_time_samples,)   [K]

    # 時間平均（台形則）
    _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
    F_out_avg_mol_s = _trapz(F_out_arr, t_samples, axis=0) / design.t_cyc
    T_out_avg = float(_trapz(T_out_arr, t_samples) / design.t_cyc)

    # mol/s → kmol/h
    F_out_avg_kmolh = {comp: float(F_out_avg_mol_s[i]) * 3600.0 / 1000.0
                       for i, comp in enumerate(_COMPS)}

    # ---- 予熱熱量 Q_preheat [GJ/h] ----
    q_w = 0.0
    for comp, key in zip(_COMPS, _COMP_KEYS):
        F_mol_s = feed.F_in.get(comp, 0.0) * 1000.0 / 3600.0
        q_w += F_mol_s * _thermo.calc_enthalpy_change(key, feed.T_feed, design.T_in)
    Q_preheat_GJh = q_w * 3600.0 / 1e9

    # ---- 装置計算 ----
    A_cross = math.pi / 4.0 * design.D ** 2
    V_cat_total = A_cross * design.z_cat * fixed.eps
    N_parallel = max(math.ceil(V_cat_total / fixed.V_cat_max_per_vessel), 1)
    N_swing_sets = math.ceil(fixed.t_regen / design.t_cyc) + 1
    N_reactors_total = N_parallel * N_swing_sets

    V_vessel_actual = (V_cat_total / N_parallel) / fixed.eps  # [m³]
    catalyst_weight_total = V_cat_total * N_swing_sets * fixed.rho_p  # [kg]

    # CAPEX: コスト推算モジュール未実装のためダミー値
    # 実装後は外部コスト関数の呼び出しに差し替えること。
    reactor_capex = float('nan')

    # ---- パフォーマンス指標 ----
    F_A_in  = feed.F_in.get('C3H8', 0.0)
    F_A_out = F_out_avg_kmolh['C3H8']
    F_B_in  = feed.F_in.get('C3H6', 0.0)
    F_B_out = F_out_avg_kmolh['C3H6']

    conversion  = (F_A_in - F_A_out) / F_A_in * 100.0 if F_A_in > 0 else 0.0
    delta_A     = F_A_in - F_A_out
    selectivity = (F_B_out - F_B_in) / delta_A * 100.0 if delta_A > 0 else 0.0

    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg=F_out_avg_kmolh,
            T_out_avg=T_out_avg,
            Q_preheat=Q_preheat_GJh,
            P_out=feed.P_in,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=V_vessel_actual,
            N_parallel=N_parallel,
            N_swing_sets=N_swing_sets,
            N_reactors_total=N_reactors_total,
            Catalyst_Weight_Total=catalyst_weight_total,
            Reactor_CAPEX=reactor_capex,
        ),
        performance=PerformanceMetrics(
            Conversion=conversion,
            Selectivity=selectivity,
        ),
    )
