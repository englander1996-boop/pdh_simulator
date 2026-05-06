"""
ダミー蒸留塔 共通エンジン

VLE や厳密な物質・エネルギー収支計算は行わない。
分配率(split_fracs)で成分を塔頂/塔底に振り分け、
簡易近似で温度・熱量・塔寸法・CAPEX を推算する。

各塔固有の設定は fake_column1/2/3.py で定義する。
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from units.utils.process_stream import ProcessStream
from src.cost_calculator import calc_cp0, calc_fp
from src.cost_parameters import (
    B1, B2, FM, CEPCI_BASE, CEPCI_CURRENT,
    USD_TO_JPY, PLANT_INDIRECT_FACTOR,
)

_R_GAS = 8.314     # [J/(mol·K)]
_P_ATM = 101325.0  # [Pa]

# 各成分の大気圧沸点 [K] と蒸発潜熱 [kJ/mol]
_T_BOIL_ATM = {
    'A': 231.1, 'B': 225.5, 'C':  20.3,
    'D': 169.4, 'E': 111.7, 'F': 184.6, 'Z': 272.7,
}
_LAMBDA_KJ = {
    'A': 19.1, 'B': 18.4, 'C': 0.9,
    'D': 13.5, 'E':  8.2, 'F': 14.7, 'Z': 22.0,
}
_T_BOIL_DEFAULT = 250.0
_LAMBDA_DEFAULT  = 15.0


@dataclass
class DistDesignVars:
    """蒸留塔設計変数"""
    P_col:        float             # 塔操作圧力 [Pa]
    N_stages:     int               # 理論段数 [-]
    reflux_ratio: float             # 還流比 L/D [-]
    split_fracs:  Dict[str, float]  # 各成分の塔頂回収率 [0〜1]


@dataclass
class DistFixedParams:
    """蒸留塔固定パラメータ"""
    tray_spacing_m: float = 0.6   # トレイ間隔 [m]
    sump_height_m:  float = 3.0   # 缶部・塔頂部付加長さ [m]
    u_vapor_ms:     float = 0.3   # 代表蒸気線速度 [m/s]（塔径推算用）


@dataclass
class DistEquipment:
    """蒸留塔の装置情報"""
    D_col:  float   # 塔径 [m]
    H_col:  float   # 塔高さ [m]
    CAPEX:  float   # 設備費 [億円]
    Q_cond: float   # コンデンサー熱量 [kW]（系外へ放出、正値）
    Q_reb:  float   # リボイラー熱量 [kW]（系外から受取、正値）


@dataclass
class DistResult:
    """蒸留塔シミュレーション結果"""
    top:       ProcessStream
    bottom:    ProcessStream
    equipment: DistEquipment


_T_COND_MIN = 313.15  # 40°C — 冷却水使用時の凝縮器最低温度


def _boil_cc(T1: float, lam_kJ: float, P_col: float) -> float:
    """Clausius-Clapeyron 式で圧力 P_col [Pa] における沸点 [K] を返す。

    ln(P2/P1) = ΔHvap/R × (1/T1 - 1/T2)
    → T2 = 1 / (1/T1 - R/ΔHvap × ln(P2/P1))
    """
    lam_J = lam_kJ * 1000.0
    inv_T2 = 1.0 / T1 - _R_GAS / lam_J * math.log(P_col / _P_ATM)
    if inv_T2 <= 0.0:
        return T1 * 5.0   # 超臨界域のガード値
    return 1.0 / inv_T2


def _weighted_boil(F_dict: dict, P_col: float) -> float:
    """流量加重平均の沸点 [K] を Clausius-Clapeyron 式で圧力補正する。"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return 298.15
    t_boil = sum(
        F_dict.get(k, 0.0) * _boil_cc(
            _T_BOIL_ATM.get(k, _T_BOIL_DEFAULT),
            _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT),
            P_col,
        )
        for k in F_dict
    ) / F_total
    return t_boil


def _weighted_lambda(F_dict: dict) -> float:
    """流量加重平均蒸発潜熱 [kJ/mol]"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return _LAMBDA_DEFAULT
    return sum(
        F_dict.get(k, 0.0) * _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT)
        for k in F_dict
    ) / F_total


def _column_capex_okuyen(V_m3: float, P_pa: float, D_m: float) -> float:
    """蒸留塔（縦型容器）の総建設費 [億円]。

    反応器用 calc_reactor_capex_okuyen と異なり K_SWING を適用しない。
    """
    P_bar = P_pa / 1.0e5 - 1.01325
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        cp0 = calc_cp0(V_m3)
    fp  = calc_fp(P_bar, D_m)
    fbm = B1 + B2 * fp * FM
    cbm = cp0 * fbm
    usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * usd * USD_TO_JPY / 1.0e8


def simulate_distillation_column(
    design: DistDesignVars,
    feed:   ProcessStream,
    fixed:  DistFixedParams | None = None,
) -> DistResult:
    """ダミー蒸留塔をシミュレーションする。

    Parameters
    ----------
    design : DistDesignVars  設計変数
    feed   : ProcessStream   入口ストリーム
    fixed  : DistFixedParams（None のときデフォルト値を使用）
    """
    if fixed is None:
        fixed = DistFixedParams()

    # ---- 物質収支 ----
    F_top:    Dict[str, float] = {}
    F_bottom: Dict[str, float] = {}
    for k, F in feed.F_in.items():
        alpha        = max(0.0, min(1.0, design.split_fracs.get(k, 0.5)))
        F_top[k]    = F * alpha
        F_bottom[k] = F * (1.0 - alpha)

    F_top_total    = sum(F_top.values())
    F_bottom_total = sum(F_bottom.values())

    # ---- 温度推算 ----
    # H2 等の非凝縮ガスが多い場合は冷却水温度下限で制約
    T_cond = max(_weighted_boil(F_top,    design.P_col) + 5.0, _T_COND_MIN)
    T_reb  = _weighted_boil(F_bottom, design.P_col) + 5.0

    # ---- 熱量推算 ----
    lam_top   = _weighted_lambda(F_top)
    Q_cond_kW = (F_top_total * (design.reflux_ratio + 1.0)
                 * lam_top * 1000.0 / 3600.0)    # [kW]
    Q_reb_kW  = Q_cond_kW * 1.05                  # ≈ 5% 損失分を加算

    # ---- 塔径推算 ----
    # V̇_vap [m³/s] = F_vap_mol_s × R × T_avg / P_col
    T_avg       = (T_cond + T_reb) / 2.0
    F_vap_mol_s = F_top_total * (design.reflux_ratio + 1.0) * 1000.0 / 3600.0
    Q_vap = F_vap_mol_s * _R_GAS * T_avg / design.P_col  # [m³/s]
    D_col = math.sqrt(4.0 * Q_vap / (math.pi * fixed.u_vapor_ms))
    D_col = max(D_col, 0.3)    # 最小直径 0.3 m

    # ---- 塔高さ推算 ----
    H_col = design.N_stages * fixed.tray_spacing_m + fixed.sump_height_m
    H_col = max(H_col, 5.0)    # 最小高さ 5 m

    # ---- CAPEX ----
    V_col = math.pi / 4.0 * D_col ** 2 * H_col
    capex = _column_capex_okuyen(V_col, design.P_col, D_col)

    return DistResult(
        top    = ProcessStream(F_in=F_top,    T_in=T_cond, P_in=design.P_col),
        bottom = ProcessStream(F_in=F_bottom, T_in=T_reb,  P_in=design.P_col),
        equipment=DistEquipment(
            D_col=D_col,
            H_col=H_col,
            CAPEX=capex,
            Q_cond=Q_cond_kW,
            Q_reb=Q_reb_kW,
        ),
    )
