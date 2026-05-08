"""
ダミー圧縮機

ポリトロピック圧縮式で出口温度と軸動力を推算する。

  W_s  = n/(n-1) × R × T_in × [(P_out/P_in)^((n-1)/n) - 1]  [J/mol]
  T_out = T_in × (P_out/P_in)^((n-1)/n)
  n     = γ / (γ - (γ-1)/η_p)  (ポリトロピック指数)

混合ガスの γ は流量加重平均で近似する。

CAPEX: Bare Module Cost 法（遠心式圧縮機）を使用。
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.cost_calculator import calc_comp_capex_okuyen

_R_GAS = 8.314  # [J/(mol·K)]

# 成分ごとの比熱比 γ = Cp/Cv (粗近似値)
_GAMMA = {
    'A': 1.13,   # C3H8
    'B': 1.15,   # C3H6
    'C': 1.40,   # H2
    'D': 1.18,   # C2H4
    'E': 1.32,   # CH4
    'F': 1.13,   # C2H6
    'Z': 1.10,   # C4H10
}
_GAMMA_DEFAULT = 1.20


@dataclass
class CompressorEquipment:
    """圧縮機の装置情報"""
    W_kW:  float   # 軸動力 [kW]
    T_out: float   # 出口温度 [K]
    CAPEX: float   # 設備費 [億円]


@dataclass
class CompressorResult:
    """圧縮機シミュレーション結果"""
    outlet:    ProcessStream
    equipment: CompressorEquipment


def simulate_compressor(
    stream:       ProcessStream,
    P_out_target: float,
    eta_poly:     float = 0.75,
) -> CompressorResult:
    """圧縮機をシミュレーションする。

    Parameters
    ----------
    stream       : 入口ストリーム
    P_out_target : 出口圧力 [Pa]（入口圧力より大きくなければならない）
    eta_poly     : ポリトロピック効率 [-]（デフォルト 0.75）
    """
    if P_out_target <= stream.P_in:
        raise ValueError(
            f"simulate_compressor: P_out ({P_out_target:.0f} Pa) は "
            f"P_in ({stream.P_in:.0f} Pa) より大きくなければなりません。"
        )

    F_total = sum(stream.F_in.values())   # [kmol/h]
    if F_total <= 0.0:
        raise ValueError("simulate_compressor: 総モル流量がゼロです。")

    # 混合ガスの加重平均 γ
    gamma = sum(
        stream.F_in.get(k, 0.0) * _GAMMA.get(k, _GAMMA_DEFAULT)
        for k in stream.F_in
    ) / F_total

    n     = gamma / (gamma - (gamma - 1.0) / eta_poly)
    ratio = P_out_target / stream.P_in
    exp   = (n - 1.0) / n
    W_mol = n / (n - 1.0) * _R_GAS * stream.T_in * (ratio ** exp - 1.0)  # [J/mol]
    T_out = stream.T_in * ratio ** exp                                      # [K]

    F_mol_s = F_total * 1000.0 / 3600.0   # kmol/h → mol/s
    W_kW    = W_mol * F_mol_s / 1000.0    # W → kW

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        capex = calc_comp_capex_okuyen(W_kW)

    outlet = ProcessStream(F_in=dict(stream.F_in), T_in=T_out, P_in=P_out_target)
    equip  = CompressorEquipment(W_kW=W_kW, T_out=T_out, CAPEX=capex)
    return CompressorResult(outlet=outlet, equipment=equip)
