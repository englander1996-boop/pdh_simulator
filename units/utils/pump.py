"""
液体ポンプ (Centrifugal pump)

液→液の昇圧。所要動力は流体動力 / 効率で算出する:
  W_kW = V_dot [m³/s] × dP [Pa] / η_pump / 1000

非圧縮性流体仮定により、出口温度は入口とほぼ等しい (粘性発熱は無視)。
contest §4-5-1 に従い、コンプレッサ (断熱効率 0.80) と区別して扱う。

CAPEX: Bare Module Cost 法 (遠心式ポンプ、炭素鋼)。
"""

import os
import sys
import warnings
from dataclasses import dataclass

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.component_data import MW, liquid_density_mix
from src.cost_calculator import calc_pump_capex_okuyen


@dataclass
class PumpEquipment:
    """ポンプの装置情報"""
    W_kW:     float   # 軸動力 [kW] (流体動力 / 効率)
    rho_liq:  float   # 入口液密度 [kg/m³] (検算用)
    V_dot:    float   # 体積流量 [m³/s] (検算用)
    CAPEX:    float   # 設備費 [億円]


@dataclass
class PumpResult:
    """ポンプシミュレーション結果"""
    outlet:    ProcessStream
    equipment: PumpEquipment


def simulate_pump(
    stream:       ProcessStream,
    P_out_target: float,
    eta_pump:     float = 0.70,
) -> PumpResult:
    """液体ポンプをシミュレーションする。

    Parameters
    ----------
    stream       : 入口ストリーム (液相想定)
    P_out_target : 出口圧力 [Pa] (入口より大きいこと)
    eta_pump     : ポンプ効率 [-] (デフォルト 0.70)
                   Ref: 化工便覧 改訂六版 5·6·4 項【例題 5·8】(η=0.7 で軸馬力算出)
    """
    if P_out_target <= stream.P_in:
        raise ValueError(
            f"simulate_pump: P_out ({P_out_target:.0f} Pa) は "
            f"P_in ({stream.P_in:.0f} Pa) より大きくなければなりません。"
        )

    F_total = sum(stream.F_in.values())   # [kmol/h]
    if F_total <= 0.0:
        raise ValueError("simulate_pump: 総モル流量がゼロです。")

    # 質量流量・体積流量
    mass_flow_kg_h = sum(stream.F_in.get(k, 0.0) * MW.get(k, 50.0) for k in stream.F_in)
    rho_liq = liquid_density_mix(stream.F_in)               # [kg/m³]
    V_dot   = mass_flow_kg_h / rho_liq / 3600.0             # [m³/s]

    dP    = P_out_target - stream.P_in                      # [Pa]
    W_kW  = V_dot * dP / eta_pump / 1000.0                  # [kW]

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        capex = calc_pump_capex_okuyen(W_kW, P_out_target)

    # 非圧縮性流体: 出口温度はほぼ不変 (粘性発熱無視)
    outlet = ProcessStream(F_in=dict(stream.F_in), T_in=stream.T_in, P_in=P_out_target)
    equip  = PumpEquipment(W_kW=W_kW, rho_liq=rho_liq, V_dot=V_dot, CAPEX=capex)
    return PumpResult(outlet=outlet, equipment=equip)
