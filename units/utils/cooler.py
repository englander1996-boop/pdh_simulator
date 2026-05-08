"""
ダミー冷却器 / 加熱器

物理的な相変化計算は行わず、感熱熱量の粗近似で熱量・伝熱面積・CAPEX を推算する。
出口ストリームは入口と同一組成（気相のまま冷却 / 加熱）。

CAPEX: Bare Module Cost 法（固定管板式熱交換器）を使用。
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
from src.cost_calculator import calc_he_capex_okuyen

# 成分ごとの代表定圧比熱 [J/(mol·K)] (300-600 K 範囲の粗近似)
_CP_APPROX = {
    'A': 100.0,   # C3H8
    'B':  90.0,   # C3H6
    'C':  29.0,   # H2
    'D':  55.0,   # C2H4
    'E':  38.0,   # CH4
    'F':  70.0,   # C2H6
    'Z': 130.0,   # C4H10
}
_CP_DEFAULT = 60.0


@dataclass
class CoolerEquipment:
    """冷却器/加熱器の装置情報"""
    Q_duty_kW: float   # 熱量 [kW]  (負=冷却, 正=加熱)
    A_est_m2:  float   # 推算伝熱面積 [m²]
    CAPEX:     float   # 設備費 [億円]


@dataclass
class CoolerResult:
    """冷却器/加熱器シミュレーション結果"""
    outlet:    ProcessStream
    equipment: CoolerEquipment


def simulate_cooler(
    stream:       ProcessStream,
    T_out_target: float,
    P_out:        float | None = None,
    U_Wm2K:       float = 200.0,
    dT_lm:        float = 30.0,
) -> CoolerResult:
    """冷却器/加熱器をシミュレーションする。

    Parameters
    ----------
    stream       : 入口ストリーム
    T_out_target : 出口目標温度 [K]
    P_out        : 出口圧力 [Pa]。None のとき入口圧力と同じ（圧損なし）
    U_Wm2K       : 総括熱伝達係数 [W/(m²·K)]（デフォルト 200）
    dT_lm        : 対数平均温度差の代替値 [K]（デフォルト 30）
    """
    if P_out is None:
        P_out = stream.P_in

    # 熱量推算 [kW]
    Q_kW = 0.0
    for k, F in stream.F_in.items():
        cp      = _CP_APPROX.get(k, _CP_DEFAULT)
        F_mol_s = F * 1000.0 / 3600.0          # kmol/h → mol/s
        Q_kW   += F_mol_s * cp * (T_out_target - stream.T_in) / 1000.0

    # 伝熱面積推算 [m²]
    A_m2 = max(abs(Q_kW) * 1000.0 / (U_Wm2K * max(dT_lm, 1.0)), 10.0)

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        capex = calc_he_capex_okuyen(A_m2)

    outlet = ProcessStream(F_in=dict(stream.F_in), T_in=T_out_target, P_in=P_out)
    equip  = CoolerEquipment(Q_duty_kW=Q_kW, A_est_m2=A_m2, CAPEX=capex)
    return CoolerResult(outlet=outlet, equipment=equip)
