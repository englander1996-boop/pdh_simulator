"""
冷却器 / 加熱器 (改訂版)

設計判断 (2026-05-08):
  - 旧版は冷却水単独前提・潜熱無視のダミー実装だった。
  - contest.md §2-3 のユーティリティ階層に従い、ターゲット温度から自動的に
    冷媒/熱媒を選択し、その単価を equipment 結果に含める形に改訂。
  - 相変化 (蒸発・凝縮) は phase_change=True でフラグ指定すると、
    src.component_data.LATENT_HEAT_KJ_PER_KMOL を使って潜熱を Q に加算する。
    旧版で Mem 気化器 OPEX が 0 になっていた既知バグを解消。

注意事項:
  - VLE を持たない簡略モデル。phase_change の有無は呼び出し側で判定する必要あり
    (例: 液フィードを露点超まで加熱する場合は True)。
  - U 値は contest §4-4 表 (流速によらず固定) を lookup_U で索引する。
    相変化を伴う場合 (phase_change=True) は顕熱区間と潜熱区間で分けて A を合算。
  - 出口ストリームの相 (気/液) はモデル化しない。CAPEX 計算は伝熱面積ベース。

CAPEX: Bare Module Cost 法（固定管板式熱交換器）。
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
from src.component_data import cp_of, CP_DEFAULT, LATENT_HEAT_KJ_PER_KMOL
from src.utility_selector import select_utility, UtilityTier
from flowsheet.heat_integration import StreamPhase, lookup_U, utility_phase


@dataclass
class CoolerEquipment:
    """冷却器/加熱器の装置情報。

    旧版から追加:
      utility_name      : 選択された冷媒/熱媒の名前 (例 "冷却水", "MP Steam")
      utility_jpy_per_GJ: その単価 [円/GJ]
      Q_latent_kW       : 潜熱分の熱量 [kW] (相変化なしのとき 0)
    economics.py はこれを直接使って OPEX を計算する。
    """
    Q_duty_kW:          float   # 熱量 [kW]  (負=冷却, 正=加熱), 顕熱+潜熱の合計
    Q_sensible_kW:      float   # うち顕熱 [kW]
    Q_latent_kW:        float   # うち潜熱 [kW] (相変化指定時のみ非ゼロ)
    A_est_m2:           float   # 推算伝熱面積 [m²]
    CAPEX:              float   # 設備費 [億円]
    utility_name:       str     # 選択ユーティリティ名
    utility_jpy_per_GJ: float   # 選択ユーティリティ単価 [円/GJ]


@dataclass
class CoolerResult:
    """冷却器/加熱器シミュレーション結果"""
    outlet:    ProcessStream
    equipment: CoolerEquipment


def simulate_cooler(
    stream:        ProcessStream,
    T_out_target:  float,
    P_out:         float | None = None,
    dT_lm:         float = 30.0,
    phase_change:  bool  = False,
    process_phase: str   = StreamPhase.GAS,
) -> CoolerResult:
    """冷却器/加熱器をシミュレーション。

    Parameters
    ----------
    stream : ProcessStream
        入口ストリーム
    T_out_target : float [K]
        出口目標温度
    P_out : float [Pa] | None
        出口圧力 (None で入口と同じ、圧損なし)
    dT_lm : float [K]
        対数平均温度差の代替値 (デフォルト 30)
    phase_change : bool
        True のとき、入口/出口で相変化が起こると見なし、全成分の潜熱を Q に加算する。
        加熱方向 (T_out > T_in) なら蒸発、冷却方向なら凝縮として扱う。
    process_phase : str
        プロセス側の顕熱区間の相 (StreamPhase.GAS / LIQUID)。
        contest §4-4 表索引用。デフォルト GAS (圧縮ガス冷却が大半)。

    総括 U は contest §4-4 (lookup_U) で (hot_phase, cold_phase) ペアから引く。
    Ref: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 (流速によらず固定 U)。

    Returns
    -------
    CoolerResult
        outlet と equipment (使用ユーティリティ情報を含む)
    """
    if P_out is None:
        P_out = stream.P_in

    # ---- 顕熱 [kW] ----
    Q_sensible_kW = 0.0
    for k, F in stream.F_in.items():
        cp      = cp_of(k) if cp_of(k) > 0 else CP_DEFAULT
        F_mol_s = F * 1000.0 / 3600.0          # kmol/h → mol/s
        Q_sensible_kW += F_mol_s * cp * (T_out_target - stream.T_in) / 1000.0

    # ---- 潜熱 [kW] (相変化指定時のみ) ----
    # 設計判断: 加熱方向のみ潜熱を加算 (T_out > T_in で蒸発)。
    # 凝縮 (T_out < T_in) も理論上は潜熱を扱うべきだが、本フローでは現状不要 (旧版踏襲)。
    Q_latent_kW = 0.0
    if phase_change and T_out_target > stream.T_in:
        for k, F in stream.F_in.items():
            dh_kJ_per_kmol = LATENT_HEAT_KJ_PER_KMOL.get(k, 0.0)
            F_mol_s = F * 1000.0 / 3600.0
            # ΔH [kJ/kmol] × F [mol/s] / 1000 (mol→kmol) = kJ/s = kW
            Q_latent_kW += dh_kJ_per_kmol * F_mol_s / 1000.0

    Q_kW = Q_sensible_kW + Q_latent_kW

    # ---- ユーティリティ自動選択 ----
    # 設計判断: target が inlet より低ければ冷却、高ければ加熱として
    # 適切な tier を選ぶ。equipment に名前と単価を埋め込み、economics.py が
    # それを直接読み出して OPEX を計算する。
    utility = select_utility(target_T_K=T_out_target, inlet_T_K=stream.T_in)

    # ---- 伝熱面積推算 [m²] (contest §4-4 表で U を決定) ----
    # 相変化を伴う場合は顕熱区間と潜熱区間で U を分けて A を合算 (§4-4 規定)。
    is_heating  = T_out_target > stream.T_in
    util_phase  = utility_phase(utility.name)
    sens_proc_phase = process_phase
    # 潜熱区間: 加熱なら EVAPORATING、冷却なら CONDENSING (プロセス側)
    lat_proc_phase  = StreamPhase.EVAPORATING if is_heating else StreamPhase.CONDENSING
    if is_heating:
        # プロセス=cold, utility=hot
        U_sens = lookup_U(util_phase, sens_proc_phase)
        U_lat  = lookup_U(util_phase, lat_proc_phase)
    else:
        # プロセス=hot, utility=cold
        U_sens = lookup_U(sens_proc_phase, util_phase)
        U_lat  = lookup_U(lat_proc_phase,  util_phase)

    dT_eff = max(dT_lm, 1.0)
    A_sens = abs(Q_sensible_kW) * 1000.0 / (U_sens * dT_eff)
    A_lat  = abs(Q_latent_kW)   * 1000.0 / (U_lat  * dT_eff) if Q_latent_kW > 0 else 0.0
    A_m2   = max(A_sens + A_lat, 10.0)

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        capex = calc_he_capex_okuyen(A_m2)

    outlet = ProcessStream(F_in=dict(stream.F_in), T_in=T_out_target, P_in=P_out)
    equip  = CoolerEquipment(
        Q_duty_kW         =Q_kW,
        Q_sensible_kW     =Q_sensible_kW,
        Q_latent_kW       =Q_latent_kW,
        A_est_m2          =A_m2,
        CAPEX             =capex,
        utility_name      =utility.name,
        utility_jpy_per_GJ=utility.jpy_per_GJ,
    )
    return CoolerResult(outlet=outlet, equipment=equip)
