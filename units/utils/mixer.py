"""
ストリームミキサー

  成分別 Cp を考慮したエンタルピーバランスで T_out を求める。単純なモル流量
  加重平均では、Cp が異なる成分が混ざる場合や温度差が大きい場合 (例: 反応器
  入口で 100°C のリサイクルと 600°C が混ざる場合) に物理的に不正確になる。

エンタルピーバランス (基準温度 T_ref を導入して両辺で消える):
  Σ_i (F_i × Cp_mix_i × (T_i - T_ref)) = (Σ F_out × Cp_mix_out) × (T_out - T_ref)
  ⇔ T_out = Σ_i (F_i × Cp_mix_i × T_i) / Σ_i (F_i × Cp_mix_i)
  Cp_mix_i は ストリーム i のモル流量加重平均 Cp。

  これにより各ストリームの「持っているエンタルピー総量」がモル流量だけでなく
  成分・Cp に応じて適切に重み付けされる。

  圧力: 最低圧力を採用 (高圧側は制御弁でドロップ)。
"""

import os
import sys
from typing import List

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.component_data import cp_of, CP_DEFAULT


def mix_streams(streams: List[ProcessStream]) -> ProcessStream:
    """複数のストリームをエンタルピーバランスで合流させる。

    Parameters
    ----------
    streams : list[ProcessStream]
        混合するストリームのリスト (1本以上)

    Returns
    -------
    ProcessStream
        混合後ストリーム (組成: 成分別加算、温度: エンタルピー保存)
    """
    if not streams:
        raise ValueError("mix_streams: 1本以上のストリームが必要です。")

    F_out: dict = {}
    P_out: float = min(s.P_in for s in streams)

    # 成分別の流量加算
    for s in streams:
        for k, v in s.F_in.items():
            F_out[k] = F_out.get(k, 0.0) + v

    # ---- T_out: エンタルピーバランス ----
    # H_total = Σ_i (Σ_k F_ik × Cp_k × T_i) と書ける (T_ref を引いて消す前提)
    # T_out = H_total / Σ_k (F_out_k × Cp_k)
    H_num = 0.0    # Σ Σ F_ik × Cp_k × T_i
    for s in streams:
        for k, v in s.F_in.items():
            H_num += v * cp_of(k) * s.T_in

    # 出口側の F × Cp 合計 (T_out の coefficient)
    H_denom = sum(F_out.get(k, 0.0) * cp_of(k) for k in F_out)

    if H_denom <= 0.0:
        # 全流量ゼロの異常ケース。298.15K (基準温度) を返すが、下流の cooler/
        # compressor が「P_in=0 流量=0」状態を受け取って ValueError を投げる場合
        # がある。silent fallback だと原因追跡が難しいため warning を出す。
        # flowsheet/run_one_pass.py の _capture_warnings 経由で BO log に残る。
        import warnings as _warnings
        _warnings.warn(
            f"mix_streams: 全流量ゼロ (H_denom={H_denom:.3e})。"
            f"T_out=298.15K (基準温度) で fallback。"
            f"入力ストリーム数={len(streams)}、上流装置の penalty 状態の可能性。",
            UserWarning, stacklevel=2,
        )
        T_out = 298.15
    else:
        T_out = H_num / H_denom

    return ProcessStream(F_in=F_out, T_in=T_out, P_in=P_out)
