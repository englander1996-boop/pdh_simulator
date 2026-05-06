"""
ストリームミキサー

複数の ProcessStream を1本に合流させる。
  物質収支 : モル流量を成分ごとに加算
  エネルギー収支 : 総モル流量加重平均温度（簡略近似）
  圧力     : 最低圧力を採用（高圧側は制御弁でドロップ）
"""

import os
import sys
from typing import List

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from units.utils.process_stream import ProcessStream


def mix_streams(streams: List[ProcessStream]) -> ProcessStream:
    """複数のストリームを混合する。

    Parameters
    ----------
    streams : list[ProcessStream]
        混合するストリームのリスト（1本以上）

    Returns
    -------
    ProcessStream
        混合後ストリーム
    """
    if not streams:
        raise ValueError("mix_streams: 1本以上のストリームが必要です。")

    F_out: dict = {}
    H_num  = 0.0   # Σ(F_i × T_i) [kmol/h × K]
    F_total = 0.0
    P_out  = min(s.P_in for s in streams)

    for s in streams:
        for k, v in s.F_in.items():
            F_out[k] = F_out.get(k, 0.0) + v
        f_i     = sum(s.F_in.values())
        H_num  += f_i * s.T_in
        F_total += f_i

    T_out = H_num / F_total if F_total > 0.0 else 298.15
    return ProcessStream(F_in=F_out, T_in=T_out, P_in=P_out)
