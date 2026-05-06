"""
全ユニット共通の入出力ストリーム型

ProcessStream は個別ユニット (Swing, PSA, Membrane, Distillation, Mixer, Cooler,
Compressor) の FeedStream / OutStream 間のブリッジとして使用する。
"""

from dataclasses import dataclass
from typing import Dict


@dataclass
class ProcessStream:
    """汎用プロセスストリーム。

    Attributes
    ----------
    F_in : Dict[str, float]
        各成分モル流量 [kmol/h]。
        keys: 'A'(C3H8), 'B'(C3H6), 'C'(H2), 'D'(C2H4),
               'E'(CH4), 'F'(C2H6), 'Z'(C4H10)
    T_in : float
        温度 [K]
    P_in : float
        圧力 [Pa]
    """
    F_in: Dict[str, float]
    T_in: float
    P_in: float

    def total_flow(self) -> float:
        """総モル流量 [kmol/h]"""
        return sum(self.F_in.values())
