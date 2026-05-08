"""
最適化対象の設計変数バンドル。

最適化器は FlowsheetDesignVars インスタンスを生成して runner に渡す。
設計変数の追加 (たとえば膜入口温度を最適化するなど) はこの dataclass を
拡張する形で行う。
"""

from dataclasses import dataclass

from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars


@dataclass(frozen=True)
class FlowsheetDesignVars:
    """フローシート全体の設計変数。最適化器はこれを生成して runner へ渡す。"""
    swing: SwingDesign
    psa:   PSADesignVars
    mem:   MemDesignVars
