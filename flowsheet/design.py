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
from src.distillation_core import ColumnTunables


@dataclass(frozen=True)
class FlowsheetDesignVars:
    """フローシート全体の設計変数。最適化器はこれを生成して runner へ渡す。

    蒸留塔の物理セマンティクス (LK/HK/回収率/K_method/q) は塔別ラッパで固定し、
    本バンドルでは BO で振る P/N/R の組のみを保持する (ColumnTunables)。
    """
    swing: SwingDesign
    psa:   PSADesignVars
    mem:   MemDesignVars
    dist1: ColumnTunables       # 脱ブタン塔
    dist2: ColumnTunables       # 脱エタン塔
    dist3: ColumnTunables       # C3 スプリッタ
