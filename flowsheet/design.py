"""
最適化対象の設計変数バンドル。

最適化器は FlowsheetDesignVars インスタンスを生成して runner に渡す。
設計変数の追加 (たとえば膜入口温度を最適化するなど) はこの dataclass を
拡張する形で行う。
"""

from dataclasses import dataclass
from typing import Union

from units.reactors.swing import DesignVars as SwingDesign
from units.reactors.radial_flow import RadialDesignVars
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from src.distillation_core import ColumnTunables


@dataclass(frozen=True)
class FlowsheetDesignVars:
    """フローシート全体の設計変数。最適化器はこれを生成して runner へ渡す。

    蒸留塔の物理セマンティクス (LK/HK/回収率/K_method/q) は塔別ラッパで固定し、
    本バンドルでは BO で振る P/N/R の組のみを保持する (ColumnTunables)。

    反応器モデル (2026-05-30): `swing` フィールドには 2 種の反応器設計変数を入れられる:
      - SwingDesign(T_in, z_cat, t_cyc, D)              … 軸流固定床 (units/reactors/swing.py)
      - RadialDesignVars(T_in, t_cyc, D_inner, bed_thickness, H) … 径方向流 (units/reactors/radial_flow.py)
    run_one_pass が型でディスパッチして対応する反応器を実行する。フィールド名は後方互換のため
    `swing` のまま (= 既存の軸流コンストラクタを無改修で維持)。両者とも T_in/t_cyc を持ち、
    runner の HI (design.swing.T_in) はどちらでも動く。
    """
    swing: Union[SwingDesign, RadialDesignVars]
    psa:   PSADesignVars
    mem:   MemDesignVars
    dist1: ColumnTunables       # 脱ブタン塔
    dist2: ColumnTunables       # 脱エタン塔
    dist3: ColumnTunables       # C3 スプリッタ
