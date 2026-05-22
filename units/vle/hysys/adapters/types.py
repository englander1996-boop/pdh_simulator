"""アダプタ層の入出力契約 (dataclass)。"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class _ColumnInputCommon:
    """全塔共通の入力。

    feed_composition は PDH キー (A/B/C/D/E/F/Z) → kmol/h の辞書。
    上流から流れてきた実組成を HYSYS feed に毎回上書きする想定。
    """
    feed_composition: Dict[str, float]
    feed_temperature_c: Optional[float] = None    # None なら HSC の値を維持
    feed_pressure_kpa: Optional[float] = None     # None なら col_p_kpa + 50 を採用


@dataclass(frozen=True)
class Column1Input(_ColumnInputCommon):
    """Dist1 (脱ブタン塔) の入力契約。

    lhs_column1.py の DOE 変数群 + 上流流量・組成:
      - col_p_kpa, feed_flow_kmolh, comp_frac_2, feed_stage
    feed_composition は上流 (Pump1 後 Fresh LPG 等) からの実組成。
    """
    col_p_kpa: float = 0.0
    feed_flow_kmolh: float = 0.0
    comp_frac_2: float = 0.99   # Comp Fraction - 2 スペック目標値 [-]
    feed_stage: int = 1


@dataclass(frozen=True)
class Column2Input(_ColumnInputCommon):
    """Dist2 (脱エタン塔) の入力契約。

    lhs_column2.py の DOE 変数群:
      - col_p_kpa, feed_flow_kmolh, reflux_ratio, feed_stage
    組成は feed_composition (PDH 側から渡される) を使う。
    """
    col_p_kpa: float = 0.0
    feed_flow_kmolh: float = 0.0
    reflux_ratio: float = 1.0
    feed_stage: int = 1


@dataclass(frozen=True)
class Column3Input(_ColumnInputCommon):
    """Dist3 (C3 スプリッタ) の入力契約。

    lhs_column3.py の DOE 変数群:
      - col_p_kpa, feed_flow_kgmols (kgmol/s 直書き), propane_frac, draw_rate_kgmols, feed_stage
    propane_frac は lhs_column3.py の運用に合わせて Propane/Propene 2 成分のみ扱う
    モードを残しているが、feed_composition (全成分) が渡されれば組成丸ごと上書きする。
    """
    col_p_kpa: float = 0.0
    feed_flow_kgmols: float = 0.0
    propane_frac: Optional[float] = None    # None なら feed_composition で全成分上書き
    draw_rate_kgmols: float = 0.0
    feed_stage: int = 1


@dataclass
class ColumnResult:
    """全塔共通の出力契約。

    HYSYS から取れる出口値だけを保持する。CAPEX 等は上位層 (HysysVLEProvider)
    で cost_calculator を呼んで後計算する。
    """
    success: bool = False
    message: str = ""
    # 熱量 [MW]
    qc_mw: float = 0.0
    qr_mw: float = 0.0
    # 塔頂・塔底 (流量 [kmol/h], 温度 [°C], 圧力 [kPa], 組成は PDH キー → mol frac)
    top_flow_kmolh: float = 0.0
    top_temp_c: float = 0.0
    top_pressure_kpa: float = 0.0
    top_composition: Dict[str, float] = field(default_factory=dict)
    bot_flow_kmolh: float = 0.0
    bot_temp_c: float = 0.0
    bot_pressure_kpa: float = 0.0
    bot_composition: Dict[str, float] = field(default_factory=dict)
    # 還流比 (取れたら格納)
    reflux_ratio: Optional[float] = None
    # フィード・各製品のエンタルピー流量 [MW] (heat integration 用、取れたら格納)
    feed_enthalpy_mw: Optional[float] = None
    top_enthalpy_mw: Optional[float] = None
    bot_enthalpy_mw: Optional[float] = None
    # ポップアップ警告の集計
    warning_count: int = 0
    warning_details: str = ""
