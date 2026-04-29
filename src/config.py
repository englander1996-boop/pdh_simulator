"""
PDH (Propane DeHydrogenation) シミュレータ — 定数・パラメータ定義

Component mapping:
    A : Propane   (C3H8)
    B : Propylene (C3H6)
    C : Hydrogen  (H2)
    D : Ethylene  (C2H4)
    E : Methane   (CH4)
    F : Ethane    (C2H6)
    Z : n-Butane  (C4H10)
"""

from dataclasses import dataclass
from typing import Dict

# ---------------------------------------------------------------------------
# 物理定数
# ---------------------------------------------------------------------------

R: float = 8.31446   # [J K⁻¹ mol⁻¹]  理想気体定数
T0: float = 793.15   # [K]              基準温度 (= 520 °C)


# ---------------------------------------------------------------------------
# 反応パラメータ — 各反応ごとに独立したデータクラスで管理
# kJ → J 変換は乗算 (* 1_000) をここで完結させる。
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Reaction1Params:
    """
    反応1  脱水素 (Dehydrogenation): A (C3H8) → B (C3H6) + C (H2)

    r1 = a * k1 * (P_A − P_B·P_C / K_eq) / (1 + P_B / K_B)

    Attributes
    ----------
    k01    : 頻度因子          [mol m⁻³ s⁻¹ Pa⁻¹]
    Ea1    : 活性化エネルギー  [J mol⁻¹]   (元値: 34.57 kJ/mol)
    deltaH : 吸着エンタルピー  [J mol⁻¹]   (元値: −85.817 kJ/mol)
    K0     : K_B の前指数因子  [Pa]
    """

    k01: float = 9.787e-5          # [mol m⁻³ s⁻¹ Pa⁻¹]
    Ea1: float = 34.57 * 1_000     # [J mol⁻¹]   34.57 kJ/mol → J
    deltaH: float = -85.817 * 1_000  # [J mol⁻¹] −85.817 kJ/mol → J
    K0: float = 3.46e5             # [Pa]


@dataclass(frozen=True)
class Reaction2Params:
    """
    反応2  クラッキング (Cracking): A (C3H8) → D (C2H4) + E (CH4)

    r2 = k2 * P_A

    Attributes
    ----------
    k02 : 頻度因子          [mol m⁻³ s⁻¹ Pa⁻¹]
    Ea2 : 活性化エネルギー  [J mol⁻¹]   (元値: 137.31 kJ/mol)
    """

    k02: float = 8.682e-7          # [mol m⁻³ s⁻¹ Pa⁻¹]
    Ea2: float = 137.31 * 1_000    # [J mol⁻¹]   137.31 kJ/mol → J


@dataclass(frozen=True)
class Reaction3Params:
    """
    反応3  水素化 (Hydrogenation): D (C2H4) + C (H2) → F (C2H6)

    r3 = k3 * P_D * P_C

    Attributes
    ----------
    k03 : 頻度因子          [mol m⁻³ s⁻¹ Pa⁻²]
    Ea3 : 活性化エネルギー  [J mol⁻¹]   (元値: 154.54 kJ/mol)
    """

    k03: float = 4.406e-8          # [mol m⁻³ s⁻¹ Pa⁻²]
    Ea3: float = 154.54 * 1_000    # [J mol⁻¹]   154.54 kJ/mol → J


@dataclass(frozen=True)
class PDHConfig:
    """
    PDH シミュレータ全体の設定をまとめるルートデータクラス。

    Attributes
    ----------
    rxn1 : 反応1（脱水素）パラメータ
    rxn2 : 反応2（クラッキング）パラメータ
    rxn3 : 反応3（水素化）パラメータ
    """

    rxn1: Reaction1Params = Reaction1Params()
    rxn2: Reaction2Params = Reaction2Params()
    rxn3: Reaction3Params = Reaction3Params()


# ---------------------------------------------------------------------------
# 熱力学パラメータ
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ThermoParams:
    """
    1 成分の熱力学データを保持するデータクラス。

    定圧比熱の多項式:
        Cp = a + b*T + c*T**2 + d*T**3   [J/(K mol)]
    ※ フィールド名 a, b, c, d は Cp 多項式の係数であり、
    成分記号 A〜F、Z とは無関係。

    Attributes
    ----------
    dHf_298 : 標準生成エンタルピー (298 K) [J mol⁻¹]
              (入力値は kJ/mol → * 1_000 で J に変換済み)
    boiling_point : 沸点 [°C]
    a       : Cp 多項式 定数項       [J K⁻¹ mol⁻¹]
    b       : Cp 多項式 T の係数     [J K⁻² mol⁻¹]
    c       : Cp 多項式 T² の係数    [J K⁻³ mol⁻¹]
    d       : Cp 多項式 T³ の係数    [J K⁻⁴ mol⁻¹]
    """

    dHf_298: float  # [J mol⁻¹]
    boiling_point: float  # [°C]
    a: float        # [J K⁻¹ mol⁻¹]
    b: float        # [J K⁻² mol⁻¹]
    c: float        # [J K⁻³ mol⁻¹]
    d: float        # [J K⁻⁴ mol⁻¹]


# 成分 A〜F、Z の熱力学データ辞書
# キー : 成分記号 ('A'〜'F', 'Z')
# 値   : ThermoParams インスタンス
# dHf_298 は kJ/mol で与えられた値を * 1_000 により J/mol に変換
THERMO_DATA: Dict[str, ThermoParams] = {
    "A": ThermoParams(          # Propane   (C3H8)
        dHf_298 = -103.85 * 1_000,  # [J mol⁻¹]  (-103.85 kJ/mol → J)
        boiling_point = -42.05,
        a =  -4.224,
        b =   3.063e-1,
        c =  -1.586e-4,
        d =   3.215e-8,
    ),
    "B": ThermoParams(          # Propylene (C3H6)
        dHf_298 =  20.41 * 1_000,  # [J mol⁻¹]  (20.41 kJ/mol → J)
        boiling_point = -47.75,
        a =   3.710,
        b =   2.345e-1,
        c =  -1.16e-4,
        d =   2.205e-8,
    ),
    "C": ThermoParams(          # Hydrogen  (H2)
        dHf_298 =   0.0 * 1_000,   # [J mol⁻¹]  (0.0 kJ/mol → J)
        boiling_point = -252.75,
        a =  27.140,
        b =   9.274e-3,
        c =  -1.381e-5,
        d =   7.645e-9,
    ),
    "D": ThermoParams(          # Ethylene  (C2H4)
        dHf_298 =  52.28 * 1_000,  # [J mol⁻¹]  (52.28 kJ/mol → J)
        boiling_point = -103.75,
        a =   3.806,
        b =   1.566e-1,
        c =  -8.348e-5,
        d =   1.755e-8,
    ),
    "E": ThermoParams(          # Methane   (CH4)
        dHf_298 = -74.85 * 1_000,  # [J mol⁻¹]  (-74.85 kJ/mol → J)
        boiling_point = -161.45,
        a =  19.250,
        b =   5.213e-2,
        c =   1.197e-5,
        d =  -1.132e-8,
    ),
    "F": ThermoParams(          # Ethane    (C2H6)
        dHf_298 = -84.67 * 1_000,  # [J mol⁻¹]  (-84.67 kJ/mol → J)
        boiling_point = -88.65,
        a =   5.410,
        b =   1.781e-1,
        c =  -6.938e-5,
        d =   8.713e-9,
    ),
    "Z": ThermoParams(          # n-Butane  (C4H10)
        dHf_298 = -126.15 * 1_000,  # [J mol⁻¹]  (-126.15 kJ/mol → J)
        boiling_point = -0.45,
        a =  11.160,
        b =   3.310e-1,
        c =  -1.110e-4,
        d =  -2.820e-9,
    ),
}
