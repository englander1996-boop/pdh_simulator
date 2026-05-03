"""
PDH (Propane DeHydrogenation) シミュレータ — 反応速度計算モジュール

Component mapping:
    A : Propane   (C3H8)  — プロパン
    B : Propylene (C3H6)  — プロピレン
    C : Hydrogen  (H2)    — 水素
    D : Ethylene  (C2H4)  — エチレン
    E : Methane   (CH4)   — メタン  ※速度式には陽に現れない
    F : Ethane    (C2H6)  — エタン  ※速度式には陽に現れない

Reactions:
    r1 : A → B + C   脱水素 (Dehydrogenation)
    r2 : A → D + E   クラッキング (Cracking)
    r3 : D + C → F   水素化 (Hydrogenation)
"""

import math
from typing import Dict

from .config import PDHConfig, R, T0


class PDHKinetics:
    """
    PDH 反応系における各反応速度を計算するクラス。

    Component mapping
    -----------------
    A : Propane   (C3H8)  — プロパン
    B : Propylene (C3H6)  — プロピレン
    C : Hydrogen  (H2)    — 水素
    D : Ethylene  (C2H4)  — エチレン
    E : Methane   (CH4)   — メタン
    F : Ethane    (C2H6)  — エタン

    Parameters
    ----------
    config : PDHConfig
        反応パラメータをまとめた設定オブジェクト。
        省略時はデフォルト値（config.py に定義された文献値）を使用する。

    Examples
    --------
    >>> kin = PDHKinetics()
    >>> rates = kin.calculate(
    ...     P_A=50_000, P_B=10_000, P_C=10_000, P_D=1_000,
    ...     T=873.15, a=1.0, K_eq=2.0e10,
    ... )
    >>> rates["r1"], rates["r2"], rates["r3"]
    """

    def __init__(self, config: PDHConfig = PDHConfig()) -> None:
        self._cfg = config

    # ------------------------------------------------------------------
    # 速度定数・吸着平衡定数（温度依存部）
    # ------------------------------------------------------------------

    def _k1(self, T: float) -> float:
        """
        反応1 速度定数 k1 [mol m⁻³ s⁻¹ Pa⁻¹]

        k1 = k01 * exp( −Ea1/R * (1/T − 1/T0) )

        Parameters
        ----------
        T : float
            温度 [K]
        """
        p = self._cfg.rxn1
        return p.k01 * math.exp(-p.Ea1 / R * (1.0 / T - 1.0 / T0))

    def _K_B(self, T: float) -> float:
        """
        プロピレン吸着平衡定数 K_B [Pa]

        K_B = K0 * exp( −deltaH/R * (1/T − 1/T0) )

        Parameters
        ----------
        T : float
            温度 [K]
        """
        p = self._cfg.rxn1
        return p.K0 * math.exp(-p.deltaH / R * (1.0 / T - 1.0 / T0))

    def _k2(self, T: float) -> float:
        """
        反応2 速度定数 k2 [mol m⁻³ s⁻¹ Pa⁻¹]

        k2 = k02 * exp( −Ea2/R * (1/T − 1/T0) )

        Parameters
        ----------
        T : float
            温度 [K]
        """
        p = self._cfg.rxn2
        return p.k02 * math.exp(-p.Ea2 / R * (1.0 / T - 1.0 / T0))

    def _k3(self, T: float) -> float:
        """
        反応3 速度定数 k3 [mol m⁻³ s⁻¹ Pa⁻²]

        k3 = k03 * exp( −Ea3/R * (1/T − 1/T0) )

        Parameters
        ----------
        T : float
            温度 [K]
        """
        p = self._cfg.rxn3
        return p.k03 * math.exp(-p.Ea3 / R * (1.0 / T - 1.0 / T0))

    # ------------------------------------------------------------------
    # 各反応速度
    # ------------------------------------------------------------------

    def _r1(
        self,
        P_A: float,
        P_B: float,
        P_C: float,
        T: float,
        a: float,
        K_eq: float,
    ) -> float:
        """
        反応1 速度 r1 [mol m⁻³-cat s⁻¹]

        r1 = a * k1 * (P_A − P_B·P_C / K_eq) / (1 + P_B / K_B)

        単位は触媒体積基準 (m³-cat)。コンテスト要項 §3-2 参照。
        活性係数 a は脱水素反応 r1 のみに適用。
        コンテスト要項 §3-2 の速度式定義による仕様であり、r2・r3 には a が含まれない。

        Parameters
        ----------
        P_A   : 分圧 プロパン   (A: C3H8)  [Pa]
        P_B   : 分圧 プロピレン (B: C3H6)  [Pa]
        P_C   : 分圧 水素       (C: H2)    [Pa]
        T     : 温度                        [K]
        a     : 触媒劣化パラメータ（活性）  [-]
        K_eq  : 反応1 平衡定数              [Pa²]
        """
        driving_force = P_A - (P_B * P_C / K_eq)
        adsorption = 1.0 + P_B / self._K_B(T)
        return a * self._k1(T) * driving_force / adsorption

    def _r2(self, P_A: float, T: float) -> float:
        """
        反応2 速度 r2 [mol m⁻³-cat s⁻¹]

        r2 = k2 * P_A   ← 活性係数 a なし

        単位は触媒体積基準 (m³-cat)。コンテスト要項 §3-2 参照。
        コンテスト要項 §3-2 でクラッキング反応 r2 の式に a は含まれない。

        Parameters
        ----------
        P_A : 分圧 プロパン (A: C3H8) [Pa]
        T   : 温度                     [K]
        """
        return self._k2(T) * P_A

    def _r3(self, P_D: float, P_C: float, T: float) -> float:
        """
        反応3 速度 r3 [mol m⁻³-cat s⁻¹]

        r3 = k3 * P_D * P_C   ← 活性係数 a なし

        単位は触媒体積基準 (m³-cat)。コンテスト要項 §3-2 参照。
        コンテスト要項 §3-2 でエチレン水素化反応 r3 の式に a は含まれない。

        Parameters
        ----------
        P_D : 分圧 エチレン (D: C2H4) [Pa]
        P_C : 分圧 水素     (C: H2)   [Pa]
        T   : 温度                     [K]
        """
        return self._k3(T) * P_D * P_C

    # ------------------------------------------------------------------
    # 公開インターフェース
    # ------------------------------------------------------------------

    def calculate(
        self,
        P_A: float,
        P_B: float,
        P_C: float,
        P_D: float,
        T: float,
        a: float,
        K_eq: float,
    ) -> Dict[str, float]:
        """
        全反応速度を一括計算して辞書形式で返す。

        Component mapping
        -----------------
        A : Propane   (C3H8)  — プロパン
        B : Propylene (C3H6)  — プロピレン
        C : Hydrogen  (H2)    — 水素
        D : Ethylene  (C2H4)  — エチレン

        Parameters
        ----------
        P_A  : 分圧 プロパン   (A: C3H8)  [Pa]
        P_B  : 分圧 プロピレン (B: C3H6)  [Pa]
        P_C  : 分圧 水素       (C: H2)    [Pa]
        P_D  : 分圧 エチレン   (D: C2H4)  [Pa]
        T    : 温度                        [K]
        a    : 触媒劣化パラメータ（0〜1）  [-]
              1.0 = 新鮮な触媒、0.0 = 完全失活
        K_eq : 反応1（脱水素）の平衡定数   [Pa²]
              外部（熱力学モデル等）から供給すること。

        Returns
        -------
        Dict[str, float]
            キー "r1", "r2", "r3" に対応する反応速度 [mol m⁻³ s⁻¹]

            r1 : A → B + C   脱水素反応速度
            r2 : A → D + E   クラッキング反応速度
            r3 : D + C → F   水素化反応速度

        Examples
        --------
        >>> kin = PDHKinetics()
        >>> rates = kin.calculate(
        ...     P_A=50_000, P_B=10_000, P_C=10_000, P_D=500,
        ...     T=873.15, a=1.0, K_eq=2.0e10,
        ... )
        >>> print(rates)
        {'r1': ..., 'r2': ..., 'r3': ...}
        """
        return {
            "r1": self._r1(P_A, P_B, P_C, T, a, K_eq),
            "r2": self._r2(P_A, T),
            "r3": self._r3(P_D, P_C, T),
        }
