"""
PDH (Propane DeHydrogenation) シミュレータ — 熱力学データモジュール

Component mapping:
    A : Propane   (C3H8)  — プロパン
    B : Propylene (C3H6)  — プロピレン
    C : Hydrogen  (H2)    — 水素
    D : Ethylene  (C2H4)  — エチレン
    E : Methane   (CH4)   — メタン
    F : Ethane    (C2H6)  — エタン
    Z : n-Butane  (C4H10)  — ノルマルブタン
"""

import math
from typing import Dict

from .config import ThermoParams, THERMO_DATA, R

# 標準圧力 [Pa]（K_eq の単位変換用: K_eq[Pa] = P_STD * Kp_dimensionless）
_P_STD: float = 101325.0


class PDHThermo:
    """
    PDH 反応系における熱力学量（Cp・エンタルピー変化）を計算するクラス。

    定圧比熱の多項式:
        Cp(T) = a + b*T + c*T**2 + d*T**3   [J/(K mol)]

    エンタルピー変化の解析的積分:
        ΔH = ∫_{T_start}^{T_end} Cp(T) dT
           = a*(T_end - T_start)
           + b/2 * (T_end**2 - T_start**2)
           + c/3 * (T_end**3 - T_start**3)
           + d/4 * (T_end**4 - T_start**4)

    Component mapping
    -----------------
    'A' : Propane   (C3H8)  — プロパン
    'B' : Propylene (C3H6)  — プロピレン
    'C' : Hydrogen  (H2)    — 水素
    'D' : Ethylene  (C2H4)  — エチレン
    'E' : Methane   (CH4)   — メタン
    'F' : Ethane    (C2H6)  — エタン
    'Z' : n-Butane  (C4H10)  — ノルマルブタン

    Parameters
    ----------
    data : Dict[str, ThermoParams], optional
        成分ごとの熱力学パラメータ辞書。
        省略時は config.py に定義された THERMO_DATA を使用する。

    Examples
    --------
    >>> thermo = PDHThermo()
    >>> cp = thermo.calc_cp("A", T=873.15)
    >>> dH = thermo.calc_enthalpy_change("A", T_start=298.15, T_end=873.15)
    """

    def __init__(self, data: Dict[str, ThermoParams] = THERMO_DATA) -> None:
        self._data = data

    # ------------------------------------------------------------------
    # 内部ヘルパー
    # ------------------------------------------------------------------

    def _get(self, component: str) -> ThermoParams:
        """
        成分キーに対応する ThermoParams を返す。

        Parameters
        ----------
        component : str
            成分記号 ('A'〜'F', 'Z')

        Raises
        ------
        KeyError
            未定義の成分記号が渡された場合。
        """
        if component not in self._data:
            valid = ", ".join(sorted(self._data.keys()))
            raise KeyError(
                f"未知の成分記号: '{component}'。有効な記号: {valid}"
            )
        return self._data[component]

    # ------------------------------------------------------------------
    # 公開インターフェース
    # ------------------------------------------------------------------

    def calc_cp(self, component: str, T: float) -> float:
        """
        指定成分の温度 T における定圧比熱 Cp を返す。

        Cp(T) = a + b*T + c*T**2 + d*T**3   [J/(K mol)]

        Component mapping
        -----------------
        'A' : Propane   (C3H8)  'D' : Ethylene  (C2H4)  'Z' : n-Butane  (C4H10)
        'B' : Propylene (C3H6)  'E' : Methane   (CH4)
        'C' : Hydrogen  (H2)    'F' : Ethane    (C2H6)

        Parameters
        ----------
        component : str
            成分記号 ('A'〜'F', 'Z')
        T : float
            温度 [K]

        Returns
        -------
        float
            定圧比熱 Cp [J/(K mol)]

        Examples
        --------
        >>> thermo = PDHThermo()
        >>> thermo.calc_cp("A", T=873.15)
        """
        p = self._get(component)
        return p.a + p.b * T + p.c * T**2 + p.d * T**3

    def calc_enthalpy_change(
        self,
        component: str,
        T_start: float,
        T_end: float,
    ) -> float:
        """
        指定成分の T_start → T_end におけるエンタルピー変化を返す。

        Cp の多項式を解析的に積分することで高速に計算する:

            ΔH = ∫_{T_start}^{T_end} Cp(T) dT
               = a * (T1 - T0)
               + b/2  * (T1**2 - T0**2)
               + c/3  * (T1**3 - T0**3)
               + d/4  * (T1**4 - T0**4)

        ここで T0 = T_start, T1 = T_end。

        Component mapping
        -----------------
        'A' : Propane   (C3H8)  'D' : Ethylene  (C2H4)  'Z' : n-Butane  (C4H10)
        'B' : Propylene (C3H6)  'E' : Methane   (CH4)
        'C' : Hydrogen  (H2)    'F' : Ethane    (C2H6)

        Parameters
        ----------
        component : str
            成分記号 ('A'〜'F', 'Z')
        T_start : float
            積分開始温度 [K]
        T_end : float
            積分終了温度 [K]

        Returns
        -------
        float
            エンタルピー変化 ΔH [J mol⁻¹]
            T_end > T_start のとき正値（吸熱方向を正）

        Examples
        --------
        >>> thermo = PDHThermo()
        >>> thermo.calc_enthalpy_change("A", T_start=298.15, T_end=873.15)
        """
        p = self._get(component)
        T0, T1 = T_start, T_end

        return (
            p.a         * (T1    - T0)
            + p.b / 2.0 * (T1**2 - T0**2)
            + p.c / 3.0 * (T1**3 - T0**3)
            + p.d / 4.0 * (T1**4 - T0**4)
        )

    def calc_keq(self, T: float) -> float:
        """
        反応1 (C3H8 → C3H6 + H2) の平衡定数 K_eq [Pa] を返す。

        Kirchhoff の法則と Gibbs-Helmholtz 式に基づく厳密計算。

        計算経路
        --------
        1. THERMO_DATA の dHf_298, dGf_298 から ΔH°_r(298), ΔG°_r(298) を算出
        2. ΔS°_r(298) = (ΔH° - ΔG°) / 298.15
        3. Cp 差分係数 Δa, Δb, Δc, Δd を計算（化学量論: νA=-1, νB=+1, νC=+1）
        4. ΔH°(T) を Kirchhoff で解析的に積分
        5. ΔS°(T) = ΔS°(298) + ∫(ΔCp/T) dT を解析的に積分
        6. ΔG°(T) = ΔH°(T) - T·ΔS°(T)
        7. K_eq [Pa] = P_STD · exp(-ΔG°(T) / (R·T))

        Parameters
        ----------
        T : float
            温度 [K]

        Returns
        -------
        float
            平衡定数 K_eq [Pa]

        Examples
        --------
        >>> thermo = PDHThermo()
        >>> thermo.calc_keq(T=873.15)
        """
        T0 = 298.15  # [K] 基準温度

        pA = self._get('A')  # C3H8
        pB = self._get('B')  # C3H6
        pC = self._get('C')  # H2

        # ΔH°_r(298), ΔG°_r(298) [J/mol]
        dH298 = pB.dHf_298 + pC.dHf_298 - pA.dHf_298
        dG298 = pB.dGf_298 + pC.dGf_298 - pA.dGf_298

        # ΔS°_r(298) [J/(mol·K)]
        dS298 = (dH298 - dG298) / T0

        # Cp 差分係数
        da = pB.a + pC.a - pA.a
        db = pB.b + pC.b - pA.b
        dc = pB.c + pC.c - pA.c
        dd = pB.d + pC.d - pA.d

        # ΔH°(T) [J/mol] — Kirchhoff の法則
        dH_T = (dH298
                + da       * (T    - T0)
                + db / 2.0 * (T**2 - T0**2)
                + dc / 3.0 * (T**3 - T0**3)
                + dd / 4.0 * (T**4 - T0**4))

        # ΔS°(T) [J/(mol·K)] — ∫(ΔCp/T) dT の解析的積分
        dS_T = (dS298
                + da        * math.log(T / T0)
                + db        * (T    - T0)
                + dc / 2.0  * (T**2 - T0**2)
                + dd / 3.0  * (T**3 - T0**3))

        # ΔG°(T) [J/mol]
        dG_T = dH_T - T * dS_T

        # K_eq [Pa] = P_STD · Kp (無次元)
        # exp 引数を上限クランプ (math.exp の最大引数 ≈ 709 を超えると OverflowError)
        exp_arg = min(-dG_T / (R * T), 700.0)
        return _P_STD * math.exp(exp_arg)
