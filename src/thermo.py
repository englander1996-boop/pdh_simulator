"""
PDH (Propane DeHydrogenation) シミュレータ — 熱力学データモジュール

Component mapping:
    A : Propane   (C3H8)  — プロパン
    B : Propylene (C3H6)  — プロピレン
    C : Hydrogen  (H2)    — 水素
    D : Ethylene  (C2H4)  — エチレン
    E : Methane   (CH4)   — メタン
    F : Ethane    (C2H6)  — エタン
"""

from typing import Dict

from .config import ThermoParams, THERMO_DATA


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
            成分記号 ('A'〜'F')

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
        'A' : Propane   (C3H8)  'D' : Ethylene  (C2H4)
        'B' : Propylene (C3H6)  'E' : Methane   (CH4)
        'C' : Hydrogen  (H2)    'F' : Ethane    (C2H6)

        Parameters
        ----------
        component : str
            成分記号 ('A'〜'F')
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
        'A' : Propane   (C3H8)  'D' : Ethylene  (C2H4)
        'B' : Propylene (C3H6)  'E' : Methane   (CH4)
        'C' : Hydrogen  (H2)    'F' : Ethane    (C2H6)

        Parameters
        ----------
        component : str
            成分記号 ('A'〜'F')
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
