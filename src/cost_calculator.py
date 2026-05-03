"""
Bare Module Cost法の計算ロジック

対象装置: 縦型プロセス容器（Vertical process vessel）
出典: 授業資料 プロセス設計R08-3.pdf

計算フロー
----------
1. 基本コスト    Cp0 = 10^(K1 + K2*log10(A) + K3*(log10(A))^2)
2. 圧力係数      Fp  = f(Pg [bar], D [m])
3. モジュール係数 FBM = B1 + B2*Fp*FM
4. ベアモジュール CBM = Cp0 * FBM * K_swing
5. 現在価値      Current_Cost = CBM * (CEPCI_CURRENT / CEPCI_BASE)
6. 総建設費      C_TM = 1.18 * Current_Cost * N_reactors_total
7. TAC          TAC  = C_TM / 減価償却年数 + Annual_OPEX
"""

import math
import warnings

from .cost_parameters import (
    K1, K2, K3, A_RANGE_MIN, A_RANGE_MAX,
    B1, B2, FM,
    CEPCI_BASE, CEPCI_CURRENT,
    K_SWING,
    USD_TO_JPY, DEPRECIATION_YEARS, PLANT_INDIRECT_FACTOR,
    K1_HE, K2_HE, K3_HE, A_HE_MIN, A_HE_MAX, B1_HE, B2_HE, FM_HE, FP_HE_DEFAULT,
    K1_COMP, K2_COMP, K3_COMP, W_COMP_MIN, W_COMP_MAX, FBM_COMP,
    MEM_UNIT_PRICE_USD_PER_M2,
)


def calc_cp0(A: float) -> float:
    """
    基本コスト Cp0 [USD（2001年基準）]

    Cp0 = 10^(K1 + K2*log10(A) + K3*(log10(A))^2)

    Parameters
    ----------
    A : float
        容器体積 [m³]（推算式適用範囲: 0.3〜520 m³）

    Returns
    -------
    float
        基本コスト Cp0 [USD]
    """
    if A <= 0:
        raise ValueError(f"calc_cp0: A={A} は正値でなければなりません。")
    if not (A_RANGE_MIN <= A <= A_RANGE_MAX):
        warnings.warn(
            f"calc_cp0: A={A:.3f} m³ は推算式の適用範囲 "
            f"[{A_RANGE_MIN}, {A_RANGE_MAX}] m³ 外です。",
            UserWarning,
            stacklevel=2,
        )
    log_A = math.log10(A)
    return 10.0 ** (K1 + K2 * log_A + K3 * log_A ** 2)


def calc_fp(P_gauge_bar: float, D: float) -> float:
    """
    圧力係数 Fp [-]

    Pg > -0.5 の場合:
        Fp = max((Pg+1)*D / (10.71 - 0.00756*(Pg+1)) + 0.5, 1.0)
    Pg ≤ -0.5 の場合（真空操作）:
        Fp = 1.25

    Parameters
    ----------
    P_gauge_bar : float
        ゲージ圧力 [bar]（= 絶対圧 [bar] - 1.01325）
    D : float
        容器内径 [m]

    Returns
    -------
    float
        圧力係数 Fp [-]（常に ≥ 1.0 または 1.25）
    """
    if P_gauge_bar > -0.5:
        numerator = (P_gauge_bar + 1.0) * D
        denominator = 10.71 - 0.00756 * (P_gauge_bar + 1.0)
        if denominator <= 0:
            # 現実的でない超高圧（>1400 bar）: 大きな Fp を返して最適化を遠ざける
            return 10.0
        return max(numerator / denominator + 0.5, 1.0)
    else:
        return 1.25


def calc_reactor_capex_okuyen(
    V_vessel_m3: float,
    P_abs_pa: float,
    D_m: float,
    N_reactors_total: int,
) -> float:
    """
    反応器システム全体の総建設費 C_TM [億円] を計算する。

    Parameters
    ----------
    V_vessel_m3 : float
        1基あたりのプロセス容器体積 [m³]
    P_abs_pa : float
        絶対圧力 [Pa]（ゲージ圧変換に使用）
    D_m : float
        容器内径 [m]
    N_reactors_total : int
        システム全体の総反応器基数

    Returns
    -------
    float
        総建設費 C_TM [億円]
    """
    # Pa → bar ゲージ圧変換
    P_abs_bar = P_abs_pa / 1.0e5
    P_gauge_bar = P_abs_bar - 1.01325

    cp0 = calc_cp0(V_vessel_m3)
    fp  = calc_fp(P_gauge_bar, D_m)
    fbm = B1 + B2 * fp * FM
    cbm = cp0 * fbm * K_SWING

    current_cost_usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    c_tm_usd = PLANT_INDIRECT_FACTOR * current_cost_usd * N_reactors_total
    return c_tm_usd * USD_TO_JPY / 1.0e8


# ---------------------------------------------------------------------------
# 熱交換器（固定管板式）
# 出典: R08-3.pdf p.9 Table A.1 / p.10 Table A.4
# ---------------------------------------------------------------------------

def calc_he_capex_okuyen(
    A_m2: float,
    Fp: float = FP_HE_DEFAULT,
    FM: float = FM_HE,
) -> float:
    """
    熱交換器（固定管板式）の総建設費 [億円]。

    C_BM = Cp0 × (B1_HE + B2_HE × Fp × FM)
    C_TM = PLANT_INDIRECT_FACTOR × C_BM × CEPCI補正

    Parameters
    ----------
    A_m2 : 伝熱面積 [m²]
    Fp   : 圧力補正係数 [-]（デフォルト 1.0、低〜中圧仮定）
    FM   : 材質係数 [-]（デフォルト 1.0、炭素鋼）
    """
    if A_m2 <= 0:
        raise ValueError(f"calc_he_capex_okuyen: A={A_m2} は正値でなければなりません。")
    if not (A_HE_MIN <= A_m2 <= A_HE_MAX):
        warnings.warn(
            f"calc_he_capex_okuyen: A={A_m2:.2f} m² は適用範囲 "
            f"[{A_HE_MIN}, {A_HE_MAX}] m² 外です。",
            UserWarning, stacklevel=2,
        )

    log_A  = math.log10(A_m2)
    cp0    = 10.0 ** (K1_HE + K2_HE * log_A + K3_HE * log_A ** 2)
    cbm    = cp0 * (B1_HE + B2_HE * Fp * FM)
    current_usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * current_usd * USD_TO_JPY / 1.0e8


# ---------------------------------------------------------------------------
# 圧縮機（遠心式、炭素鋼）
# 出典: R08-3.pdf p.9 Table A.1 / p.13–14 Table A.6
# ---------------------------------------------------------------------------

def calc_comp_capex_okuyen(W_kW: float) -> float:
    """
    遠心式圧縮機の総建設費 [億円]。

    C_BM = Cp0 × FBM_COMP  (Fp = 1.0、圧力補正なし)
    C_TM = PLANT_INDIRECT_FACTOR × C_BM × CEPCI補正

    Parameters
    ----------
    W_kW : 流体動力（実圧縮仕事） [kW]
    """
    if W_kW <= 0:
        raise ValueError(f"calc_comp_capex_okuyen: W={W_kW} は正値でなければなりません。")
    if not (W_COMP_MIN <= W_kW <= W_COMP_MAX):
        warnings.warn(
            f"calc_comp_capex_okuyen: W={W_kW:.1f} kW は適用範囲 "
            f"[{W_COMP_MIN}, {W_COMP_MAX}] kW 外です。",
            UserWarning, stacklevel=2,
        )

    log_W  = math.log10(W_kW)
    cp0    = 10.0 ** (K1_COMP + K2_COMP * log_W + K3_COMP * log_W ** 2)
    cbm    = cp0 * FBM_COMP
    current_usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * current_usd * USD_TO_JPY / 1.0e8


# ---------------------------------------------------------------------------
# 膜モジュール
# ★★ MEM_UNIT_PRICE_USD_PER_M2 が仮置き値のため、出力も暫定値 ★★
# ---------------------------------------------------------------------------

def calc_mem_capex_okuyen(A_mem_m2: float) -> float:
    """
    膜モジュールの設備費 [億円]。

    C_mem = MEM_UNIT_PRICE_USD_PER_M2 × A_mem × CEPCI補正 × 据付間接費 × 円換算

    ★★ 仮置き値を使用中 ★★
    MEM_UNIT_PRICE_USD_PER_M2 は根拠文献未確定の暫定値。
    ZIF-8 膜 TEA 論文（Hua et al. 2024 等）で単価を確認後に
    cost_parameters.MEM_UNIT_PRICE_USD_PER_M2 を更新すること。

    Parameters
    ----------
    A_mem_m2 : 総膜面積 [m²]
    """
    warnings.warn(
        f"calc_mem_capex_okuyen: 膜モジュール単価 MEM_UNIT_PRICE_USD_PER_M2 ="
        f" {MEM_UNIT_PRICE_USD_PER_M2} USD/m² は仮置き値です。"
        " ZIF-8 膜 TEA 論文等で根拠値を確認後に cost_parameters を更新してください。",
        UserWarning,
        stacklevel=2,
    )
    if A_mem_m2 <= 0:
        raise ValueError(f"calc_mem_capex_okuyen: A={A_mem_m2} は正値でなければなりません。")
    cost_usd = MEM_UNIT_PRICE_USD_PER_M2 * A_mem_m2 * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * cost_usd * USD_TO_JPY / 1.0e8

