"""
膜分離システム シミュレーター

5ユニット構成:
  1. 気化器 (Vaporizer)         前段蒸留塔底液 → 全量ガス化
  2. フィード圧縮機 (FeedComp)  P_in → P_H
  3. 膜分離モジュール (Membrane) クロスフロー ODE
  4. 製品圧縮機 (ProductComp)   P_L → P_dist
  5. 製品冷却器 (Condenser)     → 飽和液（後段蒸留塔フィード）

接続:
  [前段蒸留塔底液] --MemFeedStream--> simulate_membrane_system
  --> MemSimulationResult.product   --> [後段蒸留塔フィード（飽和液）]
  --> MemSimulationResult.retentate --> [リサイクルまたは排出]

使用方法:
    from units.separators.membrane_system import (
        MemDesignVars, MemFeedStream, MemFixedParams,
        simulate_membrane_system,
    )
    result = simulate_membrane_system(design, feed, fixed)
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.eos import (
    z_factor, residual_enthalpy,
    bubble_point_T, dew_point_T,
    compress_isentropic,
)
from src.thermo import PDHThermo

# 成分順序: index 0 = C3H6 (EOS キー 'B'), index 1 = C3H8 (EOS キー 'A')
_KEYS    = ['B', 'A']
_T_REF   = 298.15        # [K]  エンタルピー基準温度
_GPU_SI  = 3.35e-10      # [mol/(m²·s·Pa)]  1 GPU の SI 換算値
_ATM_BAR = 1.01325       # [bar] 大気圧（ゲージ圧変換用）

_thermo = PDHThermo()


# ---------------------------------------------------------------------------
# データクラス: 入力
# ---------------------------------------------------------------------------

@dataclass
class MemDesignVars:
    """最適化アルゴリズムが操作する設計変数"""
    P_H:   float  # 膜供給側（高圧）圧力 [Pa]
    P_L:   float  # 膜透過側（低圧）圧力 [Pa]
    A_mem: float  # 総膜面積 [m²]


@dataclass
class MemFeedStream:
    """前段蒸留塔底液 → 気化器への入力ストリーム"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_in:   float  # 液温度 [K]（飽和液または過冷却液）
    P_in:   float  # 圧力 [Pa]


@dataclass
class MemFixedParams:
    """プロセス固定パラメータ"""
    # 膜性能
    Q_A_GPU:         float = 40.0    # C3H6 透過度 [GPU]
    alpha:           float = 90.0    # C3H6/C3H8 選択性 [-]
    A_per_module:    float = 500.0   # モジュール 1 本あたり有効膜面積 [m²]
    # 後段蒸留塔接続
    P_dist:          float = 15.0e5  # 後段蒸留塔操作圧力 [Pa]
    # 気化器
    T_vap_superheat: float = 5.0     # 露点超過の過熱度 [K]
    U_vap:           float = 1.5     # 気化器総括伝熱係数 [kW/(m²·K)]
    T_hot:           float = 423.15  # 熱媒（低圧蒸気）温度 [K] ≈ 150°C
    # 製品冷却器
    U_cond:          float = 1.0     # 冷却器総括伝熱係数 [kW/(m²·K)]
    T_cold_in:       float = 303.15  # 冷却水入口温度 [K] (30°C)
    T_cold_out:      float = 318.15  # 冷却水出口温度 [K] (45°C)
    # 圧縮機
    eta_comp:        float = 0.75    # 断熱効率 [-]

    def __post_init__(self) -> None:
        if self.P_dist <= 0:
            raise ValueError("P_dist は正値でなければなりません。")
        if not 0 < self.eta_comp <= 1.0:
            raise ValueError("eta_comp は (0, 1] でなければなりません。")
        if self.T_hot <= 273.15:
            raise ValueError("T_hot は 273K 超でなければなりません。")


# ---------------------------------------------------------------------------
# データクラス: 出力
# ---------------------------------------------------------------------------

@dataclass
class MemRetentateStream:
    """膜非透過ストリーム（C3H8 富化; リサイクルまたは排出）"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_out:  float  # 温度 [K]（フィード圧縮機出口温度 = 膜等温仮定）
    P_out:  float  # 圧力 [Pa] = P_H


@dataclass
class MemProductStream:
    """製品ストリーム → 後段蒸留塔フィード（飽和液）"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_out:  float  # 泡点温度 [K] @ P_dist（飽和液）
    P_out:  float  # 圧力 [Pa] = P_dist


@dataclass
class MemEquipmentData:
    """機器サイズ・コスト推算テーブル（授業資料 R08-3 形式）

    CAPEX は未実装（float('nan')）。事後に cost_calculator で計算予定。

    各フィールドの意味:
      A_vap / A_cond   : 熱交換器 伝熱面積 [m²]    → 機器タイプ: 熱交換器
      W_feed / W_prod  : 圧縮機 所要動力 [kW]       → 機器タイプ: 圧縮機
      A_mem / n_modules: 膜モジュール面積・本数      → 機器タイプ: 特殊機器
      Pg_*             : ゲージ圧 [barg]             → 圧力補正係数 Fp 計算用
      Q_vap / Q_cond   : OPEX 用熱量 [kW]
    """
    # 気化器
    A_vap:        float  # 伝熱面積 [m²]
    Pg_vap:       float  # ゲージ圧 [barg]
    Q_vap_kW:     float  # 加熱量 [kW]
    # フィード圧縮機
    W_feed_kW:    float  # 所要動力 [kW]
    Pg_feed:      float  # ゲージ圧 [barg]
    # 膜モジュール
    A_mem:        float  # 総膜面積 [m²]
    n_modules:    int    # モジュール本数
    Pg_mem:       float  # ゲージ圧 [barg]
    # 製品圧縮機
    W_prod_kW:    float  # 所要動力 [kW]
    Pg_prod:      float  # ゲージ圧 [barg]
    # 製品冷却器
    A_cond:       float  # 伝熱面積 [m²]
    Pg_cond:      float  # ゲージ圧 [barg]
    Q_cond_kW:    float  # 冷却量 [kW]
    # CAPEX
    CAPEX_total:  float = float('nan')  # [億円]（後で実装）


@dataclass
class MemSimulationResult:
    """膜分離システム シミュレーション全出力"""
    retentate:   MemRetentateStream
    product:     MemProductStream
    equipment:   MemEquipmentData
    stage_cut:   float  # θ = F_perm / F_feed [-]
    perm_purity: float  # 透過ガス C3H6 モル分率 [-]
    ret_purity:  float  # 非透過ガス C3H6 モル分率 [-]


# ---------------------------------------------------------------------------
# ペナルティ結果
# ---------------------------------------------------------------------------

_PENALTY = 1e9


def _penalty_result() -> MemSimulationResult:
    """計算不能な条件のときに返すペナルティ結果。"""
    zero_stream = MemRetentateStream(0.0, 0.0, 0.0, 0.0)
    zero_prod   = MemProductStream(0.0, 0.0, 0.0, 0.0)
    eq = MemEquipmentData(
        A_vap=0.0, Pg_vap=0.0, Q_vap_kW=0.0,
        W_feed_kW=0.0, Pg_feed=0.0,
        A_mem=0.0, n_modules=0, Pg_mem=0.0,
        W_prod_kW=0.0, Pg_prod=0.0,
        A_cond=0.0, Pg_cond=0.0, Q_cond_kW=0.0,
        CAPEX_total=_PENALTY,
    )
    return MemSimulationResult(zero_stream, zero_prod, eq, 0.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# 内部ヘルパー
# ---------------------------------------------------------------------------

def _h_mol(T: float, P: float, z_C3H6: float, phase: str) -> float:
    """
    C3H6/C3H8 混合物のモルエンタルピー [J/mol]。
    基準: T_ref=298.15K の理想気体。

    H = H_ig(T_ref→T) + H^r(T, P, x, Z_phase)
    """
    x = [z_C3H6, 1.0 - z_C3H6]
    H_ig = sum(
        x[i] * _thermo.calc_enthalpy_change(_KEYS[i], _T_REF, T)
        for i in range(2)
    )
    Z   = z_factor(T, P, x, _KEYS, phase)
    H_r = residual_enthalpy(T, P, x, _KEYS, Z)
    return H_ig + H_r


def _lmtd(dT1: float, dT2: float) -> float:
    """対数平均温度差 LMTD [K]。dT1 ≠ dT2 の場合のみ厳密計算。"""
    if dT1 <= 0 or dT2 <= 0:
        return float('nan')
    if abs(dT1 - dT2) < 1e-3:
        return (dT1 + dT2) / 2.0
    return (dT1 - dT2) / math.log(dT1 / dT2)


# ---------------------------------------------------------------------------
# ユニット 1: 気化器
# ---------------------------------------------------------------------------

def _vaporizer(F_feed_mols: float, z_C3H6: float,
               T_in: float, P_in: float,
               fixed: MemFixedParams):
    """
    気化器: 液フィードを全量ガス化。

    Returns
    -------
    T_vap_out : 出口ガス温度 [K]（露点 + 過熱度）
    Q_vap_kW  : 必要加熱量 [kW]
    A_vap     : 伝熱面積 [m²]
    """
    T_dew = dew_point_T(P_in, [z_C3H6, 1.0 - z_C3H6], _KEYS)
    T_vap_out = T_dew + fixed.T_vap_superheat

    if T_vap_out >= fixed.T_hot:
        warnings.warn(
            f"気化器: T_vap_out={T_vap_out:.1f}K が T_hot={fixed.T_hot:.1f}K 以上です。"
            " T_hot を上げてください。",
            UserWarning, stacklevel=3,
        )
        return T_vap_out, float('nan'), float('nan')

    H_liq = _h_mol(T_in,      P_in, z_C3H6, 'liquid')
    H_vap = _h_mol(T_vap_out, P_in, z_C3H6, 'vapor')
    Q_vap = F_feed_mols * (H_vap - H_liq)    # [J/s = W]
    Q_vap_kW = Q_vap / 1e3

    dT1 = fixed.T_hot - T_in
    dT2 = fixed.T_hot - T_vap_out
    lmtd = _lmtd(dT1, dT2)
    A_vap = Q_vap_kW / (fixed.U_vap * lmtd) if lmtd > 0 else float('nan')

    return T_vap_out, Q_vap_kW, A_vap


# ---------------------------------------------------------------------------
# ユニット 3: 膜モジュール（クロスフロー ODE）
# ---------------------------------------------------------------------------

def _y_local(x: float, alpha: float, gamma: float) -> float:
    """
    クロスフローモデルにおける局所透過組成 y_local [-] を返す。

    x × P_H と y_local × P_L の分圧差が推進力。
    (1-alpha)*gamma * y² + [(alpha-1)*(x+gamma)+1] * y - alpha*x = 0
    の物理根（正かつ ≤ 1 の根）を返す。
    """
    a = (1.0 - alpha) * gamma
    b = (alpha - 1.0) * (x + gamma) + 1.0
    c = -alpha * x
    disc = max(0.0, b**2 - 4.0*a*c)
    denom = -b - math.sqrt(disc)
    if abs(denom) < 1e-30:
        return x  # フォールバック
    return max(0.0, min(1.0, (2.0*c) / denom))


def _membrane_ode(F_C3H6_feed: float, F_C3H8_feed: float,
                  P_H: float, P_L: float, A_mem: float,
                  Q_A_SI: float, alpha: float):
    """
    クロスフロー膜モジュールの ODE 積分。

    dF_C3H6/dA = −Q_A × (x × P_H − y_local × P_L)
    dF_C3H8/dA = −Q_B × ((1−x) × P_H − (1−y_local) × P_L)

    Returns
    -------
    F_ret_C3H6, F_ret_C3H8 : 非透過ガス流量 [mol/s]
    None, None              : 積分失敗時
    """
    gamma = P_L / P_H
    Q_B   = Q_A_SI / alpha

    def ode(A, F):
        fc  = max(F[0], 1e-12)
        fa  = max(F[1], 1e-12)
        x   = fc / (fc + fa)
        y   = _y_local(x, alpha, gamma)
        J_c = Q_A_SI * (x * P_H - y * P_L)
        J_a = Q_B    * ((1.0-x) * P_H - (1.0-y) * P_L)
        if J_c <= 0.0 or J_a <= 0.0:
            return [0.0, 0.0]
        return [-J_c, -J_a]

    try:
        sol = solve_ivp(
            ode,
            t_span=(0.0, A_mem),
            y0=[F_C3H6_feed, F_C3H8_feed],
            method='Radau',
            rtol=1e-5,
            atol=1e-8,
        )
    except Exception:
        return None, None

    if not sol.success:
        return None, None

    return float(sol.y[0, -1]), float(sol.y[1, -1])


# ---------------------------------------------------------------------------
# ユニット 5: 製品冷却器
# ---------------------------------------------------------------------------

def _condenser(F_perm_mols: float, y_C3H6: float,
               T_in: float, P_dist: float,
               fixed: MemFixedParams):
    """
    製品冷却器: 圧縮後の透過ガスを飽和液まで冷却・凝縮。

    Returns
    -------
    T_bp  : 出口（泡点）温度 [K]
    Q_cond_kW : 必要冷却量 [kW]
    A_cond    : 伝熱面積 [m²]
    """
    T_bp = bubble_point_T(P_dist, [y_C3H6, 1.0 - y_C3H6], _KEYS)

    H_gas_in  = _h_mol(T_in, P_dist, y_C3H6, 'vapor')
    H_liq_out = _h_mol(T_bp, P_dist, y_C3H6, 'liquid')
    Q_cond    = F_perm_mols * (H_gas_in - H_liq_out)   # [W]
    Q_cond_kW = Q_cond / 1e3

    # 向流熱交換器 LMTD
    # ガス側: T_in → T_bp,  冷却水側: T_cold_in → T_cold_out
    dT1 = T_in - fixed.T_cold_out     # ガス入口端
    dT2 = T_bp - fixed.T_cold_in      # ガス出口端（液出口端）
    lmtd = _lmtd(dT1, dT2)
    A_cond = Q_cond_kW / (fixed.U_cond * lmtd) if lmtd > 0 else float('nan')

    return T_bp, Q_cond_kW, A_cond


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_membrane_system(
    design: MemDesignVars,
    feed: MemFeedStream,
    fixed: MemFixedParams,
) -> MemSimulationResult:
    """
    膜分離システム全体をシミュレーションする。

    Parameters
    ----------
    design : MemDesignVars   設計変数 (P_H, P_L, A_mem)
    feed   : MemFeedStream   入力ストリーム（前段蒸留塔底液）
    fixed  : MemFixedParams  固定パラメータ

    Returns
    -------
    MemSimulationResult
    """
    # ---- 入力バリデーション ----
    if design.P_H <= design.P_L:
        warnings.warn("P_H <= P_L: 膜の駆動力がありません。")
        return _penalty_result()
    if design.A_mem <= 0 or design.P_H <= 0 or design.P_L <= 0:
        return _penalty_result()
    if feed.F_C3H6 < 0 or feed.F_C3H8 < 0:
        return _penalty_result()

    F_total_feed = feed.F_C3H6 + feed.F_C3H8
    if F_total_feed <= 0:
        return _penalty_result()

    z_C3H6_feed = feed.F_C3H6 / F_total_feed  # 供給液中 C3H6 分率

    # mol/s 換算（内部計算用）
    F_feed_mols = F_total_feed * 1000.0 / 3600.0   # [mol/s]
    Q_A_SI = fixed.Q_A_GPU * _GPU_SI                # [mol/(m²·s·Pa)]

    # ---- ユニット 1: 気化器 ----
    try:
        T_vap_out, Q_vap_kW, A_vap = _vaporizer(
            F_feed_mols, z_C3H6_feed, feed.T_in, feed.P_in, fixed
        )
    except Exception:
        return _penalty_result()
    if math.isnan(Q_vap_kW):
        return _penalty_result()

    # ---- ユニット 2: フィード圧縮機 ----
    try:
        T_feed_comp_out, W_feed_per_mol = compress_isentropic(
            T_vap_out, feed.P_in, design.P_H,
            [z_C3H6_feed, 1.0 - z_C3H6_feed], _KEYS,
            eta=fixed.eta_comp,
        )
    except Exception:
        return _penalty_result()
    W_feed_kW = F_feed_mols * W_feed_per_mol / 1e3   # [kW]

    # ---- ユニット 3: 膜モジュール ----
    F_C3H6_feed_mols = z_C3H6_feed * F_feed_mols
    F_C3H8_feed_mols = (1.0 - z_C3H6_feed) * F_feed_mols

    F_ret_C3H6_mols, F_ret_C3H8_mols = _membrane_ode(
        F_C3H6_feed_mols, F_C3H8_feed_mols,
        design.P_H, design.P_L, design.A_mem,
        Q_A_SI, fixed.alpha,
    )
    if F_ret_C3H6_mols is None:
        return _penalty_result()

    # 非透過・透過流量（mol/s）
    F_ret_C3H6_mols = max(F_ret_C3H6_mols, 0.0)
    F_ret_C3H8_mols = max(F_ret_C3H8_mols, 0.0)
    F_perm_C3H6_mols = F_C3H6_feed_mols - F_ret_C3H6_mols
    F_perm_C3H8_mols = F_C3H8_feed_mols - F_ret_C3H8_mols
    F_perm_total_mols = max(F_perm_C3H6_mols + F_perm_C3H8_mols, 1e-12)

    # 透過・非透過組成
    y_C3H6 = F_perm_C3H6_mols / F_perm_total_mols
    F_ret_total_mols = F_ret_C3H6_mols + F_ret_C3H8_mols
    x_ret_C3H6 = (F_ret_C3H6_mols / max(F_ret_total_mols, 1e-12))
    stage_cut = (F_perm_C3H6_mols + F_perm_C3H8_mols) / F_feed_mols

    # ---- ユニット 4: 製品圧縮機 ----
    # 膜は等温操作と仮定: 透過ガスは T_feed_comp_out, P_L で出る
    try:
        T_prod_comp_out, W_prod_per_mol = compress_isentropic(
            T_feed_comp_out, design.P_L, fixed.P_dist,
            [y_C3H6, 1.0 - y_C3H6], _KEYS,
            eta=fixed.eta_comp,
        )
    except Exception:
        return _penalty_result()
    W_prod_kW = F_perm_total_mols * W_prod_per_mol / 1e3   # [kW]

    # ---- ユニット 5: 製品冷却器 ----
    try:
        T_bp_perm, Q_cond_kW, A_cond = _condenser(
            F_perm_total_mols, y_C3H6,
            T_prod_comp_out, fixed.P_dist,
            fixed,
        )
    except Exception:
        return _penalty_result()

    # ---- kmol/h 換算（出力用）----
    to_kmolh = 3600.0 / 1000.0

    # ---- ゲージ圧変換 [barg] ----
    def _pg(P_pa: float) -> float:
        return P_pa / 1e5 - _ATM_BAR

    # ---- モジュール本数 ----
    n_modules = math.ceil(design.A_mem / fixed.A_per_module)

    return MemSimulationResult(
        retentate=MemRetentateStream(
            F_C3H6 = F_ret_C3H6_mols * to_kmolh,
            F_C3H8 = F_ret_C3H8_mols * to_kmolh,
            T_out  = T_feed_comp_out,
            P_out  = design.P_H,
        ),
        product=MemProductStream(
            F_C3H6 = F_perm_C3H6_mols * to_kmolh,
            F_C3H8 = F_perm_C3H8_mols * to_kmolh,
            T_out  = T_bp_perm,
            P_out  = fixed.P_dist,
        ),
        equipment=MemEquipmentData(
            A_vap     = A_vap,
            Pg_vap    = _pg(feed.P_in),
            Q_vap_kW  = Q_vap_kW,
            W_feed_kW = W_feed_kW,
            Pg_feed   = _pg(design.P_H),
            A_mem     = design.A_mem,
            n_modules = n_modules,
            Pg_mem    = _pg(design.P_H),
            W_prod_kW = W_prod_kW,
            Pg_prod   = _pg(fixed.P_dist),
            A_cond    = A_cond,
            Pg_cond   = _pg(fixed.P_dist),
            Q_cond_kW = Q_cond_kW,
        ),
        stage_cut   = float(np.clip(stage_cut,   0.0, 1.0)),
        perm_purity = float(np.clip(y_C3H6,      0.0, 1.0)),
        ret_purity  = float(np.clip(x_ret_C3H6,  0.0, 1.0)),
    )
