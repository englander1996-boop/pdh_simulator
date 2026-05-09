"""
簡易 FUG 蒸留塔モデル (共通エンジン)

設計判断 (2026-05-08):
  fake_column (split_fracs ベースのダミー) を、Fenske-Underwood-Gilliland
  (FUG) 法ベースの簡易モデルに置き換える。物理的により正しい温度・熱量を
  与え、BO による P_col / N_stages / N_feed / reflux_ratio の最適化を可能にする。

各塔の固有設定 (LK/HK 成分、デフォルト P/N/R) は units/separators/columnX/columnX.py
で定義する (Option C ハイブリッド構造)。

参考文献:
  [1] Fenske M.R., "Fractionation of straight-run Pennsylvania gasoline,"
      Ind. Eng. Chem. 24 (1932) 482-485. — N_min 最小段数式の出典
  [2] Underwood A.J.V., "Fractional distillation of multicomponent
      mixtures," Chem. Eng. Prog. 44 (8) (1948) 603-614. — R_min 最小還流比式
  [3] Eduljee H.E., "Equations replace Gilliland plot,"
      Hydrocarbon Processing 54 (9) (1975) 120. — Gilliland 相関 (式形)
  [4] Kirkbride C.G., "Process design procedure for multicomponent
      fractionators," Petroleum Refiner 23 (9) (1944) 321. — フィード段位置
  [5] プロセス設計授業資料 R08-3.pdf 付録 A
      Table A.1 (Sieve trays K1, K2, K3), 圧力係数 Fp, 据付間接費 1.18,
      CEPCI 補正など。

物性値の出典:
  - 大気圧沸点 T_b [K], 蒸発潜熱 λ [kJ/mol] : 既存 _T_BOIL_ATM, _LAMBDA_KJ
  - PR EOS パラメータ Tc, Pc, ω : src/config.py THERMO_DATA
    (n-Butane Z は化学工学便覧 改訂六版 表 1.3 物質No.181)
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from scipy.optimize import brentq

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.cost_calculator import (
    calc_cp0, calc_fp, calc_tray_capex_okuyen, calc_he_capex_okuyen,
)
from src.cost_parameters import (
    B1, B2, FM, CEPCI_BASE, CEPCI_CURRENT,
    USD_TO_JPY, PLANT_INDIRECT_FACTOR,
)
from src.component_data import cp_of, CP_DEFAULT
from src.config import THERMO_DATA
from src.eos import (
    z_factor, fugacity_coeff, residual_enthalpy,
    bubble_point_T, dew_point_T,
)
from src.utility_selector import (
    select_cooling_utility, select_heating_utility,
)
from flowsheet.heat_integration import StreamPhase, lookup_U, utility_phase


# ---------------------------------------------------------------------------
# 物性データ (大気圧沸点と蒸発潜熱、Clausius-Clapeyron 用)
# ---------------------------------------------------------------------------
_R_GAS = 8.314     # [J/(mol·K)]
_P_ATM = 101325.0  # [Pa]

# 各成分の大気圧沸点 [K] と蒸発潜熱 [kJ/mol]
# (旧版から流用、PDH の主要成分の代表値)
_T_BOIL_ATM = {
    'A': 231.1, 'B': 225.5, 'C':  20.3,
    'D': 169.4, 'E': 111.7, 'F': 184.6, 'Z': 272.7,
}
_LAMBDA_KJ = {
    'A': 19.1, 'B': 18.4, 'C': 0.9,
    'D': 13.5, 'E':  8.2, 'F': 14.7, 'Z': 22.0,
}
_T_BOIL_DEFAULT = 250.0
_LAMBDA_DEFAULT = 15.0

_T_COND_MIN = 313.15  # 40°C — 冷却水使用時の凝縮器最低温度
_PENALTY = 1e9

# PR 泡点フラッシュの探索範囲 [K]。これを外れた解は CC へフォールバック。
# 下限: H2/CH4 が混じる x_top で「全成分の和」泡点が cryogenic に張り付くケースを除外。
# 上限: 探索失敗時に hi で打ち切られて非物理になるのを除外。
_T_BUBBLE_MIN = 200.0
_T_BUBBLE_MAX = 600.0


# ===========================================================================
# データクラス: 入力・出力
# ===========================================================================

@dataclass(frozen=True)
class ColumnTunables:
    """蒸留塔の最適化対象パラメータ (P/N/R)。

    LK/HK・回収率・q・K_method といった「塔別の物理セマンティクス」は
    column1/2/3.py 側で固定。本データクラスは BO や exp で振る部分のみ保持。

    FlowsheetDesignVars.dist1/dist2/dist3 で使う。
    """
    P_col:        float    # 塔操作圧力 [Pa]
    N_stages:     int      # 理論段数
    N_feed:       int      # フィード段位置 (現状 Kirkbride 推奨値が優先、本値は記録用)
    reflux_ratio: float    # 還流比 R = L/D


# 分流 (partial condenser) で凝縮させない成分セット
# 設計判断 (2026-05-09): H2 (Tc=33K) と CH4 (Tc=190K) は工業的に到達可能な
# 冷媒温度 (-100°C エチレン = 173K) でも凝縮できないため、partial condenser では
# 凝縮させず vapor distillate でそのまま下流へパススルーする。
# C2H4 (Tc=282K) 以上は propylene/ethylene 冷媒で凝縮可能なので reflux 寄与あり。
NON_CONDENSABLE_COMPS = frozenset(['C', 'E'])  # H2, CH4


@dataclass(frozen=True)
class DistDesignVars:
    """蒸留塔設計変数 (BO 探索対象 + 仕様)。

    BO で最適化する変数:
      P_col, N_stages, N_feed, reflux_ratio

    仕様 (塔別ラッパーで固定):
      LK, HK, recovery_LK_top, recovery_HK_bot, K_method, q
    """
    # ---- BO 設計変数 ----
    P_col:            float            # 塔操作圧力 [Pa]
    N_stages:         int              # 理論段数
    N_feed:           int              # フィード段位置 (1=塔頂, N_stages=塔底)
    reflux_ratio:     float            # 還流比 R = L/D
    # ---- 分離仕様 ----
    LK:               str              # 軽キー成分 ('A'-'F','Z')
    HK:               str              # 重キー成分
    recovery_LK_top:  float = 0.99     # LK 塔頂回収率
    recovery_HK_bot:  float = 0.99     # HK 塔底回収率
    # ---- 物性計算オプション ----
    K_method:         str   = 'pr'     # 'pr' (PR EOS) または 'cc' (Clausius-Clapeyron)
    q:                float = 1.0      # フィード状態 (1=飽和液, 0=飽和気)
    # 設計判断 (2026-05-09): partial condenser フラグ。
    # True のとき塔頂は分流型で、軽質ガス (NON_CONDENSABLE_COMPS = H2/CH4) は
    # 凝縮させず vapor distillate でパススルー、reflux 用に重質分のみ凝縮する。
    # Dist2 (脱エタン塔) のように H2/CH4 主体の塔頂で必須。Dist1/Dist3 (C3 主体) は False。
    partial_condenser: bool = False


@dataclass(frozen=True)
class DistFixedParams:
    """蒸留塔固定パラメータ (寸法・流体力学、contest §4-2 準拠)。

    出典: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-2 蒸留塔
      段塔仮定、塔径はフラッディングを起こさない許容蒸気質量速度 G* で決定:
        G* = SF · K · √(ρ_v · (ρ_l - ρ_v))   [kg/(m²·s)]
      シーブトレー仮定で SF=0.8、K=0.05 m/s、段間隔 0.6m を採用。
      段効率 80% で実段数を算出、塔頂 2m + 塔底 4m を加算。
    """
    tray_spacing_m:  float = 0.6   # トレイ間隔 [m]、contest §4-2
    top_section_m:   float = 2.0   # 塔頂付加長さ [m] (還流供給・気液分離)、contest §4-2
    bot_section_m:   float = 4.0   # 塔底付加長さ [m] (液ホールドアップ)、contest §4-2
    tray_efficiency: float = 0.8   # 段効率 [-]、contest §4-2 (実段数 = 理論段数 / η)
    SF:              float = 0.8   # G* の系補正係数 [-]、contest §4-2
    K_factor:        float = 0.05  # G* の許容蒸気速度係数 [m/s]、contest §4-2


@dataclass
class DistEquipment:
    """蒸留塔の装置・診断情報。"""
    # ---- 塔寸法 ----
    D_col:               float          # 塔径 [m]
    H_col:               float          # 塔高さ [m]
    V_col:               float          # 塔体積 [m³]
    # ---- CAPEX 内訳 ----
    CAPEX_vessel:        float          # 塔本体 [億円]
    CAPEX_trays:         float          # トレイ [億円]
    CAPEX:               float          # 合計 [億円] (vessel + trays + cond + reb)
    # ---- 熱量 ----
    Q_cond:              float          # コンデンサ熱量 [kW] (放出、正値)
    Q_reb:               float          # リボイラ熱量 [kW]   (吸熱、正値)
    Q_feed_preheat_kW:   float = 0.0    # フィード予熱 [kW]
    # ---- FUG 診断 ----
    N_min:               float = 0.0    # 最小段数 (Fenske)
    R_min:               float = 0.0    # 最小還流比 (Underwood)
    N_feed_kirkbride:    int   = 0      # Kirkbride 推奨フィード段
    feasible:            bool  = True   # N >= N_min かつ R >= R_min
    message:             str   = ""
    # ---- Condenser / Reboiler 装置 (2026-05-09 追加) ----
    A_cond_m2:               float = 0.0    # コンデンサ伝熱面積 [m²]
    A_reb_m2:                float = 0.0    # リボイラ伝熱面積 [m²]
    CAPEX_cond:              float = 0.0    # コンデンサ HE CAPEX [億円]
    CAPEX_reb:               float = 0.0    # リボイラ HE CAPEX [億円]
    cond_utility_name:       str   = ""     # 例 "冷却水≒", "プロピレン冷媒-25C≒"
    cond_utility_jpy_per_GJ: float = 0.0
    reb_utility_name:        str   = ""     # 例 "LP Steam≒", "MP Steam≒"
    reb_utility_jpy_per_GJ:  float = 0.0
    T_top:                   float = 0.0    # 塔頂温度 [K] (clamp なし)
    T_bot:                   float = 0.0    # 塔底温度 [K]


@dataclass
class DistResult:
    """蒸留塔シミュレーション結果。"""
    top:       ProcessStream
    bottom:    ProcessStream
    equipment: DistEquipment


# ===========================================================================
# 物性ヘルパー (温度・K 値・α)
# ===========================================================================

def _boil_cc(T_b_atm: float, lam_kJ: float, P_col: float) -> float:
    """Clausius-Clapeyron 式で圧力 P_col [Pa] における沸点 [K] を返す。

    ln(P2/P1) = ΔHvap/R × (1/T1 - 1/T2)
    → T2 = 1 / (1/T1 - R/ΔHvap × ln(P2/P1))
    """
    lam_J = lam_kJ * 1000.0
    inv_T2 = 1.0 / T_b_atm - _R_GAS / lam_J * math.log(P_col / _P_ATM)
    if inv_T2 <= 0.0:
        return T_b_atm * 5.0
    return 1.0 / inv_T2


def _weighted_boil(F_dict: Dict[str, float], P_col: float) -> float:
    """流量加重平均の沸点 [K] (Clausius-Clapeyron で圧力補正)。"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return 298.15
    t = sum(
        F_dict.get(k, 0.0) * _boil_cc(
            _T_BOIL_ATM.get(k, _T_BOIL_DEFAULT),
            _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT),
            P_col,
        )
        for k in F_dict
    ) / F_total
    return t


def _weighted_lambda(F_dict: Dict[str, float]) -> float:
    """流量加重平均蒸発潜熱 [kJ/mol]。"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return _LAMBDA_DEFAULT
    return sum(
        F_dict.get(k, 0.0) * _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT)
        for k in F_dict
    ) / F_total


def _K_cc(comp: str, T: float, P_col: float) -> float:
    """Clausius-Clapeyron で K_i = P_sat,i(T) / P_col。"""
    T_b = _T_BOIL_ATM.get(comp, _T_BOIL_DEFAULT)
    lam = _LAMBDA_KJ.get(comp, _LAMBDA_DEFAULT) * 1000.0  # J/mol
    # P_sat = P_atm × exp(λ/R × (1/T_b - 1/T))
    P_sat = _P_ATM * math.exp(lam / _R_GAS * (1.0 / T_b - 1.0 / T))
    return P_sat / P_col


def _K_pr(comp_idx: int, comps: List[str], x: List[float],
          T: float, P_col: float) -> float:
    """PR EOS で K_i = phi_i^L / phi_i^V を計算。

    両相 (液・気) でフガシティ係数を計算し比を取る。
    """
    try:
        Z_L = z_factor(T, P_col, x, comps, 'liquid')
        Z_V = z_factor(T, P_col, x, comps, 'vapor')
        phi_L = fugacity_coeff(comp_idx, T, P_col, x, comps, Z_L)
        phi_V = fugacity_coeff(comp_idx, T, P_col, x, comps, Z_V)
        if phi_V <= 1e-30:
            return float('inf')
        return phi_L / phi_V
    except Exception:
        # 範囲外・収束失敗時は CC でフォールバック
        return _K_cc(comps[comp_idx], T, P_col)


def _K_dict(comps: List[str], x: List[float], T: float, P_col: float,
            method: str) -> Dict[str, float]:
    """全成分の K 値を計算して辞書で返す。"""
    if method == 'pr':
        # PR EOS は混合系。Tc/Pc/ω のない成分は CC にフォールバック
        K = {}
        for i, c in enumerate(comps):
            p = THERMO_DATA.get(c)
            if p is None or math.isnan(p.Tc):
                K[c] = _K_cc(c, T, P_col)
            else:
                K[c] = _K_pr(i, comps, x, T, P_col)
        return K
    # 'cc' モード
    return {c: _K_cc(c, T, P_col) for c in comps}


def _alpha_dict(K: Dict[str, float], HK: str) -> Dict[str, float]:
    """HK 基準の相対揮発度 α_i = K_i / K_HK。"""
    K_HK = K.get(HK, 1.0)
    if K_HK <= 1e-30:
        K_HK = 1e-30
    return {c: K[c] / K_HK for c in K}


def _bubble_T_K(
    x_list: List[float],
    comps:  List[str],
    P_col:  float,
    method: str,
) -> Tuple[float, List[float]]:
    """組成 x (mol fraction) の泡点温度と各成分 K 値を返す。

    method='pr' のとき: PR EOS で泡点フラッシュ (eos.bubble_point_T) を実行し、
    収束 T で K_i = φ_L/φ_V を取る。FUG 法で「平均 T で K=phi_L/phi_V」を計算
    すると、その T-P-x が単相領域に入り Z 根が 1 本で K_i ≈ 1 になる病理を
    回避できる (泡点では飽和包絡線上なので必ず 2 根存在)。

    PR が _T_BUBBLE_MIN..MAX を外れる解を返したり例外を出した場合 (例: H2 や CH4
    が x_top の主成分で泡点が cryogenic に張り付くケース) は CC にフォールバック。
    成分単位の失敗は当該成分のみ CC へ。

    method='cc' のとき: 加重平均沸点 + 純成分 P_sat ベース (旧来の挙動)。
    """
    s = sum(max(xi, 0.0) for xi in x_list)
    if s <= 1e-12:
        T_def = 298.15
        return T_def, [_K_cc(c, T_def, P_col) for c in comps]
    x_norm = [max(xi, 0.0) / s for xi in x_list]

    if method == 'cc':
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    # ---- PR 経路 ----
    pr_capable = [
        c in THERMO_DATA and not math.isnan(THERMO_DATA[c].Tc)
        for c in comps
    ]
    if sum(pr_capable) < 2:
        # PR パラメータを持つ成分が 2 つ未満 → 飽和挙動を計算できない。CC へ。
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            T_bub = bubble_point_T(P_col, x_norm, comps,
                                   T_lo=_T_BUBBLE_MIN, T_hi=_T_BUBBLE_MAX)
        if (math.isnan(T_bub)
            or T_bub <= _T_BUBBLE_MIN + 1.0
            or T_bub >= _T_BUBBLE_MAX - 1.0):
            raise ValueError("bubble T pinned to bound or nan")
    except Exception:
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    # 収束 T で K = phi_L/phi_V を成分別に取得 (個別失敗は CC へ)。
    K = []
    for i, c in enumerate(comps):
        if not pr_capable[i]:
            K.append(_K_cc(c, T_bub, P_col))
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                K_i = _K_pr(i, comps, x_norm, T_bub, P_col)
            if math.isfinite(K_i) and K_i > 0.0:
                K.append(K_i)
            else:
                K.append(_K_cc(c, T_bub, P_col))
        except Exception:
            K.append(_K_cc(c, T_bub, P_col))

    # 単相 root 検出: K の spread が小さいと単相領域に張り付いており物理的でない。
    # 例: x に H2/CH4 (Tc << 室温) が混じる x_top で bubble_point_T が cryogenic な
    # 偽解に収束し、cubic が 1 根しかないため Z_L=Z_V → φ_L/φ_V=1 → K_i ≈ 1 全成分。
    # この状態では α 計算が破綻するため CC に退避する。
    K_pos = [k for k in K if k > 0.0]
    if K_pos:
        K_max = max(K_pos)
        K_min = min(K_pos)
        if K_max / max(K_min, 1e-30) < 1.1:
            F = {comps[i]: x_norm[i] for i in range(len(comps))}
            T = _weighted_boil(F, P_col)
            return T, [_K_cc(c, T, P_col) for c in comps]
    return T_bub, K


# ===========================================================================
# FUG 計算式
# ===========================================================================

def _fenske_N_min(alpha_LK: float,
                  recovery_LK_top: float,
                  recovery_HK_bot: float) -> float:
    """Fenske 式: 最小段数 N_min。

    N_min = log[(r_LK/(1-r_LK)) × (r_HK/(1-r_HK))] / log(α_LK,HK)

    r_LK : LK の塔頂回収率
    r_HK : HK の塔底回収率
    """
    if alpha_LK <= 1.0 + 1e-9:
        # α≈1 で分離不能
        return float('inf')
    rLK = max(min(recovery_LK_top, 1.0 - 1e-9), 1e-9)
    rHK = max(min(recovery_HK_bot, 1.0 - 1e-9), 1e-9)
    ratio = (rLK / (1.0 - rLK)) * (rHK / (1.0 - rHK))
    if ratio <= 1.0:
        return 0.0
    return math.log(ratio) / math.log(alpha_LK)


def _underwood_R_min(alpha: Dict[str, float], z: Dict[str, float],
                     x_top: Dict[str, float], q: float,
                     LK: str, HK: str) -> float:
    """Underwood 式: 最小還流比 R_min。

    手順:
      1. Σ α_i × z_i / (α_i - θ) = 1 - q を θ について解く
         (θ は α_HK < θ < α_LK の範囲に存在)
      2. R_min = Σ α_i × x_top_i / (α_i - θ) - 1
    """
    alpha_LK = alpha[LK]
    alpha_HK = alpha[HK]
    if alpha_LK <= alpha_HK + 1e-9:
        # α 序列が逆転している → R_min 不定
        return float('inf')

    def eq1(theta: float) -> float:
        return sum(alpha[c] * z[c] / (alpha[c] - theta) for c in alpha) - (1.0 - q)

    # θ は α_HK と α_LK の間で 1 個の根
    eps = 1e-6
    lo = alpha_HK + eps * (alpha_LK - alpha_HK)
    hi = alpha_LK - eps * (alpha_LK - alpha_HK)
    try:
        theta = brentq(eq1, lo, hi, xtol=1e-7, maxiter=200)
    except (ValueError, RuntimeError):
        return float('inf')

    # R_min を塔頂組成から計算
    R_min = sum(alpha[c] * x_top[c] / (alpha[c] - theta) for c in alpha) - 1.0
    return max(0.0, R_min)


def _gilliland_eduljee(R: float, R_min: float, N_min: float) -> float:
    """Gilliland 相関 (Eduljee 1975 形式) で実段数 N を計算。

    X = (R - R_min) / (R + 1)
    Y = 0.75 × (1 - X^0.5668)
    N = (Y + N_min) / (1 - Y)
    """
    if R <= R_min + 1e-9 or R + 1.0 <= 0:
        return float('inf')
    X = (R - R_min) / (R + 1.0)
    if X <= 0 or X >= 1:
        return float('inf')
    Y = 0.75 * (1.0 - X ** 0.5668)
    if Y >= 1.0 - 1e-9:
        return float('inf')
    return (Y + N_min) / (1.0 - Y)


def _kirkbride_feed_stage(N: int, F_top: float, F_bot: float,
                          z: Dict[str, float],
                          x_top: Dict[str, float], x_bot: Dict[str, float],
                          LK: str, HK: str) -> int:
    """Kirkbride 式で推奨フィード段位置を計算。

    log10(N_above/N_below) = 0.206 × log10[
        (z_HK/z_LK) × (x_LK_bot/x_HK_top)^2 × (B/D)
    ]

    返り値: 塔頂から数えた段番号 (1=塔頂、N=塔底)。
    """
    z_LK = max(z.get(LK, 1e-9), 1e-9)
    z_HK = max(z.get(HK, 1e-9), 1e-9)
    x_LK_bot = max(x_bot.get(LK, 1e-9), 1e-9)
    x_HK_top = max(x_top.get(HK, 1e-9), 1e-9)
    if F_top <= 0:
        return 1
    arg = (z_HK / z_LK) * (x_LK_bot / x_HK_top) ** 2 * (F_bot / F_top)
    if arg <= 0:
        return max(1, N // 2)
    ratio = 10.0 ** (0.206 * math.log10(arg))   # N_above / N_below
    # N_above + N_below = N, N_above / N_below = ratio → N_above = N × ratio/(1+ratio)
    N_above = int(round(N * ratio / (1.0 + ratio)))
    N_above = max(1, min(N - 1, N_above))
    return N_above   # 塔頂から N_above 段目がフィード段


# ===========================================================================
# 物質収支 (Fenske による全成分分配)
# ===========================================================================

def _split_streams(F_in: Dict[str, float],
                   alpha: Dict[str, float],
                   N_min: float, HK: str,
                   recovery_LK_top: float,
                   recovery_HK_bot: float,
                   LK: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Fenske 関係を使って全成分の塔頂・塔底分配を計算。

    アプローチ:
      Key 成分は recovery で固定:
        F_top[LK] = F_in[LK] × recovery_LK_top
        F_top[HK] = F_in[HK] × (1 - recovery_HK_bot)
      Non-key 成分は Fenske の比から個別に分配:
        (x_top/x_bot)_i = α_i^N_min × (x_top/x_bot)_HK_ref
      → 個別物質収支で塔頂分率を導出。

    HK 基準量 (HK の塔頂/塔底比):
      ratio_HK = recovery_HK_bot ベースで決定
                = (F_HK_top / F_HK_bot) = (1 - rec_HK) / rec_HK
    """
    rec_HK = max(min(recovery_HK_bot, 1.0 - 1e-9), 1e-9)
    ratio_HK = (1.0 - rec_HK) / rec_HK   # HK の (top/bot) 流量比

    F_top: Dict[str, float] = {}
    F_bot: Dict[str, float] = {}
    for c, F_c in F_in.items():
        if F_c <= 0:
            F_top[c] = 0.0
            F_bot[c] = 0.0
            continue
        if c == LK:
            frac_top = recovery_LK_top
        elif c == HK:
            frac_top = 1.0 - rec_HK
        else:
            # Fenske: (top/bot)_c = α_c^N_min × (top/bot)_HK
            a_c = alpha.get(c, 1.0)
            if N_min == float('inf') or N_min > 100:
                # α > 1 なら塔頂、α < 1 なら塔底に近づく
                if a_c > 1.0:
                    frac_top = 1.0 - 1e-9
                elif a_c < 1.0:
                    frac_top = 1e-9
                else:
                    frac_top = 0.5
            else:
                ratio_c = (a_c ** N_min) * ratio_HK
                # frac_top = ratio_c / (1 + ratio_c)
                if ratio_c > 1e30:
                    frac_top = 1.0 - 1e-9
                else:
                    frac_top = ratio_c / (1.0 + ratio_c)
        F_top[c] = F_c * frac_top
        F_bot[c] = F_c * (1.0 - frac_top)
    return F_top, F_bot


# ===========================================================================
# 塔本体 CAPEX (Vessel: 既存ロジック踏襲)
# ===========================================================================

def _vessel_capex_okuyen(V_m3: float, P_pa: float, D_m: float) -> float:
    """塔本体 (vertical vessel) の CAPEX [億円]。

    出典: プロセス設計授業資料 R08-3.pdf 付録 A
      C_p0 = 10^(K1 + K2*log10(V) + K3*(log10(V))^2)
      F_p  = max((Pg+1)*D / (10.71 - 0.00756(Pg+1)) + 0.5, 1)  (Pg > -0.5 bar)
      C_BM = C_p0 × (B1 + B2 × Fp × FM)
      C_TM = 1.18 × C_BM × CEPCI 補正
    """
    P_bar = P_pa / 1.0e5 - 1.01325
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        cp0 = calc_cp0(V_m3)
    fp  = calc_fp(P_bar, D_m)
    fbm = B1 + B2 * fp * FM
    cbm = cp0 * fbm
    usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * usd * USD_TO_JPY / 1.0e8


# ===========================================================================
# メイン関数
# ===========================================================================

def simulate_distillation_column(
    design: DistDesignVars,
    feed:   ProcessStream,
    fixed:  Optional[DistFixedParams] = None,
) -> DistResult:
    """FUG ベースの蒸留塔シミュレーション。

    手順:
      1. 動作温度推定 (T_top, T_bot を Clausius-Clapeyron で初期推定)
      2. K 値・α 計算 (PR EOS or CC)
      3. Fenske: N_min
      4. 物質収支: 全成分の塔頂・塔底分配
      5. Underwood: R_min (塔頂組成必要のため物質収支後)
      6. Gilliland: feasibility 確認 (N >= N_min, R >= R_min)
      7. Kirkbride: 推奨フィード段
      8. Q_cond, Q_reb, Q_feed_preheat 計算
      9. 塔径・塔高計算
     10. CAPEX (Vessel + Trays)
     11. infeasible なら ペナルティ返却
    """
    if fixed is None:
        fixed = DistFixedParams()

    comps = list(feed.F_in.keys())
    F_total = sum(max(F, 0.0) for F in feed.F_in.values())
    if F_total <= 0.0:
        return _penalty_result(design, "feed flow ≤ 0")

    # フィード組成 z (mol fraction)
    z = {c: max(feed.F_in[c], 0.0) / F_total for c in comps}

    # ---- Step 1: 動作温度・K の初期推定 (CC、x=z で粗推定) ----
    # 塔頂は LK 主体、塔底は HK 主体と仮定して純成分沸点で T 範囲を見立てる。
    T_top = _boil_cc(_T_BOIL_ATM.get(design.LK, _T_BOIL_DEFAULT),
                     _LAMBDA_KJ.get(design.LK, _LAMBDA_DEFAULT),
                     design.P_col)
    T_bot = _boil_cc(_T_BOIL_ATM.get(design.HK, _T_BOIL_DEFAULT),
                     _LAMBDA_KJ.get(design.HK, _LAMBDA_DEFAULT),
                     design.P_col)
    # 設計判断 (2026-05-09): T_top の _T_COND_MIN clamp を解除。
    # 旧版は「冷却水到達下限 40°C」で上書きしていたが、これだと Dist2 の x_top
    # (H2/CH4 主体) で実際は cryogenic な凝縮器が必要なケースを隠蔽していた。
    # 本版は実温度を保持し、utility_selector で適切な冷媒を選ぶ。

    # CC で初期 α・初期 split を 1 回作る (PR 反復の出発点)。
    K_init = _K_dict(comps, [z[c] for c in comps],
                     0.5 * (T_top + T_bot), design.P_col, 'cc')
    alpha = _alpha_dict(K_init, design.HK)
    N_min = _fenske_N_min(
        alpha.get(design.LK, 1.0),
        design.recovery_LK_top, design.recovery_HK_bot,
    )
    F_top, F_bot = _split_streams(
        feed.F_in, alpha, N_min, design.HK,
        design.recovery_LK_top, design.recovery_HK_bot, design.LK,
    )

    # ---- Step 2-4: 塔頂/塔底それぞれの泡点で α を評価して反復 ----
    # 設計判断 (2026-05-09): 旧版は K_avg を「塔平均 T で x=z」で計算していたが、
    # この T-P-x は単相領域に入りやすく、PR の 3 次方程式が 1 根しか返さない結果
    # K_i = phi_L/phi_V ≈ 1 → α ≈ 1 → N_min = ∞ で全塔 infeasible になる病理が
    # あった (column1.py 旧コメント参照)。本版では塔頂/塔底それぞれの組成で
    # 泡点フラッシュを実行し、必ず 2 相が成立する点で K を取る。
    # α_geom = sqrt(α_top × α_bot) は FUG の標準的な平均化 (Sinnott 17.5.3)。
    for _ in range(5):
        F_top_total_it = sum(F_top.values())
        F_bot_total_it = sum(F_bot.values())
        if F_top_total_it <= 0.0 or F_bot_total_it <= 0.0:
            return _penalty_result(design, "split failed (top or bottom flow ≤ 0)")
        x_top_list = [F_top[c] / F_top_total_it for c in comps]
        x_bot_list = [F_bot[c] / F_bot_total_it for c in comps]

        T_top_new, K_top_list = _bubble_T_K(
            x_top_list, comps, design.P_col, design.K_method,
        )
        T_bot_new, K_bot_list = _bubble_T_K(
            x_bot_list, comps, design.P_col, design.K_method,
        )
        # 設計判断 (2026-05-09): _T_COND_MIN clamp 解除。実温度を保持して下流で
        # utility_selector に適切な冷媒を選ばせる。

        HK_idx = comps.index(design.HK)
        K_top_HK = max(K_top_list[HK_idx], 1e-30)
        K_bot_HK = max(K_bot_list[HK_idx], 1e-30)
        alpha_new: Dict[str, float] = {}
        for i, c in enumerate(comps):
            a_top = K_top_list[i] / K_top_HK
            a_bot = K_bot_list[i] / K_bot_HK
            alpha_new[c] = math.sqrt(max(a_top * a_bot, 1e-30))

        N_min_new = _fenske_N_min(
            alpha_new.get(design.LK, 1.0),
            design.recovery_LK_top, design.recovery_HK_bot,
        )
        F_top_new, F_bot_new = _split_streams(
            feed.F_in, alpha_new, N_min_new, design.HK,
            design.recovery_LK_top, design.recovery_HK_bot, design.LK,
        )

        # 収束判定: 塔頂流量の相対変化 (フィードに対して 0.1% 未満で打ち切り)
        diff = max(
            (abs(F_top_new[c] - F_top[c]) / max(feed.F_in.get(c, 0.0), 1e-9))
            for c in comps if feed.F_in.get(c, 0.0) > 0.0
        )
        F_top, F_bot = F_top_new, F_bot_new
        alpha, N_min = alpha_new, N_min_new
        T_top, T_bot = T_top_new, T_bot_new
        if diff < 1e-3:
            break

    # ---- 流量・組成 (mol fraction) ----
    F_top_total = sum(F_top.values())
    F_bot_total = sum(F_bot.values())
    if F_top_total <= 0 or F_bot_total <= 0:
        return _penalty_result(design, "split failed (top or bottom flow ≤ 0)")
    x_top = {c: F_top[c] / F_top_total for c in comps}
    x_bot = {c: F_bot[c] / F_bot_total for c in comps}

    # ---- Step 5: Underwood で R_min ----
    R_min = _underwood_R_min(alpha, z, x_top, design.q, design.LK, design.HK)

    # ---- Step 6: feasibility ----
    feasible = True
    msg = ""
    if N_min == float('inf'):
        feasible = False
        msg = f"alpha_LK ≈ 1 で分離不能"
    elif design.N_stages < N_min:
        feasible = False
        msg = (f"N_stages={design.N_stages} < N_min={N_min:.1f} "
               f"(LK={design.LK}, HK={design.HK})")
    elif R_min == float('inf'):
        feasible = False
        msg = "Underwood 収束失敗 (alpha 序列が逆)"
    elif design.reflux_ratio < R_min:
        feasible = False
        msg = (f"R={design.reflux_ratio:.2f} < R_min={R_min:.2f} "
               f"(LK={design.LK}, HK={design.HK})")

    if not feasible:
        return _penalty_result(design, msg, N_min=N_min, R_min=R_min)

    # ---- Step 7: Kirkbride で推奨 N_feed ----
    N_feed_kirkbride = _kirkbride_feed_stage(
        design.N_stages, F_top_total, F_bot_total,
        z, x_top, x_bot, design.LK, design.HK,
    )

    # ---- Step 8: 熱量計算 ----
    # 設計判断 (2026-05-09): partial condenser 分岐。
    # Total condenser: V = (R+1)D を全量凝縮 → Q = (R+1) × D × λ_avg、T_cond = bubble of x_top
    # Partial condenser: 軽質ガス (NON_CONDENSABLE) は気相パススルー。
    #                   reflux に使う condensable 成分のみ凝縮 → Q = R × D_C × λ_C、
    #                   T_cond = bubble of x_top_condensable (実機では暖かい温度で済む)
    if design.partial_condenser:
        F_top_C = {c: v for c, v in F_top.items()
                   if c not in NON_CONDENSABLE_COMPS and v > 0.0}
        F_C_total = sum(F_top_C.values())
        if F_C_total <= 0:
            return _penalty_result(
                design,
                "partial condenser: 凝縮可能成分が無い (全成分が H2/CH4)",
                N_min=N_min, R_min=R_min,
            )
        x_C_list = [F_top_C.get(c, 0.0) / F_C_total for c in comps]
        T_cond, _ = _bubble_T_K(x_C_list, comps, design.P_col, design.K_method)
        lam_top = _weighted_lambda(F_top_C)        # 凝縮分の λ
        # reflux liquid のみ凝縮 (vapor distillate は凝縮せず): Q = R × D_C × λ_C
        Q_cond_kW = (design.reflux_ratio * F_C_total
                     * lam_top * 1000.0 / 3600.0)
    else:
        # Total condenser (Dist1/Dist3): V = (R+1)D を全量凝縮
        T_cond  = T_top
        lam_top = _weighted_lambda(F_top)
        Q_cond_kW = (F_top_total * (design.reflux_ratio + 1.0)
                     * lam_top * 1000.0 / 3600.0)

    # リボイラ熱量 Q_reb = V' × λ_bot
    # 設計判断 (2026-05-09): 旧版 Q_reb = 1.05 × Q_cond は q=1 (飽和液フィード)
    # 専用の経験則で、q=0 (気相フィード) の Dist2 では物理的に不正確だった。
    # 標準 MESH 式 V' = (R+1)D - (1-q)F (q=0 で V' は (R+1)D - F に減る) に置換。
    # q=1 (Dist1, Dist3): V' = (R+1)D となり旧式とほぼ同等。
    # q=0 (Dist2): V' = (R+1)D - F、partial condenser 後の Q_cond と独立に決まる。
    F_feed_total  = sum(max(F, 0.0) for F in feed.F_in.values())
    V_prime_kmolh = (F_top_total * (design.reflux_ratio + 1.0)
                    - (1.0 - design.q) * F_feed_total)
    if V_prime_kmolh > 0:
        lam_bot = _weighted_lambda(F_bot)
        Q_reb_kW = V_prime_kmolh * lam_bot * 1000.0 / 3600.0
    else:
        # V' ≤ 0 は通常 q=0 で feed が大きすぎるとき。塔として動作不能なので
        # 物理的には infeasible だが、ここでは保守的に Q_cond と同等値で埋める。
        Q_reb_kW = Q_cond_kW

    # ProcessStream の T_top には実際の condenser 出口温度を入れる
    # (partial の場合は T_cond > T_top_bubble、total の場合は T_cond = T_top_bubble)
    T_top = T_cond

    # フィード予熱 (顕熱、feed.T_in → T_bot、簡易液仮定)
    Q_feed_preheat_kW = 0.0
    if feed.T_in < T_bot:
        for c, F_kmol_h in feed.F_in.items():
            if F_kmol_h <= 0:
                continue
            cp = cp_of(c) if cp_of(c) > 0 else CP_DEFAULT  # [J/(mol·K)]
            F_mol_s = F_kmol_h * 1000.0 / 3600.0
            Q_feed_preheat_kW += F_mol_s * cp * (T_bot - feed.T_in) / 1000.0

    # ---- Step 9: 塔径・塔高 (contest §4-2) ----
    # 設計判断 (2026-05-09): 旧版は u_vapor=0.3 m/s 単一固定で塔径を出していたが、
    # contest §4-2 仕様 (G* = SF·K·√(ρ_v(ρ_l-ρ_v))) に置換して物理ベースに。
    # 段塔・シーブトレー仮定。塔径は塔頂条件で評価 (蒸気量最大、最も flooding しやすい)。
    from src.component_data import MW, liquid_density_mix
    MW_top_avg = sum(F_top.get(c, 0) * MW.get(c, 50.0) for c in F_top) / max(F_top_total, 1e-9)
    # 蒸気密度: 理想気体近似 ρ_v = P·MW/(R·T) [kg/m³]
    rho_v = design.P_col * (MW_top_avg * 1e-3) / (_R_GAS * T_top)
    # 液密度: 既存 liquid_density_mix (kg/m³)、x_top を流量代用
    rho_l = liquid_density_mix(F_top)
    if rho_l <= rho_v + 1.0:
        # 物性異常 (蒸気と液の密度差が小さすぎる) → 保守側で D を大きめに
        rho_l = rho_v + 100.0
    # G* = SF · K · √(ρ_v · (ρ_l - ρ_v))   [kg/(m²·s)]
    G_star = fixed.SF * fixed.K_factor * math.sqrt(rho_v * (rho_l - rho_v))
    if G_star <= 1e-6:
        return _penalty_result(design, f"G*={G_star:.3e} <= 0 (塔径計算不能)",
                               N_min=N_min, R_min=R_min)
    # 蒸気の質量流量 [kg/s] = (R+1)·D·MW_v / 3600
    mass_flow_vap = F_top_total * (design.reflux_ratio + 1.0) * MW_top_avg / 3600.0
    A_cross = mass_flow_vap / G_star            # [m²]
    D_col = math.sqrt(4.0 * A_cross / math.pi)
    D_col = max(D_col, 0.3)

    # 塔高: 実段数 = 理論段 / 段効率、塔頂・塔底 section を加算
    N_real = math.ceil(design.N_stages / fixed.tray_efficiency)
    H_col = (N_real * fixed.tray_spacing_m
             + fixed.top_section_m + fixed.bot_section_m)
    H_col = max(H_col, 5.0)
    V_col = math.pi / 4.0 * D_col ** 2 * H_col

    # ---- Step 10a: Condenser / Reboiler の utility 選択と装置サイジング ----
    # 設計判断 (2026-05-09): 旧版は塔本体 (vessel + trays) しか CAPEX に含めず、
    # コンデンサも一律 COOLING_WATER 単価で OPEX 計上していた。実機では condenser・
    # reboiler は独立 HE で大きな CAPEX を持ち、温度依存で utility が変わる
    # (Dist2 のように T_top 氷点下ならエチレン冷媒が必要)。
    # ここで utility_selector に T_top/T_bot を渡して適切な tier を自動選択する。
    #
    # T_top が冷媒到達下限を下回るケース (Dist2 で x_top に H2/CH4 主体の場合) は、
    # 実機では partial condenser で対応するが本シミュレータは total condenser モデル。
    # 最低冷媒 (-100°C エチレン) で押し切るコストを保守見積もりとして採用する。
    _T_COOLING_FLOOR = 173.15 + 10.0 + 1.0   # -89.0°C: -100°C 冷媒 + 10K approach + 1K 余裕
    T_top_for_util = max(T_top, _T_COOLING_FLOOR)
    try:
        cond_utility = select_cooling_utility(T_top_for_util, mode='continuous')
    except ValueError:
        return _penalty_result(
            design,
            f"condenser T_top={T_top-273.15:.1f}°C は冷媒で到達不能 (要再設計)",
            N_min=N_min, R_min=R_min,
        )
    try:
        reb_utility = select_heating_utility(T_bot, mode='continuous')
    except ValueError:
        return _penalty_result(
            design,
            f"reboiler T_bot={T_bot-273.15:.1f}°C は熱媒の上限超過 (要再設計)",
            N_min=N_min, R_min=R_min,
        )

    # 伝熱面積: contest §4-4 (lookup_U) で (hot_phase, cold_phase) ペアから U 索引。
    # condenser: 蒸気凝縮 (CONDENSING、constant T_top_for_util) ⇄ 冷媒 (utility_phase)
    cond_T_in  = cond_utility.supply_T_K
    cond_T_out = cond_utility.supply_T_K + 10.0
    dT1 = T_top_for_util - cond_T_in
    dT2 = T_top_for_util - cond_T_out
    if dT1 <= 0.0 or dT2 <= 0.0:
        return _penalty_result(
            design,
            f"condenser ΔT 不成立 (T_top={T_top-273.15:.1f}°C, "
            f"supply={cond_utility.supply_T_K-273.15:.1f}°C)",
            N_min=N_min, R_min=R_min,
        )
    if abs(dT1 - dT2) < 1e-3:
        dT_lm_cond = 0.5 * (dT1 + dT2)
    else:
        dT_lm_cond = (dT1 - dT2) / math.log(dT1 / dT2)
    U_cond  = lookup_U(StreamPhase.CONDENSING, utility_phase(cond_utility.name))
    A_cond  = max(Q_cond_kW * 1000.0 / (U_cond * dT_lm_cond), 10.0)

    # reboiler: 熱媒 (utility_phase、通常 CONDENSING の蒸気か GAS の燃料炉) ⇄ 液沸騰 (EVAPORATING)
    dT_lm_reb = reb_utility.supply_T_K - T_bot
    if dT_lm_reb <= 0.0:
        return _penalty_result(
            design,
            f"reboiler ΔT 不成立 (T_bot={T_bot-273.15:.1f}°C, "
            f"supply={reb_utility.supply_T_K-273.15:.1f}°C)",
            N_min=N_min, R_min=R_min,
        )
    U_reb  = lookup_U(utility_phase(reb_utility.name), StreamPhase.EVAPORATING)
    A_reb  = max(Q_reb_kW * 1000.0 / (U_reb * dT_lm_reb), 10.0)

    # ---- Step 10b: CAPEX (Vessel + Trays + Condenser + Reboiler) ----
    capex_vessel = _vessel_capex_okuyen(V_col, design.P_col, D_col)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        capex_trays = calc_tray_capex_okuyen(D_col, design.N_stages)
        capex_cond  = calc_he_capex_okuyen(A_cond)
        capex_reb   = calc_he_capex_okuyen(A_reb)
    capex_total = capex_vessel + capex_trays + capex_cond + capex_reb

    # ---- Step 11: 結果 ----
    top_stream = ProcessStream(
        F_in=F_top, T_in=T_top, P_in=design.P_col,
    )
    bottom_stream = ProcessStream(
        F_in=F_bot, T_in=T_bot, P_in=design.P_col,
    )
    equipment = DistEquipment(
        D_col=D_col, H_col=H_col, V_col=V_col,
        CAPEX_vessel=capex_vessel,
        CAPEX_trays=capex_trays,
        CAPEX=capex_total,
        Q_cond=Q_cond_kW, Q_reb=Q_reb_kW,
        Q_feed_preheat_kW=Q_feed_preheat_kW,
        N_min=N_min, R_min=R_min,
        N_feed_kirkbride=N_feed_kirkbride,
        feasible=True, message="",
        A_cond_m2=A_cond, A_reb_m2=A_reb,
        CAPEX_cond=capex_cond, CAPEX_reb=capex_reb,
        cond_utility_name=cond_utility.name,
        cond_utility_jpy_per_GJ=cond_utility.jpy_per_GJ,
        reb_utility_name=reb_utility.name,
        reb_utility_jpy_per_GJ=reb_utility.jpy_per_GJ,
        T_top=T_top, T_bot=T_bot,
    )
    return DistResult(top=top_stream, bottom=bottom_stream, equipment=equipment)


def _penalty_result(design: DistDesignVars, msg: str,
                    N_min: float = 0.0, R_min: float = 0.0) -> DistResult:
    """infeasibility または計算失敗時のペナルティ結果。"""
    warnings.warn(
        f"simulate_distillation_column: infeasible — {msg}",
        UserWarning, stacklevel=2,
    )
    F_zero = {c: 0.0 for c in ('A', 'B', 'C', 'D', 'E', 'F', 'Z')}
    zero_top = ProcessStream(F_in=dict(F_zero), T_in=298.15, P_in=design.P_col)
    zero_bot = ProcessStream(F_in=dict(F_zero), T_in=298.15, P_in=design.P_col)
    eq = DistEquipment(
        D_col=0.3, H_col=5.0, V_col=0.4,
        CAPEX_vessel=_PENALTY,
        CAPEX_trays=0.0,
        CAPEX=_PENALTY,
        Q_cond=0.0, Q_reb=0.0, Q_feed_preheat_kW=0.0,
        N_min=N_min, R_min=R_min,
        N_feed_kirkbride=0,
        feasible=False, message=msg,
    )
    return DistResult(top=zero_top, bottom=zero_bot, equipment=eq)
