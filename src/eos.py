"""
Peng-Robinson 状態方程式モジュール

C3H6 (キー 'B') / C3H8 (キー 'A') の膜分離システム設計向け。
config.py の ThermoParams に追加した Tc, Pc, omega を使用する。

公開関数
--------
z_factor(T, P, x, keys, phase)               → Z 因子 [-]
fugacity_coeff(i, T, P, x, keys, Z)          → フガシティー係数 φ_i [-]
residual_enthalpy(T, P, x, keys, Z)          → H^r [J/mol]
residual_entropy(T, P, x, keys, Z)           → S^r [J/(mol·K)]
bubble_point_T(P, x, keys)                   → 泡点温度 T_bp [K]
dew_point_T(P, y, keys)                      → 露点温度 T_dp [K]
compress_isentropic(T1, P1, P2, x, keys, eta)→ (T2_actual [K], W_actual [J/mol])

単位系
------
温度: K,  圧力: Pa,  エンタルピー/エントロピー: J/mol または J/(mol·K)
モル流量を受け取る側で [kmol/h] → [mol/s] 変換を行うこと。
"""

import math
import warnings
from typing import List, Tuple

import numpy as np
from scipy.optimize import brentq

from .config import THERMO_DATA, R

_SQRT2: float = math.sqrt(2.0)
# Ref: Peng D.-Y. & Robinson D.B. (1976) "A New Two-Constant Equation of State",
#      Ind. Eng. Chem. Fundam. 15(1), 59-64. Eq. (5)-(6) で Ω_a/Ω_b が臨界点条件
#      (∂P/∂V)_Tc = (∂²P/∂V²)_Tc = 0 の閉形式解として導出される標準値。
_OA:    float = 0.45724   # PR Ω_a
_OB:    float = 0.07780   # PR Ω_b


# ---------------------------------------------------------------------------
# 内部ヘルパー: 単成分 PR パラメータ
# ---------------------------------------------------------------------------

def _kappa_pr(omega: float) -> float:
    """ソアベ補正係数 κ = 0.37464 + 1.54226ω − 0.26992ω²"""
    return 0.37464 + 1.54226 * omega - 0.26992 * omega**2


def _wilson_K(key: str, T: float, P: float) -> float:
    """Wilson 相関による初期 K 値推算（泡点/露点の初期値用）"""
    p = THERMO_DATA[key]
    return (p.Pc / P) * math.exp(5.373 * (1.0 + _kappa_pr(p.omega)) * (1.0 - p.Tc / T))


def _pr_single(key: str, T: float) -> Tuple[float, float, float, float]:
    """
    単成分 PR パラメータを返す。

    Returns
    -------
    a_i         : [Pa·m⁶/mol²]
    b_i         : [m³/mol]
    d_sqrt_a_dT : d(√aᵢ)/dT  [Pa^0.5·m³/(mol·K)]
    sqrt_a_i    : √aᵢ        [Pa^0.5·m³/mol]
    """
    p = THERMO_DATA[key]
    if math.isnan(p.Tc):
        raise ValueError(f"成分 '{key}' に PR パラメータ (Tc, Pc, omega) がありません。")

    ac   = _OA * R**2 * p.Tc**2 / p.Pc
    b    = _OB * R * p.Tc / p.Pc
    kap  = _kappa_pr(p.omega)
    sqrt_Tr = math.sqrt(T / p.Tc)

    sqrt_a = math.sqrt(ac) * (1.0 + kap * (1.0 - sqrt_Tr))
    # d(√a)/dT = −κ·√ac / (2·√(T·Tc))
    d_sqrt_a_dT = -kap * math.sqrt(ac) / (2.0 * math.sqrt(T * p.Tc))

    return sqrt_a**2, b, d_sqrt_a_dT, sqrt_a


# ---------------------------------------------------------------------------
# 内部ヘルパー: 混合則パラメータ
# ---------------------------------------------------------------------------

def _mix(T: float, P: float, x: List[float], keys: List[str]):
    """
    van der Waals 混合則 (kij=0) による PR パラメータを計算する。

    a_m  = (Σ xᵢ √aᵢ)²
    b_m  = Σ xᵢ bᵢ
    da_m/dT = 2·S·Ṡ  ただし S = Σ xᵢ √aᵢ,  Ṡ = Σ xᵢ d(√aᵢ)/dT

    Returns
    -------
    A      : a_m P / (R²T²)
    B      : b_m P / (RT)
    a_m, b_m, da_m_dT
    sq_list  : [√a₀, √a₁, ...]
    b_list   : [b₀, b₁, ...]
    S        : Σ xᵢ √aᵢ
    """
    n = len(keys)
    sq_list, b_list, dsq_list = [], [], []

    for key in keys:
        _, b_i, dsq_i, sq_i = _pr_single(key, T)
        sq_list.append(sq_i)
        b_list.append(b_i)
        dsq_list.append(dsq_i)

    S    = sum(x[i] * sq_list[i]  for i in range(n))
    dS   = sum(x[i] * dsq_list[i] for i in range(n))
    a_m  = S**2
    b_m  = sum(x[i] * b_list[i]   for i in range(n))
    da_m_dT = 2.0 * S * dS

    A = a_m * P / (R**2 * T**2)
    B = b_m * P / (R * T)

    return A, B, a_m, b_m, da_m_dT, sq_list, b_list, S


# ---------------------------------------------------------------------------
# 内部ヘルパー: 3 次方程式の実根
# ---------------------------------------------------------------------------

def _cubic_z(A: float, B: float) -> List[float]:
    """
    Z³ − (1−B)Z² + (A−3B²−2B)Z − (AB−B²−B³) = 0 の実根を返す。
    物理的下限 Z > B のみ採用。
    """
    coeffs = [
        1.0,
        -(1.0 - B),
        A - 3.0*B**2 - 2.0*B,
        -(A*B - B**2 - B**3),
    ]
    roots = np.roots(coeffs)
    return sorted([
        r.real for r in roots
        if abs(r.imag) < 1e-8 * max(abs(r.real), 1.0) and r.real > B + 1e-10
    ])


# ---------------------------------------------------------------------------
# 内部ヘルパー: 理想気体エンタルピー・エントロピー差分
# ---------------------------------------------------------------------------

def _dh_ig(T1: float, T2: float, x: List[float], keys: List[str]) -> float:
    """H^ig(T2) − H^ig(T1) [J/mol]"""
    dH = 0.0
    for i, key in enumerate(keys):
        p = THERMO_DATA[key]
        dH += x[i] * (
            p.a         * (T2    - T1)
            + p.b / 2.0 * (T2**2 - T1**2)
            + p.c / 3.0 * (T2**3 - T1**3)
            + p.d / 4.0 * (T2**4 - T1**4)
        )
    return dH


def _ds_ig(T1: float, T2: float, P1: float, P2: float,
           x: List[float], keys: List[str]) -> float:
    """S^ig(T2,P2) − S^ig(T1,P1) [J/(mol·K)]"""
    dS = 0.0
    for i, key in enumerate(keys):
        p = THERMO_DATA[key]
        dS += x[i] * (
            p.a         * math.log(T2 / T1)
            + p.b       * (T2    - T1)
            + p.c / 2.0 * (T2**2 - T1**2)
            + p.d / 3.0 * (T2**3 - T1**3)
        )
    dS -= R * math.log(P2 / P1)
    return dS


# ---------------------------------------------------------------------------
# 公開関数: Z 因子
# ---------------------------------------------------------------------------

def z_factor(T: float, P: float, x: List[float], keys: List[str],
             phase: str = 'vapor') -> float:
    """
    圧縮率因子 Z を返す。

    Parameters
    ----------
    T, P   : 温度 [K], 圧力 [Pa]
    x      : モル分率リスト（Σ=1）
    keys   : 成分キーリスト（例: ['B','A'] = [C3H6, C3H8]）
    phase  : 'vapor' → 最大根,  'liquid' → 最小正根
    """
    A, B, *_ = _mix(T, P, x, keys)
    roots = _cubic_z(A, B)
    if not roots:
        warnings.warn(f"z_factor: 実根なし (T={T:.1f}K, P={P/1e5:.2f}bar)。Z=1 を返します。")
        return 1.0
    return max(roots) if phase == 'vapor' else min(roots)


# ---------------------------------------------------------------------------
# 公開関数: フガシティー係数
# ---------------------------------------------------------------------------

def fugacity_coeff(comp_idx: int,
                   T: float, P: float, x: List[float], keys: List[str],
                   Z: float) -> float:
    """
    成分 comp_idx のフガシティー係数 φᵢ を返す。

    ln φᵢ = (bᵢ/b_m)(Z−1) − ln(Z−B)
            − A/(2√2 B) · (2√aᵢ/S − bᵢ/b_m) · ln L

    L = [Z+(1+√2)B] / [Z+(1−√2)B]
    """
    A, B, _, b_m, _, sq_list, b_list, S = _mix(T, P, x, keys)

    b_i  = b_list[comp_idx]
    sq_i = sq_list[comp_idx]

    ln_ZmB = math.log(max(Z - B, 1e-30))
    ln_L   = math.log(max((Z + (1.0 + _SQRT2)*B) / (Z + (1.0 - _SQRT2)*B), 1e-30))

    term1 = (b_i / b_m) * (Z - 1.0)
    term2 = -ln_ZmB
    term3 = -(A / (2.0 * _SQRT2 * B)) * (2.0 * sq_i / S - b_i / b_m) * ln_L

    return math.exp(term1 + term2 + term3)


# ---------------------------------------------------------------------------
# 公開関数: 残差熱力学量
# ---------------------------------------------------------------------------

def residual_enthalpy(T: float, P: float, x: List[float], keys: List[str],
                      Z: float) -> float:
    """
    残差エンタルピー H^r [J/mol]

    H^r = RT(Z−1) + (T·da_m/dT − a_m) / (2√2 b_m) · ln L
    """
    _, B, a_m, b_m, da_m_dT, *_ = _mix(T, P, x, keys)
    ln_L = math.log(max((Z + (1.0 + _SQRT2)*B) / (Z + (1.0 - _SQRT2)*B), 1e-30))
    return R*T*(Z - 1.0) + (T*da_m_dT - a_m) / (2.0*_SQRT2*b_m) * ln_L


def residual_entropy(T: float, P: float, x: List[float], keys: List[str],
                     Z: float) -> float:
    """
    残差エントロピー S^r [J/(mol·K)]

    S^r = R·ln(Z−B) + (da_m/dT) / (2√2 b_m) · ln L

    導出: S^r = (H^r − G^r) / T,  G^r = RT·Σxᵢ ln φᵢ
    """
    _, B, _, b_m, da_m_dT, *_ = _mix(T, P, x, keys)
    ln_L   = math.log(max((Z + (1.0 + _SQRT2)*B) / (Z + (1.0 - _SQRT2)*B), 1e-30))
    ln_ZmB = math.log(max(Z - B, 1e-30))
    return R*ln_ZmB + da_m_dT / (2.0*_SQRT2*b_m) * ln_L


# ---------------------------------------------------------------------------
# 公開関数: 泡点温度
# ---------------------------------------------------------------------------

def bubble_point_T(P: float, x: List[float], keys: List[str],
                   T_lo: float = 150.0, T_hi: float = 500.0) -> float:
    """
    泡点温度 T_bp [K] を返す（圧力 P [Pa], 液相組成 x）。

    収束条件: Σ xᵢ Kᵢ = 1,  Kᵢ = φᵢ^L / φᵢ^V
    外側ループ: T を brentq で探索
    内側ループ: thermo の PRMIX でフガシティー係数取得 (= 単相→二相遷移を正しく扱う)

    --- 実装の経緯 (2026-05-10) ---
    旧版は手作り PR EOS (`_mix`, `_cubic_z`, `fugacity_coeff` 等) を使っていたが、
    PR EOS の単相→二相遷移境界 (Z_V = Z_L 縮退点) で brentq が偽根を返す問題があり、
    Dist2 stage 1 で T = -92°C (f = -0.30) を返す不具合が発生した。
    対症療法として post-validate (|f| > 0.1 なら NaN) を入れていたが原則的に悪手。

    本版では thermo (CalebBell/thermo, MIT, 0.6.0 pin) の PRMIX を使う。
    thermo は cubic root 切替を正しく扱い:
      - 両相成立: phis_l と phis_g 両方が valid な値
      - 単相: AttributeError (= 「この T では液相 (or 気相) しか存在しない」signal)
    これで偽根問題が根本解決し、post-validate も不要。

    src/eos.py の他関数 (z_factor, fugacity_coeff 等) は手作り版のまま (= 影響なし、
    呼び出し側は src/eos.py の API を unchanged で利用)。

    Note:
      探索範囲 [150K, 500K] のデフォルトは C3H6/C3H8 等の混合を想定。極低温成分
      (H2 等) のみを含む組成は単相領域に張り付くため適用外。
    """
    # 局所 import: 起動オーバーヘッド削減 + 循環 import 回避
    from thermo.eos_mix import PRMIX as _PRMIX
    from .config import THERMO_DATA as _THERMO_DATA

    n = len(keys)
    Tcs    = [_THERMO_DATA[k].Tc    for k in keys]
    Pcs    = [_THERMO_DATA[k].Pc    for k in keys]
    omegas = [_THERMO_DATA[k].omega for k in keys]

    def obj(T: float) -> float:
        # x 固定で 液相 phi_L を取得 (thermo の PRMIX に zs=x を渡す)
        # thermo は cubic root 切替境界も正しく処理する
        try:
            eos_x = _PRMIX(T=T, P=P, zs=x, Tcs=Tcs, Pcs=Pcs, omegas=omegas)
            phi_L = eos_x.phis_l
        except (AttributeError, ValueError):
            # 単相領域 (液相 root が無い) → bubble より上、f は正の方向に発散扱い
            # brentq に sign を与えて bracketing できるようにする
            return 2.0
        if phi_L is None:
            return 2.0

        # Wilson 相関で y を初期化 (低圧近似だが反復で補正される)
        K = [_wilson_K(keys[i], T, P) for i in range(n)]
        y = [x[i] * K[i] for i in range(n)]
        s = sum(y)
        if s < 1e-30:
            return -2.0
        y = [yi / s for yi in y]

        # 逐次置換: y を更新しながら phi_V を再計算
        for _ in range(50):
            try:
                eos_y = _PRMIX(T=T, P=P, zs=y, Tcs=Tcs, Pcs=Pcs, omegas=omegas)
                phi_V = eos_y.phis_g
            except (AttributeError, ValueError):
                # 単相 (気相 root が無い) → bubble より下、f は負の方向に発散扱い
                return -2.0
            if phi_V is None:
                return -2.0
            K_new = [phi_L[i] / max(phi_V[i], 1e-30) for i in range(n)]
            y_new = [x[i] * K_new[i] for i in range(n)]
            s = sum(y_new)
            if s < 1e-30:
                return -2.0
            y_new = [yi / s for yi in y_new]
            if max(abs(y_new[i] - y[i]) for i in range(n)) < 1e-7:
                K = K_new
                break
            K = K_new
            y = y_new

        return sum(x[i] * K[i] for i in range(n)) - 1.0

    try:
        return brentq(obj, T_lo, T_hi, xtol=0.05, maxiter=200)
    except ValueError:
        warnings.warn(
            f"bubble_point_T: brentq 収束失敗 [{T_lo:.1f}, {T_hi:.1f}] K。nan を返します。",
            UserWarning, stacklevel=2,
        )
        return float('nan')


# ---------------------------------------------------------------------------
# 公開関数: 露点温度
# ---------------------------------------------------------------------------

def dew_point_T(P: float, y: List[float], keys: List[str],
                T_lo: float = 150.0, T_hi: float = 500.0) -> float:
    """
    露点温度 T_dp [K] を返す（圧力 P [Pa], 気相組成 y）。

    収束条件: Σ yᵢ / Kᵢ = 1
    外側ループ: T を brentq で探索
    内側ループ: y 固定 → φ^V 固定、逐次置換で x・φ^L を収束

    Note: デフォルト探索範囲 [150K, 500K] は C3H6/C3H8 混合を対象として設定。
          H2・CH4 など沸点が極低温の成分を含む混合物には適用不可。
    """
    n = len(keys)

    def obj(T: float) -> float:
        A, B, *_ = _mix(T, P, y, keys)
        roots = _cubic_z(A, B)

        # 実根が 1 本のとき: 単相領域
        # 露点条件 = Σ y_i/K_i − 1: T < T_dp (液相) → obj > 0, T > T_dp (気相) → obj < 0
        if len(roots) < 2:
            Z_s = roots[0] if roots else 1.0
            return +2.0 if Z_s < 0.5 else -2.0  # 泡点とは符号が逆

        Z_V = max(roots)
        phi_V = [fugacity_coeff(i, T, P, y, keys, Z_V) for i in range(n)]

        # Wilson 相関で x を初期化
        K = [_wilson_K(keys[i], T, P) for i in range(n)]
        x = [y[i] / max(K[i], 1e-30) for i in range(n)]
        s = sum(x)
        x = [xi / max(s, 1e-30) for xi in x]

        # 逐次置換: φ^L を x に合わせて更新
        for _ in range(50):
            Z_L   = z_factor(T, P, x, keys, 'liquid')
            phi_L = [fugacity_coeff(i, T, P, x, keys, Z_L) for i in range(n)]
            K_new = [phi_L[i] / max(phi_V[i], 1e-30) for i in range(n)]
            x_new = [y[i] / max(K_new[i], 1e-30) for i in range(n)]
            s     = sum(x_new)
            x_new = [xi / max(s, 1e-30) for xi in x_new]
            if max(abs(x_new[i] - x[i]) for i in range(n)) < 1e-7:
                K = K_new; x = x_new; break
            K = K_new; x = x_new

        return sum(y[i] / max(K[i], 1e-30) for i in range(n)) - 1.0

    try:
        return brentq(obj, T_lo, T_hi, xtol=0.05, maxiter=200)
    except ValueError:
        warnings.warn(
            f"dew_point_T: brentq 収束失敗 [{T_lo:.1f}, {T_hi:.1f}] K。nan を返します。",
            UserWarning, stacklevel=2,
        )
        return float('nan')


# ---------------------------------------------------------------------------
# 公開関数: 断熱圧縮
# ---------------------------------------------------------------------------

def compress_isentropic(
    T1: float, P1: float, P2: float,
    x: List[float], keys: List[str],
    eta: float = 0.80,
) -> Tuple[float, float]:
    """
    断熱圧縮シミュレーション（Peng-Robinson 補正付き）。

    手順
    ----
    1. 等エントロピー条件 ΔS_ig + ΔS^r = 0 から T2s を brentq で逆算
    2. W_isen = ΔH_ig(T1→T2s) + ΔH^r(T2s,P2) − ΔH^r(T1,P1)
    3. W_actual = W_isen / η
    4. エンタルピー収支 ΔH_ig + ΔH^r = W_actual から T2_actual を逆算

    Parameters
    ----------
    T1, P1  : 入口温度 [K], 圧力 [Pa]
    P2      : 出口圧力 [Pa]
    x       : モル分率リスト
    keys    : 成分キーリスト
    eta     : 断熱効率 [-] (デフォルト 0.80)

    Returns
    -------
    T2_actual : 実出口温度 [K]
    W_actual  : 実圧縮仕事 [J/mol]（正値 = 圧縮機への入力）
    """
    Z1   = z_factor(T1, P1, x, keys, 'vapor')
    Sr1  = residual_entropy(T1, P1, x, keys, Z1)
    Hr1  = residual_enthalpy(T1, P1, x, keys, Z1)

    # ---- 等エントロピー出口温度 T2s ----
    def entropy_balance(T2: float) -> float:
        Z2  = z_factor(T2, P2, x, keys, 'vapor')
        Sr2 = residual_entropy(T2, P2, x, keys, Z2)
        return _ds_ig(T1, T2, P1, P2, x, keys) + (Sr2 - Sr1)

    # T1 における混合 Cp から比熱比 κ = Cp / (Cp - R) を計算（理想気体近似）
    Cp_mix = sum(
        x[i] * (THERMO_DATA[keys[i]].a + THERMO_DATA[keys[i]].b * T1
                + THERMO_DATA[keys[i]].c * T1**2 + THERMO_DATA[keys[i]].d * T1**3)
        for i in range(len(keys))
    )
    kappa_approx = Cp_mix / (Cp_mix - R)
    T2s_ig = T1 * (P2 / P1) ** ((kappa_approx - 1.0) / kappa_approx)
    T_lo = max(T1 + 0.5, T2s_ig * 0.5)
    T_hi = min(T2s_ig * 2.5, 1200.0)

    try:
        T2s = brentq(entropy_balance, T_lo, T_hi, xtol=0.1, maxiter=200)
    except ValueError:
        try:
            T2s = brentq(entropy_balance, T1 + 0.5, 1200.0, xtol=0.1, maxiter=300)
        except ValueError:
            warnings.warn(
                "compress_isentropic: entropy_balance brentq 収束失敗。理想気体近似で T2s を推算。",
                UserWarning, stacklevel=2,
            )
            T2s = T1 * (P2 / P1) ** ((kappa_approx - 1.0) / kappa_approx)

    Z2s  = z_factor(T2s, P2, x, keys, 'vapor')
    Hr2s = residual_enthalpy(T2s, P2, x, keys, Z2s)
    W_isen  = _dh_ig(T1, T2s, x, keys) + (Hr2s - Hr1)
    W_actual = W_isen / eta

    # ---- 実出口温度 T2_actual ----
    def enthalpy_balance(T2: float) -> float:
        Z2  = z_factor(T2, P2, x, keys, 'vapor')
        Hr2 = residual_enthalpy(T2, P2, x, keys, Z2)
        return _dh_ig(T1, T2, x, keys) + (Hr2 - Hr1) - W_actual

    T2_lo = T1
    T2_hi = max(T2s * 2.0, T1 + 300.0)
    T2_hi = min(T2_hi, 1200.0)

    try:
        T2_actual = brentq(enthalpy_balance, T2_lo, T2_hi, xtol=0.1, maxiter=200)
    except ValueError:
        # フォールバック: 理想気体近似
        Cp_avg = sum(
            x[i] * (THERMO_DATA[keys[i]].a + THERMO_DATA[keys[i]].b * T1)
            for i in range(len(keys))
        )
        T2_actual = T1 + W_actual / max(Cp_avg, 1.0)

    return T2_actual, W_actual
