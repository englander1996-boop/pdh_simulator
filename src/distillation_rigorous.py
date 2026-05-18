"""
蒸留塔 rigorous VLE 求解 (Wang-Henke 法、MESH 方程式 tray-by-tray)。

設計判断 (2026-05-09): FUG (shortcut) との切替 BO 用 rigorous モード本体。
  ColumnTunables.solver_method='rigorous' のとき distillation_core.py の
  _simulate_rigorous() からこのモジュールが呼ばれる。

アルゴリズム: Wang-Henke (bubble-point method)
  Ref: Seader, Henley & Roper, "Separation Process Principles", 3rd ed., Ch.10.4
       (the standard tray-by-tray multi-component distillation method)

主要仮定:
  - CMO (Constant Molar Overflow): 各セクション (feed の上下) で V, L 一定
  - K 値: PR EOS (src/eos.py の bubble_point_T を流用、thermo と 0.02% 整合確認済)
  - 数値ガード: stage 毎の T は _T_BUBBLE_MIN/MAX (200/600 K) でクランプ
  - 収束: damping factor で T プロファイルを更新、max 200 iter で打ち切り

partial condenser 対応:
  - Stage 1 を「凝縮器 = 平衡段」として扱い、x_1 と y_1 で VLE が成立
  - Vapor distillate D_V = V_1 (= D in 設計)、reflux L = R × D_V
  - 通常段と同じ tridiagonal 行列で扱える (T_1 は bubble-point of x_1)

total condenser 対応:
  - Stage 1 は単純凝縮 (平衡なし): x_{i,1} = y_{i,2} = K_{i,2} x_{i,2}
  - tridiagonal 行列の最上段で B_1 = -1, C_1 = K_{i,2} を入れて表現

narrow-α (Dist3, α≈1.07) 対応:
  - damping factor (デフォルト 0.5) で T 振動を抑制
  - 200 iter で収束しなければ NotImplementedError 様の例外を投げて
    distillation_core.py 側で FUG にフォールバック (warning つき)

外部依存:
  - thermo (CalebBell/thermo, MIT License) は将来 PT flash 用にインストール済みだが、
    本モジュールでは src/eos.py の K 値を流用するため直接利用しない。
    (両者が PR EOS で 0.02% 一致を確認済、依存を最小化するため)
  - scipy.linalg.solve_banded で tridiagonal 系を高速求解
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import solve_banded

from src.eos import bubble_point_T, fugacity_coeff, z_factor


# =============================================================================
# 定数
# =============================================================================

# 数値ガード
_T_MIN_K = 200.0      # stage T の下限
_T_MAX_K = 600.0      # stage T の上限
_EPS     = 1e-12

# 収束パラメータ (デフォルト)
_DEFAULT_MAX_ITER  = 200
_DEFAULT_T_TOL_K   = 0.05    # ストレージ毎の T 変化が < 0.05 K で収束
_DEFAULT_DAMPING   = 0.5     # T プロファイル更新の重み (0=無更新, 1=全更新、Wegstein 初回フォールバック)

# 物理定数 (PR EOS と整合)
_R_GAS = 8.314    # [J/(mol·K)]

# Clausius-Clapeyron 用の蒸発潜熱 [kJ/mol]
# Ref: distillation_core.py の _LAMBDA_KJ と同値、モジュール独立性のため再定義
# Newton 1 ステップ bubble-point の df/dT 計算用 (∂lnK/∂T ≈ ΔH/RT²)
_LAMBDA_KJ_PER_MOL: Dict[str, float] = {
    'A': 19.1,    # C3H8
    'B': 18.4,    # C3H6
    'C':  0.9,    # H2
    'D': 13.5,    # C2H4
    'E':  8.2,    # CH4
    'F': 14.7,    # C2H6
    'Z': 22.0,    # C4H10
}
_LAMBDA_DEFAULT = 15.0    # [kJ/mol] フォールバック

# Newton 1 ステップ bubble-point の安全装置
_NEWTON_DT_MAX_K   = 20.0   # 1 ステップで動かす T の最大幅 [K]
_NEWTON_DT_MIN_TOL = 0.01   # T 変化がこれ未満なら収束扱い

# Wegstein 加速のクランプ範囲 (q ∈ [q_min, q_max])
# q < 0: 加速、q = 0: pure substitution、0 < q < 1: damping
# 教科書標準: 安定のため q ∈ [-5, 0] (加速のみ、damping は除く)
# Ref: Seader, Henley & Roper (2010) §6.2.3
_WEGSTEIN_Q_MIN = -5.0
_WEGSTEIN_Q_MAX = 0.0


# =============================================================================
# データクラス
# =============================================================================

@dataclass
class RigorousResult:
    """Wang-Henke 求解結果。distillation_core.py 側で DistResult に詰め直される。

    always-on validation フィールド (収束時に毎回計算):
      mesh_residual_max    : 中間段の MESH 方程式残差最大値 (正規化済 / F_total)
      mesh_residual_mean   : 同 平均
      component_balance_max: 成分マスバランス最大相対誤差 (= |feed - top - bot| / feed)
    値が大きい場合は内部解が不整合の可能性 (B_bottoms バグのような問題の早期検知)。
    """
    converged:    bool                # 収束したか
    n_iter:       int                 # 反復回数
    message:      str                 # 説明 (失敗理由など)
    T_profile_K:  List[float]         # 各段の温度 [K] (stage 1 = 塔頂)
    x_profile:    List[Dict[str, float]]  # 各段の液相組成 (mol fraction)
    y_profile:    List[Dict[str, float]]  # 各段の気相組成 (mol fraction)
    K_profile:    List[Dict[str, float]]  # 各段の K 値
    F_top:        Dict[str, float]    # 塔頂 (distillate) 流量 [kmol/h]
    F_bot:        Dict[str, float]    # 塔底 (bottoms) 流量 [kmol/h]
    V_top_kmolh:  float               # 塔頂 vapor 流量 [kmol/h]
    V_bot_kmolh:  float               # 塔底 (reboiler) vapor 流量 [kmol/h]
    L_top_kmolh:  float               # 塔頂 reflux liquid 流量 [kmol/h]
    L_bot_kmolh:  float               # 塔底 (リボイラ入口) liquid 流量 [kmol/h]
    # always-on validation (収束時のみ意味あり)
    mesh_residual_max:      float = 0.0
    mesh_residual_mean:     float = 0.0
    component_balance_max:  float = 0.0


# =============================================================================
# Wang-Henke 本体
# =============================================================================

def wang_henke_solve(
    feed_F:           Dict[str, float],     # フィード流量 [kmol/h]
    comps:            List[str],
    P_col:            float,                # 塔操作圧力 [Pa]
    N_stages:         int,                  # 全段数 (理論段数)
    N_feed:           int,                  # フィード段位置 (1=塔頂寄り)
    reflux_ratio:     float,                # R = L_reflux / D
    D_total:          float,                # 塔頂産物量 [kmol/h] (FUG 推算値を使用)
    q_feed:           float,                # フィード状態 (1=飽和液, 0=飽和気)
    partial_condenser:bool,                 # 分流型 (Dist2 用)
    K_method:         str   = 'pr',
    max_iter:         int   = _DEFAULT_MAX_ITER,
    T_tol_K:          float = _DEFAULT_T_TOL_K,
    damping:          float = _DEFAULT_DAMPING,
    T_top_init_K:     Optional[float] = None,
    T_bot_init_K:     Optional[float] = None,
) -> RigorousResult:
    """Wang-Henke で多成分蒸留塔の MESH を解く。

    入力は通常 FUG ソルバから受け取る:
      D_total (塔頂産物流量) は FUG の Fenske split 結果を使う
      T_top_init_K, T_bot_init_K も FUG の bubble-point 結果を初期値に

    出力 RigorousResult を distillation_core.py 側で DistResult に詰め直す。

    数値挙動:
      - narrow-α 系 (Dist3) は damping を低めに設定すると振動抑制
      - 単相領域に逸脱した場合は T を _T_MIN/MAX_K にクランプして次反復へ
      - max_iter で打ち切り、converged=False で返す (呼び出し側で FUG フォールバック)
    """
    # ---- 1. 設計判断 (基本パラメータ) ----
    F_total = sum(max(F, 0.0) for F in feed_F.values())
    if F_total <= _EPS:
        return _failure_result("feed flow ≤ 0", N_stages, comps)

    z_feed = {c: max(feed_F.get(c, 0.0), 0.0) / F_total for c in comps}

    # 段番号: 1 = 塔頂 (condenser, partial=平衡段 / total=単純凝縮),
    #        N = 塔底 (reboiler, 平衡段)
    if N_stages < 2:
        return _failure_result("N_stages < 2", N_stages, comps)
    if not (1 <= N_feed <= N_stages):
        return _failure_result(f"N_feed={N_feed} out of [1, {N_stages}]", N_stages, comps)

    # ---- 2. CMO による V, L プロファイル ----
    # 塔頂段の上には流出 D_total。L_top = R × D_total が下に向かう。
    # 塔底段からは B = F - D_total が出る。
    D = max(D_total, _EPS)
    B = F_total - D
    if B <= _EPS:
        return _failure_result(f"B = F - D = {B:.3f} ≤ 0 (D too large?)", N_stages, comps)

    L_top = reflux_ratio * D            # 塔頂部の液流量 [kmol/h]
    V_top = (reflux_ratio + 1.0) * D    # 塔頂部の気流量 [kmol/h]
    L_bot = L_top + q_feed * F_total    # 塔底部の液流量 (feed 後)
    V_bot = V_top - (1.0 - q_feed) * F_total  # 塔底部の気流量

    if V_bot <= _EPS:
        return _failure_result(
            f"V_bot = {V_bot:.3f} ≤ 0 (R too small or q too low for given D/F)",
            N_stages, comps,
        )

    # 段毎の V, L 配列 (1-indexed: index 0 は使わない)
    # 設計判断 (2026-05-09): Wang-Henke の正しい CMO 規約:
    #   L[j] = stage j → stage j+1 (下方向の液流量)
    #   V[j] = stage j → stage j-1 (上方向の気流量)
    # フィード段 j=N_feed では:
    #   liquid IN from above = L[N_feed-1] (= L_top)
    #   liquid IN from feed = q × F → 加算
    #   liquid OUT down      = L[N_feed]   (= L_bot = L_top + qF) ✓
    # 一方 vapor は:
    #   vapor IN from below  = V[N_feed+1] (= V_bot)
    #   vapor IN from feed   = (1-q) × F → 加算
    #   vapor OUT up         = V[N_feed]   (= V_top = V_bot + (1-q)F) ★
    # つまり V[N_feed] は **V_top** であって V_bot ではない。
    L = np.zeros(N_stages + 2)
    V = np.zeros(N_stages + 2)
    for j in range(1, N_stages + 1):
        # L[j] (下方向): N_feed 以上で L_bot、それ未満で L_top
        if j < N_feed:
            L[j] = L_top
        else:
            L[j] = L_bot
        # V[j] (上方向): フィード段 N_feed では V_top (上向きに出る方)、
        #               N_feed より下 (= N_feed+1 以降) で V_bot
        if j <= N_feed:
            V[j] = V_top
        else:
            V[j] = V_bot
    # 塔頂・塔底の境界
    V[N_stages + 1] = 0.0    # 塔底の下にはガスは無い
    L[0]            = 0.0    # 塔頂の上には液は無い (還流は L[1] として下りる)
    # partial_condenser: V[1] = D (vapor distillate)
    # total_condenser:   V[1] = 0 (no vapor leaves, all condenses then split)
    if partial_condenser:
        V[1] = D
        # L[1] = L_top (reflux) はすでにセット済み
    else:
        V[1] = 0.0
        # L[1] = L_top (reflux) はすでにセット済み

    # ---- 3. 初期 T プロファイル (線形) ----
    # 設計判断: 塔頂 (LK 主体) と塔底 (HK 主体) の純成分沸点を初期推定。
    # FUG から渡された値があればそれを使う、無ければ z 全成分で粗推定。
    if T_top_init_K is None or T_bot_init_K is None:
        # フィード組成で泡点 T を粗推定 (どちらの sec でも良い、線形 init の参考だけ)
        try:
            x_init = [z_feed[c] for c in comps]
            T_mid = bubble_point_T(P_col, x_init, comps)
        except Exception:
            T_mid = 350.0
        T_top_init_K = T_top_init_K if T_top_init_K is not None else T_mid - 30.0
        T_bot_init_K = T_bot_init_K if T_bot_init_K is not None else T_mid + 30.0
    T_top_init_K = max(_T_MIN_K, min(_T_MAX_K, T_top_init_K))
    T_bot_init_K = max(_T_MIN_K, min(_T_MAX_K, T_bot_init_K))
    if T_bot_init_K <= T_top_init_K:
        T_bot_init_K = T_top_init_K + 5.0

    T = np.linspace(T_top_init_K, T_bot_init_K, N_stages)   # 0-indexed → stage j+1

    # ---- 4. 初期 x プロファイル (フィード組成で全段一様) ----
    n_comp = len(comps)
    x = np.zeros((N_stages, n_comp))
    for j in range(N_stages):
        for i, c in enumerate(comps):
            x[j, i] = z_feed[c]

    # ---- 5. 反復ループ (Newton 1ステップ bubble-point + Wegstein 加速) ----
    # Wegstein 用の前回履歴 (各段の T と Newton-計算 T_calc を覚えておく)
    T_prev      = None     # T (= 前回の更新後)
    T_calc_prev = None     # T_calc (= 前回 Newton step が出した raw 値)
    n_iter_done = 0
    dT_max      = float('inf')

    for it in range(max_iter):
        n_iter_done = it + 1

        # 5a. 現在の T, x で K 値を計算 (各段で PR EOS + phi_L/phi_V)
        # Newton step では K 値を内部で再計算するので、ここの K_profile は
        # tridiagonal 用 (= MESH の係数行列) のみに使う。
        K = _compute_K_profile(T, x, comps, P_col, K_method)
        if K is None:
            return _failure_result(
                f"iter {it}: K 値計算失敗 (単相領域逸脱)",
                N_stages, comps,
            )

        # 5b. 各成分について tridiagonal 系を解く
        x_new = np.zeros_like(x)
        for i in range(n_comp):
            try:
                x_new[:, i] = _solve_tridiagonal(
                    i, K, V, L, F_total, z_feed, comps,
                    N_stages, N_feed, partial_condenser,
                    B_bottoms=B,
                )
            except np.linalg.LinAlgError as e:
                return _failure_result(
                    f"iter {it}: tridiagonal LinAlg error ({e})",
                    N_stages, comps,
                )

        # 5c. x の正規化 (各段で sum_i x = 1 を強制)
        for j in range(N_stages):
            s = x_new[j].sum()
            if s > _EPS:
                x_new[j] /= s
            else:
                for i, c in enumerate(comps):
                    x_new[j, i] = z_feed[c]

        # 5d. 各段で Newton 1 ステップ → T_calc プロファイル
        # 旧版は scipy.brentq で完全収束させていたが、Wang-Henke 外側で
        # 何度も呼ばれるので各段は粗い更新で十分。1 ステップ Newton で十分速い。
        # partial_condenser stage 1 では bubble-point が単相領域に落ちて
        # 真値から ±20K ズレることがある (Dist2 既知問題、KNOWN_PLACEHOLDERS.md
        # 参照)。ただし T を凍結すると x との整合が崩れて MESH 残差が爆発する
        # ため、T の更新は通常通り行う。outlier は MESH の自然な収束点。
        T_calc = np.copy(T)
        for j in range(N_stages):
            T_calc[j], _ = _bubble_T_newton_step(
                T_old=T[j], x=list(x_new[j]),
                P=P_col, comps=comps,
            )

        # 5e. T プロファイル更新: Wegstein 加速 (固定 damping より速い)
        # Wegstein: T_new[j] = q[j] T[j] + (1-q[j]) T_calc[j]
        #   q = s/(s-1), s = (T_calc[j] - T_calc_prev[j]) / (T[j] - T_prev[j])
        #   q ∈ [-5, 0] にクランプ (加速のみ、振動方向への発散を防止)
        # 初回 2 反復は履歴がないので固定 damping にフォールバック。
        if T_prev is None or T_calc_prev is None:
            T_damped = (1.0 - damping) * T + damping * T_calc
        else:
            T_damped = np.empty_like(T)
            for j in range(N_stages):
                dT_old   = T[j]      - T_prev[j]
                dT_calc  = T_calc[j] - T_calc_prev[j]
                # 分母ゼロ・極小回避
                if abs(dT_old) < _EPS:
                    q = -damping  # 履歴更新が無い段はフォールバック
                else:
                    s = dT_calc / dT_old
                    if abs(s - 1.0) < _EPS:
                        q = 0.0
                    else:
                        q = s / (s - 1.0)
                # クランプ: 加速のみ許容、damping > 0 と発散側への極端 q は避ける
                q = max(_WEGSTEIN_Q_MIN, min(_WEGSTEIN_Q_MAX, q))
                T_damped[j] = q * T[j] + (1.0 - q) * T_calc[j]
                # 物理範囲内へクランプ
                T_damped[j] = max(_T_MIN_K, min(_T_MAX_K, T_damped[j]))

        # 履歴更新
        T_prev      = T.copy()
        T_calc_prev = T_calc.copy()

        # 5f. 収束判定 (T プロファイルの最大変化)
        dT_max = float(np.max(np.abs(T_damped - T)))
        T = T_damped
        x = x_new

        if dT_max < T_tol_K:
            break
    else:
        return _failure_result(
            f"収束失敗: max_iter={max_iter}, dT_max={dT_max:.3f}K (tol={T_tol_K}K)",
            N_stages, comps,
        )

    # ---- 5g. 最終 consistency パス (MESH 残差を縮める) ----
    # 設計判断 (2026-05-10): Wegstein は T が静止する点で停止するが、その時点では
    # x が「1反復前の K」で解かれた状態のため、MESH 方程式に微小残差が残る。
    # この残差は成分別マスバランスを破ることがある (Dist1 で C3H8/C4H10 ≈ 5% 入れ替わり等)。
    # 修正: 収束後に「最終 T での K を計算 → x を tridiagonal で再求解」を 2-3 回繰り返し
    # x と K の整合を取る。T は固定 (= 静止点維持)、x のみ磨く。
    # 経験的に 2-3 回で MESH 残差は 1e-3 → 1e-6 オーダーまで縮む。
    for polish_iter in range(3):
        K_polish = _compute_K_profile(T, x, comps, P_col, K_method)
        if K_polish is None:
            break
        x_polish = np.zeros_like(x)
        polish_failed = False
        for i in range(n_comp):
            try:
                x_polish[:, i] = _solve_tridiagonal(
                    i, K_polish, V, L, F_total, z_feed, comps,
                    N_stages, N_feed, partial_condenser,
                    B_bottoms=B,
                )
            except np.linalg.LinAlgError:
                polish_failed = True
                break
        if polish_failed:
            break
        # 正規化
        for j in range(N_stages):
            s = x_polish[j].sum()
            if s > _EPS:
                x_polish[j] /= s
            else:
                x_polish[j] = x[j]
        dx_max = float(np.max(np.abs(x_polish - x)))
        x = x_polish
        if dx_max < 1e-6:
            break

    # ---- 6. 結果整理 ----
    K_final = _compute_K_profile(T, x, comps, P_col, K_method)
    y = np.zeros_like(x)
    for j in range(N_stages):
        for i in range(n_comp):
            y[j, i] = K_final[j, i] * x[j, i]
        s = y[j].sum()
        if s > _EPS:
            y[j] /= s

    # 塔頂産物・塔底産物の組成・流量
    if partial_condenser:
        # 塔頂は y_1 (vapor distillate composition × D)
        F_top_kmolh = {c: D * y[0, i] for i, c in enumerate(comps)}
    else:
        # total cond: 塔頂は x_1 (= y_2、liquid distillate composition × D)
        F_top_kmolh = {c: D * x[0, i] for i, c in enumerate(comps)}
    # 塔底は x_N (bottoms liquid composition × B)
    F_bot_kmolh = {c: B * x[N_stages - 1, i] for i, c in enumerate(comps)}

    # ---- always-on validation (内部不整合の早期検知) ----
    # 設計判断 (2026-05-10): B_bottoms バグの再発防止に、収束時に毎回
    # MESH 残差と成分マスバランスを計算して結果に含める。
    # 大きな値が出れば solver 内部の不整合 (B_bot 取り違え、index off-by-one,
    # K=1 trivial 罠など) を呼び出し側で気付ける。
    mesh_max, mesh_mean, n_res = 0.0, 0.0, 0
    for j_chk in range(2, N_stages):  # 中間段のみ
        x_jm1 = x[j_chk - 2]
        x_j   = x[j_chk - 1]
        x_jp1 = x[j_chk]
        K_j   = K_final[j_chk - 1]
        K_jp1 = K_final[j_chk]
        L_jm1 = L[j_chk - 1]
        L_j_  = L[j_chk]
        V_j_  = V[j_chk]
        V_jp1 = V[j_chk + 1]
        for i, c in enumerate(comps):
            feed_term = feed_F.get(c, 0.0) if j_chk == N_feed else 0.0
            res = (L_jm1 * x_jm1[i] + V_jp1 * K_jp1[i] * x_jp1[i]
                   - (L_j_ + V_j_ * K_j[i]) * x_j[i] + feed_term)
            ar = abs(res) / max(F_total, _EPS)
            if ar > mesh_max:
                mesh_max = ar
            mesh_mean += ar
            n_res += 1
    mesh_mean = mesh_mean / max(n_res, 1)

    # 成分マスバランス検査 (= feed - top - bot per component)
    comp_max_err = 0.0
    for c in comps:
        f_c = feed_F.get(c, 0.0)
        if f_c < _EPS:
            continue
        balance = abs(f_c - F_top_kmolh.get(c, 0.0) - F_bot_kmolh.get(c, 0.0)) / f_c
        if balance > comp_max_err:
            comp_max_err = balance

    return RigorousResult(
        converged    = True,
        n_iter       = it + 1,
        message      = f"converged in {it+1} iter (dT_max={dT_max:.4f}K)",
        T_profile_K  = list(T),
        x_profile    = [{c: x[j, i] for i, c in enumerate(comps)} for j in range(N_stages)],
        y_profile    = [{c: y[j, i] for i, c in enumerate(comps)} for j in range(N_stages)],
        K_profile    = [{c: K_final[j, i] for i, c in enumerate(comps)} for j in range(N_stages)],
        F_top        = F_top_kmolh,
        F_bot        = F_bot_kmolh,
        V_top_kmolh  = V_top,
        V_bot_kmolh  = V_bot,
        L_top_kmolh  = L_top,
        L_bot_kmolh  = L_bot,
        mesh_residual_max     = mesh_max,
        mesh_residual_mean    = mesh_mean,
        component_balance_max = comp_max_err,
    )


# =============================================================================
# 内部ヘルパ
# =============================================================================

def _bubble_T_newton_step(
    T_old: float,
    x:     List[float],
    P:     float,
    comps: List[str],
) -> Tuple[float, List[float]]:
    """各段で 1 回だけ Newton ステップを当てて bubble-point T を更新する。

    bubble-point 条件:
        f(T) = Σ K_i(T) x_i - 1 = 0
    Newton 微分は Clausius-Clapeyron 近似:
        ∂lnK_i/∂T ≈ ΔH_vap_i / (R T²)
        ∂K_i /∂T  ≈ K_i × ΔH_vap_i / (R T²)
    で 1 step 更新:
        T_new = T_old - f(T_old) / (df/dT)

    工夫:
      - 完全収束ではなく 1 ステップだけ。Wang-Henke 外側反復で何度も呼ばれて
        全体としては収束する (= "粗い内側 + Wegstein 加速外側" の戦略)
      - Newton が暴れるケース (極小 df/dT、df/dT 符号反転) に備えて
        T 変化を _NEWTON_DT_MAX_K [K] でクランプする

    Returns
    -------
    (T_new, K_at_T_old)
        K は呼び出し側で再利用できるように T_old での値を返す
    """
    try:
        Z_V = z_factor(T_old, P, x, comps, phase='vapor')
        Z_L = z_factor(T_old, P, x, comps, phase='liquid')
    except Exception:
        return T_old, [1.0] * len(comps)

    # K 値計算
    K_values: List[float] = []
    for i in range(len(comps)):
        try:
            phi_V = fugacity_coeff(i, T_old, P, x, comps, Z_V)
            phi_L = fugacity_coeff(i, T_old, P, x, comps, Z_L)
            K_values.append(phi_L / max(phi_V, _EPS))
        except Exception:
            K_values.append(1.0)

    # f(T) = Σ K x - 1
    f_T = sum(K_values[i] * x[i] for i in range(len(x))) - 1.0

    # df/dT ≈ Σ K × ΔH_vap / (R T²) × x  (Clausius-Clapeyron 微分)
    df_dT = 0.0
    for i, c in enumerate(comps):
        dH_J = _LAMBDA_KJ_PER_MOL.get(c, _LAMBDA_DEFAULT) * 1000.0  # J/mol
        df_dT += K_values[i] * dH_J / (_R_GAS * T_old**2) * x[i]

    if abs(df_dT) < _EPS:
        return T_old, K_values

    dT = -f_T / df_dT
    # 安定化: 1 ステップで大きく動きすぎないようクランプ
    dT = max(-_NEWTON_DT_MAX_K, min(_NEWTON_DT_MAX_K, dT))
    T_new = T_old + dT
    T_new = max(_T_MIN_K, min(_T_MAX_K, T_new))
    return T_new, K_values


def _compute_K_profile(
    T: np.ndarray, x: np.ndarray,
    comps: List[str], P_col: float, K_method: str,
) -> Optional[np.ndarray]:
    """各段で K_{i,j} = phi_L_i / phi_V_i を計算。

    各段の x_{:,j} は液相組成、T_j は仮の段温度。両相が共存する点 (= 泡点近傍)
    では Z_L, Z_V とも有効な根が得られる。仮 T で単相に張り付いた場合は
    None を返して呼び出し側でフォールバック。
    """
    N_stages, n_comp = x.shape
    K = np.zeros_like(x)
    for j in range(N_stages):
        x_j = list(x[j])
        try:
            Z_V = z_factor(T[j], P_col, x_j, comps, phase='vapor')
            Z_L = z_factor(T[j], P_col, x_j, comps, phase='liquid')
        except Exception:
            return None
        for i in range(n_comp):
            try:
                phi_V = fugacity_coeff(i, T[j], P_col, x_j, comps, Z_V)
                phi_L = fugacity_coeff(i, T[j], P_col, x_j, comps, Z_L)
                K[j, i] = phi_L / phi_V if phi_V > _EPS else 0.0
            except Exception:
                return None
    return K


def _solve_tridiagonal(
    i_comp:           int,
    K:                np.ndarray,    # shape (N_stages, n_comp)
    V:                np.ndarray,    # 1-indexed, length N_stages+2
    L:                np.ndarray,    # 1-indexed, length N_stages+2
    F_total:          float,
    z_feed:           Dict[str, float],
    comps:            List[str],
    N_stages:         int,
    N_feed:           int,
    partial_condenser:bool,
    B_bottoms:        float,
) -> np.ndarray:
    """Wang-Henke の段別物質収支を tridiagonal で解く。

    成分 i_comp について、各段の液相組成 x_{i,j} を求める。
    平衡関係 y_{i,j} = K_{i,j} x_{i,j} を代入後の標準形:

      A_j x_{i,j-1} + B_j x_{i,j} + C_j x_{i,j+1} = D_j

    境界条件:
      total condenser (stage 1):  x_{i,1} = K_{i,2} x_{i,2}
        → A_1=0, B_1=-1, C_1=K_{i,2}, D_1=0
      partial condenser (stage 1): mass balance with vapor distillate
        → A_1=0, B_1=-(D_V (R+K_{i,1})), C_1=V_2 K_{i,2}, D_1=0
        (D_V = V[1], 反復中の D は固定)
      reboiler (stage N): L_{N-1} x_{i,N-1} - (B_bottoms + V_N K_{i,N}) x_{i,N} = 0
        → A_N=L_{N-1}, B_N=-(B_bottoms + V_N K_{i,N}), C_N=0, D_N=0
      feed stage:                 D_j = -F z_{i,F}

    Parameters
    ----------
    B_bottoms : float
        塔底産物の流量 [kmol/h] = F_total - D_total。
        L[N_stages] (= L_bot = R·D + q·F) とは別物。
        Bug fix (2026-05-10): 旧版では `B_bot = L[N_stages]` としていたが、
        これは L_bot (= 塔底段に下りる液量) で、bottom product flow B とは
        全く異なる値 (Dist1 で ~20×)。リボイラ方程式が壊れて Dist1 の
        成分マスバランスが破綻していた。
    """
    A = np.zeros(N_stages)   # 0-indexed
    B = np.zeros(N_stages)
    C = np.zeros(N_stages)
    D = np.zeros(N_stages)

    z_i = z_feed[comps[i_comp]]
    B_bot = B_bottoms   # 塔底産物流量 [kmol/h]

    for j in range(1, N_stages + 1):
        idx = j - 1   # 0-indexed for arrays
        K_ij    = K[idx, i_comp] if idx < K.shape[0] else 1.0
        K_ijp1  = K[idx + 1, i_comp] if idx + 1 < K.shape[0] else 1.0

        if j == 1:
            # 塔頂段
            A[idx] = 0.0
            if partial_condenser:
                # 分流型: V[1] = D, L[1] = L_top
                D_V = V[1]
                R   = L[1] / max(D_V, _EPS)
                B[idx] = -(D_V * (R + K_ij))
                C[idx] = V[2] * K_ijp1 if N_stages > 1 else 0.0
                D[idx] = 0.0
            else:
                # total cond: x_1 = K_2 x_2
                B[idx] = -1.0
                C[idx] = K_ijp1
                D[idx] = 0.0
        elif j == N_stages:
            # リボイラ段
            A[idx] = L[j - 1]
            B[idx] = -(B_bot + V[j] * K_ij)
            C[idx] = 0.0
            D[idx] = 0.0
        else:
            # 中間段
            A[idx] = L[j - 1]
            B[idx] = -(L[j] + V[j] * K_ij)
            C[idx] = V[j + 1] * K_ijp1
            D[idx] = 0.0
            if j == N_feed:
                D[idx] = -F_total * z_i

        # フィード段が塔頂/塔底にある場合の補正
        if j == N_feed and (j == 1 or j == N_stages):
            D[idx] = -F_total * z_i

    # tridiagonal 行列を banded 形式に詰める (scipy.linalg.solve_banded)
    # ab[0] = upper diagonal (C), ab[1] = main (B), ab[2] = lower (A)
    ab = np.zeros((3, N_stages))
    ab[0, 1:]   = C[:-1]   # upper diag (C_j → 上方向シフト)
    ab[1, :]    = B
    ab[2, :-1]  = A[1:]    # lower diag (A_j → 下方向シフト)

    x_i = solve_banded((1, 1), ab, D)
    # 物理的に非負を強制
    x_i = np.maximum(x_i, 0.0)
    return x_i


def _failure_result(msg: str, N_stages: int, comps: List[str]) -> RigorousResult:
    """収束失敗時のダミー結果。distillation_core.py 側で警告 + FUG フォールバック対象。

    設計判断 (2026-05-18 確認): RigorousResult dataclass は array フィールド
    (T_profile_K, x/y/K_profile, F_top, F_bot, V/L) を必須としているため、failure 時も
    全フィールドを埋める必要がある。flag-based の signal は `converged=False` で表現済み。
    caller (distillation_core.py:680) は `if not rig.converged: raise RuntimeError(...)`
    で明示的に判定するため、nan/0.0 配列が下流計算に伝播することはない。
    本構造は問題ないため、追加リファクタは不要 (KNOWN_PLACEHOLDERS §4 から削除)。

    T_profile_K に nan、x/y に 0.0、K に 1.0 を埋めるのは「無意味な値」を強調するため。
    """
    return RigorousResult(
        converged    = False,
        n_iter       = 0,
        message      = msg,
        T_profile_K  = [float('nan')] * max(N_stages, 1),
        x_profile    = [{c: 0.0 for c in comps}] * max(N_stages, 1),
        y_profile    = [{c: 0.0 for c in comps}] * max(N_stages, 1),
        K_profile    = [{c: 1.0 for c in comps}] * max(N_stages, 1),
        F_top        = {c: 0.0 for c in comps},
        F_bot        = {c: 0.0 for c in comps},
        V_top_kmolh  = 0.0, V_bot_kmolh = 0.0,
        L_top_kmolh  = 0.0, L_bot_kmolh = 0.0,
    )
