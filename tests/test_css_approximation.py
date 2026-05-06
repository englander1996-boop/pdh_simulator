"""
CSS 近似精度検証スクリプト (R-3)

目的:
    psa_system.py の CSS 乗算近似
        t_abs_css = t_abs_clean × (1 - desorption_target)
    が、実際のサイクルループで得られる CSS 定常破過時間とどの程度一致するかを確認する。

検証手順:
    1. 現行の近似 (use_css_approximation=True) で t_abs_approx を取得
    2. 清浄床 (q=0) からスタートし、「吸着→脱着→吸着→…」を最大 5 サイクル繰り返す
       - 各サイクルの吸着: 前サイクルの残留固相負荷 q_residual を初期値として ODE を解く
       - 各サイクルの脱着: q_residual = q_final × exp(-KFa × t_des) を空間全点に適用
       - t_abs の変化が 1% 未満になったら収束とみなす
    3. 近似値と収束値の誤差 [%] を算出し、10% 以内かどうかで OK/NG を判定する

使い方:
    python tests/test_css_approximation.py

検証結果メモ (2026-05-06, 代表条件 D_col=1.0m, L_bed=5.0m, desorption_target=0.35):
    Langmuir: Choi et al., J. Chem. Eng. Data, 48, 603-607 (2003)
    KFa:      Rufford et al., Ind. Eng. Chem. Res. (2013) ベース
    近似 t_abs = 126.7 s、サイクルループ収束値 (5サイクル) = 97.8 s、誤差 +29.6%。
    旧仮置きパラメータでの誤差 +31.5% とほぼ同じ → 誤差は近似手法に起因し
    パラメータ値によらない。

    t_des ≈ 62.6 s、t_des/t_abs_true ≈ 0.64 → この代表条件では塔数ともに 2 塔。
    ただし N_total_columns = ceil(t_des/t_abs)+1 は CAPEX に直接影響するため、
    近似が t_abs を約 30% 過大推算することで t_abs_approx が t_des をギリギリ上回る
    境界付近では「近似: 2 塔、真値: 3 塔」の逆転が発生し最適解が誤る可能性がある。
    → 塔数境界付近の条件を重点的に確認すること（今後の課題）。

    補足 (設計空間): 新 KFa + 新 Langmuir では D=1.0m・L=5m 以上が有効域。
    D=0.5m や D=1.0m・L=3m は t_abs_css < T_ABS_MIN でペナルティ返却される。

注意:
    - 代表条件は _DESIGN / _FEED で定義。最適化の想定探索域に合わせて適宜変更すること。
    - このスクリプトは本体コード (psa_system.py) を変更せず、内部 ODE を複製して使用する。
"""

import math
import os
import sys
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.config import R
from src.cost_parameters import PSA_LANGMUIR_PARAMS, PSA_KFA
from units.separators.psa.psa_system import (
    PSADesignVars, PSAFeedStream, PSAFixedParams,
    simulate_psa_system,
    _ADS_ORDER, _N_Z, _N_ADS,
    _calc_feed_state,
)


# ---------------------------------------------------------------------------
# 代表条件 — 最適化の想定探索域中央付近。必要に応じて書き換えること。
# ---------------------------------------------------------------------------
_DESIGN = PSADesignVars(D_col=1.0, L_bed=5.0, desorption_target=0.35)
_FEED = PSAFeedStream(
    F_in={'A': 0.0, 'B': 0.0, 'C': 100.0, 'D': 2.0, 'E': 5.0, 'F': 1.0},
    T_in=350.0,
    P_in=20.0e5,
)


# ---------------------------------------------------------------------------
# ODE ヘルパー — _run_adsorption と同一の ODE を q_init 付きで再実装
# ---------------------------------------------------------------------------

def _run_adsorption_with_q_init(
    C_feed, u_0, L_bed, rho_b, eps,
    kfa, q_s, a_lang,
    breakthrough_ratio, t_max,
    q_init,
):
    """
    初期固相負荷 q_init [mol/kg] shape (N_z, N_ADS) から吸着 PDE を積分する。

    Returns
    -------
    t_abs     : 破過時刻 [s]
    q_final   : 固相最終分布 shape (N_z, N_ADS) [mol/kg]
    converged : True = 破過イベントで停止
    """
    dz    = L_bed / _N_Z
    u_eps = u_0 / eps

    def rhs(t, y):
        C = np.maximum(y[:_N_Z * _N_ADS].reshape(_N_Z, _N_ADS), 0.0)
        q = np.maximum(y[_N_Z * _N_ADS:].reshape(_N_Z, _N_ADS), 0.0)
        aC    = a_lang * C
        denom = 1.0 + aC.sum(axis=1, keepdims=True)
        q_eq  = q_s * aC / denom
        dq_dt = kfa * (q_eq - q)
        C_up       = np.empty_like(C)
        C_up[0, :] = C_feed
        C_up[1:, :] = C[:-1, :]
        dC_dt = -u_eps / dz * (C - C_up) - rho_b / eps * dq_dt
        return np.concatenate([dC_dt.ravel(), dq_dt.ravel()])

    def breakthrough_event(t, y):
        return y[(_N_Z - 1) * _N_ADS] - breakthrough_ratio * C_feed[0]
    breakthrough_event.terminal  = True
    breakthrough_event.direction = 1

    y0 = np.zeros(2 * _N_Z * _N_ADS)
    y0[_N_Z * _N_ADS:] = q_init.ravel()

    sol = solve_ivp(
        rhs, [0.0, t_max], y0,
        method='LSODA', events=breakthrough_event,
        dense_output=False, rtol=1e-4, atol=1e-7,
    )
    converged = len(sol.t_events[0]) > 0
    t_abs     = float(sol.t_events[0][0]) if converged else sol.t[-1]
    q_final   = sol.y[_N_Z * _N_ADS:, -1].reshape(_N_Z, _N_ADS)
    return t_abs, q_final, converged


def _calc_t_des(q_final, kfa, desorption_target):
    """q_final から脱着時間を逆算する（psa_system._calc_desorption_time と同一ロジック）。"""
    q_avg     = q_final.mean(axis=0)
    q_total_0 = q_avg.sum()
    if q_total_0 <= 0.0:
        return 0.0

    def residual(t):
        return (q_avg * np.exp(-kfa * t)).sum() / q_total_0 - desorption_target

    t_max = 1e5
    if residual(t_max) >= 0.0:
        return t_max
    return brentq(residual, 0.0, t_max, xtol=1.0, maxiter=200)


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------

def main():
    print("=" * 62)
    print("CSS 近似精度検証 (R-3)")
    print("=" * 62)
    print(f"  D_col={_DESIGN.D_col} m, L_bed={_DESIGN.L_bed} m, "
          f"desorption_target={_DESIGN.desorption_target}")
    print(f"  F_in={_FEED.F_in}")
    print(f"  P_in={_FEED.P_in/1e5:.1f} bar, T_in={_FEED.T_in:.1f} K")

    # 1. 現行近似 (use_css_approximation=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_approx = simulate_psa_system(
            _DESIGN, _FEED, PSAFixedParams(use_css_approximation=True)
        )
    t_abs_approx = res_approx.equipment.t_abs_sec
    if t_abs_approx == 0.0:
        print("\n[ERROR] 近似シミュレーションがペナルティを返しました。条件を見直してください。")
        return

    # 2. 清浄床 t_abs_clean を取得（参照値）
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_clean = simulate_psa_system(
            _DESIGN, _FEED, PSAFixedParams(use_css_approximation=False)
        )
    t_abs_clean = res_clean.equipment.t_abs_sec
    t_des_clean = res_clean.equipment.t_des_sec

    print(f"\n[近似]  t_abs_approx = {t_abs_approx:.1f} s")
    print(f"[清浄床] t_abs_clean  = {t_abs_clean:.1f} s  "
          f"(t_des/t_abs_clean = {t_des_clean/t_abs_clean:.2f})")

    # 3. フィード物性・パラメータ取得
    T_abs = 298.15
    C_feed_ads, _, F_non_C3_mol_s, Z = _calc_feed_state(
        _FEED.F_in, T_abs, _FEED.P_in
    )
    A_col = math.pi / 4.0 * _DESIGN.D_col ** 2
    u_0   = F_non_C3_mol_s * Z * R * T_abs / (_FEED.P_in * A_col)

    q_s    = np.array([PSA_LANGMUIR_PARAMS[k]['q_s'] for k in _ADS_ORDER])
    a_lang = np.array([PSA_LANGMUIR_PARAMS[k]['a']   for k in _ADS_ORDER])
    kfa    = np.array([PSA_KFA[k]                    for k in _ADS_ORDER])
    rho_b  = 600.0
    eps    = 0.4
    br     = 0.001
    t_max  = 7200.0

    # 4. サイクルループ（最大 5 サイクル、収束判定: 前サイクル比 1% 未満）
    print(f"\n[サイクルループ]")
    q_init = np.zeros((_N_Z, _N_ADS))
    t_abs_history = []

    for i in range(1, 6):
        t_abs_i, q_final_i, converged = _run_adsorption_with_q_init(
            C_feed_ads, u_0, _DESIGN.L_bed, rho_b, eps,
            kfa, q_s, a_lang, br, t_max, q_init,
        )
        if not converged:
            print(f"  サイクル {i}: 破過未検出（t_max 到達）、ループ中断")
            break

        t_des_i = _calc_t_des(q_final_i, kfa, _DESIGN.desorption_target)
        q_init  = q_final_i * np.exp(-kfa * t_des_i)

        t_abs_history.append(t_abs_i)
        print(f"  サイクル {i}: t_abs = {t_abs_i:.1f} s  (t_des = {t_des_i:.1f} s)")

        if i >= 2:
            delta = abs(t_abs_history[-1] - t_abs_history[-2]) / t_abs_history[-2]
            if delta < 0.01:
                print(f"  → {i} サイクルで収束（前サイクル比変化: {delta*100:.2f}%）")
                break

    if not t_abs_history:
        print("\n[ERROR] サイクルループが 1 サイクルも完了しませんでした。")
        return

    t_abs_true = t_abs_history[-1]

    # 5. 結果比較
    error_pct = (t_abs_approx - t_abs_true) / t_abs_true * 100.0
    print(f"\n{'─' * 62}")
    print(f"  近似 t_abs   : {t_abs_approx:.1f} s")
    print(f"  真の t_abs   : {t_abs_true:.1f} s  (サイクルループ収束値)")
    print(f"  誤差         : {error_pct:+.1f}%  "
          f"({'近似が短め' if error_pct < 0 else '近似が長め'})")

    if abs(error_pct) <= 10.0:
        print("  判定: OK（10% 以内） → 近似のまま最適化に進んで問題なし")
    else:
        n_approx = math.ceil(t_des_clean / t_abs_approx) + 1
        n_true   = math.ceil(t_des_clean / t_abs_true)   + 1
        print("  判定: NG（10% 超）   → 近似手法の見直しを検討すること")
        if n_approx != n_true:
            print(f"  ！ 必要塔数が変化: 近似={n_approx} 塔 vs 真値={n_true} 塔")
        else:
            print(f"  （必要塔数は両方 {n_approx} 塔で変化なし → 最適化への影響は限定的）")
    print(f"{'─' * 62}")


if __name__ == "__main__":
    main()
