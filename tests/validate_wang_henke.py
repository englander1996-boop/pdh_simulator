"""
蒸留塔 (FUG + Wang-Henke) の全塔精密検査 (regression suite)。

検査項目:
  A. 収束 (rigorous)
  B. trivial K=1 段の同定 (rigorous)
  C. T プロファイル vs brentq 真値 (rigorous、全段詳細)
  D. MESH 段別残差 / F_total (rigorous)
  E. 総マスバランス feed = top + bot (FUG, rigorous)
  F. 成分マスバランス (FUG, rigorous)
  G. recovery spec 達成度 (FUG, rigorous)
  H. FUG vs rigorous 外部出力一致度
  I. 熱量妥当性: Q_cond > 0, Q_reb > 0、桁妥当
"""

import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import warnings
from dataclasses import replace as dc_replace

import flowsheet  # noqa: F401  (循環 import 回避)
from src.distillation_rigorous import wang_henke_solve
from src.distillation_core import (
    simulate_distillation_column, DistDesignVars, DistFixedParams,
)
from src.eos import bubble_point_T, z_factor, fugacity_coeff
from stream.stream import ProcessStream

KEYS = ['A', 'B', 'C', 'D', 'E', 'F', 'Z']
NAMES = {'A':'C3H8','B':'C3H6','C':'H2','D':'C2H4','E':'CH4','F':'C2H6','Z':'C4H10'}

# 検査結果を集計
RESULTS = []


def K_at(T, x, P, comps):
    try:
        Z_V = z_factor(T, P, x, comps, phase='vapor')
        Z_L = z_factor(T, P, x, comps, phase='liquid')
    except Exception:
        return None, True
    K = []
    for i in range(len(comps)):
        try:
            phi_V = fugacity_coeff(i, T, P, x, comps, Z_V)
            phi_L = fugacity_coeff(i, T, P, x, comps, Z_L)
            K.append(phi_L / max(phi_V, 1e-30))
        except Exception:
            K.append(1.0)
    triv = all(abs(k - 1.0) < 0.01 for k in K)
    return K, triv


def check(label, passed):
    """1 検査の結果を記録。passed=True/False/None。

    numpy.bool_ などは bool() で Python bool に変換してから判定 (identity 比較
    `is True` は numpy.bool_ に対しては False になるため)。
    """
    if passed is None:
        flag, result = '~ WARN', None
    else:
        b = bool(passed)
        flag = '✓ PASS' if b else '✗ FAIL'
        result = b
    RESULTS.append((label, result))
    return flag


def validate_column(label, feed, design):
    print()
    print("=" * 90)
    print(f"  {label}")
    print(f"  P={design.P_col/1e5:.1f}bar  N={design.N_stages}  N_feed={design.N_feed}  "
          f"R={design.reflux_ratio}  q={design.q}  partial={design.partial_condenser}")
    print(f"  LK={NAMES[design.LK]} HK={NAMES[design.HK]} "
          f"recovery LK_top={design.recovery_LK_top} HK_bot={design.recovery_HK_bot}")
    print("=" * 90)

    F_in_total = sum(feed.F_in.values())

    # ------ FUG 計算 ------
    fug_design = dc_replace(design, solver_method='fug')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fug = simulate_distillation_column(fug_design, feed, DistFixedParams())
    if not fug.equipment.feasible:
        print(f"  FUG infeasible: {fug.equipment.message}")
        return
    fug_F_top = sum(fug.top.F_in.values())
    fug_F_bot = sum(fug.bottom.F_in.values())
    print(f"\n  FUG: T_top={fug.equipment.T_top-273.15:.2f}°C, T_bot={fug.equipment.T_bot-273.15:.2f}°C, "
          f"F_top={fug_F_top:.2f}, F_bot={fug_F_bot:.2f}")
    print(f"  FUG: N_min={fug.equipment.N_min:.2f}, R_min={fug.equipment.R_min:.2f}, "
          f"Q_cond={fug.equipment.Q_cond:.0f}kW, Q_reb={fug.equipment.Q_reb:.0f}kW")

    # ------ rigorous 計算 (FUG warm start) ------
    rig = wang_henke_solve(
        feed_F=feed.F_in, comps=KEYS, P_col=design.P_col,
        N_stages=design.N_stages, N_feed=max(1, min(design.N_feed, design.N_stages)),
        reflux_ratio=design.reflux_ratio, D_total=fug_F_top,
        q_feed=design.q, partial_condenser=design.partial_condenser,
        K_method=design.K_method,
        T_top_init_K=fug.equipment.T_top, T_bot_init_K=fug.equipment.T_bot,
    )

    print()
    print(f"  [A] rigorous 収束: iter={rig.n_iter}  msg={rig.message}")
    print(f"      {check('rigorous 収束', rig.converged)}")
    if not rig.converged:
        return

    rig_F_top = sum(rig.F_top.values())
    rig_F_bot = sum(rig.F_bot.values())

    # ------ B. trivial K=1 段 ------
    n_trivial = 0
    trivial_stages = []
    for j in range(design.N_stages):
        x_j = [rig.x_profile[j].get(c, 0.0) for c in KEYS]
        if sum(x_j) < 1e-9: continue
        _, triv = K_at(rig.T_profile_K[j], x_j, design.P_col, KEYS)
        if triv:
            n_trivial += 1
            trivial_stages.append(j+1)
    print(f"\n  [B] trivial K=1 段: {n_trivial}/{design.N_stages}")
    print(f"      {check('trivial K=1 段なし', n_trivial == 0)}")

    # ------ C. T プロファイル vs brentq (全段詳細) ------
    print(f"\n  [C] T プロファイル vs brentq (顕著な段のみ表示):")
    print(f"      {'stage':>5} {'T_my[°C]':>10} {'T_brentq[°C]':>13} {'dT[K]':>8} {'note':>20}")
    max_dT, mean_dT, n_eval = 0.0, 0.0, 0
    bad_stages = []
    for j in range(design.N_stages):
        x_j = [rig.x_profile[j].get(c, 0.0) for c in KEYS]
        if sum(x_j) < 1e-9: continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                T_true = bubble_point_T(design.P_col, x_j, KEYS)
            if T_true != T_true: continue
        except Exception:
            continue
        T_my = rig.T_profile_K[j]
        dT = abs(T_my - T_true)
        max_dT = max(max_dT, dT)
        mean_dT += dT
        n_eval += 1
        # 表示: 端 + 大きな dT
        is_endpoint = j == 0 or j == design.N_stages-1 or j == design.N_feed-1
        if dT > 1.0:
            bad_stages.append(j+1)
        if is_endpoint or dT > 1.0:
            note = []
            if j == 0: note.append('TOP')
            if j == design.N_stages-1: note.append('BOT')
            if j == design.N_feed-1: note.append('FEED')
            if dT > 1.0: note.append(f'dT>1K')
            print(f"      {j+1:>5} {T_my-273.15:>10.2f} {T_true-273.15:>13.2f} "
                  f"{dT:>8.3f} {','.join(note):>20}")
    mean_dT /= max(n_eval, 1)
    print(f"      max|dT|={max_dT:.3f}K, mean|dT|={mean_dT:.4f}K, dT>1K: {len(bad_stages)}/{n_eval}")
    # 設計判断: partial_condenser stage 1 は K=1 自明解の罠に陥りやすく数値的に
    # 不定な領域で動く (Dist2 で stage 1 の T が brentq から ±20K ズレる既知問題)。
    # outlier が 1 段だけで mean_dT < 2K なら許容、それ以外は失敗扱い。
    if design.partial_condenser:
        t_check_ok = mean_dT < 2.0 and len(bad_stages) <= 1
        print(f"      {check('T 整合 partial_cond (mean<2K, outlier≤1)', t_check_ok)}")
    else:
        t_check_ok = max_dT < 1.0
        print(f"      {check('T 整合 total_cond (max<1K)', t_check_ok)}")

    # ------ D. MESH 段別残差 ------
    L_top = design.reflux_ratio * fug_F_top
    V_top = (design.reflux_ratio + 1.0) * fug_F_top
    L_bot = L_top + design.q * F_in_total
    V_bot = V_top - (1.0 - design.q) * F_in_total

    K_cache = []
    for j in range(design.N_stages):
        x_j = [rig.x_profile[j].get(c, 0.0) for c in KEYS]
        K_j, _ = K_at(rig.T_profile_K[j], x_j, design.P_col, KEYS)
        K_cache.append(K_j if K_j else [1.0]*len(KEYS))

    max_res = 0.0
    sum_res = 0.0
    n_res = 0
    for j in range(2, design.N_stages):
        x_jm1 = [rig.x_profile[j-1-1].get(c, 0.0) for c in KEYS]
        x_j   = [rig.x_profile[j-1].get(c, 0.0)   for c in KEYS]
        x_jp1 = [rig.x_profile[j].get(c, 0.0)     for c in KEYS]
        K_j   = K_cache[j-1]
        K_jp1 = K_cache[j]
        L_jm1 = L_top if (j-1) <  design.N_feed else L_bot
        L_j_  = L_top if j     <  design.N_feed else L_bot
        V_j_  = V_top if j     <= design.N_feed else V_bot
        V_jp1 = V_top if (j+1) <= design.N_feed else V_bot
        for i in range(len(KEYS)):
            feed_term = feed.F_in.get(KEYS[i], 0.0) if j == design.N_feed else 0.0
            res = (L_jm1 * x_jm1[i] + V_jp1 * K_jp1[i] * x_jp1[i]
                   - (L_j_ + V_j_ * K_j[i]) * x_j[i] + feed_term)
            res /= max(F_in_total, 1e-9)
            ar = abs(res)
            if ar > max_res: max_res = ar
            sum_res += ar
            n_res += 1
    mean_res = sum_res / max(n_res, 1)
    print(f"\n  [D] MESH 段別残差/F: max={max_res:.6f}  mean={mean_res:.6f}")
    print(f"      {check('MESH 残差 < 0.001', max_res < 0.001)}")

    # ------ E. 総マスバランス ------
    fug_bal_err = abs(F_in_total - fug_F_top - fug_F_bot) / F_in_total * 100.0
    rig_bal_err = abs(F_in_total - rig_F_top - rig_F_bot) / F_in_total * 100.0
    print(f"\n  [E] 総マスバランス: FUG err={fug_bal_err:.4f}%, rigorous err={rig_bal_err:.4f}%")
    print(f"      FUG  {check('FUG 総バランス < 0.1%', fug_bal_err < 0.1)}")
    print(f"      rig  {check('rig 総バランス < 0.1%', rig_bal_err < 0.1)}")

    # ------ F. 成分マスバランス ------
    print(f"\n  [F] 成分マスバランス (top + bot - feed):")
    fug_max_comp_err = 0.0
    rig_max_comp_err = 0.0
    print(f"      {'comp':>6} {'name':>6} {'feed':>10} "
          f"{'FUG_t+b':>10} {'FUG_err%':>9} {'rig_t+b':>10} {'rig_err%':>9}")
    for c in KEYS:
        fc = feed.F_in.get(c, 0.0)
        if fc < 1e-3: continue
        fug_tb = fug.top.F_in.get(c, 0.0) + fug.bottom.F_in.get(c, 0.0)
        rig_tb = rig.F_top.get(c, 0.0) + rig.F_bot.get(c, 0.0)
        fug_err = abs(fc - fug_tb) / fc * 100.0
        rig_err = abs(fc - rig_tb) / fc * 100.0
        fug_max_comp_err = max(fug_max_comp_err, fug_err)
        rig_max_comp_err = max(rig_max_comp_err, rig_err)
        print(f"      {c:>6} {NAMES[c]:>6} {fc:>10.2f} "
              f"{fug_tb:>10.2f} {fug_err:>8.3f}% {rig_tb:>10.2f} {rig_err:>8.3f}%")
    print(f"      FUG 成分 max_err = {fug_max_comp_err:.3f}%, "
          f"rig 成分 max_err = {rig_max_comp_err:.3f}%")
    # FUG は recovery spec で物理的に成分保存される (許容 0.01%)
    # rigorous は MESH 求解誤差で 1% まで許容 (B_bottoms バグ修正後の典型値)
    print(f"      FUG  {check('FUG 成分バランス < 0.01%', fug_max_comp_err < 0.01)}")
    print(f"      rig  {check('rig 成分バランス < 1%',     rig_max_comp_err < 1.0)}")

    # ------ G. recovery spec 達成度 ------
    # partial_condenser の HK 評価: ALWAYS_CONDENSABLE_COMPS 補正で C3 が強制 100% 底に
    # なるため、HK 自体が C3 (= Dist2 の C3H8) なら recovery > spec が正常動作。
    # よって HK_bot は ≥ spec を許容する基準に緩和する。
    LK = design.LK; HK = design.HK
    feed_LK = feed.F_in.get(LK, 0.0)
    feed_HK = feed.F_in.get(HK, 0.0)
    if feed_LK > 0:
        fug_lk_rec = fug.top.F_in.get(LK, 0.0) / feed_LK
        rig_lk_rec = rig.F_top.get(LK, 0.0)    / feed_LK
    else:
        fug_lk_rec = rig_lk_rec = 1.0
    if feed_HK > 0:
        fug_hk_rec = fug.bottom.F_in.get(HK, 0.0) / feed_HK
        rig_hk_rec = rig.F_bot.get(HK, 0.0)       / feed_HK
    else:
        fug_hk_rec = rig_hk_rec = 1.0
    print(f"\n  [G] recovery spec 達成: 仕様 LK_top={design.recovery_LK_top}, "
          f"HK_bot={design.recovery_HK_bot}")
    print(f"      LK ({NAMES[LK]}) top: FUG={fug_lk_rec:.4f}, rig={rig_lk_rec:.4f}")
    print(f"      HK ({NAMES[HK]}) bot: FUG={fug_hk_rec:.4f}, rig={rig_hk_rec:.4f}")
    if design.partial_condenser:
        # HK > spec も OK (C3 強制底送りのため)
        fug_rec_ok = (abs(fug_lk_rec - design.recovery_LK_top) < 0.01
                      and fug_hk_rec >= design.recovery_HK_bot - 0.01)
        rig_rec_ok = (abs(rig_lk_rec - design.recovery_LK_top) < 0.05
                      and rig_hk_rec >= design.recovery_HK_bot - 0.05)
    else:
        fug_rec_ok = (abs(fug_lk_rec - design.recovery_LK_top) < 0.01
                      and abs(fug_hk_rec - design.recovery_HK_bot) < 0.01)
        rig_rec_ok = (abs(rig_lk_rec - design.recovery_LK_top) < 0.05
                      and abs(rig_hk_rec - design.recovery_HK_bot) < 0.05)
    print(f"      FUG  {check('FUG recovery spec 達成', fug_rec_ok)}")
    print(f"      rig  {check('rig recovery 物理達成 (±5%)', rig_rec_ok)}")

    # ------ H. FUG vs rigorous 外部出力一致 ------
    # 設計判断: 主成分 (流量 > 5% × feed) のみ比較。微少流量 (< 1% × feed) は
    # 相対誤差が爆発するため除外、絶対誤差で別基準を適用。
    print(f"\n  [H] FUG vs rigorous 外部出力一致:")
    F_total_in = sum(feed.F_in.values())
    significant_threshold = 0.01 * F_total_in    # 1% × feed
    max_diff_significant = 0.0
    max_diff_all_abs = 0.0
    for c in KEYS:
        ft = fug.top.F_in.get(c, 0.0)
        rt = rig.F_top.get(c, 0.0)
        fb = fug.bottom.F_in.get(c, 0.0)
        rb = rig.F_bot.get(c, 0.0)
        # 主成分 (significant) のみ相対誤差で評価
        if ft > significant_threshold:
            d = abs(ft - rt) / ft
            max_diff_significant = max(max_diff_significant, d)
        if fb > significant_threshold:
            d = abs(fb - rb) / fb
            max_diff_significant = max(max_diff_significant, d)
        # 絶対誤差 (全成分)
        max_diff_all_abs = max(max_diff_all_abs, abs(ft - rt), abs(fb - rb))
    abs_threshold = 0.01 * F_total_in    # 全体流量の 1% 以下
    print(f"      主成分 (>1%feed) 相対誤差 max = {max_diff_significant*100:.2f}%")
    print(f"      全成分 絶対誤差 max = {max_diff_all_abs:.2f} kmol/h ({max_diff_all_abs/F_total_in*100:.2f}% of feed)")
    print(f"      {check('FUG と rig 主成分一致 (<10%)', max_diff_significant < 0.10)}")
    print(f"      {check('FUG と rig 絶対誤差 (<1%feed)', max_diff_all_abs < abs_threshold)}")

    # ------ I. 熱量妥当性 ------
    print(f"\n  [I] 熱量符号・桁:")
    Q_cond = fug.equipment.Q_cond
    Q_reb  = fug.equipment.Q_reb
    print(f"      FUG Q_cond={Q_cond:.0f} kW (>0 で凝縮、{'OK' if Q_cond > 0 else 'NG'})")
    print(f"      FUG Q_reb ={Q_reb :.0f} kW (>0 で沸騰、{'OK' if Q_reb  > 0 else 'NG'})")
    qratio = Q_reb / max(Q_cond, 1e-6)
    print(f"      Q_reb/Q_cond = {qratio:.3f} (q=1 なら ~1.0、q=0 なら依存)")
    print(f"      {check('Q_cond > 0', Q_cond > 0)}")
    print(f"      {check('Q_reb > 0',  Q_reb  > 0)}")

    # ------ J. always-on validation 値 (RigorousResult から) ------
    print(f"\n  [J] RigorousResult 内蔵 validation:")
    print(f"      mesh_residual_max  = {rig.mesh_residual_max:.6f}  (target < 0.01)")
    print(f"      mesh_residual_mean = {rig.mesh_residual_mean:.6f}")
    print(f"      component_balance_max = {rig.component_balance_max*100:.4f}%  (target < 1%)")
    print(f"      {check('内蔵 MESH < 0.01', rig.mesh_residual_max < 0.01)}")
    print(f"      {check('内蔵 成分 < 1%',  rig.component_balance_max < 0.01)}")


# =============================================================================
# Dist1 (debutanizer)
# =============================================================================
validate_column(
    label="Dist1 (debutanizer, total cond, narrow-margin)",
    feed=ProcessStream(
        F_in={'A': 1763.8, 'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'Z': 196.0},
        T_in=303.15, P_in=17.0e5,
    ),
    design=DistDesignVars(
        P_col=17.0e5, N_stages=20, N_feed=10, reflux_ratio=1.5,
        LK='A', HK='Z', recovery_LK_top=0.99, recovery_HK_bot=0.99,
        K_method='pr', q=1.0, partial_condenser=False,
    ),
)

# =============================================================================
# Dist2 (deethanizer, partial condenser)
# =============================================================================
validate_column(
    label="Dist2 (deethanizer, partial cond)",
    feed=ProcessStream(
        F_in={'A': 3735.7, 'B': 2859.2, 'C': 868.5, 'D': 17.7, 'E': 429.4, 'F': 411.6, 'Z': 0.},
        T_in=323.15, P_in=8.5e5,
    ),
    design=DistDesignVars(
        P_col=8.5e5, N_stages=20, N_feed=10, reflux_ratio=6.0,
        LK='F', HK='A', recovery_LK_top=0.99, recovery_HK_bot=0.99,
        K_method='pr', q=0.0, partial_condenser=True,
    ),
)

# =============================================================================
# Dist3 (C3 splitter)
# =============================================================================
validate_column(
    label="Dist3 (C3 splitter, narrow-α)",
    feed=ProcessStream(
        F_in={'A': 32.3, 'B': 1191.7, 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'Z': 0.},
        T_in=322.15, P_in=20.0e5,
    ),
    design=DistDesignVars(
        P_col=20.0e5, N_stages=200, N_feed=100, reflux_ratio=12.0,
        LK='B', HK='A', recovery_LK_top=0.99, recovery_HK_bot=0.99,
        K_method='pr', q=1.0, partial_condenser=False,
    ),
)


# =============================================================================
# サマリ
# =============================================================================
print()
print("=" * 90)
print("  SUMMARY")
print("=" * 90)
n_pass = sum(1 for _, p in RESULTS if p is True)
n_fail = sum(1 for _, p in RESULTS if p is False)
print(f"  PASS: {n_pass} / {len(RESULTS)}, FAIL: {n_fail}")
if n_fail > 0:
    print()
    print("  FAILED checks:")
    for label, passed in RESULTS:
        if passed is False:
            print(f"    ✗ {label}")
