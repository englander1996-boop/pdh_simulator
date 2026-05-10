"""
exp2.py — 蒸留塔単独比較 (FUG vs rigorous VLE)

各塔を **フローシートから切り離して** 単独で動かし、FUG と rigorous の結果を
並べて比較する。exp1 のようにリサイクルループを回す必要がないので、
塔単位の検査が速い (~10 秒程度で 1 塔)。

使い方:
  - 下の "実験で振る設計変数" セクションを編集
  - .\.venv\Scripts\python.exe exp/exp2.py で実行
  - 各塔について以下が表示される:
      * FUG / rigorous の T_top, T_bot, F_top, F_bot
      * 段別 T プロファイル (rigorous の収束プロファイル)
      * 成分マスバランス (feed = top + bot per component)
      * recovery 達成度 (LK_top, HK_bot)
      * Q_cond, Q_reb, A_cond, A_reb
      * always-on validation 値 (MESH 残差、成分バランス)

調査対象を絞りたい場合は RUN_DIST1/2/3 の True/False を切替。
"""

import os
import sys

# Windows コンソール対応
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

import warnings
from dataclasses import replace as dc_replace

# 循環 import 回避
import flowsheet  # noqa

from src.distillation_core import (
    simulate_distillation_column, DistDesignVars, DistFixedParams,
)
from stream.stream import ProcessStream

# ===========================================================================
#  どの塔を比較するか (True/False で切替)
# ===========================================================================
RUN_DIST1 = True       # 脱ブタン塔
RUN_DIST2 = True       # 脱エタン塔 (partial cond)
RUN_DIST3 = True       # C3 スプリッタ (narrow-α)


# ===========================================================================
#  実験で振る設計変数 (ここを書き換えて再実行)
# ===========================================================================

# === Dist1 ====================================================================
D1_P_col        = 17.0e5      # Pa
D1_N_stages     = 20
D1_N_feed       = 10
D1_R            = 1.5
D1_LK           = 'A'         # C3H8
D1_HK           = 'Z'         # C4H10
D1_q            = 1.0
D1_partial_cond = False
D1_feed         = ProcessStream(
    F_in={'A': 1763.8, 'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'Z': 196.0},
    T_in=303.15, P_in=17.0e5,
)

# === Dist2 ====================================================================
D2_P_col        = 8.5e5
D2_N_stages     = 20
D2_N_feed       = 10
D2_R            = 6.0
D2_LK           = 'F'         # C2H6
D2_HK           = 'A'         # C3H8
D2_q            = 0.0         # 飽和気
D2_partial_cond = True
D2_feed         = ProcessStream(
    F_in={'A': 3735.7, 'B': 2859.2, 'C': 868.5, 'D': 17.7, 'E': 429.4, 'F': 411.6, 'Z': 0.},
    T_in=323.15, P_in=8.5e5,
)

# === Dist3 ====================================================================
D3_P_col        = 20.0e5
D3_N_stages     = 200
D3_N_feed       = 100
D3_R            = 12.0
D3_LK           = 'B'         # C3H6
D3_HK           = 'A'         # C3H8
D3_q            = 1.0
D3_partial_cond = False
D3_feed         = ProcessStream(
    F_in={'A': 32.3, 'B': 1191.7, 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'Z': 0.},
    T_in=322.15, P_in=20.0e5,
)


# ===========================================================================
#  比較ロジック
# ===========================================================================

KEYS = ['A', 'B', 'C', 'D', 'E', 'F', 'Z']
NAMES = {'A':'C3H8','B':'C3H6','C':'H2','D':'C2H4','E':'CH4','F':'C2H6','Z':'C4H10'}


def compare(label, feed, P_col, N_stages, N_feed, R, LK, HK, q, partial_cond):
    """1 塔について FUG と rigorous を実行して比較表示。"""
    print()
    print("=" * 84)
    print(f"  {label}")
    print(f"  P={P_col/1e5:.1f}bar  N={N_stages}  N_feed={N_feed}  R={R}  q={q}  "
          f"partial={partial_cond}  LK={NAMES[LK]} HK={NAMES[HK]}")
    print("=" * 84)
    print(f"  feed: {{ {', '.join(f'{NAMES[c]}:{feed.F_in[c]:.1f}' for c in KEYS if feed.F_in.get(c,0)>1e-6)} }}")
    print(f"  feed total: {sum(feed.F_in.values()):.2f} kmol/h, T={feed.T_in-273.15:.1f}°C")

    # 共通 design (LK/HK/recovery)
    base_design = DistDesignVars(
        P_col=P_col, N_stages=N_stages, N_feed=N_feed, reflux_ratio=R,
        LK=LK, HK=HK, recovery_LK_top=0.99, recovery_HK_bot=0.99,
        K_method='pr', q=q, partial_condenser=partial_cond,
    )

    # --- FUG ---
    import time
    fug_design = dc_replace(base_design, solver_method='fug')
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fug = simulate_distillation_column(fug_design, feed, DistFixedParams())
    fug_time = time.perf_counter() - t0
    if not fug.equipment.feasible:
        print(f"  FUG infeasible: {fug.equipment.message}")
        return

    # --- rigorous ---
    rig_design = dc_replace(base_design, solver_method='rigorous')
    t0 = time.perf_counter()
    rig = simulate_distillation_column(rig_design, feed, DistFixedParams())
    rig_time = time.perf_counter() - t0

    # ---- 表示 ----
    print()
    print(f"  {'metric':<24} {'FUG':>14} {'rigorous':>14} {'diff':>12}")
    print(f"  {'-'*24} {'-'*14} {'-'*14} {'-'*12}")

    def row(name, fug_v, rig_v, fmt='.3f'):
        diff = rig_v - fug_v
        print(f"  {name:<24} {fug_v:>14{fmt}} {rig_v:>14{fmt}} {diff:>12{fmt}}")

    row('実行時間 [秒]',           fug_time, rig_time, '.3f')
    row('T_top [°C]',              fug.equipment.T_top - 273.15, rig.equipment.T_top - 273.15, '.2f')
    row('T_bot [°C]',              fug.equipment.T_bot - 273.15, rig.equipment.T_bot - 273.15, '.2f')
    row('F_top [kmol/h]',          sum(fug.top.F_in.values()),    sum(rig.top.F_in.values()),    '.2f')
    row('F_bot [kmol/h]',          sum(fug.bottom.F_in.values()), sum(rig.bottom.F_in.values()), '.2f')
    row('Q_cond [kW]',             fug.equipment.Q_cond, rig.equipment.Q_cond, '.0f')
    row('Q_reb [kW]',              fug.equipment.Q_reb,  rig.equipment.Q_reb,  '.0f')
    row('A_cond [m²]',             fug.equipment.A_cond_m2, rig.equipment.A_cond_m2, '.0f')
    row('A_reb [m²]',              fug.equipment.A_reb_m2,  rig.equipment.A_reb_m2,  '.0f')
    row('CAPEX 計 [億円]',         fug.equipment.CAPEX, rig.equipment.CAPEX, '.3f')

    # 成分流量比較
    print()
    print(f"  成分別 top/bot 流量 [kmol/h]:")
    print(f"  {'comp':<6} {'feed':>10} {'FUG_top':>10} {'rig_top':>10} {'FUG_bot':>10} {'rig_bot':>10}")
    for c in KEYS:
        if feed.F_in.get(c, 0.0) < 1e-3 \
           and fug.top.F_in.get(c, 0.0) < 1e-3 \
           and rig.top.F_in.get(c, 0.0) < 1e-3 \
           and fug.bottom.F_in.get(c, 0.0) < 1e-3 \
           and rig.bottom.F_in.get(c, 0.0) < 1e-3:
            continue
        print(f"  {NAMES[c]:<6} {feed.F_in.get(c, 0.0):>10.3f} "
              f"{fug.top.F_in.get(c, 0.0):>10.3f} {rig.top.F_in.get(c, 0.0):>10.3f} "
              f"{fug.bottom.F_in.get(c, 0.0):>10.3f} {rig.bottom.F_in.get(c, 0.0):>10.3f}")

    # 成分マスバランス検査
    print()
    print(f"  成分マスバランス (feed - top - bot):")
    for c in KEYS:
        f_c = feed.F_in.get(c, 0.0)
        if f_c < 1e-3: continue
        fug_err = abs(f_c - fug.top.F_in.get(c, 0.0) - fug.bottom.F_in.get(c, 0.0)) / f_c * 100
        rig_err = abs(f_c - rig.top.F_in.get(c, 0.0) - rig.bottom.F_in.get(c, 0.0)) / f_c * 100
        print(f"  {NAMES[c]:<6} FUG err={fug_err:>7.4f}%  rig err={rig_err:>7.4f}%")

    # recovery
    feed_LK = feed.F_in.get(LK, 0.0)
    feed_HK = feed.F_in.get(HK, 0.0)
    if feed_LK > 0 and feed_HK > 0:
        fug_lk = fug.top.F_in.get(LK, 0.0) / feed_LK
        fug_hk = fug.bottom.F_in.get(HK, 0.0) / feed_HK
        rig_lk = rig.top.F_in.get(LK, 0.0) / feed_LK
        rig_hk = rig.bottom.F_in.get(HK, 0.0) / feed_HK
        print()
        print(f"  recovery (target spec=0.99):")
        print(f"  LK ({NAMES[LK]}) top: FUG={fug_lk:.4f}  rig={rig_lk:.4f}")
        print(f"  HK ({NAMES[HK]}) bot: FUG={fug_hk:.4f}  rig={rig_hk:.4f}")

    # FUG の bubble_point_T 偽根チェック (T で f(T) ≈ 0 か)
    # 注意: partial_cond の T_top は x_top_C (= F_top excluding non-condensable)
    # の bubble point。FUG 内部と同じ x を使わないと誤判定する。
    print()
    print(f"  FUG bubble_point_T 偽根チェック (T で f(T) = sum(K x) - 1 を計算):")
    from src.eos import z_factor as _zf, fugacity_coeff as _fc
    from src.distillation_core import NON_CONDENSABLE_COMPS as _NON_CONDENSABLE

    def check_bubble(T_K, x_dict, P, label):
        x_list = [x_dict.get(c, 0.0) for c in KEYS]
        s_x = sum(x_list)
        if s_x < 1e-9:
            print(f"    {label}: x が空、スキップ")
            return
        x_norm = [v/s_x for v in x_list]
        try:
            Z_V = _zf(T_K, P, x_norm, KEYS, phase='vapor')
            Z_L = _zf(T_K, P, x_norm, KEYS, phase='liquid')
            K = []
            for i in range(7):
                phi_V = _fc(i, T_K, P, x_norm, KEYS, Z_V)
                phi_L = _fc(i, T_K, P, x_norm, KEYS, Z_L)
                K.append(phi_L / max(phi_V, 1e-30))
            f = sum(K[i] * x_norm[i] for i in range(7)) - 1.0
            n_triv = sum(1 for k in K if abs(k - 1.0) < 0.01)
            tag = ""
            if abs(f) > 0.01:
                tag = "  ⚠ 偽根の疑い (f が 0 から大きく離れてる)"
            elif n_triv == 7:
                tag = "  ⚠ trivial K=1 領域"
            print(f"    {label}: T={T_K-273.15:.2f}°C, f(T)={f:+.4f}, "
                  f"K=1 段数={n_triv}/7{tag}")
        except Exception as e:
            print(f"    {label}: 計算失敗 ({e})")

    # T_top: partial_cond なら x_top_C (excl. non-condensable)、total_cond なら x_top
    if partial_cond:
        x_top_C = {c: v for c, v in fug.top.F_in.items() if c not in _NON_CONDENSABLE}
        check_bubble(fug.equipment.T_top, x_top_C, P_col, "FUG T_top, x=x_top_C (partial)")
    else:
        check_bubble(fug.equipment.T_top, fug.top.F_in, P_col, "FUG T_top, x=F_top")
    check_bubble(fug.equipment.T_bot, fug.bottom.F_in, P_col, "FUG T_bot, x=F_bot")

    # 熱バランス (定常): Q_reb - Q_cond + Q_preheat ≈ 産物 H - Feed H
    # 簡易チェック: Q_cond ≈ Q_reb (q=1, 飽和液 feed の場合) または依存
    Q_cond_diff = abs(rig.equipment.Q_cond - fug.equipment.Q_cond)
    Q_reb_diff = abs(rig.equipment.Q_reb - fug.equipment.Q_reb)
    print(f"\n  熱量比較: Q_cond diff={Q_cond_diff:.0f} kW ({Q_cond_diff/max(fug.equipment.Q_cond,1)*100:.1f}%), "
          f"Q_reb diff={Q_reb_diff:.0f} kW ({Q_reb_diff/max(fug.equipment.Q_reb,1)*100:.1f}%)")


# =============================================================================
# 実行
# =============================================================================
if __name__ == '__main__':
    if RUN_DIST1:
        compare("Dist1 (debutanizer)",
                D1_feed, D1_P_col, D1_N_stages, D1_N_feed, D1_R,
                D1_LK, D1_HK, D1_q, D1_partial_cond)

    if RUN_DIST2:
        compare("Dist2 (deethanizer, partial cond)",
                D2_feed, D2_P_col, D2_N_stages, D2_N_feed, D2_R,
                D2_LK, D2_HK, D2_q, D2_partial_cond)

    if RUN_DIST3:
        compare("Dist3 (C3 splitter, narrow-α)",
                D3_feed, D3_P_col, D3_N_stages, D3_N_feed, D3_R,
                D3_LK, D3_HK, D3_q, D3_partial_cond)
