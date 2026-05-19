"""Dist2 (rigorous, partial cond) を exp1 入口で単独実行し、
stage 1 (塔頂) の y/x/K と F_top の精密値を引き出して点検する。

問題意識:
  - display は v > 0.01 kmol/h で打ち切り → C3H6 in top が真に 0 か微小残量か区別不能。
  - rigorous は D_total を FUG から固定値で受け取っており、FUG では C3 ALWAYS_CONDENSABLE
    強制移動が掛かっている。rigorous で内部的にどう振り分けたかを確認したい。
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings

# 循環 import 回避: flowsheet を先に読ませる (exp1.py の流儀)
import flowsheet  # noqa
from src.distillation_core import ColumnTunables, DistFixedParams
from src.distillation_rigorous import wang_henke_solve
from stream.stream import ProcessStream
from units.separators.column2.column2 import simulate_column2

# exp1 の Dist2 フィード (Desuper → Dist2):
F_feed_kmolh = {
    'A': 4483.2,   # C3H8
    'B': 2753.0,   # C3H6
    'C': 1014.4,   # H2
    'D': 16.8,     # C2H4
    'E': 209.3,    # CH4
    'F': 192.4,    # C2H6
}
feed = ProcessStream(
    F_in=F_feed_kmolh,
    T_in=50.0 + 273.15,
    P_in=7.5e5,
)

tunables = ColumnTunables(
    P_col        = 7.4812e5,
    N_stages     = 23,
    N_feed       = 12,          # rigorous では無視されるはず
    reflux_ratio = 9.37202,
    solver_method='rigorous',
    recovery_LK_top=None,
    recovery_HK_bot=None,
)

with warnings.catch_warnings():
    warnings.simplefilter("default")
    res = simulate_column2(feed, tunables=tunables)

print("=" * 70)
print("  simulate_column2 (rigorous) 出力")
print("=" * 70)
print(f"T_top = {res.equipment.T_top - 273.15:.2f} °C")
print(f"T_bot = {res.equipment.T_bot - 273.15:.2f} °C")
print(f"N_min = {res.equipment.N_min:.2f}")
print(f"R_min = {res.equipment.R_min:.3f}")
print(f"N_feed_kirkbride = {res.equipment.N_feed_kirkbride}")
print(f"cond utility = {res.equipment.cond_utility_name} ({res.equipment.cond_utility_jpy_per_GJ} 円/GJ)")
print(f"message: {res.equipment.message}")
print()

print("--- F_top (塔頂、partial cond で vapor distillate) ---")
F_top_total = sum(res.top.F_in.values())
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    v = res.top.F_in.get(c, 0.0)
    pct = 100 * v / F_top_total if F_top_total > 0 else 0
    name = {'A': 'C3H8', 'B': 'C3H6', 'C': 'H2', 'D': 'C2H4', 'E': 'CH4', 'F': 'C2H6'}[c]
    feed_v = F_feed_kmolh.get(c, 0.0)
    rec_top = 100 * v / feed_v if feed_v > 0 else float('nan')
    print(f"  {c}({name:5s}): F_top={v:12.6f} kmol/h ({pct:6.3f}%)   "
          f"recovery_to_top={rec_top:.4f}%")
print(f"  TOTAL F_top = {F_top_total:.4f} kmol/h")

print()
print("--- F_bot (塔底) ---")
F_bot_total = sum(res.bottom.F_in.values())
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    v = res.bottom.F_in.get(c, 0.0)
    pct = 100 * v / F_bot_total if F_bot_total > 0 else 0
    name = {'A': 'C3H8', 'B': 'C3H6', 'C': 'H2', 'D': 'C2H4', 'E': 'CH4', 'F': 'C2H6'}[c]
    feed_v = F_feed_kmolh.get(c, 0.0)
    rec_bot = 100 * v / feed_v if feed_v > 0 else float('nan')
    print(f"  {c}({name:5s}): F_bot={v:12.6f} kmol/h ({pct:6.3f}%)   "
          f"recovery_to_bot={rec_bot:.4f}%")
print(f"  TOTAL F_bot = {F_bot_total:.4f} kmol/h")

print()
print("=" * 70)
print("  独立に wang_henke_solve を呼んで stage 1 の y/x/K を見る")
print("=" * 70)

rig = wang_henke_solve(
    feed_F            = F_feed_kmolh,
    comps             = ['A', 'B', 'C', 'D', 'E', 'F'],
    P_col             = 7.4812e5,
    N_stages          = 23,
    N_feed            = res.equipment.N_feed_kirkbride or 12,
    reflux_ratio      = 9.37202,
    D_total           = F_top_total,   # FUG が決めた D を流用
    q_feed            = 0.0,
    partial_condenser = True,
    K_method          = 'pr',
)

print(f"converged = {rig.converged} ({rig.n_iter} iter)")
print(f"message = {rig.message}")
print(f"mesh_residual_max = {rig.mesh_residual_max:.3e}")
print(f"mesh_residual_mean = {rig.mesh_residual_mean:.3e}")
print(f"component_balance_max = {rig.component_balance_max:.3e}")
print()
print("--- stage 1 (塔頂 partial cond) ---")
print(f"  T_1 = {rig.T_profile_K[0] - 273.15:.2f} °C")
print(f"  V_top = {rig.V_top_kmolh:.2f} kmol/h, L_top (reflux) = {rig.L_top_kmolh:.2f} kmol/h")
print(f"  {'comp':>5} | {'x_1':>12} | {'y_1':>12} | {'K_1':>10} | "
      f"{'D*y_1 [kmol/h]':>16}")
D_used = sum(rig.F_top.values())
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    x1 = rig.x_profile[0].get(c, 0.0)
    y1 = rig.y_profile[0].get(c, 0.0)
    K1 = rig.K_profile[0].get(c, 0.0)
    F1 = D_used * y1
    print(f"  {c:>5} | {x1:12.6e} | {y1:12.6e} | {K1:10.4f} | {F1:16.6f}")

print()
print("--- stage 23 (塔底 reboiler) ---")
N = 23
print(f"  T_N = {rig.T_profile_K[-1] - 273.15:.2f} °C")
print(f"  V_bot = {rig.V_bot_kmolh:.2f} kmol/h, L_bot (downflow) = {rig.L_bot_kmolh:.2f} kmol/h")
print(f"  {'comp':>5} | {'x_N':>12} | {'y_N':>12} | {'K_N':>10}")
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    xn = rig.x_profile[-1].get(c, 0.0)
    yn = rig.y_profile[-1].get(c, 0.0)
    Kn = rig.K_profile[-1].get(c, 0.0)
    print(f"  {c:>5} | {xn:12.6e} | {yn:12.6e} | {Kn:10.4f}")

print()
print("--- T プロファイル (全段) ---")
for j, T_K in enumerate(rig.T_profile_K, start=1):
    marker = ""
    if j == 1:
        marker = " <- 塔頂 (partial cond)"
    elif j == (res.equipment.N_feed_kirkbride or 12):
        marker = f" <- フィード段 ({j})"
    elif j == N:
        marker = " <- 塔底 (reboiler)"
    print(f"  stage {j:2d}: T = {T_K - 273.15:7.2f} °C{marker}")

print()
print("=" * 70)
print("  C3H6 マスバランス")
print("=" * 70)
F_feed_B = F_feed_kmolh['B']
F_top_B  = res.top.F_in.get('B', 0.0)
F_bot_B  = res.bottom.F_in.get('B', 0.0)
print(f"  Feed  C3H6: {F_feed_B:.4f} kmol/h")
print(f"  Top   C3H6: {F_top_B:.4e}  kmol/h ({100 * F_top_B / F_feed_B:.4f}%)")
print(f"  Bot   C3H6: {F_bot_B:.4f} kmol/h ({100 * F_bot_B / F_feed_B:.4f}%)")
print(f"  Sum (top+bot) = {F_top_B + F_bot_B:.4f}  (Δ={F_feed_B - F_top_B - F_bot_B:+.4e})")
