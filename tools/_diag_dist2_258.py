"""Dist2 trial #258 (P=7.55bar, R=5.22, N=26) を rigorous で診断。

BO best (trial 258) は FUG で TAC=455.4 だが rigorous 再評価で Wang-Henke が
dT_max=7K で stall (max_iter=200)。原因が
  (i) Wegstein 振動 (= 数値) なのか
  (ii) FUG の D_total が物理的に過大で MESH に解なし (= 物理) なのか
を切り分けるため、収束ログを段階的に出す。
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings
import numpy as np

import flowsheet  # noqa: 循環 import 回避
from src.distillation_core import ColumnTunables
from src.distillation_rigorous import wang_henke_solve
from stream.stream import ProcessStream
from units.separators.column2.column2 import simulate_column2

# trial #258 設計値
P_dist2_Pa   = 7.54785e5
N_dist2      = 26
reflux_dist2 = 5.21708

# Dist2 への現実的なフィードを構築するには反応器系を回す必要があるため、
# exp1 で出た値の近似を使う。trial #258 は F_C3H8_fresh=1515.7 kmol/h なので
# 反応器出口・Comp2 後の Dist2 入口は大体次のオーダー (FUG 楽観値):
F_feed = {
    'A': 4750.0,   # C3H8 (recycle 多めで増)
    'B': 2900.0,   # C3H6
    'C': 1075.0,   # H2
    'D': 18.0,     # C2H4
    'E': 220.0,    # CH4
    'F': 203.0,    # C2H6
}
feed = ProcessStream(F_in=F_feed, T_in=50.0 + 273.15, P_in=P_dist2_Pa)

# まず simulate_column2 (rigorous) を試す
print("=" * 70)
print(f"  Dist2 #258: P={P_dist2_Pa/1e5:.2f}bar, N={N_dist2}, R={reflux_dist2:.3f}")
print("=" * 70)
tunables = ColumnTunables(
    P_col=P_dist2_Pa, N_stages=N_dist2,
    N_feed=1, reflux_ratio=reflux_dist2,
    solver_method='rigorous',
)
print("[1] simulate_column2 (rigorous, strict=on)")
try:
    res = simulate_column2(feed, tunables=tunables)
    print(f"  ✓ 収束: T_top={res.equipment.T_top-273.15:.2f}°C, "
          f"T_bot={res.equipment.T_bot-273.15:.2f}°C")
    print(f"  C3H6 top: {res.top.F_in.get('B', 0):.2f} kmol/h "
          f"({100*res.top.F_in.get('B', 0)/F_feed['B']:.2f}%)")
except Exception as e:
    print(f"  ✗ 失敗: {type(e).__name__}: {e}")
print()

# 同じ系を fug でも試す (FUG が出す D_total を確認)
print("[2] simulate_column2 (fug)")
tunables_fug = ColumnTunables(
    P_col=P_dist2_Pa, N_stages=N_dist2,
    N_feed=1, reflux_ratio=reflux_dist2,
    solver_method='fug',
)
with warnings.catch_warnings(record=True) as w_list:
    warnings.simplefilter("always")
    res_fug = simulate_column2(feed, tunables=tunables_fug)
    for w in w_list:
        print(f"  warning: {w.message}")
print(f"  T_top={res_fug.equipment.T_top-273.15:.2f}°C, "
      f"T_bot={res_fug.equipment.T_bot-273.15:.2f}°C")
print(f"  N_min={res_fug.equipment.N_min:.2f}, R_min={res_fug.equipment.R_min:.3f}")
print(f"  R/R_min = {reflux_dist2 / res_fug.equipment.R_min:.3f}")
print(f"  N/N_min = {N_dist2 / res_fug.equipment.N_min:.3f}")
D_fug = sum(res_fug.top.F_in.values())
print(f"  D_fug = {D_fug:.2f} kmol/h")
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    v = res_fug.top.F_in.get(c, 0)
    f_c = F_feed.get(c, 0)
    if f_c > 0:
        print(f"    {c} top: {v:9.2f} kmol/h ({100*v/f_c:6.2f}% of feed)")
print()

# 異なる R で rigorous がどう挙動するか
print("[3] R をスイープ (R_min × margin) — rigorous 単独直呼び")
R_min = res_fug.equipment.R_min
for R_factor in [1.05, 1.10, 1.20, 1.40, 1.60, 2.00, 3.00]:
    R_try = R_min * R_factor
    rig = wang_henke_solve(
        feed_F=F_feed, comps=['A','B','C','D','E','F'],
        P_col=P_dist2_Pa, N_stages=N_dist2,
        N_feed=max(1, min(res_fug.equipment.N_feed_kirkbride, N_dist2)),
        reflux_ratio=R_try, D_total=D_fug,  # FUG D_total を流用
        q_feed=0.0, partial_condenser=True,
        T_top_init_K=res_fug.equipment.T_top,
        T_bot_init_K=res_fug.equipment.T_bot,
        max_iter=500,
    )
    cb = rig.component_balance_max
    mesh = rig.mesh_residual_max
    flag = "✓" if rig.converged and cb < 0.01 and mesh < 0.01 else "✗"
    print(f"  R={R_try:6.3f} (×{R_factor}): {flag} conv={rig.converged} "
          f"iter={rig.n_iter:3d} cb={cb:.2e} mesh={mesh:.2e} | {rig.message[:60]}")

# 違う D_total で試す (FUG が D 過大なら適切な D で rigorous が解けるか)
print()
print("[4] D_total をスイープ (R=5.22 固定)")
for D_factor in [0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]:
    D_try = D_fug * D_factor
    rig = wang_henke_solve(
        feed_F=F_feed, comps=['A','B','C','D','E','F'],
        P_col=P_dist2_Pa, N_stages=N_dist2,
        N_feed=max(1, min(res_fug.equipment.N_feed_kirkbride, N_dist2)),
        reflux_ratio=reflux_dist2, D_total=D_try,
        q_feed=0.0, partial_condenser=True,
        T_top_init_K=res_fug.equipment.T_top,
        T_bot_init_K=res_fug.equipment.T_bot,
        max_iter=500,
    )
    cb = rig.component_balance_max
    mesh = rig.mesh_residual_max
    flag = "✓" if rig.converged and cb < 0.01 and mesh < 0.01 else "✗"
    print(f"  D={D_try:7.2f} (×{D_factor}): {flag} conv={rig.converged} "
          f"iter={rig.n_iter:3d} cb={cb:.2e} mesh={mesh:.2e}")
