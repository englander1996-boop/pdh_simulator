"""HYSYS 条件 (N=60, P=8.5 bar, R=10) を当シミュレータで実走させて比較。
HYSYS 報告: 99% 分離 (C3H6 漏れ ~1%)
当方目標: 同条件で何 % が出るかを正確に見て、コードバグの可能性を切り分ける。
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings
import math

import flowsheet  # noqa
from src.distillation_core import ColumnTunables
from src.distillation_rigorous import wang_henke_solve
from stream.stream import ProcessStream
from units.separators.column2.column2 import simulate_column2

# HYSYS 同等のフィード
F_feed = {
    'A': 4750.0, 'B': 2900.0, 'C': 1075.0,
    'D': 18.0,   'E': 220.0,  'F': 203.0,
}
F_total = sum(F_feed.values())

print("=" * 100)
print("  HYSYS 同条件で当シミュレータが返す結果 (N=60, P=8.5 bar, R=10, partial cond)")
print("=" * 100)
print(f"  feed: {[(c, F_feed[c]) for c in ['A','B','C','D','E','F']]}")
print(f"  total feed = {F_total:.1f} kmol/h")
print()

# (1) simulate_column2 経由 (rigorous + 全パイプライン)
print("[1] simulate_column2 (rigorous, full pipeline)")
feed = ProcessStream(F_in=F_feed, T_in=50.0 + 273.15, P_in=8.5e5)
tun = ColumnTunables(
    P_col=8.5e5, N_stages=60,
    N_feed=1, reflux_ratio=10.0,
    solver_method='rigorous',
)
with warnings.catch_warnings():
    warnings.simplefilter("default")
    res = simulate_column2(feed, tunables=tun)

print(f"  T_top={res.equipment.T_top-273.15:.2f}°C, T_bot={res.equipment.T_bot-273.15:.2f}°C")
print(f"  N_min={res.equipment.N_min:.2f}, R_min={res.equipment.R_min:.3f}")
print(f"  R/R_min={10/res.equipment.R_min:.3f}, N/N_min={60/res.equipment.N_min:.3f}")
print(f"  N_feed_kirkbride = {res.equipment.N_feed_kirkbride}")
print()
print(f"  --- F_top (塔頂, partial cond で vapor distillate) ---")
F_top_total = sum(res.top.F_in.values())
names = {'A': 'C3H8', 'B': 'C3H6', 'C': 'H2', 'D': 'C2H4', 'E': 'CH4', 'F': 'C2H6'}
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    v = res.top.F_in.get(c, 0.0)
    f_c = F_feed[c]
    rec_top = 100 * v / f_c if f_c > 0 else float('nan')
    print(f"    {c}({names[c]:5s}): {v:10.4f} kmol/h ({100*v/F_top_total:6.2f}%)   "
          f"to_top={rec_top:7.4f}%")
print(f"    TOTAL F_top = {F_top_total:.2f}")
print()
print(f"  --- F_bot (塔底) ---")
F_bot_total = sum(res.bottom.F_in.values())
for c in ['A', 'B', 'C', 'D', 'E', 'F']:
    v = res.bottom.F_in.get(c, 0.0)
    f_c = F_feed[c]
    rec_bot = 100 * v / f_c if f_c > 0 else float('nan')
    print(f"    {c}({names[c]:5s}): {v:10.4f} kmol/h ({100*v/F_bot_total:6.2f}%)   "
          f"to_bot={rec_bot:7.4f}%")
print(f"    TOTAL F_bot = {F_bot_total:.2f}")
print()
print(f"  ★ C3H6 漏れ率 = {100*res.top.F_in.get('B',0)/F_feed['B']:.4f}%   (HYSYS 報告: ~1%)")
print(f"  ★ C3H8 漏れ率 = {100*res.top.F_in.get('A',0)/F_feed['A']:.4f}%")
print(f"  ★ C2H6 塔頂回収率 = {100*res.top.F_in.get('F',0)/F_feed['F']:.4f}%   (LK、99% 行きたい)")
print()

# (2) Wang-Henke 直呼び (内部状態詳細)
print("=" * 100)
print("[2] wang_henke_solve 直呼び (詳細トレース)")
print("=" * 100)

# まず FUG で D_total 推定
tun_fug = ColumnTunables(
    P_col=8.5e5, N_stages=60, N_feed=1, reflux_ratio=10.0,
    solver_method='fug',
)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    res_fug = simulate_column2(feed, tunables=tun_fug)
D_fug = sum(res_fug.top.F_in.values())
print(f"  FUG D_total = {D_fug:.2f} kmol/h")
print(f"  FUG F_top: {[(c, round(res_fug.top.F_in[c],2)) for c in 'ABCDEF']}")
print(f"  FUG T_top = {res_fug.equipment.T_top-273.15:.2f}°C")
print()

# rigorous 直呼び
rig = wang_henke_solve(
    feed_F=F_feed, comps=['A','B','C','D','E','F'],
    P_col=8.5e5, N_stages=60,
    N_feed=max(1, min(res_fug.equipment.N_feed_kirkbride, 60)),
    reflux_ratio=10.0, D_total=D_fug,
    q_feed=0.0, partial_condenser=True,
    T_top_init_K=res_fug.equipment.T_top,
    T_bot_init_K=res_fug.equipment.T_bot,
    max_iter=1000,
)

print(f"  rigorous: converged={rig.converged}, iter={rig.n_iter}")
print(f"  msg: {rig.message}")
print(f"  mesh_residual_max = {rig.mesh_residual_max:.3e}")
print(f"  component_balance_max = {rig.component_balance_max:.3e}")
print()
print(f"  --- T プロファイル (10 段おき + 塔頂/塔底/フィード段) ---")
for j in [1, 2, 5, 10, 15, 20, res_fug.equipment.N_feed_kirkbride, 30, 40, 50, 55, 58, 59, 60]:
    if 1 <= j <= 60:
        T_K = rig.T_profile_K[j-1]
        print(f"    stage {j:3d}: T = {T_K - 273.15:7.2f} °C")
print()

# stage 1 と stage N の組成
print(f"  --- stage 1 (塔頂 partial cond) ---")
print(f"  {'c':>3} | {'x_1':>10} | {'y_1':>10} | {'K_1':>10}")
for c in ['A','B','C','D','E','F']:
    x1 = rig.x_profile[0].get(c, 0)
    y1 = rig.y_profile[0].get(c, 0)
    K1 = rig.K_profile[0].get(c, 0)
    print(f"  {c}({names[c]:4s}) | {x1:10.5e} | {y1:10.5e} | {K1:10.4f}")
print()
print(f"  --- stage 60 (塔底 reboiler) ---")
print(f"  {'c':>3} | {'x_N':>10} | {'y_N':>10} | {'K_N':>10}")
for c in ['A','B','C','D','E','F']:
    xn = rig.x_profile[-1].get(c, 0)
    yn = rig.y_profile[-1].get(c, 0)
    Kn = rig.K_profile[-1].get(c, 0)
    print(f"  {c}({names[c]:4s}) | {xn:10.5e} | {yn:10.5e} | {Kn:10.4f}")
print()

# 中間段で B (C3H6) y がどう変わるか
print(f"  --- B (C3H6) y 値の塔内変化 (stage 1=塔頂、60=塔底) ---")
for j in range(1, 61, 3):
    yj_B = rig.y_profile[j-1].get('B', 0)
    xj_B = rig.x_profile[j-1].get('B', 0)
    Kj_B = rig.K_profile[j-1].get('B', 0)
    print(f"    stage {j:3d}: x(B)={xj_B:.4f}, y(B)={yj_B:.4f}, K(B)={Kj_B:.4f}")

# C3H6 漏れ
F_top_B = D_fug * rig.y_profile[0].get('B', 0)
print()
print(f"  ★ wang_henke 直呼び C3H6 漏れ = D × y_1(B) = {D_fug:.0f} × {rig.y_profile[0].get('B',0):.4f}")
print(f"     = {F_top_B:.2f} kmol/h = {100*F_top_B/F_feed['B']:.4f}% of feed")

# 比較 (simulate_column2 経由は内部で proxy 罰則の閾値丸めをかける)
print()
print("=" * 100)
print("  まとめ")
print("=" * 100)
print(f"  HYSYS 報告 (N=60, P=8.5bar, R=10): C3H6 漏れ ≈ 1% (99% 分離)")
print(f"  当シミュ simulate_column2: C3H6 漏れ = {100*res.top.F_in.get('B',0)/F_feed['B']:.4f}%")
print(f"  当シミュ wang_henke 直呼び: C3H6 漏れ = {100*F_top_B/F_feed['B']:.4f}%")
