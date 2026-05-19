"""Dist2 (partial cond) の C3 漏れが物理的にどこまで小さくできるか P/N/R 空間でスキャン。
exp1 trial #258 相当のフィードで rigorous を回し、C3H6 top 漏れ % を表で示す。

目的: 「そもそも漏れない設計」が探索空間内に存在するかをユーザに見せて、
search_space.py / 設計判断を物理ベースで決めるための材料にする。
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings

import flowsheet  # noqa: 循環 import 回避
from src.distillation_core import ColumnTunables
from stream.stream import ProcessStream
from units.separators.column2.column2 import simulate_column2

# trial #258 相当のフィード (Dist2 入口)
F_feed = {
    'A': 4750.0,   # C3H8
    'B': 2900.0,   # C3H6
    'C': 1075.0,   # H2
    'D': 18.0,     # C2H4
    'E': 220.0,    # CH4
    'F': 203.0,    # C2H6
}
F_C3H6 = F_feed['B']

# スキャン (P, N, R) — 探索範囲は SEARCH_SPACE と合わせる
P_list = [5.0, 7.5, 10.0, 12.5, 15.0, 18.0, 22.0]   # bar
N_list = [20, 26, 32, 40]
R_list = [4.0, 6.0, 8.0, 10.0]

print("=" * 110)
print(f"  Dist2 partial cond: C3H6 top 漏れ % スキャン (rigorous)")
print(f"  feed: C3H6={F_C3H6} kmol/h, total={sum(F_feed.values()):.1f}")
print(f"  目的: 物理的に漏れ ≪ 1% を実現できる (P, N, R) があるか")
print("=" * 110)
print()

rows = []
for P_bar in P_list:
    for N in N_list:
        for R in R_list:
            feed = ProcessStream(F_in=F_feed, T_in=50.0 + 273.15, P_in=P_bar * 1e5)
            tun = ColumnTunables(
                P_col=P_bar * 1e5, N_stages=N,
                N_feed=1, reflux_ratio=R,
                solver_method='rigorous',
            )
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    res = simulate_column2(feed, tunables=tun)
                if not res.equipment.feasible:
                    rows.append((P_bar, N, R, None, None, None, None, "infeasible"))
                    continue
                B_top = res.top.F_in.get('B', 0.0)
                A_top = res.top.F_in.get('A', 0.0)
                leak_B = 100.0 * B_top / F_C3H6
                leak_A = 100.0 * A_top / F_feed['A']
                T_top_C = res.equipment.T_top - 273.15
                util = res.equipment.cond_utility_name
                rows.append((P_bar, N, R, leak_B, leak_A, T_top_C, util, "OK"))
            except Exception as e:
                rows.append((P_bar, N, R, None, None, None, None, f"ERR: {type(e).__name__}"))

# 出力
print(f"{'P[bar]':>7} {'N':>3} {'R':>5} | {'B漏%':>6} {'A漏%':>6} {'T_top':>7} | {'冷媒':<24} | status")
print("-" * 110)
for r in rows:
    P, N, R, lb, la, t, u, st = r
    if lb is not None:
        print(f"{P:7.1f} {N:3d} {R:5.1f} | {lb:6.3f} {la:6.3f} {t:7.1f} | {u:<24} | {st}")
    else:
        print(f"{P:7.1f} {N:3d} {R:5.1f} | {'-':>6} {'-':>6} {'-':>7} | {'-':<24} | {st}")

# ベスト 10 (B漏れ最小)
print()
print("=" * 80)
print("  B (C3H6) 漏れが小さい設計 top 10")
print("=" * 80)
ok = [r for r in rows if r[3] is not None]
ok.sort(key=lambda x: x[3])
for r in ok[:10]:
    P, N, R, lb, la, t, u, st = r
    print(f"  P={P:5.1f}bar N={N:3d} R={R:5.1f} → B漏れ {lb:.3f}%, "
          f"A漏れ {la:.3f}%, T_top {t:.1f}°C, 冷媒={u}")
