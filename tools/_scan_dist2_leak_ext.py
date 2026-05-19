"""Dist2 拡張スキャン (低圧×多段×高還流)。
ユーザー (HYSYS で確認済) によると「段数↑+圧力↓」で漏れ <1% は実現可能。
当該領域が当シミュレータで再現できるかを確認する。
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings

import flowsheet  # noqa
from src.distillation_core import ColumnTunables
from stream.stream import ProcessStream
from units.separators.column2.column2 import simulate_column2

# trial #258 相当のフィード
F_feed = {
    'A': 4750.0, 'B': 2900.0, 'C': 1075.0,
    'D': 18.0,   'E': 220.0,  'F': 203.0,
}
F_C3H6 = F_feed['B']
F_total = sum(F_feed.values())

# 拡張スキャン軸
P_list = [2.0, 3.0, 4.0, 5.0, 7.0]      # bar (低圧側)
N_list = [25, 35, 50, 70]                # 多段側
R_list = [6.0, 8.0, 12.0, 16.0]          # 高還流側

print("=" * 110)
print(f"  Dist2 拡張スキャン: 低 P + 多 N + 高 R 領域 (rigorous)")
print(f"  feed: C3H6={F_C3H6} kmol/h, total={F_total:.1f}")
print("=" * 110)
print()

rows = []
total_cases = len(P_list) * len(N_list) * len(R_list)
done = 0

for P_bar in P_list:
    for N in N_list:
        for R in R_list:
            done += 1
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
                    rows.append((P_bar, N, R, None, None, None, None, None, "infeasible"))
                    continue
                B_top = res.top.F_in.get('B', 0.0)
                A_top = res.top.F_in.get('A', 0.0)
                leak_B = 100.0 * B_top / F_C3H6
                leak_A = 100.0 * A_top / F_feed['A']
                T_top_C = res.equipment.T_top - 273.15
                util = res.equipment.cond_utility_name
                N_min_v = res.equipment.N_min
                R_min_v = res.equipment.R_min
                rows.append((P_bar, N, R, leak_B, leak_A, T_top_C, util, (N_min_v, R_min_v), "OK"))
            except Exception as e:
                rows.append((P_bar, N, R, None, None, None, None, None, f"ERR: {type(e).__name__}"))
            # 経過表示 (stderr で flush)
            print(f"  [{done}/{total_cases}] P={P_bar} N={N} R={R} → {rows[-1][-1]}",
                  file=sys.stderr, flush=True)

print()
print(f"{'P[bar]':>7} {'N':>3} {'R':>5} | {'B漏%':>7} {'A漏%':>7} {'T_top':>7} | {'N_min':>5} {'R_min':>6} | {'冷媒':<24}")
print("-" * 110)
for r in rows:
    P, N, R, lb, la, t, u, fug, st = r
    if lb is not None:
        nm, rm = fug
        print(f"{P:7.1f} {N:3d} {R:5.1f} | {lb:7.4f} {la:7.4f} {t:7.1f} | {nm:5.2f} {rm:6.3f} | {u:<24}")
    else:
        print(f"{P:7.1f} {N:3d} {R:5.1f} | {st:>7}")

print()
print("=" * 80)
print("  B (C3H6) 漏れが小さい設計 top 15")
print("=" * 80)
ok = [r for r in rows if r[3] is not None]
ok.sort(key=lambda x: x[3])
for r in ok[:15]:
    P, N, R, lb, la, t, u, fug, st = r
    nm, rm = fug
    print(f"  P={P:5.1f}bar N={N:3d} R={R:5.1f} (R/Rmin={R/rm:.2f}, N/Nmin={N/nm:.1f})"
          f" → B漏れ {lb:.4f}%, A漏れ {la:.4f}%, T_top {t:.1f}°C, 冷媒={u}")
