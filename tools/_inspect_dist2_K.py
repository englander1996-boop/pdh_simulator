"""Dist2 塔頂条件で PR EOS から K_i を直接計算する一時診断。
回収率99%が物理的に可能か確認するためのスクリプト。"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.eos import z_factor, fugacity_coeff, bubble_point_T

comps = ['A', 'B', 'C', 'D', 'E', 'F']  # C3H8, C3H6, H2, C2H4, CH4, C2H6

# Dist2 塔頂 vapor 組成 (exp1 から: H2/CH4/C2H4/C2H6 のみ)
F_top = {'A': 0.0, 'B': 0.0, 'C': 1014.4, 'D': 16.8, 'E': 209.3, 'F': 190.5}
F_top_tot = sum(F_top.values())
y_top = [F_top[c] / F_top_tot for c in comps]

# Dist2 塔底 liquid 組成 (exp1 から: C3H8 + C3H6 主、C2H6 微量)
F_bot = {'A': 4483.2, 'B': 2753.0, 'C': 0.0, 'D': 0.0, 'E': 0.0, 'F': 1.9}
F_bot_tot = sum(F_bot.values())
x_bot = [F_bot[c] / F_bot_tot for c in comps]

# Feed 組成 (Desuper → Dist2)
F_feed = {'A': 4483.2, 'B': 2753.0, 'C': 1014.4, 'D': 16.8, 'E': 209.3, 'F': 192.4}
F_feed_tot = sum(F_feed.values())
z_feed = [F_feed[c] / F_feed_tot for c in comps]

P_col = 7.5e5

print("=" * 70)
print(f"  Dist2 物理点検 (P = {P_col/1e5:.1f} bar)")
print("=" * 70)
print(f"feed   z = {[f'{c}:{z:.4f}' for c, z in zip(comps, z_feed)]}")
print(f"top    y = {[f'{c}:{y:.4f}' for c, y in zip(comps, y_top)]}")
print(f"bot    x = {[f'{c}:{x:.4f}' for c, x in zip(comps, x_bot)]}")
print()

# 各温度で K_i を計算 (top stage 想定: vapor 組成=y_top, liquid 組成=未知の reflux x_1)
# ここでは便宜上、両相とも y_top を使った場合と x_bot を使った場合を比較。
for T_C in [-43.2, -60.0, -30.0, -20.0, 13.0]:
    T = T_C + 273.15
    print(f"--- T = {T_C}°C ({T:.2f} K) ---")
    # 液相を x_bot (塔底組成)、気相を y_top (塔頂組成) で評価
    try:
        Z_L = z_factor(T, P_col, x_bot, comps, phase='liquid')
        Z_V = z_factor(T, P_col, y_top, comps, phase='vapor')
        print(f"  Z_L(x_bot) = {Z_L:.5f},  Z_V(y_top) = {Z_V:.5f}")
        print(f"  {'comp':>5} | {'phi_L(x_bot)':>14} | {'phi_V(y_top)':>14} | {'K_i':>10}")
        for i, c in enumerate(comps):
            try:
                phi_L = fugacity_coeff(i, T, P_col, x_bot, comps, Z_L)
                phi_V = fugacity_coeff(i, T, P_col, y_top, comps, Z_V)
                K_i = phi_L / phi_V if phi_V > 1e-30 else float('inf')
                print(f"  {c:>5} | {phi_L:14.5f} | {phi_V:14.5f} | {K_i:10.5f}")
            except Exception as e:
                print(f"  {c}: error {e}")
    except Exception as e:
        print(f"  z_factor 失敗: {e}")
    print()

# 泡点温度参照
print("=" * 70)
print("  bubble_point_T 参照")
print("=" * 70)
for name, x_use in [("x_bot (塔底液)", x_bot), ("z_feed", z_feed)]:
    try:
        T_bp = bubble_point_T(P_col, x_use, comps)
        print(f"  T_bp({name}, 7.5 bar) = {T_bp - 273.15:.2f} °C")
    except Exception as e:
        print(f"  T_bp({name}): 失敗 {e}")

# 別圧力での bubble_point_T も比較 (HYSYS との対照用)
print()
print("=" * 70)
print("  圧力依存 (塔底液で泡点 / 塔頂で露点)")
print("=" * 70)
from src.eos import dew_point_T  # noqa
for P_bar in [7.5, 15.0, 22.0, 25.0]:
    P_use = P_bar * 1e5
    try:
        Tbp_b = bubble_point_T(P_use, x_bot, comps)
    except Exception as e:
        Tbp_b = float('nan')
    try:
        Tdp_t = dew_point_T(P_use, y_top, comps)
    except Exception as e:
        Tdp_t = float('nan')
    print(f"  P = {P_bar:5.1f} bar : T_bot(bp) = {Tbp_b - 273.15:7.2f} °C, "
          f"T_top(dp) = {Tdp_t - 273.15:7.2f} °C")
