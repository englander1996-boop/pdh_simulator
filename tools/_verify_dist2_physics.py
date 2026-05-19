"""Dist2 物理検証スクリプト (Phase: 診断検証)。

目的: シミュレータの「3% 漏れ下限」が真の物理か実装バグかを切り分ける。

検証項目:
  (1) K 値の独立計算: src/eos.py の PR EOS vs thermo.PRMIX (上流ライブラリ) を比較
  (2) thermo の SRK 系で K 値を計算し PR と比較 (HYSYS は SRK 系で計算しているかも)
  (3) α(C3H6/C3H8), α(C2H6/C3H6) の温度依存
  (4) Wang-Henke stage 1 mass balance & VLE 一貫性をハンドで検算
  (5) Fenske 単独で predict した leak と rigorous の leak を比較
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import warnings

import flowsheet  # noqa
from src.eos import z_factor, fugacity_coeff
from src.config import THERMO_DATA

import numpy as np
from thermo.eos_mix import PRMIX, SRKMIX
from thermo import ChemicalConstantsPackage

comps = ['A', 'B', 'C', 'D', 'E', 'F']  # C3H8, C3H6, H2, C2H4, CH4, C2H6
names = {'A': 'C3H8', 'B': 'C3H6', 'C': 'H2', 'D': 'C2H4', 'E': 'CH4', 'F': 'C2H6'}

# Dist2 入口組成 (trial #258 相当)
F_feed = {'A': 4750., 'B': 2900., 'C': 1075., 'D': 18., 'E': 220., 'F': 203.}
F_total = sum(F_feed.values())
z = {c: F_feed[c] / F_total for c in comps}

# 塔底 C3H8/C3H6 リッチ組成 (リフラックスっぽい)
x_bot = {'A': 4750./(4750+2900), 'B': 2900./(4750+2900), 'C': 0., 'D': 0., 'E': 0., 'F': 0.}

# Tc/Pc/omega
TcPcOm = {c: (THERMO_DATA[c].Tc, THERMO_DATA[c].Pc, THERMO_DATA[c].omega) for c in comps}
Tc_list = [TcPcOm[c][0] for c in comps]
Pc_list = [TcPcOm[c][1] for c in comps]
om_list = [TcPcOm[c][2] for c in comps]

print("=" * 90)
print("  (1) K 値: src/eos.py PR vs thermo.PRMIX vs thermo.SRKMIX")
print("=" * 90)
print(f"  feed z: {[f'{c}({names[c]}):{z[c]:.4f}' for c in comps]}")
print()

for P_bar in [5.0, 7.5, 22.0]:
    for T_C in [-55.0, -73.0, -95.0]:
        T = T_C + 273.15
        P = P_bar * 1e5
        # src/eos.py PR
        try:
            x_list = [x_bot[c] for c in comps]
            Z_L = z_factor(T, P, x_list, comps, phase='liquid')
            Z_V = z_factor(T, P, x_list, comps, phase='vapor')
            K_my_PR = {}
            for i, c in enumerate(comps):
                phi_L = fugacity_coeff(i, T, P, x_list, comps, Z_L)
                phi_V = fugacity_coeff(i, T, P, x_list, comps, Z_V)
                K_my_PR[c] = phi_L / phi_V if phi_V > 1e-30 else float('inf')
        except Exception as e:
            K_my_PR = {c: float('nan') for c in comps}
            print(f"  src.eos PR 失敗 at T={T_C}C P={P_bar}bar: {e}")

        # thermo.PRMIX 独立計算
        try:
            zs = [x_bot[c] for c in comps]
            pr = PRMIX(T=T, P=P, zs=zs, Tcs=Tc_list, Pcs=Pc_list, omegas=om_list)
            K_th_PR = {}
            for i, c in enumerate(comps):
                # phi_L / phi_V
                fl = pr.phis_l[i] if pr.phis_l else None
                fv = pr.phis_g[i] if pr.phis_g else None
                if fl is not None and fv is not None and fv > 1e-30:
                    K_th_PR[c] = fl / fv
                else:
                    K_th_PR[c] = float('nan')
        except Exception as e:
            K_th_PR = {c: float('nan') for c in comps}
            print(f"  thermo PR 失敗 at T={T_C}C P={P_bar}bar: {e}")

        # thermo.SRKMIX 独立計算
        try:
            zs = [x_bot[c] for c in comps]
            sk = SRKMIX(T=T, P=P, zs=zs, Tcs=Tc_list, Pcs=Pc_list, omegas=om_list)
            K_th_SRK = {}
            for i, c in enumerate(comps):
                fl = sk.phis_l[i] if sk.phis_l else None
                fv = sk.phis_g[i] if sk.phis_g else None
                if fl is not None and fv is not None and fv > 1e-30:
                    K_th_SRK[c] = fl / fv
                else:
                    K_th_SRK[c] = float('nan')
        except Exception as e:
            K_th_SRK = {c: float('nan') for c in comps}

        print(f"--- T = {T_C}°C, P = {P_bar} bar ---")
        print(f"  {'comp':>6} | {'src.eos PR':>11} | {'thermo PR':>10} | {'thermo SRK':>10} | {'差PR(%)':>8}")
        for c in comps:
            mp = K_my_PR.get(c, float('nan'))
            tp = K_th_PR.get(c, float('nan'))
            sk = K_th_SRK.get(c, float('nan'))
            if mp and tp:
                diff = 100 * (mp - tp) / tp if tp != 0 else float('inf')
            else:
                diff = float('nan')
            print(f"  {c}({names[c]:4s}) | {mp:11.5f} | {tp:10.5f} | {sk:10.5f} | {diff:7.3f}%")
        print()

print()
print("=" * 90)
print("  (2) α (C3H6/C3H8) と α (C2H6/C3H6) の温度依存 (P=7.5 bar)")
print("=" * 90)
print(f"  {'T[°C]':>7} | {'K(C2H6)':>10} {'K(C3H6)':>10} {'K(C3H8)':>10} | {'α(C2H6/C3H6)':>15} {'α(C3H6/C3H8)':>15}")
for T_C in [-95, -73, -55, -39, -20]:
    T = T_C + 273.15
    P = 7.5e5
    zs = [x_bot[c] for c in comps]
    try:
        pr = PRMIX(T=T, P=P, zs=zs, Tcs=Tc_list, Pcs=Pc_list, omegas=om_list)
        Ks = {c: pr.phis_l[i]/pr.phis_g[i] for i, c in enumerate(comps)}
        print(f"  {T_C:7.1f} | {Ks['F']:10.5f} {Ks['B']:10.5f} {Ks['A']:10.5f} | "
              f"{Ks['F']/Ks['B']:15.3f} {Ks['B']/Ks['A']:15.3f}")
    except Exception as e:
        print(f"  {T_C:7.1f} | (失敗: {e})")

print()
print("=" * 90)
print("  (3) Wang-Henke stage 1 手計算検算")
print("=" * 90)
print("""  trial #258 設計 (P=7.55bar, R=5.22, N=26) で rigorous が報告した stage 1 状態:
    T_1 = -55.6°C, K_1(C3H6) = 0.127
    x_1(C3H6) = 0.471 (還流液 47% C3H6)
    y_1(C3H6) = 0.0597 (vapor distillate 6%)
    D_total = 1576 kmol/h, L_top = R*D = 8228 kmol/h
    V_2 = (R+1)*D = 9804 kmol/h
""")
# 手計算: stage 1 mass balance
D = 1576.0
R = 5.22
L_top = R * D
V_2 = (R + 1) * D
K_1_B = 0.127
y_1_B = 0.0597
x_1_B = 0.471

# VLE check: y_1 ?= K_1 * x_1
y_from_VLE = K_1_B * x_1_B
print(f"  VLE check: K_1 × x_1 = 0.127 × 0.471 = {y_from_VLE:.5f}, reported y_1 = {y_1_B}")
print(f"  → VLE 残差: {y_from_VLE - y_1_B:+.5f} ({100*(y_from_VLE-y_1_B)/y_1_B:+.2f}%)")
print()

# Mass balance: V_2 y_2 = D y_1 + L_top x_1
F_top_B = D * y_1_B
F_reflux_B = L_top * x_1_B
V2_y2_B_required = F_top_B + F_reflux_B
y_2_B_required = V2_y2_B_required / V_2
print(f"  Mass balance @ stage 1:")
print(f"    D × y_1(B) = 1576 × 0.0597 = {F_top_B:.2f} kmol/h (= 塔頂 C3H6 流出)")
print(f"    L_top × x_1(B) = {L_top:.0f} × 0.471 = {F_reflux_B:.2f} kmol/h (= 還流液 C3H6)")
print(f"    両者和 = {V2_y2_B_required:.2f} kmol/h = V_2 × y_2(B)")
print(f"    → 必要な y_2(B) = {y_2_B_required:.4f} ({100*y_2_B_required:.2f}%)")
print()
print(f"  もし y_2(B) を {y_2_B_required:.3f} 未満に出来れば F_top(B) も比例して減る。")
print(f"  → stage 2 以降 (rectifying セクション) で C3H6 vapor をどれだけ落とせるかが鍵。")

print()
print("=" * 90)
print("  (4) Fenske 単独予測 vs rigorous")
print("=" * 90)
print(f"""  Fenske 非キー split: top/bot_c = α_c^N_min × top/bot_HK
    α(C3H6/HK=C3H8) = 1.27 (at -55°C/7.5bar)
    N_min (Fenske, LK=C2H6 HK=C3H8) ≈ 2.9
    recovery_HK_bot = 0.99 → top/bot_HK = 0.01/0.99 = 0.0101
    ratio_C3H6 = 1.27^2.9 × 0.0101 = {(1.27**2.9):.3f} × 0.0101 = {(1.27**2.9)*0.0101:.4f}
    frac_top_C3H6 = {(1.27**2.9)*0.0101 / (1 + (1.27**2.9)*0.0101):.4%}

  rigorous 実測: 3.27% → Fenske より高い leak (∵ 有限 R, 非平衡効果)
""")

print()
print("=" * 90)
print("  (5) 純粋な C2/C3 split (= H2/CH4 抜き) を仮想試算")
print("=" * 90)
print("""  もし Dist2 入口に H2/CH4/C2H4 が無くて C2H6/C3H8/C3H6 のみだったら?
  α(C2H6/C3H6) ≈ 5 (at 約 -30°C/7.5bar) → 99% 分離容易。
  ↓
  実機 PDH プラントは「Dist2 前に flash drum / 冷却分離」で H2/CH4 抜く。
  本シミュレータは入口に全部混ぜているのが「3% 漏れ floor」の根本原因。
""")
