import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

# ==========================================
# 1. 計算関数 (プラグフローモデル)
# ==========================================
def calculate_cross_flow_membrane(theta, z, y_target, alpha, Q_A, P_H, F_feed):
    """プラグフロー(クロスフロー)モデルによる膜分離の数値計算"""
    Q_B = Q_A / alpha
    F_in = [F_feed * z, F_feed * (1 - z)]
    
    def simulate_membrane(gamma):
        P_L = gamma * P_H
        
        def ode_sys(A, F):
            F_C3H6 = max(F[0], 1e-12)
            F_C3H8 = max(F[1], 1e-12)
            F_tot = F_C3H6 + F_C3H8
            x = F_C3H6 / F_tot
            
            a = (1 - alpha) * gamma
            b = (alpha - 1) * (x + gamma) + 1
            c = -alpha * x
            
            discriminant = max(0, b**2 - 4*a*c)
            y_local = (2 * c) / (-b - np.sqrt(discriminant))
            
            J_C3H6 = Q_A * (x * P_H - y_local * P_L)
            J_C3H8 = Q_B * ((1 - x) * P_H - (1 - y_local) * P_L)
            
            if J_C3H6 < 0 or J_C3H8 < 0:
                return [0, 0]
            return [-J_C3H6, -J_C3H8]
            
        def event_cut(A, F):
            return (F[0] + F[1]) - F_feed * (1 - theta)
        event_cut.terminal = True
        
        # 大流量対応 (積分上限を10億に拡張)
        sol = solve_ivp(ode_sys, [0, 1e9], F_in, events=event_cut, method='Radau')
        
        if sol.status != 1:
            return 0.0, np.nan
            
        F_out = sol.y[:, -1]
        A_req = sol.t[-1]
        
        F_perm_C3H6 = F_in[0] - F_out[0]
        F_perm_tot = F_feed - (F_out[0] + F_out[1])
        y_avg = F_perm_C3H6 / F_perm_tot if F_perm_tot > 0 else 0
        
        return y_avg, A_req

    def obj_func(gamma):
        y_avg, _ = simulate_membrane(gamma)
        return y_avg - y_target

    try:
        if obj_func(0.001) < 0:
            return np.nan, np.nan, np.nan, np.nan
            
        gamma_opt = brentq(obj_func, 0.001, 0.99)
        _, A_req = simulate_membrane(gamma_opt)
        
        x_out = (z - theta * y_target) / (1 - theta)
        
        AREA_PER_MODULE = 500.0
        modules = np.ceil(A_req / AREA_PER_MODULE)
        
        return x_out, gamma_opt, A_req, modules
        
    except ValueError:
        return np.nan, np.nan, np.nan, np.nan

# ==========================================
# 2. 模式図の描画関数
# ==========================================
def draw_schematic(F_feed, z_feed, P_H, theta, x_out, y_target, gamma, Area, Modules):
    """膜分離プロセスの物質収支と設備要件を模式図として描画する"""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis('off')  # 枠線を消す
    
    # 流量と圧力の計算
    F_perm = F_feed * theta
    F_ret = F_feed * (1 - theta)
    P_L = P_H * gamma
    
    # 膜モジュール（中央の四角）
    rect = patches.Rectangle((0.3, 0.4), 0.4, 0.4, linewidth=2, edgecolor='navy', facecolor='aliceblue', zorder=2)
    ax.add_patch(rect)
    ax.text(0.5, 0.7, 'Membrane Module\n(ZIF-8 Cross-Flow)', ha='center', va='center', fontsize=12, fontweight='bold', color='navy')
    ax.text(0.5, 0.5, f"Required Area: {Area:,.0f} m$^2$\nModules (8-inch): {int(Modules):,} units", ha='center', va='center', fontsize=11, color='darkred', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    # 膜の点線
    ax.plot([0.35, 0.65], [0.6, 0.6], 'k--', lw=2, zorder=3)

    # 供給ガス (Feed) の矢印とテキスト
    ax.arrow(0.05, 0.65, 0.2, 0, head_width=0.03, head_length=0.03, fc='black', ec='black')
    ax.text(0.15, 0.68, 'Feed Gas\n'
                        f'Flow: {F_feed:,.0f} kmol/h\n'
                        f'C3H6: {z_feed*100:.1f} %\n'
                        f'Press: {P_H:.1f} bar', 
            ha='center', va='bottom', fontsize=10)

    # 非透過ガス (Retentate) の矢印とテキスト
    ax.arrow(0.7, 0.65, 0.2, 0, head_width=0.03, head_length=0.03, fc='black', ec='black')
    ax.text(0.82, 0.68, 'Retentate (Non-Permeate)\n'
                        f'Flow: {F_ret:,.1f} kmol/h\n'
                        f'C3H6: {x_out*100:.1f} %\n'
                        f'Press: ~{P_H:.1f} bar', 
            ha='center', va='bottom', fontsize=10)

    # 透過ガス (Permeate) の矢印とテキスト
    ax.arrow(0.5, 0.4, 0, -0.2, head_width=0.03, head_length=0.03, fc='darkgreen', ec='darkgreen')
    ax.text(0.53, 0.25, 'Permeate (Product)\n'
                        f'Flow: {F_perm:,.1f} kmol/h\n'
                        f'C3H6: {y_target*100:.1f} %\n'
                        f'Press: {P_L:.2f} bar\n'
                        f'(Ratio $\gamma$: {gamma:.3f})', 
            ha='left', va='center', fontsize=10, color='darkgreen')

    plt.title(f'Process Flow Schematic (Stage Cut $\\theta$ = {theta:.2f})', fontsize=14, fontweight='bold')
    fig.tight_layout()
    plt.show()

# ==========================================
# 3. 実行部 (メイン処理)
# ==========================================
z_feed = 0.50          # 供給ガス プロピレン分率
y_target = 0.90        # 透過側 目標プロピレン純度
alpha_val = 85.0       # 膜の選択性
Q_A_GPU = 10.0         # プロピレン透過度 [GPU]
P_H_bar = 5.0          # 供給側圧力 [bar]
F_feed_kmol_h = 3000.0 # 供給ガス流量 [kmol/h]

Q_A_SI = Q_A_GPU * 3.35e-10
P_H_SI = P_H_bar * 1e5
F_feed_SI = F_feed_kmol_h * 1000 / 3600

# 計算の実行
theta_array = np.linspace(0.05, 0.40, 15)
results = []

for th in theta_array:
    x, gm, area, mod = calculate_cross_flow_membrane(th, z_feed, y_target, alpha_val, Q_A_SI, P_H_SI, F_feed_SI)
    results.append({
        'カット率 (θ)': th,
        '必要圧力比 (γ)': gm,
        '非透過側 C3H6純度': x,
        '必要膜面積 [m2]': area,
        '管(モジュール)本数': mod
    })

df = pd.DataFrame(results)

# --- (A) 一覧表のコンソール出力 ---
print("\n" + "="*70)
print(f"【計算条件】供給流量: {F_feed_kmol_h:,.0f} kmol/h, 目標純度: {y_target*100}%, 選択性: {alpha_val}, 供給圧: {P_H_bar} bar")
print("="*70)
df_display = df.copy()
df_display['必要圧力比 (γ)'] = df_display['必要圧力比 (γ)'].map(lambda val: f"{val:.4f}" if pd.notnull(val) else "不可")
df_display['非透過側 C3H6純度'] = df_display['非透過側 C3H6純度'].map(lambda val: f"{val*100:.1f}%" if pd.notnull(val) else "-")
df_display['必要膜面積 [m2]'] = df_display['必要膜面積 [m2]'].map(lambda val: f"{val:,.0f}" if pd.notnull(val) else "-")
df_display['管(モジュール)本数'] = df_display['管(モジュール)本数'].map(lambda val: f"{int(val):,} 本" if pd.notnull(val) else "-")
print(df_display.to_string(index=False))
print("="*70 + "\n")

# --- (B) 本数付きトレンドグラフの描画 ---
df_valid = df.dropna()

fig, ax1 = plt.subplots(figsize=(11, 6))
color1 = 'tab:blue'
ax1.set_xlabel('Stage Cut (θ)', fontsize=12)
ax1.set_ylabel('Required Pressure Ratio (P_L / P_H)', color=color1, fontsize=12)
ax1.plot(df_valid['カット率 (θ)'], df_valid['必要圧力比 (γ)'], marker='o', color=color1, linewidth=2, label='Pressure Ratio')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, linestyle='--', alpha=0.7)

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel('Required Membrane Area [m²]', color=color2, fontsize=12)
ax2.plot(df_valid['カット率 (θ)'], df_valid['必要膜面積 [m2]'], marker='s', color=color2, linestyle='--', linewidth=2, label='Area')
ax2.tick_params(axis='y', labelcolor=color2)

for i, row in df_valid.iterrows():
    ax2.annotate(f"{int(row['管(モジュール)本数']):,}本", 
                 (row['カット率 (θ)'], row['必要膜面積 [m2]']), 
                 textcoords="offset points", 
                 xytext=(0, 12), 
                 ha='center', fontsize=10, color='darkred', fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

plt.title('Cross-Flow Membrane Sizing & Module Count', fontsize=14, fontweight='bold')
fig.tight_layout()
plt.show()

# --- (C) 代表的な設計点の模式図を描画 ---
design_theta = 0.20  # ここで模式図を描画したいカット率を指定します
x_des, gm_des, area_des, mod_des = calculate_cross_flow_membrane(design_theta, z_feed, y_target, alpha_val, Q_A_SI, P_H_SI, F_feed_SI)

if not np.isnan(gm_des):
    draw_schematic(F_feed_kmol_h, z_feed, P_H_bar, design_theta, x_des, y_target, gm_des, area_des, mod_des)
else:
    print(f"カット率 {design_theta} での分離は物理的に不可能なため、模式図を描画できません。")


