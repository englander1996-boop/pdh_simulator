"""Dist2 (partial condenser) のコンデンサ熱量を「簡易式」vs「厳密エネルギー収支」で比較する診断。

目的 (2026-05-26): main(自作 rigorous) の Dist2 コンデンサ Q_cond が HYSYS(62MW) より大幅に
小さい(8MW)のは「FUG vs rigorous」ではなく、自作モデルが partial-cond の Q_cond を
  Q_cond = R × (塔頂の凝縮性成分) × λ   (潜熱だけ、非凝縮ガスの顕熱を無視)
という簡易式で計算しているため、という仮説の検証。

ここでは収束した rigorous プロファイル(V/L 流量・各段組成・温度)から凝縮器の厳密
エネルギー収支 Q_cond = V_2·Hv(y_2,T_2) − D·Hv(y_1,T_1) − L_1·Hl(x_1,T_1) を計算し、
現行の簡易式と比較する。main #282 の Dist2 フィード(収束時)を入力に使う。
"""
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

import flowsheet  # noqa: 循環 import 回避 (exp2.py と同じ)

from stream.stream import ProcessStream
from src.distillation_core import (
    DistDesignVars, DistFixedParams, simulate_distillation_column,
    _weighted_lambda, NON_CONDENSABLE_COMPS,
)
from src.distillation_rigorous import wang_henke_solve
from src.eos import z_factor, residual_enthalpy
from src.thermo import PDHThermo

_thermo = PDHThermo()
_T_REF = 298.15


def H_mol(T, P, comp_frac: dict, keys, phase: str) -> float:
    """混合物のモルエンタルピー [J/mol] = 理想気体(T_ref→T) + 残差(PR EOS)。"""
    s = sum(max(comp_frac.get(k, 0.0), 0.0) for k in keys)
    if s <= 0:
        return 0.0
    x = [max(comp_frac.get(k, 0.0), 0.0) / s for k in keys]
    H_ig = sum(x[i] * _thermo.calc_enthalpy_change(keys[i], _T_REF, T) for i in range(len(keys)))
    try:
        Z = z_factor(T, P, x, keys, phase)
        H_r = residual_enthalpy(T, P, x, keys, Z)
    except Exception:
        H_r = 0.0
    return H_ig + H_r


def main():
    # main #282 収束時の Dist2 フィード (exp1_202605261419.txt の [Desuper → Dist2] より)
    feed = ProcessStream(
        F_in={'A': 3846.6, 'B': 2277.8, 'C': 882.7, 'D': 14.6, 'E': 268.8, 'F': 254.2, 'Z': 0.0},
        T_in=323.15, P_in=558426.0,
    )
    P_col = 558426.0
    R = 8.0627
    N = 38
    keys = ['A', 'B', 'C', 'D', 'E', 'F', 'Z']

    # Dist2 = partial condenser, LK=F(C2H6), HK=A(C3H8), q=0 (column2.py 準拠)
    design = DistDesignVars(
        P_col=P_col, N_stages=N, N_feed=1, reflux_ratio=R,
        LK='F', HK='A', recovery_LK_top=0.9675, recovery_HK_bot=0.9986,
        K_method='pr', q=0.0, partial_condenser=True, solver_method='rigorous',
    )

    # --- 現行実装の rigorous 結果 (簡易式 Q_cond) ---
    res = simulate_distillation_column(design, feed, DistFixedParams())
    if not res.equipment.feasible:
        print("rigorous infeasible:", res.equipment.message); return
    Q_simpl = res.equipment.Q_cond
    T_top = res.equipment.T_top

    # --- 厳密エネルギー収支のためにプロファイルを取得 (FUG→wang_henke を再現) ---
    from dataclasses import replace as dcr
    fug = simulate_distillation_column(dcr(design, solver_method='fug'), feed, DistFixedParams())
    D_total = sum(fug.top.F_in.values())
    N_feed_kirk = max(1, min(fug.equipment.N_feed_kirkbride, N))
    kw = dict(
        feed_F=feed.F_in, comps=keys, P_col=P_col, N_stages=N, N_feed=N_feed_kirk,
        reflux_ratio=R, D_total=D_total, q_feed=0.0, partial_condenser=True,
        K_method='pr', T_top_init_K=fug.equipment.T_top, T_bot_init_K=fug.equipment.T_bot,
    )
    rig = wang_henke_solve(**kw)
    if not rig.converged or rig.mesh_residual_max > 0.01 or rig.component_balance_max > 0.01:
        # _simulate_rigorous と同じ retry (Wegstein off + 低 damping)
        from src.distillation_rigorous import _RETRY_DAMPING
        rig2 = wang_henke_solve(**kw, max_iter=1000, damping=_RETRY_DAMPING, use_wegstein=False)
        if rig2.converged:
            rig = rig2
    if not rig.converged:
        print("wang_henke 未収束:", rig.message); return

    # 凝縮器(stage1=平衡段) のエネルギー収支
    #   V_2 (stage2 から上昇する蒸気) = V_top, comp=y_profile[1], T=T_profile[1]
    #   D (vapor distillate) = V_top - L_top, comp=y_profile[0], T=T_profile[0]
    #   L_1 (reflux liquid)  = L_top, comp=x_profile[0], T=T_profile[0]
    V2 = rig.V_top_kmolh
    L1 = rig.L_top_kmolh
    D = V2 - L1
    T1 = rig.T_profile_K[0]
    T2 = rig.T_profile_K[1] if len(rig.T_profile_K) > 1 else T1
    y2 = rig.y_profile[1] if len(rig.y_profile) > 1 else rig.y_profile[0]
    y1 = rig.y_profile[0]
    x1 = rig.x_profile[0]

    Hv2 = H_mol(T2, P_col, y2, keys, 'vapor')
    Hv1 = H_mol(T1, P_col, y1, keys, 'vapor')
    Hl1 = H_mol(T1, P_col, x1, keys, 'liquid')
    # kmol/h × J/mol = kJ/h → /3600 = kJ/s = kW
    Q_EB = (V2 * Hv2 - D * Hv1 - L1 * Hl1) / 3600.0

    # 参考: 全塔頂蒸気を T_top まで顕熱冷却する分 (非凝縮ガス込み) のオーダー
    print("=" * 70)
    print("Dist2 (partial cond) コンデンサ熱量: 簡易式 vs 厳密エネルギー収支")
    print("=" * 70)
    print(f"  rigorous T_top = {T_top-273.15:.1f} °C,  R = {R},  D_total = {D_total:.1f} kmol/h")
    print(f"  V_2(塔頂蒸気) = {V2:.0f},  L_1(還流) = {L1:.0f},  D(留出) = {D:.0f} kmol/h")
    print(f"  T_stage1 = {T1-273.15:.1f}°C,  T_stage2 = {T2-273.15:.1f}°C")
    print("-" * 70)
    print(f"  [現行] 簡易式 Q_cond = R×(凝縮性成分)×λ      = {Q_simpl:10.0f} kW = {Q_simpl/1000:6.1f} MW")
    print(f"  [厳密] エネルギー収支 Q_cond                  = {Q_EB:10.0f} kW = {Q_EB/1000:6.1f} MW")
    if Q_simpl > 0:
        print(f"  → 厳密 / 簡易 = {Q_EB/Q_simpl:.2f} 倍  (簡易式の過小評価率 = {(1-Q_simpl/Q_EB)*100:.0f}%)")
    print("=" * 70)
    print("  参考: HYSYS(special同等条件) は ~62 MW。本診断は main の Dist2 feed/設計での値。")


if __name__ == '__main__':
    main()
