"""
exp2 設計変数を振って全体収率の感度を見るための一時スクリプト。
ベースライン (T_in=950, z_cat=30, t_cyc=15, A_mem=1e5) と
複数バリエーションを連続実行して比較表を出す。
"""
import os
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

config = load_operating_config()

PSA_FIXED = PSADesignVars(D_col=3.0, L_bed=20.0, desorption_target=0.35)


def run(label, T_in, z_cat, t_cyc, A_mem, D=7.0, P_H=9.5e5):
    swing = SwingDesign(T_in=T_in, z_cat=z_cat, t_cyc=t_cyc, D=D)
    mem = MemDesignVars(P_H=P_H, P_L=1.0e5, A_mem=A_mem, P_dist=20.0e5)
    design = FlowsheetDesignVars(swing=swing, psa=PSA_FIXED, mem=mem)
    res = evaluate(design, config, verbose=False)
    if res.economics is None:
        print(f"{label}: NG ({res.failure_reason})")
        return None
    R = res.solver.one_pass
    F_fresh = res.solver.fresh_C3H8
    F_prod = R['r3'].top.F_in.get('B', 0.0)
    yld = F_prod / F_fresh * 100.0
    conv = R['r_rx'].performance.Conversion
    sel = R['r_rx'].performance.Selectivity
    tac = res.effective_TAC
    print(f"{label:30s} | T_in={T_in:.0f} z={z_cat:.0f} t={t_cyc:.0f} A={A_mem:.0e} "
          f"| 収率={yld:5.2f}% 転化={conv:5.1f}% 選択={sel:5.1f}% TAC={tac:6.1f}")
    return yld, conv, sel, tac


print("=" * 110)
print(f"{'ケース':30s} | 設計変数                              | 性能")
print("=" * 110)

# 前回までのベスト (参考のため 1 ケースだけ実行)
# --- z=15 周辺をさらに絞り込み ---
run("T=890, z=15, t=2, A=5e5 (ref)", T_in=890.0, z_cat=15.0, t_cyc=2.0, A_mem=5.0e5)
run("T=895, z=15, t=2, A=5e5", T_in=895.0, z_cat=15.0, t_cyc=2.0, A_mem=5.0e5)
run("T=900, z=15, t=2, A=5e5", T_in=900.0, z_cat=15.0, t_cyc=2.0, A_mem=5.0e5)
run("T=905, z=15, t=2, A=5e5", T_in=905.0, z_cat=15.0, t_cyc=2.0, A_mem=5.0e5)
run("T=890, z=12, t=2, A=5e5", T_in=890.0, z_cat=12.0, t_cyc=2.0, A_mem=5.0e5)
run("T=890, z=18, t=2, A=5e5", T_in=890.0, z_cat=18.0, t_cyc=2.0, A_mem=5.0e5)
# t_cyc 極限
run("T=890, z=15, t=1, A=5e5", T_in=890.0, z_cat=15.0, t_cyc=1.0, A_mem=5.0e5)
