"""show_stage2_synthesis (HEN 組み合わせ・移動熱量表示) のレンダリング検証。

実フローシートを回さず、ダミー HENResult で表示関数だけを叩いて、新表が
正しく描画される (= 完了後のレポート再生成が動く) ことを確認する。
"""
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
import flowsheet  # noqa 循環import回避
from types import SimpleNamespace
from optimization.hen_synthesis import HENResult, HEMatch
from simulation.display import show_stage2_synthesis

m1 = HEMatch(name='HE_above_1', hot_label='H1_rx_cooler', cold_label='C2_rx_preheat',
             Q_kW=12300, T_h_in_K=797.15, T_h_out_K=480.15, T_c_in_K=320.15, T_c_out_K=760.15,
             LMTD_K=100, U_W_m2K=150, A_m2=820, CAPEX_okuyen=1.234, side='above')
m2 = HEMatch(name='HE_below_1', hot_label='H_Dist2_cond', cold_label='C_Dist3_reb',
             Q_kW=8000, T_h_in_K=355.15, T_h_out_K=355.15, T_c_in_K=318.15, T_c_out_K=318.15,
             LMTD_K=20, U_W_m2K=1000, A_m2=350, CAPEX_okuyen=0.512, side='below')
m3 = HEMatch(name='HE_below_2', hot_label='H_mem_cooler_gas', cold_label='C1_mem_precool_sens',
             Q_kW=4500, T_h_in_K=330.15, T_h_out_K=315.15, T_c_in_K=300.15, T_c_out_K=320.15,
             LMTD_K=12, U_W_m2K=200, A_m2=1875, CAPEX_okuyen=0.98, side='below')
hen = HENResult(matches=[m1, m2, m3], n_process_HE=3, n_utility_HE=0,
                Q_recovered_kW=24800, Q_hot_utility_kW=120000, Q_cold_utility_kW=70000,
                CAPEX_added_okuyen=2.726, OPEX_utility_okuyen={}, feasible=True)
res = SimpleNamespace(hen_result=hen, economics_synth=object())
show_stage2_synthesis(res)

print("\n--- infeasible ケース ---")
hen2 = HENResult(matches=[m1], n_process_HE=1, n_utility_HE=0, Q_recovered_kW=12300,
                 Q_hot_utility_kW=150000, Q_cold_utility_kW=80000, CAPEX_added_okuyen=1.234,
                 OPEX_utility_okuyen={}, feasible=False, message='greedy 88% しか回収できず')
res2 = SimpleNamespace(hen_result=hen2, economics_synth=object())
show_stage2_synthesis(res2)
