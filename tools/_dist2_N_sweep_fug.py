"""Ceteris paribus sweep of N_dist2 in FUG mode.

Use the BO best (#246 thin) as base, vary only N_dist2 = 22, 30, 38.
If FUG path is "N_stages doesn't change splits", we expect:
  - F_top, F_bot of Dist2 to be invariant
  - Mem feed, Dist3 feed invariant
  - yield invariant
  - Only CAPEX (vessel) and proxy_penalty differ
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'z:\pdh_simulator')

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

# Use trial #246 (thin, best feasible) as base
BASE = dict(
    T_in_K=919.1624, z_cat_m=16.8035, t_cyc_min=19.6916, D_reactor_m=9.5175,
    D_psa_col_m=3.9616, L_psa_bed_m=25.6594, desorption_target=0.2724,
    P_H_Pa=8.4341e5, P_L_Pa=1.0e5, A_mem_m2=2.8397e5,
    P_dist1_Pa=20.4932e5, N_dist1=26, reflux_dist1=2.1953,
    P_dist2_Pa=5.4111e5, N_dist2=22, reflux_dist2=8.0297,
    P_dist3_Pa=16.7200e5, N_dist3=93, reflux_dist3=11.4975,
    rec_LK_top_dist2=0.9634, rec_HK_bot_dist2=0.9997,
    F_C3H8_fresh_kmol_h=1424.135,
)

def make_design(N_dist2):
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=BASE['T_in_K'], z_cat=BASE['z_cat_m'],
                           t_cyc=BASE['t_cyc_min'], D=BASE['D_reactor_m']),
        psa=PSADesignVars(D_col=BASE['D_psa_col_m'], L_bed=BASE['L_psa_bed_m'],
                           desorption_target=BASE['desorption_target']),
        mem=MemDesignVars(P_H=BASE['P_H_Pa'], P_L=BASE['P_L_Pa'],
                           A_mem=BASE['A_mem_m2'], P_dist=BASE['P_dist3_Pa']),
        dist1=ColumnTunables(P_col=BASE['P_dist1_Pa'], N_stages=BASE['N_dist1'],
                              N_feed=1, reflux_ratio=BASE['reflux_dist1'],
                              solver_method='fug'),
        dist2=ColumnTunables(P_col=BASE['P_dist2_Pa'], N_stages=N_dist2,
                              N_feed=1, reflux_ratio=BASE['reflux_dist2'],
                              solver_method='fug',
                              recovery_LK_top=BASE['rec_LK_top_dist2'],
                              recovery_HK_bot=BASE['rec_HK_bot_dist2']),
        dist3=ColumnTunables(P_col=BASE['P_dist3_Pa'], N_stages=BASE['N_dist3'],
                              N_feed=1, reflux_ratio=BASE['reflux_dist3'],
                              solver_method='fug'),
    )

config = load_operating_config()
results = []
for N in (22, 30, 38):
    print(f"\n========== Sweep N_dist2 = {N} ==========")
    design = make_design(N)
    res = evaluate(design, config, verbose=False,
                   apply_hi=True, apply_stage2=True, hi_dT_min_K=10.0,
                   F_C3H8_override=BASE['F_C3H8_fresh_kmol_h'])
    # parse
    sr = res.solver
    one = sr.one_pass
    r2 = one['r2']
    r_mem = one['r_mem']
    r3 = one['r3']
    print(f"  Dist2 top stream F_in:   {dict((k,round(v,3)) for k,v in r2.top.F_in.items() if v>0.001)}")
    print(f"  Dist2 bot stream F_in:   {dict((k,round(v,3)) for k,v in r2.bottom.F_in.items() if v>0.001)}")
    print(f"  Mem feed (B,A):          C3H6={r_mem.product.F_C3H6:.2f} + ret C3H6={r_mem.retentate.F_C3H6:.2f}")
    print(f"  Mem stage_cut:           {r_mem.stage_cut:.4f},  perm_purity={r_mem.perm_purity:.4f}")
    print(f"  Mem product to Dist3:    C3H6={r_mem.product.F_C3H6:.2f}, C3H8={r_mem.product.F_C3H8:.2f}")
    print(f"  Dist3 top (product):     {dict((k,round(v,3)) for k,v in r3.top.F_in.items() if v>0.001)}")
    print(f"  C3H6 production:         {r3.top.F_in.get('B', 0):.2f} kmol/h")
    print(f"  Dist2 N_min:             {r2.equipment.N_min:.2f}")
    print(f"  Dist2 R_min:             {r2.equipment.R_min:.2f}")
    print(f"  Dist2 proxy_penalty:     {r2.equipment.proxy_penalty_okuyen:.2f} 億円/年")
    print(f"  Dist2 CAPEX (vessel+trays+cond+reb): {r2.equipment.CAPEX:.4f} 億円")
    print(f"  Dist2 H_col:             {r2.equipment.H_col:.2f} m")
    print(f"  feasible:               {res.is_feasible}, effective_TAC={res.effective_TAC:.2f}, failure={res.failure_reason}")
    results.append((N, r3.top.F_in.get('B',0), r2.equipment.CAPEX, r2.equipment.proxy_penalty_okuyen,
                    r_mem.stage_cut, r_mem.perm_purity, res.effective_TAC))

print("\n========== SUMMARY ==========")
print(f"{'N':>4} {'prod':>10} {'D2_CAPEX':>10} {'proxy':>10} {'stage_cut':>10} {'perm_pur':>10} {'TAC':>10}")
for r in results:
    print(f"{r[0]:>4} {r[1]:>10.3f} {r[2]:>10.4f} {r[3]:>10.4f} {r[4]:>10.4f} {r[5]:>10.4f} {r[6]:>10.3f}")
