r"""exp3 の run_one_pass を直接呼んで各ユニットの結果を 1 つずつ print する。

目的: フローシート全体評価で「どのユニットで penalty が出たか」を特定する。
"""
import os
import sys

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

# 循環 import 回避
import flowsheet  # noqa: F401

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars
from flowsheet.run_one_pass import run_one_pass
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from units.vle.hysys.provider import shutdown_default_provider


# exp3 のデフォルト値をコピー
T_in_K        = 955.6260
z_cat_m       = 21.0550
t_cyc_min     = 14.2111
D_reactor_m   = 8.4363
D_psa_col_m       = 3.3187
L_psa_bed_m       = 25.6086
desorption_target = 0.2805
P_H_Pa     = 7.5403e5
P_L_Pa     = 1.0e5
A_mem_m2   = 1.2085e5
P_dist1_kPa            = 1700.0
N_dist1                = 60
FEED_STAGE_dist1       = 30
COMP_FRAC_2_dist1      = 0.99
P_dist2_kPa            = 500.0
N_dist2                = 80
FEED_STAGE_dist2       = 40
REFLUX_RATIO_dist2     = 12.0
P_dist3_kPa            = 1800.0
N_dist3                = 150
FEED_STAGE_dist3       = 75
DRAW_RATE_dist3_kmolh  = 0.99   # 純度モード (0-1)、動的計算ロジック検証
F_C3H8_fresh_kmol_h = 1647.5800


def _dump_stream(label, s):
    """ProcessStream / PSAFeed 等の流量・T・P を 1 行で表示。"""
    try:
        F = s.F_in if hasattr(s, "F_in") else s.F_out
        total = sum(F.values())
        comps = ", ".join(f"{k}={v:.1f}" for k, v in F.items() if v > 0.01)
        T_K = getattr(s, "T_in", getattr(s, "T_out", None))
        P_Pa = getattr(s, "P_in", getattr(s, "P_out", None))
        T_c = (T_K - 273.15) if T_K else None
        P_bar = (P_Pa / 1e5) if P_Pa else None
        print(f"  [{label}] F={total:.1f} kmol/h ({comps})"
              f"  T={T_c:.1f}C P={P_bar:.2f}bar"
              if T_c is not None else f"  [{label}] F={total:.1f} ({comps})")
    except Exception as e:
        print(f"  [{label}] dump failed: {e}")


def _dump_equipment(label, eq):
    """equipment オブジェクトの主要フィールドを 1 行で表示。"""
    feasible = getattr(eq, "feasible", None)
    capex = getattr(eq, "CAPEX", None)
    reactor_capex = getattr(eq, "Reactor_CAPEX", None)
    capex_total = getattr(eq, "CAPEX_total", None)
    msg = getattr(eq, "message", "") or getattr(eq, "penalty_reason", "")
    print(f"  [{label}.equipment] feasible={feasible} CAPEX={capex} "
          f"Reactor_CAPEX={reactor_capex} CAPEX_total={capex_total}")
    if msg:
        print(f"               message: {str(msg)[:200]}")


def main():
    design = FlowsheetDesignVars(
        swing=SwingDesign(T_in=T_in_K, z_cat=z_cat_m, t_cyc=t_cyc_min, D=D_reactor_m),
        psa=PSADesignVars(D_col=D_psa_col_m, L_bed=L_psa_bed_m,
                          desorption_target=desorption_target),
        mem=MemDesignVars(P_H=P_H_Pa, P_L=P_L_Pa, A_mem=A_mem_m2,
                          P_dist=P_dist3_kPa * 1000.0),
        dist1=ColumnTunables(
            P_col=P_dist1_kPa * 1000.0, N_stages=N_dist1, N_feed=1, reflux_ratio=2.0,
            solver_method='hysys',
            hysys_spec_value=COMP_FRAC_2_dist1, hysys_feed_stage=FEED_STAGE_dist1,
        ),
        dist2=ColumnTunables(
            P_col=P_dist2_kPa * 1000.0, N_stages=N_dist2, N_feed=1, reflux_ratio=REFLUX_RATIO_dist2,
            solver_method='hysys',
            hysys_spec_value=REFLUX_RATIO_dist2, hysys_feed_stage=FEED_STAGE_dist2,
        ),
        dist3=ColumnTunables(
            P_col=P_dist3_kPa * 1000.0, N_stages=N_dist3, N_feed=1, reflux_ratio=12.0,
            solver_method='hysys',
            hysys_spec_value=DRAW_RATE_dist3_kmolh / 3600.0, hysys_feed_stage=FEED_STAGE_dist3,
        ),
    )
    config = load_operating_config()

    # tear streams を妥当値で初期化 (solver の run_recycle_convergence と同等)。
    # exp1 best 状態を参考に Dist3 塔底 (C3H8 リッチ) と Mem 保留 (C3H6 リッチ) の
    # ざっくり値を入れる。0 だと JT 膨張弁が ValueError で死ぬ。
    tear_dist3 = {'A': 500.0, 'B': 50.0}   # Dist3 塔底: C3H8 主体 + 残 C3H6
    tear_mem   = {'A': 100.0, 'B': 50.0}   # Mem 保留: C3H8/C3H6
    T_d3, T_mem = 322.15, 320.15           # Mem precool 後の典型温度

    # Fresh override → 外側ループ skip と同じ条件で 1 パスだけ実行
    print("="*70)
    print("  run_one_pass 直接呼出し (初回 iter, tear=0)")
    print("="*70)
    print(f"  Fresh: C3H8={F_C3H8_fresh_kmol_h:.2f}, C4H10={F_C3H8_fresh_kmol_h * 0.1111:.2f} kmol/h")

    try:
        result = run_one_pass(
            tear_dist3=tear_dist3, tear_mem=tear_mem,
            T_d3=T_d3, T_mem=T_mem,
            F_C3H8_feed=F_C3H8_fresh_kmol_h,
            F_C4H10_feed=F_C3H8_fresh_kmol_h * 0.1111,
            design=design, config=config,
        )
    except Exception as e:
        import traceback
        print(f"\n!!! run_one_pass で例外: {type(e).__name__}: {e}")
        print(traceback.format_exc())
        shutdown_default_provider()
        return 1

    print("\n" + "="*70)
    print("  各ユニットの結果")
    print("="*70)

    # 各ユニットを順に dump
    for key in ('fresh', 'pump1', 'r1', 'dist1_top_rx',
                'recycle_dist3', 'recycle_mem', 'reactor_inlet',
                'r_rx', 'rx_out',
                'cooled', 'comp2a', 'intercool', 'comp2b', 'desuper',
                'r2', 'r_psa', 'mem_precool', 'r_mem', 'r3'):
        if key not in result:
            continue
        obj = result[key]
        print(f"\n--- {key} ---")
        # obj が stream か result-like かで分岐
        # ProcessStream は F_in 持つ
        if hasattr(obj, "F_in"):
            _dump_stream(key, obj)
        elif hasattr(obj, "outlet"):
            _dump_stream(f"{key}.outlet", obj.outlet)
            if hasattr(obj, "equipment"):
                _dump_equipment(key, obj.equipment)
        elif hasattr(obj, "top") and hasattr(obj, "bottom"):
            _dump_stream(f"{key}.top", obj.top)
            _dump_stream(f"{key}.bottom", obj.bottom)
            if hasattr(obj, "equipment"):
                _dump_equipment(key, obj.equipment)
        elif hasattr(obj, "equipment"):
            _dump_equipment(key, obj.equipment)
        else:
            print(f"  [{key}] type={type(obj).__name__}, value={obj}")

    # Mem product = Dist3 feed の詳細
    if 'r_mem' in result:
        rm = result['r_mem']
        if hasattr(rm, 'product'):
            p = rm.product
            print(f"\n--- Mem product (= Dist3 feed) ---")
            print(f"  C3H8 = {getattr(p, 'F_C3H8', 'N/A')}")
            print(f"  C3H6 = {getattr(p, 'F_C3H6', 'N/A')}")
            print(f"  T_out = {getattr(p, 'T_out', 'N/A')} K")
            print(f"  P_out = {getattr(p, 'P_out', 'N/A')} Pa")
        else:
            print(f"\n--- r_mem attrs ---")
            for attr in dir(rm):
                if not attr.startswith('_'):
                    print(f"  {attr}")

    # warnings / shortfalls
    if 'warnings_captured' in result:
        ws = result['warnings_captured']
        if ws:
            print(f"\n--- warnings_captured ({len(ws)}) ---")
            for w in ws[:10]:
                print(f"  {w}")
    for key in ('dist1_N_shortfall', 'dist2_N_shortfall', 'dist3_N_shortfall',
                'dist1_dT_shortfall', 'dist2_dT_shortfall', 'dist3_dT_shortfall',
                'reactor_sv_shortfall', 'reactor_other_shortfall',
                'psa_t_abs_shortfall', 'psa_u_0_shortfall', 'psa_feed_shortfall',
                'mem_ph_shortfall', 'mem_bp_shortfall', 'mem_phase_shortfall', 'mem_other_shortfall',
                'trace_bypass_psa_excess', 'trace_bypass_mem_excess'):
        v = result.get(key)
        if v is not None and v != 0.0:
            print(f"  shortfall {key} = {v}")

    shutdown_default_provider()
    return 0


if __name__ == "__main__":
    sys.exit(main())
