"""Trial 0 (main_20260521_131507) を rigorous Dist2 で再実行し、Wang-Henke の
失敗箇所を verbose=True と warnings.simplefilter('always') で炙り出す診断スクリプト。

実行例:
    py -3 tools\_diag_trial0_rigorous.py

設計判断 (2026-05-21): 一回限りの調査スクリプト。原因特定後は削除して良い。
"""
import os
import sys
import warnings
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.load import load_operating_config
from flowsheet import evaluate
from optimization.search_space import build_design, extract_F_fresh_override

# main.py と同じ SEARCH_SPACE / SOLVER_BO を import するため一時的に main を load
# (実行はしない、定数だけ取り出す)
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_main_mod",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "main.py"),
)
_main_mod = importlib.util.module_from_spec(_spec)
# __name__ != '__main__' なので run_pipeline は走らない
_spec.loader.exec_module(_main_mod)
SEARCH_SPACE = _main_mod.SEARCH_SPACE
SOLVER_BO = _main_mod.SOLVER_BO


TRIAL0_PARAMS = {
    "A_mem_m2": 119739.17635685026,
    "D_psa_col_m": 2.8900466011060915,
    "D_reactor_m": 8.79597545259111,
    "F_C3H8_fresh_kmol_h": 1415.972509321058,
    "L_psa_bed_m": 17.33991780504304,
    "N_dist1": 25,
    "N_dist2": 24,
    "N_dist3": 132,
    "P_H_Pa": 923235.229154987,
    "P_dist1_Pa": 2120494.351134859,
    "P_dist2_Pa": 666488.5281600844,
    "P_dist3_Pa": 1683404.5098534338,
    "T_in_K": 926.2178083193154,
    "desorption_target": 0.16452090304204986,
    "rec_HK_bot_dist2": 0.9991625204999726,
    "rec_LK_top_dist2": 0.9642702278697041,
    "reflux_dist1": 2.9548647782429915,
    "reflux_dist2": 6.7272998688284025,
    "reflux_dist3": 15.722807884690141,
    "t_cyc_min": 21.515921243548267,
    "z_cat_m": 38.767857660247905,
}


def main() -> None:
    print(f"SOLVER_BO = {SOLVER_BO}")
    print(f"params 数 = {len(TRIAL0_PARAMS)}")

    config = load_operating_config()
    design = build_design(TRIAL0_PARAMS, SOLVER_BO, baseline=None)
    F_fresh_override = extract_F_fresh_override(TRIAL0_PARAMS, baseline=None)
    print(f"F_fresh_override = {F_fresh_override}")

    # 全 warning を表示 (stderr に流す)
    warnings.simplefilter("always")

    def _showwarning(msg, cat, fn, ln, file=None, line=None):
        sys.stderr.write(f"[WARN {cat.__name__} {os.path.basename(fn)}:{ln}] {msg}\n")
    warnings.showwarning = _showwarning

    try:
        result = evaluate(
            design, config,
            verbose=True,
            apply_hi=True,
            apply_stage2=False,
            hi_dT_min_K=10.0,
            strict_recovery_check=False,
            recovery_tolerance=0.10,
            F_C3H8_override=F_fresh_override,
        )
    except Exception as e:
        print(f"\n!! evaluate() 内で未捕捉例外: {type(e).__name__}: {e}")
        traceback.print_exc()
        return

    # one_pass dict から PSA shortfall を抽出 (run_one_pass.py の _compute_psa_shortfall 結果)
    if result.solver is not None and result.solver.one_pass is not None:
        op = result.solver.one_pass
        print("\nrun_one_pass dict 内の PSA shortfall:")
        for k in ('psa_t_abs_shortfall', 'psa_u_0_shortfall', 'psa_feed_shortfall'):
            print(f"  {k} = {op.get(k)!r}")

    print("\n" + "=" * 70)
    print("evaluate() 結果")
    print("=" * 70)
    print(f"  is_feasible    : {result.is_feasible}")
    print(f"  failure_reason : {result.failure_reason}")
    print(f"  effective_TAC  : {result.effective_TAC:.2f} 億円/年")

    if result.solver is not None and result.solver.one_pass is not None:
        op = result.solver.one_pass
        print("\n各装置の CAPEX / 状態:")
        for key in ["r1", "r2", "r3", "r_psa", "r_mem", "r_rx"]:
            r = op.get(key)
            if r is None:
                print(f"  {key}: None")
                continue
            eq = getattr(r, "equipment", None)
            if eq is None:
                print(f"  {key}: equipment=None")
                continue
            capex = getattr(eq, "CAPEX_total", getattr(eq, "CAPEX", None))
            reactor_capex = getattr(eq, "Reactor_CAPEX", None)
            feas = getattr(eq, "feasible", None)
            note = []
            if capex is not None and capex >= 1e8:
                note.append(f"!! penalty sentinel CAPEX={capex:.2e}")
            if reactor_capex is not None and reactor_capex >= 1e8:
                note.append(f"!! Reactor_CAPEX={reactor_capex:.2e}")
            note_str = "  ".join(note)
            print(f"  {key}: CAPEX={capex}  Reactor_CAPEX={reactor_capex}  "
                  f"feasible={feas}  {note_str}")

        # Dist1 (r1) の詳細 (real DistEquipment であるはず)
        r1 = op.get("r1")
        if r1 is not None and r1.equipment is not None:
            eq = r1.equipment
            print("\nr1 詳細:")
            for attr in (
                "feasible", "N_stages", "N_min", "N_needed", "R", "R_min",
                "proxy_penalty_okuyen", "proxy_penalty_reason",
                "dT_max_rigorous", "CAPEX_total", "failure_reason",
            ):
                if hasattr(eq, attr):
                    print(f"    {attr} = {getattr(eq, attr)!r}")
            # 全 attr 列挙
            print("\nr1.equipment all attrs:")
            for k in dir(eq):
                if k.startswith("_"):
                    continue
                try:
                    v = getattr(eq, k)
                except Exception:
                    continue
                if callable(v):
                    continue
                print(f"    {k} = {v!r}")

        # ---- r2.top (PSA への入口) を確認 ----
        r2_res = op.get("r2")
        if r2_res is not None and getattr(r2_res, "top", None) is not None:
            print("\nr2.top (PSA inlet 候補) F_in:")
            for c, f in r2_res.top.F_in.items():
                print(f"    {c}: {f:.4f} kmol/h")
            print(f"    T_in={r2_res.top.T_in:.2f} K  P_in={r2_res.top.P_in/1e5:.3f} bar")
        # ---- PSA を直接呼んでみる ----
        try:
            print("\n" + "=" * 70)
            print("PSA 単独再現")
            print("=" * 70)
            from units.separators.psa.psa_system import simulate_psa_system, PSAFeedStream
            from flowsheet.run_one_pass import _apply_trace_bypass, _PSA_TRACE_COMPS, _TRACE_BYPASS_FRAC
            psa_in_cleaned, psa_bypass, psa_trace_excess = _apply_trace_bypass(
                r2_res.top.F_in, _PSA_TRACE_COMPS, _TRACE_BYPASS_FRAC, label="PSA",
            )
            print(f"trace bypass 適用後 F_in (PSA):")
            total = 0.0
            for c, f in psa_in_cleaned.items():
                print(f"    {c}: {f:.4f} kmol/h")
                total += f
            print(f"    合計: {total:.4f} kmol/h")
            print(f"    excess (閾値超過): {psa_trace_excess}")
            psa_feed = PSAFeedStream(
                F_in=psa_in_cleaned, T_in=r2_res.top.T_in, P_in=r2_res.top.P_in,
            )
            # PSA 内部の penalty 経路を特定するため _calc_feed_state を直接呼ぶ
            from units.separators.psa.psa_system import (
                _calc_feed_state, PSAFixedParams, _U0_MAX, _T_ABS_MIN,
                PSA_LANGMUIR_PARAMS, PSA_KFA, _ADS_ORDER, _run_adsorption,
            )
            import numpy as np
            from src.config import R as _R
            fixed = PSAFixedParams()
            print(f"\nfixed.T_abs = {fixed.T_abs} K")
            C_feed_ads, C_H2, F_non_C3_mol_s, Z = _calc_feed_state(
                psa_in_cleaned, fixed.T_abs, r2_res.top.P_in,
            )
            print(f"C_feed_ads = {C_feed_ads}")
            print(f"C_H2 = {C_H2}")
            print(f"F_non_C3_mol_s = {F_non_C3_mol_s}")
            print(f"Z = {Z}")
            import math
            A_col = math.pi/4 * design.psa.D_col**2
            u_0 = F_non_C3_mol_s * Z * _R * fixed.T_abs / (r2_res.top.P_in * A_col)
            print(f"u_0 = {u_0:.4f} m/s   (_U0_MAX={_U0_MAX})")
            print(f"C_feed_ads[0] (CH4) = {C_feed_ads[0]}")
            print(f"_T_ABS_MIN = {_T_ABS_MIN}")

            # _run_adsorption を直接呼んで t_abs と converged を確認
            q_s    = np.array([PSA_LANGMUIR_PARAMS[k]['q_s'] for k in _ADS_ORDER])
            a_lang = np.array([PSA_LANGMUIR_PARAMS[k]['a']   for k in _ADS_ORDER])
            kfa    = np.array([PSA_KFA[k]                    for k in _ADS_ORDER])
            t_abs_clean, q_final, sol_t, C_outlet_t, converged = _run_adsorption(
                C_feed=C_feed_ads,
                u_0=u_0,
                L_bed=design.psa.L_bed,
                rho_b=fixed.rho_b,
                eps=fixed.eps,
                kfa=kfa,
                q_s=q_s,
                a_lang=a_lang,
                breakthrough_ratio=fixed.breakthrough_ratio,
                t_max=fixed.t_ads_max,
            )
            print(f"\n_run_adsorption 結果:")
            print(f"  t_abs_clean = {t_abs_clean} s")
            print(f"  converged   = {converged}")
            print(f"  sol_t.size  = {sol_t.size}")
            print(f"  fixed.t_ads_max = {fixed.t_ads_max}")
            print(f"  fixed.breakthrough_ratio = {fixed.breakthrough_ratio}")
            print(f"  desorption_target = {design.psa.desorption_target}")
            t_abs_css = t_abs_clean * (1.0 - design.psa.desorption_target) if fixed.use_css_approximation else t_abs_clean
            print(f"  t_abs (CSS 補正後) = {t_abs_css} s  (vs _T_ABS_MIN={_T_ABS_MIN})")

            psa_result = simulate_psa_system(design.psa, psa_feed)
            print("\nPSA result:")
            eq = psa_result.equipment
            for attr in ("CAPEX_total", "feasible", "t_abs_s", "t_des_s",
                         "N_total_columns", "D_col_m", "L_bed_m"):
                if hasattr(eq, attr):
                    print(f"    {attr} = {getattr(eq, attr)!r}")
        except Exception as e:
            print(f"PSA 直接呼出で例外: {type(e).__name__}: {e}")
            traceback.print_exc()



if __name__ == "__main__":
    main()
