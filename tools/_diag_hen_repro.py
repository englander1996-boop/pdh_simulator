# -*- coding: utf-8 -*-
r"""#294 の Stage2 HEN 再現性診断 (2026-05-26)

問い: 同じ #294 が BO実行時 TAC=1305 / レポート再評価 TAC=1391 と揺れた。なぜか。
仮説: フローシート評価でストリームが微妙にズレ (リサイクル/HYSYS warm 経路)、
      greedy HEN がそれを回収率 100%↔78% の段差に増幅している。

本スクリプトの切り分け:
  Phase 1: #294 を N 回 evaluate (HYSYS) → Stage1/Stage2 TAC・HEN回収率と
           「ストリーム指紋」を記録し、run 間で変わるかを見る (フローシート経路依存)。
  Phase 2: Phase1 の 1 回分のストリームを捕捉し synthesize_hen を M 回 +
           微小摂動で回す (HYSYS 不要)。greedy が決定論か・崖感度を見る。

使い方: .\.venv\Scripts\python.exe tools\_diag_hen_repro.py [N]
"""
import os
import sys
import dataclasses as dc

os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '1800')
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from flowsheet.heat_integration import (
    extract_streams, pinch_analysis, get_default_utility_tiers,
)
from optimization.hen_synthesis import synthesize_hen


# === #294 設計 (outputs/special_20260526_172608_best.json) ===
def build_294() -> FlowsheetDesignVars:
    p3_kpa = 1677.1573
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=936.9962277, z_cat=29.1345707,
                          t_cyc=19.0206554, D=9.4751160),
        psa=PSADesignVars(D_col=4.4154197, L_bed=27.8226005,
                          desorption_target=0.2974322),
        mem=MemDesignVars(P_H=8.0711387e5, P_L=1.0e5, A_mem=1.1854758e5,
                          P_dist=p3_kpa * 1000.0),
        dist1=ColumnTunables(P_col=1996.1098 * 1000.0, N_stages=31, N_feed=1,
                             reflux_ratio=2.0, solver_method='sm',
                             hysys_spec_value=0.9083715, hysys_feed_stage=28),
        dist2=ColumnTunables(P_col=645.4137 * 1000.0, N_stages=67, N_feed=1,
                             reflux_ratio=8.9916122, solver_method='hysys',
                             hysys_spec_value=8.9916122, hysys_feed_stage=37),
        dist3=ColumnTunables(P_col=p3_kpa * 1000.0, N_stages=156, N_feed=1,
                             reflux_ratio=12.0, solver_method='sm',
                             hysys_spec_value=0.99, hysys_feed_stage=112),
    )


F_FRESH = 1523.3995
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3

config = load_operating_config()
config = dc.replace(config, spec=dc.replace(config.spec, c3h6_min_wtfrac=0.9945))
design = build_294()


def stream_fingerprint(streams):
    """ストリーム集合を (name -> (T_in,T_out,FCp,Q_lat,phase)) でソート表現。"""
    rows = []
    for s in sorted(streams, key=lambda x: x.name):
        rows.append((s.name, s.T_in_K, s.T_out_K, s.F_Cp_kW_per_K,
                     s.Q_latent_kW, s.phase))
    return rows


def fp_maxdiff(fp_a, fp_b):
    """2 つの指紋の数値最大差 (名前集合が同じ前提)。違う名前があれば inf。"""
    da = {r[0]: r[1:5] for r in fp_a}
    db = {r[0]: r[1:5] for r in fp_b}
    if set(da) != set(db):
        return float('inf'), (set(da) ^ set(db))
    md = 0.0
    for k in da:
        for x, y in zip(da[k], db[k]):
            md = max(md, abs(float(x) - float(y)))
    return md, None


print("=" * 78)
print(f"  Phase 1: #294 を {N} 回 evaluate (HYSYS)。Stage1/Stage2/回収率/ストリーム指紋")
print("=" * 78)

heat_t, cool_t = get_default_utility_tiers()
records = []
for i in range(N):
    res = evaluate(design, config, verbose=False, apply_hi=True,
                   hi_dT_min_K=10.0, apply_stage2=True, F_C3H8_override=F_FRESH)
    s1 = res.economics_hi.TAC if res.economics_hi else float('nan')
    s2 = res.economics_synth.TAC if res.economics_synth else float('nan')
    hr = res.hen_result
    streams = extract_streams(res.solver.one_pass, design.swing.T_in)
    hot_total = sum(s.Q_total_kW for s in streams if s.is_hot)
    hi = res.hi_result
    target = hot_total - (hi.Q_C_min_kW if hi else 0.0)
    rec = hr.Q_recovered_kW if hr else float('nan')
    pct = 100.0 * rec / target if target > 0 else float('nan')
    fp = stream_fingerprint(streams)
    records.append({'s1': s1, 's2': s2, 'rec': rec, 'target': target,
                    'pct': pct, 'feas': hr.feasible if hr else None,
                    'msg': hr.message if hr else '', 'fp': fp,
                    'one_pass': res.solver.one_pass})
    print(f"  run {i+1}: Stage1={s1:8.2f}  Stage2={s2:8.2f}  "
          f"HEN回収={rec:8.0f}/{target:8.0f}kW ({pct:5.1f}%)  feas={hr.feasible if hr else None}")
    if hr and hr.message:
        print(f"         msg: {hr.message}")

print("\n  --- run 間のストリーム指紋差 (run1 基準の最大数値差) ---")
for i in range(1, len(records)):
    md, missing = fp_maxdiff(records[0]['fp'], records[i]['fp'])
    if missing:
        print(f"    run1 vs run{i+1}: 名前集合が違う! 差分={missing}")
    else:
        print(f"    run1 vs run{i+1}: 最大|Δ| = {md:.6g}  "
              f"(Stage2差={records[i]['s2']-records[0]['s2']:+.2f}億, "
              f"回収率差={records[i]['pct']-records[0]['pct']:+.1f}pp)")

# Phase 2: 捕捉した run1 のストリームで greedy 決定論 + 崖感度
print("\n" + "=" * 78)
print("  Phase 2: run1 のストリームで synthesize_hen を反復 (HYSYS 不要)")
print("=" * 78)
op0 = records[0]['one_pass']
streams0 = extract_streams(op0, design.swing.T_in)
hi0 = pinch_analysis(streams0, dT_min_K=10.0, heating_tiers=heat_t, cooling_tiers=cool_t)

print("  (a) 同一ストリームで 3 回 → greedy 決定論チェック")
for j in range(3):
    hen = synthesize_hen(streams0, hi0, dT_min_K=10.0,
                         heating_tiers=heat_t, cooling_tiers=cool_t)
    print(f"      回 {j+1}: Q_recovered={hen.Q_recovered_kW:10.2f}kW  "
          f"n_HE={hen.n_process_HE}  feasible={hen.feasible}")

print("\n  (b) 全ストリームの T を ±εK 摂動 → 回収率の崖感度")
print("      (修正: T_in/T_out に加え T_phase_K も同量シフトし潜熱流の整合を保つ)")
def _shift(s, eps):
    tp = (s.T_phase_K + eps) if s.T_phase_K is not None else None
    return dc.replace(s, T_in_K=s.T_in_K + eps, T_out_K=s.T_out_K + eps, T_phase_K=tp)
for eps in (0.0, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0,
            -0.001, -0.01, -0.05, -0.1, -0.5, -1.0):
    pert = [_shift(s, eps) for s in streams0]
    hi_p = pinch_analysis(pert, dT_min_K=10.0, heating_tiers=heat_t, cooling_tiers=cool_t)
    hen = synthesize_hen(pert, hi_p, dT_min_K=10.0,
                         heating_tiers=heat_t, cooling_tiers=cool_t)
    hot_total = sum(s.Q_total_kW for s in pert if s.is_hot)
    target = hot_total - hi_p.Q_C_min_kW
    pct = 100.0 * hen.Q_recovered_kW / target if target > 0 else float('nan')
    print(f"      ΔT={eps:+6.3f}K: Q_recovered={hen.Q_recovered_kW:10.2f}kW "
          f"({pct:5.1f}%)  target={target:9.0f}  n_HE={hen.n_process_HE}  feasible={hen.feasible}")

# Phase 3: 別設計を 1 回挟んでから #294 を再評価 (BO loop の warm 状態を模擬)
print("\n" + "=" * 78)
print("  Phase 3: warm 状態依存テスト [#294 → 別設計B → #294]")
print("  (BO loop / report 再評価で別 trial 後に #294 を解く状況を模擬)")
print("=" * 78)

def _eval_294():
    res = evaluate(design, config, verbose=False, apply_hi=True,
                   hi_dT_min_K=10.0, apply_stage2=True, F_C3H8_override=F_FRESH)
    s2 = res.economics_synth.TAC if res.economics_synth else float('nan')
    hr = res.hen_result
    streams = extract_streams(res.solver.one_pass, design.swing.T_in)
    hot_total = sum(s.Q_total_kW for s in streams if s.is_hot)
    target = hot_total - (res.hi_result.Q_C_min_kW if res.hi_result else 0.0)
    pct = 100.0 * hr.Q_recovered_kW / target if (hr and target > 0) else float('nan')
    return s2, (hr.Q_recovered_kW if hr else float('nan')), pct, (hr.feasible if hr else None)

# 別設計B: #294 の Dist2 だけ変えた点 (HYSYS を別 case で上書き warm 状態を変える)
design_B = dc.replace(
    design,
    dist2=dc.replace(design.dist2, P_col=600.0 * 1000.0,
                     reflux_ratio=10.0, hysys_spec_value=10.0),
)
s2a, reca, pcta, fa = _eval_294()
print(f"  #294 (A, warm from #294): Stage2={s2a:8.2f}  回収={reca:9.0f}kW ({pcta:5.1f}%)  feas={fa}")
resB = evaluate(design_B, config, verbose=False, apply_hi=True, hi_dT_min_K=10.0,
                apply_stage2=True, F_C3H8_override=F_FRESH)
s2B = resB.economics_synth.TAC if resB.economics_synth else float('nan')
print(f"  別設計B (Dist2 P=600/R=10):     Stage2={s2B if s2B==s2B else 'NA'}  feas={resB.is_feasible}")
s2c, recc, pctc, fc = _eval_294()
print(f"  #294 (C, warm from B):     Stage2={s2c:8.2f}  回収={recc:9.0f}kW ({pctc:5.1f}%)  feas={fc}")
print(f"\n  >>> #294 A vs C: Stage2差={s2c-s2a:+.2f}億  回収率差={pctc-pcta:+.1f}pp")
print(f"      (差が大 → warm 状態依存で同一設計が別 TAC = 観測された 1305/1391 の正体)")

# クリーンアップ
try:
    from units.vle.hysys.provider import shutdown_default_provider
    shutdown_default_provider()
except Exception:
    pass
