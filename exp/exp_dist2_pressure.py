# -*- coding: utf-8 -*-
r"""exp_dist2_pressure.py — Dist2 圧力↑で極低温冷媒(278億/TAC23%)を削れるか検証

final.py best #359 を起点に **col2_p_kpa と P_H_Pa を一緒に上げる**(制約 col2_p≤P_H を維持、
best のヘッドルーム ~130kPa を踏襲)。他19変数+F_fresh は #359 固定。狙い: Dist2 塔頂温度↑
→ -100℃エチレン冷媒が安い冷媒 tier に移り、Dist2 コンデンサ OPEX(best 277.9億)が下がるか。
trials.csv の corr(col2_p,TAC)=-0.357 で方向は確実、高圧側は recycle が速く高速検証可。

ログ方針 ([[feedback-detailed-logging]]): python -u + 各点 flush。1点目で econ.opex 全キーを
ダンプ(前回 Dist2冷媒列が0=キー抽出バグだったため、実キー名を確認して substring で確実に拾う)。

使い方:  .\.venv\Scripts\python.exe -u exp\exp_dist2_pressure.py
"""
import os, sys, json, time, datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)

log("import final ...")
import final
log("import 完了")

BEST = json.load(open('outputs/final_20260527_091743/best.json', encoding='utf-8'))
BASE = dict(BEST['params'])
HEADROOM_KPA = 130.0  # P_H = col2_p + headroom (best: 902-776≈126)
COL2P_GRID = [776.5, 950.0, 1150.0, 1350.0, 1550.0]  # kPa (776.5=best)

def dist2_cond_opex(opex: dict):
    """Dist2 コンデンサの OPEX [億円] と該当キーを返す (substring で確実に)。"""
    for k, v in opex.items():
        if 'Dist2' in k and ('コンデンサ' in k or '冷媒' in k or 'cond' in k.lower()):
            return float(v), k
    return 0.0, None

N = len(COL2P_GRID)
log(f"基準 = final best #{BEST['number']}  col2_p={BASE['col2_p_kpa']:.0f}kPa "
    f"P_H={BASE['P_H_Pa']/1e3:.0f}kPa  F_fresh={BASE['F_C3H8_fresh_kmol_h']:.0f} (固定)")
log(f"グリッド {N}点 col2_p[kPa]={COL2P_GRID}  (P_H=col2_p+{HEADROOM_KPA:.0f}kPa)")
print("=" * 116, flush=True)
hdr = (f"{'col2_p':>7} {'P_H':>6} {'T_top℃':>7} {'冷媒tier':>16} "
       f"{'Dist2冷媒':>9} {'TAC':>8} {'原単価':>8} {'prod':>6} {'feas':>5} {'秒':>5}")
print(hdr, flush=True); print("-" * 116, flush=True)

rows = []; durs = []; t_all = time.perf_counter()
for i, col2p in enumerate(COL2P_GRID, 1):
    P_H = col2p + HEADROOM_KPA
    eta = (sum(durs)/len(durs)*(N-i+1)) if durs else 0.0
    log(f"[{i}/{N}] col2_p={col2p:.0f}kPa P_H={P_H:.0f}kPa 計算開始 ... (残り推定 {eta:.0f}s)")
    p = dict(BASE); p['col2_p_kpa'] = col2p; p['P_H_Pa'] = P_H * 1e3
    design = final._build_design(p)
    t0 = time.perf_counter()
    try:
        res = final.evaluate(
            design, final._CONFIG, verbose=False,
            apply_hi=final.APPLY_HI, hi_dT_min_K=final.HI_DT_MIN_K,
            apply_stage2=final.APPLY_STAGE2,
            F_C3H8_override=float(p['F_C3H8_fresh_kmol_h']),
        )
    except Exception as e:
        dt = time.perf_counter()-t0; durs.append(dt)
        log(f"[{i}/{N}] col2_p={col2p:.0f} EXC {type(e).__name__}: {e}  ({dt:.0f}s)")
        continue
    dt = time.perf_counter()-t0; durs.append(dt)
    econ = res.economics_synth or res.economics_hi or res.economics
    if econ is None:
        log(f"[{i}/{N}] col2_p={col2p:.0f} solver失敗 unit={res.failure_unit}  ({dt:.0f}s)")
        continue
    if i == 1 or not rows:  # 初回: opex 全キーをダンプ (抽出バグ診断)
        log("  econ.opex キー一覧: " + " | ".join(econ.opex.keys()))
    d2c, d2key = dist2_cond_opex(econ.opex)
    try:
        eq2 = res.solver.one_pass['r2'].equipment
        ttop = eq2.T_top - 273.15; tier = eq2.cond_utility_name
    except Exception:
        ttop = float('nan'); tier = (d2key or '').replace('Dist2コンデンサ', '').strip('() ')
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    unit = econ.unit_jpy_per_t/1000.0
    feas = 'yes' if res.is_feasible else 'NO'
    print(f"{col2p:>7.0f} {P_H:>6.0f} {ttop:>7.1f} {tier[:16]:>16} "
          f"{d2c:>9.1f} {econ.TAC:>8.0f} {unit:>8.1f} {prod:>6.0f} {feas:>5} {dt:>5.0f}",
          flush=True)
    rows.append((col2p, P_H, ttop, tier, d2c, econ.TAC, unit, prod, res.is_feasible))

print("=" * 116, flush=True)
if rows:
    base = next((r for r in rows if abs(r[0]-776.5) < 1), rows[0])
    best_tac = min(rows, key=lambda r: r[5])
    best_unit = min(rows, key=lambda r: r[6])
    log(f"基準 col2_p=776kPa: Dist2冷媒={base[4]:.1f}億 TAC={base[5]:.0f} 原単価={base[6]:.1f}")
    log(f"TAC最小: col2_p={best_tac[0]:.0f}kPa → TAC={best_tac[5]:.0f} (Dist2冷媒={best_tac[4]:.1f}億, vs基準Δ{base[5]-best_tac[5]:+.0f})")
    log(f"原単価最小: col2_p={best_unit[0]:.0f}kPa → {best_unit[6]:.1f}円/kg (vs基準Δ{base[6]-best_unit[6]:+.1f})")
else:
    log("有効データ点ゼロ")
log(f"完了。総所要 {time.perf_counter()-t_all:.0f}s")
print("=== DONE ===", flush=True)
