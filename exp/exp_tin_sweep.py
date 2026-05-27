# -*- coding: utf-8 -*-
r"""exp_tin_sweep.py — 反応器 T_in 感度スイープ (転化率↔選択率トレードオフの分離評価)

final.py best #359 を起点に T_in_K だけを振り、他21変数+F_fresh は #359 値で固定。
生産規模差を打ち消す **原単価 [円/kg]** で公平比較する。

狙い (2026-05-27): trials.csv は corr(T_in,TAC)=-0.345 で「T_in↑(高転化)が net で TAC↓、
BO は上限940に張付き」を示す。よって意味があるのは **未探索の上振れ側 (T_in>940)**。
高 T_in=高転化=recycle が速く収束=高速 (best #359 wallclock 53s)。940 超で選択率↓により
TAC が下げ止まる/反転する点を探す。下振れ側 (920-930) は recycle 発散で 200s solver-abort
を量産する無駄計算なのでグリッドから除外。

ログ方針 ([[feedback-detailed-logging]]): 各点で開始時刻・No.・所要秒・状態・ETA を flush 出力。
起動は `python -u`。

使い方:  .\.venv\Scripts\python.exe -u exp\exp_tin_sweep.py
"""
import os, sys, json, time, datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)

log("import final (config/モデル ロード中)...")
import final  # _build_design / _CONFIG / evaluate / 評価オプションを流用
log("import 完了")

BEST = json.load(open('outputs/final_20260527_091743/best.json', encoding='utf-8'))
BASE = dict(BEST['params'])
# 上振れ中心: best=939.28、+1点だけ僅下(935)で傾き確認、あとは上へ
T_GRID = [935.0, 939.28, 943.0, 947.0, 951.0]

def opex_get(opex: dict, prefix: str) -> float:
    for k, v in opex.items():
        if k.startswith(prefix):
            return float(v)
    return 0.0

N = len(T_GRID)
log(f"基準 = final best #{BEST['number']}  T_in={BASE['T_in_K']:.2f}K  "
    f"F_fresh={BASE['F_C3H8_fresh_kmol_h']:.1f}kmol/h (固定)")
log(f"グリッド {N}点: {T_GRID}  (各点 evaluate: apply_hi+apply_stage2, 予算{os.environ.get('PDH_TRIAL_TIME_BUDGET_SEC')}s)")
print("=" * 112, flush=True)
hdr = (f"{'T_in':>7} {'X%':>6} {'S%':>6} {'prod':>6} {'純度wt%':>7} "
       f"{'TAC':>8} {'FreshLPG':>9} {'Dist2冷媒':>9} {'原単価':>8} {'feas':>5} {'秒':>5}")
print(hdr, flush=True); print("-" * 112, flush=True)

rows = []
durs = []
t_all = time.perf_counter()
for i, T in enumerate(T_GRID, 1):
    eta = (sum(durs) / len(durs) * (N - i + 1)) if durs else 0.0
    log(f"[{i}/{N}] T_in={T:.2f}K 計算開始 ... (これまで平均"
        f"{(sum(durs)/len(durs)) if durs else 0:.0f}s/点, 残り推定 {eta:.0f}s)")
    p = dict(BASE); p['T_in_K'] = T
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
        dt = time.perf_counter() - t0; durs.append(dt)
        log(f"[{i}/{N}] T_in={T:.2f}K EXC {type(e).__name__}: {e}  ({dt:.0f}s)")
        continue
    dt = time.perf_counter() - t0; durs.append(dt)
    econ = res.economics_synth or res.economics_hi or res.economics
    if econ is None:
        try: perf = res.solver.one_pass['r_rx'].performance; xs = f"X={perf.Conversion:.1f}% S={perf.Selectivity:.1f}%"
        except Exception: xs = ""
        log(f"[{i}/{N}] T_in={T:.2f}K solver失敗 unit={res.failure_unit} {xs}  ({dt:.0f}s)")
        continue
    perf = res.solver.one_pass['r_rx'].performance
    X, S = perf.Conversion, perf.Selectivity
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    pur = res.specs.c3h6_purity_wtfrac * 100 if res.specs else float('nan')
    fresh = opex_get(econ.opex, 'Fresh LPG')
    d2cond = opex_get(econ.opex, 'Dist2コンデンサ')
    unit = econ.unit_jpy_per_t / 1000.0  # 円/kg
    feas = 'yes' if res.is_feasible else 'NO'
    print(f"{T:>7.2f} {X:>6.1f} {S:>6.1f} {prod:>6.0f} {pur:>7.2f} "
          f"{econ.TAC:>8.0f} {fresh:>9.1f} {d2cond:>9.1f} {unit:>8.1f} {feas:>5} {dt:>5.0f}",
          flush=True)
    rows.append((T, X, S, prod, econ.TAC, fresh, d2cond, unit, res.is_feasible))

print("=" * 112, flush=True)
if rows:
    best_unit = min(rows, key=lambda r: r[7])
    log(f"原単価 最小: T_in={best_unit[0]:.2f}K → {best_unit[7]:.1f}円/kg "
        f"(S={best_unit[2]:.1f}% X={best_unit[1]:.1f}% prod={best_unit[3]:.0f} TAC={best_unit[4]:.0f})")
    base_row = next((r for r in rows if abs(r[0] - 939.28) < 0.5), None)
    if base_row:
        if best_unit[0] != base_row[0]:
            log(f"  vs best(939.28K={base_row[7]:.1f}円/kg) → 改善 {base_row[7]-best_unit[7]:.1f}円/kg")
        else:
            log("  best(939.28K) が最小 = 上振れ拡張に旨味なし")
else:
    log("有効データ点ゼロ (全点 失敗/abort)")
log(f"完了。総所要 {time.perf_counter()-t_all:.0f}s")
print("=== DONE ===", flush=True)
