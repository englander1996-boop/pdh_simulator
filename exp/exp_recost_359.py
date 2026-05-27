# -*- coding: utf-8 -*-
r"""exp_recost_359.py — best #359 を現在のコストパラメータで再評価 (TAC 比較用)

cost_parameters.py を変更した後、best #359 (固定設計) を再評価して TAC/原単価/
主要 OPEX 内訳を出す。設計は不変なので duty は同じ、価格変更ぶんだけ TAC が動く。
基準 (変更前): TAC=1225.3, 原単価=304.8円/kg, OPEX=1149.2 (final_20260527_091743)。

ログ方針 ([[feedback-detailed-logging]]): python -u + flush。
使い方:  .\.venv\Scripts\python.exe -u exp\exp_recost_359.py
"""
import os, sys, json, time, datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def log(m): print(f"[{datetime.datetime.now():%H:%M:%S}] {m}", flush=True)

# 確認: 現在の蒸気/冷媒価格
from src import cost_parameters as cp
log(f"steam LP/MP/HP = {cp.LP_STEAM_JPY_PER_GJ}/{cp.MP_STEAM_JPY_PER_GJ}/{cp.HP_STEAM_JPY_PER_GJ} 円/GJ "
    f"(日本補正={cp.JAPAN_STEAM_FUEL_CORRECTION})")
log(f"refrig -100C = {cp.ETHYLENE_REFRIG_M100C_JPY_PER_GJ} 円/GJ,  fuel = {cp.FUEL_JPY_PER_GJ} 円/GJ")

log("import final ...")
import final
log("import 完了")

BEST = json.load(open('outputs/final_20260527_091743/best.json', encoding='utf-8'))
p = dict(BEST['params'])
design = final._build_design(p)
log(f"best #{BEST['number']} を再評価 (設計固定, F_fresh={p['F_C3H8_fresh_kmol_h']:.0f})...")
t0 = time.perf_counter()
res = final.evaluate(design, final._CONFIG, verbose=False,
                     apply_hi=final.APPLY_HI, hi_dT_min_K=final.HI_DT_MIN_K,
                     apply_stage2=final.APPLY_STAGE2,
                     F_C3H8_override=float(p['F_C3H8_fresh_kmol_h']))
dt = time.perf_counter() - t0
econ = res.economics_synth or res.economics_hi or res.economics
log(f"完了 ({dt:.0f}s)  feasible={res.is_feasible}")
print("=" * 60, flush=True)
print(f"  TAC          = {econ.TAC:8.1f} 億円/年   (基準 1225.3, Δ{econ.TAC-1225.3:+.1f})", flush=True)
print(f"  OPEX         = {econ.total_opex:8.1f} 億円/年   (基準 1149.2, Δ{econ.total_opex-1149.2:+.1f})", flush=True)
print(f"  CAPEX/年     = {econ.total_capex/8:8.1f} 億円/年", flush=True)
print(f"  原単価       = {econ.unit_jpy_per_t/1000:8.1f} 円/kg   (基準 304.8, Δ{econ.unit_jpy_per_t/1000-304.8:+.1f})", flush=True)
print(f"  Revenue      = {econ.total_revenue:8.1f} 億円/年", flush=True)
print(f"  Profit       = {econ.profit:8.1f} 億円/年   (基準 -392.0)", flush=True)
print("-" * 60, flush=True)
print("  蒸気/冷媒 関連 OPEX 内訳:", flush=True)
for k, v in sorted(econ.opex.items(), key=lambda kv: -kv[1]):
    if any(s in k for s in ('Steam', '蒸気', '冷媒', 'Reactor予熱', '燃料', '予熱')):
        print(f"    {k:<36} {v:8.1f}", flush=True)
print("=== DONE ===", flush=True)
