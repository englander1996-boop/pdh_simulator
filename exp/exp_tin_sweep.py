# -*- coding: utf-8 -*-
r"""exp_tin_sweep.py — 反応器 T_in 感度スイープ (選択率↔転化率トレードオフの分離評価)

final.py best #359 を起点に T_in_K だけを振り、他21変数は固定。
F_fresh も #359 値で固定 (override) し、生産規模差を打ち消す **原単価 [円/kg]** で
公平比較する。狙い: 「T_in を下げる(選択率↑)と Fresh LPG が減り、転化率↓による
recycle 増を上回って TAC(原単価)が下がるスイートスポットがあるか」を直接確認。

使い方:  .\.venv\Scripts\python.exe exp\exp_tin_sweep.py
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

import final  # _build_design / _CONFIG / evaluate / 評価オプションを流用

BEST = json.load(open('outputs/final_20260527_091743/best.json', encoding='utf-8'))
BASE = dict(BEST['params'])
T_GRID = [900.0, 910.0, 920.0, 925.0, 930.0, 935.0, 939.28, 942.0]  # 939.28=best

def opex_get(opex: dict, prefix: str) -> float:
    for k, v in opex.items():
        if k.startswith(prefix):
            return float(v)
    return 0.0

print(f"基準 = final best #{BEST['number']}  T_in={BASE['T_in_K']:.2f}K  "
      f"F_fresh={BASE['F_C3H8_fresh_kmol_h']:.1f}kmol/h (固定)")
print("=" * 108)
hdr = (f"{'T_in':>7} {'X%':>6} {'S%':>6} {'prod':>6} {'純度wt%':>7} "
       f"{'TAC':>8} {'FreshLPG':>9} {'Dist2冷媒':>9} {'原単価':>8} {'feas':>5}")
print(hdr); print("-" * 108)

rows = []
for T in T_GRID:
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
        print(f"{T:>7.1f}  EXC {type(e).__name__}: {e}")
        continue
    dt = time.perf_counter() - t0
    econ = res.economics_synth or res.economics_hi or res.economics
    if econ is None:
        try: perf = res.solver.one_pass['r_rx'].performance; xs = f"X={perf.Conversion:.1f} S={perf.Selectivity:.1f}"
        except Exception: xs = ""
        print(f"{T:>7.1f}  solver失敗 unit={res.failure_unit} {xs}  ({dt:.0f}s)")
        continue
    perf = res.solver.one_pass['r_rx'].performance
    X, S = perf.Conversion, perf.Selectivity
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    pur = res.specs.c3h6_purity_wtfrac * 100 if res.specs else float('nan')
    fresh = opex_get(econ.opex, 'Fresh LPG')
    d2cond = opex_get(econ.opex, 'Dist2コンデンサ')
    unit = econ.unit_jpy_per_t / 1000.0  # 円/kg
    feas = 'yes' if res.is_feasible else 'NO'
    print(f"{T:>7.1f} {X:>6.1f} {S:>6.1f} {prod:>6.0f} {pur:>7.2f} "
          f"{econ.TAC:>8.0f} {fresh:>9.1f} {d2cond:>9.1f} {unit:>8.1f} {feas:>5}  ({dt:.0f}s)")
    rows.append((T, X, S, prod, econ.TAC, fresh, d2cond, unit, res.is_feasible))

print("=" * 108)
if rows:
    best_unit = min(rows, key=lambda r: r[7])
    print(f"原単価 最小: T_in={best_unit[0]:.1f}K → {best_unit[7]:.1f}円/kg "
          f"(S={best_unit[2]:.1f}% X={best_unit[1]:.1f}% prod={best_unit[3]:.0f})")
    base_row = next((r for r in rows if abs(r[0] - 939.28) < 0.5), None)
    if base_row and best_unit[0] != base_row[0]:
        print(f"  vs best(939.28K) 原単価 {base_row[7]:.1f} → 改善 {base_row[7]-best_unit[7]:.1f}円/kg")
