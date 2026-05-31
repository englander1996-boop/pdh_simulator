# -*- coding: utf-8 -*-
r"""exp_membrane_recycle_coupling.py — 膜 stage cut と反応器の連成感度 (レビュー指摘 #4)

狙い (2026-05-31):
  膜の透過量 (stage cut) は「製品側へ抜く C3H6」と「非透過側=反応器へ戻す C3H6」の
  配分を決める。膜面積 A_mem を増やすほど stage cut が上がり製品回収が増えるが、
  逆に小さいと多くの C3H6 が非透過側リサイクルで反応器へ戻り、反応器入口 C3H6 分率を
  押し上げる。C3H6 は脱水素反応 r1 (A→B+C) の **生成物** なので、入口 C3H6 が高いと
  平衡駆動力 (P_A − P_B·P_C/K_eq) が下がり転化率・選択率に効く。

  本スイープは **膜面積 A_mem だけを振り**、他は best 設計で固定して
  「A_mem → stage cut → 反応器入口 C3H6 → 転化率/選択率 → TAC」の連鎖を 1 枚に並べる。
  膜単体の純度だけでなく **系統全体の相互作用** を見るのが目的 (個別最適 ≠ 大域最適)。

バックエンドは main.py と同一 (Dist1/3=SM, Dist2=HYSYS, 反応器=径方向流)。
**HYSYS を使うため .venv (HYSYS 有効環境) で実行すること。**

ログ方針 ([[feedback-detailed-logging]]): 各点で開始時刻・No.・所要秒・状態・ETA を flush 出力。

使い方:
  .\.venv\Scripts\python.exe -u exp\exp_membrane_recycle_coupling.py [best.json のパス]
"""
import os, sys, json, glob, time, datetime

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass


def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


def _find_best_json():
    if len(sys.argv) > 1:
        return sys.argv[1]
    cands = sorted(glob.glob(os.path.join(_REPO, 'outputs', 'main_*', 'best.json')))
    if not cands:
        log("ERROR: best.json が見つからない。引数でパスを渡すか main.py を先に実行してください。")
        sys.exit(1)
    return cands[-1]


log("import main (config/SM モデル ロード中)...")
import main
from flowsheet import evaluate
log("import 完了")

BEST_PATH = _find_best_json()
BEST = json.load(open(BEST_PATH, encoding='utf-8'))
BASE = dict(BEST['params'])
F_FRESH = float(BASE['F_C3H8_fresh_kmol_h'])
A_BASE = float(BASE['A_mem_m2'])
log(f"基準 = {BEST_PATH}  (trial #{BEST.get('number')})  A_mem(best)={A_BASE:.3e} m^2")

# 膜面積グリッド (探索範囲 [5e4, 3e5] を跨ぐ。best も含めるよう自動挿入)
A_GRID = sorted(set([5.0e4, 1.0e5, 1.5e5, 2.0e5, 3.0e5, round(A_BASE)]))

N = len(A_GRID)
print("=" * 112, flush=True)
hdr = (f"{'A_mem[m2]':>10} {'stagecut':>8} {'透過純度%':>8} {'反応器入口C3H6%':>14} "
       f"{'X%':>5} {'S%':>5} {'prod':>6} {'純度wt%':>7} {'TAC':>8} {'feas':>5} {'秒':>5}")
print(hdr, flush=True); print("-" * 112, flush=True)

durs = []
rows = []
t_all = time.perf_counter()
for i, A in enumerate(A_GRID, 1):
    eta = (sum(durs) / len(durs) * (N - i + 1)) if durs else 0.0
    log(f"[{i}/{N}] A_mem={A:.3e} m^2 評価開始 ... 残り推定 {eta:.0f}s")
    p = dict(BASE); p['A_mem_m2'] = float(A)
    design = main._build_design(p)
    t0 = time.perf_counter()
    try:
        res = evaluate(design, main._CONFIG, verbose=False,
                       apply_hi=main.APPLY_HI, hi_dT_min_K=main.HI_DT_MIN_K,
                       apply_stage2=main.APPLY_STAGE2, F_C3H8_override=F_FRESH)
    except Exception as e:
        dt = time.perf_counter() - t0; durs.append(dt)
        log(f"[{i}/{N}] A_mem={A:.3e} EXC {type(e).__name__}: {e}  ({dt:.0f}s)")
        continue
    dt = time.perf_counter() - t0; durs.append(dt)
    op = res.solver.one_pass if res.solver else None
    if op is None:
        log(f"[{i}/{N}] A_mem={A:.3e} solver失敗 unit={res.failure_unit}  ({dt:.0f}s)")
        continue
    mem = op.get('r_mem')
    scut = getattr(mem, 'stage_cut', float('nan'))
    perm_pur = getattr(mem, 'perm_purity', float('nan')) * 100
    # 反応器入口 C3H6 mol%
    rin = op.get('reactor_inlet')
    if rin is not None:
        Ftot = sum(v for v in rin.F_in.values())
        c3h6_in = rin.F_in.get('B', 0.0) / Ftot * 100 if Ftot > 0 else float('nan')
    else:
        c3h6_in = float('nan')
    perf = op['r_rx'].performance
    X, S = perf.Conversion, perf.Selectivity
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    pur = res.specs.c3h6_purity_wtfrac * 100 if res.specs else float('nan')
    econ = res.economics_hi or res.economics
    tac = econ.TAC if econ else res.effective_TAC
    feas = 'yes' if res.is_feasible else 'NO'
    print(f"{A:>10.3e} {scut:>8.3f} {perm_pur:>8.2f} {c3h6_in:>14.2f} "
          f"{X:>5.1f} {S:>5.1f} {prod:>6.0f} {pur:>7.2f} {tac:>8.0f} {feas:>5} {dt:>5.0f}", flush=True)
    rows.append((A, scut, c3h6_in, X, S, tac, res.is_feasible))

print("=" * 112, flush=True)
if rows:
    feas_rows = [r for r in rows if r[6]]
    if feas_rows:
        best_tac = min(feas_rows, key=lambda r: r[5])
        log(f"feasible 中 TAC 最小: A_mem={best_tac[0]:.3e} (stagecut={best_tac[1]:.3f}, "
            f"反応器入口C3H6={best_tac[2]:.2f}%, X={best_tac[3]:.1f}% S={best_tac[4]:.1f}% TAC={best_tac[5]:.0f})")
    log("傾向: A_mem↓ → stagecut↓ → 反応器入口C3H6↑ → (平衡駆動力↓) の連鎖を上表で確認すること。")
else:
    log("有効データ点ゼロ")
log(f"完了。総所要 {time.perf_counter()-t_all:.0f}s")
print("=== DONE ===", flush=True)
