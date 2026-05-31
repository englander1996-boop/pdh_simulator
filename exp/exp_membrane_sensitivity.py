# -*- coding: utf-8 -*-
r"""exp_membrane_sensitivity.py — 膜性能の TAC 感度スイープ (レビュー指摘 #1)

狙い (2026-05-31):
  本プロセスの成立性は膜 (ZIF-8 系 MMM) の性能に最も強く依存する。文献代表値
  (Q_A=40 GPU, α=90) は室温・大気圧・単体膜の値であり、混合ガス・高圧・経時劣化
  では低下しうる。そこで **ベスト設計を固定したまま膜性能だけを劣化させ**、TAC と
  feasibility がどれだけ頑健か (= 結論「膜あり優位」がどこで崩れるか) を定量化する。

  振るもの (one-factor-at-a-time, ベスト設計起点):
    - 透過度    Q_A_factor   ∈ {1.0, 0.5, 0.25}   (= 40 / 20 / 10 GPU 相当)
    - 選択性    alpha_factor ∈ {1.0, 0.667, 0.333} (= 90 / 60 / 30 相当)
    - 膜寿命    τ_mem [年]   ∈ {5, 3, 1}          (膜交換 OPEX = CAPEX_mem/τ。後処理で算出)
  Q_A/α は env (PDH_MEM_QA_FACTOR / PDH_MEM_ALPHA_FACTOR) でフローシートに注入する
  (flowsheet/run_one_pass.py が読む。未設定なら 1.0 = 文献代表値そのまま)。

バックエンドは main.py と同一 (Dist1/3=SM, Dist2=HYSYS, 反応器=径方向流) なので、
headline の TAC と直接比較できる。**HYSYS を使うため .venv (HYSYS 有効環境) で実行すること。**

ログ方針 ([[feedback-detailed-logging]]): 各点で開始時刻・No.・所要秒・状態・ETA を flush 出力。

使い方:
  .\.venv\Scripts\python.exe -u exp\exp_membrane_sensitivity.py [best.json のパス]
  (パス省略時は outputs/main_*/best.json の最新を自動採用)
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
import main  # main.py の _build_design / _CONFIG / 評価オプションを流用 (HYSYS+SM)
from flowsheet import evaluate
log("import 完了")

BEST_PATH = _find_best_json()
BEST = json.load(open(BEST_PATH, encoding='utf-8'))
BASE = dict(BEST['params'])
F_FRESH = float(BASE['F_C3H8_fresh_kmol_h'])
log(f"基準 = {BEST_PATH}  (trial #{BEST.get('number')})  F_fresh={F_FRESH:.1f}kmol/h")


def opex_get(opex, prefix):
    return sum(float(v) for k, v in opex.items() if k.startswith(prefix))


def run_once(qa_factor, alpha_factor):
    """env に膜劣化係数をセットして 1 回評価。(res, dt) を返す。"""
    if qa_factor == 1.0:
        os.environ.pop('PDH_MEM_QA_FACTOR', None)
    else:
        os.environ['PDH_MEM_QA_FACTOR'] = repr(qa_factor)
    if alpha_factor == 1.0:
        os.environ.pop('PDH_MEM_ALPHA_FACTOR', None)
    else:
        os.environ['PDH_MEM_ALPHA_FACTOR'] = repr(alpha_factor)
    design = main._build_design(BASE)
    t0 = time.perf_counter()
    res = evaluate(design, main._CONFIG, verbose=False,
                   apply_hi=main.APPLY_HI, hi_dT_min_K=main.HI_DT_MIN_K,
                   apply_stage2=main.APPLY_STAGE2, F_C3H8_override=F_FRESH)
    return res, time.perf_counter() - t0


# one-factor-at-a-time のグリッド (ベスト = (1.0, 1.0))
CASES = [
    ('base',       1.0,   1.0),
    ('Q_A x0.5',   0.5,   1.0),
    ('Q_A x0.25',  0.25,  1.0),
    ('alpha x0.67', 1.0,  0.6667),
    ('alpha x0.33', 1.0,  0.3333),
]
LIFE_GRID = [5.0, 3.0, 1.0]   # 膜寿命 [年] (膜交換 OPEX への影響を base 結果から後処理算出)

N = len(CASES)
print("=" * 104, flush=True)
hdr = (f"{'case':>12} {'Q_A[GPU]':>8} {'alpha':>6} {'透過純度%':>8} {'stagecut':>8} "
       f"{'prod':>6} {'純度wt%':>7} {'TAC':>8} {'feas':>5} {'秒':>5}")
print(hdr, flush=True); print("-" * 104, flush=True)

durs = []
base_res = None
t_all = time.perf_counter()
for i, (tag, fq, fa) in enumerate(CASES, 1):
    eta = (sum(durs) / len(durs) * (N - i + 1)) if durs else 0.0
    log(f"[{i}/{N}] {tag} (Q_A x{fq}, alpha x{fa}) 評価開始 ... 残り推定 {eta:.0f}s")
    try:
        res, dt = run_once(fq, fa)
    except Exception as e:
        log(f"[{i}/{N}] {tag} EXC {type(e).__name__}: {e}")
        continue
    durs.append(dt)
    if tag == 'base':
        base_res = res
    op = res.solver.one_pass if res.solver else None
    mem = op['r_mem'] if op else None
    perm_pur = mem.perm_purity * 100 if mem else float('nan')
    scut = mem.stage_cut if mem else float('nan')
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    pur = res.specs.c3h6_purity_wtfrac * 100 if res.specs else float('nan')
    econ = res.economics_hi or res.economics
    tac = econ.TAC if econ else res.effective_TAC
    feas = 'yes' if res.is_feasible else 'NO'
    print(f"{tag:>12} {40*fq:>8.1f} {90*fa:>6.1f} {perm_pur:>8.2f} {scut:>8.3f} "
          f"{prod:>6.0f} {pur:>7.2f} {tac:>8.0f} {feas:>5} {dt:>5.0f}", flush=True)

print("=" * 104, flush=True)

# 膜寿命の TAC 影響 (base 結果から後処理。膜交換 OPEX = CAPEX_mem/τ のみ τ に依存)
if base_res is not None and base_res.economics is not None:
    capex_mem = base_res.economics.capex.get('Mem膜本体', 0.0)
    base_tac = (base_res.economics_hi or base_res.economics).TAC
    base_life = float(os.environ.get('PDH_MEM_LIFETIME_YEARS', '3.0'))
    log(f"膜寿命感度 (base, 膜本体CAPEX={capex_mem:.2f}億円, 既定τ={base_life:.0f}年):")
    for tau in LIFE_GRID:
        d_opex = capex_mem * (1.0 / tau - 1.0 / base_life)   # 膜交換 OPEX の増減
        log(f"  τ={tau:.0f}年: 膜交換OPEX={capex_mem/tau:6.2f}億円/年  TAC≈{base_tac + d_opex:8.0f}億円/年 "
            f"(基準比 {d_opex:+.1f})")

log(f"完了。総所要 {time.perf_counter()-t_all:.0f}s")
print("=== DONE ===", flush=True)
