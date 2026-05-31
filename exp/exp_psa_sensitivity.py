# -*- coding: utf-8 -*-
r"""exp_psa_sensitivity.py — PSA 吸着材データ不確かさの感度スイープ (PSA設計レビュー対応)

狙い (2026-05-31):
  PSA の塔数・H2 回収率・CAPEX は吸着材データ (Langmuir 定数 q_s/a, 物質移動 KFa,
  嵩密度 ρ_b) に強く依存するが，これらはいずれも `!仮置き` でベンダーデータ未確定。
  特に KFa が 2 倍ずれると破過時間・脱着時間・必要塔数が大きく変わる。そこで
  **ベスト設計を固定したまま吸着材データを係数で振り**，塔数・H2 回収率・TAC の頑健性を
  評価する。これらが大きく揺れるなら PSA は確定設計に入れない (= 暫定扱い) と判断する。

  振るもの (one-factor-at-a-time, ベスト設計起点。env で評価器に注入):
    - KFa   (物質移動)   PDH_PSA_KFA_FACTOR  ∈ {0.5, 1.0, 2.0}
    - q_s   (飽和吸着量)  PDH_PSA_QS_FACTOR   ∈ {0.7, 1.0, 1.3}   (±30%)
    - a     (Langmuir 親和) PDH_PSA_A_FACTOR  ∈ {0.5, 1.0, 1.5}   (±50%)
    - ρ_b   (嵩密度)      PDH_PSA_RHOB_FACTOR ∈ {0.67, 1.0, 1.17} (≈400/600/700 kg/m³)

バックエンドは main.py と同一 (Dist1/3=SM, Dist2=HYSYS, 反応器=径方向流)。
**HYSYS を使うため .venv (HYSYS 有効環境) で実行すること。**

ログ方針 ([[feedback-detailed-logging]]): 各点で開始時刻・No.・所要秒・状態・ETA を flush 出力。

使い方:
  .\.venv\Scripts\python.exe -u exp\exp_psa_sensitivity.py [best.json のパス]
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
log(f"基準 = {BEST_PATH}  (trial #{BEST.get('number')})  F_fresh={F_FRESH:.1f}kmol/h")

_ENV_KEYS = ('PDH_PSA_QS_FACTOR', 'PDH_PSA_A_FACTOR', 'PDH_PSA_KFA_FACTOR', 'PDH_PSA_RHOB_FACTOR')


def run_once(env):
    for k in _ENV_KEYS:
        os.environ.pop(k, None)
    for k, v in env.items():
        os.environ[k] = repr(v)
    design = main._build_design(BASE)
    t0 = time.perf_counter()
    res = evaluate(design, main._CONFIG, verbose=False,
                   apply_hi=main.APPLY_HI, hi_dT_min_K=main.HI_DT_MIN_K,
                   apply_stage2=main.APPLY_STAGE2, F_C3H8_override=F_FRESH)
    return res, time.perf_counter() - t0


# one-factor-at-a-time (base = 全係数 1.0)
CASES = [
    ('base',        {}),
    ('KFa x0.5',    {'PDH_PSA_KFA_FACTOR': 0.5}),
    ('KFa x2.0',    {'PDH_PSA_KFA_FACTOR': 2.0}),
    ('q_s x0.7',    {'PDH_PSA_QS_FACTOR': 0.7}),
    ('q_s x1.3',    {'PDH_PSA_QS_FACTOR': 1.3}),
    ('a x0.5',      {'PDH_PSA_A_FACTOR': 0.5}),
    ('a x1.5',      {'PDH_PSA_A_FACTOR': 1.5}),
    ('rho_b 400',   {'PDH_PSA_RHOB_FACTOR': 400.0 / 600.0}),
    ('rho_b 700',   {'PDH_PSA_RHOB_FACTOR': 700.0 / 600.0}),
]

N = len(CASES)
print("=" * 110, flush=True)
hdr = (f"{'case':>12} {'N_col':>6} {'H2回収%':>7} {'CH4捕捉%':>8} {'t_abs[s]':>9} "
       f"{'t_des[s]':>9} {'ΔP[bar]':>8} {'prod':>6} {'TAC':>8} {'feas':>5} {'秒':>5}")
print(hdr, flush=True); print("-" * 110, flush=True)

durs = []
t_all = time.perf_counter()
for i, (tag, env) in enumerate(CASES, 1):
    eta = (sum(durs) / len(durs) * (N - i + 1)) if durs else 0.0
    log(f"[{i}/{N}] {tag} 評価開始 ... 残り推定 {eta:.0f}s")
    try:
        res, dt = run_once(env)
    except Exception as e:
        log(f"[{i}/{N}] {tag} EXC {type(e).__name__}: {e}")
        continue
    durs.append(dt)
    op = res.solver.one_pass if res.solver else None
    psa = op.get('r_psa') if op else None
    eq = psa.equipment if psa else None
    ncol = eq.N_total_columns if eq else -1
    h2r = psa.H2_recovery * 100 if psa else float('nan')
    ch4 = psa.CH4_capture * 100 if psa else float('nan')
    tabs = eq.t_abs_sec if eq else float('nan')
    tdes = eq.t_des_sec if eq else float('nan')
    dp = getattr(eq, 'dP_bar_actual', float('nan')) if eq else float('nan')
    prod = res.specs.production_kmol_h if res.specs else float('nan')
    econ = res.economics_hi or res.economics
    tac = econ.TAC if econ else res.effective_TAC
    feas = 'yes' if res.is_feasible else 'NO'
    reason = (eq.penalty_reason if eq and getattr(eq, 'penalty_reason', '') else '')
    note = f" ({reason})" if reason else ""
    print(f"{tag:>12} {ncol:>6} {h2r:>7.1f} {ch4:>8.2f} {tabs:>9.1f} "
          f"{tdes:>9.1f} {dp:>8.3f} {prod:>6.0f} {tac:>8.0f} {feas:>5} {dt:>5.0f}{note}", flush=True)

print("=" * 110, flush=True)
log("傾向: KFa・q_s・a・ρ_b の係数で N_col / H2回収率 / TAC が大きく揺れるなら，")
log("      PSA の塔数・回収率・CAPEX はベンダーデータ確定まで暫定扱いとする。")
log(f"完了。総所要 {time.perf_counter()-t_all:.0f}s")
print("=== DONE ===", flush=True)
