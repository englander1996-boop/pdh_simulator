r"""exp4.py — (SM, rigorous, SM) 構成の 1 点評価。final.py の BO を回す前の収束 smoke test。

final.py の _build_design / _CONFIG / 評価オプションを**そのまま流用**するので、
本スクリプトの 1 評価は final.py の BO 1 trial と**同一の計算**(Dist1=SM, Dist2=rigorous,
Dist3=SM, Stage2(HEN)=実行)になる。狙いは「1 点が高速に収束するか」「rigorous Dist2 +
リサイクルが per-trial 予算 120s に収まるか」を、360 trial 並列 BO を回す前に確かめること。

下の PARAMS (final.SEARCH_SPACE のキー、物理単位) を編集して再実行:
  .\.venv\Scripts\python.exe exp\exp4.py
"""

import os
import sys
import time

# 設計判断 (2026-05-27): 収束の素の挙動を見たいので時間予算は広め (abort させない)。
# 実測 wall-clock を見て final.py の 120s 予算で足りるか判断する。
# ※ final を import する前に設定する (final 側も setdefault するため先勝ち)。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '600')
os.environ.setdefault('PDH_PER_UNIT_LOG', '1')   # 各 iter の penalty 発火ユニットを stderr へ

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

import final as F                       # ← final.py の _build_design / _CONFIG を流用
from flowsheet import evaluate
from simulation import display_full_results, show_input_snapshot, hdr


# ===========================================================================
#  1 点の設計変数 (final.SEARCH_SPACE のキー、物理単位)。ここを編集して再実行。
#  - Dist1/Dist3 は SM 学習域内 (col1_n 30-60, col3_n 115-160) であること。
#  - Dist2 は rigorous (col2_n 20-40, reflux 6-10, recovery spec)。
# ===========================================================================
PARAMS = {
    # ----- 反応器 (Swing) -----
    'T_in_K':              937.0,
    'z_cat_m':             29.0,
    't_cyc_min':           19.0,
    'D_reactor_m':         9.5,
    # ----- PSA -----
    'D_psa_col_m':         4.4,
    'L_psa_bed_m':         27.8,
    'desorption_target':   0.30,
    # ----- 膜 -----
    'P_H_Pa':              8.1e5,
    'A_mem_m2':            1.2e5,
    # ----- 原料 -----
    'F_C3H8_fresh_kmol_h': 1600.0,
    # ----- Dist1 (SM, N30-60) -----
    'col1_p_kpa':          1900.0,
    'col1_n_stages':       35,
    'col1_feed_stage':     26,
    'col1_comp_frac_2':    0.95,
    # ----- Dist2 (rigorous) — main #240 近傍 (収束実績ある点) -----
    'col2_p_kpa':          720.0,
    'col2_n_stages':       39,
    'col2_reflux_ratio':   9.4,
    'col2_rec_LK_top':     0.951,
    'col2_rec_HK_bot':     0.999,
    # ----- Dist3 (SM, N115-160) -----
    'col3_p_kpa':          1680.0,
    'col3_n_stages':       150,
    'col3_feed_ratio':     0.72,
}


def _conv_line(label: str, status) -> str:
    if status is None:
        return f"  {label:18}: (情報なし)"
    conv = getattr(status, 'converged', '?')
    nit = getattr(status, 'n_iter', '?')
    return f"  {label:18}: converged={conv}  n_iter={nit}"


def main():
    design = F._build_design(PARAMS)
    config = F._CONFIG
    eval_kwargs = {'apply_hi': F.APPLY_HI, 'apply_stage2': F.APPLY_STAGE2,
                   'hi_dT_min_K': F.HI_DT_MIN_K}

    show_input_snapshot(design, config, eval_kwargs)
    hdr("exp4: (SM, rigorous, SM) 1 点評価 — 収束速度チェック")

    t0 = time.perf_counter()
    res = evaluate(
        design, config, verbose=True,
        apply_hi=F.APPLY_HI, hi_dT_min_K=F.HI_DT_MIN_K,
        apply_stage2=F.APPLY_STAGE2,
        F_C3H8_override=float(PARAMS['F_C3H8_fresh_kmol_h']),
    )
    dt = time.perf_counter() - t0

    display_full_results(res, design, config)

    # ---- 収束サマリ (本スクリプトの主目的) ----
    sv = getattr(res, 'solver', None)
    inner = getattr(sv, 'inner_status', None) if sv is not None else None
    outer = getattr(sv, 'outer_status', None) if sv is not None else None
    fu = getattr(res, 'failure_unit', None) or 'success'
    budget = float(os.environ.get('PDH_TRIAL_TIME_BUDGET_SEC', '120'))

    print("\n" + "=" * 72)
    print("  exp4 収束サマリ — (SM, rigorous, SM)")
    print("=" * 72)
    print(f"  wall-clock         : {dt:.1f} s")
    print(_conv_line("inner (recycle)", inner))
    print(_conv_line("outer (Fresh調整)", outer))
    print(f"  effective_TAC      : {res.effective_TAC:.2f} 億円/年")
    print(f"  failure_unit       : {fu}")
    try:
        a = {}  # 経済が取れれば収率等も
        prod = getattr(getattr(res, 'specs', None), 'production_kmol_h', None)
    except Exception:
        prod = None
    # final の per-trial 予算 (120s) に収まるかの判定 (本来の予算は exp4 では 600 に拡大している)
    if dt < 120:
        verdict = "OK ✅ 高速収束 (120s 予算に余裕)"
    elif dt < 300:
        verdict = "△ やや重い (120s 予算をやや超過 — 予算引上げ or Dist2 設定見直し検討)"
    else:
        verdict = "✗ 遅い (要確認 — 収束設定/設計点を見直す)"
    print(f"  → final.py per-trial 予算 120s 比: {dt:.1f}s … {verdict}")
    print("=" * 72)

    return res


if __name__ == '__main__':
    main()
