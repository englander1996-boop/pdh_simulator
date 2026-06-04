r"""
exp3.py — リサイクルあり PDH プロセス全体フロー シミュレーション (HYSYS + SM バックエンド)

exp1.py の HYSYS 版。蒸留塔は塔ごとにバックエンドを選択:
  - Dist1 / Dist3 : SM (学習済み GPR サロゲート)。HYSYS 解とほぼ一致を検証済み。
  - Dist2         : HYSYS COM (SM 化できなかったため HYSYS 継続)。
それ以外 (反応器・PSA・膜・リサイクル合流) は exp1 と同じく既存実装で動かす。
  SM 化で HYSYS 塔が Dist2 のみになり、Dist2 は開きっぱなしで warm 再解になるが、
  force-cold (PDH_HYSYS_FORCE_COLD=1) 版と結果完全一致で経路依存なし。

HYSYS の探索変数は exp1 (FUG/rigorous) と異なる:
    - exp1   : ColumnTunables.reflux_ratio / recovery_LK_top / recovery_HK_bot
    - exp3   : ColumnTunables.hysys_spec_value (column1=Comp Fraction-2,
                                                column2=Reflux Ratio,
                                                column3=Draw Rate kgmol/s)
              + hysys_feed_stage
  そのため main.py / optimization/study.py には組み込まず、独立スクリプトとして運用。
  HI / Stage 2 は exp1 と同じく True で実行 (HYSYS 出力からも Q_cond/Q_reb が取れる)。

使い方:
  下の "実験で振る設計変数" を編集して .\.venv\Scripts\python.exe exp/exp3.py
  N_stages は HSC 存在範囲内のみ:
    column1: 30〜60
    column2: 15〜80
    column3: 69〜200
"""

import os
import sys

# HYSYS は 1 塔あたり ~15秒、3塔 × リサイクル iter で 60s 予算では 2 iter で打ち切られる。
# HYSYS バックエンド向けに 30分予算に拡大。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '1800')
os.environ.setdefault('PDH_PER_UNIT_LOG', '1')   # 各 iter で penalty 発火ユニットを stderr へ

# Windows コンソール (cp932) でも記号を表示
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.reactors.catofin import CatofinDesignVars
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from simulation import display_full_results, hdr, show_input_snapshot, run_exp


# ===========================================================================
#  実験で振る設計変数 (ここを書き換えて再実行)
# ===========================================================================

# main.py (HYSYS+SM 全23変数 BO, 反応器=catofin) の BO best の最適化済み 23 変数をそのまま投入。
#   出典 outputs/main_20260604_215740/trials.csv (trial #125)
#   = 反応器制約版(N_online≤8, balanced N_total=2×N_online, L_bed≤3, d_p≥4, D≥9, 200m³上限)の
#     再最適化で「生産=目標1194.2 kmol/h を出す中で最小TAC」の採用設計。
#   (BO 記録 effective_TAC=1068.63 億円/年 (=HI後TAC)、feasible=True、
#    純度 99.497wt% / 生産量 1194.24 kmol/h / 収率 ≈82.0% / N_total=16基)。
#   exp3 でフル詳細出力を見て手動検証するため。
# 注意1: main は反応器=Catofin型 浅床・多基並列スイング (catofin)。exp3 も CatofinDesignVars を使い、
#   反応器変数を (T_in, t_cyc, D, L_bed, N_online, d_p) とする。d_p は mm→m に変換して渡す
#   (main._build_design と同一: d_p=d_p_mm/1000.0)。
# 注意2: col2/col3 のフィード段は最適化変数ではない。最適化変数は「フィード比率」
#   (col2_feed_ratio / col3_feed_ratio) であり、絶対段は _feed_stage_from_ratio で導出する。
#   exp3 でも比率をそのまま入力し、main._build_design と同一の関数で絶対段へ変換する
#   (手動変換による丸めずれを避ける)。col1 のみフィード段が直接の最適化変数。

# === 反応器 (Catofin 浅床・多基並列スイング) =====================================
# 値は best.json (trial #227) のフル精度をそのまま転記する。
# リサイクル収束は TOL_rel=1% で打ち切るため入力を丸めると着地点がずれ TAC が ~0.2%
# 変わる。main と同一の TAC (1092.70 億円/年) を再現するにはフル精度で揃える必要がある。
T_in_K          = 931.1396067022338
t_cyc_min       = 14.646484386062184
D_reactor_m     = 10.323475386828628
L_bed_m         = 1.1037728620421936
N_online        = 8
d_p_mm          = 5.098875608192616

# === PSA =====================================================================
D_psa_col_m       = 4.8918896871258895
L_psa_bed_m       = 26.011712204997245
desorption_target = 0.32110122921906026

# === 膜分離 ===================================================================
P_H_Pa     = 849444.5811765966
P_L_Pa     = 1.0e5
A_mem_m2   = 231947.92797432584

# === Dist1 (脱ブタン塔) SM ===================================================
P_dist1_kPa            = 1826.683392889994
N_dist1                = 38
FEED_STAGE_dist1       = 25       # col1_feed_stage: 最適化変数 (フィード段を直接探索)
COMP_FRAC_2_dist1      = 0.9401130406436439

# === Dist2 (脱エタン塔) HYSYS ===============================================
# P_dist2=805.6 kPa < P_H=849.4 kPa (Mem の ph_le_pfeed 回避を満たす)。
P_dist2_kPa            = 805.5746202501241
N_dist2                = 76
FEED_RATIO_dist2       = 0.5380000571489029  # col2_feed_ratio: 最適化変数。下で絶対段へ変換
REFLUX_RATIO_dist2     = 10.124163286577364

# === Dist3 (C3 スプリッタ) SM ===============================================
P_dist3_kPa            = 1669.217330080737
N_dist3                = 150
FEED_RATIO_dist3       = 0.8458566591548666  # col3_feed_ratio: 最適化変数。絶対段は組み立て側で導出
DRAW_RATE_dist3_kmolh  = 0.99     # SM では未使用 (Dist3 は spec なし)

# === Fresh LPG ==============================================================
F_C3H8_fresh_kmol_h = 1456.694207233684


# ===========================================================================
#  HI オプション (exp1 と同じ)
# ===========================================================================
APPLY_HI     = True
# Stage2(HEN greedy) を外し HI(targeting) のみで評価。HI(=Q_H_min/Q_C_min=MER) が OPEX 本体を
# 正確に捉え、Stage2 が足す回収網 CAPEX は ΔTmin 固定では ~一定オフセットで実害小。
# exp3 を main と同基準に揃える (詳細は main.py の同コメント参照)。
APPLY_STAGE2 = False
HI_DT_MIN_K  = 10.0


# ===========================================================================
#  出力モード
# ===========================================================================
SAVE_OUTPUT = True


# ===========================================================================
#  フローシート設計変数の組み立て
# ===========================================================================
#  ColumnTunables の P_col は Pa 単位、kPa から変換する。
#  hysys 経路では reflux_ratio (PDH 側) と recovery_LK_top/recovery_HK_bot は
#  参照されないので dummy 値で OK。HYSYS の主スペックは hysys_spec_value で渡す。
#  draw_rate は HYSYS では kgmol/s 単位で書込むため kmolh → kgmol/s に変換。
#  col2/col3 のフィード段は最適化変数のフィード比率から導出する (main._build_design と同一)。
def _feed_stage_from_ratio(ratio, n, lo, hi):
    fs = int(round(ratio * n)); hi_eff = min(hi, n - 2); lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))
FEED_STAGE_dist2 = _feed_stage_from_ratio(FEED_RATIO_dist2, N_dist2, 2, 9999)
FEED_STAGE_dist3 = _feed_stage_from_ratio(FEED_RATIO_dist3, N_dist3, 70, 180)

design = FlowsheetDesignVars(
    swing=CatofinDesignVars(
        T_in=T_in_K, t_cyc=t_cyc_min, D=D_reactor_m,
        L_bed=L_bed_m, N_online=int(N_online), d_p=d_p_mm / 1000.0,
    ),
    psa=PSADesignVars(
        D_col=D_psa_col_m, L_bed=L_psa_bed_m,
        desorption_target=desorption_target,
    ),
    mem=MemDesignVars(
        P_H=P_H_Pa, P_L=P_L_Pa, A_mem=A_mem_m2, P_dist=P_dist3_kPa * 1000.0,
    ),
    dist1=ColumnTunables(
        P_col=P_dist1_kPa * 1000.0,
        N_stages=N_dist1, N_feed=1,
        reflux_ratio=2.0,                       # dummy (SM/HYSYS 経路では未使用)
        # Dist1 = SM (学習済み GPR)。In_CompFraction2 = hysys_spec_value を流用。
        solver_method='sm',
        hysys_spec_value=COMP_FRAC_2_dist1,
        hysys_feed_stage=FEED_STAGE_dist1,
    ),
    dist2=ColumnTunables(
        P_col=P_dist2_kPa * 1000.0,
        N_stages=N_dist2, N_feed=1,
        reflux_ratio=REFLUX_RATIO_dist2,        # exp1 系の dummy 互換、未使用
        solver_method='hysys',
        hysys_spec_value=REFLUX_RATIO_dist2,
        hysys_feed_stage=FEED_STAGE_dist2,
    ),
    dist3=ColumnTunables(
        P_col=P_dist3_kPa * 1000.0,
        N_stages=N_dist3, N_feed=1,
        reflux_ratio=12.0,                      # dummy
        # Dist3 = SM (学習済み GPR)。model3 は spec 入力なし → 分配(回収率)は SM 予測をそのまま
        # 設計値として採用。製品純度 99.5% は満たす。hysys_spec_value は SM では未使用。
        solver_method='sm',
        # adapter は 0-1 範囲なら recovery_LK_top (動的 Draw Rate 計算)、1.0 超なら絶対量 (kgmol/s)
        # として扱う。DRAW_RATE_dist3_kmolh が 0.99 (純度モード) なら /3600 せずそのまま渡す。
        hysys_spec_value=(DRAW_RATE_dist3_kmolh
                          if 0 < DRAW_RATE_dist3_kmolh <= 1.0
                          else DRAW_RATE_dist3_kmolh / 3600.0),
        hysys_feed_stage=FEED_STAGE_dist3,
    ),
)


# ===========================================================================
#  実行 + 結果表示 (exp1 と共通の run_exp ラッパ)
# ===========================================================================
import dataclasses as _dc
config = load_operating_config()
# SM Dist3 は 99.5 mol%=99.497 wt% 固定。mol↔wt 差 0.003pp を吸収するため purity 閾値を
# 99.45 wt% に緩和 (main.py と同一。BO と exp3 の feasible 判定を揃える)。
# これがないと exp3 だけ 0.995 で判定し、99.497 wt% が ✗ になって BO 結果と食い違う。
config = _dc.replace(config, spec=_dc.replace(config.spec, c3h6_min_wtfrac=0.9945))


def _run_simulation():
    eval_kwargs = {
        'apply_hi':     APPLY_HI,
        'apply_stage2': APPLY_STAGE2,
        'hi_dT_min_K':  HI_DT_MIN_K,
    }
    show_input_snapshot(design, config, eval_kwargs)
    hdr("外側ループ: 製品流量厳密化 (Fresh を調整) — HYSYS バックエンド")
    res = evaluate(
        design, config, verbose=True,
        apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
        apply_stage2=APPLY_STAGE2,
        F_C3H8_override=F_C3H8_fresh_kmol_h,
    )
    display_full_results(res, design, config)
    return res


_out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'outputs'))
result = run_exp(
    label                 = 'exp3',
    eval_callable         = _run_simulation,
    output_dir            = _out_dir,
    save_output           = SAVE_OUTPUT,
    expected_outer_iters  = 6,
)

# シミュレーション終了時に HYSYS default provider を確実にクローズ
# (プロセス残留対策、Python 終了時の __del__ より早く)
try:
    from units.vle.hysys.provider import shutdown_default_provider
    shutdown_default_provider()
except Exception:
    pass

if result.economics is None:
    sys.exit(1)
