"""
exp1.py — リサイクルあり PDH プロセス全体フロー シミュレーション (実験用)

設計変数を直接編集して実験する。表示部・config ロード・evaluate 呼び出しは
simulation モジュールに移譲しているので、本ファイルは「いじる対象だけ」に集中。

例:
  反応器温度を上げてみたい  → T_in_K = 970.0 に変更してこのファイルを実行
  膜面積を増やしてみたい    → A_mem_m2 = 1.5e5 に変更
  Dist3 の還流比を見たい    → reflux_dist3 = 13.0 に変更

リサイクル構成:
  - Membrane 保留側 (C3H8 富化, 残留 C3H6 含)        ─┐
  - Dist3 塔底       (未透過 C3H8, 残留 C3H6 含)     ─┴→ Reactor 直前で合流
  Dist1 (脱ブタン) には戻さない。
"""

import os
import sys

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
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from simulation import display_full_results, hdr, show_input_snapshot, run_exp


# ===========================================================================
#  実験で振る設計変数 (ここを書き換えて再実行)
# ===========================================================================

# === 反応器 (Swing) ============================================================
#  per-pass X が 0.5 bar 反応条件で 30% 前後を確保できる組み合わせを既定値に。
#  T_in を上げる (940→970K) と X↑ だが副反応 (cracking) も増えて選択率↓。
#  z_cat を下げる (30→20m) と W_cat 比例で OPEX↓ だが X↓ でリサイクル↑。
#  ↓ 2026-05-17: main.py BO v4 ベスト Trial #294 の値を採用 (Profit -98 億円/年、baseline 比 +52 改善)
T_in_K        = 917.5081    # K       入口温度 (BO ベスト、baseline 950 から低温化で選択率↑)
z_cat_m       = 21.4408     # m       触媒層長さ
t_cyc_min     = 14.0898     # min     1 サイクル運転時間
D_reactor_m   = 6.7392      # m       反応器内径


# === PSA (Dist2 塔頂から H2 回収) ==============================================
#  D_col / L_bed を増やすと容積拡大で破過時間 t_abs が伸びて N_total_columns↓。
#  desorption_target を下げると脱着時間が伸びる。
D_psa_col_m       = 4.5648    # m       塔径 (BO ベスト)
L_psa_bed_m       = 11.9958   # m       吸着層高さ
desorption_target = 0.5258    # -       q が初期値の 52.6% まで脱着


# === 膜分離 (Dist2 塔底から C3H6/C3H8 を分離) ===================================
#  A_mem を増やすと透過量↑、stage cut↑ で C3H6 回収率↑、保留側流量↓ だが CAPEX↑。
#  P_H 上げると駆動力↑ だが 9.5 bar 上限 (Hua et al. 2024)。
P_H_Pa     = 5.5150e5    # Pa      膜供給側圧力 (BO ベスト、5.5 bar)
P_L_Pa     = 1.0e5       # Pa      透過側圧力 (大気圧固定)
A_mem_m2   = 2.1897e5    # m²      総膜面積 (BO ベスト、baseline 1e5 の 2.2倍)


# === Dist1 (脱ブタン塔: C3 ←→ C4 分離) =========================================
#  PR R_min ≈ 0.95、P_col 上げると α↓ で R_min↑ だが冷却水で凝縮しやすくなる。
#  N_stages 減らすと CAPEX↓ だが N_min=12 を切ると infeasible。
P_dist1_Pa     = 17.5864e5   # Pa      操作圧力 (BO ベスト 17.6 bar)
N_dist1        = 27          # -       理論段数 (BO ベスト)
N_feed_dist1   = 14          # -       フィード段 (Kirkbride 推奨が自動採用、本値は無視)
reflux_dist1   = 1.6975      # -       還流比 (BO ベスト)


# === Dist2 (脱エタン塔: 軽質ガス ←→ C3) ========================================
#  z_LK=C2H4 が 0.26〜数 mol% で振れるため R_min が運転状態依存 (1.5〜4.8)。
#  P_col は膜の P_H ≤ 9.5 bar 制約から 8.5 bar が上限近傍。
P_dist2_Pa     = 5.2169e5    # Pa      操作圧力 (BO ベスト 5.2 bar)
N_dist2        = 27          # -       理論段数 (rigorous で 99% 達成可能)
N_feed_dist2   = 14          # -       フィード段 (Kirkbride 自動採用)
reflux_dist2   = 6.0639      # -       還流比 (BO ベスト)


# === Dist3 (C3 スプリッタ: C3H6 製品精製) ======================================
#  α 極小 (1.07) で OPEX 支配的 (Q_reb ~80MW)。R 下げる効果大、下限 11 まで。
#  N_stages 200 は実機並み、N_min ≈ 81。P_col は冷却水で凝縮可能な下限近く 20 bar。
P_dist3_Pa     = 19.5662e5   # Pa      操作圧力 (= mem.P_dist と同期、BO ベスト)
N_dist3        = 127         # -       理論段数 (BO ベスト、N_min ≈ 81 から margin 56%)
N_feed_dist3   = 63          # -       フィード段 (Kirkbride 自動採用)
reflux_dist3   = 11.4597     # -       還流比 (R_min ≈ 10 近く、Q_reb 最小狙い)


# === Fresh LPG (BO 直接指定、外側ループ skip) ==================================
#  None  → 外側ループで自動調整 (従来動作)、target 1188 kmol/h C3H6 に張り付かせる
#  float → 指定値を使って 1 回内側ループ実行、生産量はそれに応じて変動
#  ↓ 2026-05-17: BO ベスト = 1366.7 (baseline 1666 から -18%、yield 71→86% 改善)
F_C3H8_fresh_kmol_h = 1366.7087


# === 蒸留塔 recovery (None = 0.99 既定値、float = 上書き) ======================
#  recovery_LK_top: 軽質キーが塔頂で回収される割合 (high = 損失↓ 但し N/R 大に)
#  recovery_HK_bot: 重質キーが塔底で回収される割合 (high = 純度↑ 但し N/R 大に)
#  None で 0.99 (旧 hardcode)、float (例 0.95-0.999) で上書き可能
rec_LK_top_dist1 = None    # Dist1: C3H8 in top
rec_HK_bot_dist1 = None    # Dist1: C4H10 in bottom
rec_LK_top_dist2 = None    # Dist2: C2H6 in top
rec_HK_bot_dist2 = None    # Dist2: C3H8 in bottom (PSA への C3 漏洩抑制)
rec_LK_top_dist3 = None    # Dist3: C3H6 in top
rec_HK_bot_dist3 = None    # Dist3: C3H8 in bottom (C3H6 純度に直結)


# ===========================================================================
#  蒸留塔ソルバ選択 (塔ごと個別指定)
# ===========================================================================
#  'fug'      : Fenske-Underwood-Gilliland shortcut (高速、BO 用、デフォルト)
#               recovery spec を強制達成するため出口流量が確定的。
#               narrow-margin 設計 (Dist1: N_min=12 vs N=20) では楽観的になる。
#  'rigorous' : VLE tray-by-tray (Wang-Henke、CMO 仮定、厳密)
#               recovery が物理的に決まるため narrow-margin で実態が露呈する。
#               計算重め (Dist3 N=200 narrow-α は 1 回 ~10 分かかる)。
#
# 速度コスト分析:
#   Dist1 (N=20): rigorous でも数秒
#   Dist2 (N=20): rigorous でも数秒
#   Dist3 (N=200, α=1.07): rigorous で塔単独 ~9-10 分 (narrow-α × 大量段数)
#
# 推奨運用:
#   - BO ループ: 全塔 'fug' (高速)
#   - 物理検証: Dist1/Dist2 を 'rigorous'、Dist3 は 'fug' (margin 豊富で FUG で十分、
#     実機相当 N で recovery spec も達成される)
#   - 完全厳密: 全塔 'rigorous' (~12 分/評価、デバッグ用)
SOLVER_DIST1 = 'rigorous'    # 脱ブタン塔 (narrow-margin、rigorous で現実が見える)
SOLVER_DIST2 = 'rigorous'    # 脱エタン塔 (partial cond、rigorous で物理が正しい)
SOLVER_DIST3 = 'rigorous'         # C3 スプリッタ (margin 豊富、FUG で十分、Dist3 rigorous は重すぎ)


# ===========================================================================
#  フローシート設計変数の組み立て (通常ここは触らない)
# ===========================================================================
design = FlowsheetDesignVars(
    swing=SwingDesign(
        T_in=T_in_K, z_cat=z_cat_m, t_cyc=t_cyc_min, D=D_reactor_m,
    ),
    psa=PSADesignVars(
        D_col=D_psa_col_m, L_bed=L_psa_bed_m,
        desorption_target=desorption_target,
    ),
    mem=MemDesignVars(
        # mem.P_dist は Dist3 操作圧力と一致させる必要があるため P_dist3_Pa を共有
        P_H=P_H_Pa, P_L=P_L_Pa, A_mem=A_mem_m2, P_dist=P_dist3_Pa,
    ),
    dist1=ColumnTunables(
        P_col=P_dist1_Pa, N_stages=N_dist1,
        N_feed=N_feed_dist1, reflux_ratio=reflux_dist1,
        solver_method=SOLVER_DIST1,
        recovery_LK_top=rec_LK_top_dist1,
        recovery_HK_bot=rec_HK_bot_dist1,
    ),
    dist2=ColumnTunables(
        P_col=P_dist2_Pa, N_stages=N_dist2,
        N_feed=N_feed_dist2, reflux_ratio=reflux_dist2,
        solver_method=SOLVER_DIST2,
        recovery_LK_top=rec_LK_top_dist2,
        recovery_HK_bot=rec_HK_bot_dist2,
    ),
    dist3=ColumnTunables(
        P_col=P_dist3_Pa, N_stages=N_dist3,
        N_feed=N_feed_dist3, reflux_ratio=reflux_dist3,
        solver_method=SOLVER_DIST3,
        recovery_LK_top=rec_LK_top_dist3,
        recovery_HK_bot=rec_HK_bot_dist3,
    ),
)


# ===========================================================================
#  Heat Integration (HI) オプション
# ===========================================================================
#  Stage 1 (Pinch Targeting) は BO/exp1 共通でデフォルト適用 (ms オーダーで軽量)。
#    Q_H_min/Q_C_min/A_total/N_HE_min を targeting で計算し、tier 配分で OPEX 算出。
#
#  Stage 2 (Network Synthesis = top-k 用) は離散・組合せで重め。
#    APPLY_STAGE2=True で実 HEN 構成 (greedy + tick-off) を合成し、追加 HE CAPEX
#    と実 utility OPEX を計算。top-k re-evaluation の挙動を exp1 で擬似再現する。
#    通常 BO ループには含めず、top-k 候補に対してのみ実行する想定。
#
#  False にすると該当 stage を無効化する (HI 効果・synthesis 効果のデバッグ用)。
APPLY_HI     = True
APPLY_STAGE2 = True       # top-k 候補の re-evaluation を擬似再現
HI_DT_MIN_K  = 10.0       # ピンチ最小接近温度差 (textbook 標準、BO 変数にせず固定)


# ===========================================================================
#  出力モード
# ===========================================================================
#  False: ターミナルに出力 (デフォルト、デバッグ用)
#  True : outputs/exp1_<YYYYMMDDHHMM>.txt に保存 (実験管理用)
#         ファイル名のタイムスタンプで複数実行を時系列管理可能。
#         outputs/ 配下は .gitignore で git 対象外。
SAVE_OUTPUT = True


# ===========================================================================
#  実行 + 結果表示 (ここは通常触らない、simulation.run_exp が capture/ticker/save を担当)
# ===========================================================================
config = load_operating_config()


def _run_simulation():
    """exp1 のメイン処理: 入力スナップショット → evaluate → display_full_results。"""
    eval_kwargs = {
        'apply_hi':     APPLY_HI,
        'apply_stage2': APPLY_STAGE2,
        'hi_dT_min_K':  HI_DT_MIN_K,
    }
    show_input_snapshot(design, config, eval_kwargs)
    hdr("外側ループ: 製品流量厳密化 (Fresh を調整)")
    res = evaluate(design, config, verbose=True,
                   apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
                   apply_stage2=APPLY_STAGE2,
                   F_C3H8_override=F_C3H8_fresh_kmol_h)
    display_full_results(res, design, config)
    return res


_out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'outputs'))
result = run_exp(
    label                 = 'exp1',
    eval_callable         = _run_simulation,
    output_dir            = _out_dir,
    save_output           = SAVE_OUTPUT,
    expected_outer_iters  = 6,
)

if result.economics is None:
    sys.exit(1)
