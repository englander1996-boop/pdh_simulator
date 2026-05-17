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
#  ↓ 2026-05-17: main.py BO v7-third (両側spec±2% + spec_base=50) ベスト Trial #257
#    (Profit Stage 2 ≈ -225 億円/年、yield 84.5%、生産量 target ぴったり 1178 kmol/h)
#    Note: yield は v6 (77%) より改善 ✓ だが overshoot 戦略封じで Profit は悪化
T_in_K        = 908.1870    # K       入口温度
z_cat_m       = 20.3246     # m       触媒層長さ
t_cyc_min     = 15.8635     # min     1 サイクル運転時間
D_reactor_m   = 9.5686      # m       反応器内径


# === PSA (Dist2 塔頂から H2 回収) ==============================================
D_psa_col_m       = 3.8947    # m       塔径
L_psa_bed_m       = 20.3596   # m       吸着層高さ
desorption_target = 0.1724    # -       q が初期値の 17% まで脱着


# === 膜分離 (Dist2 塔底から C3H6/C3H8 を分離) ===================================
P_H_Pa     = 6.3178e5    # Pa      膜供給側圧力 (6.3 bar)
P_L_Pa     = 1.0e5       # Pa      透過側圧力 (大気圧固定)
A_mem_m2   = 9.3262e4    # m²      総膜面積 (baseline 1e5 と同水準、小さめ)


# === Dist1 (脱ブタン塔: C3 ←→ C4 分離) =========================================
P_dist1_Pa     = 12.8345e5   # Pa      操作圧力 (12.8 bar)
N_dist1        = 24          # -       理論段数
N_feed_dist1   = 12          # -       フィード段 (Kirkbride 自動採用、本値無視)
reflux_dist1   = 1.9445      # -       還流比


# === Dist2 (脱エタン塔: 軽質ガス ←→ C3) ========================================
P_dist2_Pa     = 6.0408e5    # Pa      操作圧力 (6 bar)
N_dist2        = 27          # -       理論段数
N_feed_dist2   = 14          # -       フィード段 (Kirkbride 自動採用)
reflux_dist2   = 9.2288      # -       還流比


# === Dist3 (C3 スプリッタ: C3H6 製品精製) ======================================
P_dist3_Pa     = 20.4920e5   # Pa      操作圧力 (= mem.P_dist と同期、20.5 bar)
N_dist3        = 178         # -       理論段数
N_feed_dist3   = 89          # -       フィード段 (Kirkbride 自動採用)
reflux_dist3   = 17.5116     # -       還流比


# === Fresh LPG (BO 直接指定、外側ループ skip) ==================================
#  ↓ 2026-05-17: BO v7 ベスト = 1394.9 (両側spec で target 厳守、yield 84.5% 達成)
F_C3H8_fresh_kmol_h = 1394.9415


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
