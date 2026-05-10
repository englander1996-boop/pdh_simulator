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
from simulation import display_full_results, hdr, show_input_snapshot


# ===========================================================================
#  実験で振る設計変数 (ここを書き換えて再実行)
# ===========================================================================

# === 反応器 (Swing) ============================================================
#  per-pass X が 0.5 bar 反応条件で 30% 前後を確保できる組み合わせを既定値に。
#  T_in を上げる (940→970K) と X↑ だが副反応 (cracking) も増えて選択率↓。
#  z_cat を下げる (30→20m) と W_cat 比例で OPEX↓ だが X↓ でリサイクル↑。
T_in_K        = 950.0    # K       入口温度 (上限 700°C 内)
z_cat_m       = 30.0     # m       触媒層長さ
t_cyc_min     = 15.0     # min     1 サイクル運転時間
D_reactor_m   = 7.0      # m       反応器内径


# === PSA (Dist2 塔頂から H2 回収) ==============================================
#  D_col / L_bed を増やすと容積拡大で破過時間 t_abs が伸びて N_total_columns↓。
#  desorption_target を下げると脱着時間が伸びる。
D_psa_col_m       = 3.0     # m       塔径
L_psa_bed_m       = 20.0    # m       吸着層高さ
desorption_target = 0.35    # -       q が初期値の 35% まで脱着


# === 膜分離 (Dist2 塔底から C3H6/C3H8 を分離) ===================================
#  A_mem を増やすと透過量↑、stage cut↑ で C3H6 回収率↑、保留側流量↓ だが CAPEX↑。
#  P_H 上げると駆動力↑ だが 9.5 bar 上限 (Hua et al. 2024)。
P_H_Pa     = 9.5e5     # Pa      膜供給側圧力
P_L_Pa     = 1.0e5     # Pa      透過側圧力 (大気圧)
A_mem_m2   = 1.0e5     # m²      総膜面積


# === Dist1 (脱ブタン塔: C3 ←→ C4 分離) =========================================
#  PR R_min ≈ 0.95、P_col 上げると α↓ で R_min↑ だが冷却水で凝縮しやすくなる。
#  N_stages 減らすと CAPEX↓ だが N_min=12 を切ると infeasible。
P_dist1_Pa     = 17.0e5    # Pa      操作圧力 (= pump1 出口圧力に同期)
N_dist1        = 20        # -       理論段数 (N_min ≈ 12)
N_feed_dist1   = 10        # -       フィード段 (※現状 Kirkbride 推奨が優先、記録用)
reflux_dist1   = 1.5       # -       還流比 R = L/D


# === Dist2 (脱エタン塔: 軽質ガス ←→ C3) ========================================
#  z_LK=C2H4 が 0.26〜数 mol% で振れるため R_min が運転状態依存 (1.5〜4.8)。
#  P_col は膜の P_H ≤ 9.5 bar 制約から 8.5 bar が上限近傍。
P_dist2_Pa     = 8.5e5     # Pa      操作圧力 (= comp2b 出口圧力に同期、≤ 9.5 bar)
N_dist2        = 20        # -       理論段数 (N_min ≈ 1〜1.4)
N_feed_dist2   = 10        # -       フィード段
reflux_dist2   = 6.0       # -       還流比


# === Dist3 (C3 スプリッタ: C3H6 製品精製) ======================================
#  α 極小 (1.07) で OPEX 支配的 (Q_reb ~80MW)。R 下げる効果大、下限 11 まで。
#  N_stages 200 は実機並み、N_min ≈ 81。P_col は冷却水で凝縮可能な下限近く 20 bar。
P_dist3_Pa     = 20.0e5    # Pa      操作圧力 (= mem.P_dist と同期)
N_dist3        = 200       # -       理論段数 (N_min ≈ 81)
N_feed_dist3   = 100       # -       フィード段
reflux_dist3   = 12.0      # -       還流比 (R_min ≈ 10、下限 11 まで)


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
SOLVER_DIST1 = 'fug'    # 脱ブタン塔 (narrow-margin、rigorous で現実が見える)
SOLVER_DIST2 = 'fug'    # 脱エタン塔 (partial cond、rigorous で物理が正しい)
SOLVER_DIST3 = 'fug'         # C3 スプリッタ (margin 豊富、FUG で十分、Dist3 rigorous は重すぎ)


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
    ),
    dist2=ColumnTunables(
        P_col=P_dist2_Pa, N_stages=N_dist2,
        N_feed=N_feed_dist2, reflux_ratio=reflux_dist2,
        solver_method=SOLVER_DIST2,
    ),
    dist3=ColumnTunables(
        P_col=P_dist3_Pa, N_stages=N_dist3,
        N_feed=N_feed_dist3, reflux_ratio=reflux_dist3,
        solver_method=SOLVER_DIST3,
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
#  実行 + 結果表示 (ここも通常触らない)
# ===========================================================================
import io
import contextlib
from datetime import datetime
from pathlib import Path

config = load_operating_config()


def _run_simulation():
    """exp1 のメイン処理を関数化 (stdout redirect で PDF 出力対応)。"""
    eval_kwargs = {
        'apply_hi':     APPLY_HI,
        'apply_stage2': APPLY_STAGE2,
        'hi_dT_min_K':  HI_DT_MIN_K,
    }
    show_input_snapshot(design, config, eval_kwargs)
    hdr("外側ループ: 製品流量厳密化 (Fresh を調整)")
    res = evaluate(design, config, verbose=True,
                   apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
                   apply_stage2=APPLY_STAGE2)
    display_full_results(res, design, config)
    return res


if SAVE_OUTPUT:
    # 出力を文字列にキャプチャしてから txt ファイルに保存。
    # stdout を redirect する間、ターミナルには「動いてる」progress を stderr で表示。
    import threading

    timestamp = datetime.now().strftime("%Y%m%d%H%M")
    print(f"[exp1] 計算中... (出力先 outputs/exp1_{timestamp}.txt)", file=sys.stderr, flush=True)

    _stop_progress = threading.Event()
    # 外側 iter の進行を buf から読み取って % 推定
    # 通常 5-7 iter で収束するので 6 を分母に。実際の iter 数を超えても 95% で頭打ち。
    EXPECTED_OUTER_ITERS = 6

    def _progress_ticker():
        elapsed = 0
        while not _stop_progress.wait(5.0):    # 5 秒ごと更新
            elapsed += 5
            text = buf.getvalue()
            n_outer = text.count('[外側 iter ')
            # 外側 iter 数で粗く % を推定 (収束まで 1〜n の進行)
            pct = min(95, max(5, int(n_outer * 100 / EXPECTED_OUTER_ITERS)))
            phase = f"外側 iter {n_outer}/~{EXPECTED_OUTER_ITERS}" if n_outer > 0 else "初期化中"
            print(f"  ... {elapsed}s 経過、{phase} (推定 ~{pct}%)",
                  file=sys.stderr, flush=True)

    _ticker = threading.Thread(target=_progress_ticker, daemon=True)
    _ticker.start()

    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            result = _run_simulation()
    finally:
        _stop_progress.set()
        _ticker.join(timeout=1.0)

    captured = buf.getvalue()
    out_dir = Path(__file__).resolve().parent.parent / 'outputs'
    out_dir.mkdir(parents=True, exist_ok=True)
    txt_path = out_dir / f'exp1_{timestamp}.txt'
    txt_path.write_text(captured, encoding='utf-8')
    print(f"[exp1] 完了 → {txt_path}  ({len(captured.splitlines())} 行、"
          f"{txt_path.stat().st_size / 1024:.1f} KB)", file=sys.stderr, flush=True)
else:
    result = _run_simulation()

if result.economics is None:
    sys.exit(1)
