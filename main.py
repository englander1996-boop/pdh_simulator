r"""
main.py — PDH プロセスの多変数最適化 (Optuna ベース)

使い方:
  1. 下の § 1〜5 ブロックを編集 (試行数・ソルバ・探索範囲・出力設定)
  2. `.\.venv\Scripts\python.exe main.py` で実行
  3. outputs/main_<timestamp>_*.csv / *.json / *.txt に結果が出力される

設計判断 (2026-05-14, 相談時の合意):
  - BO ループ: 全塔 FUG で高速化 (1 eval ~3 秒)
  - top-k 再評価: 上位 k 候補だけ rigorous + Stage 2 (HEN synthesis) で精密評価
  - 18 設計変数 (整数 N_stages × 3 含む)。N_feed は core 内で Kirkbride 自動採用、
    P_L は 1 atm 固定、P_dist3 と mem.P_dist は同期で 1 変数扱い
  - SEARCH_SPACE の行をコメントアウトすると、その変数は baseline 固定になる
  - 'sm' 蒸留塔モデル追加時は SOLVER_BO/TOPK の値を 'sm' に書き換えるだけで動く

進捗:
  - Optuna の show_progress_bar=True で tqdm 表示
  - SQLite storage に履歴保存 (中断・再開可、optuna-dashboard で可視化可)

出力:
  - outputs/main_<ts>.db        : Optuna SQLite (全 trial 履歴)
  - outputs/main_<ts>_trials.csv: 全 trial の params + 診断情報
  - outputs/main_<ts>_best.json : ベスト trial の要約
  - outputs/main_<ts>_topk.txt  : top-k 候補の BO vs 再評価 比較レポート
"""

import os
import sys

# Windows コンソール (cp932) で記号も出せるように
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ===========================================================================
# § 1. 最適化ハイパラ
# ===========================================================================
N_TRIALS    = 300            # Optuna 試行回数 (FUG なら ~15 分目安)
N_STARTUP   = 50             # TPE/CMAES の冒頭ランダム前置探索 (~ n_trials / 6 目安)
N_TOPK      = 10             # top-k 再評価候補数
SEED        = 42             # 乱数シード (再現用)
SAMPLER     = 'tpe'          # 'tpe' | 'cmaes' | 'random'


# ===========================================================================
# § 2. ソルバ選択 (各塔独立、'fug' | 'rigorous' | 'sm')
# ===========================================================================
#   - BO ループは速度優先 → 全塔 FUG 推奨
#   - top-k 再評価は精度優先 → rigorous 推奨
#   - 'sm' は近日実装予定。実装後は文字列を 'sm' に変えるだけで切替可
SOLVER_BO   = {
    'dist1': 'fug',
    'dist2': 'fug',          # 設計判断 (2026-05-17): FUG path に Gilliland check 追加で narrow-margin 弾けるため戻す
    'dist3': 'fug',
}
SOLVER_TOPK = {
    'dist1': 'rigorous',
    'dist2': 'rigorous',
    'dist3': 'rigorous',
}


# ===========================================================================
# § 3. 評価オプション
# ===========================================================================
APPLY_HI                 = True    # pinch targeting (BO・top-k 共通、軽量、ms オーダー)
APPLY_STAGE2_TOPK        = True    # HEN synthesis (top-k のみ、greedy + tick-off)
HI_DT_MIN_K              = 10.0    # ピンチ最小接近温度差 [K] (textbook 標準)
STRICT_RECOVERY_BO       = False   # 全塔 FUG (Gilliland check 入り) なので無効化が高速 — top-k で rigorous 検査
STRICT_RECOVERY_TOPK     = True    # top-k rigorous の non-spec 解を catch
RECOVERY_TOLERANCE       = 0.10    # spec ±10% 許容 (top-k で使用)


# ===========================================================================
# § 4. 探索空間 — 19 変数 (元 18 + F_fresh 1)
#       形式: (low, high, scale, type)
#         scale: 'linear' | 'log'
#         type:  'float'  | 'int'
#       行をコメントアウトすればその変数は baseline 固定に
#
# Bounds 出典の凡例:
#   出典: <文献/規定>      → 物理/規定値に基づく
#   制約:                   → コードや物理から導出される境界
#   !仮置き                  → 経験的に置いた範囲 (要見直し対象)
# ===========================================================================
SEARCH_SPACE = {
    # ----- 反応器 (Swing) -----
    'T_in_K':            (900.0,  970.0,  'linear', 'float'),  # K  制約: swing.py 活性式の有効範囲 400-700°C、上限 970K (≈696.85°C, 3°C 安全マージン)
    'z_cat_m':           (15.0,   40.0,   'linear', 'float'),  # m  !仮置き (反応器最大容積 200 m³/基 制約と併せて経験範囲)
    't_cyc_min':         (12.0,   25.0,   'linear', 'float'),  # min 下限/上限: 再生 30min との比 0.4-0.83 の経験範囲、極端値は SV/cyc 整合に不利
    # 設計判断 (2026-05-19): D_reactor 下限を 7m に。4-7m だと feed 流量で
    # SV > 3 m/s に必ず違反、旧 BO 良設計はすべて 8-10m に集中していた。
    'D_reactor_m':       (7.0,    10.0,   'linear', 'float'),  # m  下限: SV 制約 ≤ 3m/s から逆算 / 上限: !仮置き

    # ----- PSA -----
    # 設計判断 (2026-05-19): D/L/desorption の物理整合領域に絞り込み。
    # 旧 BO 良設計 (trial #235, #258) はすべて D 3-5m, L 19-29m, desorption 0.25 周辺。
    # それ以外の範囲では t_abs < _T_ABS_MIN 等で penalty 発火頻発。
    'D_psa_col_m':       (2.5,    5.0,    'linear', 'float'),  # m  下限: D 小で空塔速度過大→ t_abs 極小 / 上限: 化工便覧
    'L_psa_bed_m':       (15.0,   30.0,   'linear', 'float'),  # m  下限: L 短で t_abs 極小 / 上限: !仮置き
    'desorption_target': (0.15,   0.40,   'linear', 'float'),  # -  上限: target↑で吸着時間圧縮→矛盾 / 下限: !仮置き

    # ----- 膜 (P_L は 1 atm 固定、P_dist は Dist3 と同期) -----
    # 設計判断 (2026-05-19): P_H 下限を 7.5 bar に引き上げ。Mem 圧縮機が
    # 必ず正方向 (P_H > P_dist2 = 5-7 bar) になるよう構造的に保証する。
    # 旧 5-9.5 bar だと P_dist2 (5-9.5) と重なり、50% で圧縮機逆向き penalty。
    'P_H_Pa':            (7.5e5,  9.5e5,  'linear', 'float'),  # Pa 上限: Hua et al. 9.5 bar / 下限: P_dist2 上限 7 bar + 0.5 bar margin
    'A_mem_m2':          (3.0e4,  3.0e5,  'log',    'float'),  # m² !仮置き (CAPEX 支配、log scale)

    # ----- Dist1 (脱ブタン塔) -----
    'P_dist1_Pa':        (12.0e5, 25.0e5, 'linear', 'float'),  # Pa !仮置き (pump1 出口圧と同期)
    'N_dist1':           (16,     30,     'linear', 'int'  ),  # -  !仮置き (下限: 旧14→16、N_min ≈ 12 から margin 33% / 上限: 経験値)
    'reflux_dist1':      (1.3,    3.0,    'linear', 'float'),  # -  下限: R_min ≈ 1.23 を確実に上回る (Gilliland feasible 保証) / 上限: 経験値

    # ----- Dist2 (脱エタン塔, partial cond) -----
    # 設計判断 (2026-05-19): P_dist2 上限を 7 bar に引き下げ。Mem 圧縮機が
    # 必ず正方向 (P_H ≥ 7.5 > P_dist2 ≤ 7) になるよう構造的に保証する。
    # 旧 9.5 bar だと Mem P_H レンジと重複し 50% で penalty。
    'P_dist2_Pa':        (5.0e5,  7.0e5,  'linear', 'float'),  # Pa 上限: Mem P_H 下限 7.5 - 0.5 bar margin / 下限: !仮置き
    'N_dist2':           (20,     40,     'linear', 'int'  ),  # -  !仮置き (下限: 旧10→20、rigorous で 99% recovery 物理達成可能領域 / 上限: 経験値)
    'reflux_dist2':      (5.0,    10.0,   'linear', 'float'),  # -  下限: R_min ≈ 4.5 から R/R_min ≥ 1.1 で proxy_penalty 発火多発、下限 5 で margin 1.5× 確保

    # ----- Dist3 (C3 スプリッタ, narrow-α) -----
    'P_dist3_Pa':        (15.0e5, 25.0e5, 'linear', 'float'),  # Pa !仮置き (mem.P_dist と同期、冷却水凝縮可能下限近傍)
    'N_dist3':           (120,    250,    'linear', 'int'  ),  # -  下限: N_min ≈ 60-80、N/N_min ≥ 1.3 確保 (proxy_penalty 回避) / 上限: !仮置き
    'reflux_dist3':      (11.0,   20.0,   'linear', 'float'),  # -  下限: R_min ≈ 10 + margin / 上限: !仮置き

    # ----- Fresh LPG (BO 直接指定、外側ループ skip) -----
    # 設計判断 (2026-05-17): yield 0.7-0.95 領域全体を探索可能な範囲に。
    # production target = 1188 kmol/h、両側 ±2% spec 想定:
    #   F_fresh 1200 + yield 92% = 1104 (undershoot OK 範囲外、ペナルティ ~3pp)
    #   F_fresh 1700 + yield 71% = 1207 (overshoot OK 範囲)
    #   F_fresh 1300 + yield 90% = 1170 (target 近傍 ✓)
    # 上限 1700 で yield ≥ 71% (baseline 同水準) も探索可、BO が高 yield を選好するはず。
    'F_C3H8_fresh_kmol_h': (1200.0, 1700.0, 'linear', 'float'),  # kmol/h !仮置き (yield 想定からの逆算範囲、BO 結果見ながら要調整)

    # ----- 蒸留塔 recovery -----
    # 設計判断 (2026-05-19 確定): Dist2 の C3 漏れを物理的に <1% に保証する
    # ため recovery_HK_bot_dist2 の下限を 0.998 にタイトニング。BO は「0.01
    # まで分離が保証された設計領域」内で経済最適を探す。残り 0.01 は
    # PSA/Mem の trace bypass (= 1% 閾値) が吸収する。
    # rec_LK_top_dist2 は柔軟に 0.95-0.999 とし、BO に C2H6 のリサイクル比を
    # 経済性で選ばせる。
    # 'rec_LK_top_dist1':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist1':  (0.90, 0.999, 'linear', 'float'),
    'rec_LK_top_dist2':  (0.95, 0.999, 'linear', 'float'),   # C2H6 → top (柔軟)
    'rec_HK_bot_dist2':  (0.998, 0.9999, 'linear', 'float'), # C3H8 → bot ≥ 99.8% で C3 漏れ <1% 保証
    # 'rec_LK_top_dist3':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist3':  (0.95, 0.999, 'linear', 'float'),
}


# ===========================================================================
# § 5. 出力 / 保存
# ===========================================================================
OUTPUT_DIR        = 'outputs'         # 出力先ディレクトリ (リポジトリ root 直下)
SAVE_SQLITE       = True              # Optuna SQLite (中断・再開・dashboard 用)
SAVE_TRIALS_CSV   = True              # 全 trial の履歴 CSV
SAVE_BEST_JSON    = True              # ベスト trial の要約 JSON
SAVE_TOPK_REPORT  = True              # top-k 比較レポート txt
SHOW_PROGRESS     = True              # Optuna の tqdm 進捗バー
DISPLAY_BEST_FULL = True              # top-k ベスト候補について exp1 と同じ詳細レポート出力

# ----- L1: Feasibility 分類解析 (BO 終了後の post-hoc 解析) -----
RUN_FEASIBILITY_ANALYSIS = True       # False で無効化 (sklearn 未インストール時も自動無効)
FEASIBILITY_TARGET       = 'convergence'  # 'convergence' | 'spec' | 'both'
FEASIBILITY_MODEL        = 'rf'       # 'rf' | 'logreg'


# ===========================================================================
# ↑↑↑ 編集領域はここまで。以下はパイプライン呼び出し (通常触らない) ↑↑↑
# ===========================================================================

from optimization import PipelineConfig, run_pipeline


if __name__ == '__main__':
    run_pipeline(PipelineConfig(
        # § 1
        n_trials   = N_TRIALS,
        n_startup  = N_STARTUP,
        n_topk     = N_TOPK,
        seed       = SEED,
        sampler    = SAMPLER,
        # § 2
        solver_bo   = SOLVER_BO,
        solver_topk = SOLVER_TOPK,
        # § 3
        apply_hi              = APPLY_HI,
        apply_stage2_topk     = APPLY_STAGE2_TOPK,
        hi_dT_min_K           = HI_DT_MIN_K,
        strict_recovery_bo    = STRICT_RECOVERY_BO,
        strict_recovery_topk  = STRICT_RECOVERY_TOPK,
        recovery_tolerance    = RECOVERY_TOLERANCE,
        # § 4
        search_space = SEARCH_SPACE,
        # § 5
        output_dir       = OUTPUT_DIR,
        save_sqlite      = SAVE_SQLITE,
        save_trials_csv  = SAVE_TRIALS_CSV,
        save_best_json   = SAVE_BEST_JSON,
        save_topk_report = SAVE_TOPK_REPORT,
        show_progress    = SHOW_PROGRESS,
        display_best_full = DISPLAY_BEST_FULL,
        # L1
        run_feasibility_analysis = RUN_FEASIBILITY_ANALYSIS,
        feasibility_target       = FEASIBILITY_TARGET,
        feasibility_model        = FEASIBILITY_MODEL,
    ))
