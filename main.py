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

# 設計判断 (2026-05-21): per-trial 時間予算を 60s → 120s に引き上げ。
# 旧 60s では rigorous Dist2 + recycle iter が回るボーダーラインの trial が
# わずかなオーバーヘッドで timeout 打ち切り → CAPEX 扱いになり、本来 soft fail
# (生産量未達等) として TPE に方向シグナルが渡るはずの trial を取りこぼしていた
# (例: main_20260521_160951 trial 3 は 55.1s で完走 (TAC=390.81) → 同 params の
# 173003 trial 3 は 68.7s で timeout → CAPEX hit に劣化)。
# setdefault なのでユーザーが環境変数で別途上書き可能。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '120')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ===========================================================================
# § 1. 最適化ハイパラ
# ===========================================================================
N_TRIALS    = 300            # Optuna 試行回数 (FUG なら ~15 分目安)
N_STARTUP   = 50             # TPE startup 数 (Sobol QMC で広域カバレッジ → TPE 切替)
N_TOPK      = 10             # top-k 再評価候補数
SEED        = 42             # 乱数シード (再現用)
SAMPLER     = 'tpe'          # 'tpe' | 'cmaes' | 'random'
# 設計判断 (2026-05-21): n_jobs=1 デフォルト。並列化したい場合 2-4 (SQLite ロック注意)。
# penalty_scale が thread-local でないため最終評価は 1 推奨、中間探索は 2-4 で時短可。
N_JOBS      = 1


# ===========================================================================
# § 2. ソルバ選択 (各塔独立、'fug' | 'rigorous' | 'sm')
# ===========================================================================
#   - BO ループは速度優先 → 全塔 FUG 推奨 (旧方針)
#   - top-k 再評価は精度優先 → rigorous 推奨
#   - 'sm' は近日実装予定。実装後は文字列を 'sm' に変えるだけで切替可
#
# 設計判断 (2026-05-20): Dist2 のみ BO ループでも rigorous に切替。
#   FUG path では _split_streams が N_stages を流量計算に使わず、recovery 仕様
#   から直接 split を決定する (= N が CAPEX しか影響しない)。これにより BO が
#   常に N_dist2 下限 (=20) を選び、partial cond の C3H6 漏れが rigorous で
#   発覚しても BO ループ中は見えない構造になっていた (main_20260520_124508
#   BO best #246 = N=22 で top-10 全 rigorous infeasible が直接症状)。
#   Dist2 は N=20-40 / R=5-10 の探索範囲で rigorous でも軽量 (数秒/eval) のため
#   BO ループに組み込んでも全体時間は許容範囲。これで N が物理的に流量に効く
#   ようになり、Dist2 厚化の経済合理性が BO に見えるようになる。
#   Dist1/Dist3 は今のところ FUG 維持 (Dist1 は narrow-margin で FUG check 機能、
#   Dist3 は N が大きいため rigorous は重い)。
SOLVER_BO   = {
    'dist1': 'fug',
    'dist2': 'rigorous',     # 設計判断 (2026-05-20): FUG-rigorous gap の根源だったため切替
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
    # 設計判断 (2026-05-21 D-plan): 過去 run 履歴解析で「5/20 午前の bounds 縮小で
    # best が 167→641 に劣化」が判明。124508 run (best=167.74) と同等の b1712d1
    # bounds + N_dist3 下限 80 (124508 best は N=93) に巻き戻し、過去 top feas を
    # warm-start で先頭注入する戦略へ移行。
    'T_in_K':            (900.0,  970.0,  'linear', 'float'),  # K
    'z_cat_m':           (15.0,   40.0,   'linear', 'float'),  # m
    't_cyc_min':         (12.0,   25.0,   'linear', 'float'),  # min
    'D_reactor_m':       (7.0,    10.0,   'linear', 'float'),  # m

    # ----- PSA -----
    'D_psa_col_m':       (2.5,    5.0,    'linear', 'float'),  # m
    'L_psa_bed_m':       (15.0,   30.0,   'linear', 'float'),  # m
    'desorption_target': (0.15,   0.40,   'linear', 'float'),  # -

    # ----- 膜 (P_L は 1 atm 固定、P_dist は Dist3 と同期) -----
    'P_H_Pa':            (7.5e5,  9.5e5,  'linear', 'float'),  # Pa
    # 設計判断 (2026-05-21 D-plan revert): A_mem を b1712d1 と同じ (3e4, 3e5) log に戻す。
    # 124508 best (trial 247, TAC=167.74) は A_mem=2.84e5、5/20 縮小の上限 2.0e5 で切られてた。
    'A_mem_m2':          (3.0e4,  3.0e5,  'log',    'float'),  # m² (124508 best は 2.84e5)

    # ----- Dist1 (脱ブタン塔) -----
    'P_dist1_Pa':        (12.0e5, 25.0e5, 'linear', 'float'),  # Pa
    'N_dist1':           (16,     30,     'linear', 'int'  ),  # -
    'reflux_dist1':      (1.5,    3.0,    'linear', 'float'),  # -  維持 (dist1_N_shortfall 中央 R=1.68、yield 中立)

    # ----- Dist2 (脱エタン塔, partial cond) -----
    # 設計判断 (2026-05-21 D-plan revert): P_dist2 を b1712d1 (5.0, 7.0)e5 に戻す。
    # 124508 best は P_dist2=5.41e5、5/20 縮小の下限 5.5e5 で切られてた。
    'P_dist2_Pa':        (5.0e5,  7.0e5,  'linear', 'float'),  # Pa (124508 best は 5.41e5)
    # 設計判断 (2026-05-21 D-plan revert): N_dist2 を b1712d1 (20, 40) に戻す。
    # 124508 best は N=22、5/20 縮小の下限 30 で切られてた。
    'N_dist2':           (20,     40,     'linear', 'int'  ),  # - (124508 best は 22)
    'reflux_dist2':      (6.0,    10.0,   'linear', 'float'),  # -  維持 (Wang-Henke 収束、yield 中立)

    # ----- Dist3 (C3 スプリッタ, narrow-α) -----
    'P_dist3_Pa':        (15.0e5, 25.0e5, 'linear', 'float'),  # Pa
    # 設計判断 (2026-05-21 D-plan): N_dist3 下限を 80 に。124508 best (TAC=167.74) は N=93。
    'N_dist3':           (80,    250,    'linear', 'int'  ),  # - (124508 best は 93)
    # 設計判断 (2026-05-21 D-plan revert): reflux_dist3 を b1712d1 (11, 20) に戻す。
    # 124508 best は R=11.5、5/20 縮小の下限 14 で切られてた。
    'reflux_dist3':      (11.0,   20.0,   'linear', 'float'),  # - (124508 best は 11.5)

    # ----- Fresh LPG (BO 直接指定、外側ループ skip) -----
    'F_C3H8_fresh_kmol_h': (1200.0, 1700.0, 'linear', 'float'),  # kmol/h

    # ----- 蒸留塔 recovery -----
    # 設計判断 (2026-05-20): rec_HK_bot_dist2 の下限を 0.998 → 0.9995 に再タイト化。
    # 根拠: main_20260520_003551 trials.csv 分析より、
    #   - topk infeasible 3/3 (#32, #258, #208) は rec_HK_bot < 0.9995
    #   - topk feasible    2/2 (#115, #231)        は rec_HK_bot ≥ 0.9995
    #   - バケット 0.9995-0.9999 の中央 TAC は他バケット比 -25%
    # 「鋭利な feasibility 境界」が 0.9995 にあり、下限引き上げで sweet spot
    # 集中探索ができる。失う最小 TAC (#32 = 294.87) は実は infeasible なので実害なし。
    # 補完施策: PSA/Mem trace bypass の閾値超過に連続 penalty (runner.py
    # _TRACE_BYPASS_PENALTY_COEF_OKUYEN) を導入し、Dist2 の C3 漏れ自体に BO 中の
    # シグナルを与える。recovery 制約 + 連続 penalty の二段構えで「漏れない設計」へ誘導。
    # rec_LK_top_dist2 は柔軟に 0.95-0.999 とし、BO に C2H6 のリサイクル比を
    # 経済性で選ばせる。
    # 'rec_LK_top_dist1':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist1':  (0.90, 0.999, 'linear', 'float'),
    'rec_LK_top_dist2':  (0.95, 0.999, 'linear', 'float'),    # C2H6 → top (柔軟)
    # 設計判断 (2026-05-21 D-plan revert): 0.9995 → 0.998 に戻す。
    # 124508 best (trial 247, TAC=167.74) は rec_HK_bot=0.9997、5/20 縮小の下限 0.9995 で
    # 切られてないものの、TAC 247.37 の 050049 trial 294 (in-bounds 最良) は 0.998 ≤ rec ≤ 0.9999 領域。
    'rec_HK_bot_dist2':  (0.998, 0.9999, 'linear', 'float'),  # C3H8 → bot (124508 best は 0.9997)
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
# § 6. Warm-start trials (BO 開始前に先頭注入する既知良 params)
# ===========================================================================
# 設計判断 (2026-05-21): warm-start は無効化 (空リスト)。
#  理由 (ユーザー判断): 既知 best 周辺を最初に注入するのは「筋が悪い」=
#   - TPE が anchor された設計近傍に張り付く局所最適化バイアス
#   - 真に広い探索ができてるか不明 (warm-start が無いと出ない解は埋もれる)
#   - bounds + constraints_func + 適切な penalty 係数で十分なはずという信頼
#  以前の値を残す場合は下の cfg.warm_start_trials に list を渡せば動作する。
WARM_START_TRIALS: list = []

# 過去の warm-start 候補は参考までに残す (= 必要時に WARM_START_TRIALS に追加可)
_HISTORICAL_TOP_FEAS = [
    # main_20260520_124508/trial 247: TAC=167.74 (in-bounds best)
    {
        'T_in_K': 919.162410801484, 'z_cat_m': 16.803549423523734, 't_cyc_min': 19.69164234748162,
        'D_reactor_m': 9.517486109487328,
        'D_psa_col_m': 3.9615839767333703, 'L_psa_bed_m': 25.659406494276993,
        'desorption_target': 0.27240098624226927,
        'P_H_Pa': 843411.8987304909, 'A_mem_m2': 283966.3274966522,
        'P_dist1_Pa': 2049317.0549986823, 'N_dist1': 26, 'reflux_dist1': 2.1953437279306693,
        'P_dist2_Pa': 541106.4588241655, 'N_dist2': 22, 'reflux_dist2': 8.029650060703085,
        'P_dist3_Pa': 1672002.5591798534, 'N_dist3': 93, 'reflux_dist3': 11.497473089123655,
        'F_C3H8_fresh_kmol_h': 1424.1352972956554,
        'rec_LK_top_dist2': 0.9633748276747837, 'rec_HK_bot_dist2': 0.9997484430398813,
    },
    # main_20260520_124508/trial 248: TAC=190.98 (近傍 backup)
    {
        'T_in_K': 918.6227596868657, 'z_cat_m': 32.35969190835653, 't_cyc_min': 19.537051530398376,
        'D_reactor_m': 7.300460542156422,
        'D_psa_col_m': 4.001546050784177, 'L_psa_bed_m': 25.106352415496023,
        'desorption_target': 0.2729969985836635,
        'P_H_Pa': 843053.5687630774, 'A_mem_m2': 288953.37568181066,
        'P_dist1_Pa': 2085871.3610249786, 'N_dist1': 26, 'reflux_dist1': 2.2209824419850914,
        'P_dist2_Pa': 542883.3064569046, 'N_dist2': 22, 'reflux_dist2': 7.982540367890286,
        'P_dist3_Pa': 1818604.3843922198, 'N_dist3': 97, 'reflux_dist3': 11.553720803472697,
        'F_C3H8_fresh_kmol_h': 1414.0869419679552,
        'rec_LK_top_dist2': 0.9615480667075005, 'rec_HK_bot_dist2': 0.9997669495172944,
    },
    # main_20260520_124508/trial 250: TAC=195.77 (多様性、D_reactor=9.6 系)
    {
        'T_in_K': 917.228952528507, 'z_cat_m': 16.967494368089916, 't_cyc_min': 19.565596891185617,
        'D_reactor_m': 9.621697682208614,
        'D_psa_col_m': 3.920435945695369, 'L_psa_bed_m': 26.21306488402171,
        'desorption_target': 0.27244944212922584,
        'P_H_Pa': 849094.1932881174, 'A_mem_m2': 283383.2669744636,
        'P_dist1_Pa': 2068252.044189943, 'N_dist1': 26, 'reflux_dist1': 2.214558942500414,
        'P_dist2_Pa': 538758.100140692, 'N_dist2': 22, 'reflux_dist2': 8.060007316066612,
        'P_dist3_Pa': 1814383.9318691527, 'N_dist3': 101, 'reflux_dist3': 11.662553018381193,
        'F_C3H8_fresh_kmol_h': 1433.1511420033776,
        'rec_LK_top_dist2': 0.961580680788127, 'rec_HK_bot_dist2': 0.9997596686478493,
    },
    # main_20260520_050049/trial 294: TAC=247.37 (異なる cluster、T_in 高め)
    {
        'T_in_K': 940.6098665828874, 'z_cat_m': 33.22636519308083, 't_cyc_min': 21.2270426215022,
        'D_reactor_m': 8.223882145618532,
        'D_psa_col_m': 3.8719518101730284, 'L_psa_bed_m': 25.05409417641251,
        'desorption_target': 0.16377284557754657,
        'P_H_Pa': 803588.3583963995, 'A_mem_m2': 267750.35312782833,
        'P_dist1_Pa': 1256273.724291345, 'N_dist1': 19, 'reflux_dist1': 1.577742825591651,
        'P_dist2_Pa': 652996.4549965357, 'N_dist2': 30, 'reflux_dist2': 8.796035492833141,
        'P_dist3_Pa': 1680012.7124771352, 'N_dist3': 132, 'reflux_dist3': 11.520174810668923,
        'F_C3H8_fresh_kmol_h': 1617.146158989106,
        'rec_LK_top_dist2': 0.9747687346169066, 'rec_HK_bot_dist2': 0.9995878730438948,
    },
    # main_20260520_225416/trial 12: TAC=641.00 (現 in-bounds best、reflux_dist3=21.71 → clamp to 20.0)
    # 安全装置: 上記 4 件がすべて infeasible でも 641 を確実に再現できる
    {
        'T_in_K': 935.1845965173634, 'z_cat_m': 29.422597115658977, 't_cyc_min': 18.40273001964523,
        'D_reactor_m': 7.585728963394134,
        'D_psa_col_m': 4.3061302881537635, 'L_psa_bed_m': 19.211585436612836,
        'desorption_target': 0.15607899160786345,
        'P_H_Pa': 879094.4591814335, 'A_mem_m2': 63914.77412250522,
        'P_dist1_Pa': 2422596.1596587887, 'N_dist1': 30, 'reflux_dist1': 2.8552694633747624,
        'P_dist2_Pa': 605523.8050383166, 'N_dist2': 30, 'reflux_dist2': 9.641592812938626,
        'P_dist3_Pa': 1928184.1483173142, 'N_dist3': 196, 'reflux_dist3': 20.0,  # clamped from 21.71
        'F_C3H8_fresh_kmol_h': 1626.50472773368,
        'rec_LK_top_dist2': 0.9644279957114097, 'rec_HK_bot_dist2': 0.9996540390914408,
    },
]


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
        n_jobs     = N_JOBS,
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
        # § 6
        warm_start_trials        = WARM_START_TRIALS,
    ))
