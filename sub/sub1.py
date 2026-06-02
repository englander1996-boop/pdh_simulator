r"""
sub/sub1.py — PDH プロセスの多変数最適化 (Optuna ベース、FUG/rigorous)

FUG/rigorous 全フローシート版のアーカイブ (並列対応、parallel kind='sub1')。
現行 BO の本丸は ../main.py (HYSYS+SM)。

使い方:
  1. 下の § 1〜5 ブロックを編集 (試行数・ソルバ・探索範囲・出力設定)
  2. `.\.venv\Scripts\python.exe sub\sub1.py` で実行
  3. outputs/main_<timestamp>/ に結果が出力される

構成:
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

# per-trial 時間予算 120s。短いと rigorous Dist2 + recycle iter が回るボーダーラインの trial が
# timeout 打ち切り → CAPEX 扱いになり、本来 soft fail (生産量未達等) として TPE に方向シグナルが
# 渡るはずの trial を取りこぼす。setdefault なので環境変数で上書き可能。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '120')

# repo root は 1 階層上。
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ===========================================================================
# § 1. 最適化ハイパラ
# ===========================================================================
# N_WORKERS=6 と対。並列 6 worker = async TPE 実効バッチ 6 で僅かに sample 効率が落ちる
# (staleness) のを試行数 +20% で過補償する。QMC startup(=50) は不変。
N_TRIALS    = 360            # Optuna 試行回数
# QMC で広域を網羅して TPE を起動する素材を確保する。特に A_mem 上限近傍 (>2e5) の Sobol 点を
# 十分出さないと TPE が「A_mem 大」を学習素材不足で異常値扱いするため、50 で feasible ~13 件を狙う。
N_STARTUP   = 50             # TPE startup 数
N_TOPK      = 10             # top-k 再評価候補数
SEED        = 42             # 乱数シード (再現用)
SAMPLER     = 'tpe'          # 'tpe' | 'cmaes' | 'random'
# n_jobs=1 デフォルト。penalty_scale が thread-local でないため最終評価は 1 推奨、
# 中間探索は 2-4 で時短可 (SQLite ロック注意)。
N_JOBS      = 1
# マルチプロセス並列 worker 数。>1 で N プロセスが共有 SQLite study を分担 (要 SAVE_SQLITE=True)。
# 各 worker 単スレッドで penalty_scale/GIL 問題なし、constant_liar=True で冗長サンプリング抑制。
# 1 で単一プロセス。N_JOBS(スレッド)とは別物 — 並列は必ず N_WORKERS を使う。
# 品質劣化の唯一原因は「TPE フェーズの同時 worker 数 = async バッチサイズ」= async BO 固有の
# staleness。constant_liar 込みでバッチ 6 は「ほぼ逐次同等」圏の上端 (parallel.py 設計ノート)。
# 8+ は E-core(約半速)に worker が載り straggler 律速 + staleness 増 の二重劣化なので不可。
# QMC startup フェーズは Sobol が結果非依存のため何 worker でも品質ゼロ劣化 (study.py の
# _PhaseSwitchSampler が共有 study の総完了数で切替)。残る TPE staleness は N_TRIALS 側で過補償。
N_WORKERS   = 6


# ===========================================================================
# § 2. ソルバ選択 (各塔独立、'fug' | 'rigorous' | 'sm')
# ===========================================================================
#   - BO ループは速度優先 → 全塔 FUG 推奨 (旧方針)
#   - top-k 再評価は精度優先 → rigorous 推奨
#   - 'sm' は近日実装予定。実装後は文字列を 'sm' に変えるだけで切替可
#
# Dist2 のみ BO ループでも rigorous に切替。
#   FUG path では _split_streams が N_stages を流量計算に使わず、recovery 仕様から直接 split を
#   決定する (= N が CAPEX しか影響しない)。これにより BO が常に N_dist2 下限を選び、partial cond
#   の C3H6 漏れが rigorous で発覚しても BO ループ中は見えない構造になる。
#   Dist2 は N=20-40 / R=5-10 の探索範囲で rigorous でも軽量 (数秒/eval) のため BO ループに
#   組み込んでも全体時間は許容範囲。これで N が物理的に流量に効き、Dist2 厚化の経済合理性が
#   BO に見えるようになる。Dist1/Dist3 は FUG 維持 (Dist1 は narrow-margin で FUG check 機能、
#   Dist3 は N が大きいため rigorous は重い)。
SOLVER_BO   = {
    'dist1': 'fug',
    'dist2': 'rigorous',     # FUG-rigorous gap の根源だったため切替
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
    # 上限 940K: 高 T_in (940-970K) は F×D clip で実 F が大幅に下がり production を維持できず
    #   spec_production_under の主因になる。低 T_in 領域 (880-920K) が黄金帯。
    'T_in_K':            (880.0,  940.0,  'linear', 'float'),  # K
    'z_cat_m':           (15.0,   40.0,   'linear', 'float'),  # m
    't_cyc_min':         (12.0,   25.0,   'linear', 'float'),  # min
    # SV violation は bounds で切らず constraints_func (reactor_sv_shortfall) の信号で学習させる。
    #   D 下限を上げると TPE が小 V_cat 象限 (低 z_cat / 低 t_cyc) と組み合わせて SV>3 を連発し、
    #   D=9.5-10.0 の SV 安全帯という逃げ場を失う。
    'D_reactor_m':       (7.0,    10.0,   'linear', 'float'),  # m

    # ----- PSA -----
    # 下限は既知 feasible 下限 + 余裕に設定。これより小さい領域は PSA 規模不足で
    #   r_psa.t_abs_below_min を連発するため物理的に削除する。
    'D_psa_col_m':       (2.9,    5.0,    'linear', 'float'),  # m
    'L_psa_bed_m':       (22.0,   30.0,   'linear', 'float'),  # m
    'desorption_target': (0.22,   0.40,   'linear', 'float'),  # -

    # ----- 膜 (P_L は 1 atm 固定、P_dist は Dist3 と同期) -----
    'P_H_Pa':            (7.5e5,  9.5e5,  'linear', 'float'),  # Pa
    # A_mem 下限 5e4: これより小面積だと膜 ODE が発散 or 透過量不足で Dist3 starvation 経路に
    #   入り PSA/Mem CAPEX sentinel hit (silent penalty) になる。
    'A_mem_m2':          (5.0e4,  3.0e5,  'log',    'float'),  # m²

    # ----- Dist1 (脱ブタン塔) -----
    # 高 P → α 縮小 → R_min↑ → 与えた R が R_min ぎりぎり → Gilliland N_needed が爆発し
    #   N_dist1 上限で届かない r1 fail 罠。上限 20 bar で右端を切除する。
    'P_dist1_Pa':        (12.0e5, 20.0e5, 'linear', 'float'),  # Pa
    # 下限 20: N_dist1≤18 は Gilliland N_needed 不足で r1 fail。これが低 T_in 帯に集中すると
    #   TPE が「低 T_in = 悪い」と誤学習して黄金帯を避けるため、構造的に削減する。
    'N_dist1':           (20,     30,     'linear', 'int'  ),  # -
    # 下限 2.0: R≤2.0 は r1 fail (Gilliland N 不足) を連発するため低 R を切除する。
    'reflux_dist1':      (2.0,    3.0,    'linear', 'float'),  # -

    # ----- Dist2 (脱エタン塔, partial cond) -----
    # 上限 8e5: Dist2 partial-cond コンデンサ熱量を厳密エネルギー収支にしたことで「塔頂を浅冷化して
    #   安い冷媒に乗せる」のが強い TAC レバー。高 P → 塔頂温度↑ → −100°C エチレン(14373円/GJ) より
    #   暖かい tier に移れる。膜依存 (P_H≥P_dist2+0.5e5, P_H 上限 9.5e5) より 8e5 が上限実用値。
    'P_dist2_Pa':        (5.0e5,  8.0e5,  'linear', 'float'),  # Pa
    # 上限 40: N≥45 領域は Wang-Henke 不収束 (dT_max≈20K) で dist2_dT_shortfall になるため、
    #   物理的に収束不可な領域を閉鎖する。
    'N_dist2':           (20,     40,     'linear', 'int'  ),  # -
    'reflux_dist2':      (6.0,    10.0,   'linear', 'float'),  # -  (Wang-Henke 収束、yield 中立)

    # ----- Dist3 (C3 スプリッタ, narrow-α) -----
    # 上限 19 bar: 高 P_dist3 は α 悪化で N_dist3 が爆増し Dist3 CAPEX が全 CAPEX の ~70% に
    #   膨れる罠なので縛る。mem_bp_shortfall は 16bar 帯でも shortfall シグナルが TPE に渡る。
    'P_dist3_Pa':        (16.0e5, 19.0e5, 'linear', 'float'),  # Pa
    # 下限 90: 低 N + 低 rec_HK 象限は Gilliland N_needed が爆発して N が届かない罠。
    # 上限 150: N=200+ は物理現実 (C3 splitter 標準 N=80-150) と乖離し塔高 60m 超は工業的に
    #   許容不可。上限を縛らないと TPE が「N 大 → 純度安定 → feasible」と N 上限へ逃げる。
    'N_dist3':           (90,    150,    'linear', 'int'  ),  # -
    # 下限 10: R ぎりぎり (R≈R_min×1.07) 設計は N 上限 150 を超える Gilliland N_needed を要して
    #   r3 N不足罠を踏みやすい。下限 10 で Gilliland safe margin (≈R_min×1.18) を確保する。
    #   R↑→ OPEX 増だが N 上限縮小と組合せで Dist3 CAPEX 純減方向。R/N トレードオフは TPE に委ねる。
    'reflux_dist3':      (10.0,   20.0,   'linear', 'float'),  # -

    # ----- Fresh LPG (BO 直接指定、外側ループ skip) -----
    # 1380-1500: 範囲端では生産量が prod_under / prod_over 両側に逸脱して死ぬため、
    #   feasible 生産量帯に対応する F の窓に絞る。
    'F_C3H8_fresh_kmol_h': (1380.0, 1500.0, 'linear', 'float'),  # kmol/h

    # ----- 蒸留塔 recovery -----
    # rec_HK_bot_dist2 には 0.9995 付近に鋭利な feasibility 境界があり、下限を上げて sweet spot を
    #   集中探索する。補完として PSA/Mem trace bypass の閾値超過に連続 penalty (runner.py
    #   _TRACE_BYPASS_PENALTY_COEF_OKUYEN) を与え、Dist2 の C3 漏れ自体に BO 中のシグナルを渡す。
    #   rec_LK_top_dist2 は柔軟にして C2H6 のリサイクル比を経済性で選ばせる。
    # 'rec_LK_top_dist1':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist1':  (0.90, 0.999, 'linear', 'float'),
    'rec_LK_top_dist2':  (0.95, 0.999, 'linear', 'float'),    # C2H6 → top (柔軟)
    # 上限 0.9995: rec_HK_bot_dist2 が高すぎると Dist2 が C3H8 を bot に過剰回収して Dist3 入側 C3
    #   ストリームが numeric 縮退 (0.9999 は浮動小数 epsilon 帯) し r3 (Dist3 feed flow ≤ 0) になる。
    'rec_HK_bot_dist2':  (0.998, 0.9995, 'linear', 'float'),  # C3H8 → bot
    # Dist3 recovery を BO 変数化し、TPE が「purity 99.5% ギリギリ → N_dist3↓, R_dist3↓ → CAPEX↓」を
    #   探索できるようにする (ハードコード 0.99 では over-purity 設計に固定されていた)。
    #   rec_LK_top_dist3 (C3H6 → top): 0.985-0.999、C3H6 の製品側回収率。
    #     下限 0.985 = C3H6 1.5% が bot へ漏れて recycle (production loss) が増える上限。
    'rec_LK_top_dist3':  (0.985, 0.999, 'linear', 'float'),  # - C3H6 → top
    # 下限 0.97: rec_HK_bot_dist3 (C3H8 → bot) = purity 直接制御だが、低 rec × 中 R では Gilliland
    #   N_needed が急増する (R_min は下がるが必要 N が伸びる)。下限 0.97 は C3H6 top 純度
    #   ~99.5-99.8% を達成しつつ N_needed 爆発を回避できる物理的妥協点。
    'rec_HK_bot_dist3':  (0.97,  0.999, 'linear', 'float'),  # - C3H8 → bot
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
# warm-start は無効化 (空リスト)。既知 best 周辺を最初に注入すると:
#   - TPE が anchor された設計近傍に張り付く局所最適化バイアス
#   - 真に広い探索ができてるか不明 (warm-start が無いと出ない解は埋もれる)
#   - bounds + constraints_func + 適切な penalty 係数で十分なはず
#  値を注入する場合は下の cfg.warm_start_trials に list を渡せば動作する。
WARM_START_TRIALS: list = []

# 既知 feasible 候補 (必要時に WARM_START_TRIALS に追加可)
_HISTORICAL_TOP_FEAS = [
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
    # 異なる cluster (T_in 高め)
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
    # 安全装置: 上記 4 件がすべて infeasible でも確実に再現できる候補 (reflux_dist3 は上限 20.0 へ clamp)
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
        n_workers  = N_WORKERS,
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
