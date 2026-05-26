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
# 設計判断 (2026-05-22 forensic, 施策 3b): 50 → 25 に縮小。
# 旧 50 は「QMC で広域 1/6 を網羅 → TPE 起動」の意図だったが、214750 run で 50 trial
# 全部 feas=0 のまま消費 (= QMC は constraints_func を使わない)。25 に下げて TPE を
# 早期起動、QMC で得た shortfall 情報を TPE が早く活用する。
# 設計判断 (2026-05-23 forensic, main_20260523_190747, 施策 R): 25 → 50 に拡大。
# 旧 25 への縮小判断は 214750 時点 (feas 率 0%) の状況、現状は施策 A/B/C/E/F で
# 190747 run feas 率 26% (9/35) に改善済み。QMC 50 trial なら feasible ~13 件期待。
# 真の狙い: QMC 25 trial では A_mem 上限近傍 (>2e5) の Sobol 点が 1-2 件しか出ず、
# TPE が「A_mem 大」を学習素材不足で異常値扱いしていた (190747 best #18 A_mem=2.84e5 を
# TPE が再現できず TAC 1007→1068+ に劣化した構造的原因)。
# 50 trial に拡張で A_mem ∈ [2e5, 3e5] の feasible が 3-5 件期待、TPE 起動時の素材確保。
# 副作用: best 収束が QMC 段階の長期化でやや遅れる (TPE 学習開始が trial 50 から)。
N_STARTUP   = 50             # TPE startup 数 (旧 25→50、A_mem 高領域の Sobol 素材確保)
N_TOPK      = 10             # top-k 再評価候補数
SEED        = 42             # 乱数シード (再現用)
SAMPLER     = 'tpe'          # 'tpe' | 'cmaes' | 'random'
# 設計判断 (2026-05-21): n_jobs=1 デフォルト。並列化したい場合 2-4 (SQLite ロック注意)。
# penalty_scale が thread-local でないため最終評価は 1 推奨、中間探索は 2-4 で時短可。
N_JOBS      = 1
# 設計判断 (2026-05-26): マルチプロセス並列 worker 数。>1 で N プロセスが共有 SQLite study
# を分担 (要 SAVE_SQLITE=True)。各 worker 単スレッドで penalty_scale/GIL 問題なし、
# constant_liar=True で冗長サンプリング抑制。8コアなら 4 が目安 (sample 効率と速度の両立)。
# 1 で従来の単一プロセス。N_JOBS(スレッド)とは別物 — 並列は必ず N_WORKERS を使う。
N_WORKERS   = 4


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
    # 設計判断 (2026-05-23 forensic, main_20260523_085918): 上限 970→940K に縮小。
    # 085918 で TPE が高 T_in (940-970K) に 95% 集中 + F×D clip で実 F=927 (= F_min 1380 を
    # 大幅に下回る) → spec_production_under 50 件の主因に。
    # 高 T_in dead zone を物理的に切って TPE を低 T_in 領域 (880-920K = 黄金帯) へ強制誘導。
    # 黄金帯エビデンス: QMC で T_in<920K の 4 件が TAC 1054-1333 (#12, #26, #30, #5)。
    # 旧 970K 上限は「過去 best 940K 等を含める」目的だったが、085918 forensic で
    # 高 T_in は production 維持できないと判明、上限切り詰めが正解。
    'T_in_K':            (880.0,  940.0,  'linear', 'float'),  # K (上限 970→940, 高 T_in dead zone 排除)
    'z_cat_m':           (15.0,   40.0,   'linear', 'float'),  # m
    't_cyc_min':         (12.0,   25.0,   'linear', 'float'),  # min
    # 設計判断 (2026-05-23 forensic, main_20260523_172800 撤退): bounds 縮小 (7.5, 9.5) は失敗。
    # 試行: 165334 で r_rx=4件 (D≤7.6 / D=9.81) を切るため (7.0,10.0)→(7.5,9.5) に縮小。
    # 結果: 172800 trial 0-44 で r_rx=13件に悪化 (3.25倍)。TPE が新下限 D=7.5-8.2 に張り付き、
    #   小 V_cat 象限 (z_cat=15-22m, t_cyc=12-15min) と組み合わさり SV>3 を連発。
    #   旧 bounds の「D=9.5-10.0 = SV 安全帯」を切ったことで TPE の逃げ場が消えた。
    # 教訓: bounds で切るより constraint (reactor_sv_shortfall) で学習させるべき罠。
    #   F×D SV カップリング廃止 (5/23 085918) と同じく「過剰介入」のパターン。
    # 撤退: 旧 (7.0, 10.0) に巻き戻し、SV violation は constraints_func の信号に委ねる。
    'D_reactor_m':       (7.0,    10.0,   'linear', 'float'),  # m (撤退、SV は constraint 側で学習)

    # ----- PSA -----
    # 設計判断 (2026-05-23 forensic, main_20260523_001237): PSA bounds 縮小。
    # 根拠: 45 trial run の直近 20 trial で r_psa.t_abs_below_min が 12 件 (60%)、
    # TPE が「D_psa≈3.2 + L≈22 + des≈0.23」の局所最適に張り付き、PSA 規模が
    # 物理的に足りない領域から抜けられなかった (#43/#44 でほぼ同じ params 再試行)。
    # feas 2 件は L≥17 + des≥0.22。下限を「既知 feas の下限 + 余裕」に上げ、
    # 規模不足ゾーンを物理的に削除する。
    'D_psa_col_m':       (2.9,    5.0,    'linear', 'float'),  # m (下限 2.5 → 2.9, feas #15 の 2.97 を救う)
    'L_psa_bed_m':       (22.0,   30.0,   'linear', 'float'),  # m (下限 15 → 22, feas #15 の 25.3 内側)
    'desorption_target': (0.22,   0.40,   'linear', 'float'),  # - (下限 0.15 → 0.22, feas 2件は 0.23-0.24)

    # ----- 膜 (P_L は 1 atm 固定、P_dist は Dist3 と同期) -----
    'P_H_Pa':            (7.5e5,  9.5e5,  'linear', 'float'),  # Pa
    # 設計判断 (2026-05-22 E-plan): A_mem 下限を 3e4 → 5e4 に引上げ。
    # 根拠: main_20260522_005631 forensic 解析で A_mem ∈ [3e4, 5e4) の 42 trial が
    #   97.6% (41 件) PSA/Mem CAPEX sentinel hit。膜 ODE が小面積で発散 or 透過量
    #   不足で Dist3 starvation 経路に入って silent penalty 化していた。
    # 履歴 best (TAC=167.74, A_mem=2.84e5) は新範囲内、TAC=641 安全装置 (A_mem=6.39e4) も維持。
    'A_mem_m2':          (5.0e4,  3.0e5,  'log',    'float'),  # m² (124508 best は 2.84e5)

    # ----- Dist1 (脱ブタン塔) -----
    # 設計判断 (2026-05-23 forensic, main_20260523_165334 trial 0-38):
    # r1 fail 12件 (= 全体 31%) は「P_dist1 ∈ [20.5, 25.0] bar × R_dist1 ∈ [1.55, 2.17]」
    # の組合せに全件集中。高 P → α 縮小 → R_min↑ → 与えた R が R_min ぎりぎり →
    # Gilliland N_needed = 26-37 が爆発、N_dist1 上限 30 で届かない罠。
    # 物理的解釈: 「P↑したら R↑も必要」という相関が独立 suggest に反映されない構造。
    # 過去 best (124508 P=20.5/R=2.20, 050049 P=12.6/R=1.58) は P×R 平面の
    # 「低P低R or 中P高R」帯。高P低R 象限は経済的にも非報酬域。
    # 上限 25→20 bar で右端切除、R 下限 1.5→2.0 で低 R 切除して罠 A を物理排除。
    'P_dist1_Pa':        (12.0e5, 20.0e5, 'linear', 'float'),  # Pa (上限 25→20, 高P罠切除)
    # 設計判断 (2026-05-23 forensic, main_20260523_085918): 下限 16→20 に引上げ。
    # r1 fail 11 件中 7 件が N_dist1≤18 (= Gilliland N_needed 不足)。これが低 T_in 帯に
    # 集中したため、TPE は「低 T_in = r1 fail = 悪い」と spurious correlation で誤学習し、
    # 黄金帯 (低 T_in) を 95% 避ける挙動に。N_dist1≥20 にすることで r1 fail を構造的に
    # 削減、TPE の誤学習素材を除去。過去 best (#12 N=29, #272 N=28) は新範囲内。
    'N_dist1':           (20,     30,     'linear', 'int'  ),  # - (下限 16→20, r1 fail 削減)
    # 設計判断 (2026-05-23 forensic, main_20260523_165334): 下限 1.5→2.0 に引上げ。
    # 165334 run で R≤2.0 の r1 fail が連発 (12件中 9件が R<2.0)。R=1.5-2.0 は
    # 過去 best (124508 R=2.20, 050049 R=1.58) のうち 050049 の 1.58 だけ切るが、
    # 050049 系は Dist1 narrow margin で rigorous 再現性が怪しい設計のため切捨て妥協。
    'reflux_dist1':      (2.0,    3.0,    'linear', 'float'),  # -  (下限 1.5→2.0, r1 fail 削減)

    # ----- Dist2 (脱エタン塔, partial cond) -----
    # 設計判断 (2026-05-21 D-plan revert): P_dist2 を b1712d1 (5.0, 7.0)e5 に戻す。
    # 124508 best は P_dist2=5.41e5、5/20 縮小の下限 5.5e5 で切られてた。
    # 設計判断 (2026-05-26): 上限 7e5 → 8e5 に拡大。Dist2 partial-cond コンデンサ熱量を
    # 厳密エネルギー収支に改修 (src/distillation_core.py、+190億の深冷冷凍費) したことで、
    # 「塔頂を浅冷化して安い冷媒に乗せる」のが強い TAC レバーになった。高 P → 塔頂温度↑
    # → −100°C エチレン(14373円/GJ) より暖かい tier に移れる。膜依存 (P_H≥P_dist2+0.5e5,
    # P_H 上限 9.5e5) より P_dist2 上限は 8e5 が上限実用値。BO に浅冷化の余地を与える。
    'P_dist2_Pa':        (5.0e5,  8.0e5,  'linear', 'float'),  # Pa (上限 7→8, 浅冷化レバー)
    # 設計判断 (2026-05-22 forensic 反映、214750 run): N_dist2 上限を 50 → 40 に戻す。
    # 根拠: 5/22 改良 3 で 40→50 拡張したが、main_20260522_214750 (n=130 stopped)
    # で TPE phase 80 trial の 45% が N_dist2 ≥ 45 領域に張り付き、dist2_dT_shortfall
    # NZ 14 件全てが N=45-50 で Wang-Henke 不収束 (dT_max≈20K)。物理的に収束不可な
    # 領域に TPE が誤誘導された罠。上限 40 に戻して罠を物理的に閉鎖。
    # 124508 best (N=22) は新範囲内、下限 20 維持。
    'N_dist2':           (20,     40,     'linear', 'int'  ),  # - (124508 best は 22)
    'reflux_dist2':      (6.0,    10.0,   'linear', 'float'),  # -  維持 (Wang-Henke 収束、yield 中立)

    # ----- Dist3 (C3 スプリッタ, narrow-α) -----
    # 設計判断 (2026-05-23 forensic, main_20260523_004243): P_dist3 上限 25→19 + 下限 17→16。
    # 根拠: 004243 run で TPE が P_dist3 中央値 22-23 bar に向かい、α 悪化で N_dist3 が
    # 171 まで爆増、Dist3 CAPEX 414 億円 (全 CAPEX の 70%) になった。
    # 過去 best (TAC=167.74, P_dist3=16.72bar, N_dist3=93) を含めるため下限 17→16、
    # かつ「TPE が高 P_dist3 罠に再び向かう」のを防ぐため上限 25→19 に縛る。
    # 補完: mem_bp_shortfall 対策は 16bar 帯でも shortfall シグナルが TPE に渡るので
    # 自動学習を期待 (P_dist3 を 16→17 bar 上げると Mem が楽になる方向シグナルは残る)。
    'P_dist3_Pa':        (16.0e5, 19.0e5, 'linear', 'float'),  # Pa (過去 best 16.72 を含む)
    # 設計判断 (2026-05-21 D-plan): N_dist3 下限を 80 に。124508 best (TAC=167.74) は N=93。
    # 設計判断 (2026-05-23 forensic, main_20260523_172800): 下限 80→90 に引上げ。
    # 根拠: 172800 run trial #26 で N=80 + rec_HK_bot_dist3=0.982 (下限近傍) の組合せで
    # Gilliland N_needed=103 が爆発、N=80 では届かない罠を踏んだ。
    # rec_HK_bot_dist3 変数化 (5/23 004243) の副作用で「N 下限 + 低 rec」象限に
    # 新たな N不足罠が生まれていた。下限を 124508 best (N=93) より少し下 (90) まで
    # 詰めることで罠象限を排除。過去 best 93 は内側、TAC=247(050049) 系の N=132 等も影響なし。
    # 設計判断 (2026-05-23 forensic, main_20260523_182527): 上限 250→150 に縮小。
    # 根拠: 182527 run TPE phase (#26-29) が N_dist3=246-249 (上限張付き) で TAC=1156-1196 に
    # 集中、同 run の QMC best #18 (N=135, TAC=1011.78) を捨てて N 上限方向に学習。
    # 「N 大 → 純度安定 → feasible」を施策 A (spec ±5% 緩和) と組み合わせて強く学習した結果。
    # 過去 best (124508 N=93, 050049 N=132, 182527 #18 N=135) は全て N≤135 で、N=200+ は
    # 物理現実 (C3 splitter 標準 N=80-150) と乖離。塔高 60m 超は工業的に許容不可。
    # 上限 150 で過去 best 全部内側、TPE の「N 上限逃げ」経路を物理的に閉鎖。
    # 副作用: N>150 必要な低 P + 低 rec 設計は r3 (Gilliland N不足) で infeas 化。
    # これは施策 F (reflux_dist3 下限 ↑) で軽減。
    'N_dist3':           (90,    150,    'linear', 'int'  ),  # - (上限 250→150, N 上限逃げ閉鎖)
    # 設計判断 (2026-05-23 forensic, main_20260523_004243): 下限 11→9 に拡張。
    # 根拠: 004243 best (#272) は R=11.17 で R_min=8.46 → margin 32%。下限 11 だと
    # 「R ぎりぎり」の小型 Dist3 設計を TPE が試せない構造になっていた。
    # 9.0 ≈ R_min×1.07 まで下げれば、TPE が Gilliland の物理境界近傍を探索可能。
    # 副作用は Wang-Henke 不収束 (= dist3_dT_shortfall) だが信号化済みで TPE 学習可能。
    # 期待効果: Dist3 N↓ → 全 CAPEX 414→200 億円台。
    # 設計判断 (2026-05-23 forensic, main_20260523_182527): 下限 9.0→10.0 に引上げ。
    # 根拠: N_dist3 上限を 250→150 に縮小したため、R ぎりぎり (R≈R_min×1.07) 設計で
    # Gilliland N_needed が 150 を超えると r3 N不足罠を踏みやすい。
    # 過去 best (124508 R=11.50, 050049 R=11.52, 182527 #18 R=10.03) は全て R≥10、
    # 下限 10.0 で過去 best 全部内側、かつ Gilliland safe margin (≈R_min×1.18) を確保。
    # 副作用: R↑→ OPEX 増 (リボイラー熱負荷↑) だが、N 上限縮小と組合せで Dist3 CAPEX
    # 純減方向 (= TAC 純減方向)。R/N トレードオフは TPE に委ねる。
    'reflux_dist3':      (10.0,   20.0,   'linear', 'float'),  # - (下限 9→10, N=150 上限と整合)

    # ----- Fresh LPG (BO 直接指定、外側ループ skip) -----
    # 設計判断 (2026-05-22 forensic, 214750 run): 範囲を 1200-1700 → 1380-1500 に縮小。
    # 根拠: 214750 run で prod_under (中央 F=1620) と prod_over (中央 F=1480) が二極化、
    # 21 件 (16%) が F_fresh 範囲端での生産量逸脱で死亡。1380-1500 は履歴 best
    # (F=1424.1, 1414.1, 1433.2) を全部含み、prod_over/under 両側の dead zone を切る。
    # 5/22 forensic 解析の F_fresh 別ヒストグラム支持。
    'F_C3H8_fresh_kmol_h': (1380.0, 1500.0, 'linear', 'float'),  # kmol/h

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
    # 設計判断 (2026-05-23 forensic, main_20260523_165334): 上限 0.9999→0.9995 に縮小。
    # r3 (Dist3 feed flow ≤ 0) が 4件 (10%) 発生、rec_HK_bot_dist2 が高い設計で Dist2 が
    # C3H8 を bot に過剰回収 → Dist3 入側 C3 ストリームが numeric 縮退。
    # 上限 0.9999 は実装上 numeric 縮退を引き起こしやすい (浮動小数 epsilon 帯)。
    # 過去 best (124508 0.9997) は新 bounds 内、TAC=247(050049) 系も 0.9995 系。
    'rec_HK_bot_dist2':  (0.998, 0.9995, 'linear', 'float'),  # C3H8 → bot (上限 0.9999→0.9995, r3 fail 削減)
    # 設計判断 (2026-05-23 forensic, main_20260523_004243): Dist3 recovery を BO 変数化。
    # 根拠: 004243 best (#272 rigorous, TAC=1055) は purity 99.88wt% (spec 99.5%+0.38pp over-design)
    # で Dist3 CAPEX が 414 億円 (全 CAPEX の 69.6%)。ハードコード 0.99 で TPE が
    # 「purity ぎりぎり 99.5%」の Dist3 縮小設計を探れない構造になっていた。
    # 変数化することで TPE が「purity 99.5% ギリギリ → N_dist3↓, R_dist3↓ → CAPEX↓」を
    # 探索できる。期待効果: Dist3 CAPEX 414 → 200 億円台、TAC 1055 → 800 程度。
    # 範囲:
    #   rec_LK_top_dist3 (C3H6 → top): 0.985-0.999、C3H6 の製品側回収率
    #     - 下限 0.985: C3H6 1.5% が bot へ漏れて recycle = production loss が増える上限
    #     - 上限 0.999: 現状ハードコード相当
    #   rec_HK_bot_dist3 (C3H8 → bot): 0.95-0.999、C3H8 の bot 側回収率 = purity 直接制御
    #     - 下限 0.95: C3H8 5% が top へ漏れて purity ~95% = spec 違反確定
    #     - 上限 0.999: 現状ハードコード相当 (= over-purity 設計)
    #     - sweet spot 推定: 0.99 付近 (= purity 99.5% spec ギリギリ)
    'rec_LK_top_dist3':  (0.985, 0.999, 'linear', 'float'),  # - C3H6 → top
    # 設計判断 (2026-05-23 forensic, main_20260523_172800): 下限 0.95→0.97 に引上げ。
    # 根拠: 172800 run trial #8 (rec_HK=0.981, N=101, R=10.4 で N_needed=121) と trial #26
    # (rec_HK=0.982, N=80, R=12.5 で N_needed=103) で Dist3 N_needed 爆発罠を観測。
    # rec_HK 低下は purity 制御として有用だが、低 rec × 中 R では Gilliland N_needed が
    # 急増 (Underwood R_min 自体は低 rec で下がるものの、Gilliland 必要 N が伸びる構造)。
    # 下限 0.97 は purity ~97% bot で C3H6 top 純度 ~99.5-99.8% を達成できる帯、
    # かつ N_needed 爆発を回避できる物理的妥協点。
    # 副作用: 「purity ぎりぎり 99.5%」狙いの極小型 Dist3 設計 (rec_HK=0.95-0.97 帯) を
    # 探索から外す。期待効果は元の Dist3 CAPEX 削減狙い (004243 forensic) より小さくなるが、
    # N不足罠を許容するよりは妥当。
    'rec_HK_bot_dist3':  (0.97,  0.999, 'linear', 'float'),  # - C3H8 → bot (下限 0.95→0.97, N不足罠切除)
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
