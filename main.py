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
    't_cyc_min':         (10.0,   30.0,   'linear', 'float'),  # min !仮置き (触媒再生 30 min との比から経験範囲)
    'D_reactor_m':       (4.0,    10.0,   'linear', 'float'),  # m  !仮置き (大型固定床の経験範囲)

    # ----- PSA -----
    'D_psa_col_m':       (1.5,    5.0,    'linear', 'float'),  # m  制約: 空塔速度 ≤ 1 m/s (化工便覧 §13-31, 但し除湿用)
    'L_psa_bed_m':       (10.0,   30.0,   'linear', 'float'),  # m  !仮置き
    'desorption_target': (0.15,   0.55,   'linear', 'float'),  # -  !仮置き

    # ----- 膜 (P_L は 1 atm 固定、P_dist は Dist3 と同期) -----
    'P_H_Pa':            (5.0e5,  9.5e5,  'linear', 'float'),  # Pa 上限: Hua et al. (2024) 9.5 bar / 下限: !仮置き
    'A_mem_m2':          (3.0e4,  3.0e5,  'log',    'float'),  # m² !仮置き (CAPEX 支配、log scale)

    # ----- Dist1 (脱ブタン塔) -----
    'P_dist1_Pa':        (12.0e5, 25.0e5, 'linear', 'float'),  # Pa !仮置き (pump1 出口圧と同期)
    'N_dist1':           (16,     30,     'linear', 'int'  ),  # -  下限引き上げ (旧14→16): N_min ≈ 12 から margin 33% 確保
    'reflux_dist1':      (1.1,    3.0,    'linear', 'float'),  # -  下限引き上げ (旧1.0→1.1): R_min ≈ 0.95 から margin 16% 確保

    # ----- Dist2 (脱エタン塔, partial cond) -----
    'P_dist2_Pa':        (5.0e5,  9.5e5,  'linear', 'float'),  # Pa 上限: P_H 9.5 bar と整合 / 下限: !仮置き
    'N_dist2':           (20,     40,     'linear', 'int'  ),  # -  下限引き上げ (旧10→20): rigorous で 99% recovery 物理達成可能な領域に
    'reflux_dist2':      (3.0,    10.0,   'linear', 'float'),  # -  R_min 1.5-4.8 運転依存、下限引き上げ (旧2→3): narrow-margin 回避

    # ----- Dist3 (C3 スプリッタ, narrow-α) -----
    'P_dist3_Pa':        (15.0e5, 25.0e5, 'linear', 'float'),  # Pa !仮置き (mem.P_dist と同期、冷却水凝縮可能下限近傍)
    'N_dist3':           (100,    250,    'linear', 'int'  ),  # -  下限: N_min ≈ 81 + margin / 上限: !仮置き
    'reflux_dist3':      (11.0,   20.0,   'linear', 'float'),  # -  下限: R_min ≈ 10 + margin / 上限: !仮置き

    # ----- Fresh LPG (BO 直接指定、外側ループ skip) -----
    # !仮置き: baseline exp1 で外側ループ実測 F_fresh ≈ 1666 kmol/h を中心に -50% / +50% 程度。
    # production_min spec (1188 ± 1% kmol/h) は soft penalty として効くので、
    # F_fresh が低すぎる = 生産量不足 = penalty で自然に上向きに収束する想定。
    'F_C3H8_fresh_kmol_h': (800.0, 2500.0, 'linear', 'float'),  # kmol/h !仮置き

    # ----- 蒸留塔 recovery は最適化対象外 (0.99 固定、変更したい場合は ColumnTunables 経由で指定可) -----
    # 'rec_LK_top_dist1':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist1':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_LK_top_dist2':  (0.90, 0.999, 'linear', 'float'),
    # 'rec_HK_bot_dist2':  (0.90, 0.999, 'linear', 'float'),
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

# ----- L1: Feasibility 分類解析 (BO 終了後の post-hoc 解析) -----
RUN_FEASIBILITY_ANALYSIS = True       # False で無効化 (sklearn 未インストール時も自動無効)
FEASIBILITY_TARGET       = 'convergence'  # 'convergence' | 'spec' | 'both'
FEASIBILITY_MODEL        = 'rf'       # 'rf' | 'logreg'


# ===========================================================================
# ↑↑↑ 編集領域はここまで。以下はパイプライン本体 (通常触らない) ↑↑↑
# ===========================================================================

from datetime import datetime
from pathlib import Path

from config.load import load_operating_config
from optimization import (
    validate_search_space,
    make_objective,
    create_study,
    run_optimization,
    reevaluate_topk,
    best_entry,
    save_trials_csv,
    save_best_json,
    save_topk_report,
)
try:
    from optimization import analyze_feasibility
    _HAS_FEASIBILITY = True
except ImportError:
    _HAS_FEASIBILITY = False


def main():
    # ---- 入力検査 ----
    validate_search_space(SEARCH_SPACE)
    for tag in ('dist1', 'dist2', 'dist3'):
        if tag not in SOLVER_BO or tag not in SOLVER_TOPK:
            raise ValueError(f"SOLVER_BO / SOLVER_TOPK に {tag} のキーが必要")
        for s in (SOLVER_BO[tag], SOLVER_TOPK[tag]):
            if s not in ('fug', 'rigorous', 'sm'):
                raise ValueError(f"solver {s!r} は 'fug' | 'rigorous' | 'sm' のみ許容")

    # ---- パス準備 ----
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir    = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    study_name = f'pdh_{timestamp}'
    db_path        = out_dir / f'main_{timestamp}.db'
    trials_csv     = out_dir / f'main_{timestamp}_trials.csv'
    best_json      = out_dir / f'main_{timestamp}_best.json'
    topk_report    = out_dir / f'main_{timestamp}_topk.txt'
    feasibility_prefix = f'main_{timestamp}_feasibility'

    storage_url = f'sqlite:///{db_path.as_posix()}' if SAVE_SQLITE else None

    # ---- 設定スナップショット ----
    print("=" * 72)
    print(f"PDH 多変数最適化 — {timestamp}")
    print("=" * 72)
    print(f"  N_TRIALS         = {N_TRIALS}")
    print(f"  N_STARTUP        = {N_STARTUP}")
    print(f"  N_TOPK           = {N_TOPK}")
    print(f"  SAMPLER          = {SAMPLER}")
    print(f"  SEED             = {SEED}")
    print(f"  SOLVER_BO        = {SOLVER_BO}")
    print(f"  SOLVER_TOPK      = {SOLVER_TOPK}")
    print(f"  APPLY_HI         = {APPLY_HI}")
    print(f"  APPLY_STAGE2_TOPK= {APPLY_STAGE2_TOPK}")
    print(f"  探索変数数        = {len(SEARCH_SPACE)} / 19")
    print(f"  storage          = {storage_url or '(in-memory)'}")
    print("-" * 72)

    # ---- config / objective / study 準備 ----
    config = load_operating_config()
    objective = make_objective(
        search_space          = SEARCH_SPACE,
        solver_assignment     = SOLVER_BO,
        config                = config,
        apply_hi              = APPLY_HI,
        apply_stage2          = False,                # BO ループでは Stage 2 を回さない
        hi_dT_min_K           = HI_DT_MIN_K,
        strict_recovery_check = STRICT_RECOVERY_BO,
        recovery_tolerance    = RECOVERY_TOLERANCE,
    )
    study = create_study(
        study_name   = study_name,
        sampler_name = SAMPLER,
        seed         = SEED,
        n_startup    = N_STARTUP,
        storage_url  = storage_url,
    )

    # ---- BO ループ実行 (KeyboardInterrupt / 致命的例外でも部分結果を保存) ----
    print(f"[BO] {N_TRIALS} trial を実行中 ...")
    t_start = datetime.now()
    bo_interrupted  = False
    bo_fatal_error  = None
    try:
        run_optimization(
            study, objective,
            n_trials          = N_TRIALS,
            show_progress_bar = SHOW_PROGRESS,
        )
    except KeyboardInterrupt:
        bo_interrupted = True
        print("\n[BO] Ctrl+C で中断。これまでの結果を保存して終了します。")
    except Exception as e:
        bo_fatal_error = e
        print(f"\n[BO] 致命的例外 ({type(e).__name__}: {e})。"
              f" SQLite に保存された分は再開可能 (study_name='{study_name}')。")
    t_bo = (datetime.now() - t_start).total_seconds()
    n_completed = sum(1 for t in study.trials
                      if t.state.name == 'COMPLETE')
    print(f"[BO] 経過 {t_bo:.1f} 秒、完了 trial = {n_completed} / 試行 {len(study.trials)}")

    try:
        best = study.best_trial
        print(f"[BO] ベスト trial #{best.number}: "
              f"effective_TAC = {best.value:.4f} 億円/年")
    except ValueError:
        print("[BO] 完了 trial がありません (全失敗)")

    # ---- top-k 再評価 (BO が中断/失敗してもこれまでの trial で再評価試行) ----
    entries = []
    if N_TOPK > 0 and n_completed > 0 and bo_fatal_error is None:
        print(f"[top-k] 上位 {min(N_TOPK, n_completed)} 候補を再評価中 "
              f"(solver={SOLVER_TOPK}, stage2={APPLY_STAGE2_TOPK}) ...")
        t_start = datetime.now()
        try:
            entries = reevaluate_topk(
                study              = study,
                k                  = N_TOPK,
                solver_assignment  = SOLVER_TOPK,
                config             = config,
                apply_hi           = APPLY_HI,
                apply_stage2       = APPLY_STAGE2_TOPK,
                hi_dT_min_K        = HI_DT_MIN_K,
                strict_recovery_check = STRICT_RECOVERY_TOPK,
                recovery_tolerance = RECOVERY_TOLERANCE,
                verbose            = False,
            )
        except KeyboardInterrupt:
            print("\n[top-k] Ctrl+C で中断。これまでに完了した再評価分のみ保存します。")
        except Exception as e:
            print(f"\n[top-k] 再評価で例外 ({type(e).__name__}: {e})。"
                  f" 部分結果を保存します。")
        t_topk = (datetime.now() - t_start).total_seconds()
        print(f"[top-k] 経過 {t_topk:.1f} 秒、完了 {len(entries)} 候補")
        be = best_entry(entries)
        if be is not None:
            tag = "feasible" if be.is_feasible_re else "infeasible"
            print(f"[top-k] 再評価ベスト: rank {be.rank} (trial #{be.trial_number}), "
                  f"effective_TAC = {be.effective_TAC_re:.4f} 億円/年 ({tag})")
    elif N_TOPK > 0 and n_completed == 0:
        print("[top-k] BO 完了 trial が 0 件、再評価スキップ。")
    elif bo_fatal_error is not None:
        print("[top-k] BO で致命的例外発生、再評価スキップ。")

    # ---- 出力 (例外が出ても可能な限り部分結果を残す) ----
    print("[出力] 保存中 ...")
    if SAVE_TRIALS_CSV:
        try:
            n = save_trials_csv(study, trials_csv)
            print(f"[出力] trials CSV → {trials_csv} ({n} 行)")
        except Exception as e:
            print(f"[出力] trials CSV 失敗: {type(e).__name__}: {e}")
    if SAVE_BEST_JSON:
        try:
            save_best_json(study, best_json)
            print(f"[出力] best JSON → {best_json}")
        except Exception as e:
            print(f"[出力] best JSON 失敗: {type(e).__name__}: {e}")
    if SAVE_TOPK_REPORT and entries:
        try:
            save_topk_report(entries, topk_report)
            print(f"[出力] top-k 比較 → {topk_report}")
        except Exception as e:
            print(f"[出力] top-k report 失敗: {type(e).__name__}: {e}")

    # ---- L1: Feasibility 分類解析 (post-hoc、副作用ゼロ) ----
    if RUN_FEASIBILITY_ANALYSIS and _HAS_FEASIBILITY and n_completed > 0:
        try:
            print("[feasibility] 分類解析を実行中 ...")
            analyze_feasibility(
                study       = study,
                output_dir  = out_dir,
                prefix      = feasibility_prefix,
                target_type = FEASIBILITY_TARGET,
                model       = FEASIBILITY_MODEL,
            )
        except Exception as e:
            print(f"[feasibility] 解析失敗: {type(e).__name__}: {e}")
    elif RUN_FEASIBILITY_ANALYSIS and not _HAS_FEASIBILITY:
        print("[feasibility] スキップ: scikit-learn が未インストール")

    print("=" * 72)
    if bo_interrupted:
        print(f"中断完了。再開コマンド: 同じ main.py を再実行 (storage URL が同じ)")
        print(f"  → storage: {storage_url}")
    elif bo_fatal_error is not None:
        print(f"異常終了。エラー: {type(bo_fatal_error).__name__}: {bo_fatal_error}")
        raise bo_fatal_error
    else:
        print("完了。")


if __name__ == '__main__':
    main()
