"""
optimization/pipeline.py — BO + top-k + L1 + 詳細結果表示の一括オーケストレーション

main.py の main() を関数化し、設定 dataclass で全パラメータを受ける形に。
top-k 完了後、ベスト候補について exp1 と同じ詳細レポート (show_input_snapshot +
display_full_results) を出力する機能も追加。

設計判断:
  - main.py は § 1-5 の編集領域 + run_pipeline(...) 呼び出し 1 行のみに集中
  - 例外/KeyboardInterrupt の救出ロジックは pipeline 内に集約
  - 詳細表示は top-k 再評価の result を再利用 (= 追加計算なし)
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Any

from config.load import load_operating_config
from optimization.search_space import (
    VarSpec, validate_search_space, build_design,
)
from optimization.objective import make_objective
from optimization.study import create_study, run_optimization
from optimization.topk import reevaluate_topk, best_entry, TopKEntry
from optimization.reporting import (
    save_trials_csv, save_best_json, save_topk_report,
)
from simulation import (
    display_full_results, hdr, show_input_snapshot,
)


# ---------------------------------------------------------------------------
# 設定 dataclass
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    """`run_pipeline()` に渡す全ハイパラ。main.py の § 1-5 を 1 オブジェクトに束ねる。"""

    # § 1. 最適化ハイパラ
    n_trials:    int             = 300
    n_startup:   int             = 50
    n_topk:      int             = 10
    seed:        int             = 42
    sampler:     str             = 'tpe'

    # § 2. ソルバ選択
    solver_bo:   Dict[str, str]  = field(default_factory=lambda: {
        'dist1': 'fug', 'dist2': 'fug', 'dist3': 'fug',
    })
    solver_topk: Dict[str, str]  = field(default_factory=lambda: {
        'dist1': 'rigorous', 'dist2': 'rigorous', 'dist3': 'rigorous',
    })

    # § 3. 評価オプション
    apply_hi:              bool  = True
    apply_stage2_topk:     bool  = True
    hi_dT_min_K:           float = 10.0
    strict_recovery_bo:    bool  = False
    strict_recovery_topk:  bool  = True
    recovery_tolerance:    float = 0.10

    # § 4. 探索空間
    search_space:          Dict[str, VarSpec] = field(default_factory=dict)

    # § 5. 出力
    output_dir:            str   = 'outputs'
    save_sqlite:           bool  = True
    save_trials_csv:       bool  = True
    save_best_json:        bool  = True
    save_topk_report:      bool  = True
    show_progress:         bool  = True

    # L1: feasibility 分類解析
    run_feasibility_analysis:  bool = True
    feasibility_target:        str  = 'convergence'
    feasibility_model:         str  = 'rf'

    # ベスト候補の詳細表示 (exp1 と同じ display_full_results)
    display_best_full:         bool = True


# ---------------------------------------------------------------------------
# 入力検査
# ---------------------------------------------------------------------------

def _validate(cfg: PipelineConfig) -> None:
    """SEARCH_SPACE と SOLVER_* を検査。問題あれば ValueError。"""
    validate_search_space(cfg.search_space)
    for tag in ('dist1', 'dist2', 'dist3'):
        if tag not in cfg.solver_bo or tag not in cfg.solver_topk:
            raise ValueError(f"solver_bo / solver_topk に {tag} のキーが必要")
        for s in (cfg.solver_bo[tag], cfg.solver_topk[tag]):
            if s not in ('fug', 'rigorous', 'sm'):
                raise ValueError(f"solver {s!r} は 'fug' | 'rigorous' | 'sm' のみ許容")


# ---------------------------------------------------------------------------
# 設定スナップショット
# ---------------------------------------------------------------------------

def _print_snapshot(cfg: PipelineConfig, timestamp: str, storage_url: Optional[str]) -> None:
    print("=" * 72)
    print(f"PDH 多変数最適化 — {timestamp}")
    print("=" * 72)
    print(f"  N_TRIALS         = {cfg.n_trials}")
    print(f"  N_STARTUP        = {cfg.n_startup}")
    print(f"  N_TOPK           = {cfg.n_topk}")
    print(f"  SAMPLER          = {cfg.sampler}")
    print(f"  SEED             = {cfg.seed}")
    print(f"  SOLVER_BO        = {cfg.solver_bo}")
    print(f"  SOLVER_TOPK      = {cfg.solver_topk}")
    print(f"  APPLY_HI         = {cfg.apply_hi}")
    print(f"  APPLY_STAGE2_TOPK= {cfg.apply_stage2_topk}")
    print(f"  探索変数数        = {len(cfg.search_space)} (うち最大 19)")
    print(f"  storage          = {storage_url or '(in-memory)'}")
    print("-" * 72)


# ---------------------------------------------------------------------------
# ベスト候補の詳細表示
# ---------------------------------------------------------------------------

def _write_readme(out_dir: Path, cfg, paths, study, be, timestamp: str) -> None:
    """run subdir に README.md を出力 (結果の見方ガイド)。"""
    lines = []
    lines.append(f"# PDH 多変数最適化 run — {timestamp}")
    lines.append("")
    lines.append("## まず見るべきファイル (推奨順)")
    lines.append("")
    lines.append("1. **`topk.txt`** ─ ★最終結果。top-k 候補の rigorous + Stage 2 再評価詳細。")
    lines.append("   - BO ベスト ≠ 真のベストの場合あり、本ファイルの `rank 1` が実際の最良。")
    lines.append("   - `feas_re=True` の中で `effective_TAC_re` 最小が「採用すべき設計」。")
    lines.append("2. **`best.json`** ─ BO 単体ベスト trial の params + 経済値。簡易確認用。")
    lines.append("3. **`feasibility.txt`** ─ 収束分類器の AUC + 特徴量重要度。")
    lines.append("   - 「どの設計変数が収束/非収束を左右するか」のヒント。")
    lines.append("")
    lines.append("## 詳細解析用")
    lines.append("")
    lines.append("4. **`trials.csv`** ─ 全 trial 履歴。Excel/pandas で散布図・統計解析。")
    lines.append("5. **`feasibility_2d.png`** ─ 最重要 2 変数で feasible/infeasible 散布図。")
    lines.append("6. **`optuna.db`** ─ Optuna SQLite。")
    lines.append("   - 中断時の再開: 同じ main.py で再実行 (storage 自動検出)。")
    lines.append("   - 可視化: `optuna-dashboard sqlite:///optuna.db`")
    lines.append("")
    lines.append("## この run の設定")
    lines.append("")
    lines.append(f"- N_TRIALS = {cfg.n_trials}, N_STARTUP = {cfg.n_startup}, N_TOPK = {cfg.n_topk}")
    lines.append(f"- SAMPLER = {cfg.sampler}, SEED = {cfg.seed}")
    lines.append(f"- SOLVER_BO = {cfg.solver_bo}")
    lines.append(f"- SOLVER_TOPK = {cfg.solver_topk}")
    lines.append(f"- 探索変数数 = {len(cfg.search_space)}")
    lines.append("")
    lines.append("## ベスト要約")
    lines.append("")
    if be is not None:
        tag = "feasible ✓" if be.is_feasible_re else "infeasible ✗"
        lines.append(f"- top-k 再評価ベスト: **trial #{be.trial_number}** (rank {be.rank})")
        lines.append(f"- effective_TAC (再評価) = **{be.effective_TAC_re:.3f}** 億円/年 ({tag})")
        if be.result.economics_synth is not None:
            lines.append(f"- Profit (Stage 2 後) = **{be.result.economics_synth.profit:+.3f}** 億円/年")
    else:
        try:
            best = study.best_trial
            lines.append(f"- BO ベスト trial #{best.number}: effective_TAC = {best.value:.3f} 億円/年")
            lines.append(f"- top-k 再評価未実施 or 全 infeasible")
        except ValueError:
            lines.append(f"- 完了 trial なし")
    lines.append("")
    (out_dir / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _display_best_full(entry: TopKEntry, cfg: PipelineConfig, config) -> None:
    """top-k のベスト候補について exp1 と同じ詳細レポートを出力。

    既存の re-eval 結果 (entry.result) を再利用するので追加計算なし。
    show_input_snapshot + display_full_results を呼ぶ。
    """
    print()
    hdr(f"top-k ベスト候補の詳細 (trial #{entry.trial_number}, rank {entry.rank})")

    # design を rebuild (params + SOLVER_TOPK で entry.result と整合)
    design = build_design(entry.params, cfg.solver_topk)
    eval_kwargs = {
        'apply_hi':     cfg.apply_hi,
        'apply_stage2': cfg.apply_stage2_topk,
        'hi_dT_min_K':  cfg.hi_dT_min_K,
    }
    show_input_snapshot(design, config, eval_kwargs)
    display_full_results(entry.result, design, config)


# ---------------------------------------------------------------------------
# パイプライン本体
# ---------------------------------------------------------------------------

def run_pipeline(cfg: PipelineConfig) -> dict:
    """BO + top-k + L1 + 詳細表示を一括実行。

    Returns
    -------
    dict
      {'study': ..., 'topk_entries': [...], 'best_entry': ...,
       'paths': {db, trials_csv, best_json, topk_report, feasibility_prefix}}
    """
    _validate(cfg)

    # ---- パス準備 ----
    # 設計判断 (2026-05-17): 各 run の成果物 6 ファイルを 1 subdir にまとめる。
    # outputs/main_<timestamp>/{optuna.db, trials.csv, best.json, topk.txt,
    #   feasibility.txt, feasibility_2d.png}
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_root   = Path(cfg.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    out_dir    = out_root / f'main_{timestamp}'   # 本 run 専用 subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    study_name = f'pdh_{timestamp}'
    paths = {
        'db':                  out_dir / 'optuna.db',
        'trials_csv':          out_dir / 'trials.csv',
        'best_json':           out_dir / 'best.json',
        'topk_report':         out_dir / 'topk.txt',
        'feasibility_prefix':  'feasibility',   # → out_dir/feasibility.txt, feasibility_2d.png
    }
    storage_url = f'sqlite:///{paths["db"].as_posix()}' if cfg.save_sqlite else None

    _print_snapshot(cfg, timestamp, storage_url)

    # ---- config / objective / study 準備 ----
    config = load_operating_config()
    objective = make_objective(
        search_space          = cfg.search_space,
        solver_assignment     = cfg.solver_bo,
        config                = config,
        apply_hi              = cfg.apply_hi,
        apply_stage2          = False,                # BO ループでは Stage 2 を回さない
        hi_dT_min_K           = cfg.hi_dT_min_K,
        strict_recovery_check = cfg.strict_recovery_bo,
        recovery_tolerance    = cfg.recovery_tolerance,
    )
    study = create_study(
        study_name   = study_name,
        sampler_name = cfg.sampler,
        seed         = cfg.seed,
        n_startup    = cfg.n_startup,
        storage_url  = storage_url,
    )

    # ---- BO ループ実行 (KeyboardInterrupt / 致命的例外でも部分結果を保存) ----
    # 設計判断 (2026-05-20): Optuna 標準 logger を WARNING に絞り、自前 compact callback で
    # 1 trial = 構造化 5 行表示 (status + vars 3 行 + progress/ETA) に置換。可読性向上。
    import optuna as _optuna_for_log
    _optuna_for_log.logging.set_verbosity(_optuna_for_log.logging.WARNING)
    from optimization.callbacks import make_compact_callback
    _compact_cb = make_compact_callback(n_trials_total=cfg.n_trials)

    print(f"[BO] {cfg.n_trials} trial を実行中 ...")
    t_start = datetime.now()
    bo_interrupted = False
    bo_fatal_error = None
    try:
        run_optimization(
            study, objective,
            n_trials          = cfg.n_trials,
            show_progress_bar = False,           # tqdm は自前 ETA と競合するため無効
            callbacks         = [_compact_cb],
        )
    except KeyboardInterrupt:
        bo_interrupted = True
        print("\n[BO] Ctrl+C で中断。これまでの結果を保存して終了します。")
    except Exception as e:
        bo_fatal_error = e
        print(f"\n[BO] 致命的例外 ({type(e).__name__}: {e})。"
              f" SQLite に保存された分は再開可能 (study_name='{study_name}')。")
    t_bo = (datetime.now() - t_start).total_seconds()
    n_completed = sum(1 for t in study.trials if t.state.name == 'COMPLETE')
    print(f"[BO] 経過 {t_bo:.1f} 秒、完了 trial = {n_completed} / 試行 {len(study.trials)}")

    try:
        best = study.best_trial
        print(f"[BO] ベスト trial #{best.number}: "
              f"effective_TAC = {best.value:.4f} 億円/年")
    except ValueError:
        print("[BO] 完了 trial がありません (全失敗)")

    # ---- top-k 再評価 ----
    entries = []
    if cfg.n_topk > 0 and n_completed > 0 and bo_fatal_error is None:
        print(f"[top-k] 上位 {min(cfg.n_topk, n_completed)} 候補を再評価中 "
              f"(solver={cfg.solver_topk}, stage2={cfg.apply_stage2_topk}) ...")
        t_start = datetime.now()
        try:
            entries = reevaluate_topk(
                study              = study,
                k                  = cfg.n_topk,
                solver_assignment  = cfg.solver_topk,
                config             = config,
                apply_hi           = cfg.apply_hi,
                apply_stage2       = cfg.apply_stage2_topk,
                hi_dT_min_K        = cfg.hi_dT_min_K,
                strict_recovery_check = cfg.strict_recovery_topk,
                recovery_tolerance = cfg.recovery_tolerance,
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
    elif cfg.n_topk > 0 and n_completed == 0:
        print("[top-k] BO 完了 trial が 0 件、再評価スキップ。")
    elif bo_fatal_error is not None:
        print("[top-k] BO で致命的例外発生、再評価スキップ。")

    # ---- 出力ファイル保存 ----
    print("[出力] 保存中 ...")
    if cfg.save_trials_csv:
        try:
            n = save_trials_csv(study, paths['trials_csv'])
            print(f"[出力] trials CSV → {paths['trials_csv']} ({n} 行)")
        except Exception as e:
            print(f"[出力] trials CSV 失敗: {type(e).__name__}: {e}")
    if cfg.save_best_json:
        try:
            save_best_json(study, paths['best_json'])
            print(f"[出力] best JSON → {paths['best_json']}")
        except Exception as e:
            print(f"[出力] best JSON 失敗: {type(e).__name__}: {e}")
    if cfg.save_topk_report and entries:
        try:
            save_topk_report(entries, paths['topk_report'])
            print(f"[出力] top-k 比較 → {paths['topk_report']}")
        except Exception as e:
            print(f"[出力] top-k report 失敗: {type(e).__name__}: {e}")

    # ---- L1: Feasibility 分類解析 ----
    if cfg.run_feasibility_analysis and n_completed > 0:
        try:
            from optimization.feasibility import analyze_feasibility
            print("[feasibility] 分類解析を実行中 ...")
            analyze_feasibility(
                study       = study,
                output_dir  = out_dir,
                prefix      = paths['feasibility_prefix'],
                target_type = cfg.feasibility_target,
                model       = cfg.feasibility_model,
            )
        except ImportError:
            print("[feasibility] スキップ: scikit-learn が未インストール")
        except Exception as e:
            print(f"[feasibility] 解析失敗: {type(e).__name__}: {e}")

    # ---- ベスト候補の詳細表示 (exp1 と同じレポート) ----
    be = best_entry(entries)
    if cfg.display_best_full and be is not None:
        try:
            _display_best_full(be, cfg, config)
        except Exception as e:
            print(f"[詳細表示] 失敗: {type(e).__name__}: {e}")

    # ---- 出力ガイド + README 生成 ----
    _write_readme(out_dir, cfg, paths, study, be, timestamp)
    print()
    print("=" * 72)
    print(f"成果物: {out_dir.resolve()}/")
    print("=" * 72)
    print(f"  📌 README.md       … 結果の見方ガイド (最初に開いて)")
    print(f"  ★ topk.txt        … top-{cfg.n_topk} 候補の詳細比較 (★最終結果はここ)")
    print(f"  ・ best.json       … BO ベスト trial (JSON、ベスト 1 件のみ簡易)")
    print(f"  ・ trials.csv      … 全 {cfg.n_trials} trial 履歴 (Excel/pandas で解析)")
    if cfg.run_feasibility_analysis:
        print(f"  ・ feasibility.txt   … 収束分類学習レポート (特徴量重要度)")
        print(f"  ・ feasibility_2d.png … feasible/infeasible 2D 散布図")
    print(f"  ・ optuna.db       … SQLite (中断・再開・dashboard 用)")

    # ---- 終了メッセージ ----
    print("=" * 72)
    if bo_interrupted:
        print(f"中断完了。再開コマンド: 同じ main.py を再実行 (storage URL が同じ)")
        print(f"  → storage: {storage_url}")
    elif bo_fatal_error is not None:
        print(f"異常終了。エラー: {type(bo_fatal_error).__name__}: {bo_fatal_error}")
        raise bo_fatal_error
    else:
        print("完了。")

    return {
        'study':         study,
        'topk_entries':  entries,
        'best_entry':    be,
        'paths':         paths,
        'timestamp':     timestamp,
    }
