# optimization/ — ベイズ最適化 (Optuna TPE) 層

## 役割

PDH フローシートの設計変数 (反応器・PSA・膜・蒸留塔 Dist1/2/3・原料流量・
各塔 recovery) を、Optuna の TPE (Tree-structured Parzen Estimator) ベイズ
最適化で探索し、`effective_TAC` (= TAC + ソフトペナルティ) を最小化する層。

最適化は 2 段構成:

1. **BO ループ** … 全塔を高速な FUG shortcut で評価し、多数の trial を回す。
2. **top-k 再評価** … BO 上位候補だけを rigorous (+任意で Stage 2 HEN) で
   再評価し、FUG の楽観性バイアスを是正して最終解を確定する。

`main.py` / `sub/sub1.py` / `sub/sub2.py` は本モジュールを薄く
オーケストレーションする。最適化アルゴリズムの理論・式は
`optimization/SPEC_bo.md` を参照。

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `__init__.py` | 公開シンボルの集約。`optuna` / `sklearn` の有無で条件付き import。 |
| `search_space.py` | `SEARCH_SPACE` スキーマの検証・suggest、params → `FlowsheetDesignVars` 変換 (軸流/径方向流対応)。 |
| `objective.py` | Optuna objective ファクトリ `make_objective()`。suggest → 評価 → `effective_TAC` 返却 + 診断情報を `user_attrs` に格納。 |
| `penalty_scale.py` | adaptive penalty scaling。trial 進行に応じてソフトペナルティ係数を線形に強化 (`set_scale`/`get_scale`/`default_schedule`)。 |
| `study.py` | Optuna Study 生成・最適化ループ。2-phase sampler (Sobol QMC → TPE)、TPE 用 `constraints_func` を定義。 |
| `feasibility.py` | (L1 事後解析) study から収束/feasibility の二値分類器 (Random Forest) を学習し、特徴量重要度・2D 散布図を出力。`sklearn` 任意。 |
| `callbacks.py` | BO ループ用の compact 表示 callback (1 trial = 数行の構造化ログ + 進捗/ETA + 失敗モード集計)。 |
| `topk.py` | BO 上位 k 候補の rigorous (+Stage 2) 再評価。`TopKEntry`/`select_topk`/`reevaluate_topk`/`best_entry`。 |
| `reporting.py` | 結果出力。全 trial CSV・ベスト JSON・top-k 比較 txt (BO vs 再評価)。 |
| `parallel.py` | マルチプロセス並列最適化。共有 SQLite study を N worker で分担 (`spawn_workers` / worker CLI)。 |
| `pipeline.py` | BO + top-k + L1 + 詳細表示の一括オーケストレーション (`PipelineConfig` / `run_pipeline`)。 |

## 最適化の流れ

`run_pipeline(PipelineConfig)` (pipeline.py) が以下を実行する:

1. **検証** … `validate_search_space` で `SEARCH_SPACE` のキー・型・範囲を検査。
2. **パス準備** … `outputs/main_<timestamp>/` に成果物 (optuna.db, trials.csv,
   best.json, topk.txt, feasibility.txt, feasibility_2d.png, README.md) を集約。
3. **objective 構築** … `make_objective()` で `objective(trial) -> effective_TAC` を生成。
   trial 開始時に penalty scale を更新し、診断情報を `trial.user_attrs` に格納。
4. **study 生成** … `create_study()`。sampler は TPE (前段に Sobol QMC startup)、
   `direction='minimize'`、SQLite storage で中断・再開可。
5. **warm-start (任意)** … 既知良設計を `study.enqueue_trial()` で先頭注入。
6. **BO ループ実行** … `run_optimization()` で N trial を回す。全塔 FUG。
   `n_workers>1` のときは `parallel.spawn_workers()` でマルチプロセス並列。
   compact callback でライブ表示。Ctrl+C / 致命的例外でも部分結果を保存。
7. **top-k 再評価** … `reevaluate_topk()` で上位候補を rigorous (+任意 Stage 2)
   再評価。`best_entry()` が feasible 優先で最良候補を選ぶ。
8. **出力保存** … `save_trials_csv` / `save_best_json` / `save_topk_report`。
9. **L1 解析 (任意)** … `analyze_feasibility()` で収束分類器を学習・出力。
10. **詳細表示** … ベスト候補について `display_full_results` を出力 (再評価結果を再利用)。

## 位置づけ・依存

- 上位呼び出し: `main.py` (HYSYS/SM), `sub/sub1.py` (FUG/rigorous), `sub/sub2.py` (SM/rigorous + Stage 2)。
- 下位依存: `flowsheet.evaluate` / `FlowsheetResult`, `config.load`,
  `src.distillation_core.ColumnTunables`, 各 unit の `DesignVars`
  (`units.reactors.swing`, `units.reactors.radial_flow`,
  `units.separators.psa`, `units.separators.membrane`)。
- `optuna` は必須 (study/topk/reporting)。未インストール時は
  `search_space`/`objective`/`pipeline` のみ公開され、警告が出る。
- `sklearn` は `feasibility.py` (L1) のみ任意依存。
- 並列実行には `save_sqlite=True` (共有 SQLite storage) が必須。
- 暫定値・既知の不備の扱いは リポジトリ直下 `KNOWN_PLACEHOLDERS.md` を参照。
