# simulation/ — 結果表示・実験ランナー

## 役割

実験スクリプト (`exp/exp*.py` 等) や最適化パイプラインから共通で使う
「実験ごとに変えない部分」を集約する。具体的には (1) フローシート評価結果の
整形・表示、(2) 実験スクリプトの共通ランナー (stdout キャプチャ・進捗 ticker・
ファイル保存)。設計変数の具体値は実験スクリプト側に直書きし、本モジュールは
表示と実行の定型処理に専念する。

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `display.py` | フローシート評価結果 (`FlowsheetResult`) の各セクション表示関数群。 |
| `exp_runner.py` | 実験スクリプト共通ランナー (stdout キャプチャ + 進捗 ticker + txt 保存)。 |
| `__init__.py` | display / exp_runner の公開関数を再エクスポート。 |

### `display.py` の主な関数

| 関数 | 説明 |
|---|---|
| `hdr(title)` | 区切り線付きの見出しを表示。 |
| `show_input_snapshot(design, config, eval_kwargs)` | 実行条件 (全 unit 設計変数 + 評価オプション) のスナップショット表示。軸流/径方向流を自動判別。 |
| `show_stream(label, stream)` | 単一ストリームの組成・T・P を表示。 |
| `show_hi_summary(result)` | ヒートインテグレーション (HI) 結果の要約。 |
| `show_stage2_synthesis(result)` | Stage 2 (HEN 合成) 結果の表示。 |
| `show_final_summary_box(result, ...)` | 最終サマリ (TAC・収支等) のボックス表示。 |
| `display_full_results(result, design, config)` | 上記を含む全セクションをまとめて出力する総合表示。 |

(他に `show_streams_overview`, `show_unit_details`, `show_production`,
`show_capex`, `show_opex`, `show_revenue`, `show_specs`, `show_tac_summary`
等の内部セクション関数を持つ。)

### `exp_runner.py` の主な関数

| 関数 | 説明 |
|---|---|
| `run_with_capture(run_fn, label, output_dir, ...)` | `run_fn()` を実行し、stdout を `StringIO` にキャプチャしつつ stderr に進捗を流し、`outputs/<label>_<ts>.txt` に保存。例外時もログを残す。 |
| `outer_iter_progress(expected_iters)` | 外側ループ反復数から進捗文字列を返す `progress_phrase` ファクトリ。 |
| `run_exp(label, eval_callable, ...)` | exp 系スクリプトのデフォルト実行ラッパ (`run_with_capture` の薄い前面)。 |

## 位置づけ・依存

- `display.py` は `stream.ProcessStream` と `src.cost_parameters` (単価等の表示用) に依存。
  成分キー (A–Z) の名称対応表 `_COMP_NAMES` を持つ (`stream/README.md` 参照)。
- `exp_runner.py` は標準ライブラリのみ (`contextlib`, `io`, `threading`, `pathlib`)。
- 最適化パイプライン (`optimization/pipeline.py`) はベスト候補の詳細出力に
  `show_input_snapshot` / `display_full_results` を再利用する。
