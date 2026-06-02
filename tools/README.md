# tools/ — 実用ツール群

## 役割

最適化・解析を補助する単発実行ツール。ノートブックの生成、並列 run の進捗監視、
プロファイリング、事後 feasibility 解析など。venv の python で実行する。

## ファイル一覧

| ファイル | 用途 |
|---|---|
| `build_reactor_ceiling_nb.py` | `monitor/reactor_conversion_ceiling.ipynb` を生成 (単通転化率~28%・選択率~80% が熱力学・反応速度論で必然であることを示すノート)。 |
| `build_stage_comparison_nb.py` | `monitor/stage_comparison.ipynb` を生成 (単段断熱床→径方向流→多段化の変遷と段数 1/2/3/4 段の比較、なぜ 3 段かを示すノート、HYSYS 不要)。 |
| `monitor_main.py` | sub1/sub2 並列 run の全体進捗 (X/N trial・feasible・best・ETA) を共有 SQLite から read-only 表示。`--watch` で定期更新。 |
| `profile_flowsheet.py` | `flowsheet.evaluate()` を cProfile で関数別累積時間を出力し、収束高速化の真のボトルネックを特定する。 |
| `run_feasibility.py` | 既存 Optuna study (`outputs/main_*.db`) を後から読み込み、feasibility 分類器を再学習する事後解析 (`--target` / `--model` 切替可)。 |
| `top_designs.py` | sub1/sub2 並列 run の feasible 設計を 収率 / Dist3 段数 / 利益 軸で一覧 (`--sort tac\|profit\|yield\|n3`)。共有 SQLite を read-only 参照。 |

## 使い方・位置づけ

- 例: `python tools\monitor_main.py --watch 15` / `python tools\top_designs.py --sort profit` / `python tools\profile_flowsheet.py --top 30`。
- `monitor_main.py` / `top_designs.py` は SQLite を read-only で読むだけなので、実行中の最適化に干渉しない (WAL で並行読取り可)。
- nb ビルダーは生成後に nbconvert --execute で出力を埋める運用。
