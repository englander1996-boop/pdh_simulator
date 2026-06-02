# sub/ — 別バックエンド構成の最適化アーカイブ

## 役割

現行 BO の本丸は本体 `../main.py` (蒸留塔 backend = Dist1=SM / Dist2=HYSYS / Dist3=SM)。
`sub/` はそれとは **別のバックエンド構成を持つ最適化ドライバのアーカイブ**で、HYSYS の単一 COM 制約を
避けてマルチプロセス並列を得たり、過去構成を保存・再利用するために残してある。
いずれも Optuna ベースの全変数 制約付き最適化で、`outputs/` に結果一式を出力する。

## ファイル一覧

| ファイル | 構成 | main.py との違い |
|---|---|---|
| `sub1.py` | **FUG / rigorous 全系** (旧 main 系)。BO ループは全塔 FUG で高速 (1 eval ~3 秒)、top-k 候補だけ rigorous + Stage2 で精密再評価。parallel kind=`sub1`。 | Dist1/Dist3 が BO ループでは FUG (低精度) であり、Stage2 (HEN 合成) も top-k のみ → BO の目的関数と再評価に乖離がある。 |
| `sub2.py` | **SM / rigorous / SM 全系** (旧 final 系)。Dist1/Dist3=学習済み SM(GPR)、Dist2=in-house rigorous (Wang-Henke)。全塔 pure Python ゆえ 6 worker 並列可。Stage2 を全 trial で実行。parallel kind=`sub2`。 | Dist2 を HYSYS → rigorous に置換した main.py のフォーク。HYSYS COM を使わないので並列化でき、全 trial で Stage2 を回すため目的関数が本物の TAC にほぼ一致する。 |

構成の精度・速度トレードオフ:

- `main` (sm, hysys, sm): 精度◎だが HYSYS COM で並列不可・遅い。
- `sub1` (fug, rig, fug → top-k rig): 並列可だが BO ループの Dist1/Dist3 が FUG (低精度)、Stage2 は top-k のみ。
- `sub2` (sm, rig, sm) + Stage2 全 trial: SM が Dist1/Dist3 を HYSYS 精度・FUG 速度で置換、Dist2 は精度保護で rigorous 据え置き → main の精度構造 × sub1 の並列性。

## 使い方・位置づけ

- `sub1.py`: ファイル冒頭の §1〜5 ブロック (試行数・ソルバ・探索範囲・出力) を編集し `.\.venv\Scripts\python.exe sub\sub1.py`。
- `sub2.py`: `.\.venv\Scripts\python.exe sub\sub2.py > outputs\sub2_run.log 2>&1`。`exp/` の感度スクリプトの一部は `sub2` を import して評価器に使う。
- 並列 run の進捗・上位設計の俯瞰は別ターミナルで `tools/monitor_main.py` / `tools/top_designs.py` を使う (共有 SQLite を read-only 参照)。
