# comparing/ — 欠陥最適化手法の再現と ΔTAC 比較ケース群

## 役割

過去の学部レポートで広く使われてきた **欠陥のある最適化手法 (逐次1次元・粗グリッド・部分最適化・
後置ピンチ・検証なし 等) を本シミュレータ上で忠実に再現**し、その best 設計を本体 BO (`../main.py`,
全変数 + 熱統合を同時最適化) の best と突き合わせて **ΔTAC (損失額) を定量化**する比較ケース群。
これにより「素朴手法に対して BO がどれだけ得をするか」を本物のプラント TAC で示す。

## 一括実行 (推奨)

全 case を 1 コマンドで直列実行し、BO ベスト (`../outputs/main_*/best.json`) との
ΔTAC = TAC(本手法) − TAC(BO) を比較表 (CSV)・棒グラフ (PNG)・レポート貼付用 Markdown 表に
まとめるエントリ `run_all.py` を用意している (HYSYS のある PC で実行)。

```powershell
# 計画と import 検証だけ (HYSYS 不要)
.\.venv\Scripts\python.exe comparing\run_all.py --dry-run
# 全 22 ケース実行 (HYSYS 直列、~33h 目安。1 ケース完了ごとに results へ逐次保存)
.\.venv\Scripts\python.exe comparing\run_all.py
# 一部だけ
.\.venv\Scripts\python.exe comparing\run_all.py --only case_rep_styrene2025,case_rep_eo2025
# BO 基準を差し替えたい場合のみ
.\.venv\Scripts\python.exe comparing\run_all.py --baseline outputs\main_20260601_150117\best.json
```

成果物は `comparing/results/comparison_<ts>.{csv,md,png}`。HYSYS は単一 COM インスタンスのため
**直列実行のみ (並列不可)**。BO 基準は既定で**レポートと同一の同梱 `baseline_best.json`**
(trial #194 / TAC 1056.3 億円/年) を使う (ラボ PC でもレポートと揃う)。`--baseline` で上書き可。

## case フォルダ構成

各 case は独立スクリプト `comparing/<case>/main.py`。実行は
`.\.venv\Scripts\python.exe comparing\<case>\main.py`。共通土台は `shared/` を流用し、
結果は `comparing/results/<method>_<ts>/` に main.py と同形式 (trials.csv / best.json /
top-N 詳細 / README) で出力される。蒸留塔 backend は Dist1=SM / Dist2=HYSYS / Dist3=SM
(FUG は P10 の可視化目的でのみ使用)。**HYSYS 実走は HYSYS のある PC で** (本 PC では import 確認まで)。

命名規約 `case_<カテゴリ>_<名前>`:

| カテゴリ | 意味 |
|---|---|
| `case_p##_*` | 欠陥そのものの再現 (P## = 問題点レポートの P01〜P12)。汎用の手法部品。 |
| `case_combo_*` | 複数 C/P の意味ある組み合わせ (実在の典型ワークフローをまるごと再現)。 |
| `case_rep_*` | 実在レポートの手法をそのまま再現 (テーマ非依存で PDH sim に写す。著者名は出さず匿名 ID)。 |

## ファイル / サブ構成一覧

| パス | 内容 |
|---|---|
| `CASES.md` | 全 case の単一索引 (命名規約・欠陥対応表・実在レポート再現の成熟度一覧)。**まずここを読む。** |
| `REPORT_METHODS_ANALYSIS.md` | `case_rep_*` の手法・出典 (匿名 ID) の詳細分析。 |
| `PROGRESS_AND_PLAN.md` | 進捗と未着手項目 (P03/P07-09/P11 等) の計画。 |
| `case_p01〜p12_*/main.py` | 欠陥部品の再現 (部分最適化 / 後置ピンチ / 逐次1次元 / 粗グリッド / マルチスタート / FUG 精度問題 / 未収束)。 |
| `case_combo_typical/main.py` | 典型的な学部レポート全体 (C01 逐次 → C04 後置ピンチ、P02+P04+P05+P06+P12) の忠実再現。 |
| `case_rep_*/main.py` | 実在レポート (PDH・EO・スチレン・トルエン・ブタジエン・DME・メタノール 等 8 テーマ、成熟度 2〜7) の手法再現。 |
| `shared/harness.py` | 素朴手法を「enqueue 駆動の Optuna study」として走らせる共通エンジン (RandomSampler、各点を決定的に enqueue)。 |
| `shared/simulator.py` / `space.py` / `reporting.py` | 評価器ラッパ / 探索空間定義 / 結果レポート出力 (main.py 流用)。 |

## 使い方・位置づけ

- 比較は利用者側で行う: 各 case の `best.json` と `outputs/main_*/best.json` (BO best) を突き合わせて ΔTAC を算出する。
- リポジトリには過去レポート著者の人名を残さない (case は匿名 ID で対応付け)。
