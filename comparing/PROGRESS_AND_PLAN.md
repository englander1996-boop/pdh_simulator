# comparing/ 引継ぎ — 進捗と今後の詳細方針

最終更新: 2026-05-29 (命名統一・実在レポート14本・P10/FUG・各case内に欠陥部品節 まで) / 作成者: Claude

このドキュメントは **再開のための引き継ぎ書 (なぜ/経緯/設計判断/再開手順)**。
**全 case の最新一覧は `CASES.md`**、**実在レポートの手法分析は `REPORT_METHODS_ANALYSIS.md`** が正
(このファイルと食い違ったら CASES.md / REPORT_METHODS_ANALYSIS.md を優先)。

---

## 0. いま再開する人がまず知るべきこと (TL;DR)

- `comparing/` は **過去レポートの「問題のある最適化のやり方」を本シミュレータ上で忠実に再現し、
  その欠陥の損失 (ΔTAC = TAC(本手法) − TAC(BO)) を定量化する** パッケージ。BO (= `special.py`) の有用性を逆説的に立証するのが目的。
- **構成**: `case_<カテゴリ>_<名前>/main.py` が **独立スクリプト**として単独で走る。中央ディスパッチャ無し。共通土台のみ `shared/`。
  - `case_p##_*` = **欠陥部品** (P01〜P12 の単体実演): p01_subsystem/p02_pinch/p04_sequential/p05_grid/p06_multistart/p10_fug/p12_converge (7本)
  - `case_combo_*` = **組合せ**: combo_typical (C01逐次→C04後置の典型ワークフロー)
  - `case_rep_*` = **実在レポート再現** (14本・11テーマ・成熟度2〜7。各 case 内に「含まれる欠陥部品」節あり)
  - = 計 **22 case**。索引は CASES.md。
- **実装・import 確認は全 case 完了**。残りは **実走行 (HYSYS) と BO との突合**のみ。
- backend は special.py と同じ **Dist1=SM / Dist2=HYSYS / Dist3=SM**。**FUG は既定では使わない**が、
  例外として **P10 (case_p10_fug) は FUG vs HYSYS真値 の精度乖離を見せる目的でのみ FUG を使う**。
- **HYSYS はこの PC では動かせない** (ユーザ環境)。→ このマシンでは import/構文確認まで。実走行は HYSYS のある PC へ移してから。
- 既知バグは全て解消済み (循環 import / cp932 UnicodeError)。
- 次の一手: **HYSYS のある PC で各 `comparing\<case>\main.py` を実行 → 各 results を BO (outputs/special_*) と突合し ΔTAC**。
  内容調整 (掃引変数・点数) は各 case 先頭の設定定数 (`VARS_ORDER`/`GRID_VARS`/`BLOCKS` 等) を編集するだけ。
  テーマ追加は `REPORT_METHODS_ANALYSIS.md` の手順に従い `case_rep_*` を増やす (手法はテーマ非依存)。

---

## 1. 目的・背景 (なぜ作っているか)

ユーザの研究: PDH (プロパン脱水素) プロセス最適化。自作の `pdh_simulator` は
**ベイズ最適化 (BO) + rigorous モデル + 並列** が売り。

過去 7 年分 (2019–2025) の学生レポートを分析した 2 つのレポート (デスクトップ、§7 参照) が、
- 誰も高度な最適化 (BO/NLP/MINLP/多目的) を使っておらず C01 逐次1次元 と C04 ピンチ熱統合に集中、
- 方法論的欠陥 P01〜P12 (リサイクル不変数化・整数の連続扱い・大域性未保証・単目的・検証なし 等) が蔓延、
を指摘している。

**ユーザの狙い (確定):** 「こういう欠陥があるから BO を採用した」を、逆に過去の素朴な手法を
**同じシミュレータ・同じ目的関数の土俵で実際に再現・実行**し、BO best との差で
「この欠陥は年間 X 億円損している (+Y%)」と **定量的に**示す。

### フレーミング上の最重要点 (ユーザが繰り返し強調)
- 作るのは「素朴な最適化器」ではなく、**「問題のあるやり方そのもの」の忠実再現**。
- 「それが問題である」ことを **定量的に表現**する (ΔTAC = TAC(本手法) − TAC(BO))。
- PDH に限定しなくてよい (論理破綻が明白な手法が多く、軽い再現で十分露呈する)。
- 選定は「レポート単位」でなく **「問題 (P) 単位」** で、再現したら測れる悪い結果が出るものを狙う。

---

## 2. 確定済みの設計判断 (勝手に変えないこと)

| # | 判断 | 理由 |
|---|------|------|
| D1 | 各手法は **自分の `comparing/results/<method>_<ts>/` に結果を吐くだけ**。横断比較レポートは作らない | 横断比較はユーザ側でやる (2026-05-28 合意)。`report.py` や BO 読込は不要 |
| D2 | **`special.py` の I/O・ライブログ・結果保存・ファイル構成を忠実再現** | special のファイル構成はエラー調査・分析がしやすい。`reporting.py` がその移植 |
| D3 | 各手法 = **enqueue 駆動の Optuna study (RandomSampler、決定的)** | special の作り込んだ資産 (callback / trials.csv / top-N 詳細 / README / `_store_diagnostics`) は study/trial 前提。enqueue で探索点を与えれば全部流用できる |
| D4 | backend は **Dist1=SM / Dist2=HYSYS / Dist3=SM 固定** (special.py:225-232 と一致 = 公平比較) | BO も素朴手法も同一 backend・同一目的でないと「手法だけの差」が切り出せない |
| D5 | `comparing/` は **`special.py` を import しない**。`space.py` は SEARCH_SPACE / build_design の **手動コピー (単一の真実)** | `special.py` → `main.py` 改名・`final.py`/旧`main.py` 削除予定。改名に巻き込まれない疎結合に |
| D6 | purity は **99.45 wt% に緩和** (`simulator.CONFIG`) | special.py の決定A (2026-05-25)。SM Dist3 の 99.5mol%=99.497wt% を尊重 |
| D7 | BO ベースラインは **`special.py` の既存 run (`outputs/special_*/best.json`) を流用**。再走不要 | 同一 backend で既に走っているので公平。比較時にユーザが突合 |
| D8 | **FUG は既定では使わない**。例外: **P10 (case_p10_fug) のみ FUG を使い、その精度問題を可視化する** | ユーザ指示 (2026-05-29)。FUG を黙って既定にするのが NG (proxy_penalty が乗り無意味)。ただし「FUG を使うなら精度問題を必ず見せる」= P10 の本質。case_p10_fug は同設計を FUG と HYSYS真値で二重評価し乖離を出す |
| D9 | **case ごとにフォルダ + `<case>/main.py` を独立スクリプト化。中央ディスパッチャ無し** | ユーザ指示 (2026-05-29)。まとめて比較する作りはエラーの元 (1 手法のバグが全体を巻き込む)。各 case は単独実行 |
| D10 | **HYSYS はこの PC では動かない**。実走行は HYSYS のある PC で | ユーザ環境 (2026-05-29)。このマシンでは import/構文確認まで。FUG 代替はしない |
| D11 | **真値は HYSYS** (`rigorous` solver ですら HYSYS より精度が劣る) | ユーザ指示 (2026-05-29)。FUG 精度比較 (case_p10_fug) の「真値」側は Dist2=HYSYS を含む special backend を使う |
| D12 | **命名は `case_<カテゴリ>_<名前>` に統一** (case_p##/case_combo_/case_rep_)。全 case 索引は CASES.md | ユーザ指示 (2026-05-29)。命名散乱の解消 |
| D13 | **過去レポート著者の人名を repo (git) に残さない**。匿名 ID (case_rep_<テーマ><年><連番>) を使う。Claude メモリ (git外) は可 | ユーザ指示 (2026-05-29)。第三者氏名を公開 repo に残さない |
| D14 | **各 case_rep の docstring に「含まれる欠陥部品」節を必ず書く** (◎=sim再現→case_p##、△=検出だが sim非対応)。分析は分析 doc でなく各 case 内 | ユーザ指示 (2026-05-29)。レポート再現は欠陥部品の束なのでそれを各 case 内で明示 |
| D15 | **テーマを PDH に偏らせない**。問題点レポートの**成熟度スコア (12−検出問題数)** を軸に低〜高成熟度を横断 | ユーザ指示 (2026-05-29)。手法はテーマ非依存なのでいろいろなテーマを写す。低成熟度ほど欠陥多→ΔTAC大 |

### リポジトリ側の予定 (背景として重要)
- `Z:\pdh_simulator\final.py` と `Z:\pdh_simulator\main.py` は **削除予定**。
- `Z:\pdh_simulator\special.py` を **`main.py` に改名して最適化の本丸**にする。
- → だから `comparing/` は `special.py` に依存させていない (D5)。

---

## 2.5 会話の流れ (時系列) と却下された方向 — 文脈ゼロの後任向け

設計がどう固まったかの経緯。**特に「却下された方向」を再提案しないこと。**

1. ユーザがデスクトップの分析フォルダ (§11) を見せ「何をやろうとしているか当てて」。
2. Claude 初回推測「過去手法を再現したい」→ **ユーザ「違う」**。
3. Claude 再推測「高度な最適化を誰も使っていない空白を立証し、自作 BO の新規性 (related work) を主張」
   → **ユーザ「半分正解」**。正確には「**こういう問題があるから BO を採用した。逆に過去の最適化手法を
   ここで再現することで BO の有用性を示す**」。
4. 2 つの詳細レポート (手法分析・問題点) を読了し構図を確認。
5. 「レポート抜粋_md があれば再現可能か?」→ Claude「再現可能。ただし **📖最適化セクション原文だけを信じ、
   🎯変数表・🔁再現テンプレ等のメタデータ残骸は無視**」と評価 (§11)。
6. ユーザ補足: 「軽い分析でよい。全部やると大変。**調査をもとにどれを再現するか決める。内包される問題が
   露わになる結果が得られたらよい**」「**PDH に限らなくてよい**」「論理破綻が明白なので再現は容易」
   「**ファイルが大量なので構成も課題**」。
7. Claude「**レポート単位でなく問題 (P) 単位で選ぶ**」方針 + 露呈度の表を提示。置き場所を質問
   → **ユーザ「pdh_simulator リポジトリ内」**。
8. ユーザ「**`comparing/` 内に 10 個くらいファイルを作り、その中で pdh_simulator の構成を再現。分野ごとに
   モジュール化、main.py で最適化。流用できるものは流用**」。さらに「**`final.py`/`main.py` は削除予定、
   `special.py` を `main.py` に改名して本丸にする**」。
9. Claude が本丸 `special.py` を読み、継ぎ目 (`SEARCH_SPACE` → `_build_design` → `evaluate`) を確認。
   backend を質問 → **ユーザ「本丸と同じ Dist2=HYSYS」** (D4)。
10. ディレクトリと初期ファイルを作成。**ここで Claude が簡易自作の `io.py` を書いた → 却下** (下記)。
11. ユーザ「**入出力・ログ・結果保存は special でかなり工夫している。忠実に再現せよ。ファイル構成も見やすく、
    エラー調査・分析がしやすい**」(D2) → `io.py` 破棄、`reporting.py` に special の I/O を移植。
12. ユーザ「**時間をかけていいから丁寧に**」。
13. ユーザ「**問題のあるやり方を表現するんだよ**」「**問題であることを定量的に表現するんだよ**」
    → フレーム再構築 (下記の却下②)。
14. ユーザ「**各最適化の結果が各 outputs に出力されればいい。比較はこっちでやる**」(D1)。
15. ユーザ「**ひと段落したら `comparing/` 内に md で引継ぎを保存。再開に必要な全情報を入れて**」(= このファイル)。
16. Claude が `simulator.py` を special 流に作り直し、`reporting.py` / `harness.py` / `sequential_1d.py` /
    `main.py` を作成。**smoke テスト中に循環 import で落ち、さらに API 通信エラーでセッションが中断**。
17. 再開セッション: 循環 import を修正、Unicode バグを特定 → 修正、shared+sequential_1d を FUG smoke で
    end-to-end 完走確認 (proxy_penalty は FUG 由来の想定挙動)、このファイルを拡充。
18. ユーザ追加指示 (2026-05-29): 「**FUG は使わない**」「**まとめて比較する作りはエラーの元、やめて**」
    「**`comparing/` 下に case ごとにフォルダを作り、各 `caseN/main.py` を独立スクリプトに**」
    「蒸留塔は **SM と HYSYS**」「**HYSYS はこの PC では動かせない**」。
19. Claude が再編: 中央 `main.py` 廃止、FUG/smoke 全廃、欠陥部品 case を caseN=PN に整理、shared/ 流用維持。
20. ユーザ「**実在レポートのデータを上げてるんだからそれをモデルに再現すればいい。手法はテーマ非依存。
    PDH 以外のテーマでも再現できる (読むのが面倒なだけ)。いろんなテーマを分析するのは必須**」
    → `case_rep_*` を新設し、レポート抜粋の 📖原文から手法を抽出して PDH sim 変数に写す方針に。
21. ユーザ「**FUG を使うなら精度問題をちゃんと可視化しないと**」「**rigorous ですら HYSYS より精度が悪い**」
    → P10 を「不可」から復活させ `case_p10_fug` (FUG vs HYSYS真値 の乖離可視化) を新設 (D8/D11)。
22. ユーザ「**case 名の散乱を何とかして**」「**人名を repo に残すな (メモリは可)**」「**PDH に偏りすぎ。成熟度を参考に**」
    → 命名を `case_<カテゴリ>_<名前>` に統一 (D12)、人名を全 repo から除去し匿名 ID 化 (D13)、
    成熟度 2〜7・11 テーマを横断する `case_rep_*` を整備 (D15)。CASES.md (索引) と REPORT_METHODS_ANALYSIS.md (分析) を新設。
23. ユーザ「**レポート再現も欠陥部品を含んでるんだから、それを各 case 内で説明・分析して**」
    → 全 14 本の `case_rep_*` の docstring に「含まれる欠陥部品」節 (◎sim再現→case_p## / △検出のみ非対応) を追加 (D14)。

### ❌ 却下された方向 (再提案禁止)
1. **「素朴な最適化器を作って BO と性能比較する」フレーム** → ❌。正しくは「**欠陥手順そのものを忠実に再現し、
   その損失を定量化**」(P 単位、ΔTAC)。「いい最適化器を作る」のではなく「悪いやり方を再現する」。
2. **自前の簡易 `io.py` 出力形式** → ❌。`special.py` の I/O 機構を `reporting.py` に忠実移植 (D2)。
3. **横断比較レポート (`report.py`) や BO best の自動読込** → ❌。各手法は自分の outputs を出すだけ、比較はユーザ (D1)。
4. **30×30 のような理想化した密グリッド** → ❌。過去手法は粗い (例 3×3) ので **忠実に粗く** (HYSYS が重い事情とも一致)。
5. **`special.py` を import して流用** → ❌。改名予定なので `space.py` 等に**手動コピー** (D5、§9)。
6. **FUG バックエンド / 高速 smoke** → ❌ (D8)。蒸留塔は SM と HYSYS のみ。FUG は proxy_penalty が乗り無意味。
7. **中央ディスパッチャ (`main.py` で全手法をまとめて実行/比較)** → ❌ (D9)。case ごとに独立した `main.py`。

---

## 3. ファイル構成 (現状)

```
comparing/
  __init__.py              # パッケージ目的・構成の説明
  PROGRESS_AND_PLAN.md     # このファイル
  shared/                  # ★全 case が流用する共通土台 (ここは比較ではなく「部品」)
    __init__.py
    space.py        # 21設計変数の SEARCH_SPACE + build_design (special._build_design の手動コピー)
    simulator.py    # raw_evaluate(design,...) -> FlowsheetResult。CONFIG/EVAL_KWARGS_DEFAULT 保持
    reporting.py    # special の I/O 機構の忠実移植 (callback/trials.csv/top-N/README)
    harness.py      # enqueue 駆動 study エンジン (suggest_all/make_objective/run_batch/finalize 等)
  case_p01_subsystem/main.py    # P01 部分最適化                    [欠陥部品]
  case_p02_pinch/main.py        # P02 熱統合の後置                  [欠陥部品]
  case_p04_sequential/main.py   # P04 1次元逐次 (+P05)              [欠陥部品]
  case_p05_grid/main.py         # P05 粗グリッド (+P06)             [欠陥部品]
  case_p06_multistart/main.py   # P06 大域性未保証 (multi-start)    [欠陥部品]
  case_p10_fug/main.py          # P10 FUG精度問題 (FUG vs HYSYS真値) [欠陥部品]
  case_p12_converge/main.py     # P12 収束未検証 (複数巡)           [欠陥部品]
  case_combo_typical/main.py    # C01逐次→C04後置 典型ワークフロー  [組合せ]
  case_rep_*/main.py            # 実在レポート再現 (索引=CASES.md / 詳細=REPORT_METHODS_ANALYSIS.md)
  results/          # 各 run の出力先 (実行時生成)
```
> 注: 命名は `case_<カテゴリ>_<名前>` に統一 (case_p##/case_combo_/case_rep_)。**全 case の最新一覧は CASES.md**
> (case_rep_* 14本の個別行も含む)。実在レポートの手法・成熟度・出典は REPORT_METHODS_ANALYSIS.md。

- 各 `<case>/main.py` は **独立スクリプト** (`def run(...)` 本体 + `def main(): run()` + `__main__`)。
  実行: `.\.venv\Scripts\python.exe comparing\<case>\main.py`。中央ディスパッチャは無い (D9)。
- 各 case は冒頭に `sys.path.insert(0, repo_root)` (= `__file__` の 3 つ上) と
  `sys.stdout.reconfigure(encoding='utf-8')` を持ち、単独で起動できる。掃引変数・点数は各 case 先頭の
  設定定数 (`VARS_ORDER`/`GRID_VARS`/`BLOCKS`/`PHASE1_VARS`/`N_STARTS` 等) を編集して調整する。
- case フォルダは package ではない (`__init__.py` 無し)。`comparing` と `comparing/shared` のみ package。

### 各ファイルの役割 (詳細)
- **`shared/space.py`**: `SEARCH_SPACE` (21変数: 反応器4/PSA3/膜2/原料1/Dist1×4/Dist2×4/Dist3×3)、
  `build_design(params, backend)` → `FlowsheetDesignVars`、`midpoint_params()`、`grid_points(name,k)`、
  `clamp()`、`is_int()`、`bounds()`。`DEFAULT_BACKEND={dist1:sm,dist2:hysys,dist3:sm}` (= special、FUG 無し)。
  **import 順序が重要** (§5 の循環 import 参照)。
- **`shared/simulator.py`**: `raw_evaluate(design, *, apply_hi, hi_dT_min_K, apply_stage2, F_fresh)`
  → `flowsheet.evaluate` を呼んで `FlowsheetResult` を返すだけの薄いラッパ。`CONFIG` (purity 緩和済)、
  `EVAL_KWARGS_DEFAULT` (apply_hi=True/hi_dT_min_K=10/apply_stage2=True)、`shutdown()`。
- **`shared/reporting.py`**: `make_callback(n_total)` (ライブログ)、`summarize(study)`、
  `save_trials_csv`、`save_best_reports` (top-N 再評価 → CAPEX/OPEX/spec 内訳 txt、feasible 無しはスキップ)、
  `write_readme`。**ファイル書き込みは全て `encoding='utf-8'`**。
- **`shared/harness.py`**: `suggest_all` (21変数+制御変数を suggest)、`make_objective` (eval_opts フックあり=
  後置ピンチ系 case_p02_pinch/case_combo_typical で apply_hi/ΔTmin を差し替えるのに使用)、`new_study` (RandomSampler)、
  `run_batch` (points を enqueue→評価→batch の trial を返す)、`best_of`、`finalize` (trials.csv/best.json/top-N/README を保存)、
  `new_run_dir`、`save_table_csv`。penalty scale は固定 1.0。
- **各 case の `run()` 共通形**: 範囲中央 (or 乱数) 始点 → 探索点リストを作る → `harness.run_batch` で評価 →
  `harness.finalize` で保存 → 手法固有 CSV を `harness.save_table_csv` で出力。手法固有 CSV:
  case_p01_subsystem=block_curves / case_p02_pinch=dtmin_sweep+phase1_no_hi / case_p04_sequential=cost_curves /
  case_p05_grid=grid / case_p06_multistart=multistart / case_p10_fug=fug_vs_rigorous /
  case_p12_converge=pass_curves+pass_best / case_combo_typical=phase1_cost_curves / case_rep_*=cost_curves または grid。
- **case_p10_fug の特殊性**: harness の make_objective ではなく **自前の二重評価 objective** (FUG backend と
  HYSYS真値 backend で同じ設計を 2 回評価し、FUG が feasible/安いと誤判定→真値が覆す乖離を fug_vs_rigorous.csv に出す)。

### 出力物 (各 run dir に入るもの = special と同形式)
`trials.csv` / `best.json` / `top{1..N}_trial*.txt` (CAPEX/OPEX/spec/HI 内訳) / `README.md`
（+ 各 case 固有 CSV、上記）。

---

## 4. 完了済み

- [x] `shared/` 4 モジュール (space/simulator/reporting/harness) — special からの移植・流用
- [x] **欠陥部品 7 本** (case_p01_subsystem/p02_pinch/p04_sequential/p05_grid/p06_multistart/p10_fug/p12_converge)
- [x] **組合せ 1 本** (case_combo_typical = C01逐次→C04後置)
- [x] **実在レポート再現 14 本** (case_rep_*。11 テーマ・成熟度 2〜7。各 case 内に「含まれる欠陥部品」節)
- [x] 循環 import 修正・cp932 UnicodeError 修正 (両方解消済)
- [x] 全 22 case の import 確認 (importlib で run+main 存在を確認)・人名ゼロ確認
- [x] 命名統一 (case_<カテゴリ>_<名前>)・中央ディスパッチャ廃止・CASES.md / REPORT_METHODS_ANALYSIS.md 整備
- [ ] **実走行 (HYSYS のある PC で各 case)** ← 残り。このマシンでは HYSYS 不可 (D10)
- [ ] 各 case の results を BO (outputs/special_*) と突合し ΔTAC を算出 (ユーザ側)
- [ ] (任意) さらにテーマ追加 / P03・P07・P08・P09 (要モデル拡張・仮定) / P11 は不可

---

## 5. 既知バグ (すべて解消済み — 記録として残す)

### [解消済] 循環 import
`space.py` が `src.distillation_core` を `flowsheet` より **先に**直 import していたため、
`ImportError: cannot import name 'ColumnTunables' ... partially initialized` が出ていた。
→ **`special.py` と同じ順序 (flowsheet を先) に直して解消**。
**教訓: 新規モジュールでも `from flowsheet import ...` を `from src.distillation_core import ...` より先に書く。**

### [解消済] UnicodeEncodeError (cp932)
ライブログの `✗ ✓ ★ █ ░ Δ →` 等が Windows の cp932 標準出力で encode 不可で落ちていた。
→ **各 case 冒頭 (import 直後) に `try: sys.stdout.reconfigure(encoding='utf-8') except Exception: pass` を入れて解消**
(本家 `special.py:50` と同じ)。ファイル出力は全て `encoding='utf-8'` 指定済みなので成果物には元から影響なし。
**教訓: 新規 case を作る時もこの 1 行を冒頭に必ず入れる (テンプレに含まれている)。**

> 現時点でオープンなコードバグは無い。残りは「HYSYS のある PC での実走行」だけ (D10)。

---

## 6. 再開手順 (順番に)

前提: 実装・import 確認は全 22 case 完了。残りは **HYSYS のある PC での実走行と BO 突合**。

1. **HYSYS のある PC にリポジトリごと移す** (このマシンでは HYSYS 不可、D10)。`.venv` が壊れていれば `fix_venv.ps1` (§10)。
2. **疎通: 軽い case を 1 本走らせる**。例:
   `.\.venv\Scripts\python.exe comparing\case_p04_sequential\main.py`
   - `comparing/results/sequential_1d_*/` に trials.csv / best.json / README.md / cost_curves.csv が出るか。
   - ライブログ (★✓✗ 等) が文字化けせず出るか (utf-8 reconfigure 済)。feasible 0 でも finalize がクラッシュしないこと
     (`save_best_reports` は feasible 無しでスキップする設計)。~150s/eval × 評価数。
3. **各 case を順に走らせる** (`comparing\<case>\main.py`、一覧は CASES.md)。重い grid 系は先頭 `K_POINTS` を落として時短可。
4. **BO 突合 (ユーザ側)**: 各 `best.json` と `outputs/special_*/best.json` を突合し **ΔTAC = TAC(case) − TAC(BO)** を算出。
   case_rep_* は成熟度が低いほど ΔTAC が大きくなる想定 (各 case の「含まれる欠陥部品」節参照)。
5. **case_p10_fug** は `fug_vs_rigorous.csv` で「FUG が feasible/安いと誤判定→HYSYS真値が覆す」乖離を確認 (P10 の証拠)。

> メモリ方針: 「import が通る」≠「実際に走る」。挙動主張は必ず実コードを走らせて裏取り
> (このマシンは HYSYS 不可なので実走行検証は HYSYS PC へ持ち越し)。

---

## 7. 残り作業・今後の拡張

欠陥部品 (P01/02/04/05/06/10/12)・組合せ・実在レポート 14 本は実装済み (CASES.md)。今後やれること:

### 7.1 さらにテーマを足す (`case_rep_*` を増やす)
手順は `REPORT_METHODS_ANALYSIS.md` の「今後」節に従う:
1. レポート md の **📖原文だけ**読み、最適化因子・順序/グリッド・刻み・報告最適値を抜く。
2. PDH sim の最寄り設計変数に写す (写せない変数=パージ率/O2比/吸収塔等は明記して置換 or 省略)。
3. `case_rep_<テーマ><年>/main.py` を **既存 case_rep_* をコピー**して作る (VARS_ORDER か GRID_VARS を差し替え)。
4. docstring に「**含まれる欠陥部品**」節 (◎=sim再現→case_p## / △=検出だが sim非対応) を必ず入れる (D14)。
5. **成熟度を散らす** (D15)、**著者名を書かない** (D13、匿名 ID。対応は Claude メモリに)。

### 7.2 未実装の P (要追加対応)
- **P03** リサイクル/パージ固定 … パージ率が sim 未露出 → flowsheet 本体にパージ変数を足す改造が要る。
- **P07/P08/P09** 単目的/NPV/CO2 … CO2 係数・割引率・寿命等の**仮定を置けば作成可** (新 `case_p07/p08/p09`)。
- **P11** 定常のみ … 動的機能なし → この sim では不可。
- **C05-C10 高度手法** (SQP/Nelder-Mead/GA/MINLP/NSGA-II) … GA/MINLP/NSGA-II は数百〜千 eval で HYSYS では非現実的。
  SQP/Nelder-Mead (scipy) のみ現実的。BO との比較軸として欲しければ追加検討。

### 7.3 新しい case を作る時のテンプレ (既存 case をコピーが速い)
冒頭: `_REPO_ROOT = __file__ の 3 つ上` + `sys.path.insert` + `sys.stdout.reconfigure(encoding='utf-8')` (cp932 回避必須)、
`from comparing.shared import space, simulator, reporting, harness`、設定定数、`run()` (points 生成→
`harness.run_batch`→`harness.finalize` + `save_table_csv`)、`def main(): run()`、`if __name__=='__main__': main()`。
case フォルダは package ではない (`__init__.py` 不要)。

---

## 8. 流用している既存リポジトリ資産 (依存関係)

`comparing/` は物理・経済モデルを一切書き直さず、以下を流用:
- `flowsheet`: `evaluate`, `FlowsheetDesignVars`, `FlowsheetResult`
- `config.load.load_operating_config`
- `src.distillation_core.ColumnTunables`
- `units.reactors.swing.DesignVars`, `units.separators.psa.psa_system.PSADesignVars`,
  `units.separators.membrane.membrane_system.MemDesignVars`
- `optimization.study`: `make_sampler`, `run_optimization`
- `optimization.objective._store_diagnostics` (trial.user_attrs に診断値を詰める)
- `optimization.penalty_scale.set_scale`
- `optimization.callbacks`: `_fmt_dur`, `_fmt_reason_from_trial`, `_fmt_tally`
- `simulation`: `display_full_results`, `show_input_snapshot`

**import 順序の注意**: `flowsheet` を `src.distillation_core` より先に import すること (循環回避、§5)。

---

## 9. 同期義務 (single source of truth はリポジトリ本体、コピーは手で追従)

- **`shared/space.py` の `SEARCH_SPACE` / `build_design`** は `special.py` の `SEARCH_SPACE` /
  `_build_design` の **手動コピー**。`special.py` (→将来 `main.py`) を変えたら **手で同期する**。
- **`shared/reporting.py`** は `special.py` の `_make_special_callback` / `_save_trials_csv` /
  `_save_best_reports` / `_write_readme` の移植。special 側を変えたら同期。
- **`simulator.CONFIG` の purity 緩和 (99.45wt%)** は special.py の決定Aに追従 (D6)。

---

## 10. 環境メモ

- 実行は `.\.venv\Scripts\python.exe`。venv が「pythoncore not found」で落ちる場合は
  `fix_venv.ps1` (pyvenv.cfg を現存 Python 3.13 に向け直す) を使う。PC をまたぐと再発しやすい。
- HYSYS は単一 COM インスタンス・**並列不可**。`raw_evaluate` は必ず単一プロセス・直列で呼ぶ
  (special.py が N_WORKERS=1 の理由と同じ)。

---

## 11. 参照: 元レポート (デスクトップ)

- 手法分析: `C:\Users\yuish\OneDrive\デスクトップ\Analysis_of_previous_research\最適化手法分析\分析結果\最適化手法分析_詳細レポート.md`
  (C01〜C10 分類。C01 逐次1次元 41.8% + C04 ピンチ 54.9% で 96.7%、C05〜C10 はゼロ)
- 問題点: `...\最適化手法分析\分析結果\最適化方法論の問題点_詳細レポート.md`
  (P01〜P12 の方法論的欠陥。§4.2 に IntegratedSimulator/Optimizer/Validator の統一 API 案、§4.1 にロードマップ)
- レポート抜粋 (再現の核): `...\最適化手法分析\分析結果\レポート抜粋_md\<カテゴリ>\<年>\*.md`
  の **「📖 最適化セクション原文」** が信頼できる (変数・順序・目的関数・範囲・報告最適値)。
  「🎯 最適化変数の表」「🔁 再現テンプレート」は正規表現の残骸/ボイラープレートなので無視。

---

## 12. 用語集 (前提知識ゼロの後任向け)

### プロセス
- **PDH**: プロパン (C3H8) 脱水素 → プロピレン (C3H6)。本プロジェクトの対象プロセス。
- フローシート構成: 反応器 (Swing 触媒) → PSA → 膜分離 → 蒸留 3 塔 (Dist1/Dist2/Dist3)。
- **Dist1/2/3**: 蒸留塔。`solver_method` は `fug`(短絡) / `rigorous`(厳密) / `sm`(自作 SM) / `hysys`(HYSYS COM) の 4 択。
  既定は Dist1=SM/Dist2=HYSYS/Dist3=SM。**`fug` は P10 (case_p10_fug) の精度問題可視化でのみ使う** (D8)。
  精度序列: FUG < SM/rigorous < HYSYS (真値は HYSYS、D11)。
- **HYSYS**: 商用プロセスシミュレータ。COM 経由で呼ぶ。単一インスタンス・**並列不可**・~150s/eval。この PC では不可 (D10)。

### 評価値 (FlowsheetResult / trial.user_attrs / best.json に入る)
- **TAC** (Total Annualized Cost): 年間総コスト [億円/年]。CAPEX 年償却 + OPEX。
- **effective_TAC**: 最適化が最小化する目的値 = TAC + 制約違反ペナルティ [億円/年]。`objective` の返り値。
- **feasible / is_feasible**: 純度・生産量帯などの制約を満たすか。満たさないと penalty が乗る。
- **purity (c3h6_purity_wtfrac)**: 製品プロピレン純度 (重量分率)。spec 下限 99.45wt% (D6)。
- **production_kmol_h / yield**: 生産量と収率。`yield = production / F_fresh`。
- **economics / economics_hi / economics_synth**: HI 前 / HI 適用後 / HEN 合成後の経済内訳。
- **apply_hi / hi_dT_min_K / apply_stage2**: 熱統合を適用するか / ピンチ ΔTmin [K] / Stage2 (HEN 合成) を行うか。
  既定 True / 10.0 / True (`simulator.EVAL_KWARGS_DEFAULT`)。

### 最適化
- **BO** (Bayesian Optimization): 本丸 `special.py` が使う賢い手法。Optuna TPE。比較の「正解側」。
- **Optuna study / trial**: 最適化の実行単位 / 1 評価点。`enqueue_trial` で探索点を外から固定指定できる。
  本パッケージは学習しない RandomSampler study に点を enqueue するだけ (= 決定的)。
- **座標降下 / Gauss-Seidel**: 1 変数ずつ他を固定して最適化する素朴手法 (= P04、sequential_1d が再現)。

### 分類コード (元レポート由来)
- **C01〜C10**: 最適化「手法」の分類。実態は C01 逐次1次元 (41.8%) と C04 ピンチ熱統合 (54.9%) に集中、
  C05〜C10 (SQP/直接探索/RSM/メタヒューリスティック/数理計画/多目的) はゼロ。
- **P01〜P12**: 方法論的「欠陥」の分類。単体実演の欠陥部品 case (実装済):
  - P01 部分最適化 … `case_p01_subsystem`
  - P02 熱統合の後置 … `case_p02_pinch`
  - P04 1次元逐次最適化 … `case_p04_sequential`
  - P05 整数の連続扱い/粗い掃引 … `case_p05_grid`
  - P06 大域性未保証 (単一始点) … `case_p06_multistart`
  - P10 非理想で FUG (精度問題) … `case_p10_fug` (FUG vs HYSYS真値)
  - P12 検証なし (収束未確認) … `case_p12_converge`
  - P03 リサイクル/パージ不変数化、P07 単目的、P08 経済単純化 (NPV無視)、P09 環境/CO2、P11 定常のみ … この sim では未実装/不可 (§7.2)
- 各 case_rep (実在レポート再現) は **これら欠陥部品 P を複数束ねる**。各 case の「含まれる欠陥部品」節に
  ◎(sim再現)/△(検出だが非対応) を明記。
- 各欠陥の損失 = **ΔTAC = TAC(その欠陥手法の best) − TAC(BO best)** [億円/年、%]。これが「定量化」の中身。
- **成熟度スコア** = 12 − (検出された問題数)。低いほど欠陥が多い。case_rep は成熟度 2〜7 を横断。
