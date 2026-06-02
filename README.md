# PDH Simulator — プロパン脱水素プロセス設計・最適化シミュレータ

プロパン脱水素 (Propane Dehydrogenation, PDH) によるプロピレン製造プロセス
全体を、設計・シミュレーション・経済評価・ベイズ最適化まで一気通貫で扱う
Python シミュレータです。

対象プロセスは以下のユニットから構成されます。

- **反応器**: スイング固定床 (軸流) / 径方向流の触媒反応器 (Cr₂O₃-Al₂O₃, Catofin 相当)
- **PSA**: オフガスからの水素分離 (圧力スイング吸着, 活性炭)
- **膜分離器**: プロピレン/プロパン選択膜
- **蒸留塔 3 塔** (Dist1 / Dist2 / Dist3)
- **リサイクル**: 未反応プロパンの循環 (tear ストリームの逐次代入 / Wegstein 加速で収束)
- **熱統合 (HI / HEN)**: ピンチ解析によるユーティリティ削減
- **経済評価**: CAPEX / OPEX / TAC / Profit (長谷部式ベース)
- **ベイズ最適化**: Optuna による全 22 設計変数の制約付き最適化

---

## 1. プロジェクト概要

`main.py` を本体とし、フローシート全体を 1 回評価する関数 (`flowsheet.evaluate`)
を目的関数として、反応器・PSA・膜・原料・蒸留 3 塔の設計変数を同時に最適化します。
反応器は `REACTOR_KIND` で軸流 (21 変数) / 径方向流 (22 変数, 既定) を切り替えます。

蒸留塔は塔ごとにバックエンドを選択できます。

- **Dist1 / Dist3**: SM (学習済み GPR サロゲートモデル, `models/*.pkl`) ─ ほぼ瞬時
- **Dist2**: Aspen HYSYS COM (`main.py` / `exp/exp3.py` の既定) または in-house rigorous solver (Wang-Henke 法)

---

## 2. 動作環境

- **Python 3.13** (`.venv` は 3.13 系で構築。`.python-version` = `3.13.0`)
- **OS**: Windows (HYSYS COM 連携と `.venv` 修復スクリプトが Windows 前提)
- **Aspen HYSYS (COM)**: **Dist2 を HYSYS バックエンドで評価する場合のみ必須** (ライセンス必要)

### HYSYS 無しで動く範囲

HYSYS が無い環境でも、以下の経路は **pure Python** で完結し動作します。

- 反応器・PSA・膜・リサイクル・熱統合・経済評価の全フロー
- 蒸留塔の **SM** (Dist1/Dist3) と **FUG / rigorous** (Wang-Henke) バックエンド

HYSYS が必要なのは Dist2 を HYSYS で解く構成 (`main.py` 既定, `exp/exp3.py`) のみです。
HYSYS 無しで全塔を回したい場合は、全塔 pure Python の **`sub/sub2.py`** (Dist1/Dist3=SM,
Dist2=rigorous, 6 worker 並列対応) や **`sub/sub1.py`** (FUG/rigorous) を使用してください。

---

## 3. セットアップ

```powershell
# 1. venv 作成 (初回のみ)
py -3.13 -m venv .venv

# 2. 依存インストール
.\.venv\Scripts\python.exe -m pip install -r requirements.txt

# 3. (PC を変えた直後など) venv が起動できない場合に修復
.\fix_venv.ps1
```

- **`requirements.txt`** / **`pyproject.toml`**: 依存パッケージ定義 (numpy, scipy,
  thermo/chemicals/fluids ─ PR EOS 物性, tomli, optuna, scikit-learn)。
  物性パッケージはバージョン pin されています (理由は `requirements.txt` 内コメント参照)。
- **`fix_venv.ps1`**: `.venv` をフォルダごと別 PC へ持ち回ると `pyvenv.cfg` に
  焼き込まれた絶対パスが壊れて起動失敗します。これを現マシンの Python 3.13 へ
  向け直す修復スクリプト (site-packages は触らず、完全に可逆)。

---

## 4. 使い方

### 4.1 全変数ベイズ最適化 (本体)

```powershell
.\.venv\Scripts\python.exe main.py > outputs\main_run.log 2>&1
```

- 全 22 変数 (既定 `REACTOR_KIND='radial'`) を Optuna BO で探索 (既定 300 trials)。
- 試行数等は環境変数で短縮可: `PDH_N_TRIALS`, `PDH_N_STARTUP`, `PDH_SEED`, `PDH_SAMPLER`。
- 結果は `outputs/main_<timestamp>/` に出力 (`README.md`, `trials.csv`, `best.json`,
  上位候補の詳細レポート `top*.txt`)。flush 付きログなのでリダイレクトしてもライブに書かれます。
- Dist2 を HYSYS で評価する既定構成では HYSYS ライセンスが必要です。

### 4.2 単一設計点の詳細評価

```powershell
.\.venv\Scripts\python.exe exp\exp3.py
```

- `exp/exp3.py` 冒頭の「実験で振る設計変数」ブロックを編集して 1 つの設計点を
  詳細に評価 (CAPEX/OPEX/spec/HI 内訳)。
- Dist1/Dist3=SM, Dist2=HYSYS。`main.py` の `best.json` をそのまま投入して再現確認にも使えます。

---

## 5. ディレクトリ構成

各ディレクトリの詳細は、配下の `SPEC_*.md` / `ISSUES_*.md` / `CASES.md` 等の
ドキュメント (存在する場合) も参照してください。

| ディレクトリ | 役割 |
|---|---|
| `src/` | コア計算ライブラリ (PR EOS `eos.py`, 反応速度 `kinetics.py`, 触媒 `catalyst_model.py`, 蒸留コア/rigorous/SM, コスト・物性パラメータ, 熱力学) |
| `units/` | 装置ユニットモデル: `reactors/` (swing 軸流・radial 径方向流), `separators/` (PSA・膜・蒸留塔 column1-3), `utils/` (圧縮機・冷却器・ポンプ・膨張弁・混合器), `vle/hysys/` (HYSYS COM アダプタ) |
| `flowsheet/` | フローシート結線・リサイクル収束ソルバ (`solvers/` 逐次代入/Wegstein)・熱統合・経済評価・スペック判定 |
| `optimization/` | Optuna BO インフラ (探索空間・目的関数・制約・ペナルティスケール・並列・コールバック・feasibility 解析・top-k) |
| `simulation/` | 結果表示・実験ランナー (`display.py`, `exp_runner.py`) |
| `stream/` | プロセスストリーム (流量・組成・状態) データ構造 |
| `config/` | 運転設定 (`operating.toml`) とローダ (ソルバ収束パラメータ・ペナルティ・初期推定値) |
| `exp/` | 個別実験スクリプト (`exp3.py` 単一設計点詳細評価, 膜・PSA 感度解析) |
| `sub/` | BO 本体のアーカイブ系: `sub1.py` (FUG/rigorous), `sub2.py` (SM/rigorous/SM, 6 worker 並列。感度解析が import) |
| `tools/` | 補助ツール (モニタリング, プロファイリング, feasibility 実行, notebook 生成, 上位設計抽出) |
| `monitor/` | 解析用 Jupyter ノートブック・図 (反応器転化率上限・圧損/形状・熱統合・LPG 損益分岐 等) |
| `comparing/` | 過去の欠陥的最適化手法を本シミュで再現し ΔTAC で損失を定量化する比較ケース群 (`CASES.md` 参照) |
| `models/` | 学習済み蒸留塔サロゲート (`column1_sm.pkl`, `column3_sm.pkl`) |
| `data/` | 入力データ (触媒 a パラメータ fitting 用 Excel 等) |
| `hysys_cases/` | 段数別の HYSYS ケースファイル (`.hsc`, column1/2/3 各段数) |

---

## 6. 主要な出典・前提

- **設計課題ベース**: プロセスデザイン学生コンテスト Ver.2.0 §4 (蒸留塔・熱交換器・
  圧縮機・加熱炉の設計式)。年間稼働 8000 h/年、生産目標などの規定を採用。
- **コスト式**: 長谷部・外輪 式 (CAPEX/OPEX/TAC 集計, 労務費 C_OL, 保全費等)。
- **物性 (VLE)**: Peng-Robinson 状態方程式 (Peng & Robinson, 1976)。bubble/dew point は
  `thermo` パッケージ (CalebBell, MIT) を内部利用し in-house PR EOS と 0.02% 一致を検証済み。
- **蒸留 rigorous solver**: Wang-Henke bubble-point 法 (Seader, Henley & Roper, Ch.10.4)。
- **触媒**: Cr₂O₃-Al₂O₃ (Catofin 相当)。単価・寿命・嵩密度は文献 citation 付き。
- **膜**: Hua et al. (2024) の実測透過性能 (Q_A, α)。
- **ユーティリティ・原料・製品単価**: 日本実勢値を独自調査 (各定数に出典コメント)。

確定済み出典の一覧は `KNOWN_PLACEHOLDERS.md` §A、各ユニットの `SPEC_*.md` を参照してください。

---

## 7. 仮置き値・既知の課題

経済・物性パラメータの一部や、BO ペナルティ・ソルバ設定の一部は **仮置き値** です。
仮置き・未解決の品質課題の一覧と出典確定状況は **[`KNOWN_PLACEHOLDERS.md`](KNOWN_PLACEHOLDERS.md)**
にまとめてあります。

コード中の仮置き箇所は `!仮置き` マーカーでグレップできます。

```powershell
Select-String -Path .\**\*.py,.\**\*.md -Pattern "!仮置き"
```
