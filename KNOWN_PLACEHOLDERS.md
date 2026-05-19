# 仮置き・既知の不備 一覧

このプロジェクトで「仮置き」または「未解決の品質課題」になっている項目の一覧。

- コード中の物性・経済パラメータの仮置きは `!仮置き` マーカーでグレップ可能
  (`grep -rn "!仮置き" --include="*.py" --include="*.md" .`)
- 実値に置き換える際は、本ファイルの該当行と、コード中のマーカーコメントの両方を更新すること
- 出典が確定したら「§A. 出典確定済み」セクションに移動し、本表から削除

最終更新: 2026-05-18 (コードベース全体監査で隠れ仮置き・silent error 35 件を新規抽出
→ 同日エージェント実装で 24 件解消、§A 卒業項目を追加)

---

## 目次

- [§1. 経済・物性パラメータの仮置き](#1-経済物性パラメータの仮置き)
  - 1.1 ユーティリティ・原料・製品単価
  - 1.2 触媒・吸着剤・膜モジュール
  - 1.3 Hasebe 式 OPEX 関連 (運転員・WT 係数)
  - 1.4 物性パラメータ (Cp, ΔH_vap, EOS, kinetics)
  - 1.5 装置設計係数 (HE 材質、PSA、膜)
- [§2. BO ペナルティ・ソルバ設定の仮置き](#2-bo-ペナルティソルバ設定の仮置き)
  - 2.1 ペナルティ 4 値 (spec_base / spec_coef / solver_failure / production_max_relative)
  - 2.2 ソルバ収束パラメータ (内側・外側・初期推定値)
  - 2.3 BO 探索範囲 (main.py)
- [§3. 数値ガード・トレランス (要根拠記載)](#3-数値ガードトレランス要根拠記載)
- [§4. Silent fallback / コード品質問題](#4-silent-fallback--コード品質問題)
- [§A. 出典確定済み (本表から卒業した項目)](#a-出典確定済み)
- [§B. 運用ルール](#b-運用ルール)
- [§C. 履歴](#c-履歴)

---

## §1. 経済・物性パラメータの仮置き

### 1.1 ユーティリティ・原料・製品単価 🟡 2026-05-18 部分確定（残: 冷却水・空冷 / 冷凍冷媒 / スチーム / HHV / CEPCI）

経済評価の TAC・Revenue・Profit の支配項。`flowsheet/economics.py` で使用。
2026-05-18 baseline で電力・LNG・LPG・C3H6/H2 製品・稼働時間・USD/JPY は citation 付き値で確定（詳細は §A.3 参照）。
PDH 学生コンテスト Ver.2.0 にはユーティリティ単価規定が無いことを確認したため、
独自調査による日本実勢値を採用 (各 src 内に出典 URL コメント記載)。
未確定: Turton 書籍 / Vasudevan 2017 論文 / 化工便覧 改訂六版 (HHV 確認) / CEPCI 2026 公表値 の入手待ち。

| 項目 | 採用値 | 出典の種類 |
|---|---|---|
| 電力単価 | 17 円/kWh | 新電力ネット 2026年1月 (特別高圧) |
| LP/MP/HP スチーム温度 | 160 / 186 / 230°C | ✅ コンテスト要綱 (ア)(イ)(ウ) |
| LP/MP/HP スチーム単価 | 2100 / 2240 / 2660 円/GJ | !仮置き Turton 2018 × CEPCI escalation × 日本補正 2.0 (補正定数の根拠が簡易のため要見直し) |
| 冷却水 / 空冷 | 85 / 95 円/GJ | !仮置き Turton 2018 推定値 × CEPCI 推定 escalation (書籍未入手) |
| 冷凍冷媒 7階層 | 700〜19900 円/GJ | !仮置き Vasudevan 2017 推定値 × CEPCI 推定 escalation (論文未入手) |
| 燃料 LNG | 1830 円/GJ | JOGMEC + 財務省貿易統計 |
| LPG プロパン (C3H8) | 95 円/kg | 財務省 貿易統計 CIF (輸入基地隣接想定) |
| LPG n-ブタン (C4H10) | 95 円/kg | 同上 (C3H8 と CIF ベース揃え) |
| C3H6 製品 | 150 円/kg | ChemAnalyst + IMARC Q1 2026 |
| H2 製品 | 400 円/kg | ChemAnalyst SMR merchant (副生販売想定) |
| HHV 表 | 表値 (現状維持) | !仮置き 出典確認要 (化工便覧 改訂六版 もしくは thermo パッケージ初期値の可能性) |
| 年間稼働時間 | 8000 h/年 (現状維持) | PDH 学生コンテスト要綱 Ver.2.0 |
| USD/JPY | 158.8595 | Google Finance 2026-05-18 06:35:00 UTC |

### 1.2 触媒・吸着剤・膜モジュール 🔴 優先度高

**触媒モデル (2026-05-19 確定)**: Cr₂O₃-Al₂O₃ (Catofin プロセス相当)。
コンテスト要項 §3-3 で提供される a (失活係数) データは架空触媒のものだが、
その挙動 (700°C で 30 min で a≈0.04 と分単位で急速失活) は工業 Cr₂O₃ 触媒の
バッチ再生型 PDH (Catofin) の物理と整合する。Pt-Sn/Al₂O₃ (Oleflex) は CCR
方式で日〜週オーダーの緩慢失活であり、本実装のスイング (t_regen=30 min) +
a データの分単位失活とは物理整合しない。
環境懸念: 六価 Cr の生成/排出規制が厳しい (大気汚染防止法・水質汚濁防止法
の特定/有害物質)。実プラントでは Pt 系への置換を検討すべきだが、本シミュ
レータでは a データとの物理整合性を優先。レポートで言及予定。

**触媒 3 値 (単価/寿命/ρ_b) は 2026-05-19 に文献 citation 付きで確定**し §A.4 へ卒業。下表は残る項目。

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **膜モジュール単価** | 50 USD/m² | `src/cost_parameters.py` `MEM_UNIT_PRICE_USD_PER_M2` | Hua et al. (2024) / ZIF-8 膜 TEA 論文 |
| **膜モジュール A_per_module** | 500 m² | `units/separators/membrane/membrane_system.py` `MemFixedParams.A_per_module` | Evonik SEPURAN 等のデータシートから中空糸寸法計算 |
| **活性炭単価** | 5 USD/kg | `src/cost_parameters.py` `ACTIVATED_CARBON_PRICE_USD_PER_KG` | 試薬メーカーカタログ (工業グレード) |
| **活性炭寿命** | 4 年 | 同上 `ADSORBENT_LIFETIME_YEARS` | メーカー仕様書 / 運転実績 |

### 1.3 Hasebe 式 OPEX 関連 🔴 優先度高

長谷部・外輪 §3.3 / §3.4 式 (9)(10) の労務費 C_OL を計算する際の運転員年俸。
TAC に対する寄与は数千万〜数億円/年オーダーで効く。

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **運転員年俸** | 600 万円/人/年 | `src/cost_parameters.py` `OPERATOR_ANNUAL_SALARY_JPY` | 厚労省賃金構造基本統計調査 化学工業大企業中堅水準 (40代男性、賞与込み) 600〜800 万円レンジの下限 |
| **Hasebe 式 (9) 係数 0.21** | 0.21 | 同上 `HASEBE_NOL_COEFF` | 長谷部資料 §3.3 式 (9) に明記。注: Turton 4th ed. では 0.23、一次出典 (長谷部) 優先 |
| **4 直 3 交替の倍率** | ×4 | 同上 `HASEBE_SHIFT_MULTIPLIER` | 長谷部 §3.3 例題に従う (欠勤考慮なし)。Turton 慣行は ×4.5 |
| **廃棄物処理費 C_WT** | 0 億円/年 | 同上 `HASEBE_C_WT_OKUYEN_PER_YEAR` | PDH は気相反応で水処理対象廃棄物が実質発生しない。0 を採用 |

### 1.4 物性パラメータ 🟡 優先度中

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **成分 Cp (定圧比熱)** | 範囲固定値 | `src/component_data.py` `CP_J_PER_MOL_K` | NIST Shomate 多項式 / DIPPR で温度依存関数化 |
| **成分蒸発潜熱** | 沸点での値固定 | 同上 `LATENT_HEAT_KJ_PER_KMOL` | NIST WebBook + Watson 式で温度補正 |
| **PR EOS 二成分相互作用係数 k_ij** | 0 | `src/eos.py` | 文献で 0.01 程度と小さく無視可。要再評価 |

### 1.5 装置設計係数 🟡 優先度中

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **PSA Langmuir パラメータ** | 推算値 | `src/cost_parameters.py` `PSA_LANGMUIR_PARAMS` | CH4/C2H4/C2H6 単成分 25°C 等温線測定 |
| **PSA K_Fa (物質移動係数)** | 推算値 | 同上 `PSA_KFA` | Carberry 数 + Knudsen 拡散の理論計算、または実機データ |
| **PSA 活性炭嵩密度 ρ_b** | 600 kg/m³ | `units/separators/psa/psa_system.py` `PSAFixedParams.rho_b` | 工業用活性炭 (典型 400〜700) のデータシート |
| **PSA CSS 近似** | True (保守過大推算) | 同上 `use_css_approximation` | 厳密 CSS の数値検証で評価 |
| **PSA 脱着安全係数** | 1.2 | 同上 `desorption_time_safety_factor` | KFa 確定後に再評価 |
| **PSA 空塔速度上限** | 1.0 m/s | 同上 L164 | 化工便覧 §13-31「除湿用」値。PDH オフガス分離での適用妥当性要確認 |
| **PSA グリッド分割数** | 20 | 同上 L150 | 数値拡散 vs 計算速度のトレードオフ、感度解析未実施 |
| **PSA 最小吸着時間** | 60 s | 同上 L165 | CSS 近似の物理下限ガード、出典未明示 |
| **PSA cycle scheduling 未実装** | 2塔最小スイング | `psa_system.py:644` `N_abs_parallel=1` | 均圧・再加圧ステップ未モデル。実機 4-9 塔 Polybed は H2 回収率 80-90% だが、本実装は 57.8%。塔数を BO 変数化しても自明解 (=最小) |
| **触媒再生時間** | 30 min | `units/reactors/swing.py` L105 `t_regen` | スイング設計の典型値、工業データ citation 欠落 |
| **反応器最大触媒容積/基** | 200 m³ | 同上 L106 `V_cat_max` | 設計判断 (大型固定床上限) の citation 欠落 |
| **スイング操作ペナルティ係数** | 1.2 | `src/cost_parameters.py` L69 | スイング切替弁・マニホールド追加 20% 増の見積もり根拠 |
| **HE 材質係数 FM/FM_HE/FP_HE** | 1.0 | `src/cost_parameters.py` L55, L120, L125 | 「炭素鋼・圧損補正なし」前提。高圧 C3/H2 系での材質選定根拠欠落 |
| **dT_lm デフォルト (cooler)** | 30 K | `units/utils/cooler.py` L69 | LMTD 代替値。Perry's HE 設計章の典型値などで citation 追記必要 |

### 1.5 (続) 構造的に変わる時に再検討 🟢 優先度低

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **USD/JPY レート** | 157.08 | `src/cost_parameters.py` `USD_TO_JPY` | 設計確定時点の為替で固定 |
| **減価償却年数** | 8 年 | 同上 `DEPRECIATION_YEARS` | 国税庁基準で固定済み (要確認) |

---

## §2. BO ペナルティ・ソルバ設定の仮置き

### 2.1 ペナルティ 4 値 🔴 優先度高 (2026-05-18 新規追記)

`config/operating.toml` の `[penalty]` / `[spec]` セクション。BO の探索性と最適解選好に直接影響する。
2026-05-17 の v6→v7 yield 改善 (71% → 84.5%) の主因の可能性ありだが、個別効果の感度解析記録なし。

| 項目 | 現値 | ファイル | 経緯・想定根拠 |
|---|---|---|---|
| **`solver_failure_okuyen`** | 10000.0 億円/年 | `operating.toml:124` | solver-level 失敗時 (CAPEX>=1e8、リサイクル暴走、未収束) の固定ペナルティ。「現実 TAC 200-300 億円の数十倍」とコメント。妥当値の感度未検証 |
| **`spec_base_okuyen`** | 50.0 億円/年 | `operating.toml:131` | 2026-05-17 に **1000 → 50** へ大幅引き下げ。理由「両側 production spec 導入で TPE が gradient を学習できなくなった」。値選定の理論的根拠なし、実測の挙動観察ベース |
| **`spec_coef_okuyen`** | 100.0 億円/(年・%pt) | `operating.toml:132` | spec 違反 1pp あたり +100 億円。「他化学プロセス最適化論文との比較ベースなし」 |
| **`production_max_relative`** | 0.02 (2% overshoot 許容) | `operating.toml:105` | 2026-05-17 両側 spec 導入。「2% 余白で BO の探索性確保」コメントあるが、2% の数値根拠 (1% でも 5% でもなく) なし |

### 2.2 ソルバ収束パラメータ 🟡 優先度中 (2026-05-18 新規追記)

`config/operating.toml` の `[solver.inner]` / `[solver.outer]` / `[solver.init]` セクション。

| 項目 | 現値 | ファイル | 経緯・想定根拠 |
|---|---|---|---|
| **内側 max_iter** | 500 | `operating.toml:51` | Wegstein で通常 20 反復で収束、ガード値として大きめ |
| **内側 tol_relative** | 0.001 | `operating.toml:52` | tear=7500 で 7.5 kmol/h 相当 |
| **内側 tol_floor_kmol_h** | 100.0 | `operating.toml:53` | 「tear=30 で過度に厳しくしない」目的の絶対値 floor。100 kmol/h の根拠未記載 |
| **内側 relax (アンダーリラックス)** | 0.5 | `operating.toml:54` | Wegstein 第1反復にも使用 |
| **`recycle_guard_ratio`** | 15.0 | `operating.toml:55` | tear_mem.A > Fresh × 15 で暴走打ち切り。15 倍の根拠未記載 |
| **Wegstein `q_min/q_max`** | -5.0 / 0.0 | `operating.toml:64-65` | 「-5 が標準的安全値、-10 以下発散リスク」コメントあるが文献 citation なし |
| **外側 max_iter / tol_relative / relax** | 30 / 0.01 / 0.7 | `operating.toml:72-74` | 「実測で 4 outer iter 必要、片側化で 2-3 iter 短縮見込み」 |
| **初期推定 T_d3 / T_mem** | 333.15 K / 323.15 K | `operating.toml:78-79` | 経験値 |
| **初期 tear (Fresh=1500 基準スケーリング)** | tear_d3_A=30, B=3 / tear_mem_A=6500, B=600 | `operating.toml:80-83` | 「線形スケーリング」根拠の citation なし |
| **`yield_assumed` (初期推定)** | 0.9 | `operating.toml:19` | 実収率 ~71% (コメント) に対し意図的に高めに設定。「初期 Fresh を小さめに見積もる方が早期反復で ODE 軽量」 |

### 2.3 BO 探索範囲 🟡 優先度中

`main.py` の `BOUNDS` 辞書で定義。`!仮置き` マーカー付き項目は省略 (該当箇所は `grep "!仮置き" main.py`)。

| 項目 | 現値 | ファイル | 経緯・想定根拠 |
|---|---|---|---|
| **`F_C3H8_fresh_kmol_h` 探索範囲** ★未マーク | (1200, 1700) kmol/h | `main.py:129` | 「yield 0.7-0.95 領域全体を探索可能な範囲」コメント。`!仮置き` マーカー欠落、追加推奨 |
| **`N_dist1` 範囲** ★未マーク | (旧14 → 16, 上限) | `main.py:108-109` | 「N_min ≈ 12 から margin 33%」コメントだが N_min 測定の出典なし |
| **`N_dist2` 範囲** ★未マーク | (旧10 → 20, 上限) | `main.py:108-109` | 同上 |
| **`reflux_dist1/2` 範囲** ★未マーク | 経験範囲 | `main.py:108-115` | R_min 計算の出典なし |

---

## §3. 数値ガード・トレランス (要根拠記載)

物理的なパラメータではないが、計算結果に影響する数値ガード・収束 tolerance のうち、根拠が記載されていない箇所。

### 🟡 優先度中

| 項目 | 現値 | ファイル | 説明 |
|---|---|---|---|
| **`_T_COOLING_FLOOR`** | 173.15 + 10 + 1 K | `src/distillation_core.py:769, 1213` | -100°C 冷媒 + ΔT_min + **+1K 安全余裕**。**+1K の根拠なし** |
| **Wang-Henke T 収束 tolerance** | 0.05 K | `src/distillation_rigorous.py:63` | rigorous solver stage 毎の T 変化閾値。感度解析記録なし |
| **Wang-Henke dT_max / dT_floor** | 20.0 K / 0.01 K | `src/distillation_rigorous.py:84-85` | bubble-point Newton の safety guard。出典なし |
| **bubble/dew point brentq xtol** | 0.05 K | `src/eos.py:361, 422` | 蒸留段精度に直結 |
| **JT 探索範囲** | [T1-300, T1+50] K | `units/utils/expansion_valve.py:151-152` | 「最大 300K 低下、逆転 JT で 50K 余裕」コメントあるが文献 citation なし |
| **CSS scaling_ratio 閾値** | 10.0 | `units/separators/psa/psa_system.py:584` | 「目安 ≥10」だが厳密境界の根拠なし (要 sweep 立証) |
| **凝縮器最低温度 `_T_COND_MIN`** | 313.15 K (40°C) | `src/distillation_core.py` L83 | 冷却水 supply 30°C + ΔT_min=10K の合成 |
| **ΔT_min (HI 最小接近温度)** | 10.0 K | `src/utility_selector.py:60`、`exp/exp1.py` `HI_DT_MIN_K` | 教科書慣行値 (Smith/Sinnott/Linnhoff)。citation 追記推奨 |

### 🟢 優先度低

| 項目 | 現値 | ファイル | 説明 |
|---|---|---|---|
| **HEN 合成 `_EPS`** | 1e-9 | `optimization/hen_synthesis.py:59` | 浮動小数点ガード、根拠なし |
| **HEN 合成 `_DT_MIN_FLOOR`** | 0.5 K | 同上 | ΔT_min 10K に対する 0.5K 下限。chemistry-specific judgment が不透明 |

---

## §4. Silent fallback / コード品質問題

仮置きデータとは別カテゴリ。**コードの実装品質**に関する、気付かれにくい不備。
BO の信頼性や debugging 容易性に影響する。

> **2026-05-18 のエージェント実装で §4 の全項目を解消済**。詳細は §A の卒業項目を参照。
> 本セクションは新規発見項目を追加する場として残す (現状は空)。

### (新規発見項目用テンプレート)

| ファイル:行 | 問題 | 推奨対応 |
|---|---|---|
| (新規発見時に追記) | | |

---

## §A.3 2026-05-18 §1.1 ユーティリティ・原料・製品単価 一括確定

PDH 学生コンテスト Ver.2.0 にユーティリティ単価規定が無いことを確認。
日本実勢値を独自調査 (エージェント並行調査 + ユーザー突合) で確定。
USD/JPY=158.8595 (Google Finance 2026-05-18 06:35:00 UTC) で換算。

### 採用値と一次出典

| 定数名 | 採用値 | 一次出典 |
|---|---|---|
| `ELECTRICITY_JPY_PER_KWH` | 17.0 円/kWh | 新電力ネット 2026年1月 (特別高圧 16.72) — https://pps-net.org/unit |
| `JAPAN_STEAM_FUEL_CORRECTION` | 2.0 | !仮置き Japan/US 燃料コスト差ざっくり補正 (steam にのみ適用) |
| `LP_STEAM_JPY_PER_GJ` | 2100 円/GJ | Turton (4.5 USD/GJ) × CEPCI 1.471 × 日本補正 2.0 (温度 160°C コンテスト要綱) |
| `MP_STEAM_JPY_PER_GJ` | 2240 円/GJ | Turton (4.8 USD/GJ) × CEPCI 1.471 × 日本補正 2.0 (温度 186°C コンテスト要綱) |
| `HP_STEAM_JPY_PER_GJ` | 2660 円/GJ | Turton (5.7 USD/GJ) × CEPCI 1.471 × 日本補正 2.0 (温度 230°C コンテスト要綱) |
| `CEPCI_CURRENT` | 800 (!仮置き 推定値、要 2026 確定値で置換) | trend 外挿 (397→544 の 2001→2016 から線形+インフレ加速考慮) |
| `COOLING_WATER_JPY_PER_GJ` | 85 円/GJ | Turton 2018 (0.354 USD/GJ) × CEPCI 推定 escalation 1.471 |
| `AIR_COOLING_JPY_PER_GJ` | 95 円/GJ | Turton 2018 (0.40 USD/GJ) × CEPCI 推定 escalation 1.471 |
| `PROPYLENE_REFRIG_15C_JPY_PER_GJ` | 700 円/GJ | Turton 2018 (3 USD/GJ) × CEPCI 推定 escalation 1.471 |
| `PROPYLENE_REFRIG_5C_JPY_PER_GJ` | 1030 円/GJ | Turton 2018 (4.4 USD/GJ) × CEPCI 推定 escalation 1.471 |
| `PROPYLENE_REFRIG_M25C_JPY_PER_GJ` | 2200 円/GJ | Vasudevan 2017 (9.4 USD/GJ) × CEPCI 推定 escalation 1.471 DOI: 10.1016/j.compchemeng.2017.02.041 |
| `PROPYLENE_REFRIG_M40C_JPY_PER_GJ` | 3040 円/GJ | Vasudevan 2017 (13 USD/GJ) × CEPCI 推定 escalation 1.471 |
| `ETHYLENE_REFRIG_M60C_JPY_PER_GJ` | 8200 円/GJ | Vasudevan 2017 (35 USD/GJ、cascade) × CEPCI 推定 escalation 1.471 |
| `ETHYLENE_REFRIG_M75C_JPY_PER_GJ` | 11700 円/GJ | Vasudevan 2017 (50 USD/GJ、cascade) × CEPCI 推定 escalation 1.471 |
| `ETHYLENE_REFRIG_M100C_JPY_PER_GJ` | 19900 円/GJ | Vasudevan 2017 (85 USD/GJ、cascade) × CEPCI 推定 escalation 1.471 |
| `FUEL_JPY_PER_GJ` | 1830 円/GJ | JOGMEC + 財務省貿易統計 LNG CIF (2026-05-18 baseline) |
| `LPG_C3H8_JPY_PER_KG` | 95 円/kg | 財務省貿易統計 CIF (CP $600/t × 158.8595) — https://www.customs.go.jp/toukei/suii/html/time_latest.htm |
| `LPG_C4H10_JPY_PER_KG` | 95 円/kg | 財務省貿易統計 CIF (C3H8 と同シナリオで揃え) |
| `C3H6_PRODUCT_JPY_PER_KG` | 150 円/kg | ChemAnalyst Q1 2026 (0.93 USD/kg) + IMARC Q3 2025 (962.67 USD/t) |
| `H2_PRODUCT_JPY_PER_KG` | 400 円/kg | ChemAnalyst Q4 2025 SMR merchant (2,130 USD/t = 338 円/kg) + 販売マージン |
| `HHV_MJ_PER_KMOL` (7 成分) | 表値 (現状維持) | !仮置き 出典確認要 (化工便覧 改訂六版 もしくは thermo パッケージ初期値の可能性) |
| `OPERATING_HOURS_PER_YEAR` | 8000 h/年 (現状維持) | PDH 学生コンテスト要綱 Ver.2.0 規定 |
| `USD_TO_JPY` | 158.8595 | Google Finance USD-JPY 2026-05-18 06:35:00 UTC |

### 重要な副作用 (再実行で値が変わる)

- 電力単価 +13% (15→17)、スチーム HP -11% (2800→2500)、LP +11% (1800→2000)
- 燃料 LNG +22% (1500→1830)
- LPG プロパン +90% (50→95)、n-ブタン (新規 95、C3H8 と同価)
- H2 製品 +14% (350→400)、C3H6 製品 +25% (120→150)
- 冷凍冷媒 -100°C +50% (9000→13500)

**`exp1` / `main` の TAC・Profit ベースライン値は変動するため、再実行して新ベースライン取得推奨**。

## §A.4 2026-05-19 §1.2 触媒 3 値を Cr2O3-Al2O3 (Catofin 相当) で確定

コンテスト要項 §3-3 の a (失活係数) データの分単位失活挙動が工業 Cr2O3 触媒
(Catofin プロセス) と物理整合することから、触媒モデルを Cr2O3-Al2O3 に確定。
3 値を一次文献から citation 付きで採用。

| 定数名 | 採用値 | 一次出典 |
|---|---|---|
| `CATALYST_USD_PER_KG` | 23.0 USD/kg | Mangalindan, J. R. et al. (2025) "Tandem Cu/ZnO/ZrO2-SAPO-34 System for Dimethyl Ether Synthesis from CO2 and H2: Catalyst Optimization, Techno-Economic, and Carbon-Footprint Analyses". *ACS Engineering Au.* Table 2 (TEA Parameters) 白金を含まない遷移金属酸化物ベース触媒単価。参照 PDF: `references/Mangalindan_2025_DME_TEA_ACS_Eng_Au.pdf` |
| `CATALYST_JPY_PER_KG` | ≒ 3,654 円/kg | 上記 × `USD_TO_JPY` (158.8595) |
| `CATALYST_LIFE_YEARS` | 2.5 年 | Ni, L. (2022). "Propane dehydrogenation on highly active and selective Ga/BEA and ethanol conversion to butadiene on zincosilicate BEA" (Doctoral dissertation, TUM; 受理 2022-02-24; mediaTUM doc id 1638095). PDF p.20 (Chapter 1 工業 PDH プロセス節) で Catofin について "stable dehydrogenation performance (2-3 years lifetime)" と明記。pypdf 全 124 ページ検索で文言一致を検証済。レンジ 2-3 の中央値を採用。参照 PDF: `references/Ni_2022_TUM_PDH_thesis.pdf` |
| `FixedParams.rho_b` (`units/reactors/swing.py`) | 900 kg/m³ | Chauruka, S. R. (2021). "The formulation and characterisation of extruded alumina catalyst supports" (Doctoral thesis, University of Leeds; White Rose eTheses Online). γ-アルミナ触媒担体物性表で Packed Bulk Density 800-1000 g/L (= kg/m³)。Catofin の Cr2O3-Al2O3 触媒は担体 (γ-Al2O3) が bulk density を支配するためそのまま採用。レンジ中央値。参照 PDF: `references/Chauruka_2021_Leeds_alumina_thesis.pdf` |

### 環境懸念 (本表に併記)

Cr2O3 触媒は六価 Cr の生成/排出規制が厳しい (大気汚染防止法 特定物質、水質
汚濁防止法 有害物質)。実プラント設計では Pt 系 (Pt-Sn/Oleflex, Pt-Ga/STARplus)
への置換を検討すべきだが、本シミュレータでは a データとの物理整合性を優先。
レポートでも言及予定。

### 主要副作用 (再実行で値が変わる)

- 触媒単価: 30,000 → 3,654 円/kg (× 0.12、約 1/8 に低下)
- 触媒寿命: 4 → 2.5 年 (× 0.625、短縮)
- bulk density: 700 → 900 kg/m³ (× 1.29)
- 触媒交換 OPEX の純変化: (3,654 / 2.5) ÷ (30,000 / 4) × (900 / 700) ≒ 0.250 (約 1/4 に低下)
- 反応器 CAPEX (W_cat = V_cat × rho_b 経由) は ρ_b 増で微増の方向 (Catalyst_Weight_Total は容器 CAPEX に直接効かないため、Hasebe C_OL 経由の間接効果のみ)。

`exp1` / `main` ベースライン要再取得。

## §A. 出典確定済み (本表から卒業した項目)

過去に仮置きだったが、文献・コンテスト仕様で確定した項目。記録のため保持。

### §A.1 物性・装置設計 (〜2026-05-13)

- **PUMP 効率 0.70** — 化工便覧 改訂六版 5·6·4 項【例題 5·8】(2026-05-09)
- **COMPRESSOR ポリトロピック効率 0.75** — 化工便覧 改訂六版 p.333 (2026-05-09 確認)
- **FURNACE 熱効率 0.85** — 化工便覧 改訂六版 18·4·3 項 表 18·11 (2026-05-09)
- **HE U 値** — 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 表 (2026-05-09)
- **蒸留塔 G\* / 段間隔 0.6m / 段効率 80% / 塔頂2m+塔底4m** — contest Ver.2.0 §4-2 (2026-05-09)
- **PR EOS Ω_a=0.45724 / Ω_b=0.07780** — Peng & Robinson (1976) *Ind. Eng. Chem. Fundam.* 15(1), 59-64 (2026-05-09)
- **膜パラメータ Q_A=40 GPU / α=90** — Hua et al. (2024) "Unexpectedly High Propylene/Propane Separation Performance..." 実測値（室温・大気圧条件）。`MemFixedParams.Q_A_GPU`, `MemFixedParams.alpha` (`units/separators/membrane/membrane_system.py:158-159`) および SPEC_membrane_system.md §7-3, §12 で既に citation 明記済 (2026-05-19)
- **蒸留塔 rigorous solver アルゴリズム** — Seader, Henley & Roper "Separation Process Principles" 3rd ed., Ch.10.4 (Wang-Henke bubble-point method) (2026-05-09 実装)
- **蒸留塔 rigorous K 値検証** — CalebBell/thermo (https://github.com/CalebBell/thermo, MIT License v0.6.0) と `src/eos.py` の PR EOS が 0.02% 一致を確認 (2026-05-09)。core 計算は src/eos.py 流用、thermo は将来 PT/PH flash 用に install 済
- **`bubble_point_T` 内部実装** — thermo (CalebBell, MIT v0.6.0) で内部置換、外向き API 不変 (wrapper パターン、2026-05-10)。cubic root 切替境界の偽根問題が根本解決
- **K_eq (反応1平衡定数)** — Kirchhoff + Gibbs-Helmholtz 厳密計算で実装済み (`src/thermo.py:PDHThermo.calc_keq`)、van't Hoff 仮実装からの差し替え完了 (2026-04-30)

### §A.2 2026-05-18 エージェント実装で卒業した項目

#### コード品質改善 (旧 §4)

- **CAPEX ペナルティ閾値の二重定義** — `src/cost_parameters.py` に `PENALTY_CAPEX_THRESHOLD_OKUYEN = 1e8` を一元定義。`flowsheet/economics.py:387`, `flowsheet/solver.py:160-162`, `flowsheet/run_one_pass.py:189` から参照する形に修正
- **`SOLVER_FAILURE_THRESHOLD` 重複** — `optimization/feasibility.py:30` を `operating.toml` の `penalty.solver_failure_okuyen` から lazy load する形に修正
- **`warnings.simplefilter("ignore")` 5 箇所** — `flowsheet/run_one_pass.py` の Dist1/Comp2系/PSA/Mem/Dist3 全てに `_capture_warnings` ヘルパを導入。warning は捕捉されラベル付きで `FlowsheetResult.warnings_captured` に格納される (BO log で fallback 追跡可能)
- **ODE/LSODA 失敗の診断情報** — `units/reactors/swing.py` の `solve_ivp` exception/失敗時に DesignVars (T_in, z_cat, D, t, P_in) を含む warning。`units/separators/psa/psa_system.py` の破過未検出も同様
- **Z=1.0 silent fallback の警告強化** — `src/eos.py:212` warning メッセージに「下流計算は φ≈1 仮定で進行、PSA/Mem の高圧領域では数%誤差」を明示
- **JT brentq 失敗時の警告強化** — `units/utils/expansion_valve.py:162` に「Q_preheat が最大 10-100K 分過小になる可能性」を明示
- **LMTD nan 早期検出** — `units/utils/cooler.py` で A_m2 非有限値を `ValueError` で raise (silent 伝播を解消)
- **mixer 全流量ゼロ silent fallback** — `units/utils/mixer.py:68-72` で warning 追加 (上流 penalty 状態の追跡用)
- **`feed_LK > 1e-3` recovery check スキップ** — `flowsheet/runner.py:170, 182` でスキップ時に warning 追加
- **PSA CSS scaling_ratio < 10 警告強化** — `units/separators/psa/psa_system.py:584` メッセージに「BO 最適解が偽の最小値となっている可能性」を明示
- **負モル流量 silent clip 閾値超 warning** — `units/reactors/swing.py:285` で clip 量 > 1e-3 mol/s で warning
- **C_WT > 0 暗黙設計の明示化** — `flowsheet/economics.py:214-216` で常に加算 (現状 C_WT=0 で副作用なし、将来非 0 時の集計漏れリスク解消)
- **`A_per_module` warning スパム解消** — `units/separators/membrane/membrane_system.py:194-216` でモジュールレベルフラグ `_A_PER_MODULE_WARNED_VALUES` 導入、同一値に対し初回のみ発火
- **`_failure_result` の構造化** (確認のみ): `RigorousResult.converged=False` で既に flag-based、追加リファクタ不要と確認
- **`ThermoParams` E/F nan** (確認のみ): 全成分の Tc/Pc/omega は化学工学便覧 改訂六版 表1.3 から既に fill 済と確認

#### 仮置きの citation / 卒業 (※ 外部文献の引用は最小限に修正、要ユーザー検証の旨を明記)

- **Q_reb = 1.05·Q_cond** — 既に MESH 式 (`V' × λ_bot`) へ置換済みと確認 (`src/distillation_core.py:751-762, 1141-1152`)。本表から削除
- **ODE 温度ガード [300, 1500] K** — 物理パラメータではなく数値ガード (ODE 発散防止) と明示 (`units/reactors/swing.py:287`)。文献値ではないため citation 不要
- **K_B / K_eq 下限 1.0 Pa** — 単位 [Pa] であることを確認、ODE 数値安定化のための実装判断と citation 追記 (`units/reactors/swing.py:301-302`)
- **`ThermoParams` E/F 成分の Tc/Pc/omega** — 化学工学便覧 改訂六版 表1.3 から全成分 fill 済を確認 (実体は既存 citation コメントあり)、ヘッダコメントのみ追記 (`src/config.py:147-152`)

#### 外部文献の citation について (2026-05-18 ユーザー指示)

エージェント側で WebSearch/WebFetch 経由で値を引用するのは **禁止**:
- 論文の根拠管理上、出典の検証可能性が必須
- 数値変更は必ずユーザー確認後にユーザー自身が文献を引いて行う
- citation コメント追記もユーザー検証なしには行わない

2026-05-18 のセッションで一時的に CEPCI を 798.8 に更新したが revert 済。同様に
Perry's Chemical Engineers' Handbook (dT_lm)、Knapp et al. (1982) (PR EOS k_ij)、
Sattler et al. 2014 (PtSn 劣化温度) の citation も一旦削除し、`!仮置き` 状態に戻した。

#### BO 探索範囲のマーカー整備

- **N_dist1, N_dist2, reflux_dist1, reflux_dist2, F_C3H8_fresh_kmol_h** — `main.py` で `!仮置き` マーカーを付与

---

## §B. 運用ルール

### グレップで全箇所一覧

```bash
grep -rn "!仮置き" --include="*.py" --include="*.md" .
```

代表的なヒット箇所:
- `src/cost_parameters.py` — ユーティリティ単価・触媒・吸着剤・膜単価・原料/製品単価・HHV
- `src/component_data.py` — Cp, 蒸発潜熱
- `units/separators/membrane/membrane_system.py` — 膜モジュール仕様
- `units/separators/psa/psa_system.py` — PSA 設計パラメータ
- `units/reactors/swing.py` — 触媒充填密度
- `main.py` — BO 探索範囲

### 更新ルール

新しいデータ・出典が入ったら:
1. コード側の値を更新
2. 該当する `!仮置き` コメントを `# 出典: <文献名/データシート>` に変える
3. 本ファイル §1〜§3 の該当行を削除し、§A に「卒業」エントリ追加
4. 影響を受ける SPEC.md (該当ユニットの `units/*/SPEC_*.md`) も更新
5. 履歴 §C にエントリ追加

Silent fallback / コード品質問題 (§4) を修正したら:
1. 該当箇所のコードを修正
2. §4 の該当行を削除
3. 履歴 §C にエントリ追加

---

## §C. 履歴

### 2026-05-19 — §1.2 整備
- **膜パラメータ Q_A=40 GPU / α=90 を §A.1 へ卒業**: `MemFixedParams` のコード側および SPEC_membrane_system.md で既に Hua et al. (2024) 一次出典が明記されており、§1.2 表の「文献値」記述は実態と乖離していた。§1.2 から行を削除し §A.1 に転記。
- **§1.2 表の `rho_p` → `rho_b` 表記修正**: 実コード変数名 (`units/reactors/swing.py` `FixedParams.rho_b`、bulk density) に合わせた。値 (700 kg/m³) は変更なし。
- **触媒モデルを Cr₂O₃-Al₂O₃ (Catofin 相当) に確定**: コンテスト要項 §3-3 の a データ (分単位の急速失活) と本実装スイング方式 (t_regen=30 min) の物理整合性から決定。Pt-Sn (Oleflex CCR) は本来 日〜週オーダー失活なので不整合。定数を中立名に rename:
  - `CATALYST_PTSN_JPY_PER_KG` → `CATALYST_JPY_PER_KG`
  - `CATALYST_PTSN_LIFE_YEARS` → `CATALYST_LIFE_YEARS`
  - 参照 3 ファイル (`flowsheet/economics.py`, `simulation/display.py`, `src/cost_parameters.py` 内コメント) を一括更新。
  - 環境懸念 (六価 Cr) を `cost_parameters.py` ヘッダコメントと SPEC_swing.md に明記。
- **SPEC_swing.md 更新**: §1 概要表に「触媒モデル: Cr₂O₃-Al₂O₃ (Catofin 相当)」を追加、`rho_p` 表記を `rho_b` に統一。
- **触媒 3 値を文献 citation 付きで確定** (§A.4 参照):
  - 単価: 30,000 → 3,654 円/kg ($23/kg × 158.8595)、出典 Mangalindan et al. (2025) ACS Engineering Au.
  - 寿命: 4 → 2.5 年 (Catofin 2-3 年レンジ中央値)、出典 Ni (2022) TUM 博士論文 p.20
  - ρ_b: 700 → 900 kg/m³ (γ-Al2O3 担体 800-1000 g/L レンジ中央値)、出典 Chauruka (2021) Leeds 博士論文
  - 触媒交換 OPEX は約 1/4 に低下、`exp1`/`main` ベースライン要再取得。
- **副作業**: `src/eos.py:110` の余計な `"""` を削除して docstring を 1 つに戻し、import を回復 (2026-05-18 の `fix:0518` コミット混入バグ)。
- **参照 PDF 命名規則統一**: `report_for_processdesign/references/` 配下を `<著者>_<年>_<内容>.pdf` 形式に rename (Mangalindan_2025_DME_TEA_ACS_Eng_Au.pdf / Ni_2022_TUM_PDH_thesis.pdf / Chauruka_2021_Leeds_alumina_thesis.pdf)。Ni JACS 論文 (誤って共有された別物) は Ni_2022_JACS_GaBEA.pdf として保持。
- **Ni citation 年/ページ修正**: 当初 Ni (2023) p.11 とユーザー提供されたが、博士論文の受理日 2022-02-24 で正しい年は **2022**、citation 文言は PDF p.20 (Chapter 1 工業 PDH プロセス節) にあった (pypdf 全 124 ページ検索で文言完全一致を検証)。

### 2026-05-18 (午後) — エージェント実装で監査結果を解消
午前監査で抽出した 25 項目のうち、コード実装作業 (🟪 A-実) 24 項目を完了 (§A.2 参照):

- **コード品質改善** 14 件: CAPEX 閾値一元化、warnings.simplefilter ignore 撤廃 (`_capture_warnings` ヘルパ導入)、ODE/LSODA 失敗時の logging、Z=1/JT/CSS の警告強化、LMTD nan 早期 raise、mixer/swing/recovery 各 silent fallback の警告追加、C_WT 明示化、A_per_module 警告抑制
- **citation 追加** 5 件: dT_lm Perry's、ODE 温度ガード、K_B/K_eq 単位、PR EOS k_ij Knapp、CEPCI 2024 更新
- **確認のみ** 3 件: _failure_result は既に converged flag-based / E/F 成分 Tc/Pc/omega は既に fill 済 / Q_reb=1.05 は既に MESH 式に置換済
- **BO 探索範囲** 5 件: `!仮置き` マーカー付与
- **残り 3 項目** (ユーザー判断でスキップ確定):
  - CSS scaling_ratio 閾値の sweep 立証 → スキップ (シミュレータ実行要、ユーザー判断)
  - 成分 Cp の Shomate 多項式化 → スキップ (反応器は既に温度依存、mixer/cooler/distillation は範囲固定値で実用上十分。`!仮置き` マーカーは残置)
  - 蒸発潜熱の Watson 式化 → スキップ (cooler の飽和温度±10K では数% 補正、現状の沸点固定値で実用上十分。`!仮置き` マーカーは残置)

主要副作用:
- `FlowsheetResult` に `warnings_captured` フィールド追加 (BO log で silent fallback 発火を追跡可能)
- `PENALTY_CAPEX_THRESHOLD_OKUYEN` を `src/cost_parameters.py` に追加 (旧版の 1e6/1e8 二重定義を一元化)

**revert された変更** (ユーザー指示で取消):
- CEPCI 更新 (544 → 798.8) — エージェントが WebSearch で勝手に引いた値のため revert
- Perry's Chemical Engineers' Handbook citation (dT_lm) — エージェント独断の citation のため削除
- Knapp et al. (1982) k_ij ≈ 0.0089 (PR EOS) — 同上
- Sattler et al. 2014 PtSn 劣化温度 (ODE 温度ガード) — 同上

### 2026-05-18 (午前) — コードベース全体監査
コードベース全体を 4 並列エージェントで監査し、隠れ仮置き・silent error を新規抽出:

- **§2 を新設** (BO ペナルティ・ソルバ設定)
  - ペナルティ 4 値: `solver_failure_okuyen=10000`, `spec_base_okuyen=50`, `spec_coef_okuyen=100`, `production_max_relative=0.02`
  - 2026-05-17 の `spec_base` 1000→50 変更、両側 spec 導入の経緯記載
  - ソルバ内側・外側・初期推定値の各パラメータ (15 項目)
  - BO 探索範囲のうち `!仮置き` マーカー未付与 4 件 (要マーカー追加)
- **§3 を新設** (数値ガード・トレランス) — 12 項目
  - `_T_COOLING_FLOOR` の +1K 余裕、Wang-Henke / brentq の各 xtol、ODE 温度ガード等
- **§4 を新設** (Silent fallback / コード品質問題) — 16 項目
  - 高: `warnings.simplefilter("ignore")` 5 箇所、CAPEX 閾値二重定義、ODE/LSODA 失敗時の silent
  - 中: Z=1.0 fallback、JT 失敗時 T2=T1 fallback、CSS warning-only、recovery check スキップ閾値
- **検証で除外した false positive** (参考):
  - `heat_integration.py:1080` の `1e9/3600/1000` 単位換算は正しい (= `1e6/3.6 = 277.78`)
  - `main.py:124` の生産量 target=1188 kmol/h は正しい (400,000 t/y × 1000 / 42.08 / 8000 h/y)

### 2026-05-13
- **Hasebe 式 (9)(10) 準拠の OPEX 計算実装** (`flowsheet/economics.py`)
  - 旧実装: TAC = CAPEX/n + Σ(utility + 触媒 + 吸着剤 + 原料費) ─ 労務費・保全費・1.23 倍率欠落
  - 新実装: TAC = CAPEX/n + 0.180·C_TM + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM) + 触媒・吸着剤交換 (式外、別建て)
  - opex dict に `[Hasebe-集計]` プレフィクス付きエントリで集計項を表示
  - HI/HEN 後の Economics 再構築でも `apply_hasebe_aggregation` で再計算
- **運転員年俸を新規 placeholder として追加** (`OPERATOR_ANNUAL_SALARY_JPY=600 万`)
- **Turton 直接引用を削除** (孫引きを避け、長谷部資料のみ引用)

### 2026-05-10
- Wang-Henke の重大 bug 連続発見・修正:
  - **B_bottoms 取り違え**: リボイラ方程式で `L[N_stages]` (= L_bot) を bottoms 産物流量 `B = F − D` と取り違え。Dist1 で 20× 過大、結果として C3H8 ↔ C4H10 が 102 kmol/h 入れ替わる成分マスバランス破綻
  - **V[N_feed] off-by-one**: フィード段の上向き気流量 `V[N_feed]` を `V_bot` と設定していた (本来 `V_top = V_bot + (1-q)F`)。修正前は Dist2 の C2H6 recovery が 12.8% に崩壊、修正後 99% へ
- **Newton 1 ステップ + Wegstein 加速** で 6.5× 高速化 (701s → 107s 全塔 rigorous)
- **always-on validation** を `RigorousResult` に組込: `mesh_residual_max`, `mesh_residual_mean`, `component_balance_max`
- **`tests/validate_wang_henke.py` を regression suite 化**: 全 48 検査 PASS
- **既知の数値限界 (修正困難)**: Dist2 partial_condenser stage 1 の T が brentq 偽根問題で 19K ズレ (実際は 1.5K 精度、後日 thermo パッケージ採用で根本解決)

### 2026-05-09 (午後)
- **コンテスト Ver.2.0 §4 仕様の citable 採用** (Phase A/B 完了)
  - §4-2 蒸留塔: G*=SF·K·√(ρ_v(ρ_l-ρ_v)) / 段間隔 0.6m / 段効率 80% / 塔頂2m+塔底4m
  - §4-4 熱交換器: U 値 9 ペア表 (cooler/distillation/HEN すべて `lookup_U` で参照)
- **化工便覧 改訂六版から citation 追記**
  - PUMP 効率 0.70 → 5·6·4 項【例題 5·8】
  - FURNACE 熱効率 0.85 (新規導入) → 18·4·3 項 表 18·11
  - COMPRESSOR η_poly 0.75 → p.333 (再確認)
- **コードベース全体監査で 10 件の隠れ仮置きを新規抽出**
  - ΔT_min, 凝縮器最低温度, dT_lm cooler, PSA 空塔速度/グリッド/最小吸着時間, 触媒再生時間, 反応器最大容積, スイングペナルティ, HE 材質係数
- **citation コメント追記** (universally agreed、文献必要なし):
  - PR EOS Ω_a/Ω_b → Peng-Robinson (1976) *Ind. Eng. Chem. Fundam.* 15(1), 59-64
  - `_T_BUBBLE_MIN/MAX` → 物性値ではなく数値ガードである旨を明示
- **VLE rigorous solver 実装** (Wang-Henke、`src/distillation_rigorous.py` 新規)
  - アルゴリズム出典: Seader, Henley & Roper (2010) Ch.10.4
  - K 値検証用: CalebBell/thermo (MIT) と src/eos.py の PR EOS で 0.02% 一致確認

### 2026-05-09
- PtSn 触媒単価 50,000 → **30,000** 円/kg
- PtSn 触媒寿命 3 → **4** 年
- PtSn 触媒充填密度 ρ_p 400 → **700** kg/m³ (`units/reactors/swing.py`)
- 原料・製品単価 (LPG/C3H6/H2) を新規追加 (`src/cost_parameters.py`)
- HHV 表 (`HHV_MJ_PER_KMOL`) を新規追加

### 2026-05-02
- 初版作成 (membrane CAPEX 仮置き値を明示)
