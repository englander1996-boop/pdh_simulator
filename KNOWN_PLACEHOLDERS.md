# 仮置きデータ一覧 (要出典・要更新)

このプロジェクトで「仮置き」になっている数値・パラメータの一覧。
コード中は `!仮置き` マーカーでグレップ可能 (`grep -rn "!仮置き" .`)。

実値に置き換える際は、本ファイルの該当行と、コード中のマーカーコメントの両方を更新すること。

最終更新: 2026-05-09 (Phase A/B コンテスト仕様採用、化工便覧 citation 追記後)

---

## 🔴 優先度高 (TAC に大きく効く / 制約判定に直結)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **電力単価** | 15 円/kWh | `src/cost_parameters.py` `ELECTRICITY_JPY_PER_KWH` | コンテスト課題 Ver.2.0 サイト仕様 / 産業統計 |
| **LP/MP/HP スチーム単価** | 1800/2200/2800 円/GJ | 同上 `*_STEAM_JPY_PER_GJ` | 同上、温度ごとに別単価 |
| **冷却水単価** | 60 円/GJ | 同上 `COOLING_WATER_JPY_PER_GJ` | 同上 |
| **冷凍冷媒単価 (7階層)** | 100〜9000 円/GJ | 同上 `*_REFRIG_*` | 冷凍冷媒製造プラントの実績 / Carnot 効率推定 |
| **燃料 (LNG) 単価** | 1500 円/GJ | 同上 `FUEL_JPY_PER_GJ` | 都市ガス料金 / コンテスト仕様 |
| **PtSn 触媒単価** | **30,000 円/kg** (2026-05-09 50,000→30,000) | 同上 `CATALYST_PTSN_JPY_PER_KG` | 触媒メーカー (Clariant 等) のデータシート、Pt 含有率依存 |
| **PtSn 触媒寿命** | **4 年** (2026-05-09 3→4) | 同上 `CATALYST_PTSN_LIFE_YEARS` | 工業実績 (Catofin/Oleflex 共に 3〜5 年)、コーキング・再生回数依存 |
| **膜モジュール単価** | 50 USD/m² | 同上 `MEM_UNIT_PRICE_USD_PER_M2` | Hua et al. (2024) / ZIF-8 膜 TEA 論文 |
| **活性炭単価** | 5 USD/kg | 同上 `ACTIVATED_CARBON_PRICE_USD_PER_KG` | 試薬メーカーカタログ (工業グレード) |
| **活性炭寿命** | 4 年 | 同上 `ADSORBENT_LIFETIME_YEARS` | メーカー仕様書 / 運転実績 |

## 🔴 原料・製品単価 / 副産物クレジット (NEW: 2026-05-09 追加)

経済評価の TAC・Revenue・Profit の支配項。`flowsheet/economics.py` で使用。

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **LPG 原料単価 (C3H8 + C4H10 共通)** | 50 円/kg | `src/cost_parameters.py` `LPG_FEED_JPY_PER_KG` | 2024 LPG CIF ~600 USD/t × 157 ≒ 94 円/kg 程度。コンテスト規定値で確定 |
| **C3H6 製品売価** | 120 円/kg | 同上 `C3H6_PRODUCT_JPY_PER_KG` | 国内市況 100〜140 円/kg レンジの中央値 |
| **H2 製品売価** | 350 円/kg | 同上 `H2_PRODUCT_JPY_PER_KG` | 副生 H2 OEM 30〜50 円/Nm³ ≒ 330〜560 円/kg の下限 |
| **HHV (高位発熱量)** | 表 (`HHV_MJ_PER_KMOL`) | 同上 | NIST WebBook (25°C, H2O 液基準)。オフガス燃料クレジット計算用 |

> **注意**: コンテスト Ver.2.0 サイト仕様で原料・製品単価が指定されていれば、
> その値を最優先で適用すること。LPG・C3H6 は地域市況依存性が大きい。

## 🟡 優先度中 (物理計算精度に効く)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **成分 Cp (定圧比熱)** | 範囲固定値 | `src/component_data.py` `CP_J_PER_MOL_K` | NIST Shomate 多項式 / DIPPR で温度依存関数化 |
| **成分蒸発潜熱** | 沸点での値固定 | 同上 `LATENT_HEAT_KJ_PER_KMOL` | NIST WebBook + Watson 式で温度補正 |
| **PSA Langmuir パラメータ** | 推算値 | `src/cost_parameters.py` `PSA_LANGMUIR_PARAMS` | CH4/C2H4/C2H6 単成分 25°C 等温線測定 |
| **PSA K_Fa (物質移動係数)** | 推算値 | 同上 `PSA_KFA` | Carberry 数 + Knudsen 拡散の理論計算、または実機データ |
| **PSA 活性炭嵩密度 ρ_b** | 600 kg/m³ | `units/separators/psa/psa_system.py` `PSAFixedParams.rho_b` | 工業用活性炭 (典型 400〜700) のデータシート |
| **PSA CSS 近似** | True (保守過大推算) | 同上 `use_css_approximation` | 厳密 CSS の数値検証で評価 |
| **PSA 脱着安全係数** | 1.2 | 同上 `desorption_time_safety_factor` | KFa 確定後に再評価 |
| **膜モジュール A_per_module** | 500 m² | `units/separators/membrane/membrane_system.py` `MemFixedParams.A_per_module` | Evonik SEPURAN 等のデータシートから中空糸寸法計算 |
| **膜パラメータ (Q_A=40 GPU, α=90)** | 文献値 | 同上 `MemFixedParams` | Hua et al. (2024) 実測値、室温・大気圧 |
| **PR EOS 二成分相互作用係数 k_ij** | 0 | `src/eos.py` | 文献で 0.01 程度と小さく無視可。要再評価 |
| **PtSn 触媒充填密度 ρ_p** | **700 kg/m³** (2026-05-09 400→700) | `units/reactors/swing.py` `FixedParams.rho_p` | 実触媒 PtSn/Al2O3 ペレットの粒子密度寄り。空隙の二重控除になっている可能性は別途要確認 |

## 🟡 優先度中 — 2026-05-09 監査で新規抽出 (隠れ仮置き)

コードベース全体監査 (2026-05-09) で発見、本ファイル未収録だった項目:

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **ΔT_min (HI 最小接近温度)** | 10.0 K | `src/utility_selector.py` L60、`exp/exp1.py` `HI_DT_MIN_K` | 教科書慣行値 (Smith/Sinnott/Linnhoff)。textbook citation 追記 |
| **凝縮器最低温度** | 313.15 K (40°C) | `src/distillation_core.py` L83 `_T_COND_MIN` | 冷却水 supply 30°C + ΔT_min=10K の合成。citation 追記必要 |
| **dT_lm デフォルト (cooler)** | 30 K | `units/utils/cooler.py` L69 | LMTD 代替値。Perry's HE 設計章の典型値などで citation 追記 |
| **PSA 空塔速度上限** | 1.0 m/s | `units/separators/psa/psa_system.py` L164 | 化工便覧 §13-31 引用ありだが「除湿用」値。PDH オフガス分離での適用妥当性要確認 |
| **PSA グリッド分割数** | 20 | 同上 L150 | 数値拡散 vs 計算速度のトレードオフ、感度解析未実施 |
| **PSA 最小吸着時間** | 60 s | 同上 L165 | CSS 近似の物理下限ガード、出典未明示 |
| **触媒再生時間** | 30 min | `units/reactors/swing.py` L105 `t_regen` | スイング設計の典型値、工業データ citation 欠落 |
| **反応器最大触媒容積/基** | 200 m³ | 同上 L106 `V_cat_max` | 設計判断 (大型固定床上限) の citation 欠落 |
| **スイング操作ペナルティ係数** | 1.2 | `src/cost_parameters.py` L69 | スイング切替弁・マニホールド追加 20% 増の見積もり根拠 |
| **HE 材質係数 FM/FM_HE/FP_HE** | 1.0 | `src/cost_parameters.py` L55, L120, L125 | 「炭素鋼・圧損補正なし」前提。本フロー (高圧 C3/H2 系) での材質選定根拠欠落 |

## 🟢 優先度低 (構造的に変わる時に再検討)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **CEPCI 基準年** | 2001 (BASE) / 2016 (CURRENT) | `src/cost_parameters.py` `CEPCI_BASE`, `CEPCI_CURRENT` | 最新の CEPCI 値で更新 (現在 2024 年で ~800) |
| **USD/JPY レート** | 157.08 | 同上 `USD_TO_JPY` | 設計確定時点の為替で固定 |
| **減価償却年数** | 8 年 | 同上 `DEPRECIATION_YEARS` | 国税庁基準で固定済み (要確認) |
| **Q_reb = 1.05 × Q_cond (蒸留塔)** | 5% 損失仮定 | `src/distillation_core.py:554` | 文献根拠未確定 |

> **出典確定済み (本ファイルから卒業した項目)**:
> - PUMP 効率 0.70 — 化工便覧 改訂六版 5·6·4 項【例題 5·8】(2026-05-09)
> - COMPRESSOR ポリトロピック効率 0.75 — 化工便覧 改訂六版 p.333 (2026-05-09 確認)
> - FURNACE 熱効率 0.85 — 化工便覧 改訂六版 18·4·3 項 表 18·11 (2026-05-09)
> - HE U 値 — 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 表 (2026-05-09)
> - 蒸留塔 G* / 段間隔 0.6m / 段効率 80% / 塔頂2m+塔底4m — contest Ver.2.0 §4-2 (2026-05-09)

---

## 📋 グレップで全箇所一覧

```bash
grep -rn "!仮置き" --include="*.py" --include="*.md" .
```

代表的なヒット箇所:
- `src/cost_parameters.py` — ユーティリティ単価・触媒・吸着剤・膜単価・原料/製品単価・HHV
- `src/component_data.py` — Cp, 蒸発潜熱
- `units/separators/membrane/membrane_system.py` — 膜モジュール仕様
- `units/separators/psa/psa_system.py` — PSA 設計パラメータ
- `units/reactors/swing.py` — 触媒充填密度

---

## 更新ルール

新しいデータが入ったら:
1. コード側の値を更新
2. 該当する `!仮置き` コメントを `# 出典: <文献名/データシート>` に変える
3. 本ファイルの該当行を消す or 出典欄を埋める
4. 影響を受ける SPEC.md (該当ユニットの `units/*/SPEC_*.md`) も更新

---

## 履歴

- **2026-05-09 (午後)**:
  - **コンテスト Ver.2.0 §4 仕様の citable 採用** (Phase A/B 完了)
    - §4-2 蒸留塔: G*=SF·K·√(ρ_v(ρ_l-ρ_v)) / 段間隔 0.6m / 段効率 80% / 塔頂2m+塔底4m
    - §4-4 熱交換器: U 値 9 ペア表 (cooler/distillation/HEN すべて lookup_U で参照)
  - **化工便覧 改訂六版から citation 追記**
    - PUMP 効率 0.70 → 5·6·4 項【例題 5·8】
    - FURNACE 熱効率 0.85 (新規導入) → 18·4·3 項 表 18·11
    - COMPRESSOR η_poly 0.75 → p.333 (再確認)
  - **コードベース全体監査で 10 件の隠れ仮置きを新規抽出** (本ファイル「優先度中 — 監査で新規抽出」表)
    - ΔT_min, 凝縮器最低温度, dT_lm cooler, PSA 空塔速度/グリッド/最小吸着時間, 触媒再生時間, 反応器最大容積, スイングペナルティ, HE 材質係数
  - **citation コメント追記** (citation 文献必要なし、universally agreed):
    - PR EOS Ω_a/Ω_b → Peng-Robinson (1976) Ind. Eng. Chem. Fundam. 15(1), 59-64
    - `_T_BUBBLE_MIN/MAX` → 物性値ではなく数値ガードである旨を明示
- **2026-05-09**:
  - PtSn 触媒単価 50,000 → **30,000** 円/kg
  - PtSn 触媒寿命 3 → **4** 年
  - PtSn 触媒充填密度 ρ_p 400 → **700** kg/m³ (`units/reactors/swing.py`)
  - 原料・製品単価 (LPG/C3H6/H2) を新規追加 (`src/cost_parameters.py`)
  - HHV 表 (`HHV_MJ_PER_KMOL`) を新規追加
- **2026-05-02**: 初版作成 (membrane CAPEX 仮置き値を明示)
