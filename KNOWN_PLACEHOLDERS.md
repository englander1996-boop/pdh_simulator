# 仮置きデータ一覧 (要出典・要更新)

このプロジェクトで「仮置き」になっている数値・パラメータの一覧。
コード中は `!仮置き` マーカーでグレップ可能 (`grep -rn "!仮置き" .`)。

実値に置き換える際は、本ファイルの該当行と、コード中のマーカーコメントの両方を更新すること。

最終更新: 2026-05-09

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
| **U_Wm2K (cooler の U 値)** | 200 W/(m²·K) | `units/utils/cooler.py` | phase 組み合わせで本来変わる (contest §4-4 表) |

## 🟢 優先度低 (構造的に変わる時に再検討)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **CEPCI 基準年** | 2001 (BASE) / 2016 (CURRENT) | `src/cost_parameters.py` `CEPCI_BASE`, `CEPCI_CURRENT` | 最新の CEPCI 値で更新 (現在 2024 年で ~800) |
| **USD/JPY レート** | 157.08 | 同上 `USD_TO_JPY` | 設計確定時点の為替で固定 |
| **減価償却年数** | 8 年 | 同上 `DEPRECIATION_YEARS` | 国税庁基準で固定済み (要確認) |
| **PUMP 効率** | 0.70 | `units/utils/pump.py` `eta_pump` | 汎用遠心ポンプ典型値 |
| **COMPRESSOR ポリトロピック効率** | 0.75 | `units/utils/compressor.py` `eta_poly` | 化工便覧 改訂六版 p.333 (0.7〜0.8 中央値) |
| **Q_reb = 1.05 × Q_cond (蒸留塔)** | 5% 損失仮定 | `src/distillation_core.py:554` | 文献根拠未確定 |

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

- **2026-05-09**:
  - PtSn 触媒単価 50,000 → **30,000** 円/kg
  - PtSn 触媒寿命 3 → **4** 年
  - PtSn 触媒充填密度 ρ_p 400 → **700** kg/m³ (`units/reactors/swing.py`)
  - 原料・製品単価 (LPG/C3H6/H2) を新規追加 (`src/cost_parameters.py`)
  - HHV 表 (`HHV_MJ_PER_KMOL`) を新規追加
- **2026-05-02**: 初版作成 (membrane CAPEX 仮置き値を明示)
