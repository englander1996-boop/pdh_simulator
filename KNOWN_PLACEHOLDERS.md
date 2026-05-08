# 仮置きデータ一覧 (要出典・要更新)

このプロジェクトで「仮置き」になっている数値・パラメータの一覧。
コード中は `!仮置き` マーカーでグレップ可能 (`grep -rn "!仮置き" .`)。

実値に置き換える際は、本ファイルの該当行と、コード中のマーカーコメントの両方を更新すること。

---

## 🔴 優先度高 (TAC に大きく効く / 制約判定に直結)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **電力単価** | 15 円/kWh | `src/cost_parameters.py` `ELECTRICITY_JPY_PER_KWH` | コンテスト課題 Ver.2.0 サイト仕様 / 産業統計 |
| **LP/MP/HP スチーム単価** | 1800/2200/2800 円/GJ | 同上 `*_STEAM_JPY_PER_GJ` | 同上、温度ごとに別単価 |
| **冷却水単価** | 60 円/GJ | 同上 `COOLING_WATER_JPY_PER_GJ` | 同上 |
| **冷凍冷媒単価 (7階層)** | 100〜9000 円/GJ | 同上 `*_REFRIG_*` | 冷凍冷媒製造プラントの実績 / Carnot 効率推定 |
| **燃料 (LNG) 単価** | 1500 円/GJ | 同上 `FUEL_JPY_PER_GJ` | 都市ガス料金 / コンテスト仕様 |
| **PtSn 触媒単価** | 50000 円/kg | 同上 `CATALYST_PTSN_JPY_PER_KG` | 触媒メーカー (Clariant 等) のデータシート |
| **PtSn 触媒寿命** | 3 年 | 同上 `CATALYST_PTSN_LIFE_YEARS` | 同上 + 運転条件依存 |
| **膜モジュール単価** | 50 USD/m² | 同上 `MEM_UNIT_PRICE_USD_PER_M2` | Hua et al. (2024) / ZIF-8 膜 TEA 論文 |
| **活性炭単価** | 5 USD/kg | 同上 `ACTIVATED_CARBON_PRICE_USD_PER_KG` | 試薬メーカーカタログ (工業グレード) |
| **活性炭寿命** | 4 年 | 同上 `ADSORBENT_LIFETIME_YEARS` | メーカー仕様書 / 運転実績 |

## 🟡 優先度中 (物理計算精度に効く)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **成分 Cp (定圧比熱)** | 範囲固定値 | `src/component_data.py` `CP_J_PER_MOL_K` | NIST Shomate 多項式 / DIPPR で温度依存関数化 |
| **成分蒸発潜熱** | 沸点での値固定 | 同上 `LATENT_HEAT_KJ_PER_KMOL` | NIST WebBook + Watson 式で温度補正 |
| **PSA K_Fa (C2H4, C2H6)** | 推算値 | `src/cost_parameters.py` `PSA_KFA` | 実機データで CH4 以外も直接測定して更新 |
| **膜モジュール A_per_module** | 500 m² | `units/separators/membrane/membrane_system.py` MemFixedParams | Evonik SEPURAN 等のデータシートから中空糸寸法計算 |

## 🟢 優先度低 (構造的に変わる時に再検討)

| 項目 | 現値 | ファイル | 想定出典 |
|---|---|---|---|
| **CEPCI 基準年** | 2001 (BASE) / 2016 (CURRENT) | `src/cost_parameters.py` `CEPCI_BASE`, `CEPCI_CURRENT` | 最新の CEPCI 値で更新 (現在 2024 年で ~800) |
| **USD/JPY レート** | 157.08 | 同上 `USD_TO_JPY` | 設計確定時点の為替で固定 |
| **減価償却年数** | 8 年 | 同上 `DEPRECIATION_YEARS` | 国税庁基準で固定済み (要確認) |

---

## 📋 グレップで全箇所一覧

```bash
grep -rn "!仮置き" --include="*.py" --include="*.md" .
```

代表的なヒット箇所:
- `src/cost_parameters.py` — ユーティリティ単価・触媒・吸着剤・膜単価
- `src/component_data.py` — Cp, 蒸発潜熱
- `units/separators/membrane/membrane_system.py` — 膜モジュール仕様
- `units/separators/psa/psa_system.py` — PSA 設計パラメータ

---

## 更新ルール

新しいデータが入ったら:
1. コード側の値を更新
2. 該当する `!仮置き` コメントを `# 出典: <文献名/データシート>` に変える
3. 本ファイルの該当行を消す
