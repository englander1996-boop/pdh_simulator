# SPEC: economics.py — 経済集計 (CAPEX / OPEX / TAC / Profit)

**ファイルパス**: `flowsheet/economics.py`
**単価出典の集約先**: `src/cost_parameters.py` (最終単価チューニングはそちら 1 ファイルで完結)

---

## 1. 目的

`run_one_pass` の戻り値 (1 パス分の全ユニット結果) から、装置別 CAPEX・OPEX・
Revenue を抽出して集計し、TAC・Profit・製品単価を計算する。CAPEX は各ユニットが
Bare Module 法で算出済みの値を受け取り、OPEX は Hasebe 経験式 (10) に従って集計する。

---

## 2. モデル理論 (式)

### 2-1. CAPEX (Bare Module 法)

各装置の Bare Module Cost は **個々のユニットモジュール側** (`units/.../`) で算出され、
`equipment.CAPEX` 系属性として返る。`economics.py` はそれを収集するだけ
(`collect_capex_opex`)。集計対象 (キー / 億円):

```
Pump1, Dist1, Reactor, Cooler, Comp2a, Intercool, Comp2b, Desuper, Dist2,
PSA容器, PSA活性炭, MemPrecool, Mem気化器, Mem F圧縮機, Mem P圧縮機,
Mem冷却器, Mem段間冷却器, Mem膜本体, Dist3
```

- ペナルティ装置 (CAPEX ≥ `PENALTY_CAPEX_THRESHOLD_OKUYEN = 1e8 億円` sentinel) は
  `total_capex` 合計から除外。
- `total_capex = Σ (v for v in capex.values() if v < 1e8)`

### 2-2. OPEX (Hasebe 式 (10))

出典: 長谷部『プロセス設計 R08-3.pdf』§3.4 式 (10)。製造コストの経験式:

```
製造コスト = 0.180·C_TM + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM) + 減価償却費
```

- `C_TM` = 全プラント Bare Module コスト (= `total_capex`)
- `C_OL` = 労務費 (Hasebe 式 (9) で算出、後述)
- `C_UT` = 用役費 (電力・蒸気・冷却水・冷媒・燃料・予熱)
- `C_WT` = 廃棄物処理費 (本プロセスは `HASEBE_C_WT_OKUYEN_PER_YEAR = 0.0`)
- `C_RM` = 原料費 (Fresh LPG)

本実装での扱い (`Economics` / `apply_hasebe_aggregation`):

- 減価償却費は OPEX に入れず、TAC 側で `CAPEX / DEPRECIATION_YEARS` として別途加算。
- `opex` dict は**装置別の生 (1.00 倍) エントリ** + **Hasebe 集計項**
  (`[Hasebe-集計] ` プレフィクス) を持つ。用役費・原料費は生エントリが既に 1.00×
  計上済みなので、Hasebe 集計項には **delta 分のみ** (= 0.23×) を加える。
- 不変条件: `sum(opex.values()) == total_opex`。

集計項 (`[Hasebe-集計] `):

```
0.180·C_TM   (保全 + 諸経費 + 税保険)
2.73·C_OL    (労務費)
0.23·C_UT    (用役費 上乗せ分 = 1.23 − 1.00)
0.23·C_RM    (原料費 上乗せ分)
1.23·C_WT    (廃棄物処理、本プロセスは 0)
```

従って実効 OPEX:
```
total_opex = 0.180·C_TM + 2.73·C_OL + 1.23·C_UT + 1.23·C_WT + 1.23·C_RM
             + 触媒・吸着剤・膜交換費 (Hasebe 式枠外、別建て、1.23× なし)
```

#### 労務費 C_OL — Hasebe 式 (9)

```
1 班の人数 = ceil( sqrt(6.29 + 0.21·N_eq) )       (式 (9))
総運転員  = 1 班の人数 × 4                          (4 直 3 交替)
C_OL      = 総運転員 × OPERATOR_ANNUAL_SALARY_JPY / 1e8   [億円/年]
```

`N_eq` (主要機器数) は `_count_main_equipment`: capex dict の物理機器キー数
(`PSA活性炭` と `[Hasebe...` を除く) に、運転監視上独立した加熱炉として
**Reactor予熱炉 +1** を加算。

### 2-3. ユーティリティ換算式

```
電力:  W [kW] × ELECTRICITY_JPY_PER_KWH × OPERATING_HOURS_PER_YEAR / 1e8     [億円/年]
熱:    Q [kW] × 3.6e-3 [GJ/(kW·h)] × OPERATING_HOURS_PER_YEAR × 円/GJ / 1e8  [億円/年]
流量:  F [kmol/h] × MW [kg/kmol] × OPERATING_HOURS_PER_YEAR × 円/kg / 1e8    [億円/年]
```

- 蒸留塔のリボイラ/コンデンサ単価は `equipment` が `utility_selector` で動的選択した
  名前・単価 (`reb_utility_*`, `cond_utility_*`) を使う (例: Dist2 塔頂が氷点下なら
  エチレン冷媒)。Cooler/Intercool/Desuper/MemPrecool も同様に `equipment.utility_*`。
- 反応器予熱燃料: `Q_fuel = Q_preheat / FURNACE_EFFICIENCY` (LHV ベース) を
  `FUEL_JPY_PER_GJ` で評価。
- 触媒交換: `Catalyst_Weight_Total × CATALYST_JPY_PER_KG / CATALYST_LIFE_YEARS / 1e8`
- PSA 活性炭交換: `equipment.OPEX_adsorbent_okuyen_per_year`
- 膜交換: `CAPEX_mem / MEM_LIFETIME_YEARS` (!仮置き、env `PDH_MEM_LIFETIME_YEARS` で上書き可)

### 2-4. Revenue (全て正値で計上)

```
C3H6 製品売上         = annual(F_C3H6_prod, MW_B, C3H6_PRODUCT_JPY_PER_KG)
H2 製品売上           = annual(F_H2_prod,   MW_C, H2_PRODUCT_JPY_PER_KG)
PSA オフガス燃料CR    = HHV(offgas) [GJ/h] × hours × FUEL_JPY_PER_GJ / 1e8
Dist1 塔底 燃料CR     = HHV(r1.bottom) [GJ/h] × hours × FUEL_JPY_PER_GJ / 1e8
```

オフガス・Dist1 塔底 (C4H10 主) は反応器プリヒーター燃料として利用 → 浮いた燃料費を
クレジットとして計上。HHV は `HHV_MJ_PER_KMOL` を成分別に積算。

### 2-5. TAC / Profit / 製品単価

```
TAC            = total_capex / DEPRECIATION_YEARS + total_opex     [億円/年]
profit         = total_revenue − TAC                              [億円/年, 正=黒字]
unit_jpy_per_t = TAC × 1e8 / (annual_kg_C3H6 / 1000)              [円/ton]
```

---

## 3. 入出力

- 入力: `calculate_economics(one_pass: dict, mw_C3H6_kg_per_kmol: float)`
- 出力: `Economics` dataclass
  - `capex` (dict, 億円), `opex` (dict, 億円/年, 生 + Hasebe 集計, 全て正),
    `revenue` (dict, 億円/年, 全て正)
  - `total_capex`, `total_opex`, `total_revenue`, `TAC`, `profit`,
    `annual_kg_C3H6`, `unit_jpy_per_t`

補助関数: `collect_capex_opex`, `apply_hasebe_aggregation`,
`_classify_opex_term` (UT/RM/CATALYST_OUT/HASEBE_AGGR の区分判定),
`_count_main_equipment`, `_compute_labor_cost_okuyen`。

---

## 4. 主要パラメータ (`src/cost_parameters.py`)

| 定数 | 値 | 区分 |
|---|---|---|
| `HASEBE_COEFF_C_TM` | 0.180 | 保全+諸経費+税保険 |
| `HASEBE_COEFF_C_OL` | 2.73 | 労務費倍率 |
| `HASEBE_COEFF_C_UT_WT_RM` | 1.23 | 用役費・廃棄物・原料費の間接費込み倍率 |
| `HASEBE_NOL_COEFF` | 0.21 | 式 (9) の N_eq 係数 |
| `HASEBE_SHIFT_MULTIPLIER` | 4 | 4 直 3 交替倍率 |
| `HASEBE_C_WT_OKUYEN_PER_YEAR` | 0.0 | 廃棄物処理費 (本プロセスは 0) |
| `OPERATOR_ANNUAL_SALARY_JPY` | 6,000,000 | 運転員年俸 (!仮置き) |
| `DEPRECIATION_YEARS` | 8 | 減価償却年数 |
| `OPERATING_HOURS_PER_YEAR` | 8000 | 年間稼働時間 (contest 規定) |
| `MEM_LIFETIME_YEARS` | 3.0 | 膜寿命 (!仮置き) |
| `CATALYST_LIFE_YEARS` | 2.5 | 触媒寿命 |
| `PENALTY_CAPEX_THRESHOLD_OKUYEN` | 1e8 | penalty sentinel 判定閾値 |

---

## 5. 出典

- OPEX 経験式 (10) / 労務費式 (9): 長谷部『プロセス設計 R08-3.pdf』§3.3, §3.4。
- 係数 0.21 / ×4 は一次出典 (長谷部資料) 準拠 (Turton 4th ed. の 0.23 / ×4.5 とは
  あえて異なる、`cost_parameters.py` に注記)。
- 稼働時間 8000 h/年: PDH 学生コンテスト要綱 Ver.2.0。
- TAC に原料費を含めるのは化工標準 (Sinnott §6.5, Turton §8.2)。
- 個々の装置 CAPEX (Bare Module) の算式は各ユニット SPEC を参照。

---

## 6. 既知の限界・仮置き

- 運転員年俸・膜寿命・各種単価 (LPG/製品/冷却水/冷媒/スチーム/HHV/CEPCI) の一部は
  `!仮置き`。詳細は `KNOWN_PLACEHOLDERS.md` §1.1〜1.3。
- `N_eq` の数え方 (蒸留塔は付帯 HE 込みで 1 単位、予熱炉 +1 等) は Hasebe/Turton
  慣行に基づく判断であり、厳密な機器カウント規約ではない。
- 膜交換費は CAPEX_mem を耐用年数で割った年均等費 (簡易)。
- `C_WT = 0` 固定 (気相反応で水処理対象廃棄物が実質発生しない前提)。
