# SPEC: psa_system.py — PSA (Pressure Swing Adsorption) システム

**ファイルパス**: `units/separators/psa/psa_system.py`
**最終更新**: 2026-05-09

---

## 1. 目的

Dist2 塔頂ガス (H2 + CH4 + C2 軽質 + C3 微量) から **H2 を高純度で回収**し、
吸着成分 (CH4/C2H4/C2H6) と C3 成分をオフガスとして燃料系へ送る。

```
[Dist2 塔頂] → Preheater → 吸着塔 (PDE + LDF) → H2 製品 (≥99.9 mol%)
                                                  │
                                              脱着 (簡易指数減衰) → オフガス (燃料へ)
```

## 2. 構成

| 段 | 役割 |
|---|---|
| 1. Preheater | フィードを T_abs = 25°C に温度調整 (顕熱) |
| 2. 吸着塔 (Adsorption) | 1D 上流差分 PDE + LDF。CH4 出口/入口比が `breakthrough_ratio` (0.1%) に到達した時点で吸着終了 |
| 3. 脱着 (Desorption) | 大気圧パージ。`q(t) = q₀·exp(−KFa·t)`、`Σq/Σq₀ = desorption_target` で停止 |

## 3. 成分マッピング

| キー | 化学式 | PSA での扱い |
|---|---|---|
| A | C3H8 | **全量オフガス** (吸着層前段で完全捕捉、PDE から除外) |
| B | C3H6 | **全量オフガス** (同上) |
| C | H2 | **非吸着、全量プロダクト** (吸着量 < 0.01 mol/kg @25°C で無視) |
| D | C2H4 | 吸着 (PDE 対象) |
| E | CH4 | 吸着 (PDE 対象、**破過基準成分**) |
| F | C2H6 | 吸着 (PDE 対象) |

## 4. インターフェース

### 入力

`PSADesignVars` (BO 探索):
- `D_col` 塔径 [m]
- `L_bed` 吸着層高さ [m]
- `desorption_target` 脱着完了基準 [-] (例 0.35 → q が初期値の 35% まで)

`PSAFeedStream`: `F_in [kmol/h]`, `T_in [K]`, `P_in [Pa]`

`PSAFixedParams`:
- `T_abs = 298.15 K`、`P_des = 101325 Pa`、`rho_b = 600 kg/m³` (★仮置き)、
  `eps = 0.4`、`breakthrough_ratio = 0.001`、`t_ads_max = 7200 s`、
  `use_css_approximation = True` (★仮置き)、`desorption_time_safety_factor = 1.2` (★仮置き)

### 出力

`PSASimulationResult`:
- `product: dict` H2 リッチプロダクト [kmol/h] (キー A〜F)
- `offgas: dict` オフガス [kmol/h] (キー A〜F)
- `equipment: PSAEquipmentData` (下記)
- `H2_recovery: float`, `CH4_capture: float`

`PSAEquipmentData`:
- `N_abs_parallel = 1` (並列なし) / `N_cycle_sets = ceil(t_des/t_abs)+1` / `N_total_columns`
- `t_abs_sec`, `t_des_sec`, `u_0`
- `W_adsorbent_kg`, `Q_preheat_kW`
- `CAPEX_vessels`, `CAPEX_adsorbent`, `CAPEX_total`
- `H2_loss_blowdown_kmolh`, `H2_loss_purge_kmolh`
- `OPEX_adsorbent_okuyen_per_year` (吸着剤交換、`ADSORBENT_LIFETIME_YEARS` 依存)

## 5. 物理モデル

### 吸着 PDE (1D 上流差分)
```
∂C_i/∂t = −(u_0/ε)·∂C_i/∂z − (ρ_b/ε)·∂q_i/∂t       (ガス相)
∂q_i/∂t = KFa_i · (q_i* − q_i)                       (LDF、固相)
q_i*    = q_si · a_i · C_i / (1 + Σ a_j · C_j)       (Markham-Benton 多成分 Langmuir)
```

### 脱着 (簡易)
```
q_i(t) = q_i0 · exp(−KFa_i · t)
停止条件: Σ q_i(t) / Σ q_i(0) = desorption_target
```

### CSS (Cyclic Steady State) 近似
`use_css_approximation=True` のとき初期固相負荷を「フィード濃度の Langmuir 平衡 ×
desorption_target」で初期化 → 床全体が飽和の保守的設定 (実 CSS より過大、t_abs が
清浄床より短く出る → 必要塔数は多めに出る = コスト保守的)。

### 物理ガード
| 制約 | 違反時の挙動 |
|---|---|
| `u_0 ≤ U0_MAX = 1.0 m/s` | 超過は CAPEX = 1e9 ペナルティ (流動化・チャネリング防止) |
| `t_abs ≥ T_ABS_MIN = 60 s` | 未満はペナルティ (scale 発散防止 + 物理的に短サイクル無意味) |
| `t_ads ≤ 7200 s` | solve_ivp 打ち切り上限。超過はペナルティ |

## 6. CAPEX / OPEX

- **CAPEX_vessels**: `calc_reactor_capex_okuyen(V_col, P_in, D_col, N_total_columns)` (縦型容器)
- **CAPEX_adsorbent**: `W × ACTIVATED_CARBON_PRICE_USD_PER_KG × CEPCI 補正` (★仮置き 5 USD/kg)
- **OPEX 吸着剤交換**: `W × 単価 / ADSORBENT_LIFETIME_YEARS` (★仮置き 4 年)
- **OPEX 予熱**: `Q_preheat_kW` を `economics.py` 側でスチーム単価と掛けて計上 (※現状未経路、要追加)

## 7. 仮置き / 仮定

### 文献根拠あり
- C3 成分の完全捕捉 (Schell et al. 2012)
- H2 非吸着 (Poirier & Darriet 2001)
- Langmuir Markham-Benton 形 (Yang 1987)
- LDF (Glueckauf 1955; Ruthven 1984)

### 設計判断 (根拠あり)
- `T_abs = 25°C` (冷却水到達下限 + Langmuir 仮置き条件と整合)
- `P_des = 1 atm` (真空ポンプ不要、コスト最小化優先)
- 等速近似 (初期スクリーニング用)。最良点では吸着対象 (CH4/C2H4/C2H6) 合計が数十 mol% に達し
  厳密には無視できないが、等速近似の誤差は破過を早める方向 = 塔数を多めに見積もる保守側に出る。
  実機サイズ確定時に可変速度・非等温・温度依存 Langmuir へ拡張 (レポート「水素分離」章と整合)
- `breakthrough_ratio = 0.001` (CH4 捕捉率 > 99.9% 保証)
- `N_z = 20` (10/20 で t_abs 差 < 5% を確認済み)
- `_U0_MAX = 1.0 m/s` (化工便覧第 13 章 図 13・31)

### ★ 仮置き — 文献未確定
| 項目 | 値 | 確認方法 |
|---|---|---|
| Langmuir パラメータ (q_s, a) | 推算値 | 25°C 等温線測定 (CH4/C2H4/C2H6 単成分) |
| KFa (物質移動係数) | 推算値 | Carberry 数 + Knudsen 拡散の理論計算、または実機データ |
| `rho_b` 活性炭嵩密度 | 600 kg/m³ | 工業用活性炭 (典型 400〜700) のデータシート |
| `use_css_approximation` | True | 厳密 CSS の数値検証 |
| `desorption_time_safety_factor` | 1.2 | KFa 確定後に再評価 |
| 活性炭単価 | 5 USD/kg | メーカー (Calgon Carbon, Cabot 等) カタログ |
| 活性炭寿命 | 4 年 | 運転実績、再生回数依存 |

## 8. 依存

- `src/config.py`: `R`, `THERMO_DATA`
- `src/eos.py`: `z_factor`
- `src/thermo.py`: `PDHThermo` (エンタルピー計算)
- `src/cost_calculator.py`: `calc_reactor_capex_okuyen`
- `src/cost_parameters.py`: `PSA_LANGMUIR_PARAMS`, `PSA_KFA`, `ACTIVATED_CARBON_PRICE_USD_PER_KG`, `ADSORBENT_LIFETIME_YEARS`, etc.

## 9. 呼び出し箇所

- `flowsheet/run_one_pass.py`: Dist2 塔頂を直接 `simulate_psa_system` へ
- `exp/exp1.py`
- `monitor/psa_random_test.ipynb`: 単独感度試験

## 10. 関連ドキュメント

- `units/separators/psa/error_handling_fixes_20260505.md` (旧バグ修正履歴)
- `KNOWN_PLACEHOLDERS.md` (★仮置き項目の全体一覧。Langmuir/KFa/rho_b/活性炭単価・寿命/
  use_css_approximation/desorption_time_safety_factor の確認方法)
