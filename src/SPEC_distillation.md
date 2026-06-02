# SPEC: 蒸留計算コア — FUG / Wang-Henke 厳密 / GPR サロゲートの 3 経路

**ファイルパス**: `src/distillation_core.py`, `src/distillation_rigorous.py`, `src/distillation_sm.py`
**最終更新**: 2026-06-02 (新規作成。実コードからの式・定数・入出力の抽出)

---

## 目次

1. [目的](#1-目的)
2. [3 経路の使い分け](#2-3-経路の使い分け)
3. [共通データクラスとディスパッチ](#3-共通データクラスとディスパッチ)
4. [経路1: FUG ショートカット (distillation_core.py)](#4-経路1-fug-ショートカット-distillation_corepy)
5. [経路2: Wang-Henke 厳密 (distillation_rigorous.py)](#5-経路2-wang-henke-厳密-distillation_rigorouspy)
6. [経路3: GPR サロゲート (distillation_sm.py)](#6-経路3-gpr-サロゲート-distillation_smpy)
7. [塔径・塔高・CAPEX・ユーティリティ](#7-塔径塔高capexユーティリティ)
8. [入力・出力](#8-入力出力)
9. [主要パラメータ](#9-主要パラメータ)
10. [出典](#10-出典)
11. [既知の限界・仮置き](#11-既知の限界仮置き)

---

## 1. 目的

多成分蒸留塔の物質収支・温度・熱量・装置寸法・CAPEX を計算する。同一の入出力契約
（`DistDesignVars` / `ProcessStream` → `DistResult`）のもとで **求解アルゴリズムを 3 通り**
切り替えられる。`ColumnTunables.solver_method` / `DistDesignVars.solver_method` で選択する:

| solver_method | 経路 | 速度 | 用途 |
|---|---|---|---|
| `'fug'` | FUG ショートカット | 高速 | BO 探索のデフォルト |
| `'rigorous'` | Wang-Henke 厳密 (MESH tray-by-tray) | 中速 | top-k 設計の物理検証 |
| `'sm'` | GPR サロゲート | 最速 | HYSYS 置換（Dist1/Dist3） |
| `'hysys'` | (HYSYS COM、本ディレクトリ外) | 低速 | 真値（`units/vle/hysys`） |

成分マッピングは `src/` 共通（A=C3H8, B=C3H6, C=H2, D=C2H4, E=CH4, F=C2H6, Z=C4H10）。

---

## 2. 3 経路の使い分け

- **FUG**: Fenske-Underwood-Gilliland のショートカット式で N_min / R_min / 実段数・分配を
  代数的に解く。最も速く BO の内側ループで使う。ただし楽観的（narrow-α や partial condenser
  で feasible 判定が甘い）ため、**proxy 罰則**（`_compute_proxy_penalty`）を BO objective に加算して
  rigorous との乖離を補正する。
- **Wang-Henke 厳密**: MESH 方程式を段ごとに解く反復法。FUG では捉えられない段プロファイル・
  非鍵成分分配・partial condenser の non-condensable パススルーを物理的に扱う。FUG が出した
  D_total・T プロファイルを初期値として受け取る。narrow-α（Dist3, α≈1.07）では damping で振動抑制。
- **GPR サロゲート**: HYSYS の COM オーバーヘッド（毎パス HSC を開き直す数百 ms〜秒）を回避するため、
  Dist1/Dist3 を学習済みガウス過程回帰（GPR）で置換。入力→出力の決定的写像なので経路依存がなく、
  HYSYS 塔が Dist2 のみになり swap が原理的に消える。

---

## 3. 共通データクラスとディスパッチ

### `ColumnTunables`（BO/exp で振る P/N/R）

| フィールド | 型 | 意味 |
|---|---|---|
| `P_col` | float | 塔操作圧力 [Pa] |
| `N_stages` | int | 理論段数 |
| `N_feed` | int | フィード段位置（rigorous/sm では Kirkbride 推奨を自動採用、FUG では post-hoc 出力のみ） |
| `reflux_ratio` | float | 還流比 R = L/D |
| `solver_method` | str | 'fug' \| 'rigorous' \| 'sm' \| 'hysys' |
| `recovery_LK_top` / `recovery_HK_bot` | Optional[float] | 回収率（None でラッパ既定 0.99） |

### `DistDesignVars`（設計変数 + 分離仕様）

BO 変数（P_col, N_stages, N_feed, reflux_ratio）に加え、塔別ラッパで固定する仕様
`LK`（軽キー）, `HK`（重キー）, `recovery_LK_top`(0.99), `recovery_HK_bot`(0.99),
`K_method`('pr'|'cc'), `q`（フィード状態, 1=飽和液）, `partial_condenser`(bool),
`solver_method`, `D_override`（HYSYS の Distillate Rate spec 相当、直接 D 指定）。

### `DistFixedParams`（寸法・流体力学、コンテスト Ver.2.0 §4-2）

`tray_spacing_m=0.6`, `top_section_m=2.0`, `bot_section_m=4.0`, `tray_efficiency=0.8`,
`SF=0.8`, `K_factor=0.05`。

### ディスパッチ（`simulate_distillation_column`）

```
simulate_distillation_column(design, feed, fixed)
  ├─ solver_method == 'rigorous' → _simulate_rigorous (= wang_henke_solve)
  │     失敗時: PDH_RIGOROUS_STRICT=1(既定) なら penalty_result/RuntimeError、=0 なら warn+FUG
  ├─ solver_method == 'sm'       → _simulate_sm (未対応なら NotImplementedError → FUG)
  └─ それ以外('fug')             → FUG 本体
```

> rigorous は **fail-fast がデフォルト**。silent FUG fallback は「表示は rigorous なのに実体は FUG」
> という信頼性問題を生むため避ける。BO で退路が必要なら環境変数 `PDH_RIGOROUS_STRICT=0`。

---

## 4. 経路1: FUG ショートカット (distillation_core.py)

### 4-1. 手順

1. 動作温度推定（T_top, T_bot を Clausius-Clapeyron で初期推定）
2. K 値・相対揮発度 α 計算（PR EOS または CC）
3. **Fenske**: 最小段数 N_min
4. 物質収支: Fenske split で全成分の塔頂・塔底分配
5. **Underwood**: 最小還流比 R_min（塔頂組成が必要なので物質収支後）
6. **Gilliland (Eduljee)**: feasibility 確認（N ≥ N_min, R ≥ R_min）
7. **Kirkbride**: 推奨フィード段
8. Q_cond, Q_reb, Q_feed_preheat 計算
9. 塔径・塔高計算 → 10. CAPEX（Vessel + Trays + Cond + Reb）→ 11. infeasible ならペナルティ

### 4-2. Fenske（`_fenske_N_min`）

$$N_{min} = \frac{\log\!\left[\dfrac{r_{LK}}{1-r_{LK}}\cdot\dfrac{r_{HK}}{1-r_{HK}}\right]}{\log \alpha_{LK,HK}}$$

$r_{LK}$ = LK 塔頂回収率、$r_{HK}$ = HK 塔底回収率。$\alpha \le 1$ なら分離不能（∞）。

### 4-3. Underwood（`_underwood_R_min`）

$$\sum_i \frac{\alpha_i z_i}{\alpha_i - \theta} = 1 - q
\quad(\alpha_{HK} < \theta < \alpha_{LK}\ \text{を brentq で}),
\qquad R_{min} = \sum_i \frac{\alpha_i x_{top,i}}{\alpha_i - \theta} - 1$$

### 4-4. Gilliland-Eduljee（`_gilliland_eduljee`）

$$X = \frac{R - R_{min}}{R + 1}, \qquad Y = 0.75\left(1 - X^{0.5668}\right), \qquad N = \frac{Y + N_{min}}{1 - Y}$$

### 4-5. Kirkbride（`_kirkbride_feed_stage`）

$$\log_{10}\!\frac{N_{above}}{N_{below}} = 0.206\,\log_{10}\!\left[\frac{z_{HK}}{z_{LK}}\left(\frac{x_{LK,bot}}{x_{HK,top}}\right)^2\frac{B}{D}\right]$$

返り値は塔頂から数えたフィード段番号（1=塔頂）。

### 4-6. K 値・α

`K_method='pr'`（既定）は PR EOS（`_K_pr`、`src/eos.fugacity_coeff` 経由の $K=\varphi^L/\varphi^V$）、
`'cc'` は Clausius-Clapeyron（`_K_cc`、大気圧沸点 `_T_BOIL_ATM` と蒸発潜熱 `_LAMBDA_KJ` から）。
$\alpha_i = K_i / K_{HK}$。泡点温度は `_bubble_T_K`（PR brentq、範囲外なら CC fallback、ガード [200K, 600K]）。

### 4-7. proxy 罰則（`_compute_proxy_penalty`）

FUG の楽観性を BO objective で是正する罰則 [億円/年]:
- **(a) C3 漏れ罰則**: partial cond の F_top に C3 が閾値超で漏れる設計に超過 pp × 係数（rigorous パスは PSA/Mem trace bypass 罰則と二重計上になるためスキップ）
- **(b) margin 罰則**: R/R_min または N/N_min が閾値未満（rigorous Wang-Henke が dT_max stall する領域）

partial cond の塔（Dist2）は非鍵分配が乖離しやすいので係数を底上げ。

---

## 5. 経路2: Wang-Henke 厳密 (distillation_rigorous.py)

### 5-1. アルゴリズム

bubble-point method（Wang-Henke）で MESH 方程式を tray-by-tray に解く反復法。`wang_henke_solve(...)`。

主要仮定:
- **CMO（Constant Molar Overflow）**: feed の上下各セクションで V, L 一定
- **K 値**: PR EOS（`src/eos.bubble_point_T` を流用、thermo と 0.02% 整合）。各段で Newton 1 ステップ
  bubble-point 更新（`_bubble_T_newton_step`、∂lnK/∂T ≈ ΔH/RT² を Clausius-Clapeyron で近似）
- **数値ガード**: 段の T を [200K, 600K] にクランプ
- **収束**: damping で T プロファイル更新、stage 毎の |ΔT| < 0.05 K で収束、最大 500 iter で打ち切り

線形系は `scipy.linalg.solve_banded` で tridiagonal を高速求解（`_solve_tridiagonal`）。
収束加速に Wegstein（retry 時は無効化 + damping 0.2 へ）。

### 5-2. condenser モデル

- **partial condenser**: Stage 1 を平衡段として扱い、x_1 と y_1 で VLE 成立。vapor distillate
  $D_V = V_1$、reflux $L = R \times D_V$。non-condensable（H2='C', CH4='E'）は液側から除くため
  $K \to 10^8$（理想化 $K\to\infty$）に上書きし、$x_1 \approx 0$ を作る。
- **total condenser**: Stage 1 は単純凝縮（平衡なし）。$x_{i,1} = y_{i,2} = K_{i,2}\,x_{i,2}$。
  tridiagonal 最上段で $B_1=-1, C_1=K_{i,2}$。

### 5-3. always-on 検証（`RigorousResult`）

収束時に毎回計算し、内部解の不整合（過去の B_bottoms バグ等）を早期検知:
- `mesh_residual_max` / `mesh_residual_mean`: 中間段 MESH 残差（F_total で正規化）
- `component_balance_max`: 成分マスバランス相対誤差 |feed − top − bot| / feed

### 5-4. narrow-α と失敗時

Dist3（α≈1.07）は damping を低めにして T 振動を抑制。200〜500 iter で収束しなければ
`converged=False`（または例外）で返し、`distillation_core` 側が dT_max を抽出して FUG にフォールバック
（warning + BO の連続シグナル化）。

---

## 6. 経路3: GPR サロゲート (distillation_sm.py)

### 6-1. モデル構造

`models/column{1,3}_sm.pkl`（dict）。プロセス内キャッシュ（`_MODEL_CACHE`）で 1 回だけロード:
- `regressors`: ターゲット別 `GaussianProcessRegressor`
- `classifier`: feasibility 判定の RandomForest Pipeline
- `x_scaler` / `y_scalers`: StandardScaler

### 6-2. 予測手順（`_predict`、HYSYS 実測と突合せ確定）

```
Xs = x_scaler.transform(X)                                  # 入力スケーリング
y  = y_scalers[t].inverse_transform(regressors[t].predict(Xs))  # ターゲット別逆変換
feasible = classifier.predict(X)                            # raw X (Pipeline 内蔵 scaler)
```

入力が `bounds` 外なら clip し `clamped=True`。

### 6-3. 入力特徴量（`input_columns`）

| 塔 | 特徴量 |
|---|---|
| Dist1 | In_Total_Stages, In_Feed_Stage, In_Column_P[kPa], In_Flow[kgmol/h], In_CompFraction2(spec) |
| Dist3 | In_Total_Stages, In_Feed_Stage, In_Column_P[kPa], In_Flow[kgmol/s], In_Propane(feed C3H8 mol frac) |

### 6-4. 物質収支の復元

組成は **binary feed 前提**（Dist1=A/Z, Dist3=A/B）で `bottom = feed − top` のマスバランス
（`_binary_split_streams`）により保存則を満たす形で復元。Dist3 SM は spec 入力を持たず、
サンプリング時の固定 spec 戦略を埋め込むため、分配（回収率）は SM 準拠を「正」とする。

### 6-5. 出力の組立

`_build_column_result` → `units.vle.hysys.provider._column_result_to_dist_result` を再利用し、
ユーティリティ選択・HE 面積・CAPEX を **HYSYS 経路と完全に揃える**。エントリは
`solve_column1_via_sm` / `solve_column3_via_sm`（汎用 `simulate_sm` はスタブで NotImplementedError、
column{1,3}.py ラッパから直接呼ぶのが標準）。

---

## 7. 塔径・塔高・CAPEX・ユーティリティ

（3 経路共通の後処理。コンテスト Ver.2.0 §4-2 準拠）

### 7-1. 塔径（フラッディング許容蒸気速度 G*）

$$G^* = SF \cdot K \cdot \sqrt{\rho_v(\rho_l - \rho_v)} \quad [\mathrm{kg/(m^2 s)}]$$

塔頂蒸気質量流量を G* で割って断面積→塔径 $D_{col}$ を決める（SF=0.8, K=0.05）。

### 7-2. 塔高

$$H_{col} = \frac{N_{theoretical}}{\eta_{tray}}\cdot \text{tray\_spacing} + \text{top\_section} + \text{bot\_section}$$

段効率 80%、段間隔 0.6m、塔頂 +2m、塔底 +4m。

### 7-3. CAPEX

- 塔本体: 縦型プロセス容器の Bare Module（`_vessel_capex_okuyen` → `cost_calculator`）
- トレイ: Sieve tray（`calc_tray_capex_okuyen`、段数係数 Fq 付き）
- コンデンサ/リボイラ: HE 面積を U 値（`flowsheet.heat_integration.lookup_U`）と LMTD から算出し
  `calc_he_capex_okuyen`
- ユーティリティは `src.utility_selector.select_cooling_utility` / `select_heating_utility` で
  塔頂/塔底温度から tier と単価を選択（`Q_cond`/`Q_reb`、`*_utility_name`/`*_jpy_per_GJ`）

### 7-4. 物性

- 大気圧沸点 `_T_BOIL_ATM` [K] / 蒸発潜熱 `_LAMBDA_KJ` [kJ/mol]（CC 用、`distillation_core` と
  `distillation_rigorous` に同値を独立定義）
- PR EOS の Tc/Pc/ω は `src/config.py THERMO_DATA`

---

## 8. 入力・出力

```python
from src.distillation_core import (
    DistDesignVars, DistFixedParams, simulate_distillation_column,
)
result = simulate_distillation_column(design, feed, fixed)   # feed: ProcessStream
```

`DistResult`（抜粋）:

| フィールド | 説明 |
|---|---|
| `top` / `bottom` | 塔頂・塔底 `ProcessStream` |
| `equipment` (`DistEquipment`) | D_col, H_col, V_col, CAPEX 内訳, Q_cond/Q_reb, N_min/R_min, N_feed_kirkbride, feasible, T_top/T_bot, A_cond/A_reb, ユーティリティ名・単価, proxy_penalty_okuyen, N_needed, dT_max_rigorous |

`RigorousResult`（rigorous 経路の内部結果、`DistResult` に詰め直し）: converged, n_iter, message,
T_profile_K, x/y/K_profile, F_top/F_bot, V/L_top/bot_kmolh, mesh_residual_*、component_balance_max。

---

## 9. 主要パラメータ

| 定数 | 値 | 場所 | 説明 |
|---|---|---|---|
| Gilliland 指数 | 0.5668 / 係数 0.75 | core | Eduljee 形式 |
| Kirkbride 係数 | 0.206 | core | フィード段 |
| `_T_COND_MIN` | 313.15 K (40°C) | core | 冷却水使用時の凝縮器最低温度 |
| `_T_BUBBLE_MIN/MAX` | 200 / 600 K | core/rigorous | 泡点 brentq ガード（物性値でなく数値安全域） |
| `_PENALTY` | 1e9 億円 | core | infeasible シグナル |
| max_iter | 500 | rigorous | 反復打ち切り |
| T_tol | 0.05 K | rigorous | 段 T 収束 |
| damping | 0.5（retry 0.2） | rigorous | T 更新重み |
| `_K_NONCONDENSABLE_LARGE` | 1e8 | rigorous | partial cond の non-condensable K |
| `_NON_CONDENSABLE_KEYS` | {'C','E'} | rigorous | H2, CH4 |
| `tray_efficiency` | 0.8 | core (FixedParams) | 段効率（contest §4-2） |
| `SF` / `K_factor` | 0.8 / 0.05 | core (FixedParams) | G*（contest §4-2） |

---

## 10. 出典

- **Fenske, M.R. (1932)** *Ind. Eng. Chem.* 24, 482-485 — N_min。
- **Underwood, A.J.V. (1948)** *Chem. Eng. Prog.* 44(8), 603-614 — R_min。
- **Eduljee, H.E. (1975)** *Hydrocarbon Processing* 54(9), 120 — Gilliland 相関（式形）。
- **Kirkbride, C.G. (1944)** *Petroleum Refiner* 23(9), 321 — フィード段位置。
- **Seader, Henley & Roper, "Separation Process Principles", 3rd ed., Ch.10.4** — Wang-Henke (bubble-point method)。
- **プロセス設計授業資料 R08-3.pdf 付録 A** — Sieve tray K1/K2/K3、圧力係数 Fp、据付間接費 1.18、CEPCI 補正。
- **第17回プロセスデザイン学生コンテスト Ver.2.0 §4-2** — G*、段間隔 0.6m、段効率 80%、塔頂 2m + 塔底 4m。
- K 値検証: CalebBell/thermo（MIT v0.6.0）と `src/eos.py` の PR EOS が 0.02% 一致を確認。

---

## 11. 既知の限界・仮置き

- **成分 Cp / 蒸発潜熱は範囲固定値**（`src/component_data.py`、CC 用の `_LAMBDA_KJ`）。温度依存の
  Shomate/Watson 化は未実施（蒸留段温度域では実用上十分と判断、`!仮置き` マーカー残置）。
- **FUG は楽観的**。narrow-margin（R/R_min, N/N_min ≪ 1.5）や partial cond の C3 漏れで feasible でも
  rigorous で詰むことがある。proxy 罰則で補正するが、罰則係数・閾値の感度解析は未記録。
- **Wang-Henke 数値ガード**: T 収束 tol 0.05 K、dT_max/dT_floor 20/0.01 K、T クランプ [200,600] K は
  いずれも出典なし（数値安定化のための判断）。Dist2 partial_condenser stage 1 の T に brentq 偽根由来の
  ズレが残る既知問題（thermo 採用で大幅緩和）。
- **SM サロゲート**: binary feed 前提（Dist1=A/Z, Dist3=A/B）で組成を復元。多成分フィードや学習範囲外
  （bounds clip 発火）では精度低下。Dist3 SM は分配を spec ではなくモデル準拠で「正」とする近似。
- **CEPCI_CURRENT=800** は !仮置き 推定値（CAPEX 全体に effected）。

> 上記「!仮置き」「【確認中】」項目の最新状況・置換ルールは
> [`../KNOWN_PLACEHOLDERS.md`](../KNOWN_PLACEHOLDERS.md)（§1.4 物性、§3 数値ガード=`_T_COOLING_FLOOR`/
> Wang-Henke T tol/dT_max/dT_floor/`_T_COND_MIN`、§A 蒸留 rigorous アルゴリズム/K 値検証）を参照。
