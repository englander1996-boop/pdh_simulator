# SPEC: column1/2/3.py — FUG ベース蒸留塔 3 本

**ファイルパス**:
- `units/separators/column1/column1.py` — 脱ブタン塔
- `units/separators/column2/column2.py` — 脱エタン塔
- `units/separators/column3/column3.py` — C3 スプリッタ

**共通エンジン**: `src/distillation_core.py`（FUG/rigorous）/ `src/distillation_sm.py`（SM）/ `units/vle/hysys/`（HYSYS COM）
**最終更新**: 現行コードへ整合（solver_method 4 経路・Dist2 真正 deethanizer 構成・Dist3 動的 recovery を反映）

---

## 0. 共通設計

3 塔とも `simulate_columnX(feed, tunables, fixed)` という薄いラッパ。各塔ファイルで
「塔別固定設定」(LK/HK 成分、回収率、K_method、q、partial_condenser、デフォルト P/N/R)
のみ定義する Option C ハイブリッド構造。

### solver_method による経路分岐

`ColumnTunables.solver_method` で計算経路を選ぶ（各ラッパ冒頭で dispatch）:

| solver_method | 経路 | 実装 | 用途 |
|---|---|---|---|
| `'fug'`（既定） | Fenske-Underwood-Gilliland shortcut | `src/distillation_core.simulate_distillation_column` | 高速、BO 用 |
| `'rigorous'` | VLE tray-by-tray（Wang-Henke） | `src/distillation_core` | 厳密、top-k 再評価用 |
| `'sm'` | 学習済み GPR で HYSYS 解を近似 | `src/distillation_sm.solve_columnX_via_sm` | SM フォーク |
| `'hysys'` | HYSYS COM で段数別 HSC を実行 | `units/vle/hysys.solve_columnX_via_hysys` | 検証・special フォーク |

`'hysys'`/`'sm'` 経路では LK/HK/recovery は使われず、`hysys_spec_value`（塔別の主スペック）
と `hysys_feed_stage` を使う。以下 §0〜§3 の式・物性は主に **FUG 経路** の説明。

### 共通エンジンのフロー (distillation_core.py)
1. **動作温度初期推定** (CC pure component 沸点)
2. **初期 split** (CC、x=z で粗推定)
3. **反復** (最大 5 回): 塔頂/塔底それぞれの**泡点フラッシュ** (PR or CC) で K → α
   - α_geom = √(α_top × α_bot) を Fenske/Underwood に投入
   - 旧版は「平均 T で K=φ_L/φ_V」だったが単相 root で K≈1 病理 → 2026-05-09 改修
4. Fenske → N_min / Underwood → R_min / Gilliland → feasibility check
5. Kirkbride → 推奨フィード段
6. Q_cond, Q_reb, Q_feed_preheat 計算
7. 塔径・塔高、CAPEX (Vessel + Trays)
8. infeasible なら ペナルティ (`CAPEX = 1e9`) 返却

### PR 単相 root 問題のフォールバック
`_bubble_T_K` 内で泡点フラッシュ → K 計算 → K の spread が 1.1 倍未満 (= K≈1 全成分)
の場合は CC へフォールバック。Dist2 の x_top (H2/CH4 主体) で発生しがちだが安全。
詳細は `src/distillation_core.py:_bubble_T_K` のコメント参照。

### 共通インターフェース
```python
def simulate_columnX(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,   # None で各塔のデフォルト
    fixed:    DistFixedParams | None = None,
) -> DistResult:
```

### `ColumnTunables` (BO/exp で振る部分のみ)
ラッパが受け取るのは `ColumnTunables`。LK/HK・回収率・q・K_method・partial_condenser
は各塔ラッパ側で固定し、`ColumnTunables` には保持しない。

| フィールド | 種別 | 説明 |
|---|---|---|
| `P_col` | BO | 塔操作圧力 [Pa] |
| `N_stages` | BO | 理論段数 |
| `N_feed` | BO | フィード段位置（rigorous/sm では Kirkbride 推奨自動採用、本値は無視。FUG は post-hoc 出力のみ） |
| `reflux_ratio` | BO | R = L/D |
| `solver_method` | 仕様 | `'fug'`/`'rigorous'`/`'sm'`/`'hysys'`（既定 `'fug'`） |
| `recovery_LK_top`, `recovery_HK_bot` | BO/仕様 | None でラッパ既定 0.99、float で上書き |
| `D_override` | 検証 | None で FUG が Fenske split で D 決定、float で rigorous の D_total 強制（"Distillate Rate spec" 相当） |
| `hysys_spec_value`, `hysys_feed_stage` | 仕様 | HYSYS/SM 経路の主スペック値・フィード段（FUG/rigorous では None で可） |

ラッパが内部で組み立てる `DistDesignVars` には LK/HK・recovery・K_method・q が加わる。

### `DistResult.equipment` (主要)
- `D_col`, `H_col`, `V_col` 塔寸法
- `CAPEX_vessel`, `CAPEX_trays`, `CAPEX` 内訳
- `Q_cond`, `Q_reb`, `Q_feed_preheat_kW` 熱量
- `N_min`, `R_min`, `N_feed_kirkbride`, `feasible`, `message` 診断

---

## 1. column1.py — Dist1 (脱ブタン塔)

### 目的
Pump1 後の Fresh LPG (30°C, 17 bar 飽和液) から C4H10 を分離。
塔頂 (C3H8 + 微量 C4) → 反応器系 (膨張弁経由 0.5 bar)
塔底 (C4H10 + 微量 C3) → 廃棄

### デフォルト設計
| 項目 | 値 | 備考 |
|---|---|---|
| `P_col` | 17 bar | Pump1 後の圧力 |
| `N_stages` / `N_feed` | 20 / 10 | N_min ≈ 12 (PR) に対し margin 1.6× |
| `reflux_ratio` | **1.5** | PR R_min ≈ 0.95 に対し 1.57× margin (旧 0.6 は CC R_min=0.44 基準) |
| `LK` / `HK` | A (C3H8) / Z (C4H10) | α ≈ 2.3 @17bar (PR), CC は 3.4 と過大推定 |
| `K_method` | **'pr'** | 2026-05-09 復帰 |
| `q` | 1.0 | 飽和液フィード |

### 物性
- T_top ≈ 50°C (C3H8 17bar 沸点)
- T_bot ≈ 98°C (C4H10 17bar 沸点)

---

## 2. column2.py — Dist2 (脱エタン塔)

### 目的（真正 deethanizer 構成）
反応器出口の冷却・圧縮ガス (8.5 bar, 47°C ガス相) から軽質成分を分離。
塔頂 (H2 + CH4 + C2H4 + C2H6) → PSA 原料 (含むオフガス → 燃料)
塔底 (C3H8 + C3H6 の clean C3 のみ) → Membrane 系へ

### デフォルト設計
| 項目 | 値 | 備考 |
|---|---|---|
| `P_col` | 8.5 bar | Mem の P_H ≤ 9.5 bar 制約に合わせ低圧運転。差圧 1bar マージン |
| `N_stages` / `N_feed` | 20 / 10 | recovery 99/99 で N_min ≈ 7 (Fenske)、margin 約 3× |
| `reflux_ratio` | **7.0** | PR α(C2H6/C3H8 @8.5bar)≈3.9、R_min ≈ 5 (推算)、margin 1.4× |
| `LK` / `HK` | **F (C2H6) / A (C3H8)** | 実機 deethanizer 標準: LK = 最も重い軽キー (C2H6)。recovery_LK_top=0.99 で C2H6 を 99% 塔頂固定 |
| `recovery_LK_top` / `recovery_HK_bot` | 0.99 / 0.99 | 質量保存が塔内で自然に閉じる |
| `K_method` | **'pr'** | α 計算精度重視。x_top に H2/CH4 主体だと CC へ自動フォールバック |
| `q` | 0.0 | 気フィード (Desuperheater で 50°C まで冷やした飽和蒸気) |
| `partial_condenser` | **True** | H2/CH4 凝縮不能のため分流型 |

### 非キー成分の自動分配 (Fenske 後)
- C2H4(D): 塔頂 99.9% / CH4(E): 塔頂 100% / H2(C): 塔頂 100% (lighter than LK)
- C3H6(B): α/α_LK ≈ 0.29 → 塔底 98% (HK=A の隣) / C4H10(Z): 塔底 100% (Dist1 で除去済)

### partial condenser の扱い (distillation_core)
H2(C)/CH4(E) は工業冷媒 (-100°C エチレン) でも凝縮不能なので分流型で vapor distillate に
パススルー。C3 (A/B) は塔頂温度で大半凝縮するが PR EOS 実値で K≠0 のため、feed 比
`_PARTIAL_COND_C3_ROUND_FRAC`(1%) 未満の塔頂漏れのみ 0 に丸めて塔底集約 (閾値超は物理値
を残し警告)。FUG 楽観性に対し `_compute_proxy_penalty` で C3 漏れ・narrow margin 罰則を BO
objective に乗せる (partial cond は係数 ×2)。

### 物性
- T_top ≈ 40°C (反応器圧縮済みガスの露点) / T_bot ≈ 16°C (8.5bar 下の C3 主体液)

---

## 3. column3.py — Dist3 (C3 スプリッタ)

### 目的
Membrane 透過後 (20 bar, 飽和液) からポリマーグレード C3H6 を回収。
**PDH プロセスで最もエネルギー集約的**な塔 (C3H6/C3H8 揮発度差が小さい)。

### デフォルト設計
| 項目 | 値 | 備考 |
|---|---|---|
| `P_col` | 20 bar | Mem 透過側冷却で液化可能な下限近傍 |
| `N_stages` / `N_feed` | 200 / 100 | C3H6/C3H8 分離は 150〜250 段が典型 |
| `reflux_ratio` | **12.0** | PR R_min ≈ 10 に対し 1.19× margin (旧 7.7 は CC R_min=7.22 基準で PR では infeasible) |
| `LK` / `HK` | B (C3H6) / A (C3H8) | α ≈ 1.05〜1.10 @20bar |
| `recovery_LK_top` | 0.99 | C3H6 塔頂回収率 |
| `recovery_HK_bot` | **動的計算** | 製品純度 spec (C3H6 ≥ 99.5 wt%) から逆算 (下記) |
| `K_method` | **'pr'** | α 極小なので精度必須 |
| `q` | 1.0 | 飽和液フィード (Mem 透過後の冷却液) |

### recovery_HK_bot の動的計算 (Dist3 特有)
`recovery_HK_bot` を固定値にすると、フィード C3H8 量が極少 (~1%) のとき製品純度 100% の
overspec を生み N_stages が過大化する。そこでフィード組成から純度 99.5 wt% spec を満たす
最低 recovery を逆算する (`_dist3_dynamic_recovery_HK_bot`):

$$\mathrm{rec}_{HK,bot} \ge 1 - \mathrm{rec}_{LK,top}\cdot\frac{F_{LK,feed}}{F_{HK,feed}}\cdot\frac{1-p}{p}$$

- $p$ = `_DIST3_PURITY_SPEC_MOL` = 0.995 (C3H6/C3H8 は MW 差 <5% で mol≒wt 近似)
- Gilliland (Eduljee) の楽観性補正に `_DIST3_PURITY_SAFETY_FACTOR` = 1.2 (=spec margin 20%、
  経験的 magic constant、要 rigorous 比較で調整) を適用

### 物性
- T_top ≈ 49°C, T_bot ≈ 57°C
- Q_reb ≈ 80 MW (実プラント規模)、**OPEX の主要因** (exp1 で 42.6 億円/年)

### 設計感度
α 極小ゆえ R_min/N_min が α に超敏感。BO で R を振る場合、下限 11.0 程度
(margin 1.09) を限度にすること。

---

## 共通の依存

- `src/distillation_core.py`: `simulate_distillation_column`, `ColumnTunables`, `DistDesignVars/FixedParams/Result/Equipment`
- `src/distillation_sm.py` (SM モード時): `solve_columnX_via_sm`
- `units/vle/hysys/` (HYSYS モード時): `solve_columnX_via_hysys`（→ `provider`/`registry`/`session`/`adapters`）
- `src/eos.py` (PR モード時): `bubble_point_T`, `z_factor`, `fugacity_coeff`
- `src/cost_calculator.py`: `calc_cp0`, `calc_fp`, `calc_tray_capex_okuyen`
- `src/cost_parameters.py`: `B1`, `B2`, `FM`, `CEPCI_*`, `USD_TO_JPY`, `PLANT_INDIRECT_FACTOR`
- `src/component_data.py`: `cp_of`, `CP_DEFAULT`
- `src/config.py`: `THERMO_DATA` (PR Tc/Pc/ω)
- 上位: `flowsheet/run_one_pass.py`（3 塔を配管）

## 仮置き

| 項目 | 値 | ファイル | 備考 |
|---|---|---|---|
| 物性データ (Tc/Pc/ω) | NIST | `src/config.py` | C4H10 (Z) は化工便覧 改訂六版 表1.3 No.181 |
| `_T_BOIL_ATM`, `_LAMBDA_KJ` | 代表値固定 | `src/distillation_core.py` | CC モード時に使用、PR モード時はフォールバック専用 |
| Q_reb = 1.05 × Q_cond | 5% 損失仮定 | `src/distillation_core.py` | 文献根拠未確定、将来課題 |
| Q_feed_preheat | 顕熱のみ、液仮定 | `src/distillation_core.py` | フィード相が混在する場合に過大/過小推算 |
| `_DIST3_PURITY_SAFETY_FACTOR` | 1.2 | `column3.py` | Gilliland 楽観性補正の magic constant、要 rigorous 比較 |
| FUG proxy 罰則係数 | 各種 | `src/distillation_core.py` | C3 漏れ/narrow margin の罰則、BO 楽観性是正用 |

> `!仮置き` / `【確認中】` の詳細は `KNOWN_PLACEHOLDERS.md` を参照。
