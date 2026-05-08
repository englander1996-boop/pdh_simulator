# SPEC: column1/2/3.py — FUG ベース蒸留塔 3 本

**ファイルパス**:
- `units/separators/column1/column1.py` — 脱ブタン塔
- `units/separators/column2/column2.py` — 脱エタン塔
- `units/separators/column3/column3.py` — C3 スプリッタ

**共通エンジン**: `src/distillation_core.py`
**最終更新**: 2026-05-09

---

## 0. 共通設計

3 塔とも `src/distillation_core.simulate_distillation_column(design, feed, fixed)` を
呼ぶ薄いラッパ。各塔ファイルで「塔別固定設定」(LK/HK 成分、デフォルト P/N/R) のみ
定義する Option C ハイブリッド構造。

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
    feed:   ProcessStream,
    design: DistDesignVars | None = None,    # None で各塔のデフォルト
    fixed:  DistFixedParams | None = None,
) -> DistResult:
```

### `DistDesignVars` (BO 探索対象 + 仕様)
| フィールド | 種別 | 説明 |
|---|---|---|
| `P_col` | BO | 塔操作圧力 [Pa] |
| `N_stages` | BO | 理論段数 |
| `N_feed` | BO | フィード段位置 (塔頂=1, 塔底=N) |
| `reflux_ratio` | BO | R = L/D |
| `LK`, `HK` | 仕様 | 軽キー / 重キー成分 (塔別ラッパで固定) |
| `recovery_LK_top`, `recovery_HK_bot` | 仕様 | キー成分回収率 (デフォルト 0.99) |
| `K_method` | 仕様 | `'pr'` (PR EOS) または `'cc'` (Clausius-Clapeyron) |
| `q` | 仕様 | フィード状態 (1=飽和液, 0=飽和気) |

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

### 目的
反応器出口の冷却・圧縮ガス (8.5 bar, 47°C ガス相) から軽質ガスを分離。
塔頂 (H2 + CH4 + C2H4 + 微量 C3) → PSA
塔底 (C3H8 + C3H6 + C2H6) → Membrane

### デフォルト設計
| 項目 | 値 | 備考 |
|---|---|---|
| `P_col` | 8.5 bar | Mem の P_H ≤ 9.5 bar 制約に合わせ低圧運転。差圧 1bar マージン |
| `N_stages` / `N_feed` | 20 / 10 | N_min ≈ 1〜1.4 (α 大なので余裕大) |
| `reflux_ratio` | **6.0** | PR R_min は供給組成依存 (1.5〜4.8 のレンジ)。single-pass 想定 (z_LK=0.26%) で R_min=4.8、margin 1.25。**旧 4.5 は CC 基準で過大設計だったが、PR でも feasibility ぎりぎり** |
| `LK` / `HK` | D (C2H4) / A (C3H8) | C2H6 を key にすると Underwood 過大評価、C2H4 を LK にして R_min 現実値 |
| `K_method` | **'pr'** | 'pr' に切替 (2026-05-09)。x_top に H2/CH4 主体だと CC へ自動フォールバック |
| `q` | 0.0 | 気フィード |

### 物性
- T_top ≈ 40°C (Comp2b 出口を冷却して反応器圧縮済みガスの露点)
- T_bot ≈ 16°C (8.5bar 下の C3 主体液)

### 注意点
single-pass (z_LK=0.26%) と recycle (z_LK 大) で R_min が大きく変動するため、
最適化器は R_min を見ながら R を振る必要がある。

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
| `LK` / `HK` | B (C3H6) / A (C3H8) | α ≈ 1.07〜1.10 @20bar |
| `recovery_LK_top` / `recovery_HK_bot` | 0.99 / 0.99 | C3H6 製品 99wt% 仕様、塔底リサイクル 99% |
| `K_method` | **'pr'** | α 極小なので精度必須 |
| `q` | 1.0 | 飽和液フィード |

### 物性
- T_top ≈ 49°C, T_bot ≈ 57°C
- Q_reb ≈ 80 MW (実プラント規模)、**OPEX の主要因** (exp1 で 42.6 億円/年)

### 設計感度
α 極小ゆえ R_min/N_min が α に超敏感。BO で R を振る場合、下限 11.0 程度
(margin 1.09) を限度にすること。

---

## 共通の依存

- `src/distillation_core.py`: `simulate_distillation_column`, `DistDesignVars/FixedParams/Result/Equipment`
- `src/eos.py` (PR モード時): `bubble_point_T`, `z_factor`, `fugacity_coeff`
- `src/cost_calculator.py`: `calc_cp0`, `calc_fp`, `calc_tray_capex_okuyen`
- `src/cost_parameters.py`: `B1`, `B2`, `FM`, `CEPCI_*`, `USD_TO_JPY`, `PLANT_INDIRECT_FACTOR`
- `src/component_data.py`: `cp_of`, `CP_DEFAULT`
- `src/config.py`: `THERMO_DATA` (PR Tc/Pc/ω)

## 仮置き

| 項目 | 値 | ファイル | 備考 |
|---|---|---|---|
| 物性データ (Tc/Pc/ω) | NIST | `src/config.py` | C4H10 (Z) は化工便覧 改訂六版 表1.3 No.181 |
| `_T_BOIL_ATM`, `_LAMBDA_KJ` | 代表値固定 | `src/distillation_core.py:68-77` | CC モード時に使用、PR モード時はフォールバック専用 |
| Q_reb = 1.05 × Q_cond | 5% 損失仮定 | `src/distillation_core.py:554` | 文献根拠未確定、将来課題 |
| Q_feed_preheat | 顕熱のみ、液仮定 | `src/distillation_core.py:558-564` | フィード相が混在する場合に過大/過小推算 |
