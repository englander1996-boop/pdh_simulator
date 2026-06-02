# SPEC: heat_integration.py — ヒートインテグレーション (ピンチ解析 Stage1)

**ファイルパス**: `flowsheet/heat_integration.py`

---

## 1. 目的

フローシートの全熱交換ストリームに対し **ピンチ解析の targeting (Stage 1)** を行い、
熱力学的限界値 (最少加熱量 Q_H_min / 最少冷却量 Q_C_min = MER、ピンチ温度、理想総
伝熱面積、最少熱交換器数、ユーティリティ tier 別 Q 配分) を計算する。BO ループ内では
**離散的な流体ペアリング (HE network 合成 = Stage 2) は行わず**、滑らかで決定論的な
ターゲティング値だけを返す (GP に不向きな組合せ問題を回避)。結果は `economics.py` 側で
HI 後 OPEX に反映される。

Stage 2 (実 HEN 合成) は不採用。HI は Stage 1 (pinch targeting) のみで評価する。

---

## 2. モデル理論 (式)

### 2-1. ストリーム展開

各 `HIStream` を顕熱区間 + 潜熱区間に分解 (`_expand_stream`):

- `is_hot` は `T_in > T_out` で自動判定 (与熱流体 = 冷却される側 = hot)。
- 顕熱: `Q_sensible = F_Cp × |T_out − T_in|`
- 潜熱: 単一の相変化温度 `T_phase_K` に `Q_latent_kW` を集中させる簡易モデル
  (潜熱区間は `T_top == T_bot` の水平線、`F_Cp = 0`)。

### 2-2. Problem Table Algorithm (`_problem_table`)

ΔTmin が決まれば MER は一意。

1. 全 segment 温度をシフト: hot は `−ΔTmin/2`、cold は `+ΔTmin/2`。
2. シフト温度のユニーク値を降順ソート → 区間境界。
3. 各区間の熱不足量:
   `deficit = (ΣF_Cp_cold − ΣF_Cp_hot)·ΔT + ΣQ_latent_cold − ΣQ_latent_hot`
   (正 = 外部加熱が必要、負 = 過剰)。
4. 高温から累積 (外部加熱 0 開始) → 最大正値が **Q_H_min**。
5. Q_H_min を最高温で投入後に再積算 → GCC (グランドコンポジットカーブ、
   `(T_shifted, H_residual)`)。末尾値 = **Q_C_min**。`H_residual ≈ 0` の点 =
   ピンチ温度 (シフト)。
6. 実ピンチ温度: `T_pinch_hot = T_pinch_shifted + ΔTmin/2`,
   `T_pinch_cold = T_pinch_shifted − ΔTmin/2`。

### 2-3. 複合線 (`_build_composite`)

与熱/受熱の `(H_kW, T_K)` 昇順折れ線。潜熱区間は同じ温度で H が水平に伸びる。

### 2-4. 総伝熱面積 A_total — Bath 式 (`_calc_A_total`)

```
A_total = Σ_k (1/ΔT_LM,k) × Σ_i (q_i,k / h_i)
```

各シフト温度区間 k で向流の対数平均温度差 `ΔT_LM = (ΔT1 − ΔT2)/ln(ΔT1/ΔT2)`
(`ΔT1 = T_hot_high − T_cold_high`, `ΔT2 = T_hot_low − T_cold_low`)。区間内で
交換可能な量 `Q_exch = min(Q_hot_total, Q_cold_total)` にスケールして加算。

### 2-5. 最少熱交換器数 — Linnhoff 式 (`_calc_N_HE_min`)

```
ピンチなし (閾値問題): N_HE,min = N_streams + N_utilities − 1
ピンチあり          : 上下分割で N_streams + N_utilities − 2 (簡易)
```

### 2-6. ユーティリティ tier 別 Q 配分 (`_assign_utilities_to_tiers`)

GCC を各 tier の有効シフト温度で水平に切り、Linnhoff 標準手法で多 tier 配分:

- 加熱 (ピンチ上): 安価順 = `supply_T` 昇順 (LP < MP < HP < 燃料)。
  各 tier の `T_eff_shifted = supply_T − ΔTmin/2` まで「使えるだけ使う」。
- 冷却 (ピンチ下): 安価順 = `supply_T` 降順 (空冷 → 冷却水 → プロピレン → エチレン)。
  `T_eff_shifted = supply_T + ΔTmin/2`。
- GCC pocket による過剰計上は `_normalize` で合計を Q_H_min / Q_C_min に揃える
  (pocket は内部熱交換で消化扱い)。

### 2-7. 総括熱伝達係数 U (contest §4-4 表)

`U_TABLE_CONTEST[(hot_phase, cold_phase)]` を引く (`lookup_U`)。流速によらず固定 U
を使用。相分類 `StreamPhase` = GAS / LIQUID / CONDENSING / EVAPORATING。

| 相ペア | U [W/(m²·K)] |
|---|---|
| ガス-ガス | 150 |
| 液-ガス / ガス-液 | 200 |
| 液-液 | 300 |
| ガス凝縮-ガス / ガス-液蒸発 | 500 |
| ガス凝縮-液 / 液-液蒸発 | 1000 |
| ガス凝縮-液蒸発 | 1500 |
| 表外 (フォールバック) | 150 |

Bath 式の `h` は `U/2` を仮置き (`1/U = 1/h_h + 1/h_c`、両側等価仮定で `h ≈ 2U`)。

---

## 3. 入出力

### 公開 API

- `extract_streams(one_pass, swing_T_in) → List[HIStream]`
  run_one_pass 結果から HI 対象ストリームを抽出。各ユニットの T_in/T_out/Q を
  equipment から直接読み、新たな仮定を導入しない。Q が小さい流れは自動除外。
  抽出される代表ストリーム:
  - H1 Cooler (反応器出口冷却), H2 Comp2 段間冷却, H3 Desuper
  - C1 mem_precool (顕熱+潜熱), C2 反応器予熱
  - 各塔 condenser (hot 潜熱) / reboiler (cold 潜熱) / feed preheat (顕熱)
  - 膜 vaporizer / 製品冷却器 (顕熱+潜熱) / 圧縮機段間冷却 / PSA 予熱 (符号で hot/cold)
- `pinch_analysis(streams, dT_min_K=10.0, heating_tiers, cooling_tiers) → HIResult`
- `get_default_utility_tiers() → (heating, cooling)`
  (`src.utility_selector` の tier を `UtilityTier` 型に変換)
- `apply_hi_to_economics(economics, hi_result, ...) → Economics`
  熱系 OPEX を HI tier 別 OPEX (`HI: <tier>` キー) に置換した新 Economics を返す。
  非熱系 OPEX (電力・触媒・吸着剤・原料費) と Revenue・CAPEX は据え置き。Hasebe 集計項
  は C_UT 変化を反映して再計算。CAPEX は据え置き (targeting only)。
- `classify_heat_opex_key`, `summarize_current_heat_opex`, `calc_hi_opex_okuyen`
  (economics.opex キー ⇄ 熱種別 hot/cold の対応、二重計上防止に使用)
- 可視化: `plot_TQ` (T-Q 線図), `plot_GCC` (グランドコンポジットカーブ)

### `HIResult` の主フィールド

`Q_H_min_kW`, `Q_C_min_kW`, `T_pinch_hot_K`, `T_pinch_cold_K`, `A_total_m2`,
`N_HE_min`, `utility_breakdown {tier: Q_kW}`, `composite_hot/cold`, `GCC`,
`feasible`, `message`。

---

## 4. 主要パラメータ

| 項目 | 値 | 意味 |
|---|---|---|
| `dT_min_K` | 10.0 (既定) | 最小接近温度差 (教科書標準)。BO 設計変数に含めず固定が一般的。 |
| U 表 | contest §4-4 | 9 区分の固定 U |
| `H_*_W_m2K` | 300/600/1000 | Bath 式の h 代表値 (≈2U) |
| 加熱 tier | LP/MP/HP Steam, 燃料燃焼 | `src.utility_selector._HEATING_TIERS` |
| 冷却 tier | 空冷, 冷却水, プロピレン冷媒, エチレン冷媒 | `src.utility_selector._COOLING_TIERS` |

---

## 5. 出典

- [1] 長谷部 伸治・外輪 健一郎『プロセスシステム工学 (No.4) — 熱交換器ネットワーク
  の最適合成』京都大学講義資料 令和7年度 (2025)。§4.4 最小接近温度差と最少必要
  加熱・冷却量、§4.5 熱複合線図の分解、§4.6 グランドコンポジットカーブ。
  例題 4.2/4.3 は `tests/test_heat_integration.py` で再現検証済み。
- [2] Linnhoff B., Hindmarsh E., "The pinch design method for heat exchanger
  networks," *Chem. Eng. Sci.* 38 (1983) 745-763。多 tier 配分・最少 HE 数式
  `N_HE = N_streams + N_utilities − 1` の出典。
- U 表: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4。

---

## 6. 既知の限界・仮置き

- **潜熱は単一温度集中の簡易モデル**。VLE フル対応 (温度依存の相変化曲線) は将来課題。
- `N_HE_min` のピンチ分割は簡易式 (上下それぞれのアクティブストリーム数の厳密カウントを
  していない)。
- Bath 式の `h ≈ 2U` 仮置きは targeting 推算の精度範囲。
- tier 配分は GCC pocket を内部熱交換で消化扱いとし `_normalize` で合計補正する近似。
- HI は **targeting only** (Stage 1)。実 HE network・追加 HE CAPEX は計算せず、
  理論限界 OPEX のみ反映 (CAPEX 据え置き)。Stage 2 は不採用。
- `dT_min_K = 10` は教科書標準の固定値 (最適化変数にしていない)。
