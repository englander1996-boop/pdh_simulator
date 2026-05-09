# SPEC: units/utils/ — 共通ユーティリティユニット (5 本)

**ファイルパス**:
- `units/utils/pump.py`
- `units/utils/compressor.py`
- `units/utils/cooler.py`
- `units/utils/mixer.py`
- `units/utils/expansion_valve.py`

**最終更新**: 2026-05-09

---

## 0. 共通

- 全モジュールが `stream.stream.ProcessStream` を入出力ベース型として使う
  - `F_in: dict[str, float]`  各成分流量 [kmol/h]、キーは `'A'..'F','Z'`
  - `T_in: float [K]`、`P_in: float [Pa]`
- CAPEX 計算は `src/cost_calculator.py` の Bare Module Cost 法ヘルパに集約

成分マッピング (config.py 準拠): A=C3H8, B=C3H6, C=H2, D=C2H4, E=CH4, F=C2H6, Z=C4H10

---

## 1. pump.py — 液体ポンプ (Centrifugal Pump)

### 目的
液→液の昇圧。contest §3-3-3「加圧すべき箇所には、ポンプ(液) を入れること」。

### IF
| 入力 | 説明 |
|---|---|
| `stream: ProcessStream` | 入口 (液相想定) |
| `P_out_target: float [Pa]` | 出口圧力 (`> P_in` 必須、ValueError) |
| `eta_pump: float = 0.70` | ポンプ効率 |
| 出力 | `PumpResult(outlet, equipment.W_kW/rho_liq/V_dot/CAPEX)` |

### 物理モデル
```
ρ_liq = liquid_density_mix(F_in)        [kg/m³]   (component_data.py)
V_dot = mass_flow / ρ_liq / 3600        [m³/s]
W_kW  = V_dot × ΔP / η_pump / 1000      [kW]
T_out = T_in                            (非圧縮性、粘性発熱無視)
```

### CAPEX
`calc_pump_capex_okuyen(W_kW, P_out_target)` (FM_PUMP=1.0、CS)

### 仮定 / 仮置き
| 項目 | 値 | 備考 |
|---|---|---|
| `eta_pump` | 0.70 | 化工便覧 改訂六版 5·6·4 項【例題 5·8】(η=0.7 で軸馬力算出例) |
| 液密度 | C3H8/C4H10 のみ精緻化 | 他成分はフォールバック値 |
| 粘性発熱 | 無視 | LPG 系で ΔP=7bar 程度なら数 K 以内 |

### 呼び出し箇所
- `flowsheet/run_one_pass.py`: Pump1 (Fresh LPG 10bar → Dist1 17bar)

---

## 2. compressor.py — 遠心式圧縮機

### 目的
ガス→ガスの昇圧 (例: 反応器出口 0.5 bar → Dist2 8.5 bar)。
ポリトロピック圧縮式で W と T_out を推算。

### IF
| 入力 | 説明 |
|---|---|
| `stream: ProcessStream` | 入口 (ガス相) |
| `P_out_target: float [Pa]` | 出口圧力 (`> P_in` 必須) |
| `eta_poly: float = 0.75` | ポリトロピック効率 |
| 出力 | `CompressorResult(outlet, equipment.W_kW/T_out/CAPEX)` |

### 物理モデル
```
γ_mix  = Σ x_i × γ_i                                (流量加重平均、_GAMMA テーブル)
n      = γ / (γ - (γ-1)/η_poly)                     (ポリトロピック指数)
ratio  = P_out / P_in
W_mol  = n/(n-1) × R × T_in × (ratio^((n-1)/n) - 1) [J/mol]
T_out  = T_in × ratio^((n-1)/n)
W_kW   = W_mol × F_mol_s / 1000
```

成分別 γ (粗近似): A=1.13, B=1.15, C=1.40, D=1.18, E=1.32, F=1.13, Z=1.10

### CAPEX
`calc_comp_capex_okuyen(W_kW)` (遠心式、CS)

### 仮定 / 仮置き
| 項目 | 値 | 備考 |
|---|---|---|
| `eta_poly` | 0.75 | 化工便覧 改訂六版 p.333 (0.7〜0.8 中央値) |
| `γ_i` | 上記表 | 温度依存性無視。NIST WebBook で精緻化可能 |
| 多段化判定 | 呼び出し側 | 圧縮比 > 4 程度なら 2 段+段間冷却を呼び出し側で組む |

### 呼び出し箇所
- `flowsheet/run_one_pass.py`: Comp2a, Comp2b (反応器出口の昇圧、等圧縮比 √17 で 2 段)
- `units/separators/membrane/membrane_system.py`: 内部で `compress_isentropic` 経由

---

## 3. cooler.py — 冷却器 / 加熱器 (Heat Exchanger)

### 目的
ストリームを目標温度まで温度制御。冷却にも加熱にも使う。
ターゲット温度から自動的に冷媒/熱媒を選択し (`utility_selector.py`)、equipment に
utility 名・単価を埋め込む (economics.py が直接読む)。

### IF
| 入力 | 説明 |
|---|---|
| `stream: ProcessStream` | 入口 |
| `T_out_target: float [K]` | 出口目標温度 |
| `P_out: float [Pa] \| None` | None で圧損なし |
| `dT_lm: float = 30` | 対数平均温度差 |
| `phase_change: bool = False` | True で潜熱を Q に加算 |
| `process_phase: str = GAS` | プロセス側顕熱区間の相 (StreamPhase.GAS/LIQUID)、§4-4 U 表索引用 |
| 出力 `equipment.Q_duty_kW` | 顕熱+潜熱合計 (符号: 負=冷却, 正=加熱) |
| 出力 `equipment.Q_sensible_kW`/`Q_latent_kW` | 内訳 |
| 出力 `equipment.A_est_m2`, `CAPEX` | 伝熱面積、設備費 |
| 出力 `equipment.utility_name`, `utility_jpy_per_GJ` | 選択ユーティリティ |

### 物理モデル
```
Q_sensible = Σ_k F_k × Cp_k × (T_out - T_in)             [kW]
Q_latent   = Σ_k F_k × λ_k    (phase_change=True かつ T_out>T_in のときのみ)
A          = max(|Q| × 1000 / (U × dT_lm), 10)           [m²]
```

VLE を持たないため、相変化の有無は呼び出し側で `phase_change=True` を渡して指示
(例: `mem_precool` の液→気フィード化)。

### CAPEX / OPEX
- CAPEX: `calc_he_capex_okuyen(A_m2)` (固定管板式 HE)
- OPEX: `economics.py` で `_heat(|Q|, equipment.utility_jpy_per_GJ)` を計上

### 仮定 / 仮置き
| 項目 | 値 | 備考 |
|---|---|---|
| U 値 | contest §4-4 表 (lookup_U) で (hot_phase, cold_phase) 別に決定。相変化区間は分割 | 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 |
| `dT_lm` | 30 K | 代表値。寒冷側で温度クロスする条件は呼び出し側で確認 |
| `A_min` | 10 m² | CAPEX 下限保護 (HE 適用範囲は 10〜1000 m²) |
| 凝縮潜熱 | 加熱側のみ計上 | T_out<T_in (凝縮) は現状無視 (Mem 内部で別途計算済み) |
| `LATENT_HEAT_KJ_PER_KMOL` | 沸点での値固定 | `component_data.py`、Watson 式で温度補正可能 |

### 呼び出し箇所
- `flowsheet/run_one_pass.py`: Cooler (反応器出口冷却), Intercool (Comp 段間冷却), MemPrecool (Mem フィードガス化)

---

## 4. mixer.py — ストリームミキサー

### 目的
複数ストリームの合流。エンタルピーバランスで T_out を解く (旧版のモル流量加重平均は
Cp 異成分や温度差が大きい場合に物理的に不正確だった)。

### IF
| 入力 | 説明 |
|---|---|
| `streams: list[ProcessStream]` | 1 本以上 |
| 出力 | `ProcessStream` (組成: 成分別加算、温度: エンタルピー保存、圧力: 最低値) |

### 物理モデル
```
F_out_k = Σ_i F_i_k                                              (成分別加算)
T_out   = Σ_i Σ_k F_i_k × Cp_k × T_i  /  Σ_k F_out_k × Cp_k     (T_ref キャンセル)
P_out   = min_i P_i                                              (高圧側は制御弁でドロップ)
```

### 仮定
| 項目 | 値 | 備考 |
|---|---|---|
| Cp | `cp_of()` の固定値 | `component_data.py`、温度依存性無視 |
| 圧損 | 高圧側は弁で min まで降圧 | 弁の CAPEX/OPEX は計上しない |
| 相変化 | 無視 | 部分凝縮を伴う混合は精度低下、本フロー内 (Reactor 入口) は気相同士なので OK |

### 呼び出し箇所
- `flowsheet/run_one_pass.py`: 反応器入口 (Dist1 塔頂 + Recycle Dist3 + Recycle Mem)

---

## 5. expansion_valve.py — Joule-Thomson 等エンタルピー膨張

### 目的
配管中の減圧操作で起こる JT 効果 (温度低下) を陽に扱う。旧 run_one_pass は P を
書き換えるだけで T 維持していたため反応器入口プレヒート Q_preheat が過小評価だった。

### IF
| 入力 | 説明 |
|---|---|
| `stream: ProcessStream` | 入口 |
| `P_out: float [Pa]` | 出口圧力 (`< P_in` 必須、圧縮方向は ValueError) |
| `T_search_lower_K: float = 100` | T_out 探索下限 |
| 出力 | `ProcessStream` (F_in 維持、T を JT 等エンタルピーで更新、P=P_out) |

### 物理モデル
```
等エンタルピー条件: H(T_in, P_in) = H(T_out, P_out)
H(T,P) = H_ig(T) + H_residual(T,P)              (PR EOS、src/eos.py)
        ⇒ brentq で T_out を解く
```

### 仮定 (重要)
| 項目 | 状態 | 備考 / 限界 |
|---|---|---|
| **vapor 相のまま膨張** | 仮定 | 部分気化は扱わない。飽和近い流体 (例: 17bar dist1_top_rx) は実際は部分気化で温度低下が緩和される。**温度低下を過大評価する側に偏る**。将来 VLE フラッシュ (Rachford-Rice) を `src/eos.py` に追加して二相膨張対応 |
| 探索失敗時 | T_in 維持 | brentq 失敗時はワーニング出して理想気体仮定 (T 維持) で返す |

### 機器コスト
**なし**: 配管中の絞り弁扱い、CAPEX/OPEX に計上しない (run_one_pass.py のコメント参照)。

### 呼び出し箇所
- `flowsheet/run_one_pass.py`:
  - Dist1 塔頂 17bar → 反応器入口 0.5bar
  - Recycle Dist3 (20bar) → 0.5bar
  - Recycle Mem (9.5bar) → 0.5bar

---

## 共通の依存

- `stream/stream.py`: `ProcessStream`
- `src/cost_calculator.py`: `calc_pump/comp/he_capex_okuyen`
- `src/component_data.py`: `MW`, `cp_of`, `CP_DEFAULT`, `liquid_density_mix`, `LATENT_HEAT_KJ_PER_KMOL`
- `src/utility_selector.py`: `select_utility` (cooler のみ)
- `src/eos.py`: `z_factor`, `residual_enthalpy`, `_dh_ig` (expansion_valve のみ)
- `src/cost_parameters.py`: 単価系定数 (経由は `cost_calculator.py`)
