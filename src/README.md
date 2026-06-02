# `src/` — 基盤層（物性・EOS・蒸留計算コア・コスト/物性パラメータ）

PDH (Propane DeHydrogenation) プロセスシミュレータの **基盤ライブラリ層**。
個々の装置モデル（`units/` 配下の反応器・蒸留塔ラッパ・PSA・膜分離器）や
フローシート求解（`flowsheet/`）、最適化（`optimization/`）が共通で参照する、
**物性データ・状態方程式・蒸留計算エンジン・コスト推算・経済/物性パラメータ** を提供する。

このディレクトリ自体は装置やフローシートを「組み立てない」。あくまで
「下から支える純関数・データ・計算コア」の集合である。

---

## 役割の整理

| レイヤ | このディレクトリでの担当 |
|---|---|
| **物理定数・反応設定** | `config.py`（気体定数・基準温度・反応速度パラメータ・成分熱力学データ `THERMO_DATA`） |
| **物性データ** | `component_data.py`（MW / Cp / 沸点 / 蒸発潜熱 / 液密度） |
| **反応速度・反応熱力学** | `kinetics.py`（PDH 3 反応の速度式）、`thermo.py`（Cp 積分・反応エンタルピー・平衡定数 K_eq）、`catalyst_model.py`（触媒失活 a(T,t)） |
| **状態方程式 (EOS)** | `eos.py`（Peng-Robinson: Z 因子・フガシティー・残差熱量・泡点/露点・断熱圧縮） |
| **蒸留計算コア** | `distillation_core.py`（FUG ショートカット本体 + ディスパッチ）、`distillation_rigorous.py`（Wang-Henke 厳密法）、`distillation_sm.py`（GPR サロゲート） |
| **コスト・経済パラメータ** | `cost_parameters.py`（Bare Module 係数・ユーティリティ単価・原料/製品単価・触媒/吸着剤/膜パラメータ・Hasebe OPEX 係数）、`cost_calculator.py`（容器/HE/圧縮機/ポンプ/トレイ/膜の CAPEX 推算）、`utility_selector.py`（温度→ユーティリティ tier 選択と単価） |

---

## ファイル一覧

| ファイル | 役割（1 行説明） |
|---|---|
| `config.py` | 物理定数（R, T0）、PDH 反応速度パラメータ（3 反応のデータクラス）、成分熱力学データ辞書 `THERMO_DATA`（生成エンタルピー/Gibbs、Cp 多項式係数、PR EOS の Tc/Pc/ω）。 |
| `component_data.py` | 成分の分子量・定圧比熱（範囲固定値）・標準沸点・蒸発潜熱・液体密度のテーブルと加重平均ヘルパ（`cp_mix`, `liquid_density_mix`）。 |
| `kinetics.py` | `PDHKinetics` クラス。脱水素 r1（可逆・吸着抑制）/クラッキング r2/水素化 r3 の速度式を修正 Arrhenius で計算。 |
| `thermo.py` | `PDHThermo` クラス。Cp(T) 多項式、エンタルピー変化の解析積分、Kirchhoff+Gibbs-Helmholtz による反応1平衡定数 `K_eq`。 |
| `catalyst_model.py` | 触媒活性 `a(T_℃, t_min)` を `data/a_parameter_fitting.xlsx` 由来のグリッドから 2 次元スプライン補間（外挿はクリップ）。 |
| `eos.py` | Peng-Robinson EOS。Z 因子、フガシティー係数、残差エンタルピー/エントロピー、泡点/露点温度（泡点は `thermo` パッケージの PRMIX を内部利用）、断熱圧縮。 |
| `distillation_core.py` | 蒸留塔エントリポイント `simulate_distillation_column`。FUG（Fenske-Underwood-Gilliland）ショートカット本体 + `solver_method` による rigorous/sm へのディスパッチ、塔径・塔高・CAPEX・ユーティリティ・proxy 罰則を含む `DistResult` 組立。 |
| `distillation_rigorous.py` | `wang_henke_solve`。MESH 方程式 tray-by-tray の Wang-Henke 法（bubble-point method）。partial/total condenser 対応、always-on MESH/成分収支検証付き。 |
| `distillation_sm.py` | 学習済み GPR サロゲート（`models/column{1,3}_sm.pkl`）で Dist1/Dist3 を HYSYS の代替として高速予測。`solve_column1_via_sm` / `solve_column3_via_sm`。 |
| `cost_parameters.py` | Bare Module Cost 法の係数（容器/HE/圧縮機/ポンプ/トレイ/膜）、CEPCI、為替、ユーティリティ単価（蒸気・冷媒・電力・燃料）、原料/製品単価、触媒/吸着剤/膜パラメータ、Hasebe OPEX 係数、HHV。 |
| `cost_calculator.py` | `cost_parameters` の係数を使った CAPEX 計算関数群（縦型容器・HE・圧縮機・ポンプ・Sieve トレイ・膜モジュール、すべて億円返し）。 |
| `utility_selector.py` | ターゲット温度から冷却/加熱ユーティリティ tier を選択し単価を返す。離散モード（階段関数）と連続モード（区分線形補間、BO 向け）。 |

---

## 依存関係の概観

レイヤ内の依存は概ね一方向（下から上）：

```
config.py  ─┬─→ thermo.py ─┐
            ├─→ kinetics.py │
            └─→ eos.py ─────┤
                            ▼
component_data.py ──────────┤
cost_parameters.py ─→ cost_calculator.py
cost_parameters.py ─→ utility_selector.py
                            ▼
            distillation_core.py  (FUG 本体 + ディスパッチ)
              ├─ import → eos.py / thermo は不使用 / component_data.py
              ├─ import → cost_calculator.py / cost_parameters.py / utility_selector.py
              ├─ import → flowsheet.heat_integration (U 値・相判定)
              ├─ import → stream.stream.ProcessStream
              ├─ 'rigorous' dispatch → distillation_rigorous.wang_henke_solve
              └─ 'sm'       dispatch → distillation_sm (column{1,3} ラッパ経由が標準)

distillation_rigorous.py ─→ eos.py (bubble_point_T, fugacity_coeff, z_factor)
distillation_sm.py       ─→ distillation_core (型) / units.vle.hysys.provider (結果組立)
catalyst_model.py        ─→ (scipy のみ。reactor 側から呼ばれる)
```

外部パッケージ依存:
- `numpy`, `scipy`（補間・最適化・線形代数）
- `thermo`（CalebBell/thermo, MIT, v0.6.0 ピン留め）— `eos.py` の `bubble_point_T` が PRMIX を内部利用。未 import 時は `eos.py` ロードで fail-fast。

被参照（このレイヤを使う側、抜粋）:
`units/reactors/swing.py`・`units/reactors/radial_flow.py`（kinetics/thermo/catalyst_model/cost）、
`units/separators/column{1,2,3}/`（distillation_*）、`units/separators/psa`・`membrane`（eos/cost）、
`units/utils/`（component_data/cost/utility_selector）、
`flowsheet/economics.py`・`flowsheet/run_one_pass.py`（cost_parameters/utility_selector）。

---

## 詳細仕様（SPEC）

複雑なモジュールは別途 SPEC を用意している。理論式・入出力・既知の限界はそちらを参照:

- [`SPEC_distillation.md`](SPEC_distillation.md) — 蒸留 3 経路（FUG / Wang-Henke 厳密 / GPR サロゲート）の理論・式・入出力・使い分け
- [`SPEC_eos.md`](SPEC_eos.md) — Peng-Robinson EOS と泡点/露点、`thermo` 連携

「仮置き」「確認中」の値は各 SPEC 末尾に列挙し、詳細は
[`../KNOWN_PLACEHOLDERS.md`](../KNOWN_PLACEHOLDERS.md) に集約している。

成分マッピング（全 `src/` 共通）:

| キー | 化学式 | 名称 |
|---|---|---|
| A | C₃H₈ | プロパン |
| B | C₃H₆ | プロピレン |
| C | H₂ | 水素 |
| D | C₂H₄ | エチレン |
| E | CH₄ | メタン |
| F | C₂H₆ | エタン |
| Z | C₄H₁₀ | n-ブタン |
