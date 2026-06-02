# SPEC: radial_flow.py — PDH 径方向流スイング反応器システム シミュレーター

**ファイルパス**: `units/reactors/radial_flow.py`
**位置づけ初版**: 2026-05-30 / **多段化 (Oleflex 型)**: 2026-05-31

---

## 0. 本 SPEC の読み方

本ファイルは軸流 `swing.py` の幾何を「軸 z → 半径 r」に差し替えた**姉妹実装**である。
反応速度・熱力学・Ergun 物性・触媒失活・コスト計算・出力データクラスは **すべて
`swing.py` から import して共有**しており、独自に持つのは「幾何 (径方向)」と「径方向 ODE」
「多段直列 + 段間再加熱」だけ。したがって本 SPEC は径方向流 **固有の差分**を中心に記述し、
共通部分は `SPEC_swing.md` を参照する。

| 共有元 (`swing.py` から import) | 内容 |
|---|---|
| `calc_a`, `calc_Cp`, `calc_rate_constants`, `_reaction_enthalpies` | 失活・比熱・速度定数・反応エンタルピー |
| `_gas_viscosity`, `_STOICH`, `_COMPS`, `_R_GAS`, `_P_FLOOR_PA` | Ergun 物性・量論行列・成分順序・数値ガード |
| `FeedStream`, `FixedParams` | 入口流体条件・固定パラメータ (径方向流も同一型を使用) |
| `EffluentStream`, `EquipmentCost`, `PerformanceMetrics`, `SimulationResult` | 出力データクラス (下流 run_one_pass/economics はそのまま) |
| `_penalty_result`, `_PENALTY_CAPEX`, `calc_reactor_capex_okuyen` | ペナルティ・CAPEX 推算 |

---

## 1. 目的

### なぜ径方向流か

反応器設計レビューで「0.5 bar 低圧固定床では圧力損失が支配的になり、軸流深床
(`swing.py`) は現実的粒径 (3 mm) でほぼ全域 infeasible になる」と判明した
(詳細は `swing.py` ヘッダ + `monitor/reactor_pressure_drop_and_geometry.ipynb`)。

その物理的解が **径方向流 (radial flow) ベッド**:

- ガスを薄い環状床に「半径方向に短く」通す → 流路長 = 床厚 Δr (~0.3–1 m) と短い
- 流路断面積 = 2π r H (塔高 H を稼げば大きい) → 空塔速度を低く抑えられる
- 触媒量 (= 塔高 H と床厚 Δr で決まる) と圧損 (= 床厚 Δr のみで決まる) を **分離**できる

⇒ 体積 vs 圧損のトレードオフが解消し、3 mm 粒でも 0.5 bar で ΔP が収まる。これは実機
UOP Oleflex (径方向流) / Lummus Catofin (減圧・浅床) の設計思想そのもの。

### 設計方針

`swing.py` (軸流) は温存し、本ファイルを **新規追加**で並走させる。化学 (反応・熱・失活)
は軸流と同一なので、同じ V_cat・同じ feed なら転化率/選択率/温度プロファイルは軸流と
ほぼ一致し、**差が出るのは ΔP だけ** (= 径方向流が「圧損だけ」を救う、という主張を担保)。

---

## 2. 反応系・速度・熱力学・失活モデル

`swing.py` と完全に同一 (import 共有)。詳細は [SPEC_swing.md](SPEC_swing.md) の
§2 (反応系)・§4 (速度・熱力学モデル)・§5 (触媒失活モデル) を参照。

- 成分: A=C₃H₈, B=C₃H₆, C=H₂, D=C₂H₄, E=CH₄, F=C₂H₆
- 反応: r1 (脱水素・可逆), r2 (クラッキング), r3 (水素化)
- 量論行列 `_STOICH`、速度式 r1/r2/r3、修正 Arrhenius、Kirchhoff/Gibbs-Helmholtz K_eq

---

## 3. 物理モデルと計算手法 (径方向流 固有)

### 3-1. 幾何

環状床を「外周 r_o から中心 r_i へ」流す centripetal モデル。

| 記号 | 意味 | 関係 |
|---|---|---|
| $r_i$ | 環状床の内半径 (中心捕集管の外半径) | $= D_{inner}/2$ |
| $r_o$ | 環状床の外半径 | $= r_i + \Delta r$ |
| $\Delta r$ | 床厚 (= 流路長、圧損を決める) | `bed_thickness` |
| $H$ | 環状床の高さ | `H` |
| $A(r)$ | 局所断面積 (円筒側面積) | $= 2\pi r H$ (r 依存) |

軸流では $A_{cross} = \pi D^2/4$ が一定なのに対し、径方向流は $A(r)=2\pi r H$ が r に依存し、
内側ほど面積が小さく速度が上がる点が本質的な違い。

### 3-2. 半径方向 ODE (`_ode_radial`)

状態ベクトル: `y = [F_A, F_B, F_C, F_D, F_E, F_F, T, P]`(単位: mol/s, K, Pa)。
独立変数 $r$ を $r_o \to r_i$ へ積分 (r 減少方向 = 流れ方向)。

**物質収支**

$$\frac{dF_i}{dr} = -(1-\varepsilon)\cdot 2\pi r H \cdot \sum_j \nu_{ij}\, r_j$$

**エネルギー収支 (断熱)**

$$\frac{dT}{dr} = +\frac{(1-\varepsilon)\cdot 2\pi r H \cdot \sum_j \Delta H_{rxn,j}(T)\, r_j}{\sum_i F_i\, C_{p,i}(T)}$$

**圧力損失 (Ergun、流れ方向に圧損)**

$$\frac{dP}{dr} = +\left(150\frac{(1-\varepsilon_b)^2\,\mu\, u}{\varepsilon_b^3\,(\phi d_p)^2}
+ 1.75\frac{(1-\varepsilon_b)\,\rho\, u^2}{\varepsilon_b^3\,\phi d_p}\right)$$

ここで局所空塔速度 $u = Q_{vol}/A(r) = Q_{vol}/(2\pi r H)$ (r 内側ほど面積↓→速度↑)。
符号は「r 減少方向 = 流れ方向」で反応物消費・降温・降圧が起こるよう取る (dF/dr<0,
dT/dr>0 で r 減少につれ降温, dP/dr>0 で r 減少につれ降圧)。

| 記号 | 意味 | 単位 |
|---|---|---|
| $\varepsilon$ | 床/容器比 (`FixedParams.eps`、軸流と同義) | — |
| $\varepsilon_b$ | 粒子間空隙率 (`FixedParams.eps_bed`) | — |
| $d_p$ | 触媒粒径 (`FixedParams.d_p_m`) | m |
| $\phi$ | 形状係数 (`FixedParams.sphericity`) | — |
| $u$ | 局所空塔速度 $= Q_{vol}/(2\pi r H)$ | m/s |
| $\rho$ | 理想気体 + 組成平均 MW のガス密度 | kg/m³ |
| $\mu$ | `_gas_viscosity(T)` による温度依存近似 (!仮置き) | Pa·s |

**ソルバー設定**: Radau、`rtol=1e-4`, `atol=1e-7`(軸流 `swing.py` と同値)。

### 3-3. 並列分割の規約 (案B, 2026-05-31 確定)

設計変数 (`D_inner`, `bed_thickness`, `H`) は「作る反応器全体」を表す。触媒体積 $V_{cat}$
が $V_{cat,max}$ (200 m³) を超える場合、総断面 $2\pi r H$ を $N_{parallel}$ 本のサブカラムに
**機械分割**する。

> 重要: 総流量 ÷ 総断面 = 各サブカラムの空塔速度は分割で不変なので、**u・ΔP・SV は
> $N_{parallel}$ に依存しない**。`_ode_radial` 内で u を $N_{parallel}$ で割らないのはこのため
> (旧実装は割っており ΔP・SV を過小評価するバグだった)。$N_{parallel}$ は容器分割・コスト
> 計上 ($V_{vessel}=V_{design}/N_{parallel}$、基数 $=N_{parallel}\cdot N_{swing}$) にのみ使う。

### 3-4. 時間方向積分

`swing.py` と同じ。t=0…t_cyc を `n_time_samples`(既定 20) 点で台形則平均。各時刻 t で
$a(T_{in},t)$ を 1 回評価し r 方向に定数として使用。

---

## 4. サイジング・制約・コスト

### 4-1. 触媒体積・並列基数

$$V_{cat,total} = \pi (r_o^2 - r_i^2)\, H\, (1-\varepsilon)$$

$$N_{parallel} = \max\!\left(\left\lceil V_{cat,total}/V_{cat,max}\right\rceil,\ 1\right)$$

### 4-2. 圧力損失ハード制約 (床 Ergun × 内部品マージン)

床 ODE で得た真の床 ΔP に `dP_margin_factor`(既定 1.4、!仮置き) を掛けて総反応器 ΔP とし、
出口圧・feasibility に適用する (床内の化学計算には混入させない)。

$$\frac{\Delta P}{P_{in}} = \mathrm{clip}\!\left(\frac{(P_{in}-P_{out,床})}{P_{in}}\times f_{margin},\ 0,\ 1\right)$$

$P_{out} = P_{in}(1-\Delta P/P_{in})$ を下流へ伝播。$\Delta P/P_{in} >$ `dP_over_P_max`(既定 0.10)
で infeasible (`reason='dP_excess'`)。

### 4-3. 空塔速度チェック (SV_max のみ、内側面)

軸流の [SV_min, SV_max] のうち **SV_max のみ**を、速度最大の内側面 $r_i$ について課す:

$$SV_{inner} = \frac{Q_{vol,in}}{2\pi r_i H}$$

$SV_{inner} > $ `SV_max_m_per_s` で infeasible (`reason='sv_out_of_range'`)。
**SV_min は課さない**: 径方向流は塔高 H を稼ぐと空塔速度が下がるのが本来の利点で、
低速はディストリビュータ設計で扱う想定 (軸流の channeling 懸念とは事情が異なる)。
SV も案B により $N_{parallel}$ で割らない。

### 4-4. スイング・容器・触媒量・CAPEX

$$N_{swing,sets} = \lceil t_{regen}/t_{cyc}\rceil + 1,\quad
N_{reactors,total} = N_{parallel}\times N_{swing,sets}$$

容器体積は中心捕集管 (半径 $r_i$) の void を除外し、触媒層を含む**環状部のみ**計上
(ユーザー決定 2026-05-30)。この定義では $V_{vessel}=V_{cat}/(1-\varepsilon)$ となり軸流と整合:

$$V_{vessel,actual} = \frac{\pi (r_o^2 - r_i^2)\, H}{N_{parallel}}$$

$$W_{cat,total} = V_{cat,total}\times N_{swing,sets}\times \rho_b$$

CAPEX は縦型プロセス容器として `calc_reactor_capex_okuyen` で推算 (外径 $= 2 r_o$)。
式・係数は [SPEC_swing.md §7](SPEC_swing.md) と同一。

---

## 5. 多段反応器 (Oleflex 型: `simulate_radial_multibed_reactor_system`)

### 5-1. なぜ多段か

単段断熱の径方向流は、PDH の強い吸熱で床が降温し、出口温度での平衡に張り付いて
per-pass 転化率が ~30% で頭打ちになる (触媒量を増やしても 30% で飽和)。実機 PDH は
これを「単段で解かない」: UOP Oleflex は径方向流断熱反応器を 3〜4 基直列に並べ、各反応器
の間に加熱炉 (interstage heater) を置き、吸熱で冷えたガスを反応温度まで戻してから次段へ
送る。各段で平衡がリセットされ累積転化が上がる (実測 3 段で ~59%)。本関数はこの Oleflex
型を、既存の単段ソルバを N 回直列に呼ぶだけで表現する。

### 5-2. 時刻同期で直列計算 (2026-05-31)

スイング列は時刻 t ごとに Bed1(t)→reheat→Bed2(t)→reheat→Bed3(t) と流れる。各段の活性
$a(t)$ は同一 t (同時に新鮮再生) で決まる。反応は非線形なので「各段で時間平均してから次段
へ渡す」と $f(\text{平均})\neq\text{平均}(f)$ のズレ (特に選択率) が出る。よって本実装は
**時刻サンプル t ごとに全段を直列積分し、最終段出口を最後に時間平均する**。段間再加熱の
熱量も t ごとに評価して平均する。

### 5-3. 段間再加熱の燃料費

各段の積分は `design.T_in` から開始する (= reheat 済)。各段 t の昇温熱量「前段 per-t 出口
温度 → T_in」を合算・時間平均して `Q_preheat` に計上 → economics の「Reactor 予熱炉燃料」
に総量が乗る。

### 5-4. 装置の積算

全段同一ジオメトリ (`RadialDesignVars` の D_inner/bed_thickness/H) とし、触媒量・CAPEX・
総基数を n_beds 段ぶん積算する ($N_{reactors,total} = N_{parallel}\cdot N_{swing,sets}\cdot n\_beds$)。
`n_beds=1` のときは単段ソルバ `simulate_radial_flow_reactor_system` に委譲し完全一致。

---

## 6. 入出力

### 入力

#### `RadialDesignVars` — 設計変数 (径方向流 固有)

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `T_in` | float | K | 反応器入口温度 (各段は reheat でここまで戻す) |
| `t_cyc` | float | min | 1サイクル反応フェーズ時間 |
| `D_inner` | float | m | 環状床の内径 (中心捕集管の外径)。$r_i = D_{inner}/2$ |
| `bed_thickness` | float | m | 環状床の厚さ $\Delta r$。$r_o = r_i + \Delta r$ |
| `H` | float | m | 環状床の高さ |

プロパティ `r_i`(=D_inner/2), `r_o`(=r_i+bed_thickness) を提供。
軸流 `DesignVars` の (T_in, z_cat, t_cyc, D) に対応する幾何違いの設計変数。

#### `FeedStream` / `FixedParams`

`swing.py` と同一型を import して使用。[SPEC_swing.md §9](SPEC_swing.md) 参照
(`dP_margin_factor` 含む全パラメータが共通)。

### 出力

`SimulationResult`(= `swing.py` と同一型)。`EffluentStream`/`EquipmentCost`/
`PerformanceMetrics` の各フィールドも軸流と同じ。多段の場合 `EquipmentCost` は全段合計
(N_reactors_total・Catalyst_Weight_Total・Reactor_CAPEX)、`V_vessel_actual` は 1 基代表値。

### 関数シグネチャ

```python
def simulate_radial_flow_reactor_system(
    design: RadialDesignVars, feed: FeedStream, fixed: FixedParams,
    n_time_samples: int = 20,
) -> SimulationResult: ...

def simulate_radial_multibed_reactor_system(
    design: RadialDesignVars, feed: FeedStream, fixed: FixedParams,
    n_beds: int = 3, n_time_samples: int = 20,
) -> SimulationResult: ...
```

---

## 7. 主要パラメータ

`FixedParams`(swing と共有) の Ergun/SV/サイジング系がそのまま効く。径方向流固有で
特に効くのは:

| パラメータ | 既定 | 効き方 |
|---|---|---|
| `bed_thickness` (設計変数) | — | 流路長 = ΔP を支配。小さいほど圧損↓ |
| `H` (設計変数) | — | 断面積 $2\pi r H$ → SV を支配。大きいほど SV↓ |
| `dP_margin_factor` | 1.4 (!仮置き) | 床以外 (分配器/中心管/弁/段間配管) 圧損の一括マージン |
| `d_p_m` | 0.003 (!仮置き) | Ergun 粘性項に $1/(\phi d_p)^2$ で強く効く |
| `n_beds` (引数) | 3 | 多段化による累積転化率の引き上げ |

---

## 8. 出典

- Ergun 式: 標準的固定床圧損相関 (Ergun 1952)。係数 150 / 1.75 は文献標準値。
- 多段 + 段間 reheat 構成: UOP Oleflex / Lummus Catofin の実機設計思想。
- 反応速度・触媒・コスト出典: すべて [SPEC_swing.md](SPEC_swing.md) を参照。

---

## 9. 既知の限界・仮置き

`swing.py` と共通のモデル限界 (コード冒頭ヘッダに明記):

| # | 限界 | 備考 |
|---|---|---|
| 1 | 再生動特性なし | 各サイクル開始時に新鮮触媒 $a=1.0$ で初期化 |
| 2 | 粒内拡散 (Weisz-Prater) 未考慮 | 反応速度は粒外バルク濃度ベース |
| 3 | 触媒活性は入口温度で代表し空間一定 | $a(T_{in},t)$ を r 方向に定数使用 |
| 4 | 床 Ergun 以外の ΔP は係数一括 | 分配器/中心管/弁/段間配管を `dP_margin_factor` で粗く見込む。確定値は【確認中】 |
| 5 | 段間加熱炉 CAPEX は独立計上せず | 「反応器バンドルに内包」の現行規約。燃料 OPEX は Q_preheat 合算で計上。独立計上は今後の精緻化候補 |
| 6 | 段間加熱炉 ΔP | 圧縮機ループで回復する別項として `dP_margin_factor` に含めない |

`!仮置き` (`_gas_viscosity` 近似、`d_p_m`、`eps_bed`、`sphericity`、`dP_margin_factor`) /
`【確認中】` の詳細は **`KNOWN_PLACEHOLDERS.md`** を参照。

---

## 10. エラーハンドリング

`swing.py` と同方針 (無効入力・数値異常はすべて `_penalty_result()` で吸収しクラッシュ
しない)。`reason` ラベル: `'input_invalid'` / `'sim_failure'` / `'dP_excess'` /
`'sv_out_of_range'` / `'volume_zero'`。多段ではどこか 1 段・1 時刻でも積分が失敗したら
`sim_failure` を伝播。詳細は [SPEC_swing.md §13](SPEC_swing.md) を参照。
