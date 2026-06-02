# SPEC: membrane_system.py — PDH 膜分離システム シミュレーター

**ファイルパス**: `units/separators/membrane/membrane_system.py`  
**依存 EOS モジュール**: `src/eos.py`  
**最終更新**: 現行コードへ整合 (多段圧縮+段間冷却・vapor_feed・膜性能感度係数・P_L 1atm 固定方針を反映)

---

## 目次

1. [シミュレータ概要](#1-シミュレータ概要)
2. [成分系と分離対象](#2-成分系と分離対象)
3. [物理モデルと計算手法](#3-物理モデルと計算手法)
4. [Peng-Robinson EOS モデル (src/eos.py)](#4-peng-robinson-eos-モデル-srceospyf)
5. [気化器モデル](#5-気化器モデル)
6. [圧縮機モデル](#6-圧縮機モデル)
7. [膜分離モデル（クロスフロー ODE）](#7-膜分離モデルクロスフロー-ode)
8. [製品冷却器モデル](#8-製品冷却器モデル)
9. [OPEX/CAPEX 出力](#9-opexcapex-出力)
10. [仮定一覧](#10-仮定一覧)
11. [データクラス定義](#11-データクラス定義)
12. [定数・パラメータ一覧](#12-定数パラメータ一覧)
13. [依存モジュール](#13-依存モジュール)
14. [使い方](#14-使い方)
15. [エラーハンドリング](#15-エラーハンドリング)

---

## 1. シミュレータ概要

### 何をシミュレートするか

PDH プロセスにおける **C₃H₆/C₃H₈ 膜分離システム** を 5 ユニット一貫でシミュレートする。  
前段蒸留塔塔底液（C₃H₆/C₃H₈ 混合液）を受け取り、膜を用いて C₃H₆ を濃縮し、  
後段精製蒸留塔へ飽和液フィードとして送り出す。

| 項目 | 内容 |
|---|---|
| 分離対象 | C₃H₆ (プロピレン) / C₃H₈ (プロパン) 二成分系 |
| 分離方式 | ZIF-8 系ポリマー膜（クロスフロー、圧力駆動型） |
| 熱力学モデル | Peng-Robinson EOS（実在気体補正） |
| 設計変数 | 膜供給側圧力 $P_H$、透過側圧力 $P_L$、総膜面積 $A_{mem}$、後段蒸留塔操作圧力 $P_{dist}$ |
| 出力 | 後段蒸留塔向け飽和液ストリーム＋機器サイズ（OPEX/CAPEX 推算テーブル） |

### プロセスフロー

```
[前段蒸留塔 塔底液]
       │ MemFeedStream (F_C3H6, F_C3H8, T_in, P_in)
       ▼
  ┌──────────┐     全量ガス化 (PR 露点計算 + エンタルピー差分)
  │ 気化器    │──▶ T_vap_out, Q_vap_kW, A_vap
  └──────────┘
       │ ガス @ P_in
       ▼
  ┌──────────────┐  多段等エントロピー圧縮 + 段間冷却 (PR 補正)
  │フィード圧縮機 │──▶ T_feed_comp_out, W_feed_kW, Q/A_intercool
  └──────────────┘
       │ ガス @ P_H
       ▼
  ┌────────────────┐  クロスフロー ODE (solve_ivp Radau)
  │ 膜モジュール   │──▶ 非透過 (Retentate) / 透過 (Permeate)
  └────────────────┘
       │ Permeate @ P_L
       ▼
  ┌────────────────┐  多段等エントロピー圧縮 + 段間冷却 (PR 補正)
  │ 製品圧縮機     │──▶ T_prod_comp_out, W_prod_kW, Q/A_intercool
  └────────────────┘
       │ ガス @ P_dist
       ▼
  ┌──────────────┐  凝縮・液化 (PR 泡点計算 + エンタルピー差分)
  │ 製品冷却器   │──▶ T_bp_perm, Q_cond_kW, A_cond
  └──────────────┘
       │ MemProductStream (飽和液 @ P_dist)
       ▼
[後段蒸留塔 フィード]
```

---

## 2. 成分系と分離対象

### 成分定義

膜分離モジュールは C₃H₆ と C₃H₈ の **2 成分系** を扱う。  
EOS および膜計算内での成分インデックスは以下のとおり固定。

| インデックス | config.py キー | 化学式 | 名称 | 役割 |
|---|---|---|---|---|
| 0 | `'B'` | C₃H₆ | プロピレン | 膜透過優先成分（製品） |
| 1 | `'A'` | C₃H₈ | プロパン | 非透過成分（未反応原料） |

```python
_KEYS = ['B', 'A']  # index 0 = C3H6, index 1 = C3H8
```

### 前段・後段との接続

| 接続先 | ストリーム | 状態 |
|---|---|---|
| 前段蒸留塔 塔底 → 気化器 | `MemFeedStream` | 飽和液または過冷却液 |
| 製品冷却器 → 後段蒸留塔 フィード | `MemProductStream` | 飽和液 @ `P_dist` |
| 非透過ガス 出口 | `MemRetentateStream` | ガス @ `P_H`（リサイクルまたは排出） |

---

## 3. 物理モデルと計算手法

### 3-1. 全体フロー

```
入力 (MemDesignVars, MemFeedStream, MemFixedParams)
│
├─ バリデーション: P_H > P_L, A_mem > 0, etc.
│
├─ ユニット 1: 気化器
│     dew_point_T → T_dew
│     T_vap_out = T_dew + T_vap_superheat
│     _h_mol(liquid) → _h_mol(vapor) → Q_vap, A_vap
│
├─ ユニット 2: フィード圧縮機
│     compress_isentropic(T_vap_out, P_in, P_H) → T_feed_comp_out, W_feed_per_mol
│
├─ ユニット 3: 膜モジュール (クロスフロー ODE)
│     solve_ivp Radau: F_C3H6(A), F_C3H8(A) → F_ret_C3H6, F_ret_C3H8
│     → y_C3H6, x_ret_C3H6, stage_cut
│
├─ ユニット 4: 製品圧縮機
│     compress_isentropic(T_feed_comp_out, P_L, P_dist) → T_prod_comp_out, W_prod_per_mol
│
├─ ユニット 5: 製品冷却器
│     bubble_point_T → T_bp
│     _h_mol(vapor) → _h_mol(liquid) → Q_cond, A_cond
│
└─ 出力 (MemSimulationResult)
```

### 3-2. エンタルピー計算の統一基準

全ユニットで共通のモルエンタルピー $h$ を使用する。

$$h(T, P, z_{C3H6}, \text{phase}) = H^{ig}(T) + H^r(T, P, x, Z_{\text{phase}})$$

$$H^{ig}(T) = \sum_i x_i \int_{T_{ref}}^{T} C_{p,i}(T')\, dT'$$

基準温度 $T_{ref} = 298.15\,\text{K}$（理想気体状態）

| 記号 | 意味 | 単位 |
|---|---|---|
| $H^{ig}(T)$ | 理想気体モルエンタルピー（$T_{ref}$ 基準） | J/mol |
| $H^r(T, P, x, Z)$ | 残差エンタルピー（PR EOS による実在気体補正） | J/mol |
| $C_{p,i}(T)$ | 成分 $i$ の定圧モル比熱多項式 | J K⁻¹ mol⁻¹ |
| $x_i$ | 成分 $i$ のモル分率 | — |

---

## 4. Peng-Robinson EOS モデル (src/eos.py)

### 4-1. PR 状態方程式

$$P = \frac{RT}{V_m - b_m} - \frac{a_m(T)}{V_m(V_m + b_m) + b_m(V_m - b_m)}$$

Z 因子形式（3 次方程式）:

$$Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0$$

$$A = \frac{a_m P}{R^2 T^2}, \quad B = \frac{b_m P}{RT}$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $Z$ | 圧縮率因子 $= PV_m/(RT)$ | — |
| $A, B$ | 無次元 PR パラメータ | — |
| $a_m(T)$ | 混合物引力パラメータ（温度依存） | Pa·m⁶ mol⁻² |
| $b_m$ | 混合物斥力パラメータ（温度不依存） | m³ mol⁻¹ |

### 4-2. 単成分パラメータ

$$a_i(T) = \left[\sqrt{a_{c,i}}\left(1 + \kappa_i\left(1 - \sqrt{T/T_{c,i}}\right)\right)\right]^2$$

$$a_{c,i} = 0.45724\, \frac{R^2 T_{c,i}^2}{P_{c,i}}, \quad b_i = 0.07780\, \frac{R T_{c,i}}{P_{c,i}}$$

$$\kappa_i = 0.37464 + 1.54226\,\omega_i - 0.26992\,\omega_i^2$$

$$\frac{d\sqrt{a_i}}{dT} = -\frac{\kappa_i \sqrt{a_{c,i}}}{2\sqrt{T\, T_{c,i}}}$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $a_{c,i}$ | 臨界点における $a_i$ | Pa·m⁶ mol⁻² |
| $\kappa_i$ | ソアベ補正係数 | — |
| $T_{c,i}$ | 臨界温度 | K |
| $P_{c,i}$ | 臨界圧力 | Pa |
| $\omega_i$ | アセントリック因子 | — |

**C₃H₆/C₃H₈ の臨界定数（NIST 値）**

| 成分 | $T_c$ [K] | $P_c$ [Pa] | $\omega$ |
|---|---|---|---|
| C₃H₆ (B) | 364.90 | 4.6000×10⁶ | 0.1408 |
| C₃H₈ (A) | 369.89 | 4.2512×10⁶ | 0.1521 |

### 4-3. 混合則（van der Waals、kij=0）

$$a_m = \left(\sum_i x_i \sqrt{a_i}\right)^2 = S^2, \quad S \equiv \sum_i x_i \sqrt{a_i}$$

$$b_m = \sum_i x_i b_i$$

$$\frac{da_m}{dT} = 2 S \dot{S}, \quad \dot{S} \equiv \sum_i x_i \frac{d\sqrt{a_i}}{dT}$$

### 4-4. フガシティー係数

$$\ln \varphi_i = \frac{b_i}{b_m}(Z-1) - \ln(Z-B) - \frac{A}{2\sqrt{2}\,B}\left(\frac{2\sqrt{a_i}}{S} - \frac{b_i}{b_m}\right)\ln L$$

$$L \equiv \frac{Z + (1+\sqrt{2})\,B}{Z + (1-\sqrt{2})\,B}$$

### 4-5. 残差エンタルピー

$$H^r = RT(Z-1) + \frac{T\,\dfrac{da_m}{dT} - a_m}{2\sqrt{2}\,b_m}\ln L$$

### 4-6. 残差エントロピー

$$S^r = R\ln(Z-B) + \frac{\dfrac{da_m}{dT}}{2\sqrt{2}\,b_m}\ln L$$

### 4-7. 泡点温度 `bubble_point_T`

圧力 $P$・液相組成 $x$ を固定し、$\sum_i x_i K_i = 1$ を満たす温度 $T_{bp}$ を求める。

**外側ループ（brentq）**: $T$ を変化させ、目的関数 $f(T)=\sum_i x_i K_i - 1$ のゼロ点を探索。  
**内側ループ（逐次置換、20 回）**: $\varphi^L$ を $x$ で固定し、$y \leftarrow x K / \sum x K$ を繰り返して $\varphi^V$ を更新。

単相領域（3 次方程式の実根が 1 本）では、Z 値の大小から符号センチネルを返して brentq の符号変化を保証する:

| 状態 | $Z_s < 0.5$ (液相様) | $Z_s \geq 0.5$ (気相様) |
|---|---|---|
| 泡点 `obj` の戻り値 | $-2.0$ | $+2.0$ |

**許容誤差**: `xtol=0.05 K`,  `maxiter=200`

### 4-8. 露点温度 `dew_point_T`

圧力 $P$・気相組成 $y$ を固定し、$\sum_i y_i / K_i = 1$ を満たす温度 $T_{dp}$ を求める。  
内側ループで $x$ を更新する点が泡点と対称であり、センチネルの符号は反転する:

| 状態 | $Z_s < 0.5$ (液相様) | $Z_s \geq 0.5$ (気相様) |
|---|---|---|
| 露点 `obj` の戻り値 | $+2.0$ | $-2.0$ |

### 4-9. 断熱圧縮 `compress_isentropic`

等エントロピー出口温度 $T_{2s}$ を求め、効率 $\eta$ で実仕事・実出口温度を計算する。

**手順:**

1. $\Delta S^{ig}(T_1 \to T_{2s}, P_1 \to P_2) + S^r(T_{2s},P_2) - S^r(T_1,P_1) = 0$ を brentq で $T_{2s}$ について解く
2. $W_{isen} = \Delta H^{ig}(T_1 \to T_{2s}) + H^r(T_{2s},P_2) - H^r(T_1,P_1)$
3. $W_{actual} = W_{isen} / \eta$
4. $\Delta H^{ig}(T_1 \to T_{2,act}) + H^r(T_{2,act},P_2) - H^r(T_1,P_1) = W_{actual}$ を brentq で $T_{2,act}$ について解く

$$\Delta S^{ig} = \sum_i x_i \left[ a_i \ln\frac{T_2}{T_1} + b_i(T_2 - T_1) + \frac{c_i}{2}(T_2^2 - T_1^2) + \frac{d_i}{3}(T_2^3 - T_1^3) \right] - R \ln\frac{P_2}{P_1}$$

$$\Delta H^{ig} = \sum_i x_i \left[ a_i(T_2 - T_1) + \frac{b_i}{2}(T_2^2-T_1^2) + \frac{c_i}{3}(T_2^3-T_1^3) + \frac{d_i}{4}(T_2^4-T_1^4) \right]$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $T_{2s}$ | 等エントロピー出口温度 | K |
| $T_{2,act}$ | 実出口温度 | K |
| $W_{isen}$ | 等エントロピー仕事 | J/mol |
| $W_{actual}$ | 実圧縮仕事（= 消費エネルギー） | J/mol |
| $\eta$ | 断熱効率（デフォルト 0.75） | — |

---

## 5. 気化器モデル

### 5-1. 概要

前段蒸留塔塔底液（C₃H₆/C₃H₈ 混合液）を受け取り、低圧蒸気で全量ガス化する。

### 5-2. 計算手順

1. **露点温度**: $T_{dew} = \text{dew\_point\_T}(P_{in}, y, \text{keys})$
2. **出口温度**: $T_{vap,out} = T_{dew} + \Delta T_{superheat}$（デフォルト +5 K）
3. **加熱量**:

$$Q_{vap} = F_{feed} \cdot \left[h(T_{vap,out}, P_{in}, z_{C3H6}, \text{vapor}) - h(T_{in}, P_{in}, z_{C3H6}, \text{liquid})\right]$$

4. **LMTD**（向流）:

$$\Delta T_1 = T_{hot} - T_{in}, \quad \Delta T_2 = T_{hot} - T_{vap,out}$$

$$\mathrm{LMTD} = \frac{\Delta T_1 - \Delta T_2}{\ln(\Delta T_1 / \Delta T_2)}$$

5. **伝熱面積**:

$$A_{vap} = \frac{Q_{vap}\,[\mathrm{kW}]}{U_{vap}\,[\mathrm{kW/(m^2 K)}] \times \mathrm{LMTD}}$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $T_{dew}$ | 供給組成・供給圧力における露点温度 | K |
| $T_{vap,out}$ | 気化器出口ガス温度 | K |
| $Q_{vap}$ | 気化器加熱量（= OPEX 蒸気消費基礎） | kW |
| $A_{vap}$ | 気化器伝熱面積（= CAPEX サイズ） | m² |
| $U_{vap}$ | 気化器総括伝熱係数（1.0, 化工便覧 改訂六版 表6・18） | kW/(m²·K) |
| $T_{hot}$ | 熱媒（LP Steam）温度（433.15 K = 160°C, コンテスト仕様） | K |

---

## 6. 圧縮機モデル（多段 + 段間冷却）

フィード圧縮機・製品圧縮機の両方で内部ヘルパ `_compress_multistage`（`compress_isentropic`
を段数分ループ + 段間冷却）を共通使用する。

### 6-0. 多段化ロジック `_compress_multistage`

圧縮比が大きいと 1 段では出口温度が過大になるため、`max_compression_ratio_per_stage`
(!仮置き 4.0) を超えないよう段数 $n$ を決め、各段で等比圧縮比 $(P_{out}/P_{in})^{1/n}$ で
昇圧する。各段の間で冷却水によりガスを `intercool_T_K`(!仮置き 40°C) まで顕熱冷却し、
次段入口温度とする。

- 返り値: `(T_out, W_total_per_mol, Q_intercool_kW, A_intercool_m2, n_stages, T_intercool_in_K)`
- `n=1`(圧縮比 ≤ r_max) のときは段間冷却なし → 従来の単段モデルと一致 (感度比較用)。
- 段間冷却はガス顕熱の与熱流体なので、温度域 ($T_{intercool,in}\to$ `intercool_T_K`) も返し、
  `flowsheet/heat_integration.extract_streams` で `'H_mem_intercool'` ホットストリームとして
  熱統合 (HI) に乗せる (二重計上なし)。段間冷却器面積 $A=Q/(U_{intercool}\cdot \mathrm{LMTD})$。

> 段間冷却の値はすべて !仮置き (`max_compression_ratio_per_stage` / `intercool_T_K` /
> `U_intercool`)。詳細は `KNOWN_PLACEHOLDERS.md`。

### フィード圧縮機

$$W_{feed} = F_{feed}\,[\mathrm{mol/s}] \times W_{actual}\,[\mathrm{J/mol}] \times 10^{-3} \quad [\mathrm{kW}]$$

| 入力 | 値 |
|---|---|
| 入口温度 | $T_{vap,out}$（気化器出口） |
| 入口圧力 | $P_{in}$（前段蒸留塔底圧力） |
| 出口圧力 | $P_H$（設計変数） |
| 流量 | $F_{feed}$ [mol/s] |

### 製品圧縮機

$$W_{prod} = F_{perm}\,[\mathrm{mol/s}] \times W_{actual}\,[\mathrm{J/mol}] \times 10^{-3} \quad [\mathrm{kW}]$$

| 入力 | 値 |
|---|---|
| 入口温度 | $T_{feed,comp,out}$（膜等温仮定）|
| 入口圧力 | $P_L$（設計変数） |
| 出口圧力 | $P_{dist}$（後段蒸留塔操作圧力） |
| 流量 | $F_{perm}$ [mol/s] |

---

## 7. 膜分離モデル（クロスフロー ODE）

### 7-1. 基本仮定

- **クロスフローモデル**: 非透過側はプラグフロー、透過ガスは即時排出（混合なし）
- **膜全域で温度・圧力一定**: $P_{feed} = P_H$, $P_{perm} = P_L$（等温・等圧）
- **透過係数は組成・圧力非依存**（線形透過則）

### 7-2. 透過フラックス

$$J_{C3H6} = Q_A \cdot (x \cdot P_H - y_{local} \cdot P_L)$$

$$J_{C3H8} = Q_B \cdot \left[(1-x) \cdot P_H - (1-y_{local}) \cdot P_L\right]$$

$$Q_B = Q_A / \alpha$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $Q_A$ | C₃H₆ 実効透過度 $= Q_{A,GPU}\times Q_{A,factor}$（既定 40 GPU × 1.0 → $1.340 \times 10^{-8}$ mol/(m²·s·Pa)） | mol m⁻² s⁻¹ Pa⁻¹ |
| $Q_B$ | C₃H₈ 透過度（$= Q_A / \alpha$） | mol m⁻² s⁻¹ Pa⁻¹ |
| $\alpha$ | C₃H₆/C₃H₈ 実効膜選択性 $= \alpha\times \alpha_{factor}$（既定 90 × 1.0） | — |
| $x$ | 非透過側局所 C₃H₆ モル分率 | — |
| $y_{local}$ | 透過側局所 C₃H₆ モル分率 | — |
| $\gamma$ | 圧力比 $= P_L / P_H$ | — |

**単位換算**: $1\,\mathrm{GPU} = 3.35 \times 10^{-10}\,\mathrm{mol/(m^2 \cdot s \cdot Pa)}$

### 7-3. 局所透過組成 `_y_local`

透過側と非透過側の分圧平衡条件から次の 2 次方程式を導く。

$$(1-\alpha)\gamma\, y^2 + \left[(\alpha-1)(x+\gamma)+1\right] y - \alpha x = 0$$

物理根（$0 \le y \le 1$）の選択には **負判別根**（$y < 0$ を回避する側）を採用:

$$y_{local} = \frac{2c}{-b - \sqrt{b^2 - 4ac}}$$

$$a = (1-\alpha)\gamma, \quad b = (\alpha-1)(x+\gamma)+1, \quad c = -\alpha x$$

### 7-4. クロスフロー ODE

$$\frac{dF_{C3H6}}{dA} = -J_{C3H6}, \quad \frac{dF_{C3H8}}{dA} = -J_{C3H8}$$

**初期条件**: $A = 0$ における $F_{C3H6}$, $F_{C3H8}$ は気化器・フィード圧縮機通過後の流量 [mol/s]  
**積分範囲**: $A = 0$ から $A = A_{mem}$ まで（$A_{mem}$ は設計変数）  
**数値ソルバー**: `scipy.integrate.solve_ivp`（Radau 法）、`rtol=1e-5`, `atol=1e-8`

フラックスがゼロ以下（分圧差の逆転）になった場合は `[0, 0]` を返して積分を停止する。

### 7-5. 出力量

| 量 | 計算式 | 単位 |
|---|---|---|
| Stage cut $\theta$ | $(F_{perm,C3H6} + F_{perm,C3H8}) / F_{feed}$ | — |
| 透過純度 $y_{C3H6}$ | $F_{perm,C3H6} / (F_{perm,C3H6} + F_{perm,C3H8})$ | — |
| 非透過 C₃H₆ 分率 $x_{ret}$ | $F_{ret,C3H6} / (F_{ret,C3H6} + F_{ret,C3H8})$ | — |
| モジュール本数 $n$ | $\lceil A_{mem} / A_{per,module} \rceil$ | — |

---

## 8. 製品冷却器モデル

透過ガスを製品圧縮機出口温度から飽和液まで凝縮させる。

### 8-1. 計算手順

1. **泡点温度**: $T_{bp} = \text{bubble\_point\_T}(P_{dist}, y, \text{keys})$
2. **冷却量**:

$$Q_{cond} = F_{perm} \cdot \left[h(T_{prod,comp,out}, P_{dist}, y_{C3H6}, \text{vapor}) - h(T_{bp}, P_{dist}, y_{C3H6}, \text{liquid})\right]$$

3. **LMTD**（向流熱交換器）:

$$\Delta T_1 = T_{prod,comp,out} - T_{cold,out}, \quad \Delta T_2 = T_{bp} - T_{cold,in}$$

4. **伝熱面積**:

$$A_{cond} = \frac{Q_{cond}\,[\mathrm{kW}]}{U_{cond}\,[\mathrm{kW/(m^2 K)}] \times \mathrm{LMTD}}$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $T_{bp}$ | 透過ガス組成・$P_{dist}$ における泡点温度 | K |
| $Q_{cond}$ | 冷却量（= OPEX 冷却水消費基礎） | kW |
| $A_{cond}$ | 冷却器伝熱面積（= CAPEX サイズ） | m² |
| $U_{cond}$ | 冷却器総括伝熱係数（1.0, 化工便覧 改訂六版 表6・18） | kW/(m²·K) |
| $T_{cold,in}$ | 冷却水入口温度（303.15 K = 30°C, コンテスト仕様） | K |
| $T_{cold,out}$ | 冷却水出口温度（313.15 K = 40°C, コンテスト仕様） | K |

---

## 9. OPEX/CAPEX 出力

### 9-1. OPEX 用データ

| 機器 | OPEX 基礎データ | 意味 |
|---|---|---|
| 気化器 | `Q_vap_kW` [kW] | 蒸気消費量換算の基礎 |
| フィード圧縮機 | `W_feed_kW` [kW] | 電力消費量 |
| 製品圧縮機 | `W_prod_kW` [kW] | 電力消費量 |
| 製品冷却器 | `Q_cond_kW` [kW] | 冷却水消費量換算の基礎 |
| 膜モジュール | `n_modules` [本] | 膜交換コスト計算の基礎 |

### 9-2. CAPEX 推算（Turton Bare Module Cost 法）

**一次出典**: プロセス設計R08-3.pdf「プロセス設計(No. 3) 建設費と運転費の推算」（長谷部 伸治・外輪 健一郎）

計算フロー（`src/cost_calculator.py`）:

$$C_{p0} = 10^{K_1 + K_2 \log_{10}A + K_3 (\log_{10}A)^2}$$

$$C_{BM} = C_{p0} \times F_{BM}$$

$$C_{TM} [\text{億円}] = 1.18 \times C_{BM} \times \frac{\text{CEPCI}_{2016}}{\text{CEPCI}_{2001}} \times \frac{\text{JPY}}{\text{USD}} \times 10^{-8}$$

CEPCI: 397.0 (2001基準) → 544.0 (2016年8月)、為替: 110 JPY/USD

| 機器名 | 機器タイプ | サイズ変数 | $K_1$ / $K_2$ / $K_3$ | $F_{BM}$ | 出典ページ |
|---|---|---|---|---|---|
| 気化器（Vaporizer） | 熱交換器（固定管板式） | $A_{vap}$ [m²] | 4.3247 / −0.3030 / 0.1634 | $B_1 + B_2 F_p F_M = 1.63 + 1.66 \times 1.0 \times 1.0$ | p.9, p.10 |
| フィード圧縮機 | 遠心式圧縮機（CS） | $W_{feed}$ [kW] | 2.2897 / 1.3604 / −0.1027 | 2.15（グラフ ID=1） | p.9, p.13–14 |
| 製品圧縮機 | 遠心式圧縮機（CS） | $W_{prod}$ [kW] | 同上 | 2.15 | 同上 |
| 製品冷却器（Condenser） | 熱交換器（固定管板式） | $A_{cond}$ [m²] | 同上（気化器と同型） | 同上（気化器と同値） | p.9, p.10 |
| 膜モジュール | 特殊機器 | $A_{mem}$ [m²] | — | — | **TODO: 単価未確定** |

> **圧縮機の適用範囲注意**: Turton 相関の適用範囲は 450–3000 kW。  
> 100 kmol/h フィード・P_H=10 bar・P_dist=20 bar 条件では W_feed≈182 kW, W_prod≈108 kW となり  
> 範囲外の外挿になる。最終設計流量確定後に再確認すること。

| 出力フィールド | 単位 | 状態 |
|---|---|---|
| `CAPEX_vap` | 億円 | 実装済み（`calc_he_capex_okuyen`） |
| `CAPEX_comp_feed` | 億円 | 実装済み（`calc_comp_capex_okuyen`） |
| `CAPEX_comp_prod` | 億円 | 実装済み（`calc_comp_capex_okuyen`） |
| `CAPEX_cond` | 億円 | 実装済み（`calc_he_capex_okuyen`） |
| `CAPEX_mem` | 億円 | 実装済み（**★ 仮置き**: 単価 50 USD/m²、呼び出し時 UserWarning 発行） |
| `CAPEX_total` | 億円 | 実装済み（5 機器合算、**★ CAPEX_mem が仮置きのため暫定値**） |

---

## 10. 仮定一覧

| # | 仮定 | 根拠・影響 |
|---|---|---|
| 1 | **膜等温操作** | 透過ガスの出口温度はフィード圧縮機出口温度と等しいとする。Hua et al. (2024) の測定条件（室温・大気圧）と整合 |
| 2 | **膜等圧操作** | 非透過側 $P_H$ = 一定、透過側 $P_L$ = 一定（圧力損失なし） |
| 3 | **クロスフローモデル** | 非透過側プラグフロー、透過ガスは即時排出（逆拡散なし）。スパイラル型モジュールの標準近似 |
| 4 | **線形透過則（溶解拡散）** | フラックス = 透過度 × 分圧差（非線形効果は無視） |
| 5 | **透過度は組成・圧力非依存** | $Q_A$, $\alpha$ は膜設計時の代表値を使用 |
| 6 | **van der Waals 混合則 kij=0** | C₃H₆/C₃H₈ は分子構造が近く文献値 kij≈0.01 と小さいため初期設計では省略 |
| 7 | **2 成分系のみ** | 反応器流出ガス中の C₂H₄, CH₄, C₂H₆, H₂, n-C₄H₁₀ は前段蒸留塔で除去済みと仮定 |
| 8 | **気化器 LMTD: 熱媒温度一定** | 蒸気凝縮を熱媒と仮定し $T_{hot}$ = const（LP Steam 160°C, コンテスト仕様） |
| 9 | **冷却器 LMTD: 向流** | ガス入口端 = $T_{in} - T_{cold,out}$、液出口端 = $T_{bp} - T_{cold,in}$（冷却水 30→40°C, コンテスト仕様） |
| 10 | **冷媒不使用（Case A）** | 製品冷却器は冷却水のみ。泡点が $T_{cold,out}$(40°C) を下回る場合は温度クロスが発生しペナルティ返却。冷媒モデルは目的関数の不連続性を生むため採用しない |
| 11 | **CAPEX: 全機器実装済み** | Turton Bare Module Cost 法（プロセス設計R08-3.pdf）。膜モジュール単価（50 USD/m²）と A_per_module（500 m²）は仮置き値。確定次第 `cost_parameters.MEM_UNIT_PRICE_USD_PER_M2` と `MemFixedParams.A_per_module` を更新すること |
| 12 | **圧縮機は多段 + 段間冷却** | 圧縮比 > `max_compression_ratio_per_stage`(!仮置き 4.0) で段数分割、段間で `intercool_T_K`(!仮置き 40°C) まで冷却。段間冷却は HI ホットストリーム `'H_mem_intercool'` に登録。値は全て !仮置き |
| 13 | **vapor_feed フラグ** | 低圧 Dist2 で液フィード不可のとき True にして気化器をスキップ |
| 14 | **膜性能感度係数** | `Q_A_factor`/`alpha_factor`(既定 1.0 = 挙動不変) で混合ガス/高圧/経時劣化を感度解析。実測が出るまで 1.0 固定で文献代表値を使用 |
| 15 | **P_L は 1 atm 固定が設計思想** | 真空ポンプなし。逸脱時は simulate 内で初回のみ警告 |

> `!仮置き`（A_per_module / 膜単価 / 多段冷却 3 値 / 膜性能感度係数）の詳細は
> `KNOWN_PLACEHOLDERS.md`、レビュー課題は `ISSUES_membrane_system.md` を参照。

---

## 11. データクラス定義

### 入力

#### `MemDesignVars` — 設計変数（最適化器が操作）

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `P_H` | float | Pa | 膜供給側（高圧）圧力 |
| `P_L` | float | Pa | 膜透過側（低圧）圧力。**大気圧 1 atm 固定が設計思想**（真空ポンプなし）。全体最適化側は常に 1 atm を渡す。逸脱は simulate 内で一度だけ警告 |
| `A_mem` | float | m² | 総膜面積 |
| `P_dist` | float | Pa | 後段蒸留塔操作圧力（製品圧縮機出口圧力）。C3H6 97%組成では PR EOS 試算で ≳ 17 bar が冷却水使用の目安。P_L より大きくなければならない（`__post_init__` で検証） |

#### `MemFeedStream` — 入力ストリーム（前段蒸留塔底液）

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `F_C3H6` | float | kmol/h | C₃H₆ モル流量 |
| `F_C3H8` | float | kmol/h | C₃H₈ モル流量 |
| `T_in` | float | K | 液温度（飽和液または過冷却液） |
| `P_in` | float | Pa | 圧力 |

#### `MemFixedParams` — 固定パラメータ

| フィールド | デフォルト | 単位 | 説明・根拠 |
|---|---|---|---|
| `Q_A_GPU` | 40.0 | GPU | C₃H₆ 透過度。Hua et al. (2024) 実測値 |
| `alpha` | 90.0 | — | C₃H₆/C₃H₈ 膜選択性。Hua et al. (2024) 実測値 |
| `A_per_module` | 500.0 | m² | 1 モジュールあたり有効膜面積（**★ 仮置き**、Evonik SEPURAN 等カタログで要確認、インスタンス生成時に UserWarning 発行） |
| `T_vap_superheat` | 5.0 | K | 気化器 露点超過過熱度（設計ヒューリスティクス） |
| `U_vap` | 1.0 | kW/(m²·K) | 気化器総括伝熱係数。化工便覧 改訂六版 表6・18（範囲 0.45〜1.14 の中央〜上限値） |
| `T_hot` | 433.15 | K | 熱媒（LP Steam）温度 = 160°C。コンテスト仕様（入手可能スチームのうち最安） |
| `U_cond` | 1.0 | kW/(m²·K) | 冷却器総括伝熱係数。化工便覧 改訂六版 表6・18（範囲 0.45〜1.14 の中央〜上限値） |
| `T_cold_in` | 303.15 | K | 冷却水入口温度 = 30°C。コンテスト仕様 |
| `T_cold_out` | 313.15 | K | 冷却水出口温度 = 40°C。コンテスト仕様 |
| `eta_comp` | 0.75 | — | 圧縮機断熱効率。化工便覧 改訂六版 p.333（ポリトロープ効率 0.7〜0.8 の中央値） |
| `vapor_feed` | False | — | False: 液フィード（気化器で気化）/ True: ガスフィード（気化器スキップ、圧縮機が直接受け取る）。Hua 検証範囲 P_H≤9.5bar 確保のため Dist2 低圧運転 → 液フィード不可のとき True |
| `max_compression_ratio_per_stage` | 4.0 | — | **!仮置き** 1 段あたり最大圧縮比。遠心圧縮機慣行 3〜4 の上端。全圧縮比 ≤ これなら n=1（単段） |
| `intercool_T_K` | 313.15 | K | **!仮置き** 段間冷却の到達温度 = 40°C（冷却水で届く） |
| `U_intercool` | 0.5 | kW/(m²·K) | **!仮置き** 段間冷却器の総括伝熱係数（ガス顕熱-冷却水、ガス-液 0.2 と凝縮 1.0 の中間） |
| `Q_A_factor` | 1.0 | — | **!仮置き** 透過度劣化係数。$Q_{A,eff}=Q_{A,GPU}\times Q_{A,factor}$。可塑化/界面リーク/経時劣化の感度解析用（0<f≤1 で劣化） |
| `alpha_factor` | 1.0 | — | **!仮置き** 選択性劣化係数。$\alpha_{eff}=\alpha\times \alpha_{factor}$（0<f≤1 で劣化） |

> `MemFixedParams` の `__post_init__` は `0 < eta_comp ≤ 1.0`、`T_hot > 273 K`、
> `Q_A_GPU>0`・`alpha>0`、`Q_A_factor>0`・`alpha_factor>0`、
> `max_compression_ratio_per_stage > 1.0`、`U_intercool > 0` を検証する。
> `A_per_module` が仮置き値のとき初回のみ UserWarning。  
> `P_dist` は MemDesignVars に移動済み（最適化変数）。

### 出力

#### `MemSimulationResult` — ルートオブジェクト

| フィールド | 型 | 説明 |
|---|---|---|
| `retentate` | MemRetentateStream | 膜非透過ガス（C₃H₈ 富化） |
| `product` | MemProductStream | 製品（飽和液 @ P_dist） |
| `equipment` | MemEquipmentData | 機器サイズ・OPEX/CAPEX データ |
| `stage_cut` | float | $\theta = F_{perm} / F_{feed}$ [-] |
| `perm_purity` | float | 透過ガス C₃H₆ モル分率 [-] |
| `ret_purity` | float | 非透過ガス C₃H₆ モル分率 [-] |

#### `MemRetentateStream` — 非透過ガス

| フィールド | 単位 | 説明 |
|---|---|---|
| `F_C3H6` | kmol/h | C₃H₆ モル流量 |
| `F_C3H8` | kmol/h | C₃H₈ モル流量 |
| `T_out` | K | 温度（フィード圧縮機出口温度） |
| `P_out` | Pa | 圧力（= $P_H$） |

#### `MemProductStream` — 製品ストリーム（飽和液）

| フィールド | 単位 | 説明 |
|---|---|---|
| `F_C3H6` | kmol/h | C₃H₆ モル流量 |
| `F_C3H8` | kmol/h | C₃H₈ モル流量 |
| `T_out` | K | 泡点温度 @ $P_{dist}$ |
| `P_out` | Pa | 圧力（= $P_{dist}$） |

#### `MemEquipmentData` — 機器サイズ・コスト推算テーブル

| フィールド | 単位 | 説明 |
|---|---|---|
| `A_vap` | m² | 気化器伝熱面積 |
| `Pg_vap` | barg | 気化器ゲージ圧 |
| `Q_vap_kW` | kW | 気化器加熱量（OPEX 用） |
| `W_feed_kW` | kW | フィード圧縮機動力（OPEX 用） |
| `Pg_feed` | barg | フィード圧縮機吐出ゲージ圧 |
| `A_mem` | m² | 総膜面積（= 設計変数） |
| `n_modules` | — | モジュール本数（OPEX 膜交換費用用） |
| `Pg_mem` | barg | 膜モジュール供給側ゲージ圧 |
| `W_prod_kW` | kW | 製品圧縮機動力（OPEX 用） |
| `Pg_prod` | barg | 製品圧縮機吐出ゲージ圧 |
| `A_cond` | m² | 製品冷却器伝熱面積 |
| `Pg_cond` | barg | 製品冷却器ゲージ圧 |
| `Q_cond_kW` | kW | 製品冷却器冷却量（OPEX 用） |
| `CAPEX_vap` | 億円 | 気化器 CAPEX（実装済み） |
| `CAPEX_comp_feed` | 億円 | フィード圧縮機 CAPEX（実装済み） |
| `CAPEX_comp_prod` | 億円 | 製品圧縮機 CAPEX（実装済み） |
| `CAPEX_cond` | 億円 | 製品冷却器 CAPEX（実装済み） |
| `CAPEX_mem` | 億円 | 膜モジュール CAPEX（**★ 仮置き**: 単価 50 USD/m²） |
| `CAPEX_total` | 億円 | 合計資本費の合算（**★ CAPEX_mem が仮置きのため暫定値**） |
| `Q_intercool_kW` | kW | 段間冷却 総熱量（フィード+製品圧縮機、冷却水 OPEX 用、n=1 なら 0） |
| `A_intercool_m2` | m² | 段間冷却器 総伝熱面積 |
| `CAPEX_intercool` | 億円 | 段間冷却器 CAPEX（n=1 なら 0） |
| `T_intercool_in_K` / `T_intercool_out_K` | K | 段間冷却 入口（各段出口温度の熱量加重平均）/ 出口（= `intercool_T_K`）。HI 登録用 |

---

## 12. 定数・パラメータ一覧

### `src/config.py` の追加フィールド（`ThermoParams`）

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `Tc` | float | K | 臨界温度（デフォルト `nan`） |
| `Pc` | float | Pa | 臨界圧力（デフォルト `nan`） |
| `omega` | float | — | アセントリック因子（デフォルト `nan`） |

### モジュール内定数

| 定数 | 値 | 説明 |
|---|---|---|
| `_KEYS` | `['B', 'A']` | EOS 成分キー順序（C₃H₆, C₃H₈） |
| `_T_REF` | 298.15 K | エンタルピー基準温度 |
| `_GPU_SI` | 3.35×10⁻¹⁰ mol/(m²·s·Pa) | 1 GPU の SI 換算値 |
| `_ATM_BAR` | 1.01325 bar | ゲージ圧変換用大気圧 |
| `_PENALTY` | 1×10⁹ | ペナルティ値（最適化への無効シグナル） |

### `src/eos.py` 内定数（PR EOS 固有）

| 定数 | 値 | 説明 |
|---|---|---|
| `_OA` | 0.45724 | PR $\Omega_a$ |
| `_OB` | 0.07780 | PR $\Omega_b$ |
| `_SQRT2` | $\sqrt{2}$ | 対数項の定数 |

---

## 13. 依存モジュール

| モジュール | 使用箇所 | 役割 |
|---|---|---|
| `src.eos.z_factor` | `_h_mol`, 圧縮機 | PR Z 因子計算 |
| `src.eos.residual_enthalpy` | `_h_mol` | 残差エンタルピー計算 |
| `src.eos.bubble_point_T` | `_condenser` | 製品飽和液温度計算 |
| `src.eos.dew_point_T` | `_vaporizer` | 気化器出口温度計算 |
| `src.eos.compress_isentropic` | `_compress_multistage`（フィード/製品圧縮機） | 各段の出口温度・仕事計算 |
| `src.cost_calculator.calc_he/comp/mem_capex_okuyen` | CAPEX 計算 | 気化器/冷却器/段間冷却器・圧縮機・膜の Bare Module Cost |
| `src.thermo.PDHThermo` | `_h_mol` の理想気体項 | $C_p$ 積分（`calc_enthalpy_change`） |
| `src.config.THERMO_DATA` | `src/eos.py` 全体 | 臨界定数・$C_p$ 係数 |
| `scipy.integrate.solve_ivp` | `_membrane_ode` | クロスフロー ODE ソルバー（Radau） |
| `scipy.optimize.brentq` | `src/eos.py` 全体 | 泡点/露点/断熱圧縮の根探索 |
| `numpy` | `simulate_membrane_system` | `np.clip` によるクリッピング |

---

## 14. 使い方

### 基本的な呼び出し

```python
from units.separators.membrane.membrane_system import (
    MemDesignVars, MemFeedStream, MemFixedParams,
    simulate_membrane_system,
)
from src.eos import bubble_point_T

# 前段蒸留塔底液の泡点温度を計算して T_in に使用
P_in = 5.0e5  # 5 bar
T_in = bubble_point_T(P_in, [0.5, 0.5], ['B', 'A'])

feed = MemFeedStream(
    F_C3H6 = 500.0,   # kmol/h
    F_C3H8 = 500.0,   # kmol/h
    T_in   = T_in,    # 飽和液
    P_in   = P_in,    # Pa
)
design = MemDesignVars(
    P_H    = 15.0e5,   # Pa (15 bar)
    P_L    =  1.5e5,   # Pa (1.5 bar)
    A_mem  = 30000.0,  # m²
    P_dist = 20.0e5,   # Pa (20 bar) — 冷却水使用のため ≳ 17 bar 目安
)

result = simulate_membrane_system(design, feed, MemFixedParams())

print(f"stage_cut   = {result.stage_cut:.4f}")
print(f"perm_purity = {result.perm_purity*100:.2f}%")
print(f"W_feed_kW   = {result.equipment.W_feed_kW:.0f} kW")
print(f"Q_vap_kW    = {result.equipment.Q_vap_kW:.0f} kW")
```

### 最適化ループでの想定使用パターン

```python
def objective(x):
    P_H, P_L, A_mem, P_dist = x
    design = MemDesignVars(P_H=P_H, P_L=P_L, A_mem=A_mem, P_dist=P_dist)
    r = simulate_membrane_system(design, feed, fixed)
    purity_penalty = max(0.0, 0.99 - r.perm_purity) * 1e6
    return r.equipment.W_feed_kW + r.equipment.W_prod_kW + purity_penalty

# 無効条件では CAPEX_total = 1e9 億円が返るので最適化器に安全に渡せる
```

### パラメータ変更例

```python
# 膜性能を変更する場合
fixed = MemFixedParams(Q_A_GPU=60.0, alpha=100.0)

# 後段蒸留塔操作圧力を変更する場合（MemDesignVars で指定）
design = MemDesignVars(P_H=15.0e5, P_L=1.5e5, A_mem=30000.0, P_dist=22.0e5)

# 圧縮機効率を変更する場合
fixed = MemFixedParams(eta_comp=0.80)
```

---

## 15. エラーハンドリング

### 方針

無効入力・数値異常はすべて `_penalty_result()` を返して**クラッシュしない**。  
最適化器に対して大きなペナルティ値を送り、探索を無効領域から引き戻す。

### 入力バリデーション

`simulate_membrane_system` 冒頭で以下をチェック。条件に合致した場合は即 `_penalty_result()` を返す。

| チェック対象 | 条件 |
|---|---|
| `design.P_H` と `design.P_L` | `P_H <= P_L`（駆動力なし） |
| `design.A_mem`, `design.P_H`, `design.P_L` | ≤ 0 |
| `feed.F_C3H6`, `feed.F_C3H8` | 負値 |
| `F_total_feed` | ≤ 0（空フィード） |

### ユニットレベルの保護

| ユニット | 保護内容 |
|---|---|
| 気化器 | `T_vap_out >= T_hot` の場合に `UserWarning` を発行し NaN を返す → ペナルティ |
| 気化器・冷却器 | `_h_mol` 全体を `try/except` で囲む |
| フィード/製品圧縮機 | `compress_isentropic` 全体を `try/except` で囲む |
| 膜 ODE | `solve_ivp` 失敗時に `(None, None)` を返す → ペナルティ |
| 製品冷却器 | `_condenser` 全体を `try/except` で囲む |

### ODE 数値安定化

| 対象 | 処理 |
|---|---|
| 非透過モル流量 `F` | `max(F, 1e-12)` — 負値・ゼロ割りを防止 |
| フラックス $J$ | $J \leq 0$ のとき `[0, 0]` を返す — 逆拡散を無視 |
| $y_{local}$ の分母 | `abs(denom) < 1e-30` のとき $y = x$（フォールバック） |

### ペナルティ値

```
各 CAPEX_*  = nan  [億円]（計算例外時）
CAPEX_total  = nan  [億円]（計算例外時）
stage_cut    = 0.0
perm_purity  = 0.0
ret_purity   = 0.0
全ストリーム  = 0 kmol/h
```
