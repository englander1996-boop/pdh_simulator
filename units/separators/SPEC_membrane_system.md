# SPEC: membrane_system.py — PDH 膜分離システム シミュレーター

**ファイルパス**: `units/separators/membrane_system.py`  
**依存 EOS モジュール**: `src/eos.py`  
**最終更新**: 2026-05-01

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
| 設計変数 | 膜供給側圧力 $P_H$、透過側圧力 $P_L$、総膜面積 $A_{mem}$ |
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
  ┌──────────────┐  等エントロピー圧縮 (PR 補正)
  │フィード圧縮機 │──▶ T_feed_comp_out, W_feed_kW
  └──────────────┘
       │ ガス @ P_H
       ▼
  ┌────────────────┐  クロスフロー ODE (solve_ivp Radau)
  │ 膜モジュール   │──▶ 非透過 (Retentate) / 透過 (Permeate)
  └────────────────┘
       │ Permeate @ P_L
       ▼
  ┌────────────────┐  等エントロピー圧縮 (PR 補正)
  │ 製品圧縮機     │──▶ T_prod_comp_out, W_prod_kW
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
| $U_{vap}$ | 気化器総括伝熱係数（デフォルト 1.5） | kW/(m²·K) |
| $T_{hot}$ | 熱媒（低圧蒸気）温度（デフォルト 423.15 K = 150°C） | K |

---

## 6. 圧縮機モデル

フィード圧縮機・製品圧縮機の両方で `compress_isentropic`（[4-9節](#4-9-断熱圧縮-compress_isentropic)）を共通使用。

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
| $Q_A$ | C₃H₆ 透過度（デフォルト 40 GPU → $1.340 \times 10^{-8}$ mol/(m²·s·Pa)） | mol m⁻² s⁻¹ Pa⁻¹ |
| $Q_B$ | C₃H₈ 透過度（$= Q_A / \alpha$） | mol m⁻² s⁻¹ Pa⁻¹ |
| $\alpha$ | C₃H₆/C₃H₈ 膜選択性（デフォルト 90） | — |
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
| $U_{cond}$ | 冷却器総括伝熱係数（デフォルト 1.0） | kW/(m²·K) |
| $T_{cold,in}$ | 冷却水入口温度（デフォルト 303.15 K = 30°C） | K |
| $T_{cold,out}$ | 冷却水出口温度（デフォルト 318.15 K = 45°C） | K |

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

### 9-2. CAPEX 用データ（授業資料 R08-3 形式）

| 機器名 | 機器タイプ | 特徴サイズ $A$ | 単位 | ゲージ圧 $P_g$ [barg] |
|---|---|---|---|---|
| Vaporizer | 熱交換器 | `A_vap` | m² | `Pg_vap` = $P_{in}/10^5 - 1.01325$ |
| Feed Compressor | 圧縮機 | `W_feed_kW` | kW | `Pg_feed` = $P_H/10^5 - 1.01325$ |
| Membrane | 特殊機器 | `A_mem` | m² | `Pg_mem` = $P_H/10^5 - 1.01325$ |
| Product Compressor | 圧縮機 | `W_prod_kW` | kW | `Pg_prod` = $P_{dist}/10^5 - 1.01325$ |
| Condenser | 熱交換器 | `A_cond` | m² | `Pg_cond` = $P_{dist}/10^5 - 1.01325$ |

> **注**: ゲージ圧は絶対圧 [bar] から大気圧 1.01325 bar を引いて算出する。  
> CAPEX 相関係数（$K_1$, $K_2$, $K_3$, $B_1$, $B_2$ 等）は `src/cost_parameters.py` に実装予定（現在 `CAPEX_total = nan`）。

---

## 10. 仮定一覧

| # | 仮定 | 根拠・影響 |
|---|---|---|
| 1 | **膜等温操作** | 透過ガスの出口温度はフィード圧縮機出口温度と等しいとする |
| 2 | **膜等圧操作** | 非透過側 $P_H$ = 一定、透過側 $P_L$ = 一定（圧力損失なし） |
| 3 | **クロスフローモデル** | 非透過側プラグフロー、透過ガスは即時排出（逆拡散なし） |
| 4 | **線形透過則（溶解拡散）** | フラックス = 透過度 × 分圧差（非線形効果は無視） |
| 5 | **透過度は組成・圧力非依存** | $Q_A$, $\alpha$ は膜設計時の代表値を使用 |
| 6 | **van der Waals 混合則 kij=0** | C₃H₆/C₃H₈ 近似成分のため交差パラメータを省略 |
| 7 | **2 成分系のみ** | 反応器流出ガス中の C₂H₄, CH₄, C₂H₆, H₂, n-C₄H₁₀ は前段蒸留塔で除去済みと仮定 |
| 8 | **気化器 LMTD: 熱媒温度一定** | 蒸気凝縮を熱媒と仮定し $T_{hot}$ = const |
| 9 | **冷却器 LMTD: 向流** | ガス入口端 = $T_{in} - T_{cold,out}$、液出口端 = $T_{bp} - T_{cold,in}$ |
| 10 | **CAPEX は現在 NaN** | `cost_parameters.py` 拡張後に実装予定 |

---

## 11. データクラス定義

### 入力

#### `MemDesignVars` — 設計変数（最適化器が操作）

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `P_H` | float | Pa | 膜供給側（高圧）圧力 |
| `P_L` | float | Pa | 膜透過側（低圧）圧力 |
| `A_mem` | float | m² | 総膜面積 |

#### `MemFeedStream` — 入力ストリーム（前段蒸留塔底液）

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `F_C3H6` | float | kmol/h | C₃H₆ モル流量 |
| `F_C3H8` | float | kmol/h | C₃H₈ モル流量 |
| `T_in` | float | K | 液温度（飽和液または過冷却液） |
| `P_in` | float | Pa | 圧力 |

#### `MemFixedParams` — 固定パラメータ

| フィールド | デフォルト | 単位 | 説明 |
|---|---|---|---|
| `Q_A_GPU` | 40.0 | GPU | C₃H₆ 透過度 |
| `alpha` | 90.0 | — | C₃H₆/C₃H₈ 膜選択性 |
| `A_per_module` | 500.0 | m² | 1 モジュールあたり有効膜面積 |
| `P_dist` | 15.0×10⁵ | Pa | 後段蒸留塔操作圧力 |
| `T_vap_superheat` | 5.0 | K | 気化器 露点超過過熱度 |
| `U_vap` | 1.5 | kW/(m²·K) | 気化器総括伝熱係数 |
| `T_hot` | 423.15 | K | 熱媒（低圧蒸気）温度 |
| `U_cond` | 1.0 | kW/(m²·K) | 冷却器総括伝熱係数 |
| `T_cold_in` | 303.15 | K | 冷却水入口温度 |
| `T_cold_out` | 318.15 | K | 冷却水出口温度 |
| `eta_comp` | 0.75 | — | 圧縮機断熱効率 |

> `MemFixedParams` は `__post_init__` で `P_dist > 0`、`0 < eta_comp ≤ 1.0`、`T_hot > 273 K` を検証する。

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
| `CAPEX_total` | 億円 | 合計資本費（現在 `nan`） |

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
| `src.eos.compress_isentropic` | フィード/製品圧縮機 | 圧縮機出口温度・仕事計算 |
| `src.thermo.PDHThermo` | `_h_mol` の理想気体項 | $C_p$ 積分（`calc_enthalpy_change`） |
| `src.config.THERMO_DATA` | `src/eos.py` 全体 | 臨界定数・$C_p$ 係数 |
| `scipy.integrate.solve_ivp` | `_membrane_ode` | クロスフロー ODE ソルバー（Radau） |
| `scipy.optimize.brentq` | `src/eos.py` 全体 | 泡点/露点/断熱圧縮の根探索 |
| `numpy` | `simulate_membrane_system` | `np.clip` によるクリッピング |

---

## 14. 使い方

### 基本的な呼び出し

```python
from units.separators.membrane_system import (
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
    P_H   = 15.0e5,   # Pa (15 bar)
    P_L   =  1.5e5,   # Pa (1.5 bar)
    A_mem = 30000.0,  # m²
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
    P_H, P_L, A_mem = x
    design = MemDesignVars(P_H=P_H, P_L=P_L, A_mem=A_mem)
    r = simulate_membrane_system(design, feed, fixed)
    purity_penalty = max(0.0, 0.99 - r.perm_purity) * 1e6
    return r.equipment.W_feed_kW + r.equipment.W_prod_kW + purity_penalty

# 無効条件では CAPEX_total = 1e9 億円が返るので最適化器に安全に渡せる
```

### パラメータ変更例

```python
# 膜性能を変更する場合
fixed = MemFixedParams(Q_A_GPU=60.0, alpha=100.0)

# 後段蒸留塔操作圧力を変更する場合
fixed = MemFixedParams(P_dist=18.0e5)

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
CAPEX_total  = 1×10⁹  [億円]
stage_cut    = 0.0
perm_purity  = 0.0
ret_purity   = 0.0
全ストリーム  = 0 kmol/h
```
