# SPEC: eos.py — Peng-Robinson 状態方程式モジュール

**ファイルパス**: `src/eos.py`
**最終更新**: 2026-06-02 (新規作成。実コードからの式・定数・入出力の抽出)

---

## 目次

1. [目的](#1-目的)
2. [モデル・理論](#2-モデル理論)
   - 2-1. PR 状態方程式と混合則
   - 2-2. Z 因子（3 次方程式の実根）
   - 2-3. フガシティー係数
   - 2-4. 残差熱力学量
   - 2-5. 泡点温度（thermo 連携）
   - 2-6. 露点温度
   - 2-7. 断熱圧縮
3. [入力・出力](#3-入力出力)
4. [主要パラメータ](#4-主要パラメータ)
5. [thermo パッケージ連携](#5-thermo-パッケージ連携)
6. [出典](#6-出典)
7. [既知の限界・仮置き](#7-既知の限界仮置き)

---

## 1. 目的

Peng-Robinson (PR) 状態方程式に基づく気液平衡・熱力学量を提供する純関数モジュール。
主に **膜分離システム・圧縮機・蒸留塔 rigorous 求解** が参照する。理想気体近似では
誤差が大きくなる中〜高圧領域（PSA/膜: <25 bar、圧縮機昇圧、蒸留塔 8.5〜20 bar）の
フガシティー・K 値・圧縮仕事を扱う。

`config.py` の `THERMO_DATA` に格納された各成分の Tc / Pc / ω（アセントリック因子）と
Cp 多項式係数（理想気体エンタルピー/エントロピー差分用）を入力に使う。

**単位系**: 温度 K、圧力 Pa、エンタルピー/エントロピー J/mol または J/(mol·K)。
モル流量 [kmol/h] → [mol/s] 変換は呼び出し側の責任。

---

## 2. モデル・理論

### 2-1. PR 状態方程式と混合則

単成分パラメータ（`_pr_single`）:

$$a_c = \Omega_a \frac{R^2 T_c^2}{P_c}, \qquad b = \Omega_b \frac{R T_c}{P_c}$$

$$\kappa = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^2$$

$$\sqrt{a_i} = \sqrt{a_c}\left(1 + \kappa\left(1 - \sqrt{T/T_c}\right)\right), \qquad
\frac{d\sqrt{a_i}}{dT} = -\frac{\kappa\sqrt{a_c}}{2\sqrt{T\,T_c}}$$

van der Waals 混合則（`_mix`、**k_ij = 0**）:

$$S = \sum_i x_i \sqrt{a_i}, \quad a_m = S^2, \quad b_m = \sum_i x_i b_i, \quad
\frac{da_m}{dT} = 2 S \sum_i x_i \frac{d\sqrt{a_i}}{dT}$$

無次元化:

$$A = \frac{a_m P}{R^2 T^2}, \qquad B = \frac{b_m P}{R T}$$

| 記号 | 意味 | 単位 |
|---|---|---|
| $\Omega_a, \Omega_b$ | PR の臨界点定数（0.45724 / 0.07780） | — |
| $a_c, b$ | 単成分引力/排除体積パラメータ | Pa·m⁶/mol² / m³/mol |
| $\kappa$ | ソアベ型温度補正係数 | — |
| $\omega$ | アセントリック因子 | — |
| $T_c, P_c$ | 臨界温度・臨界圧力 | K / Pa |
| $a_m, b_m$ | 混合物パラメータ | — |
| $x_i$ | 成分 $i$ のモル分率 | — |
| $A, B$ | 無次元化 PR パラメータ | — |

### 2-2. Z 因子（3 次方程式の実根）

`z_factor` / 内部 `_cubic_z`。次の 3 次方程式の実根のうち物理下限 $Z > B$ を満たすものを採用:

$$Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0$$

- `phase='vapor'` → 最大根、`phase='liquid'` → 最小正根
- 実根なし（臨界点近傍・低温真空などで cubic root 消失）→ **Z = 1（理想気体）で fallback し UserWarning**

### 2-3. フガシティー係数

`fugacity_coeff`:

$$\ln\varphi_i = \frac{b_i}{b_m}(Z-1) - \ln(Z-B)
- \frac{A}{2\sqrt{2}\,B}\left(\frac{2\sqrt{a_i}}{S} - \frac{b_i}{b_m}\right)\ln L$$

$$L = \frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B}$$

K 値は $K_i = \varphi_i^L / \varphi_i^V$。

### 2-4. 残差熱力学量

`residual_enthalpy` / `residual_entropy`:

$$H^r = R T (Z-1) + \frac{T\,\dfrac{da_m}{dT} - a_m}{2\sqrt{2}\,b_m}\ln L \quad [\mathrm{J/mol}]$$

$$S^r = R\ln(Z-B) + \frac{\dfrac{da_m}{dT}}{2\sqrt{2}\,b_m}\ln L \quad [\mathrm{J/(mol\cdot K)}]$$

理想気体差分（`_dh_ig`, `_ds_ig`）は `THERMO_DATA` の Cp 多項式 $C_p = a + bT + cT^2 + dT^3$ を解析積分:

$$\Delta H^{ig} = \sum_i x_i\!\left[a(T_2-T_1) + \tfrac{b}{2}(T_2^2-T_1^2) + \tfrac{c}{3}(T_2^3-T_1^3) + \tfrac{d}{4}(T_2^4-T_1^4)\right]$$

$$\Delta S^{ig} = \sum_i x_i\!\left[a\ln\tfrac{T_2}{T_1} + b(T_2-T_1) + \tfrac{c}{2}(T_2^2-T_1^2) + \tfrac{d}{3}(T_2^3-T_1^3)\right] - R\ln\tfrac{P_2}{P_1}$$

### 2-5. 泡点温度（thermo 連携）

`bubble_point_T(P, x, keys, T_lo=150, T_hi=500)`。収束条件 $\sum_i x_i K_i = 1$。

- 外側ループ: 温度 $T$ を `scipy.optimize.brentq`（xtol=0.05 K, maxiter=200）で探索
- 内側ループ: 気相組成 $y$ を Wilson 相関で初期化 → 逐次置換（最大 50 回）で $\varphi^V$ を更新
- **液相/気相のフガシティーは `thermo` パッケージの `PRMIX`（`phis_l` / `phis_g`）で取得**。
  単相→二相遷移境界（$Z_V = Z_L$ 縮退点）で手作り PR の brentq が偽根を返す問題を回避するため。
  単相領域では PRMIX が `AttributeError`/`None` を返すので、それを brentq の bracketing 用の
  符号付き発散値（±2.0）に変換する。
- 収束失敗時は **nan を返し UserWarning**。

### 2-6. 露点温度

`dew_point_T(P, y, keys, T_lo=150, T_hi=500)`。収束条件 $\sum_i y_i / K_i = 1$。
こちらは **手作り PR**（`z_factor` + `fugacity_coeff`）を使用。
外側 brentq（xtol=0.05 K）、内側で $\varphi^L$ を逐次置換（最大 50 回）。実根 < 2 本（単相）の
区間では符号付き発散値を返す。

### 2-7. 断熱圧縮

`compress_isentropic(T1, P1, P2, x, keys, eta=0.80)`:

1. 等エントロピー条件 $\Delta S^{ig} + \Delta S^r = 0$ を brentq で解き $T_{2s}$ を逆算
   （初期推定は理想気体多変断熱 $T_{2s} = T_1 (P_2/P_1)^{(\kappa-1)/\kappa}$、$\kappa = C_p/(C_p - R)$）
2. $W_{isen} = \Delta H^{ig}(T_1\to T_{2s}) + H^r(T_{2s},P_2) - H^r(T_1,P_1)$
3. $W_{actual} = W_{isen} / \eta$
4. エンタルピー収支 $\Delta H^{ig} + \Delta H^r = W_{actual}$ から実出口温度 $T_2$ を逆算

返り値: $(T_{2,actual}\,[\mathrm{K}],\ W_{actual}\,[\mathrm{J/mol}])$（正値 = 圧縮機への入力仕事）。

---

## 3. 入力・出力

| 関数 | 主な入力 | 出力 |
|---|---|---|
| `z_factor(T, P, x, keys, phase)` | T[K], P[Pa], モル分率 x, 成分キー keys, phase='vapor'/'liquid' | Z [-] |
| `fugacity_coeff(i, T, P, x, keys, Z)` | 成分 index i ほか | φᵢ [-] |
| `residual_enthalpy(T, P, x, keys, Z)` | 同上 | Hʳ [J/mol] |
| `residual_entropy(T, P, x, keys, Z)` | 同上 | Sʳ [J/(mol·K)] |
| `bubble_point_T(P, x, keys, T_lo, T_hi)` | P[Pa], 液相組成 x | 泡点 T_bp [K]（失敗時 nan） |
| `dew_point_T(P, y, keys, T_lo, T_hi)` | P[Pa], 気相組成 y | 露点 T_dp [K]（失敗時 nan） |
| `compress_isentropic(T1, P1, P2, x, keys, eta)` | 入口 T1/P1, 出口 P2, η | (T2_actual [K], W_actual [J/mol]) |

`keys` は成分キーのリスト（例 `['B','A']` = [C3H6, C3H8]）。`x`/`y` は同順のモル分率（Σ=1）。

---

## 4. 主要パラメータ

| 記号/定数 | 値 | 説明 | 出典 |
|---|---|---|---|
| `_OA` (Ω_a) | 0.45724 | PR 引力項定数 | Peng-Robinson (1976) Eq.(5)-(6) |
| `_OB` (Ω_b) | 0.07780 | PR 排除体積定数 | 同上 |
| `R` | 8.31446 J K⁻¹ mol⁻¹ | 気体定数（`config.py`） | — |
| k_ij | 0 | 二成分相互作用係数（van der Waals 古典近似） | !仮置き |
| bubble/dew brentq xtol | 0.05 K | 温度探索の許容誤差 | 数値ガード |
| 内側逐次置換 | 50 回 / 収束 1e-7 | 組成収束 | 数値ガード |
| `eta` 既定 | 0.80 | 断熱圧縮効率（呼び出し側で上書き） | — |
| K_eq exp 引数クランプ | 700 | `thermo.py` 側 OverflowError 防止 | 数値ガード |

各成分の Tc / Pc / ω は `config.py` の `THERMO_DATA` を参照（H2/CH4/C2 系は化学工学便覧 改訂六版 表1.3、C3 系は NIST）。

---

## 5. thermo パッケージ連携

- `bubble_point_T` のみ **`thermo.eos_mix.PRMIX`（CalebBell/thermo, MIT, v0.6.0 ピン）** を内部利用。
  単相⇔二相切替境界での偽根問題を根本解決するため。
- 他関数（`z_factor`, `fugacity_coeff`, 残差量, `dew_point_T`, `compress_isentropic`）は
  **すべて自前 PR 実装**。両者は PR EOS で 0.02% 一致を確認済（KNOWN_PLACEHOLDERS §A.1）。
- `thermo` 未インストール時は **モジュールロード時に `ImportError` で fail-fast**。
  （関数内 lazy import だと呼び出し側 except に握り潰され、PSA/膜の偽 penalty として現れるため。）
- 蒸留 rigorous（`distillation_rigorous.py`）は K 値計算に本モジュールの `bubble_point_T` /
  `fugacity_coeff` / `z_factor` を流用する。

---

## 6. 出典

- **Peng, D.-Y. & Robinson, D.B. (1976)** "A New Two-Constant Equation of State,"
  *Ind. Eng. Chem. Fundam.* 15(1), 59-64. — Ω_a / Ω_b、PR 式形、κ(ω) 相関。
- **CalebBell/thermo** (https://github.com/CalebBell/thermo, MIT License, v0.6.0) — `bubble_point_T` の PRMIX。
- 成分臨界物性: 化学工学便覧 改訂六版 表1.3（H2/CH4/C2H4/C2H6/C4H10）、NIST（C3H8/C3H6）。
- Wilson 相関（初期 K 値推算）: 標準的な VLE 初期推定式。

---

## 7. 既知の限界・仮置き

- **!仮置き k_ij = 0**（van der Waals 古典近似）。文献では C3 系で 0.01 程度と小さく無視可と
  見られるが、出典確定はユーザー判断要。完全実装には `_mix` を二重和
  $a_m = \sum_i\sum_j x_i x_j (1-k_{ij})\sqrt{a_i a_j}$ に書き換え + k_ij テーブル追加が必要。
- **泡点/露点の探索範囲デフォルト [150K, 500K]** は C3H6/C3H8 系を想定。H2/CH4 など極低温成分
  のみを含む組成は単相領域に張り付き適用外（nan を返す）。
- **Z=1 silent fallback**: cubic root 消失領域で理想気体に退化。下流（PSA/膜の単純流量計算）では
  φ≈1 仮定の誤差が直接 CAPEX に伝播し得る（定量はユーザー検証要）。蒸留 rigorous では MESH 残差で
  検出され FUG fallback する。
- **bubble/dew point brentq xtol = 0.05 K**（数値ガード、出典なし。蒸留段精度に直結）。

> 上記「!仮置き」「【確認中】」項目の最新状況・置換ルールは
> [`../KNOWN_PLACEHOLDERS.md`](../KNOWN_PLACEHOLDERS.md)（§1.4 物性パラメータ、§3 数値ガード）を参照。
