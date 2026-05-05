# エラーハンドリング修正レポート

**作成日**: 2026-05-05  
**対象**: PSA システム最適化変数スイープにおけるフリーズ・クラッシュ・サイレントエラーの修正

---

## 修正内容一覧

### `units/separators/psa/psa_system.py`

---

#### 【修正 P-1】空塔速度 `u_0` 上限チェック追加（LSODA フリーズ防止）

**問題**  
`u_0 = F / A_col` は `D_col²` に反比例して増大する。`u_0` が大きいと PDE の移流項が
硬くなり、LSODA が極めて小さいタイムステップを繰り返し、7200 秒積分に数百万ステップを
要してフリーズする。タイムアウト機構がないため最適化ループ全体が停止する。

**修正内容**  
`u_0` 計算直後に上限チェックを追加。超過時は `_penalty_result()` を即時返却する。

```python
_U0_MAX = 2.0  # [m/s]  ← 後述「要判断事項」参照
if u_0 > _U0_MAX:
    return _penalty_result()
```

---

#### 【修正 P-2】CSS 補正後 `t_abs` 最小値チェック追加（`scale` 発散防止）

**問題**  
`t_abs = t_abs_clean * (1 - desorption_target)` において `desorption_target` が 1 に
近い場合 `t_abs ≈ 0` になり、`scale = 3600 / (t_abs * 1000)` が発散して産物流量が
入力流量を超える誤結果が警告なく返されていた。

**修正内容**  
CSS 補正直後に最小値チェックを追加。`t_abs < 1.0 s` の場合はペナルティを返す。

```python
_T_ABS_MIN = 1.0  # [s]
if t_abs < _T_ABS_MIN:
    return _penalty_result()
```

---

#### 【修正 P-3】`total_moles_out` を CSS 補正後の `t_abs` まで再積分

**問題**  
旧コードの `total_moles_out` は `[0, t_abs_clean]`（清浄床破過まで）の積分だった。
CSS 補正で `t_abs < t_abs_clean` に短縮しても積分範囲を修正していなかったため、
実際には操作しない区間 `[t_abs, t_abs_clean]` の出口成分が産物流量に含まれ、不純物を
過大推算していた。

**修正内容**  
`_run_adsorption` の返り値を `(t_abs, q_final, sol_t, C_outlet, converged)` に変更し、
呼び出し側で `[0, t_abs]` まで切り詰めて `np.trapezoid` で積分し直す。
端点 `t_abs` が `sol.t` のサンプル点の間にある場合は線形補間で補完する。

> **仮定**: 破過前の出口濃度はほぼゼロのため、この修正が産物純度の推算に与える影響は
> 小さい。C2H4・C2H6 が CH4 より先に破過するケースでは差が出る可能性がある。

---

#### 【修正 P-4】`rhs` 内の負値クランプ追加

**問題**  
LSODA は `C ≥ 0`、`q ≥ 0` を強制しない。移流支配の PDE では数値誤差で濃度が
わずかに負転することがあり、Langmuir 式で `q_eq < 0` → LDF で `q` がさらに負に
引き込まれる不安定フィードバックが起きうる。

**修正内容**  
```python
C = np.maximum(y[:_N_Z * _N_ADS].reshape(_N_Z, _N_ADS), 0.0)
q = np.maximum(y[_N_Z * _N_ADS:].reshape(_N_Z, _N_ADS), 0.0)
```

> **注意**: クランプは RHS 評価時のみ行うため、`sol.y` の状態ベクトル自体は
> 負値を含むことがある。必要であれば `_run_adsorption` の返り値 `C_outlet` にも
> `np.maximum(..., 0.0)` を適用すること。

---

#### 【修正 P-5】H2 損失の物質収支クランプ

**問題**  
`H2_loss_purge` が `F_H2_in` を超える場合、オフガス中の H2 損失が入力流量を超え
物質収支が壊れていた。`F_H2_product` は `max(0, ...)` でクランプされていたが
`offgas['C']` は未クランプのまま大きな値が残っていた。

**修正内容**  
```python
H2_loss_purge = min(H2_loss_purge, max(0.0, F_H2_in - H2_loss_blowdown))
```

---

#### 【修正 P-6】`z_factor` フォールバック時の警告追加

**問題**  
PR EOS が収束失敗した場合、`Z=1.0`（理想気体）にサイレントフォールバックしており
`u_0` と `C_feed_ads` に体系的な誤差が混入してもログに出なかった。

**修正内容**  
```python
except Exception as e:
    warnings.warn(f"_calc_feed_state: z_factor 計算失敗 ({e})。理想気体 Z=1.0 を使用。...")
```

---

### `src/eos.py`

---

#### 【修正 E-1】`_cubic_z` 虚数根フィルタを相対閾値に変更

**問題**  
`abs(r.imag) < 1e-8` の絶対閾値は、根の絶対値が大きい高圧条件で数値誤差が
`1e-8` を超え、有効な実根を除外することがあった。roots が空になると
`z_factor` が `Z=1.0` を無音で返す。

**修正内容**  
```python
abs(r.imag) < 1e-8 * max(abs(r.real), 1.0)
```

---

#### 【修正 E-2】`bubble_point_T` / `dew_point_T` の内側ループ改善

**問題**  
逐次置換ループが 20 回固定で収束チェックなし。未収束の K 値が `obj(T)` に
渡ると brentq が誤根を見つけるか ValueError を投げる可能性があった。

**修正内容**  
- 最大反復を 50 回に増加
- `max(|y_new - y|) < 1e-7` で早期終了

---

#### 【修正 E-3】`bubble_point_T` / `dew_point_T` の brentq に try/except 追加

**問題**  
探索範囲 `[T_lo=150K, T_hi=500K]` が固定のため、対象組成の泡点/露点が範囲外に
ある場合 `ValueError` が uncaught で伝播し最適化ループがクラッシュしていた。

**修正内容**  
```python
try:
    return brentq(obj, T_lo, T_hi, xtol=0.05, maxiter=200)
except ValueError:
    warnings.warn(f"bubble_point_T: brentq 収束失敗 [...] K。nan を返します。")
    return float('nan')
```

> **要判断**: `nan` を返す設計にしたため、呼び出し側でのハンドリングが必要。
> 膜分離システム等で `bubble_point_T` / `dew_point_T` を呼んでいる箇所を確認し、
> `nan` チェックを追加すること。

---

#### 【修正 E-4】`compress_isentropic` の 2 段目 brentq に try/except 追加

**問題**  
1 段目の brentq 失敗後の広範囲再試行（`[T1+0.5, 1200K]`）にも try/except がなく、
これも失敗した場合 `ValueError` が伝播していた。

**修正内容**  
```python
except ValueError:
    try:
        T2s = brentq(entropy_balance, T1 + 0.5, 1200.0, xtol=0.1, maxiter=300)
    except ValueError:
        warnings.warn("compress_isentropic: entropy_balance brentq 収束失敗。理想気体近似で T2s を推算。")
        T2s = T1 * (P2 / P1) ** ((kappa_approx - 1.0) / kappa_approx)
```

> **仮定**: フォールバックに用いた `kappa_approx = 1.13` は C3 混合の代表値。
> フィード組成が C3 から大きく外れる場合（例: H2 リッチな PSA フィード）は
> この近似が不適切になる。必要に応じて `kappa_approx` を組成から計算すること。

---

### `src/thermo.py`

---

#### 【修正 T-1】`calc_keq` の exp オーバーフロー防止

**問題**  
`math.exp(-dG_T / (R * T))` は T が非常に高い（または ΔG_T が極端に負な）
条件で `OverflowError` を投げる可能性があった（`math.exp` の最大引数 ≈ 709）。

**修正内容**  
```python
exp_arg = min(-dG_T / (R * T), 700.0)
return _P_STD * math.exp(exp_arg)
```

> **影響**: PDH 操作温度（500〜700°C）では exp_arg ≈ 3〜10 程度であり、クランプは
> 発動しない。最適化で T が極端な値を取った場合のみ保護として機能する。

---

### `src/kinetics.py`

---

#### 【修正 K-1】K_eq 単位のドキュメント誤記修正

**問題**  
`_r1` および `calculate` の docstring に `K_eq [Pa²]` と記載されていたが、
速度式 `P_A - P_B * P_C / K_eq` を単位整合させると `K_eq [Pa]` が正しい
（`thermo.calc_keq` も `[Pa]` を返す）。

**修正内容**  
該当箇所のドキュメントを `[Pa²]` → `[Pa]` に修正。機能への影響なし。

---

## 要判断事項（ユーザーが決定すること）

### J-1: `_U0_MAX = 2.0 m/s` の妥当性

現在 `_U0_MAX = 2.0 m/s` を上限として設定しているが、これは **工業用 PSA の典型的な
操作範囲（0.1〜1 m/s）の上限に余裕を持たせた仮置き値** である。

- 最適化で `D_col` の下限制約を別途設けている場合は `_U0_MAX` との整合性を確認すること
- `u_0 > _U0_MAX` の領域をペナルティ（CAPEX = 1e9）で埋めると勾配ベースのソルバーが
  崖付近で不安定になる可能性がある。`D_col` に直接下限制約を課す方が数値的に安定

**場所**: `psa_system.py` の `_U0_MAX` 定数

---

### J-2: `bubble_point_T` / `dew_point_T` が `nan` を返す場合の扱い

`eos.py` の修正で、brentq 収束失敗時に `nan` を返すようにした。
**これらの関数を呼んでいる箇所（膜分離システム等）でのエラーハンドリングを追加する必要がある。**

現状では `nan` を受け取った計算がそのまま伝播し、コスト関数が `nan` になる可能性がある。
対策案:
- 呼び出し側で `math.isnan(T_bp)` チェック → ペナルティ返却
- または brentq の探索範囲 `[T_lo, T_hi]` を PSA/膜分離の実際の操作範囲に合わせて
  絞り込む（現状 `[150K, 500K]` は C3 系専用）

---

### J-3: `_T_ABS_MIN = 1.0 s` の妥当性

`desorption_target` が高い（≥ 0.99 など）ケースでペナルティを返す閾値として
1 秒を設定した。最適化で `desorption_target` が実際に 0.9 台を探索するなら
より高い閾値（例: 60 s）に引き上げて「物理的に意味のない超短吸着サイクル」を
早期に排除する方が効率的。

**場所**: `psa_system.py` の `_T_ABS_MIN` 定数

---

## 要調査事項

### R-1: `compress_isentropic` フォールバック `kappa_approx = 1.13` の適用範囲

PSA フィードは H2 リッチ（C: H2 が主成分）であり、C3 混合を想定した
`kappa_approx = 1.13` は不適切な可能性がある。`compress_isentropic` が
PSA 上流の圧縮機計算にも使われるかを確認し、H2 リッチ組成での `kappa` を
別途計算するか、フォールバック自体が実際に発火するケースを把握すること。

### R-2: `bubble_point_T` / `dew_point_T` の探索範囲 `[150K, 500K]`

膜分離システムで使われている場合、操作温度・圧力によっては泡点/露点が
この範囲外に出る可能性がある。実際の操作条件から範囲が十分かを確認すること。

### R-3: CSS 補正と `q_final` の整合性（影響評価）

`q_final` は `t_abs_clean` 時点の固相負荷量（清浄床シミュレーションの終点）。
CSS 補正で `t_abs < t_abs_clean` に短縮しても `q_final` は修正されていない。
`t_abs_clean > t_abs` なので `q_final` は実際の CSS 状態より多くの吸着量を含み、
`t_des` を保守的（やや長め）に推算している。この誤差の大きさは
`t_des / t_abs_clean` 比によって変わるため、代表ケースで影響度を確認すること。

### R-4: H2 損失パージモデルの精度（既存の TODO を再確認）

`H2_loss_purge` の計算はパージ速度を `u_0` で固定しており保守的な過大推算。
物質収支クランプ（修正 P-5）を入れたことで最悪ケースが回避されたが、
`t_des/t_abs >> 1` のケースでは依然として H2 損失が現実より大きく出る。
最適化で H2 回収率が目的関数に入っている場合、この過大推算がどの程度
最適解に影響するかを確認すること。

---

## 修正対象外（既存設計の意図的な仮置き）

以下は今回の修正スコープ外。MEMORY.md や コード内コメントに既に記録済み。

| 項目 | 場所 | 状態 |
|------|------|------|
| Langmuir パラメータ (q_s, a) | `cost_parameters.py` | ★★★ 仮置き |
| KFa 物質移動係数 | `cost_parameters.py` | ★★★ 仮置き |
| 活性炭単価 | `cost_parameters.py` | ★★★ 仮置き |
| CSS 近似精度 (scaling_ratio < 10 警告) | `psa_system.py` | 既存警告で対応済み |
| 脱着モデル `q_final` 空間非一様性 | `psa_system.py` | 安全係数 1.2 でカバー済み |
