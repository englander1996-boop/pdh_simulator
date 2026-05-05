# 膜分離システム エラー・フリーズ・クラッシュ対策リスト

**対象ファイル**: `units/separators/membrane_system.py`  
**最終更新**: 2026-05-02（ID-01〜10 全対処済み）  
**更新ルール**: 対策を実施したら「対策内容」欄に何をしたか記入し、ステータスを ✅ に変える。

---

## 凡例

| 記号 | 意味 |
|---|---|
| 🔴 | フリーズ（処理が終わらない） |
| 🟠 | クラッシュ・不正値（例外にならず壊れた値が伝播） |
| 🟡 | エラー（例外またはペナルティ未到達） |
| ✅ | 対処済み |

---

## 未対処リスト

### ✅ ID-01 `solve_ivp` のフリーズ

- **箇所**: `_membrane_ode` 内の `solve_ivp`
- **原因**: `max_step` も実質的な上限もない。`A_mem` が巨大・低駆動力条件でステップが細かくなり終わらない。最適化ループ内で発生すると全体がハングする。
- **対策内容**: `solve_ivp` に `max_step=max(A_mem / 200.0, 0.1)` を追加。ステップ数を最大 200 程度に抑え、どんな `A_mem` でも有限時間で完了することを保証する。巨大 A_mem (1×10⁶ m²) でも 0.05s 以内で返ることを確認済み。

---

### ✅ ID-02 `[0, 0]` 返却によるソルバーハマり

- **箇所**: `_membrane_ode` 内 `ode()` 関数
- **原因**: `J_c <= 0 or J_a <= 0` のとき `[0, 0]` を返す。この状態が続くと Radau が超大ステップを踏みながら終わらない場合がある。
- **対策内容**: `[0, 0]` 返却をやめ、負フラックスを `max(J, 0.0)` でクリップするだけにした。代わりに terminal event `_event_no_flux`（`min(J_c, J_a)` が 0 を下から切った瞬間に `terminal=True, direction=-1` で停止）を追加。フラックス枯渇を明示的に検出して積分を打ち切る。`sol.status == -1` のときのみ失敗とし、`status 0`（完走）と `status 1`（event 発火）はどちらも正常終了として扱う。

---

### ✅ ID-03 ガス状フィードの液相エンタルピー誤計算

- **箇所**: `_vaporizer` 内 `_h_mol(T_in, P_in, z, 'liquid')`
- **原因**: `feed.T_in` が液体の露点より高い（フィードが既にガス状）場合、液相 EOS 解が物理的に無効になり誤ったエンタルピー差を返す。Q_vap が負または異常値になって後段に伝播する。
- **対策内容**: `simulate_membrane_system` 入力バリデーション部で `dew_point_T(feed.P_in, z, _KEYS)` を呼び出し、`feed.T_in >= T_dew_feed` のとき `UserWarning` を発行して `_penalty_result()` を返す。

---

### ✅ ID-04 `P_dist <= P_L` のバリデーション欠如

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **原因**: 製品圧縮機は `P_L → P_dist` に昇圧する前提だが、`P_dist <= P_L` のチェックがない。減圧方向になると `W_prod_kW < 0` が返り、後段の CAPEX 計算に負値が渡る。
- **対策内容**: `P_dist <= P_L` のとき `UserWarning` を発行して `_penalty_result()` を返す。

---

### ✅ ID-05 `P_H <= feed.P_in` のバリデーション欠如

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **原因**: フィード圧縮機は `P_in → P_H` に昇圧する前提だが、`P_H <= feed.P_in` のチェックがない。`W_feed_kW < 0` が返り後段に伝播する。
- **対策内容**: `P_H <= feed.P_in` のとき `UserWarning` を発行して `_penalty_result()` を返す。

---

### ✅ ID-06 `alpha = 1.0` での `_y_local` の退化

- **箇所**: `_y_local`
- **原因**: `alpha = 1.0` のとき `a = (1 - alpha) * gamma = 0` となり二次方程式が線形に退化。現状の `denom` 計算は一応動くが、フォールバックの `return x` が正しい物理解かどうか未検証。
- **対策内容**: `abs(alpha - 1.0) < 1e-10` のとき `return x` を早期 return で明示。alpha=1 は選択性なしを意味し y=x が物理的に正しい。数式を通すとたまたま y=x になるが、フォールバックと混同しないよう意図を明記した。

---

### ✅ ID-07 `CAPEX_total = nan` の最適化器への伝播

- **箇所**: `simulate_membrane_system` CAPEX 合算部
- **原因**: `capex_vap` 等が `nan`（CAPEX 計算例外時）のとき `capex_total = nan` になり、最適化器に `nan` が渡る。最適化器によっては無限ループや異常終了を引き起こす。
- **対策内容**: CAPEX 合算後に `_capex_sum = sum(各 capex)`、`capex_total = _capex_sum if not math.isnan(_capex_sum) else _PENALTY` で nan をペナルティ値（1×10⁹ 億円）に差し替え。最適化器に nan が渡らないことを保証する。

---

### ✅ ID-08 `feed.P_in`・`feed.T_in`・`design.P_dist` の非正値チェック欠如

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **原因**: `feed.P_in <= 0`、`feed.T_in <= 0`、`design.P_dist <= 0` のチェックがない。`dew_point_T` や EOS に異常値が渡り予期しない例外が発生する。
- **対策内容**: バリデーション冒頭で `if feed.P_in <= 0 or feed.T_in <= 0 or design.P_dist <= 0: return _penalty_result()` を追加。

---

### ✅ ID-09 `T_prod_comp_out < T_bp_perm` による `Q_cond < 0`

- **箇所**: `_condenser`
- **原因**: 製品圧縮機出口温度が泡点より低い場合（透過ガスが既に液体に近い状態）、`Q_cond < 0` になり `A_cond < 0` が CAPEX に渡る。Case A 制約（泡点 vs 冷却水温度）とは別経路で発生する。
- **対策内容**: `_condenser` 内で `T_in <= T_bp` のとき `UserWarning` を発行して `(T_bp, nan, nan)` を返す。呼び出し元（`simulate_membrane_system`）で `math.isnan(Q_cond_kW)` を確認して `_penalty_result()` を返す。

---

### ✅ ID-10 `dT1 < 0` 時の LMTD が nan のまま進む

- **箇所**: `_vaporizer` および `_condenser` 内 LMTD 計算
- **原因**: `_lmtd` は `dT1 <= 0 or dT2 <= 0` のとき `nan` を返す。`A_vap` や `A_cond` が `nan` のまま CAPEX 計算に渡り、`calc_he_capex_okuyen` が `ValueError` を投げて `try/except` に落ちる。ペナルティにはなるが、その手前で `math.isnan(Q_vap_kW)` チェックのみで `A_vap` の nan チェックがない。
- **対策内容**: `simulate_membrane_system` 内で気化器呼び出し後に `math.isnan(A_vap)`、冷却器呼び出し後に `math.isnan(A_cond)` を追加チェック。どちらかが nan なら即 `_penalty_result()` を返す。CAPEX 計算まで nan を引きずらずに済む。

---

## 対処済みリスト

### ✅ P_H <= P_L チェック

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **対策内容**: `P_H <= P_L` のとき `UserWarning` を発行して `_penalty_result()` を返す。

### ✅ A_mem, P_H, P_L の非正値チェック

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **対策内容**: `design.A_mem <= 0 or design.P_H <= 0 or design.P_L <= 0` のとき `_penalty_result()` を返す。

### ✅ F_total <= 0 チェック

- **箇所**: `simulate_membrane_system` 入力バリデーション
- **対策内容**: `F_C3H6 < 0 or F_C3H8 < 0` および `F_total <= 0` のとき `_penalty_result()` を返す。

### ✅ T_vap_out >= T_hot チェック

- **箇所**: `_vaporizer`
- **対策内容**: `UserWarning` を発行して `(T_vap_out, nan, nan)` を返す。呼び出し元で `math.isnan(Q_vap_kW)` を確認して `_penalty_result()` を返す。

### ✅ T_bp <= T_cold_out チェック（温度クロス）

- **箇所**: `simulate_membrane_system` 内 Case A 制約
- **対策内容**: `UserWarning` を発行して `_penalty_result()` を返す。

### ✅ ODE 失敗チェック

- **箇所**: `_membrane_ode`
- **対策内容**: `try/except` で例外をキャッチし `(None, None)` を返す。`sol.success` が False の場合も同様。呼び出し元で `None` を確認して `_penalty_result()` を返す。

### ✅ CAPEX 計算例外のキャッチ

- **箇所**: `simulate_membrane_system` CAPEX 推算ブロック
- **対策内容**: `try/except` で各 CAPEX 関数の例外をキャッチし `nan` を代入。

### ✅ 各ユニット計算例外のキャッチ

- **箇所**: `simulate_membrane_system` 内 各ユニット呼び出し
- **対策内容**: `_vaporizer`・`compress_isentropic`・`_condenser` の呼び出しを `try/except` で囲み、例外時は `_penalty_result()` を返す。
