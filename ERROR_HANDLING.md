# エラーハンドリング実装サマリー

実装日: 2026-04-30

## 対処した問題と実装箇所

---

### A. ODEソルバーのフリーズリスク

| ID | 問題 | 対処 |
|----|------|------|
| A1 | `solve_ivp` が病的な条件で極めて多くのステップを踏む | `_simulate_one_time` 内の `solve_ivp` 呼び出しを `try/except` で囲み、例外時は `(None, None)` を返す |
| A2 | `n_time_samples` 回のループが全てスタックする | A1 の保護 + 後述の F/T クリッピング（C1/C2）で発散条件そのものを排除 |

**補足**: Radau ソルバーは陰的法のため、根本的なフリーズには入力バリデーション（D1〜D4）による病的入力の排除が第一の防衛線。

---

### B. 数値演算クラッシュ

| ID | 問題 | 箇所 | 対処 |
|----|------|------|------|
| B1 | `K_eq → 0` 時に `P_B * P_C / K_eq → ∞` で駆動力が発散 | `_ode_axial` | `K_eq = max(rc['K_eq'], 1.0)` でクリップ（1 Pa が事実上ゼロ駆動力） |
| B2 | `K_B → 0` 時に吸着項 `1 + P_B/K_B → ∞` でゼロ除算相当 | `_ode_axial` | `K_B = max(rc['K_B'], 1.0)` でクリップ |
| B3 | `calc_fp` の分母 `10.71 - 0.00756*(Pg+1)` がゼロ（>1400 bar） | `cost_calculator.py` | `denominator ≤ 0` 時はペナルティ値 `10.0` を返す |
| B4 | `calc_cp0` に `A ≤ 0` が渡されると `math.log10` で `ValueError` | `cost_calculator.py` | `A ≤ 0` 時は `ValueError` を raise（呼び出し元の B5 で捕捉） |
| B5 | `V_vessel_actual ≤ 0` が B4 を連鎖トリガー | `swing_reactor_simulator.py` | コスト計算前に `if V_vessel_actual <= 0: return _penalty_result()` |

---

### C. サイレントな数値異常

| ID | 問題 | 箇所 | 対処 |
|----|------|------|------|
| C1 | ODE 数値積分の一時的アンダーシュートでモル流量が負値になる | `_ode_axial` | `F = np.maximum(y[:6], 0.0)` でクリップ |
| C2 | 負の F が分圧計算・速度計算を汚染し発散を引き起こす | `_ode_axial` | C1 の対処で連鎖を遮断 |
| C3 | C1/C2 の余波で Conversion/Selectivity が 0〜100% 範囲外になる | `simulate_swing_reactor_system` | `np.clip(value, 0.0, 100.0)` でクリップ |

**ODE 温度クリッピング（追加）**: `T_local = float(np.clip(y[6], 300.0, 1500.0))` — ODE 積分中の温度が物理的範囲を超えて発散するのを防止。

---

### D. 入力バリデーション不足

`simulate_swing_reactor_system` の冒頭で以下をチェック。条件に引っかかった場合は即座に `_penalty_result()` を返す。

| ID | チェック内容 | 対処 |
|----|------------|------|
| D1 | `design.T_in <= 0` または `feed.T_feed <= 0` | `_penalty_result()` |
| D2 | `design.D <= 0`（既存チェックを拡充） | `_penalty_result()` |
| D3 | `feed.P_in <= 0` | `_penalty_result()` |
| D4 | `feed.F_in` のいずれかの値が負 | `_penalty_result()` |
| D5 | `feed.F_in` の合計流量がゼロ（空フィード） | `_penalty_result()` |

---

### E. コスト計算の保護

`calc_reactor_capex_okuyen` の呼び出し全体を `try/except Exception` で囲み、予期しない例外が発生した場合は `Reactor_CAPEX = _PENALTY_CAPEX`（= 1e9 億円）、`TAC = _PENALTY_CAPEX / DEPRECIATION_YEARS` を返す。

---

## 変更ファイル一覧

| ファイル | 変更内容 |
|--------|---------|
| `src/cost_calculator.py` | B3（Fp 分母ガード）、B4（calc_cp0 の A ≤ 0 バリデーション） |
| `units/reactors/swing_reactor_simulator.py` | A1（solve_ivp try/except）、B1/B2（K_eq・K_B クリップ）、B5（V_vessel ガード）、C1/C2（F クリップ）、C3（Conversion/Selectivity クリップ）、温度クリップ、D1〜D5（入力バリデーション）、E（コスト計算保護） |

---

## ペナルティ値について

無効な計算条件はすべて `_penalty_result()` を返す。

```
Reactor_CAPEX = 1e9  [億円]   ← 最適化アルゴリズムへの「ここは探索しない」シグナル
TAC           = 1e9 / 8 [億円/年]
Conversion    = 0.0 [%]
Selectivity   = 0.0 [%]
```
