# flowsheet/ — 全系評価層 (設計変数 → TAC / Profit)

**役割**: PDH (プロパン脱水素) プロセスの**統合層**。個々のユニットモデル
(`units/` 配下: 反応器・蒸留塔・PSA・膜) を 1 本のフローシートに連結し、
設計変数 1 組を入力として「全系を物理的に解いて経済指標を返す」までを担う。

具体的には以下を 1 関数 `evaluate(design, config)` に束ねる:

1. **1 パス計算** — fresh 原料 + リサイクル仮定値からユニットを順に連結評価
   (`run_one_pass.py`)。
2. **リサイクル収束** — tear ストリーム (Dist3 塔底 / 膜保留側) を逐次置換 +
   Wegstein 加速で収束 (内側ループ) し、さらに生産量目標へ向けて fresh 流量を
   調整 (外側ループ) (`solver.py`, `solvers/`)。
3. **経済集計** — CAPEX (Bare Module 法) + OPEX (Hasebe 式 (10)) から TAC /
   Revenue / Profit を計算 (`economics.py`)。
4. **熱統合** — ピンチ解析 (Stage 1 ターゲティング) で最少加熱/冷却量を求め、
   tier 別ユーティリティ OPEX に置換 (`heat_integration.py`)。
5. **spec 判定** — 製品純度 (C3H6 / H2) と生産量の達成可否を独立評価
   (`specs.py`)。

最適化器 (BO/DE/GA) はこの層の `evaluate()` を目的関数として呼び、
`effective_TAC` (= TAC + 連続ソフトペナルティ) を最小化する。

---

## ファイル一覧

| ファイル | 役割 (1 行) |
|---|---|
| `runner.py` | 最上位エントリ `evaluate(design, config) → FlowsheetResult`。solver→economics→specs を束ね、2 階層ペナルティで `effective_TAC` を確定。 |
| `solver.py` | リサイクル収束 (内側) + fresh 調整 (外側) の 2 重ループ。失敗は例外でなく状態として返す。 |
| `run_one_pass.py` | リサイクル 1 反復分のユニット連結評価。tear/トレースバイパス/penalty 早期検出/shortfall 抽出を含む。詳細 → `SPEC_run_one_pass.md` |
| `economics.py` | CAPEX/OPEX/Revenue 集計、TAC/Profit/製品単価計算 (Bare Module + Hasebe 式 (10))。詳細 → `SPEC_economics.md` |
| `heat_integration.py` | ピンチ解析 Stage1 ターゲティング (MER/A_total/N_HE/tier 配分) + ストリーム抽出 + HI 後経済再計算。詳細 → `SPEC_heat_integration.md` |
| `specs.py` | 製品仕様 (C3H6 純度 / H2 純度 / 生産量) の compliance 判定と %pt 正規化違反量。 |
| `design.py` | 最適化対象の設計変数バンドル `FlowsheetDesignVars` (swing/psa/mem/dist1/dist2/dist3)。 |
| `solvers/base.py` | tear 加速法の抽象基底 `TearAccelerator` (strategy pattern)。 |
| `solvers/wegstein.py` | Wegstein 加速 (履歴 2 点から最適加速係数 q を自動計算)。推奨デフォルト。 |
| `solvers/successive_substitution.py` | 逐次置換 + アンダーリラックス (ベースライン)。 |
| `solvers/__init__.py` | `make_accelerator(method, config)` ファクトリ。 |

---

## 評価フロー

```
evaluate(design, config)                         [runner.py]
  │
  ├─ solve_flowsheet(design, config)             [solver.py]
  │     │  外側ループ: fresh 流量を生産量目標へ調整 (片側相対収束)
  │     │   └─ run_recycle_convergence(...)      内側ループ: tear 収束
  │     │        └─ run_one_pass(tear, fresh, ...) [run_one_pass.py]
  │     │             Pump1 → Dist1 → (膨張) ─┐
  │     │             recycle(Dist3,Mem) ──── Mix → Reactor
  │     │               → Cooler → Comp2(2段) → Desuper → Dist2
  │     │               → PSA → Membrane → Dist3
  │     │             tear_new (Dist3 塔底, Mem 保留) を返す
  │     │        TearAccelerator.step() で次反復推定 [solvers/]
  │     └─ SolverResult (one_pass, inner/outer status)
  │
  ├─ (a) solver-level 失敗判定 → ハード打ち切り (TAC 連続化)
  │
  ├─ calculate_economics(one_pass)               [economics.py]
  │     CAPEX (Bare Module) + OPEX (Hasebe 式 10) → TAC/Revenue/Profit
  ├─ check_specs(one_pass, config)               [specs.py]
  │     C3H6 純度 / H2 純度 / 生産量 → pass + violation_pp
  │
  ├─ (b) spec 違反 → 連続ソフトペナルティ加算
  ├─ HI (apply_hi=True): pinch_analysis → apply_hi_to_economics [heat_integration.py]
  │
  └─ FlowsheetResult(effective_TAC, is_feasible, ...)
```

### 2 階層ペナルティ構造 (`runner.py`)

- **(a) solver-level 失敗** (PSA/Mem/反応器/塔 CAPEX sentinel, リサイクル暴走,
  未収束, タイムアウト): 数値結果が信頼できないため `effective_TAC` をハード
  打ち切り。shortfall 総和に応じて `silent_base`〜`solver_failure_okuyen` で
  連続化し、TPE が infeasible 内の序列を学習できるようにする。
- **(b) spec 違反** (純度・生産量): 数値は信頼できるので連続ソフトペナルティ:
  `effective_TAC = TAC + n_violations×spec_base + Σviolation%pt×spec_coef`
  (全 spec 違反は %pt スケールに正規化済み)。

> 仮置き値 (運転員年俸、膜寿命、各種単価・ペナルティ 4 値、ソルバ設定など) は
> `!仮置き` / 【確認中】 でコード中にマークされている。詳細は
> リポジトリルートの `KNOWN_PLACEHOLDERS.md` を参照。
