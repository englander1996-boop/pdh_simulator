# SPEC_bo — ベイズ最適化の目的関数・制約・スケジュール仕様

本書は `optimization/` のベイズ最適化 (BO) 層の理論・式・パラメータを定義する。
コード本体は `objective.py` / `study.py` / `penalty_scale.py` / `topk.py`、
固定パラメータは `config/operating.toml`、目的関数の組み立ては
`flowsheet/runner.py` (= `evaluate`) にある。

---

## 1. 目的

設計変数ベクトル **x** (反応器・PSA・膜・Dist1/2/3・原料流量・各塔 recovery)
に対し、`effective_TAC(x)` を最小化する。

```
minimize_x  effective_TAC(x)        [億円/年]
```

`effective_TAC` は「年次総コスト (TAC) + 制約違反のソフトペナルティ」であり、
最適化器 (Optuna, `direction='minimize'`) が見る唯一のスカラ目的関数。
診断用に `is_feasible` / `failure_reason` / `specs` / 各種 shortfall が
`trial.user_attrs` に併記される。

---

## 2. 理論 (式)

### 2.1 目的関数 = effective_TAC (2 段ソフトペナルティ構造)

`flowsheet/runner.py` が定義する 2 階層構造:

**(a) solver-level 失敗** (PSA/膜 CAPEX のペナルティ発火、リサイクル暴走、
NaN、未収束、未処理例外)。数値結果が信頼できないため評価不能として
大きな固定値で打ち切る:

```
effective_TAC = penalty.solver_failure_okuyen        (= 10000 億円/年)
```

現実 TAC (200–300 億円/年) の数十倍。例外経路は連続化不可のため固定。

**(b) spec 違反** (純度・生産量)。数値結果は信頼できるので、連続ソフト
ペナルティで方向感を与える:

```
effective_TAC = TAC + spec_base + spec_coef × Σ(violation_%pt)
```

- `spec_base` = `penalty.spec_base_okuyen` (= 50 億円/年): 違反が 1 つでもあれば固定加算。
- `spec_coef` = `penalty.spec_coef_okuyen` (= 100 億円/(年·%pt)): 違反量に比例。
- 全 spec 違反は %ポイント (%pt) スケールに正規化済み (`flowsheet/specs.py`)。

採用理由: 化学プロセス最適化は微分不要系 (DE/GA/BO) が実用的で、それらは
feasibility 領域外でもスカラ値を要する。陽な制約よりソフトペナルティの方が
探索効率が良い (設計判断 2026-05-08)。

### 2.2 探索 sampler = TPE + Sobol QMC 前置 (2-phase)

`study.py::make_sampler('tpe', ...)`:

- **Phase 1 (完了 trial < `n_startup`)**: Sobol QMC (`QMCSampler`, 低乖離点列)
  で広域カバレッジ。高次元 (最大 22 変数) で pure random の偏りを避け、
  狭い feasible 領域 (本問題では 1% 以下) を掴む確率を底上げする。
- **Phase 2 (それ以降)**: TPE。`_PhaseSwitchSampler` が完了 trial 数で切替える
  (Optuna 標準 TPE は startup を内部 RandomSampler で生成し Sobol へ差し替える
  正規 API が無いため、薄い delegator wrapper を自作)。

TPE の主な設定:

- `multivariate=True`: Parzen estimator を多変量化し変数間相関を学習。
- `n_ei_candidates=200`: 各 trial で 200 候補から最良 EI を選ぶ (既定 24 → exploration 強化)。
- `constant_liar`: 並列時に pending trial を悲観仮値で埋め、worker の群がりを抑制。

Phase 1 (QMC) 中も `failure_reason` は記録され、Phase 2 突入時に TPE が
Phase 1 を含む全 trial 履歴から学習する。

### 2.3 制約 — `constraints_func` (TPE への連続シグナル)

`study.py::_default_constraints_func(trial)` が trial → `Sequence[float]` を返す。
各要素は **負値 = feasible / 正値 = 違反** と TPE が解釈する。値は
`objective.py::_store_diagnostics` が `user_attrs` に格納したものを読む。

目的: binary な `[proxy, feas_flag]` だけでは異種 infeasibility が潰れて
TPE が「どちらに動かすべきか」を学べない。塔別・装置別の **連続 shortfall** を
渡すことで、infeasible 領域内でも勾配を立てる。

現行は 25 要素 (順序は `study.py` の return を正とする):

| idx | キー | 意味 |
|---|---|---|
| 0 | proxy_penalty_total | rigorous プロキシ罰則合計 (×0.1 正規化, 10億円=1.0) |
| 1 | feas_violation | `is_feasible=False` で 1.0 |
| 2–4 | dist1/2/3_N_shortfall | FUG Gilliland 段数不足比 |
| 5–7 | dist2/1/3_dT_shortfall | rigorous Wang-Henke 収束不足 (log10 比) |
| 8–10 | psa_t_abs / psa_u_0 / psa_feed | PSA 吸着時間・空塔速度・feed 異常 |
| 11–12 | reactor_sv / reactor_other | 反応器 SV 範囲外 (×20)・その他 |
| 13–14 | production_under / over_pp | 生産量下限不足・上限超過 [%pt] (×0.2, 5pp=1.0) |
| 15–18 | mem_ph / mem_bp / mem_phase / mem_other | 膜の各 silent penalty 経路 |
| 19–20 | trace_bypass_psa / mem_excess | 非 C3 漏れ閾値超過 (×100, 1%=1.0) |
| 21 | unknown_failure | `is_feasible=False` かつ全 shortfall=0 の保険 |
| 22 | dist2_cond_shortfall | HYSYS Dist2 凝縮器ΔT不成立 (×0.1) |
| 23 | reactor_dp_shortfall | 反応器 Ergun 圧損 ΔP/P 超過 (×5) |
| 24 | psa_dp_shortfall | PSA 床 Ergun 圧損 上限超過 (×5) |

各 raw 値は「典型的な marginal 違反 ≈ 1.0」に正規化される (係数は上表 / `study.py` 参照)。
QMC phase の trial にも `constraints_func` を実行して `system_attrs['constraints']` に
格納し、TPE 切替後の「constraint 欠落」警告を抑制する (`_PhaseSwitchSampler.after_trial`)。

### 2.4 penalty_scale スケジュール (adaptive penalty)

`penalty_scale.py`。trial 進行に応じてソフトペナルティ係数 (trace_bypass /
proxy / spec) を動的に強化する module-level 変数 `_CURRENT_SCALE`。各 trial 開始時に
`objective.py` が `set_scale(schedule(trial.number, n_total))` で更新し、各 penalty
計算箇所が `get_scale()` を掛ける。

既定スケジュール `default_schedule` (2026-05-23, 線形):

```
scale(t) = 0.2 + (3.0 - 0.2) × (t / n_total)
```

- 序盤 (scale≈0.2): 弱ペナルティで探索を広げ、spec ぎりぎりの trial も
  TPE に「悪くない TAC」として学ばせる。
- 終盤 (scale→3.0): infeasible 領域から強制退出させ feasible へ収束。

`linear_schedule(start, end)` で任意の線形区間も指定可能。
注意: `_CURRENT_SCALE` は process-local global。`n_jobs>1` (スレッド並列) では
レース可能性があり、最終評価は `n_jobs=1` 推奨。マルチ**プロセス**並列
(`parallel.py`) は各プロセスで scale が分離されるため安全。

### 2.5 多忠実度 (multi-fidelity)

蒸留塔ソルバを忠実度別に使い分ける。`ColumnTunables.solver_method`:

| method | 内容 | 用途 |
|---|---|---|
| `fug` | Fenske-Underwood-Gilliland shortcut (高速・決定的) | BO ループ (多数 trial) |
| `rigorous` | VLE ベース tray-by-tray (Wang-Henke MESH, 厳密) | top-k 再評価 |
| `sm` | 簡易厳密モデル (special フォーク) | 中間忠実度 (sub2 等) |
| `hysys` | HYSYS COM バックエンド (段数別 HSC) | 検証・main 経路 |

**2 段最適化**:
- BO ループは全塔 `fug` で高速に広く探索 (`solver_bo`)。
- BO 終了後、Profit 上位 k 候補だけを `solver_topk` (例: 全塔 `rigorous`) で
  再評価 (`topk.py::reevaluate_topk`)。FUG は narrow-margin 設計で楽観的に
  振れるため、rigorous で is_feasibility / TAC を是正し最終解を確定する。
- BO loop には rigorous プロキシ罰則 (`proxy_penalty`) を乗せ、FUG では feasible
  だが rigorous で詰まる領域に張り付くのを抑制する (BO と top-k の乖離を縮める)。

HYSYS は単一 COM のため並列不可 (worker=1)。FUG/rigorous/SM はマルチプロセス並列可。

---

## 3. 入出力

**入力:**
- 設計変数: `SEARCH_SPACE` (`{key: (low, high, scale, type)}`, `search_space.py`)。
  suggest 対象外のキーは `DEFAULT_BASELINE` で補完。許容キーは `EXPECTED_KEYS` (最大 22)。
- 固定運転条件: `OperatingConfig` (`config.load.load_operating_config`)。
- ソルバ割当: `solver_bo` / `solver_topk` (`{'dist1'/'dist2'/'dist3': method}`)。
- ハイパラ: `PipelineConfig` (n_trials, n_startup, n_topk, seed, sampler, n_jobs, n_workers, …)。

**出力 (`outputs/main_<ts>/`):**
- `optuna.db` (SQLite, 再開/ダッシュボード用)
- `trials.csv` (全 trial 履歴: params + user_attrs)
- `best.json` (BO 単体ベスト)
- `topk.txt` (★最終結果: BO vs rigorous 再評価の比較)
- `feasibility.txt` / `feasibility_2d.png` (L1 収束分類解析, 任意)
- `README.md` (結果の見方ガイド)

---

## 4. 主要パラメータ

| パラメータ | 既定 | 出所 | 意味 |
|---|---|---|---|
| `n_trials` | 300 | `PipelineConfig` | BO trial 総数 (penalty schedule の分母) |
| `n_startup` | 50 | `PipelineConfig` | Sobol QMC startup trial 数 |
| `n_topk` | 10 | `PipelineConfig` | rigorous 再評価する上位候補数 |
| `seed` | 42 | `PipelineConfig` | 乱数シード |
| `sampler` | `'tpe'` | `PipelineConfig` | `'tpe'`/`'cmaes'`/`'random'` |
| `n_jobs` | 1 | `PipelineConfig` | スレッド並列数 (penalty_scale 競合のため 1 推奨) |
| `n_workers` | 1 | `PipelineConfig` | マルチプロセス並列 worker 数 (>1 は SQLite 必須) |
| `n_ei_candidates` | 200 | `study.py` | TPE の EI 候補数 (exploration 強化) |
| `solver_failure_okuyen` | 10000 | `operating.toml [penalty]` | solver 失敗時の固定 TAC |
| `spec_base_okuyen` | 50 | `operating.toml [penalty]` | spec 違反の固定加算 |
| `spec_coef_okuyen` | 100 | `operating.toml [penalty]` | spec 違反 %pt 当たりの係数 |
| `c3h6_min_wtfrac` | 0.995 | `operating.toml [spec]` | C3H6 製品純度下限 |
| `h2_min_molfrac` | 0.999 | `operating.toml [spec]` | H2 製品純度下限 |
| `production_min/max_relative` | 0.05 | `operating.toml [spec]` | 生産量許容不足/超過 (両側 ±5%) |
| penalty schedule | 線形 0.2→3.0 | `penalty_scale.py` | adaptive penalty scale |

---

## 5. 出典

- **目的関数 = TAC + 2 段ソフトペナルティ**: 設計判断 (2026-05-08)。
  `flowsheet/runner.py` 冒頭 docstring および `config/operating.toml [penalty]` のコメント。
- **TPE**: J. Bergstra et al., *Algorithms for Hyper-Parameter Optimization*, NeurIPS 2011
  (Optuna `TPESampler` の基礎理論)。本リポジトリでは Optuna 実装をそのまま使用。
- **Sobol QMC startup**: 低乖離点列による高次元初期サンプリング (Optuna `QMCSampler`)。
- **FUG shortcut**: Fenske-Underwood-Gilliland 法 (`src/distillation_core.py`)。
- **製品仕様 (C3H6 99.5 wt% 等)**: コンテスト課題 Ver.2.0 §2-1 を必須制約として採用、
  H2 純度はコンテスト規定外の独自設定 (`operating.toml [spec]` のコメント参照)。
- **penalty schedule / constraints_func / 正規化係数**: 本リポジトリ内の forensic
  (run ログ解析) に基づく経験的調整。日付付きの設計判断コメントを各ソースに記載。

注意: 外部出典の URL/数値の最終確定はエージェント単独では行わない方針
(検索補助は可)。上記コンテスト仕様の引用は MEMORY の contest_spec を根拠とする。

---

## 6. 既知の限界

- **狭い feasible 領域**: 本問題は feasible 領域が極めて狭く (≤1%)、startup random
  だけでは全 infeasible になり TPE が立ち上がれないことがある。Sobol QMC startup +
  warm-start で緩和しているが、根本的に探索効率はこの狭さに律速される。
- **BO ≠ 真のベスト**: BO ループは全塔 FUG の楽観評価。最終解は必ず top-k の
  rigorous 再評価 (`topk.txt` rank 1, `feas_re=True` の最小) で判断すること。
- **penalty_scale の並列安全性**: module-level global。スレッド並列 (`n_jobs>1`) では
  競合し得る。マルチプロセス並列では scale がプロセス分離され影響軽微。
- **constraints_func の正規化係数**: 経験則ベース。スケール (×0.1/×5/×20/×50/×100 等)
  が変わると TPE の制約バランスが変わる。
- **rigorous 再評価のコスト**: Dist3 N≥170 では 1 リサイクル周回が ~10 分に達し、
  内側収束 tol を 1% に緩和し max_iter=100 で打ち切る運用 (`operating.toml [solver.inner]`)。
- **多忠実度の不整合**: SM/HYSYS と FUG/rigorous の間で同一設計の評価値が乖離する
  ケースがある (例: Dist3 CAPEX)。検証は別途課題。
- **HYSYS 経路は並列不可** (単一 COM, worker=1)。
- **暫定値 (`!仮置き` / 【確認中】)**: 物性・経済パラメータ・一部の正規化係数に暫定値が
  含まれる。詳細は リポジトリ直下 `KNOWN_PLACEHOLDERS.md` を参照。
