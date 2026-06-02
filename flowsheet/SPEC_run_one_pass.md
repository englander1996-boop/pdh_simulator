# SPEC: run_one_pass.py — リサイクル 1 パス計算 (フローシート連結)

**ファイルパス**: `flowsheet/run_one_pass.py`

---

## 1. 目的

リサイクル収束の **1 反復** (= 1 パス) に相当するプロセスシミュレーションを実行する。
tear ストリーム (Dist3 塔底 / 膜保留側) の前反復値と fresh LPG 流量を入力に、
全ユニットを規定の順序で連結評価し、各ユニット結果 + tear の更新値 + 各種診断
フィールドを 1 つの `dict` で返す。外側ソルバ (`flowsheet/solver.py`) がこの関数を
反復呼び出ししてリサイクルを収束させる。

成分キー (ProcessStream 仕様): A=C3H8, B=C3H6, C=H2, D=C2H4, E=CH4, F=C2H6, Z=C4H10。

---

## 2. モデル理論 (ユニット連結順序)

### 2-1. tear ストリーム

リサイクルループを切断する 2 本の tear:

| tear | 由来 | 成分 | 元圧力 |
|---|---|---|---|
| `tear_dist3` | Dist3 (C3 スプリッタ) 塔底 | A, B | `design.dist3.P_col` |
| `tear_mem` | 膜 (Membrane) 保留側 (retentate) | A, B | `design.mem.P_H` |

各 tear は温度 (`T_d3`, `T_mem`) も前反復値として受け取り、反応器入口で合流する。

### 2-2. パイプライン順序

```
fresh LPG (30°C 飽和液, C3H8:C4H10=9:1)
  │
Step1: Pump1 (液送昇圧 → design.dist1.P_col) → Dist1 (脱ブタン塔)
  │   Dist1 塔頂 → JT 膨張弁 (→ P_rx = reactor_inlet_Pa) ──┐
  │   Dist1 塔底 (C4H10 富化) → 燃料クレジット (economics)  │
recycle_dist3 → JT 膨張 (→ P_rx) ───────────────────────────┤
recycle_mem   → JT 膨張 (→ P_rx) ───────────────────────────┤
  │                                                        ▼
       reactor_inlet = mix_streams([dist1_top, recycle_dist3, recycle_mem])
  │
Step2: Swing Reactor (型でディスパッチ: 軸流固定床 or 径方向流 多段)
  │   rx_out = 時間平均 effluent (F_out_avg, T_out_avg, P_out)
  │
Step3: Cooler → Comp2a → Intercool → Comp2b → Desuper → Dist2 (脱エタン塔)
  │   圧縮比 ~17:1 を等圧縮比 √17≈4.12 の 2 段に分割、段間冷却
  │   Comp2b 出口 (~151°C) を desuperheater で dew 直上 (50°C) まで冷却水冷却
  │   Dist2 塔頂 → PSA へ / Dist2 塔底 → 膜へ
  │
Step4: PSA (H2 精製、塔頂ガス受入)
  │   製品 H2 / offgas (燃料クレジット)
  │
Step5: Membrane (C3H6/C3H8 分離)
  │   膜前: r2.bottom.P_in > design.mem.P_H なら JT 膨張で P_H まで減圧
  │   mem_precool で mem_feed_K まで気化・過熱 (顕熱+潜熱、ガスフィード)
  │   product (C3 富化) → Dist3 / retentate → tear_mem_new (recycle)
  │
Step6: Dist3 (C3 スプリッタ)
      塔頂 → C3H6 製品 / 塔底 → tear_dist3_new (recycle)
```

反応器入口・膨張弁後の圧力 `P_rx = config.pressure.reactor_inlet_Pa`
(contest 規定 0.5 bar)。Dist1/recycle 各流は反応器圧力まで **JT (等エンタルピー)
膨張** を経由する。膨張弁本体はコストフリー (配管中の絞り弁) で、温度低下を反映して
反応器入口プレヒート量 `Q_preheat` の見積精度を上げる目的。

### 2-3. 反応器モデルのディスパッチ

`design.swing` の型で分岐 (両者とも `(DesignVars, FeedStream, FixedParams) →
SimulationResult` の同一インターフェース):

- `RadialDesignVars` → 径方向流。`n_beds ≥ 2` (既定 3、env `PDH_RADIAL_N_BEDS`
  で上書き可) のとき `simulate_radial_multibed_reactor_system` (段間再加熱付き
  断熱床直列)、単段なら `simulate_radial_flow_reactor_system`。
- それ以外 (`SwingDesign`) → 軸流固定床 `simulate_swing_reactor_system`。

### 2-4. トレースバイパス (モデル簡略化への数値補正)

PSA / 膜の design モデルは入口に「主要成分のみ」を前提する (PSA: CH4 破過計算、
膜: C3H6/C3H8 二成分透過)。上流 (Dist2 partial cond) から微量不純物が漏れ込むと
簡略モデルが破綻し penalty を返してしまう。対処として orchestration 層で
**閾値未満の微量成分を design 計算から除き、マスバランスを保ったまま下流に直接
ルーティング**する (`_apply_trace_bypass`)。これは「物理装置の追加」ではなく
「シミュレータの数値処理」であり、マスは保たれる。

- 閾値 `_TRACE_BYPASS_FRAC = 0.01` (入口総モル流量の 1%、固定)。
  意図: Dist2 を設計で詰めて物理的に漏れ <1% を達成し、残りを bypass で吸収。
  1% 超過は「設計が悪い」=本物の警告。
- PSA で許容しない成分 `_PSA_TRACE_COMPS = ('A','B')` (C3) → 除去分は offgas に合算。
- 膜で許容しない成分 `_MEM_TRACE_COMPS = ('C','D','E','F')` (非 C3) → retentate に合算。
- 閾値超過分 `max_excess_frac` は `trace_bypass_{psa,mem}_excess` として返し、
  `runner.py` が連続ペナルティに使う (係数 1000 億円/fraction)。

### 2-5. penalty 早期検出と shortfall 連続化

各ユニットが penalty (CAPEX sentinel `≥ 1e8 億円`, `feasible=False`, F=T=P=0 の
ゼロ流) を返すと、下流装置が「P_in=0」等で例外 crash し、solver の penalty 判定に
到達できない。これを避けるため各段で早期検出し、下流を stub
(`_PenaltyResult` / `_PenaltyEquipment`, CAPEX=1e9 sentinel) に差し替えた結果
dict を返す:

- Dist1 (`r1.feasible=False`) → `_build_penalty_after_column('r1', ...)`
- Reactor (`Reactor_CAPEX ≥ 1e8`) → `_build_penalty_one_pass_result(...)`
- Dist2 / Dist3 (`feasible=False`) → `_build_penalty_after_column(...)`
- PSA / Mem は早期 return せず通過させ、solver が CAPEX sentinel で後から判定。

同時に、各 penalty 経路の **連続 shortfall** を計算して結果 dict に積み、
`runner.py` 経由で TPE constraints_func に渡す (どちらへ逃げれば feasible に
近づくか BO が学習できるようにする):

| 関数 | 抽出する shortfall (代表) |
|---|---|
| `_compute_reactor_shortfall` | `reactor_sv_shortfall` (SV 範囲外の log10 距離), `reactor_dp_shortfall` (ΔP/P 超過分), `reactor_other_shortfall` |
| `_compute_psa_shortfall` | `psa_t_abs_shortfall`, `psa_u_0_shortfall`, `psa_dp_shortfall`, `psa_feed_shortfall` |
| `_compute_mem_shortfall` | `mem_ph_shortfall`, `mem_bp_shortfall`, `mem_phase_shortfall`, `mem_other_shortfall` |
| `_compute_dist_shortfalls` | `dist{1,2,3}_N_shortfall` (FUG Gilliland infeasible), `dist{1,2,3}_dT_shortfall` (Wang-Henke 収束失敗の log10 比), `dist{1,2,3}_cond_shortfall` (HYSYS cold-top) |

### 2-6. warning 集約

silent fallback (PR EOS Z=1 fallback, brentq 偽根, Wang-Henke 残差超過の FUG
fallback 等) を見逃さないため、各ユニット呼び出しを `_capture_warnings(source, log)`
(= `catch_warnings(record=True) + simplefilter("always")`) で囲む。捕捉した warning
は `_CapturedWarning` (source ラベル付き) として `warnings_captured` に集約し、
`runner.py` が `FlowsheetResult.warnings_captured` に流す。

### 2-7. 観測ラベル抽出

- `_determine_first_failed_unit` : パイプライン順 (r1→r_rx→r2→r_psa→r_mem→r3) で
  最初に penalty を起こしたユニットを返す (`first_failed_unit`)。success path でも
  PSA/Mem は早期 return しないため、Dist3 まで素通りした trial の「実は PSA で
  死んでいた」を取りこぼさない。
- `_extract_unit_diagnostics` : 各 equipment から `penalty_reason` + key actual 値
  (SV, t_abs, u_0, P_H, T_dew, N_needed, dT_max 等) を 1 dict (`_EMPTY_UNIT_DIAG`
  と同キーセット) に抽出。live 表示 / CSV groupby 用。

---

## 3. 入出力

### 入力 (`run_one_pass(...)`)

| 引数 | 型 / 単位 | 説明 |
|---|---|---|
| `tear_dist3` | dict {A,B} [kmol/h] | Dist3 塔底由来リサイクル (前反復値) |
| `tear_mem` | dict {A,B} [kmol/h] | 膜保留側由来リサイクル (前反復値) |
| `T_d3`, `T_mem` | float [K] | 各リサイクルの温度 (前反復値) |
| `F_C3H8_feed`, `F_C4H10_feed` | float [kmol/h] | Fresh LPG 中 C3H8 / C4H10 流量 |
| `design` | `FlowsheetDesignVars` | 設計変数 (swing/psa/mem/dist1-3) |
| `config` | `OperatingConfig` | 圧力・温度・製品仕様・原料状態 |

### 出力 (`dict`)

- 各ユニット結果: `pump1, r1, dist1_top_rx, reactor_inlet, r_rx, rx_out, cooled,
  comp2a, intercool, comp2b, desuper, r2, r_psa, mem_precool, r_mem, r3`
- tear 更新値: `tear_dist3_new {A,B}`, `tear_mem_new {A,B}`, `T_d3_new`, `T_mem_new`
- 診断: `warnings_captured`, `trace_bypass_{psa,mem}_excess`, 各 shortfall フィールド,
  `first_failed_unit`, および `_EMPTY_UNIT_DIAG` の各観測ラベル

---

## 4. 主要パラメータ

| 定数 | 値 | 意味 |
|---|---|---|
| `_TRACE_BYPASS_FRAC` | 0.01 | トレースバイパス微量判定閾値 (入口総量の 1%) |
| `_PSA_TRACE_COMPS` | (A,B) | PSA で許容しない成分 (C3) |
| `_MEM_TRACE_COMPS` | (C,D,E,F) | 膜で許容しない成分 (非 C3) |
| `_RIG_TOL_K` | 0.05 | Wang-Henke 収束 tol (距離計算用、`src/distillation_rigorous.py` と同期) |
| `_PSA_T_ABS_MIN_S` / `_PSA_U0_MAX_MS` / `_PSA_DP_MAX_BAR` | 60.0 / 1.0 / 0.3 | PSA shortfall 基準 (`psa_system.py` と同期) |
| `_REACTOR_SV_MIN_MS` / `_REACTOR_SV_MAX_MS` / `_REACTOR_DP_OVER_P_MAX` | 0.5 / 3.0 / 0.10 | 反応器 shortfall 基準 (`swing.py` と同期) |
| `_MEM_T_MARGIN_K` | 5.0 | 膜 bp/dew shortfall の健全マージン |
| `T_dist2_feed_K` | 323.15 (50°C) | desuperheater 出口温度 (8.5 bar dew +5K margin) |
| env `PDH_RADIAL_N_BEDS` | 3 | 径方向流多段の段数 |
| env `PDH_PER_UNIT_LOG` | 0 | penalty 発火を stderr に 1 行ログ |
| env `PDH_MEM_QA_FACTOR` / `PDH_MEM_ALPHA_FACTOR` | 1.0 | 膜性能劣化係数 (感度解析用) |

---

## 5. 出典

- ユニット連結圧力・段数: contest §3-3-3 (加圧箇所には液ポンプ)、反応器入口圧 0.5 bar
  (contest 規定)。
- 膜 P_H ≤ 9.5 bar: Hua et al. (2024) (膜検証範囲)。
- 反応器 SV / ΔP / PSA / 蒸留塔の各 penalty 基準値は各ユニット SPEC
  (`SPEC_swing.md`, `SPEC_psa_system.md`, `SPEC_columns.md`,
  `SPEC_membrane_system.md`) と同期。

---

## 6. 既知の限界・仮置き

- **トレースバイパスは数値処理**であり物理装置ではない。1% 超過は penalty 化されるが
  多成分化 (PSA/膜モデルの拡張) は「コスト対効果が悪い」として実施していない。
- 膨張弁・JT 膨張はコストフリー (CAPEX/OPEX 非計上)。
- 反応器段数のデフォルト 3 段は「収率≈選択率で頭打ち、TAC 上 3 段有利」という
  判断 (詳細はレポート反応器章)。env で感度解析可。
- shortfall の log10 / 比スケールは TPE への方向シグナル目的の正規化であり、
  物理量そのものではない。
- penalty 経路で下流装置を stub 化するため、その trial の下流ユニット結果は
  ダミー (ゼロ流 / sentinel CAPEX) になる。
