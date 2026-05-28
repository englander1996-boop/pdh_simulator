# 各ユニット(蒸留塔以外)の数学・物理・化学的妥当性 精緻検証

**作成日**: 2026-05-28
**対象**: 蒸留塔を除く全ユニットと src 実装
  `units/reactors/swing.py`, `src/kinetics.py`, `src/catalyst_model.py`, `src/thermo.py`,
  `src/eos.py`, `units/utils/{compressor,pump,cooler,expansion_valve,mixer}.py`,
  `units/separators/psa/psa_system.py`, `units/separators/membrane/membrane_system.py`,
  `src/cost_calculator.py`
**方法**: 各式を支配方程式・熱力学関係・標準相関に照合し、次元整合・符号・保存則・近似の妥当性を確認。
**凡例**: ✓ 妥当 / ⚠ 限界(文書化済・許容) / ✗ 要修正(バグ・内部不整合)

---

## 総括

| ユニット | 判定 | 要点 |
|---|---|---|
| 反応器(化学量論・速度式・エネルギー収支) | ✓ | 原子収支・LHHW 速度式・断熱 dT/dz すべて整合 |
| 反応器(SV / 並列・触媒量の幾何) | ✓ | **実行検証で自己整合を確認**(§1.5、当初の不整合指摘は撤回) |
| 反応速度定数(Arrhenius)・触媒活性 | ✓ | 参照温度形 Arrhenius は標準形と等価。活性スプラインは運転域内 |
| 熱力学(Cp, ΔH, K_eq) | ✓ | Kirchhoff+Gibbs-Helmholtz、K_eq[Pa]=P_STD·exp(−ΔG/RT) 厳密 |
| Peng-Robinson EOS | ✓ | Ω, κ, cubic, φ, 残差 H/S すべて教科書と一致 |
| 圧縮機(`compressor.py`, polytropic) | ✗ | **軸動力 = polytropic head のまま。η_p で割っておらず ~25% 過小** |
| 圧縮機(`eos.compress_isentropic`) | ✓ | W_actual = W_isen/η と正しい(膜用)。← 上と非整合 |
| ポンプ | ✓ | V·dP/η、非圧縮で T 不変、正しい |
| 冷却器/加熱器 | ⚠ | ΔT_lm=30K 固定の placeholder + 冷却時の凝縮潜熱を無視 |
| 膨張弁(JT) | ✓⚠ | PR-EOS 等エンタルピー解は厳密。単相気体仮定で T 低下を過大評価 |
| 混合器 | ⚠ | エンタルピー収支だが**定数 cp**(反応器の多項式 cp と不整合の近似) |
| PSA | ✓⚠ | PDE/LDF/Langmuir/脱着は正しい。均圧なし2塔で H2 回収を保守的過小 |
| 膜分離 | ✓ | GPU→SI 換算・クロスフロー透過 ODE・局所組成2次式すべて正しい |
| コスト(Turton Bare Module) | ✓ | Cp0, Fp, FBM, CEPCI, C_TM すべて標準 |

**最重要(1件)**: 圧縮機 `compressor.py` の η_p 欠落(主要2圧縮機の動力・電力 OPEX・CAPEX を ~25% 過小)。要修正。
**撤回事項**: 初版で「反応器 SV の幾何不整合(真SV 4.77 m/s)」を ✗ としたが、実コード実行検証(§1.5)で**完全に自己整合(SV=0.80 m/s が正しい)と確認、撤回した**。当初の指摘は出力レポートの数値からの手計算 + 並列分割の誤解釈によるもので、コードにバグはない。

---

## 1. 反応器 `units/reactors/swing.py` + `src/kinetics.py` + `src/catalyst_model.py`

### 1.1 化学量論・物質収支 ✓
3 反応 r1: C3H8→C3H6+H2、r2: C3H8→C2H4+CH4、r3: C2H4+H2→C2H6。`_STOICH`(swing.py:93)は全反応で C・H 原子収支が閉じている。`dFdz = (1-eps)·A_cross·(STOICH @ rates)` は次元 [mol/(m·s)] で整合(rate [mol/m³/s] × 触媒断面 [m²])。

### 1.2 反応エンタルピー・エネルギー収支 ✓
`_reaction_enthalpies`(swing.py:315): H_k = ΔHf_298 + ∫_{298}^{T}Cp dT で各成分の絶対エンタルピーを取り、ΔH_j = Σν H を計算 = 標準的な ΔH_rxn(T)。断熱式 `dT/dz = −Q_rxn/ΣF_i Cp_i`(swing.py:392)は、吸熱(ΔH1>0、PDH ≈ +124 kJ/mol)で dT/dz<0(温度降下)を与え符号正しい。

### 1.3 速度式(LHHW)・速度定数 ✓
- r1 = a·k1·(P_A − P_B P_C/K_eq)/(1 + P_B/K_B):プロピレン吸着阻害 + 逆反応(平衡接近)を持つ Langmuir-Hinshelwood 形。駆動力 P_A − P_B P_C/K_eq の単位は Pa(K_eq[Pa])で整合、吸着項は無次元。
- Arrhenius は参照温度形 `k = k0·exp(−Ea/R·(1/T − 1/T0))`(kinetics.py:60)。標準形 k0'·exp(−Ea/RT)(k0'=k0·exp(Ea/RT0))と数学的に等価。✓
- 次元: k1,k2 [mol/m³/s/Pa]、k3 [mol/m³/s/Pa²] で r=[mol/m³/s] に整合。
- a は r1 のみ(クラッキング r2・水素化 r3 には掛けない)= コンテスト要項 §3-2 準拠。

### 1.4 触媒活性 ✓
`calculate_activity_a`(catalyst_model.py)は要項 §3-3 の 7×7 グリッド(T 400–700℃, t 0–30min)を 3 次 RectBivariateSpline で補間、入力 clip(外挿禁止)+ 出力 [0,1] clip。運転域(T_in 880–940K = 607–667℃、t_cyc 12–25min)はグリッド内。a は T_in(= 床内最高温度)で一律決定 = 要項準拠の簡略化。✓

### 1.5 ✓ SV / 並列基数 / 触媒量の幾何 — 自己整合(実行検証済、当初の ✗ 指摘を撤回)
**経緯**: 初版で「SV 計算 `Q_vol/(A_cross·N_parallel)` が触媒量/CAPEX と矛盾、真SV 4.77 m/s」と指摘したが、これは出力レポートの数値からの手計算 + 並列分割の**誤解釈**だった。`simulate_swing_reactor_system` を best #171 設計で実行して検証した結果、**完全に自己整合**で**バグなし**。

**正しい解釈**(検証で確定):
- `A_cross = π/4·D²` は **1塔あたりの断面積**。`D` は 1塔あたり径。全流量 Q_vol は N_parallel 基に均等分配される。
- **SV(1塔あたり)= (Q_vol/N_parallel)/A_cross = Q_vol/(A_cross·N_parallel)** = swing.py:541 と一致。私は「全流量を1塔断面に流す」と誤り Q/A を計算していた。
- ODE は「全流量 × 断面 A_cross × 長さ z_cat」の等価単一 PFR。1塔(流量 Q/N・触媒 V_cat/N・長さ z_cat/N)と**空間時間 V_cat/F が一致**するため転化率が保存される(数学的に妥当な等価変換)。

**実行検証結果(best #171: D=8.336, z_cat=37.576, T_in=939.674, F=6005 kmol/h, P=0.5bar)**:
| 量 | 値 |
|---|---|
| `penalty_reason` | `''`(SV 制約通過 = feasible)|
| N_parallel | 6 |
| コード SV | **0.796 m/s**(∈[0.5,3.0] ✓)|
| 物理 (Q/N)/A = 43.4/54.6 | **0.796 m/s**(一致)|
| 1塔: V_cat / 高さ | 171 m³ / 6.26 m(z_cat=37.6 = 6×6.26)|
| 空間時間 τ: 等価PFR vs 1塔 | 3.93 s = 3.93 s(転化率保存)|
| 転化率 / 選択率 | 33.1% / 77.5%(report と一致)|

→ **SV=0.80 m/s が物理的に正しい1塔あたり空塔速度。幾何の不整合もバグもない。SV 制約は正常に機能している。**(なお、サブ大気圧 0.5bar・長さ 37.6m 床での Ergun 圧損無視は docstring に既知限界として記載済 ⚠ — これは別件の文書化済み限界。)

---

## 2. 熱力学 `src/thermo.py` ✓
- Cp(T)=a+bT+cT²+dT³、ΔH=∫Cp dT を解析積分(thermo.py:172)= 正しい。
- `calc_keq`(thermo.py:179): ΔS298=(ΔH−ΔG)/T0、Kirchhoff で ΔH(T)、∫(ΔCp/T)dT で ΔS(T)、ΔG=ΔH−TΔS、**K_eq[Pa]=P_STD·exp(−ΔG/RT)**。反応 Δn_gas=+1 の標準状態換算(P_B P_C/P_A = Kp·P°)が正しく入っており、速度式の駆動力 P_A−P_B P_C/K_eq と単位整合。厳密。

## 3. Peng-Robinson EOS `src/eos.py` ✓
- Ω_a=0.45724, Ω_b=0.07780(PR1976)、κ=0.37464+1.54226ω−0.26992ω²、a_c=Ω_a R²Tc²/Pc、b=Ω_b RTc/Pc、d√a/dT=−κ√a_c/(2√(T·Tc)) すべて標準。
- 混合則 vdW(kij=0、文書化済 ⚠)、cubic Z³−(1−B)Z²+(A−3B²−2B)Z−(AB−B²−B³)=0、φ_i・残差 H/S(departure functions)すべて教科書(Smith-Van Ness / Poling)と一致。
- bubble_point は外部 `thermo` ライブラリ PRMIX で偽根回避(自作版の単相⇄二相遷移での偽根問題を根治)。dew_point は自作 EOS で逐次置換。✓

## 4. 圧縮機

### 4.1 ✗ `units/utils/compressor.py`(polytropic、Comp2a/Comp2b が使用)— 要修正
- ポリトロピック指数 `n = γ/(γ−(γ−1)/η_p)` → (n−1)/n = (γ−1)/(γη_p) は標準(検算済)。T_out = T_in·ratio^((n−1)/n) も正しい(実 T_out、非効率込み)。
- **問題**: 軸動力 `W = n/(n−1)·R·T_in·[ratio^((n−1)/n)−1]`(compressor.py:90)は **polytropic head**。断熱 1st 則では軸動力 = ṅ·Cp·(T_out−T_in) = γ/(γ−1)·R·T_in·[…]。n/(n−1) = η_p·γ/(γ−1) なので **W = η_p × 実エンタルピー上昇**。
- **実行検証**(Comp2a 相当, 0.5→1.9bar, γ_mix=1.179): `W_code/[ṅ·Cp_γ·(T_out−T_in)] = 0.7500` = ちょうど η_p。→ モデル内エンタルピー上昇 10906 kW に対し W_code=8180 kW。GPSA 標準では Power = polytropic head / η_p。
- **影響**: 真の軸動力は報告値の **1/0.75 = 1.33 倍**。`run_one_pass.py:874,877` の Comp2a/Comp2b(best #171 で各 7996kW)+ Mem 圧縮機の電力 OPEX(~33 億)と圧縮機 CAPEX(~29 億)が ~25% 過小。`T_out` は正しいが `W` と整合しない(エネルギー収支非整合)。
- **副次(別件・軽微)**: γ 表(`_GAMMA`)は C3 等多原子ガスで粗く、Cp_implied=54.8 J/mol/K と実 Cp(~85)より低い(定数γ近似は docstring 既知)。η_p 欠落とは独立の精度課題。

### 4.2 ✓ `src/eos.compress_isentropic`(膜の圧縮機が使用)
等エントロピー T2s → W_isen → **W_actual = W_isen/η**(eos.py:513)→ エンタルピー収支で T2_actual。η で正しく割っており実気体補正も厳密。→ **同一フローシート内で 4.1 と非整合な2モデルが併存**。統一を推奨。

## 5. ポンプ `units/utils/pump.py` ✓
`W = V_dot·dP/η_pump/1000`(pump.py:73)= 流体動力/効率。非圧縮で T 不変。次元・効率の扱いとも正しい(圧縮機と対照的に η で割っている)。

## 6. 冷却器/加熱器 `units/utils/cooler.py` ⚠
- 顕熱 Q=ΣF·cp·ΔT ✓、相変化潜熱は `LATENT_HEAT` 加算。U は contest §4-4 lookup ✓。
- **限界1**: 潜熱は**加熱方向のみ**(cooler.py:114)。冷却時の凝縮潜熱を無視(コメントに「本フローでは不要」)。反応器後 Cooler(523→47℃)は 0.5bar で C3 凝縮しないため実害なしだが、用途を変えると過小評価。
- **限界2**: `dT_lm=30K` 固定の placeholder(cooler.py:84「出典未確認」)。A=Q/(U·ΔT_lm) が線形に効くため、大型 HE(Cooler A=17734 m² 等)の CAPEX が ΔT_lm 仮定に比例して感度を持つ。要文献化。

## 7. 膨張弁 `units/utils/expansion_valve.py` ✓⚠
PR-EOS で等エンタルピー H(T2,P_out)=H(T1,P1) を brentq 求解(`_dh_ig + (Hr2−Hr1)`)= JT 効果の厳密計算。✓ ただし**単相気体仮定**(二相フラッシュなし、docstring 明記)で、飽和近い流体(17bar dist1 塔頂)では T 低下を過大評価する側に偏る ⚠。

## 8. 混合器 `units/utils/mixer.py` ⚠
T_out = Σ(F_ik·cp_k·T_i)/Σ(F_out_k·cp_k)(mixer.py:60)= エンタルピー収支。ただし `cp_of(k)` は**温度非依存の定数 cp** を使用 → 反応器/EOS の多項式 cp と不整合。高温リサイクル(~600℃)混合で T_out に数%級の近似誤差。実装上は許容範囲だが厳密性では劣る。

## 9. PSA `units/separators/psa/psa_system.py` ✓⚠
- ガス相収支 ε∂C/∂t + u₀∂C/∂z + ρ_b∂q/∂t = 0(1次風上差分、interstitial u₀/ε)✓、LDF ∂q/∂t=KFa(q*−q)(Glueckauf)✓、多成分 Langmuir(Markham-Benton)✓、指数脱着 q=q₀exp(−KFa·t) ✓、塔数 ceil(t_des/t_abs)+1 ✓。u₀=F·ZRT/(P·A) 次元整合 ✓。H2 損失(ブローダウン+パージ)は物質収支保護付き ✓。
- **限界**: ① 1次風上の数値拡散(N_z=20、文書化)。② 均圧ステップなし・2塔最小 → H2 回収 74%(実機 85-90%)を**保守的に過小** → 高価な H2 製品(400円/kg)を低価値オフガス燃料へ移し revenue 過小。③ Langmuir/KFa が !仮置き。④ CSS 近似(t_abs×(1−desorption_target))。いずれも文書化済。

## 10. 膜分離 `units/separators/membrane/membrane_system.py` ✓
- **1 GPU = 3.35e-10 mol/(m²·s·Pa)**(membrane_system.py:111)を手計算で検証(1e-6 cm³STP/cm²/s/cmHg → 3.346e-10)= 正確。
- クロスフロー透過 ODE `dF_i/dA = −Q_i·(x_i P_H − y_i P_L)`(分圧差駆動)✓。局所透過組成の2次式 `(1−α)γy² + [(α−1)(x+γ)+1]y − αx = 0`(γ=P_L/P_H)を独立導出で一致確認 ✓。
- 気化器(露点+過熱、LMTD)、製品冷却器(顕熱/潜熱分離、向流 LMTD、Case A: 泡点>冷却水出口で制約)すべて妥当 ✓。圧縮機は `compress_isentropic`(§4.2、正しい)。膜等温・kij=0 は文書化済の仮定。

## 11. コスト `src/cost_calculator.py` ✓
Turton Bare Module 法: Cp0=10^(K1+K2 log10A+K3(log10A)²)、Fp=(Pg+1)D/(10.71−0.00756(Pg+1))+0.5(真空は 1.25)、FBM=B1+B2 Fp FM、C_TM=1.18·CBM·(CEPCI_cur/base)·N すべて標準。適用範囲(容器 0.3–520 m³、HE/Comp/Pump 別係数)を warning 付きでチェック。✓

---

## 推奨アクション(優先順)
1. **圧縮機 `compressor.py` を η_p で割る(§4.1)**、または膜と同じ `compress_isentropic` に統一。電力 OPEX・圧縮機 CAPEX が ~25% 過小(実行検証済 W/Δh=0.75)。あわせて γ 表の精度(§4.1 副次)も検討。
2. 冷却器 `dT_lm=30K` の文献化(§6)、混合器の T 依存 cp 化(§8)は精度向上の中位課題。
3. PSA 均圧モデル化(§9、別 PR 規模)は H2 回収 → revenue の過小評価を改善するが工数大。

> 注: 本検証は式の妥当性に対する静的レビュー + 該当箇所の実行検証。経済影響の定量は ceteris-paribus sweep で確認してから主張すること。
> **教訓(§1.5)**: コードの挙動主張は出力レポートの数値からの手計算でなく、必ず**実コードを実行**して確認する。初版の SV 指摘はこれを怠って誤った。
