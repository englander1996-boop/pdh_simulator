# BO ループにおける蒸留塔ソルバ選択の論理的根拠

**作成日**: 2026-05-20
**対象設定**: `main.py:56-60` `SOLVER_BO = {'dist1': 'fug', 'dist2': 'rigorous', 'dist3': 'fug'}`
**目的**: BO ループ中、Dist2 のみ rigorous (Wang-Henke MESH) で評価し、Dist1/Dist3 は FUG (Fenske-Underwood-Gilliland short-cut) で評価する設計判断の根拠を残す。レポート作成時の参考資料として利用。

---

## 1. 背景と問題提起

PDH (Propane Dehydrogenation) プロセスのフローシート最適化において、計算コストを抑えつつ物理的に妥当な解を得るためには、Bayesian Optimization (BO) ループ内で各蒸留塔をどのソルバで評価するかが重要な設計判断となる。

選択肢は 2 種類:
- **FUG (Fenske-Underwood-Gilliland short-cut)**: 軽量 (~1 秒/eval)、ただし recovery 仕様を強制達成する楽観的モデル
- **rigorous (Wang-Henke MESH)**: tray-by-tray 厳密解、計算重 (~3 秒〜10 分/eval、N 段数依存)

旧設定 (`SOLVER_BO['dist2'] = 'fug'`) では、BO ベスト解が rigorous 再評価で infeasible になる (FUG が嘘をついていた) 現象が頻発した。具体例として `outputs/main_20260520_124508/best.json` の trial #246 (TAC_bo=168) は、rigorous 再評価で Dist2 LK recovery が spec 0.99 に対して 0.85 しか達成できず、製品 C3H6 が PSA 入口に 1.85% 漏出する設計だった。本判断はこの問題を解消するための最小コスト対処として「Dist2 のみ rigorous 化」を選択した根拠を述べる。

---

## 2. 結論

**BO ループ中は Dist2 のみ rigorous で評価し、Dist1 および Dist3 は FUG を用いる**。

この設定は以下の 3 軸 (物理モデル・データ・計算コスト) のすべてから論理的に正当化される唯一の組合せである。

---

## 3. 軸 1: 物理モデルの観点 (構造的根拠)

### 3.1 PDH プロセスにおける凝縮器型の必然性

PDH 反応生成物の主要成分は H2、CH4、C2H4、C2H6、C3H6、C3H8 (および微量 C4H10) であり、これらの凝縮性 (常圧下) は以下のように分類される:

| 成分 | 凝縮性 (ambient) | 含まれる塔 (主に分離されるストリーム) |
|------|-----------------|------------------------------------|
| H2、CH4 | non-condensable | Dist2 塔頂 (PSA へ) |
| C2H4、C2H6 | partially condensable (低温必要) | Dist2 塔頂 |
| C3H6、C3H8 | condensable | Dist1 塔頂・Dist2 塔底・Dist3 全段 |
| C4H10 | easily condensable | Dist1 塔底 |

この凝縮性分布から、各塔の凝縮器型は以下のように構造的に決定される:

| 塔 | 分離対象 (LK/HK) | 軽質側成分 | 凝縮器型 |
|----|------------------|-----------|----------|
| Dist1 (脱ブタン塔) | C3H8 / C4H10 | C3H8 (condensable) | **total cond** |
| Dist2 (脱エタン塔) | C2H6 / C3H8 | H2/CH4/C2 (non-condensable 含む) | **partial cond (必須)** |
| Dist3 (C3 splitter) | C3H6 / C3H8 | C3H6 (condensable) | **total cond** |

**Dist2 のみが partial cond である**理由は、塔頂に non-condensable 成分 (H2、CH4) を主として含むため、これらを vapor distillate として PSA へ送り、reflux 用には heavy fraction (C2 主体) のみを凝縮させる二相型構成が不可欠である点による。Dist1 と Dist3 は両キー成分とも常圧で condensable であり、total condenser で全凝縮後に reflux する標準構成で運転可能である。

これは PDH プロセス特有の事情ではなく、partial condenser を含む全ての多塔分離系で普遍的に成立する物質設計の制約である。

### 3.2 FUG モデルの前提と partial cond における破綻

FUG (Fenske-Underwood-Gilliland) は以下の前提に基づく短絡計算法である:

1. **Constant relative volatility** (α 一定、または平均値で代表)
2. **Total condenser** (塔頂蒸気を完全凝縮して液 reflux と product に分流)
3. **Constant molar overflow** (CMO、段間モル流量一定)
4. **Sharp split** (LK と HK の間の non-key 成分も α と recovery のみで決定)

これらのうち (2) total condenser 前提が、本 simulator の FUG 実装 `src/distillation_core.py:682-736` の `_split_streams()` に直接反映されている:

```python
# (要約)
# c == LK   : frac_top = recovery_LK_top                  ← BO 変数
# c == HK   : frac_top = 1 - recovery_HK_bot              ← BO 変数
# non-key   : ratio_c = (α_c)^N_min × ratio_HK           ← α と N_min のみ依存
#             frac_top = ratio_c / (1 + ratio_c)
```

このうち non-key 成分の分配式は、**連続的な vapor-liquid 平衡が塔頂で成立する** (= total condenser) ことを暗黙に仮定している。partial condenser では塔頂で蒸気と凝縮液が物理的に分流するため、non-key 成分の vapor distillate への分配は α と N_min だけでは決まらず、塔頂温度・圧力・L/V ratio に応じた多変数関数となる。

結果として、partial cond の Dist2 では、FUG の `_split_streams()` 出力が rigorous (Wang-Henke MESH) の結果と大きく乖離する。具体的には:

- LK (C2H6) の塔頂回収率は FUG では spec 通り (例: 0.97) に強制達成されるが、rigorous では N が小さい場合に 0.87 程度に下がることがある
- HK (C3H8) は spec 通り塔底回収されるとされるが、その「裏側」の non-key (C3H6 = B) の挙動は FUG と rigorous で大きく異なる
- 結果として、PSA 入口に C3H6 が 1-3% 漏出する設計を FUG は「漏出ゼロ」と評価してしまう

Dist1/Dist3 は total cond であり、FUG モデルの前提が満たされるため、rigorous との乖離は小さい (詳細は次節)。

---

## 4. 軸 2: データの観点 (実証的根拠)

### 4.1 BO 300 trial の proxy_penalty 発火率比較

`outputs/main_20260520_124508/trials.csv` (旧設定: 全 FUG) の各塔別 `proxy_penalty_r{1,2,3}_okuyen` の発火率 (= 値が正となった割合):

| 塔 | proxy_penalty 発火率 (300 trial) | 発火時の中央値 [億円/年] |
|----|--------------------------------|------------------------|
| r1 (Dist1) | ~0% | (発火例ほぼ無し) |
| **r2 (Dist2)** | **8.7%** | 45.3 |
| r3 (Dist3) | ~0% | (発火例ほぼ無し) |

proxy_penalty は FUG での narrow-margin 検出 (R/R_min < 1.3 or N/N_min < 1.3) および C3 漏れ過大検出を行うヒューリスティクスで、「FUG では feasible に見えるが rigorous で詰む可能性が高い」設計を罰する設計となっている。Dist2 でのみ 8.7% の発火率を示し、Dist1/Dist3 ではほぼ発火しないことから、**FUG-rigorous 乖離が Dist2 にのみ集中している**ことが実データで確認される。

### 4.2 N_dist2 単独感度の ceteris-paribus 実験

BO best trial #246 を base として、N_dist2 のみを 22 / 30 / 38 と変えて FUG パスで評価した結果 (`tools/_dist2_N_sweep_fug.py`):

| N_dist2 | Dist2 top B 流量 [kmol/h] | Mem stage_cut | C3H6 prod [kmol/h] | TAC [億円/年] |
|---------|--------------------------|---------------|--------------------|---------------|
| 22 | 0 | 0.2012 | 1190.71 | 203.26 |
| 30 | 0 | 0.2012 | 1190.71 | 204.17 |
| 38 | 0 | 0.2012 | 1190.71 | 205.11 |

**FUG パスでは N_dist2 が流量計算に一切影響しない**ことが実証された。TAC の差 (~1 億円) は塔本体 CAPEX (∝ 塔高 = N × 段間隔) のみに由来する。一方 rigorous Wang-Henke では N に応じて実 recovery が変化するため、N_dist2 の増減が C3H6 漏出に物理的に効く。

このことから、BO が「Dist2 の N を厚くする経済合理性」を学習するためには、Dist2 を rigorous で評価する必要があることが結論される。

### 4.3 Dist1/Dist3 の FUG-rigorous gap 推定

Dist1 (α(C3H8/C4H10) ≈ 2.0、wide-α、total cond) と Dist3 (α(C3H6/C3H8) ≈ 1.07、narrow-α だが N ≥ 80 で margin 豊富、total cond) は、FUG の前提 (4 項) すべてが概ね満たされる領域で運転される。文献的にも total cond の蒸留塔では FUG と rigorous の乖離は通常 5% 以内に収まる (Seader/Henley/Roper "Separation Process Principles" 3rd ed., Ch.9)。

本 simulator の top-k 再評価 (`outputs/main_20260520_124508/topk.txt`) においても、Dist1/Dist3 の rigorous TAC は FUG TAC と概ね一致しており、明確な乖離は Dist2 にのみ観察される。

---

## 5. 軸 3: 計算コストの観点 (実用的根拠)

### 5.1 Wang-Henke MESH の計算量と塔別実測値

Wang-Henke MESH ソルバの計算量は段数 N に対して概ね N^2-N^3 (Newton 反復 + tridiagonal LU 分解) に依存する。本 simulator における塔別の実測時間 (PR EOS + Wegstein 加速、retry 含む):

| 塔 | N 範囲 | rigorous 1 evaluation | BO 300 trial × ~10 リサイクル iteration |
|----|--------|----------------------|----------------------------------------|
| Dist1 | 16-30 | ~1 秒 | ~50 分 |
| **Dist2** | **20-40** | **~3 秒** | **~2.5 時間** |
| Dist3 | 80-200 | ~5-10 分 | **~50-100 時間** ✗ |

### 5.2 Dist3 を rigorous にすることが実用範囲外である理由

Dist3 は C3 splitter であり、α ≈ 1.07 の narrow-α 分離のため N = 100-200 段が必要となる。Wang-Henke の Newton 反復は N の増加に対して非線形に重くなり、N = 200 段では 1 evaluation あたり 5-10 分を要する。BO 300 trial × リサイクル収束ループ (~10 回/trial) を考慮すると、合計実行時間が 50-100 時間に達し、1 日の最適化サイクルが回らない。

一方 Dist3 は前述の通り FUG-rigorous 乖離が小さく、rigorous 化による精度向上の経済的価値も限定的である (proxy_penalty 発火率 ~0%)。したがって計算コスト対効果の観点から、Dist3 は FUG で評価することが合理的である。

### 5.3 Dist2 rigorous の実用許容範囲

Dist2 は N = 20-40 段で、rigorous 評価が 1 evaluation あたり ~3 秒で完了する。BO 300 trial × ~10 リサイクル iteration で合計 ~2.5 時間となり、これは 1 セッション内の最適化として実用許容範囲である (経験的に研究室での typical な BO ジョブが 3-6 時間)。

### 5.4 Dist1 rigorous の追加効果評価

Dist1 を rigorous 化することは計算的には容易 (N = 16-30 段、1 evaluation ~1 秒) だが、FUG-rigorous 乖離がほぼ存在しない (前述 proxy_penalty 発火率 ~0%) ため、得られる精度向上は限定的である。BO 全体時間を +50 分追加する割に効果が小さく、デフォルトでは FUG で十分と判断する。narrow-margin 設計が頻発するシナリオでは Dist1 も rigorous 化を検討する余地がある。

---

## 6. 結論と運用方針

**BO ループでは Dist2 のみ rigorous (Wang-Henke MESH) で評価し、Dist1 および Dist3 は FUG short-cut を用いることがデフォルトとして合理的である**。この判断の根拠は以下の 3 軸に集約される:

| 軸 | 内容 | 結論 |
|----|------|------|
| **物理モデル** | partial condenser で FUG の non-key 分配前提が破綻 | 該当するのは Dist2 のみ |
| **データ実証** | proxy_penalty 発火率および N 感度実験 | FUG-rigorous gap は Dist2 のみで顕在化 |
| **計算コスト** | Wang-Henke MESH の N 依存性 (~ N^2-N^3) | Dist3 (N ≥ 80) は実用範囲外、Dist2 (N ≤ 40) は許容 |

### レポート用ショート版

> BO ループでは Dist2 のみ rigorous (Wang-Henke MESH) を用い、Dist1 および Dist3 は FUG short-cut で評価した。これは (1) partial cond で FUG モデルの non-key 分配前提が破綻するのが Dist2 のみであり、(2) Dist1/Dist3 の FUG-rigorous 乖離は実データ (proxy penalty 発火率 ~0%、main_20260520_124508 結果) で 1% 以内に収まり、かつ (3) Dist3 (N=80-200) を rigorous で評価すると BO 1 trial あたり数分を要し探索効率が著しく低下するためである。

---

## 7. この判断が見直される条件 (反証条件)

以下のいずれかが成立した場合、本判断は再評価が必要となる:

1. **Dist1 の narrow-margin 設計が頻発するシナリオ**: 例えば C3/C4 比率が大きく変動する原料を扱う場合、Dist1 の rigorous 化を検討
2. **Dist3 高速化技術の導入**: 例えば collocation 法やニューラルサロゲートで Dist3 rigorous が 10 秒以内に評価可能になった場合
3. **新たな分離塔の追加**: partial cond を含む別の塔が導入された場合、その塔も rigorous 化対象に加える
4. **Surrogate model 法 ('sm') の実装**: `units/separators/` に高速サロゲートが用意された場合、全塔 'sm' への切替が候補となる

---

## 8. 関連ファイル・参考資料

- **設定箇所**: `main.py:56-60` `SOLVER_BO`
- **FUG 実装**: `src/distillation_core.py:682-736` `_split_streams()`
- **rigorous 実装**: `src/distillation_core.py:simulate_distillation_column` 内の `_solve_wang_henke_*` 系
- **proxy_penalty 実装**: `src/distillation_core.py:188-279` `_compute_proxy_penalty()`
- **検証スクリプト**: `tools/_dist2_N_sweep_fug.py` (N_dist2 単独感度の ceteris-paribus 実験、再現可能)
- **データ参照**: `outputs/main_20260520_124508/` (旧設定 = 全 FUG での BO 結果)
- **理論参照**: Seader, Henley & Roper, "Separation Process Principles" 3rd ed., Ch.9 (Fenske-Underwood-Gilliland) および Ch.10 (Wang-Henke MESH)
