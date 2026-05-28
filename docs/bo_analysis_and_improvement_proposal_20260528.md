# special.py BO 挙動分析 & 改良案

**作成日**: 2026-05-28
**対象**: `special.py`(HYSYS+SM 全21変数 BO) と BO 基盤(`optimization/study.py`・`objective.py`・`penalty_scale.py`)
**データ**: `outputs/special_20260528_124401`(baseline)、`outputs/special_20260528_150021`(cond-signal 実験)
**目的**: 最新2 run の結果を精査し、BO の改良余地を実データで切り分ける。

---

## 0. 結論(先出し)

**この構成では BO は既に達成可能な最良(effective_TAC ≈ 1130 億円/年)に頑健に到達しており、BO の小細工(sampler・制約・penalty 調整)では best は動かない。** これは `docs/bo_dist2_cond_shortfall_experiment_20260528.md` の結論をデータで多角的に裏付けた。改善余地は **(A) 探索効率**(死領域の除去・pre-screen)と **(B) 設計レバー/モデル**(SM Dist3 の妥当性、収率↔Dist2 cold-top トレードオフ)にあり、後者が本丸。

---

## 1. 現状診断(最新2 run の実測)

| 指標 | 124401 (baseline) | 150021 (cond-signal) |
|---|---|---|
| 完了 / feasible | 300 / 107 (**35.7%**) | 300 / 104 (34.7%) |
| best effective_TAC | **1129.5 億円/年**(#171) | 1131.1(#187) |
| best 利益 | **−276 億円/年(赤字)** | ほぼ同 |
| best 設計性格 | F=1560, prod=1194, 99.5wt%, 収率76.5% | F=1532, 同等 |

### best #171 の経済構造(`top1_trial171.txt`)
- **CAPEX 606 億円**: Dist3 = **334 億(57%)** が断トツ(SM、D=7.36m × H94.2m × N117 のシェル)。反応器 78、Dist2 46、膜 42、Dist1 21。
- **OPEX 1054 億/年**: **原料 LPG = 599 億(49%)** + Hasebe markup 138 が支配。**Dist2 cold-top 冷媒(エチレン −100℃)= 105 億**。反応器予熱 81、触媒交換 40。
- **Revenue 853 億/年**: C3H6 製品 603(70.7%)+ PSA オフガス燃料 134 + Dist1 塔底燃料 74 + H2 製品 42。
- → **原料費(599+markup=737)> 製品売上(603)**。yield 76.5%(選択率77.5%、T_in=940K=探索上限)。**構造的赤字**。

### 失敗内訳(infeasible 193 件 @124401)
| failure_unit | 件数 | 備考 |
|---|---|---|
| **r2 (Dist2/HYSYS)** | **145 (48%)** | うち **86% が cold-top**(塔頂 < エチレン −100℃+ΔT)、14% が COM エラー |
| r3 (Dist3/SM) | 32 | |
| r_rx (反応器) | 9 | SV 範囲外等 |
| spec_production_under | 4 | F_fresh 下限不足はごく少数 |

### 収束タイムライン(best-so-far)
trial 50: 1277 → 100: 1168 → 150: 1147 → **200: 1129.5 → 300: 1129.5(横ばい)**。
**BO は trial 200 で収束、後半 100 trial は改善ゼロ。** feasible のうち best+5 以内は 1 件のみ、best+20 以内 7 件 → best はやや孤立した谷底。
総 wallclock 59 分 / 300 trial(feasible 中央 20.6s、r2 失敗は 6.1s と安価)。

---

## 2. データが否定する「効かない改良」(着手しない)

| 案 | 否定根拠(実データ) |
|---|---|
| cold-top の制約シグナル強化 | 150021 で実証済み無効(feasible 35.7→34.7%、best 悪化)。`bo_dist2_cond_shortfall_experiment` 参照 |
| **col2_p 下限引上げで cold-top 削減** | **否定**: cold-top の col2_p は 506–699 kPa に分布(median 648 ≈ feasible median 667)、**125件中59件が >650kPa**(高圧側)。全圧力域で発生する相互作用効果(リサイクル組成依存)で bounds では潰せない |
| trial 数増 | trial 200 で plateau。basin 探索済み |
| sampler 再調整(TPE↔CMA-ES、n_ei_candidates 等) | 制約構成が異なる2 run(22 vs 23 制約)が共に ~1130 へ収束 → 頑健 basin、sampler 依存ではない |

補足: `study.py:99` の制約 23 本は feasible では全て 0(= 勾配なし)、infeasible 同士の重み付けにしか効かない。Optuna 制約 TPE は本質的に feasible/infeasible の二分であり、48% を占める r2 cold-top はこの仕組みでは回避させられない。

---

## 3. BO 改良案(価値順)

### Tier 1 — 探索効率の確実な改良(best はほぼ不変、反復が速くなる)

**T1. T_in 下限を 880 → 910K に引上げ**(`special.py:112`)
- 根拠: T_in 別 feasible 率 → **(880,900]=0%(17/17 r2 失敗), (900,915]=0%(22/22 r2), (915,925]=14%, (925,940]=49%**。T_in<915K は 39/39 が infeasible で、しかも失敗は **r2 (Dist2 cold-top)**。
- これは col2_p のケースと違い**クリーンな崖**(下は完全に 0%)なので bounds で切ってよい。QMC50+TPE 予算の約13%を生領域に集中できる。
- **重要な訂正**: 探索空間コメント `special.py:107-111` は「低T=production不足」と書くが、**データ上は誤り**(production_under は 4 件のみ)。低 T_in が infeasible になる真因は **Dist2 cold-top**(低 T → 反応器でのクラッキング↓ → Dist2 feed の軽質分↓ → 塔頂組成が変わり partial cond が成立しない)。コメント訂正を推奨。

**T2. r2 cold-top の pre-screen(任意・効率のみ)**
- r2 は総 wallclock の **28%(16分/59分)**。HYSYS 収束後に物理判定で棄却している。dew-point 事前判定で時間回収可。
- ただし best は動かない(plateau 済み)。価値は **Tier 2 のレバー実験を高速に回せる**点のみ。`bo_dist2_cond_shortfall_experiment` の Rec1 と同位置づけ。

### Tier 2 — best/利益を動かしうる構造レバー(BO の外、or BO+model)

**B1.【最優先】SM Dist3 モデルの検証**
- Dist3 が CAPEX の **57%(334 億)**。SM が径 7.36m を返すが rigorous 相当は ~6.8m の疑い(memory `project_dist3_capex_verification`、未着手)。さらに SM の **N≥115 floor とスループット制約が T_in を 940K に固定**している。
- 検証で過大評価が確認されれば: CAPEX 直接 −約110 億(償却 −14 億/yr)+ T_in 解放余地。**単一で最大効果**。

**B2. T_in を Dist2 cold-top 制約から解放**
- 実測: 低T feasible(915–925K)は**収率 77.8% > 高T帯 75.9%**、TAC 中央も 1187 と悪くない。低T = 低クラッキング = 高収率 = feed 削減(収率 +2pt ≈ feed −25 億/yr)。
- だが低Tは Dist2 cold-top で大半 infeasible。cold-top は recycle 組成との相互作用で BO/bounds では潰せない(実証)→ Dist2 段数/圧力の再設計、冷媒 tier、またはフローシート側で扱う問題。

**B3. 目的関数**: 生産量下限に張り付く TAC 最小化は、赤字下では loss 最小化と一致(= 現状で正しい)。変更不要。報告では「製造原単価 281 円/kg・利益」を併記して BO 結果の意味を明示するのが有効。

### Tier 3 — 総括
利益化は BO の洗練ではなく **(a) SM Dist3 の妥当性(B1)、(b) 収率↔Dist2 cold-top のトレードオフ(B2)、(c) 価格前提**で決まる。**次の一手は B1(SM Dist3 検証)を推す。**

---

## 4. 参照
- データ: `outputs/special_20260528_124401/`、`outputs/special_20260528_150021/`
- 関連 doc: `docs/bo_dist2_cond_shortfall_experiment_20260528.md`、`docs/solver_choice_rationale.md`
- BO 機構: `optimization/study.py`(`_default_constraints_func` 23 本、`_PhaseSwitchSampler`)、`optimization/objective.py`(`_store_diagnostics`)、`optimization/penalty_scale.py`(`default_schedule` 0.2→3.0)
