# HYSYS Dist2 cold-top 失敗の連続シグナル化 — 実験記録(結論:効果なし)

**作成日**: 2026-05-28
**対象**: `special.py`(HYSYS+SM 全21変数 BO) / `flowsheet/run_one_pass.py`・`optimization/objective.py`・`optimization/study.py`
**目的**: r2(Dist2 凝縮器ΔT不成立 = 塔頂が冷えすぎて凝縮不可)失敗が全 trial の ~42% を占める問題に対し、**Dist2 SM 化を断念した上での代替手**として「凝縮器ΔTの連続制約シグナル化(Rec 2)」を試し、feasible 率・best TAC への効果を検証した記録。**結論は効果なし**。コードは保持(将来の pre-screen 土台として再利用可)。

---

## 1. 背景・問題提起

`special.py` の BO(300 trial、TPE+QMC startup50、seed42)では Dist2 のみ HYSYS で評価する。
ベースライン run `outputs/special_20260528_124401`(用役単価0527改定後)を分解したところ:

- 完了 300 / feasible 107(**35.7%**)、best #171 = **1129.51 億円/年**。
- infeasible 193 件のうち **r2(Dist2)= 145 件が最大**(= 全 trial の 48%)。
- r2 の内訳:**125件(86%)が「凝縮器ΔT不成立」**(塔頂温度 < 最冷媒エチレン -100°C + ΔT、`units/vle/hysys/provider.py:594-597`)、20件(14%)が HYSYS COM エラー。
- 重要:この cold-top 失敗は HYSYS が列を収束させ塔頂温度を出した**後**の物理判定であり、N_needed / dT_max がともに 0。よって `_compute_dist_shortfalls`(`run_one_pass.py`)は `dist2_{N,dT}_shortfall=0` を返し、`_default_constraints_func`(`study.py`)では **`unknown_failure=1.0`(方向不明 infeasible)**としてしか TPE に伝わっていなかった。

→ 仮説:**「silent な cold-top 失敗が TPE の feasible 学習のボトルネック」**。塔頂を暖める方向(col2_p↑ / reflux↓)の連続勾配を与えれば、TPE が cold-top 領域を避け、feasible 率が上がるはず。

## 2. 仮説と狙い

cold-top 失敗に「塔頂が凝縮可能下限を何K下回るか」の連続シグナルを与え、TPE に col2_p↑/reflux↓ の勾配を学習させる。**feasible 率の上昇**(と best TAC の改善)を期待値とした。is_feasible 判定は不変なので best 悪化リスクは無い設計。

## 3. 実装(3ファイル、最小侵襲)

新シグナル `dist2_cond_shortfall = max(0, (T_util + 10) − T_top)` [K](T_util+10 = utility 戻り温度 = 凝縮可能下限)。

1. **`flowsheet/run_one_pass.py` `_compute_dist_shortfalls`**: 失敗 message(`"...condenser ΔT 不成立 (T_top=-99.4°C, T_util=-100.0°C)"`)を正規表現でパースし `dist{idx}_cond_shortfall` を出力。HYSYS 経路は T_top を equipment に構造化保存しないため message パースで復元。`import re` 追加。
2. **`optimization/objective.py` `_store_diagnostics`**: 診断ループの kind に `'cond'` を追加 → `dist2_cond_shortfall` を `trial.user_attrs` へ。
3. **`optimization/study.py` `_default_constraints_func`**: 制約ベクトルに **[22] = dist2_cond_shortfall × 0.1**(10K=1.0 正規化)を追加。`raw_total` にも算入し、cold-top が `unknown_failure` 扱いにならないよう修正。ベクトル長 22 → 23。

### スモークテスト(配線確認、HYSYS 1回)
cold-top 既知失敗 trial #0(col2_p=536, T_top=-99.4°C)を再評価:
- `one_pass['dist2_cond_shortfall'] = 9.4`(=(-100+10)-(-99.4))✓
- 制約ベクトル 23 要素(feasible/infeasible とも一致)✓
- cold-top trial: [22]=0.940、unknown_failure[21]=0.0(silent 解消)✓
- feasible trial: [22]=0.0(回帰なし)✓

## 4. 実験設定

- run: `outputs/special_20260528_150021`(ログ `outputs/special_run_condfix_20260528_150019.out.log`)
- N_TRIALS=300, N_STARTUP(QMC)=50, seed=42(ベースラインと同一)。
- 比較対象ベースライン: `outputs/special_20260528_124401`(同設定、cond シグナル無し)。
- 注意:制約ベクトルが 22→23 に変わるため trial50 以降の TPE 軌道はベースラインと分岐する(同一 seed でも別ラン扱い)。

## 5. 結果

| 指標 | ベースライン (124401) | cond-signal (150021) | 差 |
|---|---|---|---|
| 完了 / feasible | 300 / 107 (**35.7%**) | 300 / 104 (**34.7%**) | −0.9pt |
| best effective_TAC | **1129.51**(#171) | **1131.14**(#187) | +1.6 億(悪化) |
| r2(Dist2)失敗 | 145 | 141 | −4 |
| best 設計性格 | F=1560 / prod=1194 / 99.5wt% | F=1532 / prod=1194 / 99.5wt% | ほぼ同一 |

## 6. 判定

**効果なし。** feasible 率はわずかに低下、best はわずかに悪化、r2 は 4 件減のみ。いずれも **run 間ばらつきの範囲内**で、cond シグナルが有意に効いたとは言えない。仮説「silent cold-top が feasible 学習のボトルネック」は**支持されなかった**。

## 7. 効かなかった理由(考察)

1. **cold-top は相互作用効果**。塔頂温度は (col2_p, reflux) だけでなく上流リサイクル組成に依存し、リサイクル収束で揺れる。TPE が直接動かせる col2 ノブと塔頂温度の関係が一意でないため、「暖める方向」の勾配が素直に立たない。
2. **Optuna 制約 TPE は本質的に feasible/infeasible の二分**で trial を分類する。cold-top は前も後も「infeasible」であり、連続値(勾配)は「どの infeasible がより近いか」の重み付けにしか効かない。23 制約の中の 1 本では希釈され、サンプリング挙動を動かすレバーにならなかった。

## 8. 副次的収穫(BO の頑健性=局所解でない傍証)

制約構成が異なる 2 ラン(22 要素 vs 23 要素、trial50 以降で軌道が分岐)が、ともに **~1130 億円・同じ設計性格(prod≈1194 / F≈1530-1560 / 99.5wt%)に収束**した。これは小さな**再現実験**であり、「~1130 は seed/制約構成依存の局所罠ではなく頑健な basin」の傍証になる(報告書での「BO 結果が局所解でない」主張の一本として利用可)。

## 9. 含意・今後

- **r2/cold-top の浪費(42%)は制約シグナルでは潰せない**(探索空間 + 相互作用に内在)。
- Rec 1(HYSYS 前の露点プレスクリーン)は**計算時間は節約できるが feasible 率・best は改善しない**(提案される trial は同じ、棄却が速くなるだけ)。本実験で「TPE に避けさせる」のが困難と判明した以上、Rec 1 の価値は純粋に高速化のみ。
- BO はこの構成では**既に達成可能な best(~1130)付近**にいる可能性が高い。さらなる改善は BO の工夫ではなく**設計レバー(収率/反応器-原料トレードオフ)や探索範囲**側で行うべき。

## 10. コードの扱い

ユーザー判断(2026-05-28)により**変更は保持**。無害(cold-top 以外は常に 0、is_feasible 判定不変)で、将来 Rec 1(露点プレスクリーン)を実装する際の `dist2_cond_shortfall` 復元ロジックを再利用できるため。
