# monitor/ — 解析・検証ノートブック群

## 役割

シミュレータ各部の挙動を可視化・検証し、レポートの主張を実コードで裏取りするための
Jupyter ノートブック置き場。各ノートは設計条件を冒頭で編集 → 全セル実行で結果を確認する。
一部はレポート貼付用の図 (PNG/PDF) も出力する。

> **反応器は現在 Catofin型 浅床軸流スイング (catofin) が採用**(径方向流 radial は一時検討後に撤去)。
> 反応器設計を扱うノート/図は catofin 基準。BO best は `outputs/main_20260605_170938/best.json` (trial #201)。

## ファイル一覧

| ノートブック | 何を可視化/検証するか |
|---|---|
| `reactor_conversion_ceiling.ipynb` | Catofin の単通転化率 (~41%) と選択率 (~80%) が熱力学+速度論+**HGM 等価熱補償**で必然であることを示す。HGM が床温を T_in−ΔT_max(=612℃)に維持するため、**転化率は断熱平衡頭打ちでなく床体積 (L_bed/N_online) でスケール**(radial 断熱床の「体積無効」と逆)。実走で X=40.6%, S=80.1%, 炭素収支 closure=100% を再現。 |
| `reactor_pressure_drop_and_geometry.ipynb` | 反応器の圧損 (Ergun) と幾何を検証。現実的 2–6mm 触媒では**軸流深床が 0.5 bar でどの床長も ΔP 過大=全 infeasible**を実証し、**Catofin 中浅床(L_bed≤3m)+多基並列(N_online)**で ΔP/P≈3.3% で成立することを示す (#201)。捨てたのは「軸流」でなく「深床」。径方向流は内部品/分配/再生温度の不確かさで撤去。 |
| `stage_comparison.ipynb` | 反応器**形式の選択**を同一入口で実走比較 (軸流深床=圧損破綻 / radial 多段=段数↑で選択率↓ / **Catofin 浅床多基**=高選択率を維持)。総収率≈選択率なので Catofin の低per-pass転化・高選択率 (S=80.1%) が radial 多段 (3段 S=68.8%) より総収率で有利、かつ段間再加熱不要で単純、を示す。Catofin の設計レバー N_online 掃引も。 (HYSYS 不要) |
| `heat_integration_targeting.ipynb` | `flowsheet/heat_integration.py` を教科書例題で動作確認。Q_H,min/Q_C,min (Problem Table)、複合線・GCC、A_total (Bath 式)、N_HE,min、ΔTmin 感度、潜熱対応を検証。 |
| `dist2_inspect.ipynb` | Dist2 (脱エタン塔, partial cond) の rigorous 入出力を詳細検証。ストリーム・マスバランス・T/K/x/y プロファイル・装置/CAPEX、HYSYS との比較、D 感度スキャン。 |
| `fitting_analysis.ipynb` | 触媒活性パラメータ `a(T, t)` の 2 次元スプライン補間 (RectBivariateSpline) を可視化。データ源は `data/a_parameter_fitting.xlsx` (0-30 min)。 |
| `lpg_price_breakeven.ipynb` | 最適設計 (catofin best #201) を固定したまま LPG 原料単価だけを振り、年間利益と損益分岐価格を算出。原料費だけが線形にスケールすることを実コードで実証。図を PNG/PDF 出力。**結果: 損益分岐 ≒64.4 円/kg、価格感度 7.32 億円/年·(円/kg)⁻¹、原料費は TAC の 66%**。 |
| `dist3_sm_domain.ipynb` | **r3(#201 run で 2 番目の失敗 ~24%、最多は反応器 r_rx ~40%)の根本** = Dist3 SM classifier の受理面(In_Flow×In_Propane)を可視化。受理は「高流量(≥0.45, propane不問)」か「低流量×超高純度(In_Propane≤0.05)」の2択で中間が谷。**#201 best は高流量(In_Flow≈0.562→学習域上限0.5を超え clamp)に着地し頑健に受理**(旧#227 の低流量コーナーから安全側へ移動)。ただし In_Flow が学習上限超=外挿の残課題。off-best trial が谷に落ちて r3 失敗。(HYSYS 不要) |
| `reactor_axial_profile.ipynb` | **反応器の床内部 z 方向プロファイルとサイクル t 方向経時**を実コード (`swing._ode_axial` を per-vessel 直接積分) から再構成 (報告書 §4.5 / スライド7 用)。**Fig1★** T(z)・X(z) の HGM補償 vs 無補償断熱 対比 (HGM は床温612℃維持で X 伸長、断熱は冷却で頭打ち) / **Fig2★** 微分・累積選択率 S_diff(z)・S_int(z) (床後半で崩落→出口へ収束) / **Fig5★** サイクル経時 a(t)・X(t)・S(t) (失活で性能低下、**時間平均 X̄=40.6% が公式と完全一致**) / Fig3 分圧+P(z) Ergun圧損 / Fig4 速度 r1,r2,r3(z) / Fig6 ΔTmax感度・Q_HGM軸方向。priority 図(1/2/5)は PNG/PDF 出力。既存反応器ノート・report図と非重複。(HYSYS 不要、`tools/build_reactor_axial_profile_nb.py` で生成) |
| `hgm_dtmax_sensitivity.ipynb` | **catofin の最大の !仮置き** = HGM 等価熱補償 ΔT_max(=50K)を 30/50/80K で実走し、単通転化率・生産量・TAC・feasibility への効きを定量化。**結果(#201): 50K が feasible 最低 TAC(eff_TAC=1055)、80K(弱 HGM)で X 40.6%→33%→TAC↑(1075)、30K は反応器制約(r_rx)で不可**。(HYSYS 必須・pdh-venv カーネル) |
| `membrane_degradation_sensitivity.ipynb` | 膜 Q_A/α(!仮置き, 文献値)を劣化方向に振り、C3H6 回収率・TAC・feasibility の頑健性を実走。**結果(#201): Q_A 30%劣化で 回収率 60→40%・TAC +3.5%(1092)、α 30%劣化はほぼ無影響(+0.2%)、いずれも feasible = 膜性能は TAC の主リスクでない(頑健)**。(HYSYS 必須・pdh-venv) |
| `bound_relax_gain.ipynb` | BO 最適が探索箱の端に張り付く変数を端の外へ 1 点動かし TAC の局所下げ代を見積り。**結果(#201): col2_N/col3_N/d_p を端の外へ動かすと ΔTAC +0.7〜+3.7(いずれも悪化)、F_fresh を下限(1450)外へ下げると生産量厳格化で infeasible = 端の張り付きは人為的箱制約でなく実トレードオフ(BO はほぼ真の最適に到達)**。単点摂動の局所見積り。(HYSYS 必須・pdh-venv) |

## 同梱の図ファイル

- `lpg_breakeven_profit.png` / `.pdf` — LPG 価格 vs 年間利益 (本命図)。
- `lpg_breakeven_unitcost.png` / `.pdf` — LPG 価格 vs C3H6 製造原単価 (補足図)。

いずれも `lpg_price_breakeven.ipynb` が `rule.md` の図表ルール (グリッドなし・四辺内向き目盛り・グレースケール+ハッチング) に沿って生成する。

## 使い方・位置づけ

- 実行は repo 直下を import パスに通す前提 (各ノート冒頭で `..` を sys.path に追加)。
- HYSYS を呼ぶノート/セル (`lpg_price_breakeven` のセル3, `hgm_dtmax_sensitivity` / `membrane_degradation_sensitivity` / `bound_relax_gain`) は HYSYS 有効環境・VPN が必要。
  - これら HYSYS ノートは **`.venv` カーネルで実行**すること(optuna/pywin32 がシステム Python に無いため)。`.venv` カーネルは下記で登録済み(表示名「PDH .venv」):
    ```
    .venv\Scripts\python.exe -m ipykernel install --user --name pdh-venv --display-name "PDH .venv"
    ```
    nbconvert で実行する場合: `... -m jupyter nbconvert --execute --inplace --ExecutePreprocessor.kernel_name=pdh-venv <nb>`。HYSYS 不要のノートはシステム Python でも可。
- `reactor_conversion_ceiling` / `reactor_pressure_drop_and_geometry` / `stage_comparison` は `tools/build_*_nb.py` で生成する(直接編集せず builder 側を編集)。**3 つとも catofin 基準で生成・実行済み**。再生成は例えば:
  ```
  .venv\Scripts\python.exe tools\build_reactor_ceiling_nb.py
  .venv\Scripts\python.exe -m jupyter nbconvert --to notebook --execute --inplace monitor\reactor_conversion_ceiling.ipynb
  ```
  (`build_stage_comparison_nb.py` / `_build_reactor_pdrop_nb.py` も同様。いずれも HYSYS 不要。)
