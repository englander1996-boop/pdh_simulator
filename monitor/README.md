# monitor/ — 解析・検証ノートブック群

## 役割

シミュレータ各部の挙動を可視化・検証し、レポートの主張を実コードで裏取りするための
Jupyter ノートブック置き場。各ノートは設計条件を冒頭で編集 → 全セル実行で結果を確認する。
一部はレポート貼付用の図 (PNG/PDF) も出力する。

## ファイル一覧

| ノートブック | 何を可視化/検証するか |
|---|---|
| `reactor_conversion_ceiling.ipynb` | 反応器の単通転化率 (~28%) と選択率 (~80%) がバグではなく熱力学 (断熱平衡) + 反応速度論で必然であることを厳密に示す。実走で X=28.3%, S=79.4% を再現。 |
| `reactor_pressure_drop_and_geometry.ipynb` | 反応器の圧損 (Ergun) と幾何 (軸流深床 vs 径方向流) を検証。軸流深床が 0.5 bar で全 infeasible となる問題と径方向流での解決を示す。 |
| `stage_comparison.ipynb` | 反応器設計の変遷 (単段断熱床→径方向流→多段化) と段数 1/2/3/4 段の比較 (累積転化率・選択率・コスト)。なぜ 3 段かを示す (HYSYS 不要)。 |
| `heat_integration_targeting.ipynb` | `flowsheet/heat_integration.py` を教科書例題で動作確認。Q_H,min/Q_C,min (Problem Table)、複合線・GCC、A_total (Bath 式)、N_HE,min、ΔTmin 感度、潜熱対応を検証。 |
| `dist2_inspect.ipynb` | Dist2 (脱エタン塔, partial cond) の rigorous 入出力を詳細検証。ストリーム・マスバランス・T/K/x/y プロファイル・装置/CAPEX、HYSYS との比較、D 感度スキャン。 |
| `fitting_analysis.ipynb` | 触媒活性パラメータ `a(T, t)` の 2 次元スプライン補間 (RectBivariateSpline) を可視化。データ源は `data/a_parameter_fitting.xlsx` (0-30 min)。 |
| `lpg_price_breakeven.ipynb` | 最適設計を固定したまま LPG 原料単価だけを振り、年間利益と損益分岐価格を算出。原料費だけが線形にスケールすることを実コードで実証。図を PNG/PDF 出力。 |

## 同梱の図ファイル

- `lpg_breakeven_profit.png` / `.pdf` — LPG 価格 vs 年間利益 (本命図)。
- `lpg_breakeven_unitcost.png` / `.pdf` — LPG 価格 vs C3H6 製造原単価 (補足図)。

いずれも `lpg_price_breakeven.ipynb` が `rule.md` の図表ルール (グリッドなし・四辺内向き目盛り・グレースケール+ハッチング) に沿って生成する。

## 使い方・位置づけ

- 実行は repo 直下を import パスに通す前提 (各ノート冒頭で `..` を sys.path に追加)。
- HYSYS を呼ぶノート/セル (例: `lpg_price_breakeven` のセル3) は HYSYS 有効環境・VPN が必要。それ以外のセルは不要。
- `reactor_conversion_ceiling` / `stage_comparison` は `tools/build_*_nb.py` で生成される (直接編集せず builder 側を編集)。
