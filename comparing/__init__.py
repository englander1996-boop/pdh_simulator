r"""
comparing/ — 過去レポートの「問題のある最適化のやり方」を本シミュレータ上で忠実に再現し、
その問題が結果に与える損失を **定量化** するパッケージ。

目的 (2026-05-28 ユーザ合意):
  「最適化方法論の問題点_詳細レポート.md」が指摘する P01〜P12 の方法論的欠陥を、
  素朴な最適化器を作ることではなく **欠陥そのものを忠実に再現** して走らせ、
  BO (main.py 旧special の既存 best) を基準に「この欠陥は年間 X 億円損している (+Y%)」と
  定量的に示す。これにより BO 採用の有用性を逆説的に立証する。

構成 (2026-05-29 ユーザ指示で確定):
  case ごとにフォルダを分け、各フォルダの main.py が **独立したスクリプト** として単独で走る。
  中央ディスパッチャ (まとめて比較する main.py) は作らない (= エラーの元なので廃止)。
  共通の土台だけ shared/ に置いて流用する。

  shared/space.py     : 21 変数の探索空間 + build_design (main.py 旧special から移植)
  shared/simulator.py : flowsheet.evaluate の薄いラッパ + CONFIG (purity 緩和)
  shared/reporting.py : main.py 旧special の I/O 機構 (ライブログ/trials.csv/top-N/README) の忠実移植
  shared/harness.py   : enqueue 駆動 Optuna study エンジン (実行・記録を main.py 旧special と揃える)

  caseN = PN (1 case = 1 problem)。実装済み (この sim で再現可能なもの):
  case_p01_subsystem/main.py  : P01 部分最適化 (サブシステム別)           [subsystem]
  case_p02_pinch/main.py  : P02 熱統合の後置 (post-hoc pinch)         [pinch_posthoc]
  case_p04_sequential/main.py  : P04 1次元逐次 (+P05 も観察)              [sequential_1d]
  case_p05_grid/main.py  : P05 整数の連続扱い (粗いグリッド)         [grid_search]
  case_p06_multistart/main.py  : P06 大域性未保証 (multi-start)           [multistart]
  case_p12_converge/main.py : P12 検証なし (収束未確認・複数巡)         [converge_check]

  C×P の意味ある組み合わせ (documented practice をまるごと再現):
  case_combo_typical/main.py : C01 1次元逐次 → C04 後置ピンチ の典型レポート全体
                         (P02+P04+P05+P06+P12 を 1 本で同時体現。コーパスの 96.7%)。

  実在レポートの手法をそのまま再現 (テーマ非依存。著者名は repo に残さず匿名 ID。詳細は
  REPORT_METHODS_ANALYSIS.md):
  case_rep_pdh2025a  : PDH 2025 — C3 分離塔を 圧力→段数→feed の逐次1次元
  case_rep_pdh2024a  : PDH 2024 — 反応器の単通転化率を 1 次元掃引
  case_rep_pdh2025b  : PDH 2025 — 1 塔の段数のみ 1 次元掃引
  case_rep_propro2019: n-プロピルプロピオネート 2019 — 反応条件の 3×3 粗グリッド (C02)
  case_rep_toluene2023: トルエン脱アルキル 2023 — 反応器の多変数フルグリッド絨毯爆撃 (C03)

  P10 = case_p10_fug で「FUG vs HYSYS真値 の精度乖離」として再現済み。
  未着手/要検討: P03 (要 flowsheet 改造: パージ変数), P07/P08/P09 (要 CO2/NPV 等の仮定),
  P11 (定常のみで不可)。全 case の索引は CASES.md。

実行: 各 case を個別に  `.\.venv\Scripts\python.exe comparing\caseN\main.py`
  結果は comparing/results/<method>_<ts>/ に出る。BO との比較はユーザ側で best.json を突合。

バックエンド: 本丸 main.py (旧 special.py) と同じ Dist1=SM / Dist2=HYSYS / Dist3=SM。**FUG は使わない**。
  HYSYS は単一 COM・並列不可・~150s/eval なので、各手法はレポートが実際に行った粗いサンプリングを
  忠実に再現し評価回数を抑える。なお HYSYS はマシン依存で、無い PC では実行できない (import は可)。
"""
