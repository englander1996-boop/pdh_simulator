# exp/ — 単一設計点の詳細評価・感度スクリプト群

## 役割

最適化器 (`main.py` / `sub/`) を回さず、**1 つの設計点を固定して全フローシートを詳細評価**したり、
ベスト設計を起点に **特定パラメータだけを振って感度を見る** ための独立スクリプト置き場。
いずれも `.\.venv\Scripts\python.exe exp\<file>.py` で単体実行する。

## ファイル一覧

| ファイル | 用途 |
|---|---|
| `exp3.py` | **主力。HYSYS + SM バックエンドでリサイクルあり PDH 全系を詳細評価** (Dist1/Dist3=SM、Dist2=HYSYS COM、反応器・PSA・膜・合流は本実装)。HYSYS 探索変数 (`hysys_spec_value` / `hysys_feed_stage`) を直接振る。N_stages は HSC 存在範囲のみ (col1=30-60 / col2=15-80 / col3=69-200)。 |
| `exp_membrane_sensitivity.py` | 膜性能 (透過度 Q_A・選択性 α・膜寿命 τ) を OFAT で劣化させ、ベスト設計の TAC と feasibility がどこで崩れるかを定量化 (レビュー指摘 #1)。 |
| `exp_membrane_recycle_coupling.py` | 膜面積 A_mem だけを振り「A_mem → stage cut → 反応器入口 C3H6 → 転化率/選択率 → TAC」の系統連成を 1 枚に並べる (レビュー指摘 #4)。 |
| `exp_psa_sensitivity.py` | PSA 吸着材データ (KFa・q_s・a・ρ_b) を係数で振り、塔数・H2 回収率・TAC の頑健性を評価。揺れが大きければ PSA を暫定扱いと判断する (PSA レビュー対応)。 |

## 使い方・位置づけ

- バックエンドは sensitivity 3 本とも `main.py` と同一 (Dist1/3=SM、Dist2=HYSYS、反応器=径方向流) なので、
  headline の TAC と直接比較できる。**HYSYS を使うため `.venv` (HYSYS 有効環境) で実行すること。**
- sensitivity スクリプトは引数にベスト設計 (`best.json`) のパスを取る (省略時は `outputs/main_*/best.json` の最新を自動採用)。
- `exp3.py` は HYSYS 探索変数が `main.py` / 最適化器と異なるため最適化器には組み込まず、独立運用する。
