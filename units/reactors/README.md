# プロパン脱水素反応：仮想複数管切り替え反応器シミュレーションモジュール仕様書

## 1. モジュールの目的

本モジュールは、触媒失活を伴うプロパン脱水素反応を行う断熱反応器システム（スイング操作）をシミュレーションし、以下の情報を出力することを目的とする：

- 後段の分離工程に渡す時間平均ストリーム情報
- 装置コスト・基数

## 2. データ構造の定義（Input）

### ① DesignVars（設計変数）

最適化アルゴリズムが操作する変数

| 変数 | 単位 | 説明 |
|------|------|------|
| $T_{in}$ | K | 反応器入口温度 |
| $z_{cat}$ | m | 触媒層の長さ（空間積分の範囲） |
| $t_{cyc}$ | min | 1サイクルあたりの反応フェーズ運転時間（時間積分の範囲） |
| $D$ | m | 反応器1基の直径（断面積$A$およびコスト計算に使用） |

### ② FeedStream（入口流体条件）

前段工程またはリサイクルから供給されるストリーム

| 変数 | 単位 | 説明 |
|------|------|------|
| $F_{in}$ | kmol/h | 各成分の入口モル流量（C3H8, C3H6, H2, C2H4, CH4, C2H6） |
| $T_{feed}$ | K | 加熱炉入口（予熱前）の原料温度 |
| $P_{in}$ | Pa | 反応器入口圧力 |

### ③ FixedParams（固定定数・制約条件）

システム上の制約および物理定数。熱力学・反応速度・触媒失活に関する物性値やパラメータは Section 6 の外部モジュール側で管理する。

| 変数 | 値 | 単位 | 説明 |
|------|------|------|------|
| $t_{regen}$ | 30.0 | min | 触媒再生時間（固定） |
| $V_{cat\_max\_per\_vessel}$ | 200.0 | m³ | 1基あたりの最大触媒量（固定） |
| $\varepsilon$ | 0.5 | — | 触媒充填率（容器体積が触媒量の2倍という制約より） |
| $\rho_p$ | — | kg/m³ | 触媒充填密度 |

## 3. データ構造の定義（Output）

### ① EffluentStream（出口流体情報）

後段の分離工程や熱統合計算に渡す

| 変数 | 単位 | 説明 |
|------|------|------|
| $F_{out\_avg}$ | kmol/h | 各成分の出口モル流量の時間平均値 |
| $T_{out\_avg}$ | K | 出口温度の時間平均値 |
| $Q_{preheat}$ | GJ/h | $T_{feed}$から$T_{in}$まで加熱するのに必要な熱量 |
| $P_{out}$ | Pa | 出口圧力（今回は$P_{in}$と同値） |

### ② EquipmentCost（装置・経済性情報）

プラント全体の投資額評価に使用する

| 変数 | 説明 |
|------|------|
| $V_{vessel\_actual}$ | 1基あたりのプロセス容器容積 [m³] |
| $N_{parallel}$ | 200m³制約により必要な並列基数 |
| $N_{swing\_sets}$ | 再生時間をカバーするための切り替えセット数：$\lceil t_{regen}/t_{cyc} \rceil + 1$ |
| $N_{reactors\_total}$ | プラント全体の総反応器基数：$N_{parallel} \times N_{swing\_sets}$ |
| $Catalyst\_Weight\_Total$ | kg | システム全体の必要触媒総量 |
| $Reactor\_CAPEX$ | 億円 | 全基分の建設コスト合計 |

### ③ PerformanceMetrics（プロセス指標）

| 変数 | 単位 | 説明 |
|------|------|------|
| Conversion | % | プロパンの単通反応率（時間平均値ベース） |
| Selectivity | % | プロピレンの選択率（時間平均値ベース） |

## 4. 変数と数式の相関（ロジック定義）

### A. 物理モデル（微分方程式レイヤー）

**入口温度（$T_{in}$）**
- 反応速度定数$k$の初期値を決定
- 熱収支式の積分初期条件となる

**フィード流量（$F_{in}$）**
- 物質収支式の積分初期条件となる
- 分圧$P_i$を通じて反応速度を支配

**直径（$D$）**
- 断面積$A$を決定
- 管内の空塔速度および単位長さあたりの触媒量を決定

**触媒層長さ（$z_{cat}$）**
- 空間方向の積分範囲を規定

**運転時間（$t_{cyc}$）**
- 時間方向の積分範囲を規定
- 触媒劣化度$a$を通じて反応速度を低下させる

### B. 装置・コストモデル（代数計算レイヤー）

**出口平均値**

反応終了時間$t_{cyc}$までの全瞬時値を時間平均して$F_{out\_avg}$、$T_{out\_avg}$を算出する

**基数計算**

$$V_{cat\_total\_phase} = A \times z_{cat} \times \varepsilon$$

$$N_{parallel} = \left\lceil \frac{V_{cat\_total\_phase}}{200} \right\rceil$$

$$N_{swing\_sets} = \left\lceil \frac{30}{t_{cyc}} \right\rceil + 1$$

**CAPEX**

1基あたりの容器体積：

$$V_{vessel} = \frac{V_{cat\_total\_phase}}{N_{parallel}} \times 2$$

全体容器容積に対して基数$N_{reactors\_total}$倍して算出する

## 5. I/O 方法の実装仕様

- **モジュール名**: `swing_reactor_simulator.py`
- **データクラス定義**: `dataclasses`モジュールを使用して上記データクラスを定義
- **メイン関数**: `simulate_swing_reactor_system(design, feed, fixed)`
  - 入力：DesignVars、FeedStream、FixedParams
  - 出力：SimulationResultオブジェクト（EffluentStream、EquipmentCost、PerformanceMetricsを内包）
- **数値積分**: `scipy.integrate.solve_ivp`を使用し、4次のルンゲ・クッタ法で空間積分を実施

## 6. 外部依存モジュール（Global Functions）

本モジュールは、`src/` モジュールをインポートして以下の関数をモジュールレベルで定義する。パラメータを引数で引き回す必要はない。

### ① `calc_a(t, T, P)`

触媒劣化度を返す関数

| 引数 | 単位 | 説明 |
|------|------|------|
| $t$ | min | 反応フェーズの経過時間 |
| $T$ | K | 反応器温度 |
| $P$ | Pa | 反応器圧力 |

| 戻り値 | 説明 |
|--------|------|
| $a$ | — | 触媒活性度（0〜1の範囲、1.0が新触媒、0に近いほど失活） |

### ② `calc_Cp(T)`

各成分のモル熱容量を返す関数
|--------|------|
| `dict` または `list` | — | 各成分（C3H8, C3H6, H2, C2H4, CH4, C2H6）のモル熱容量 [J/(mol·K)] |

### ③ `calc_rate_constants(T)`

反応速度定数および平衡定数を返す関数

| 引数 | 単位 | 説明 |
|------|------|------|
| $T$ | K | 温度 |

| 戻り値 | 説明 |
|--------|------|
| `dict` | — | 反応速度定数や平衡定数を含む辞書。キー例：`k_forward`, `k_backward`, `Keq` など |
```python
from src.catalyst_model import calculate_activity_a
from src.thermo import PDHThermo
from src.kinetics import PDHKinetics

# モジュールレベルでインスタンス化し、ラッパー関数として定義
# 微分方程式の時間・空間積分内で使用
a = calc_a(t, T_local, P_local)
Cp_dict = calc_Cp(T_local)
rate_consts = calc_rate_constants(T_local)
```
