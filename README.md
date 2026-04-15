# PDH Simulator

プロパン脱水素 (Propane DeHydrogenation) プロセスのシミュレータ。

## 成分記号

| 記号 | 物質 | 化学式 |
|------|------|--------|
| A | Propane (プロパン) | C3H8 |
| B | Propylene (プロピレン) | C3H6 |
| C | Hydrogen (水素) | H2 |
| D | Ethylene (エチレン) | C2H4 |
| E | Methane (メタン) | CH4 |
| F | Ethane (エタン) | C2H6 |

## 反応系

| 反応 | 式 | 種別 |
|------|----|------|
| r1 | A → B + C | 脱水素 (Dehydrogenation) |
| r2 | A → D + E | クラッキング (Cracking) |
| r3 | D + C → F | 水素化 (Hydrogenation) |

## ファイル構成

```
pdh_simulator/
├── README.md
├── data/
│   └── a_parameter_fitting.xlsx   # 触媒活性データ (0〜1000 min, 400〜700 °C)
├── monitor/
│   └── fitting_analysis.ipynb     # [済] 触媒活性補間解析・可視化
└── src/
    ├── config.py                  # [済] 定数・パラメータ定義
    ├── kinetics.py                # [済] 反応速度計算
    ├── thermo.py                  # [済] 熱力学計算
    └── catalyst_model.py          # [済] 触媒活性パラメータ計算
```

## 実装済みモジュール

### `src/config.py` — 定数・パラメータ定義

物理定数とすべての数値パラメータを集約するファイル。kJ → J の単位変換もここで完結させる。

| クラス / 変数 | 内容 |
|---|---|
| `R`, `T0` | ガス定数 [J/(K mol)]、基準温度 [K] |
| `Reaction1Params` | 脱水素反応パラメータ (k01, Ea1, deltaH, K0) |
| `Reaction2Params` | クラッキング反応パラメータ (k02, Ea2) |
| `Reaction3Params` | 水素化反応パラメータ (k03, Ea3) |
| `PDHConfig` | 上記3クラスをまとめるルートデータクラス |
| `ThermoParams` | 1成分の熱力学データ (dHf_298, Cp多項式係数 a〜d) |
| `THERMO_DATA` | 成分A〜F の `ThermoParams` 辞書 |

### `src/kinetics.py` — 反応速度計算

**クラス: `PDHKinetics(config: PDHConfig)`**

| メソッド | 説明 | 戻り値の単位 |
|---|---|---|
| `calculate(P_A, P_B, P_C, P_D, T, a, K_eq)` | r1, r2, r3 を一括計算して辞書で返す | mol m⁻³ s⁻¹ |

- 触媒劣化パラメータ `a` と平衡定数 `K_eq` は引数として外部から渡す設計
- 各速度定数 (`_k1`, `_K_B`, `_k2`, `_k3`) はメソッドに分離済み

### `src/thermo.py` — 熱力学計算

**クラス: `PDHThermo(data: Dict[str, ThermoParams])`**

| メソッド | 説明 | 戻り値の単位 |
|---|---|---|
| `calc_cp(component, T)` | 温度T における定圧比熱 Cp | J/(K mol) |
| `calc_enthalpy_change(component, T_start, T_end)` | Cp の解析的積分によるエンタルピー変化 | J/mol |

- `calc_enthalpy_change` は多項式の解析的積分で実装（scipy 不使用）

### `src/catalyst_model.py` — 触媒活性パラメータ計算

**関数: `calculate_activity_a(T_celsius, t_min) -> float`**

| 引数 | 単位 | 説明 |
|---|---|---|
| `T_celsius` | °C | 反応温度（有効範囲: 400〜700 °C） |
| `t_min` | min | 経過運転時間（有効範囲: 0〜30 min） |
| 戻り値 | - | 触媒活性 `a`（常に 0.0〜1.0 に収まる） |

- 手法：2次元スプライン補間（`RectBivariateSpline`, kx=3, ky=3）
- `data/a_parameter_fitting.xlsx` の 0〜30 min データをスクリプト内にハードコード（毎回ファイル読み込みなし）
- データ点上での再現誤差は機械精度レベル（≈ 1e-16）
- 入出力ともにクリッピング処理でシミュレータの発散を防止

```python
from src.catalyst_model import calculate_activity_a

a = calculate_activity_a(T_celsius=600.0, t_min=15.0)  # -> 0.24
```

### `monitor/fitting_analysis.ipynb` — 触媒活性補間解析

- `data/a_parameter_fitting.xlsx` から 0〜30 min データを読み込み
- `RectBivariateSpline` による補間曲面を構築し、全データ点での誤差ゼロを検証
- 既知温度（400〜700 °C）の補間カーブと実データの一致、および未知温度（525 °C）の補間結果を可視化

## 未実装 / 今後の作業

- [ ] 物質収支計算モジュール
- [ ] エネルギー収支計算モジュール
- [ ] 反応器モデル (PFR等)
- [ ] テストスイート
