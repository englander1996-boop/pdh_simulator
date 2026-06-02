# config/ — 固定運転条件の設定層

## 役割

「設計者が決め打ちしている運転条件」と「最適化器が触る設計変数」を明確に
分離するための設定層。原料仕様・目標生産量・主要昇圧/冷却ターゲット・
収束ソルバ設定・製品仕様 (feasibility 制約)・ペナルティ係数を
TOML ファイルに集約し、型付き `dataclass` にロードする。

設計変数 (反応器/PSA/膜/蒸留塔の dataclass) はコード側 (`optimization/`,
各 unit モジュール) に残り、optimizer が探索する。経済パラメータ
(電気代・蒸気代・触媒単価等) は本ディレクトリではなく
`src/cost_parameters.py` にある。

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `operating.toml` | 固定運転条件の本体 (値とその設計判断コメント)。 |
| `load.py` | TOML → 型付き `OperatingConfig` dataclass へのローダ + 各セクションの dataclass 定義。 |
| `__init__.py` | `OperatingConfig` 系 dataclass と `load_operating_config` を再エクスポート。 |

## 読み込み方法

```python
from config.load import load_operating_config
config = load_operating_config()          # config/operating.toml を読む
config = load_operating_config(path=...)  # 別ファイルを指定する場合
```

`load_operating_config()` は `operating.toml` を `tomllib` でパースし、
`OperatingConfig` (frozen dataclass) を返す。`optimization/` や
`flowsheet.evaluate()` はこの `config` を受け取って評価を行う。

## 設定セクション (`operating.toml` / `OperatingConfig`)

| セクション | dataclass | 主な内容 |
|---|---|---|
| `[product]` | `ProductSpec` | C3H6 目標生産量 `target_mta` [t/年]、分子量。 |
| `[feed]` | `FeedSpec` | LPG 原料組成 (C3H8:C4H10 = 9:1)、温度/圧力、外側ループ初期推定用 `yield_assumed`。 |
| `[pressure]` | `PressureSpec` | 主要昇圧ターゲット (ポンプ1出口・Comp2 出口・反応器入口) [Pa]。 |
| `[temperature]` | `TemperatureSpec` | 反応器後冷却・膜フィードの目標温度 [K]。 |
| `[solver.inner]` | `InnerSolverSpec` | リサイクル収束 (tear stream) 設定。相対許容誤差・floor・`method` (`successive_substitution`/`wegstein`)。 |
| `[solver.outer]` | `OuterSolverSpec` | 製品流量を目標に合わせる Fresh 調整ループ (片側相対 TOL)。 |
| `[solver.init]` | `InitSpec` | 反復の初期推定値 (Fresh=1500 基準の線形スケーリング)。 |
| `[spec]` | `SpecSpec` | 製品仕様の feasibility 制約 (C3H6 質量分率下限・H2 モル分率下限・生産量許容幅)。 |
| `[penalty]` | `PenaltySpec` | `effective_TAC` への加算ロジック (solver 失敗固定値・spec 違反の base/coef)。 |

`SolverSpec` は `inner`/`outer`/`init` の 3 サブ dataclass を束ねる。

## 位置づけ・依存

- 依存は標準ライブラリ (`tomllib`, `dataclasses`, `os`) のみ。
- `optimization/objective.py`, `optimization/pipeline.py`, `flowsheet.evaluate()`,
  `optimization/feasibility.py` (solver 失敗閾値の導出) などから参照される。
- `[spec]` と `[penalty]` の値は最適化の目的関数・制約に直結する
  (詳細は `optimization/SPEC_bo.md`)。
- 一部の値は設計判断の途中段階・暫定値を含む。`!仮置き` / 【確認中】 の扱いは
  リポジトリ直下の `KNOWN_PLACEHOLDERS.md` を参照。
