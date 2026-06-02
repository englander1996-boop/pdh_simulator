# stream/ — 共通プロセスストリーム型

## 役割

全ユニット (Swing/Radial 反応器, PSA, 膜, 蒸留塔, Mixer, Cooler, Compressor) の
入出力を表す汎用データ構造 `ProcessStream` を提供する。個別ユニットの
`FeedStream` / `OutStream` 間を橋渡しする「共通通貨」として使われる。

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `stream.py` | `ProcessStream` dataclass の定義 (成分モル流量・温度・圧力)。 |
| `__init__.py` | `ProcessStream` を再エクスポート。 |

## データ構造 (`ProcessStream`)

`@dataclass` で、属性は 3 つのみ。

| 属性 | 型 | 単位 | 意味 |
|---|---|---|---|
| `F_in` | `Dict[str, float]` | kmol/h | 成分別モル流量 (キーは下記の成分コード) |
| `T_in` | `float` | K | 温度 |
| `P_in` | `float` | Pa | 圧力 |

メソッド:

- `total_flow() -> float` … `F_in` の全成分の合計モル流量 [kmol/h] を返す。

### 成分キー (A–Z) の意味

`F_in` のキーは 1 文字の成分コード。`stream.py` の docstring と
`simulation/display.py` の `_COMP_NAMES` で定義され、両者は一致している。

| キー | 化学種 | 名称 |
|---|---|---|
| `A` | C3H8  | プロパン (原料 LPG の主成分・反応物) |
| `B` | C3H6  | プロピレン (目的製品) |
| `C` | H2    | 水素 (副生・製品候補) |
| `D` | C2H4  | エチレン |
| `E` | CH4   | メタン |
| `F` | C2H6  | エタン |
| `Z` | C4H10 | ブタン (原料 LPG 中の 1 割成分) |

## 位置づけ・依存

- 依存なし (標準ライブラリ `dataclasses` / `typing` のみ)。最下層の純データ型。
- 反応器・分離器・蒸留塔・経済計算など上位の全モジュールから参照される。
- 表示層 (`simulation/display.py`) はこの成分キー対応表を用いて化学種名を表示する。
