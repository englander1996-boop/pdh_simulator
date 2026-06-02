# hysys_cases/ — 段数別 HYSYS ケースファイル群

## 役割

HYSYS バックエンドが蒸留塔を解くときに読み込む **段数別の HYSYS ケースファイル (.hsc)** 置き場。
塔ごと・理論段数ごとに 1 ファイルを用意してあり、最適化器が段数を変えると対応する HSC が
切り替えて読み込まれる。

## サブ構成

| サブディレクトリ | 塔 | 収録段数 (ファイル名 = 段数.hsc) |
|---|---|---|
| `column1/` | Dist1 (塔1) | 30〜60 |
| `column2/` | Dist2 (塔2, 脱エタン塔) | 15〜80 |
| `column3/` | Dist3 (塔3, C3 スプリッタ) | 69〜200 |

- `{段数}.hsc` が HYSYS ケース本体。`.bk0` は HYSYS が作るバックアップファイル。

## 使い方・位置づけ

- `units/vle/hysys/registry.py` (`HsysRegistry`) が `hysys_cases/column{1,2,3}/` を走査し、
  要求された段数 → HSC パスを解決する。存在しない段数を要求すると `StageNotAvailableError`。
- そのため最適化器の N_stages 探索範囲は **各塔の HSC 収録範囲内に限定**する必要がある
  (col1=30-60 / col2=15-80 / col3=69-200)。
- 実際に HYSYS で塔を解くには HYSYS 有効環境 (本体 COM・ライセンス/VPN) が要る。
