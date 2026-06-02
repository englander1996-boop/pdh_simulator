# units/vle/hysys/ — HYSYS COM 連携層 (VLE バックエンド)

蒸留塔 (column1/2/3) の気液平衡計算を **Aspen HYSYS** に委譲する COM 連携バックエンド。
`units/vle/` は VLE 実装の差し替えポイントで、その HYSYS 実装が本ディレクトリ。

## 役割

PDH 側の設計変数 (`ColumnTunables`) + 入口ストリーム (`ProcessStream`) を受け取り、段数別
HSC ファイルを開いて HYSYS ソルバで厳密に蒸留塔を解き、出口値を PDH 側の `DistResult`
に変換して返す。CAPEX は HYSYS 出力 (流量・温度・組成・熱量) から `cost_calculator` で
後計算する。非収束/エラーは `feasible=False` の penalty DistResult で吸収し、BO 探索から
自動除外する。

蒸留塔の `solver_method == 'hysys'` のときに `column*.py` からディスパッチされる
(special フォーク・検証用)。構造パラメータ (段数) は COM で変更できないため、**段数別に
事前作成した HSC ファイルを切替える**方式を採る。

## アーキテクチャ (4 層)

```
column*.py (solver_method='hysys')
   │  solve_columnX_via_hysys(feed, tunables)
   ▼
provider.py   HysysVLEProvider … 高水準 API (3塔の solve)
   │            ├ ColumnTunables+ProcessStream → Column{1,2,3}Input 組立
   │            ├ ColumnResult → DistResult 変換 (utility/HE/CAPEX 後計算)
   │            ├ _SessionCache … app 1つで case を swap して高速化
   │            └ Dist1 メモ化 (入力不変時の cold 結果再利用、warm-start ではない)
   ▼
registry.py   HsysRegistry … 段数 → HSC ファイルパス解決 (hysys_cases/columnN/{N}.hsc)
   ▼
session.py    HysysSession … HYSYS 起動/HSC open-close の context manager
   │            ├ HysysPopupMonitor … モーダル警告を別スレッドで自動押下
   │            ├ set_feed_stage … SpecifyFeedLocation 方式でフィード段書込み
   │            ├ run_column_and_wait … sub-flowsheet Solver で収束待機
   │            └ swap_case … app 保持のまま case だけ Close→Open
   ▼
adapters/     塔ごとの HYSYS オブジェクト名・スペック名・流量単位を吸収
              Column{1,2,3}Adapter: Column{1,2,3}Input → 書込み → solve → ColumnResult
```

## ファイル一覧

| ファイル | 役割 |
|---|---|
| `__init__.py` | 公開 API の re-export (HysysSession, HsysRegistry, set_feed_stage 等) |
| `session.py` | HYSYS COM 低レイヤ。`HysysSession` (起動/open/close)、`HysysPopupMonitor` (モーダル自動 dismiss)、`set_feed_stage`、`run_column_and_wait`、`swap_case`。例外 `HysysConnectionError`/`HysysConvergenceError` |
| `provider.py` | 高水準 API `HysysVLEProvider`。入力組立・結果変換・セッションキャッシュ・Dist1 メモ化。process-global singleton (`get_default_provider`) |
| `registry.py` | `HsysRegistry`。`hysys_cases/column{1,2,3}/{N}.hsc` を走査し段数→パス解決。`StageNotAvailableError` |
| `adapters/__init__.py` | アダプタ層の re-export |
| `adapters/types.py` | 入出力契約 dataclass: `Column{1,2,3}Input`, `ColumnResult` |
| `adapters/base.py` | 3塔共通ヘルパ (SPR-1 圧力/流量書込み、feed 組成/T/P 書込み、出力回収、empty sentinel 判定) |
| `adapters/column1.py` | Dist1 アダプタ (主スペック = Comp Fraction-2) |
| `adapters/column2.py` | Dist2 アダプタ (主スペック = Reflux Ratio) |
| `adapters/column3.py` | Dist3 アダプタ (主スペック = Draw Rate [kgmol/s]) |
| `adapters/components.py` | PDH 成分キー (A〜F,Z) ↔ HYSYS コンポーネント名のマッピング (部分一致) |

## 設計上の要点

- **段数別 HSC 切替**: COM で段数を変えられないので registry が段数→HSC を解決。要求段数の
  HSC が無ければ `StageNotAvailableError` → penalty。
- **フィード段は最後に書く**: 他入力をすべて済ませた後に `SpecifyFeedLocation` で書かないと
  収束しない (失敗時は DeleteFeedStream + AddFeedStream フォールバック)。
- **ポップアップ自動押下**: ソルバ実行中のモーダル警告で計算が止まらないよう、別スレッドで
  全ウィンドウを poll し HYSYS/Aspen のダイアログを PostMessage で dismiss する。
- **収束判定**: 塔の sub-flowsheet Solver (`col.ColumnFlowsheet.Solver`) の IsSolving +
  Converged で判定 (メイン flowsheet の Solver は塔単独運転時に誤判定するため使わない)。
- **empty sentinel**: HYSYS が Solver 未走/未収束のとき `-32767.0` を返す。`is_hysys_empty`
  で検出して明示的に失敗扱いにし、無音の精度劣化を防ぐ。
- **warm-start 不採用**: case を開きっぱなしにすると前回プロファイルを引きずる warm-start に
  なり BO が経路依存で張り付くため、`force_cold`/`PDH_HYSYS_FORCE_COLD` で cold 解を保証。
  Dist1 メモ化は「入力完全一致時の cold 結果の厳密再利用」であって warm-start ではない。

## パイプライン内の位置づけ・依存

- **依存 (外部)**: `pythoncom`, `win32com`, `win32gui/con/process` (pywin32)、Windows 専用
- **依存 (`src/`)**: `distillation_core` (ColumnTunables/DistResult/DistEquipment),
  `cost_calculator`, `cost_parameters`, `utility_selector`
- **依存 (`flowsheet/`)**: `heat_integration` (StreamPhase/lookup_U/utility_phase)
- **HSC 配置**: `hysys_cases/column{1,2,3}/{N}.hsc` (リポジトリ未同梱の場合あり)
- **呼び出し**: `column*.py` の `solver_method='hysys'` 経路 → `solve_columnX_via_hysys`

`units/vle/__init__.py` / `units/vle/hysys/__init__.py` は本層のパッケージ宣言・公開 API
re-export のみを担う。
