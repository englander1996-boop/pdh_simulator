# units/utils/ — 補助機器モデル

フローシートの配管をつなぐ汎用補助機器 5 種。いずれも `stream.stream.ProcessStream`
(各成分流量 [kmol/h] + T [K] + P [Pa]) を入出力ベース型とする。

## 役割

昇圧・冷却/加熱・合流・減圧といった基本操作を物理モデルで推算し、出口ストリーム +
動力/熱量 + CAPEX を返す。コスト推算は `src/cost_calculator.py` の Bare Module Cost 法
ヘルパに集約する。

成分マッピング (config.py 準拠): **A=C₃H₈, B=C₃H₆, C=H₂, D=C₂H₄, E=CH₄, F=C₂H₆, Z=C₄H₁₀**

## ファイル一覧

| ファイル | 機器 | 物理モデル | 主な呼び出し箇所 |
|---|---|---|---|
| `pump.py` | 液体ポンプ (遠心式) | $W = V\!\cdot\!\Delta P/\eta$、T 不変 (非圧縮) | Pump1 (Fresh LPG 昇圧) |
| `compressor.py` | 遠心式圧縮機 | ポリトロピック圧縮 ($\gamma$ は流量加重平均) | Comp2a/Comp2b (反応器出口昇圧) |
| `cooler.py` | 冷却器/加熱器 (固定管板式 HE) | 顕熱 + 潜熱 (phase_change 指定時)。U は contest §4-4 表 | Cooler/Intercool/MemPrecool |
| `mixer.py` | ストリームミキサー | 組成加算 + エンタルピー保存で T_out、P は最低値 | 反応器入口 (Dist1塔頂 + リサイクル合流) |
| `expansion_valve.py` | Joule-Thomson 等エンタルピー膨張弁 | $H(T_{in},P_{in})=H(T_{out},P_{out})$ を PR EOS + brentq | Dist1塔頂/Dist3/Mem リサイクルの減圧 |
| `SPEC_utils.md` | 5 機器の詳細仕様 | — | — |

## 各機器の要点

- **pump**: `P_out_target > P_in` 必須 (ValueError)。`eta_pump=0.70`。液密度は
  `component_data.liquid_density_mix`。
- **compressor**: `P_out_target > P_in` 必須。`eta_poly=0.75`。多段化判定は呼び出し側。
- **cooler**: 冷媒/熱媒を `utility_selector` で自動選択し単価を equipment に埋め込む
  (economics.py が直読)。VLE を持たないため相変化は呼び出し側が `phase_change=True` で指示。
  プロセス側顕熱区間の相は `process_phase` で指定 (§4-4 U 表索引用)。
- **mixer**: 旧版のモル流量加重平均からエンタルピーバランス解法へ改修済み。
- **expansion_valve**: vapor 相のまま膨張する仮定 (部分気化は未対応、温度低下を過大評価
  する側に偏る)。配管中の絞り弁扱いで CAPEX/OPEX は計上しない。

## パイプライン内の位置づけ・依存

- **依存 (`src/`)**: `cost_calculator` (CAPEX), `component_data` (MW/Cp/密度/潜熱),
  `utility_selector` (cooler), `eos` (expansion_valve の PR 残差エンタルピー)
- **依存 (`flowsheet/`)**: `heat_integration` (cooler の StreamPhase/lookup_U)
- **依存 (`stream/`)**: `ProcessStream`
- **上位**: `flowsheet/run_one_pass.py` が全機器を配管

`!仮置き` (各効率・$\gamma$ テーブル・$dT_{lm}$ 等) の詳細は `KNOWN_PLACEHOLDERS.md` を参照。
