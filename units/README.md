# units/ — 単位操作モデル群

PDH (プロパン脱水素) プロセスシミュレータの**単位操作 (unit operation) モデル**を
集約したサブツリー。反応器・分離器・補助機器・気液平衡 (VLE) バックエンドを、それぞれ
独立した Python モジュールとして実装する。

## 役割

各単位操作を「設計変数 (DesignVars) + 入口ストリーム (FeedStream) + 固定パラメータ
(FixedParams) を受け取り、出口ストリーム + 装置サイズ/コストを返す純粋関数」として
モデル化する。フローシート (`flowsheet/run_one_pass.py`) はこれらを配管して 1 パスの
プロセス全体を構成し、最適化 (`main.py` 等) が設計変数を振る。

- 無効入力・数値異常は例外を投げず、大きな CAPEX ペナルティ (`1e9` 億円) を返して
  最適化器に「探索しない領域」を伝える方針 (各モジュール共通)。
- コスト推算は `src/cost_calculator.py` の Bare Module Cost 法ヘルパに集約。
- 成分マッピング (全モジュール共通、`src/config.py` 準拠):
  **A=C₃H₈, B=C₃H₆, C=H₂, D=C₂H₄, E=CH₄, F=C₂H₆, Z=C₄H₁₀**

## サブディレクトリ一覧

| ディレクトリ | 役割 | 詳細 README |
|---|---|---|
| `reactors/` | PDH 反応器 (軸流 swing / 径方向流 radial_flow) | [reactors/README.md](reactors/README.md) |
| `separators/` | 蒸留塔 3 本 (column1/2/3) + 膜分離 + PSA | [separators/README.md](separators/README.md) |
| `utils/` | 補助機器 (ポンプ・圧縮機・冷却器・ミキサー・膨張弁) | [utils/README.md](utils/README.md) |
| `vle/` | 気液平衡バックエンド (HYSYS COM 連携層) | [vle/hysys/README.md](vle/hysys/README.md) |

## パイプライン内の位置づけ

```
Fresh LPG
  │ Pump1 (utils)
  ▼
Dist1 (separators/column1)  脱ブタン塔 → 塔頂 C3 を反応器へ
  │ ExpansionValve (utils) で 0.5 bar へ
  ▼
Mixer (utils): Dist1塔頂 + Recycle(Dist3) + Recycle(Mem)
  ▼
Reactor (reactors): 軸流 swing または 径方向流 radial_flow
  │ Cooler + Compressor (utils) で昇圧
  ▼
Dist2 (separators/column2)  脱エタン塔
  ├─ 塔頂 (軽質) → PSA (separators/psa) → H2 製品 + オフガス
  └─ 塔底 (C3) → Membrane (separators/membrane) → Dist3
                                                    ▼
                          Dist3 (separators/column3) C3 スプリッタ → C3H6 製品
```

蒸留塔 (column1/2/3) は VLE バックエンドを切替可能で、`solver_method` に応じて
`src/distillation_core.py` の FUG/rigorous/SM 経路、または `vle/hysys/` の HYSYS COM
経路にディスパッチする。

## 依存関係

- 上位 (呼び出し側): `flowsheet/`, `main.py`, `exp/`, `monitor/`, `comparing/`
- 下位 (本サブツリーが依存): `src/` (kinetics, thermo, eos, cost_calculator,
  component_data, config, distillation_core 等), `stream/` (ProcessStream)

## ドキュメント

各単位操作の詳細仕様は同ディレクトリ内の `SPEC_*.md` を参照:

| SPEC | 対象 |
|---|---|
| `reactors/SPEC_swing.md` | 軸流スイング反応器 |
| `reactors/SPEC_radial_flow.md` | 径方向流スイング反応器 (単段 + Oleflex 型多段) |
| `separators/SPEC_columns.md` | 蒸留塔 3 本 (FUG / SM / HYSYS の 3 経路) |
| `separators/membrane/SPEC_membrane_system.md` | 膜分離システム |
| `separators/psa/SPEC_psa_system.md` | PSA (圧力スイング吸着) |
| `utils/SPEC_utils.md` | 補助機器 5 種 |

`!仮置き` / `【確認中】` マーカーの詳細一覧はリポジトリルートの
`KNOWN_PLACEHOLDERS.md` を参照。
