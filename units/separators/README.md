# units/separators/ — 分離器モデル

PDH プロセスの分離工程を構成する 5 つの分離器: 蒸留塔 3 本 (column1/2/3)・膜分離
(membrane)・圧力スイング吸着 (psa)。

## 役割

反応器出口の多成分混合ガスから、ポリマーグレード C₃H₆ 製品・H₂ 製品を回収し、未反応
プロパンをリサイクルする。各分離器は出口ストリーム + 熱量 (Q) + 装置サイズ/CAPEX を返す。

## 各分離器の役割

```
                      ┌─────────────┐
 Fresh LPG ──Pump1──▶│ Dist1 脱ブタン塔 │
                      └─────────────┘
                       塔頂 C3 → 反応器     塔底 C4H10 → 廃棄
                                  │
                              (反応器ループ)
                                  │ 反応器出口を冷却・圧縮
                                  ▼
                      ┌──────────────┐
                      │ Dist2 脱エタン塔 │
                      └──────────────┘
              塔頂 (H2/CH4/C2)        塔底 (C3H8/C3H6 clean)
                    │                        │
                    ▼                        ▼
            ┌──────────┐            ┌──────────────┐
            │   PSA     │            │  Membrane     │
            └──────────┘            └──────────────┘
        H2 製品 + オフガス(燃料)   C3H6 濃縮 (透過) → Dist3
                                     非透過 → リサイクル
                                              │
                                              ▼
                                  ┌──────────────┐
                                  │ Dist3 C3スプリッタ │
                                  └──────────────┘
                              塔頂 C3H6 製品    塔底 C3H8 → リサイクル
```

| 分離器 | 役割 | 分離原理 |
|---|---|---|
| **column1** (Dist1) | 脱ブタン塔。Fresh LPG から C₄H₁₀ を分離、塔頂 C₃ を反応器へ | 蒸留 (FUG/SM/HYSYS) |
| **column2** (Dist2) | 脱エタン塔。反応器出口から軽質ガス (H₂/CH₄/C₂) を分離 | 蒸留 (partial condenser) |
| **column3** (Dist3) | C₃ スプリッタ。ポリマーグレード C₃H₆ を回収。**最もエネルギー集約的** (α≈1.07) | 蒸留 (高還流・多段) |
| **membrane** | C₃H₆/C₃H₈ 膜分離。Dist2 塔底液を膜で C₃H₆ 濃縮し Dist3 へ | ZIF-8 系膜 (クロスフロー、圧力駆動) |
| **psa** | 圧力スイング吸着。Dist2 塔頂から H₂ を高純度回収 | 活性炭吸着 (Langmuir + LDF) |

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `column1/column1.py` | Dist1 脱ブタン塔。`src/distillation_core` の薄いラッパ (LK='A', HK='Z' 固定) |
| `column2/column2.py` | Dist2 脱エタン塔。partial condenser、LK='F'(C2H6), HK='A'(C3H8) 固定 |
| `column3/column3.py` | Dist3 C3 スプリッタ。LK='B'(C3H6), HK='A'(C3H8)、製品純度から recovery 動的計算 |
| `membrane/membrane_system.py` | 膜分離 5 ユニット (気化器→圧縮機→膜 ODE→圧縮機→冷却器) |
| `membrane/SPEC_membrane_system.md` | 膜分離の詳細仕様 |
| `membrane/ISSUES_membrane_system.md` | 膜分離のレビュー課題・既知の不備メモ |
| `psa/psa_system.py` | PSA 吸着 PDE (1D 上流差分 + LDF) + 脱着 (指数減衰近似) |
| `psa/SPEC_psa_system.md` | PSA の詳細仕様 |
| `psa/error_handling_fixes_20260505.md` | PSA エラー処理の旧バグ修正履歴 |
| `SPEC_columns.md` | 蒸留塔 3 本の詳細仕様 (FUG/SM/HYSYS 3 経路) |

## 蒸留塔の VLE バックエンド切替

column1/2/3 はいずれも `ColumnTunables.solver_method` で計算経路を切替える薄いラッパ:

| solver_method | 経路 | 用途 |
|---|---|---|
| `'fug'` (既定) | `src/distillation_core` Fenske-Underwood-Gilliland shortcut | 高速、BO 用 |
| `'rigorous'` | `src/distillation_core` VLE tray-by-tray | 厳密、top-k 再評価用 |
| `'sm'` | `src/distillation_sm` 学習済み GPR で HYSYS 解を近似 | SM フォーク |
| `'hysys'` | `units/vle/hysys/` 経由で HYSYS COM 実行 | 検証・special フォーク |

LK/HK・回収率・q・K_method・partial_condenser といった塔別の物理セマンティクスは
各ラッパ (column*.py) で固定し、`ColumnTunables` は BO/exp で振る部分 (P_col, N_stages,
N_feed, reflux_ratio + HYSYS 用追加) のみ保持する。

## パイプライン内の位置づけ・依存

- **依存 (`src/`)**: `distillation_core` (3 塔共通エンジン), `eos` (PR EOS: 膜・PSA・PR 蒸留),
  `thermo`, `cost_calculator`, `cost_parameters`, `component_data`, `config`
- **依存 (`units/`)**: 蒸留塔の HYSYS 経路は `units/vle/hysys/` を呼ぶ
- **上位**: `flowsheet/run_one_pass.py` が 5 分離器を配管、`economics.py` が OPEX 計上
