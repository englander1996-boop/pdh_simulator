# units/reactors/ — PDH 反応器モデル

プロパン脱水素 (PDH: C₃H₈ → C₃H₆ + H₂) の**断熱固定床スイング反応器**を 2 つの幾何で
実装する。化学 (反応速度・熱力学・触媒失活・コスト) は共通で、**幾何だけが異なる**
2 つのモデルが並走する。

## 役割

触媒失活を伴う PDH 反応を時間 × 空間で積分し、後段分離工程への時間平均ストリームと
装置コストを出力する。Ergun 圧損を ODE に連成し、ΔP/P が閾値超過なら infeasible 化する。

- 反応: r1 脱水素 (可逆・主反応) / r2 クラッキング / r3 水素化 (副反応)
- 操作: スイング (反応フェーズと触媒再生フェーズを複数基で切替)
- 触媒: Cr₂O₃-Al₂O₃ (Catofin 相当)。コンテスト §3-3 の分単位失活データと物理整合
- 出力: 時間平均出口流量・温度・出口圧 (Ergun 減衰後) + 反応器 CAPEX

## 軸流 (swing) と径方向流 (radial_flow) の関係

| | `swing.py` (軸流) | `radial_flow.py` (径方向流) |
|---|---|---|
| 幾何 | 円筒床を軸方向 z に通す | 環状床を半径方向 r に通す (外周→中心) |
| 断面積 | $A=\pi D^2/4$ (一定) | $A(r)=2\pi r H$ (r 依存) |
| 設計変数 | T_in, z_cat, t_cyc, D | T_in, t_cyc, D_inner, bed_thickness, H |
| 圧損の決まり方 | 床長 z_cat (= 体積と連動) | 床厚 Δr のみ (体積 = H と分離) |
| 0.5 bar での成立性 | 3 mm 粒・深床でほぼ全域 infeasible | 薄床 + 高 H で ΔP が収まる |
| 多段化 | 未実装 (単段等価長さ) | Oleflex 型 N 段直列 + 段間 reheat を実装 |

**径方向流は軸流の幾何だけを差し替えた姉妹実装**で、反応速度・熱力学・Ergun 物性・
失活・コスト・出力データクラスはすべて `swing.py` から import して共有する (二重実装回避)。
したがって同じ V_cat・同じ feed なら転化率/選択率/温度プロファイルはほぼ一致し、差が出る
のは **ΔP だけ** (= 径方向流が「圧損だけ」を救う、という設計主張を担保)。

実機対応: 軸流 ≒ Lummus Catofin (減圧・浅床)、径方向流 ≒ UOP Oleflex (径方向流・多段)。
詳細な背景は両ファイルの冒頭ヘッダコメントを参照。

## ファイル一覧

| ファイル | 説明 |
|---|---|
| `swing.py` | 軸流スイング反応器。`simulate_swing_reactor_system()`。反応速度/熱力学/Ergun 物性/コスト/出力データクラスの**定義元** (径方向流が import する) |
| `radial_flow.py` | 径方向流スイング反応器。単段 `simulate_radial_flow_reactor_system()` + 多段 (Oleflex 型) `simulate_radial_multibed_reactor_system()`。化学は swing から共有 |
| `SPEC_swing.md` | 軸流の詳細仕様 (反応系・速度式・失活・サイジング・コスト・エラー処理) |
| `SPEC_radial_flow.md` | 径方向流の差分仕様 (径方向 ODE・並列分割規約・多段) |
| `swing_timeline.png` | スイング操作タイムライン図 (SPEC_swing.md から参照) |

## パイプライン内の位置づけ・依存

- **上流**: Mixer (Dist1 塔頂 + Dist3/Mem リサイクル) → 反応器入口 (`FeedStream`)
- **下流**: 出口ストリーム → Cooler → Compressor → Dist2。出口圧 `P_out` (Ergun 減衰後)
  を下流圧縮機に伝播し、圧損→圧縮比悪化が TAC に反映される
- **依存 (`src/`)**: `kinetics.PDHKinetics`, `thermo.PDHThermo`, `catalyst_model`,
  `config.THERMO_DATA`, `cost_calculator.calc_reactor_capex_okuyen`,
  `component_data.MW`, `scipy.integrate.solve_ivp` (Radau)
- **REACTOR_KIND 配線**: `run_one_pass.py` 等で軸流/径方向流を切替 (径方向流は 22 変数)

検証ノート: `monitor/reactor_pressure_drop_and_geometry.ipynb`,
`monitor/reactor_conversion_ceiling.ipynb`, `tools/_smoke_test_ergun.py`。
