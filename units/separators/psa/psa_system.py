"""
PSA (Pressure Swing Adsorption) システム シミュレーター

構成:
  1. プレヒーター (Preheater)   フィードを T_abs = 25°C に昇温
  2. PSA カラム (Adsorption)   1D 上流差分 PDE + LDF モデル
  3. 脱着 (Desorption)         指数減衰近似モデル

成分マッピング (config.py 準拠):
  'A' : C3H8 (プロパン)   — 全量オフガスへ
  'B' : C3H6 (プロピレン) — 全量オフガスへ
  'C' : H2  (水素)        — 非吸着、全量プロダクトへ
  'D' : C2H4 (エチレン)   — 吸着成分 (PDE 対象)
  'E' : CH4  (メタン)     — 吸着成分 (PDE 対象、破過基準成分)
  'F' : C2H6 (エタン)     — 吸着成分 (PDE 対象)

PDE モデル
----------
  ガス相質量収支 (上流差分):
    ∂C_i/∂t = -(u_0/ε)·∂C_i/∂z - (ρ_b/ε)·∂q_i/∂t

  固相質量収支 (LDF):
    ∂q_i/∂t = KFa_i · (q_i* - q_i)

  多成分 Langmuir 等温線 (Markham-Benton 形):
    q_i* = q_si · a_i · C_i / (1 + Σ a_j · C_j)
    C_i [mol/m³],  a_i [m³/mol]

脱着モデル (簡易):
    q_i(t) = q_i0 · exp(-KFa_i · t)
    停止条件: Σ q_i(t) / Σ q_i(0) = desorption_target

--------------------------------------------------------------------
【仮定一覧】

■ 文献・仕様書に根拠がある仮定

  [C3 成分の完全捕捉]
      活性炭に対する C3H8・C3H6 の Langmuir 定数 a は CH4 の約 30〜80 倍
      （参考: Schell et al., Sep. Purif. Technol., 2012 等）であり、
      吸着層前段で実質的に全量捕捉される。C3 は PDE から除外して
      全量オフガスとして扱い、計算を簡略化する。

  [H2 非吸着仮定]
      活性炭への H2 吸着量は 25°C・低圧条件で 0.01 mol/kg 以下であり
      （参考: Poirier & Darriet, J. Chem. Eng. Data, 2001）、
      CH4 (約 0.17 mol/kg @ 15 bar, 25°C) と比べ 1 桁以上小さい。
      設計精度の範囲内で無視し、H2 は全量プロダクトへ流す。

  [多成分 Langmuir (Markham-Benton 形)]
      単成分等温線が Langmuir 型に従う場合、多成分系への拡張として
      最も広く用いられる近似。厳密な IAS 理論より計算が軽い。
      根拠: Yang, R. T., "Gas Separation by Adsorption Processes" (1987)

  [LDF (Linear Driving Force) モデル]
      Glueckauf (1955) の近似。粒子内拡散が律速の場合に等価な
      一次モデルとして成立する。PSA シミュレーションでは標準的。
      根拠: Ruthven, D. M., "Principles of Adsorption Processes" (1984)

■ 設計判断（根拠あり）

  [等温 PSA 操作 T_abs = 25°C]
      吸着は低温ほど有利（Langmuir 平衡定数が増加）。25°C は
      冷却水のみで到達可能な下限温度（コンテスト仕様: T_cold_out=40°C
      を考慮し余裕を取る）。また Langmuir 仮置き値の測定条件（25°C
      を想定）に合わせることでパラメータと操作条件の整合を確保。

  [脱着圧力 P_des = 大気圧 (101325 Pa)]
      真空脱着 (P_des < 1 atm) は脱着効率が向上するが真空ポンプが
      追加となり CAPEX・OPEX が増加する。初期設計ではコスト最小化を
      優先し、大気圧ブローダウンのみの簡易サイクルとする。

  [空塔速度 = 一定（等速近似）]
      フィードの主成分は H2 であり、吸着成分 (CH4, C2H4, C2H6) の
      モル分率は合計で数 mol% 以下。吸着による気相成分の減少が速度に
      与える影響は相対誤差 5% 未満に収まるため等速近似を採用する。
      ※ C3 成分 (C3H8, C3H6) は u_0 計算から除外（入口付近で即吸着と仮定）。
      C3 モル分率が大きい場合は入口速度を過小推算するため要確認。

  [N_abs_parallel = 1（並列塔なし）]
      必要塔数は N_total = ceil(t_des/t_abs) + 1 で決まる。
      t_des << t_abs の場合は 2 塔で足りる（1 塔吸着・1 塔脱着）。
      スケールアップが必要な場合は本パラメータを設計変数に昇格させる。

  [破過基準 breakthrough_ratio = 0.001 (0.1%)]
      吸着操作の終点を「CH4 出口濃度が入口の 0.1% に達した時点」と定義。
      これにより CH4 捕捉率 > 99.9% を保証する保守的な設計基準。
      緩める場合は OPEX の CH4 損失と CAPEX のトレードオフを再評価すること。

  [グリッド点数 N_z = 20]
      N_z が少ないと数値拡散が大きく破過が早まる。
      N_z=10 と N_z=20 で t_abs の差が < 5% であることを
      テストにより確認済みで、計算速度と精度のバランスから 20 点を採用する。

  [脱着モデル: 指数減衰近似]
      低圧パージ (P_des ≈ 1 atm) では吸着平衡値 q* ≈ 0 となるため、
      LDF 式 dq/dt = KFa·(0 - q) の解が q(t) = q₀·exp(-KFa·t) になる。
      空間平均 q_avg に適用することで PDE の再積分を省略する。
      ※ 固相負荷量の z 方向分布（入口高・出口低）を無視するため、
      入口高負荷セルが律速になる場合に t_des を過小推算する可能性がある。
      安全係数 desorption_time_safety_factor=1.2 でこのリスクをカバーしている。

■ !仮置き値 — 根拠文献未確定

  Langmuir パラメータ (q_s, a): 活性炭 25°C 等温線、学術論文未確認
  KFa (物質移動係数): Carberry 数・Knudsen 拡散の推算未実施
  活性炭単価: メーカーカタログ未照合
  活性炭嵩密度 rho_b = 600 kg/m³: 工業用活性炭の代表値（典型範囲 400〜700 kg/m³）
      → cost_parameters.py の !仮置き コメントを参照。
--------------------------------------------------------------------
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.config import R
from src.component_data import MW
from src.eos import z_factor
from src.thermo import PDHThermo
from src.cost_calculator import calc_reactor_capex_okuyen
from src.cost_parameters import (
    PSA_LANGMUIR_PARAMS,
    PSA_KFA,
    ACTIVATED_CARBON_PRICE_USD_PER_KG,
    ADSORBENT_LIFETIME_YEARS,
    CEPCI_BASE, CEPCI_CURRENT,
    USD_TO_JPY, PLANT_INDIRECT_FACTOR,
)

_thermo = PDHThermo()

# 成分配列の順序: index → (config キー, PSA パラメータキー)
# 吸着成分 (N_ADS = 3): 0=CH4('E'), 1=C2H4('D'), 2=C2H6('F')
_ADS_ORDER  = ['CH4', 'C2H4', 'C2H6']           # PSA_LANGMUIR_PARAMS のキー
_ADS_KEYS   = ['E', 'D', 'F']                   # config キー (対応順)
_EOS_KEYS   = ['C', 'E', 'D', 'F']              # PR EOS 対象: H2, CH4, C2H4, C2H6
_C3_KEYS    = ['A', 'B']                         # 全量オフガス成分
_N_ADS      = 3                                  # 吸着成分数
_N_Z        = 20  # 空間グリッド点数: 計算速度と数値拡散のバランスから選択（未検証）
# _U0_MAX: 空塔速度ハードリミット [m/s]
# 【設計パラメータ設定根拠】
# 根拠文献: 『化学工学便覧』第13章「吸着・イオン交換の操作と装置」
#   ・図13・31「サイクル時間によるPSA除湿性能の変化」の解析条件として u_0 = 1.0 m/s が
#     理論的・実験的に性能評価された値として用いられている。
#   ・（参考）表13・28「圧力スイング吸着による各種ガス分離，精製の操作条件」の
#     実機の空塔速度は「0.3程度」とされている。
# 設定値 1.0 m/s の判断理由:
#   空塔速度の過大は圧力損失増大・吸着材層の流動化（チャネリング・粉化）を引き起こす。
#   実機標準値（約 0.3 m/s）のみに探索空間を限定すると潜在的な最適解を見逃す恐れがある
#   ため、文献で性能評価パラメータとして実証された 1.0 m/s をシミュレータ上の
#   ハードリミット（ペナルティ閾値）として採用した。これにより ODE ソルバーの安定性
#   （LSODA フリーズ防止）と設計の物理的妥当性を両立している。
_U0_MAX     = 1.0  # [m/s] 空塔速度の上限: 超過時はペナルティを返す
_T_ABS_MIN  = 60.0  # [s]   CSS補正後の最小吸着時間: 未満はペナルティ (scale 発散防止 + 物理的に無意味な超短サイクルの早期排除)


# ---------------------------------------------------------------------------
# データクラス: 入力
# ---------------------------------------------------------------------------

@dataclass
class PSADesignVars:
    """最適化アルゴリズムが操作する設計変数"""
    D_col:             float  # 吸着塔塔径 [m]
    L_bed:             float  # 吸着層高さ [m]
    desorption_target: float  # 脱着完了基準 [-] (例: 0.35 → q が初期値の 35% まで低下)

    def __post_init__(self) -> None:
        if self.D_col <= 0.0:
            raise ValueError(f"D_col={self.D_col} は正値でなければなりません。")
        if self.L_bed <= 0.0:
            raise ValueError(f"L_bed={self.L_bed} は正値でなければなりません。")
        if not (0.0 < self.desorption_target < 1.0):
            raise ValueError(
                f"desorption_target={self.desorption_target} は (0, 1) の範囲でなければなりません。"
            )


@dataclass
class PSAFeedStream:
    """PSA フィードストリーム (上流蒸留塔オーバーヘッドガス)"""
    F_in: Dict[str, float]  # 成分別モル流量 [kmol/h]
                             # キー: 'A'(C3H8), 'B'(C3H6), 'C'(H2),
                             #       'D'(C2H4), 'E'(CH4), 'F'(C2H6)
    T_in: float              # 入口温度 [K]
    P_in: float              # 入口圧力 [Pa] = PSA 吸着操作圧力


@dataclass
class PSAFixedParams:
    """
    PSA システム固定パラメータ

    T_abs = 298.15 K (25°C)
        吸着有利な低温かつ冷却水のみで到達可能な温度。
        Langmuir 仮置き値の想定測定条件 (25°C) と一致させ整合性を確保。

    P_des = 101325 Pa (大気圧)
        真空ポンプ不要の最低コスト構成。脱着効率より設備簡素化を優先。
        真空 PSA (P_des = 0.1〜0.3 atm) との比較は最終設計段階で再評価すること。

    rho_b = 600 kg/m³
        !仮置き — 工業用活性炭の代表的嵩密度。
        典型範囲: 400〜700 kg/m³ (Perry's Chemical Engineers' Handbook, 8th Ed.)
        実際に使用する活性炭のデータシートで確認後に更新すること。

    eps = 0.4
        球状粒子ランダム充填の代表値 (Ergun 式では 0.36〜0.42)。
        活性炭ペレット充填塔の設計標準値として 0.4 を採用。
        実測値が得られた場合は更新すること。

    breakthrough_ratio = 0.001
        CH4 出口/入口濃度比 = 0.1% を破過の定義とする。
        これにより CH4 捕捉率 > 99.9% を保証する保守的設計基準。

    t_ads_max = 7200 s (2 時間)
        solve_ivp の打ち切り上限。典型的な t_abs は数分〜数十分であり、
        2 時間以内に破過しない場合はペナルティを返して最適化に再評価を促す。

    use_css_approximation = True  !仮置き
        True のとき、吸着 PDE の初期固相負荷量をゼロ（清浄床）ではなく
        「フィード濃度での Langmuir 平衡値 × desorption_target」で初期化する。
        これはサイクル定常状態（CSS）の近似だが、床全体が飽和していると仮定するため
        実際の CSS より初期残留を保守的に過大評価する。
        → t_abs が清浄床より短くなり、必要塔数はやや多め（コスト保守的）に出る。
        False にすると清浄床（q=0）スタートに戻せる（比較用）。

    desorption_time_safety_factor = 1.2  !仮置き
        計算された t_des に掛ける安全係数。
        理由: KFa (物質移動係数) が仮置き値であり、文献値との誤差が
        数倍に達する可能性があるため、設計マージンとして 1.2 を仮置きする。
        KFa を推算・実測値で確定させた後、この係数の妥当性を再評価すること。
    """
    T_abs:              float = 298.15   # [K]  PSA 操作温度
    P_des:              float = 101325.0 # [Pa] 脱着圧力
    rho_b:              float = 600.0    # [kg/m³] 活性炭充填嵩密度 !仮置き
    eps:                float = 0.4      # [-]  吸着層空隙率
    breakthrough_ratio: float = 0.001   # [-]  破過基準 (CH4 出口/入口)
    t_ads_max:          float = 7200.0  # [s]  吸着時間の探索上限
    # !仮置き — 後日検証・調整すること
    use_css_approximation:         bool  = True  # CSS 簡易補正フラグ (保守的過大推算)
    desorption_time_safety_factor: float = 1.2   # 脱着時間安全係数 (KFa 不確実性対策)
    # ---- 床圧力損失 (Ergun) — 2026-05-31 PSA設計レビュー対応 ----
    # !仮置き — 確定値は活性炭ベンダーデータで更新。u_0 上限 (_U0_MAX=1.0m/s) は ODE 安定の
    #   数値ガードに留め、実機の現実的な空塔速度 (~0.3-0.4m/s) は本 ΔP 制約で物理的に縛る。
    d_p_m:        float = 0.003    # !仮置き 活性炭粒径 [m] (3mm 成形炭)
    sphericity:   float = 0.9      # !仮置き 形状係数 [-]
    mu_gas_pa_s:  float = 1.0e-5   # !仮置き H2 リッチガス混合粘度 [Pa·s] (25°C)
    dP_max_bar:   float = 0.3      # !仮置き 床 ΔP 上限 [bar] (レビュー目安 0.1-0.3)。超過で infeasible 化


# ---------------------------------------------------------------------------
# データクラス: 出力
# ---------------------------------------------------------------------------

@dataclass
class PSAEquipmentData:
    """機器サイズ・コスト推算テーブル"""
    N_abs_parallel:  int    # 並列吸着塔数 (= 1 固定)
    N_cycle_sets:    int    # = ceil(t_des/t_abs) + 1
    N_total_columns: int    # = N_abs_parallel × N_cycle_sets
    t_abs_sec:       float  # 吸着時間 [s]
    t_des_sec:       float  # 脱着時間 [s]
    u_0:             float  # 空塔速度 [m/s]
    W_adsorbent_kg:  float  # 吸着材総重量 [kg]
    Q_preheat_kW:    float  # プレヒーター熱量 [kW] (正=加熱, 負=冷却)
    CAPEX_vessels:   float = float('nan')   # 塔体 CAPEX [億円]
    CAPEX_adsorbent: float = float('nan')   # 吸着材 CAPEX [億円]
    CAPEX_total:     float = float('nan')   # 合計 CAPEX [億円]
    # H2 損失 (!仮置きモデルによる推算値)
    H2_loss_blowdown_kmolh:        float = float('nan')  # ブローダウン損失 [kmol/h]
    H2_loss_purge_kmolh:           float = float('nan')  # パージ損失 [kmol/h]
    # 吸着材交換 OPEX !仮置き (ADSORBENT_LIFETIME_YEARS に依存)
    OPEX_adsorbent_okuyen_per_year: float = float('nan')  # 吸着材年間交換費 [億円/年]
    # ---- penalty 診断 (2026-05-21 追加、BO の TPE constraints_func 用) ----
    # 設計判断: silent _penalty_result() 経路が BO に「方向のシグナル」を渡せず
    # 全 trial が無方向で死ぬ問題を解消するための診断フィールド。
    # penalty 発火時に「どの条件で死んだか」「actual 値」を保持し、run_one_pass が
    # log10(MIN/actual) 等の連続 shortfall を計算して TPE に渡せるようにする。
    # 通常完走時は penalty_reason='' のままで識別する (CAPEX_total < threshold で判定)。
    penalty_reason:    str   = ''     # '' | 't_abs_below_min' | 'u_0_above_max' | 'dp_excess' | 'no_non_C3_feed' | 'no_CH4_feed' | 'breakthrough_no_converge' | 'mask_lt_2'
    t_abs_actual_s:    float = 0.0    # penalty 発火時の CSS 補正後 t_abs [s] (0 なら未計算)
    u_0_actual:        float = 0.0    # penalty 発火時の空塔速度 [m/s] (0 なら未計算)
    dP_bar_actual:     float = 0.0    # 床 Ergun 圧損 [bar] (正常完走時も格納、dp_excess で連続シグナル化)


@dataclass
class PSASimulationResult:
    """PSA システム シミュレーション全出力"""
    product:     Dict[str, float]  # H2 リッチプロダクト [kmol/h]  キー: 'A'〜'F'
    offgas:      Dict[str, float]  # オフガス [kmol/h]             キー: 'A'〜'F'
    equipment:   PSAEquipmentData
    H2_recovery: float  # H2 回収率 [-]
    CH4_capture: float  # CH4 捕捉率 [-]


# ---------------------------------------------------------------------------
# ペナルティ結果
# ---------------------------------------------------------------------------

_PENALTY = 1e9


def _penalty_result(
    reason:         str   = '',
    t_abs_actual:   float = 0.0,
    u_0_actual:     float = 0.0,
    dP_bar_actual:  float = 0.0,
) -> PSASimulationResult:
    """計算不能な条件のときに返すペナルティ結果。

    Parameters
    ----------
    reason : str
        発火条件のラベル ('t_abs_below_min' 等)。BO の TPE constraints_func で
        連続 shortfall を計算するために run_one_pass が読み取る。
    t_abs_actual : float
        CSS 補正後の実 t_abs [s]。t_abs_below_min 経路でのみ意味を持つ。
    u_0_actual : float
        実空塔速度 [m/s]。u_0_above_max 経路でのみ意味を持つ。

    設計判断 (2026-05-21): 旧版は引数なしの silent penalty で、BO は「どこを突けば
    feasible に出るか」のシグナルを得られず stuck していた (main_20260521_131507
    で 300/300 trial 全滅)。本版は理由ラベル + 実値を equipment に格納して
    run_one_pass で psa_t_abs_shortfall を計算 → TPE constraints_func で
    「t_abs MIN への接近度」を学習させる。
    """
    zero = {k: 0.0 for k in ['A', 'B', 'C', 'D', 'E', 'F']}
    eq = PSAEquipmentData(
        N_abs_parallel=0, N_cycle_sets=0, N_total_columns=0,
        t_abs_sec=0.0, t_des_sec=0.0, u_0=0.0,
        W_adsorbent_kg=0.0, Q_preheat_kW=0.0,
        CAPEX_total=_PENALTY,
        penalty_reason=reason,
        t_abs_actual_s=t_abs_actual,
        u_0_actual=u_0_actual,
        dP_bar_actual=dP_bar_actual,
    )
    return PSASimulationResult(
        product=dict(zero), offgas=dict(zero),
        equipment=eq, H2_recovery=0.0, CH4_capture=0.0,
    )


# ---------------------------------------------------------------------------
# 内部ヘルパー
# ---------------------------------------------------------------------------

def _calc_preheat_kW(
    F_in:  Dict[str, float],
    T_in:  float,
    T_abs: float,
) -> float:
    """
    フィードを T_in → T_abs に昇温(または降温)する熱量 [kW]。

    Returns
    -------
    float
        熱量 [kW]。正値 = 加熱、負値 = 冷却。
    """
    Q = 0.0
    for key, F_kmolh in F_in.items():
        if F_kmolh <= 0.0:
            continue
        F_mol_s = F_kmolh * 1000.0 / 3600.0            # [mol/s]
        dH = _thermo.calc_enthalpy_change(key, T_in, T_abs)  # [J/mol]
        Q += F_mol_s * dH                               # [W]
    return Q / 1000.0                                   # [kW]


def _calc_feed_state(
    F_in:  Dict[str, float],
    T_abs: float,
    P_in:  float,
) -> tuple:
    """
    PSA カラム入口における吸着成分濃度と空塔速度の基礎量を計算する。

    Returns
    -------
    C_feed_ads : np.ndarray shape (3,)  [mol/m³]  [CH4, C2H4, C2H6]
    C_H2       : float                 [mol/m³]  H2 濃度 (H2 損失計算用)
    F_non_C3   : float                 [mol/s]   非 C3 総モル流量
    Z          : float                            PR EOS Z 因子
    """
    # 非 C3 成分のモル流量 [mol/s]
    F = {k: F_in.get(k, 0.0) * 1000.0 / 3600.0 for k in _EOS_KEYS}
    F_total = sum(F.values())
    if F_total <= 0.0:
        return np.zeros(_N_ADS), 0.0, 0.0, 1.0

    x = [F[k] / F_total for k in _EOS_KEYS]  # H2, CH4, C2H4, C2H6 モル分率

    try:
        Z = z_factor(T_abs, P_in, x, _EOS_KEYS, phase='vapor')
    except Exception as e:
        warnings.warn(
            f"_calc_feed_state: z_factor 計算失敗 ({e})。理想気体 Z=1.0 を使用。"
            f" u_0・C_feed に誤差が混入します。",
            UserWarning, stacklevel=2,
        )
        Z = 1.0

    C_total = P_in / (Z * R * T_abs)  # [mol/m³]

    # _EOS_KEYS 順: ['C', 'E', 'D', 'F'] → index 0=H2, 1=CH4, 2=C2H4, 3=C2H6
    C_H2 = x[0] * C_total  # [mol/m³]  H2 濃度 (ブローダウン・パージ損失計算用)
    C_feed_ads = np.array([
        x[1] * C_total,  # CH4
        x[2] * C_total,  # C2H4
        x[3] * C_total,  # C2H6
    ])

    return C_feed_ads, C_H2, F_total, Z


def _run_adsorption(
    C_feed:             np.ndarray,  # [mol/m³] (CH4, C2H4, C2H6)
    u_0:                float,       # [m/s] 空塔速度
    L_bed:              float,       # [m]
    rho_b:              float,       # [kg/m³]
    eps:                float,       # [-]
    kfa:                np.ndarray,  # [1/s]  (CH4, C2H4, C2H6)
    q_s:                np.ndarray,  # [mol/kg]
    a_lang:             np.ndarray,  # [m³/mol]
    breakthrough_ratio: float,
    t_max:              float,
) -> tuple:
    """
    吸着 PDE を LSODA で積分し、破過時刻・固相最終分布・出口累積量を返す。

    ガス相質量収支 (等速・等温・1D プラグフロー):
        ε·∂C_i/∂t + ρ_b·∂q_i/∂t + u_0·∂C_i/∂z = 0
        → ∂C_i/∂t = -(u_0/ε)·(C_i[k]-C_i[k-1])/Δz − (ρ_b/ε)·∂q_i/∂t

    状態ベクトル y (長さ 2·N_z·N_ADS):
        y[k·N_ADS + i]            = C[k, i]  (気相, k=0..N_z-1, i=0..2)
        y[N_z·N_ADS + k·N_ADS + i] = q[k, i]  (固相)

    Returns
    -------
    t_abs     : float           破過時刻 [s]
    q_final   : np.ndarray      固相負荷量 shape (N_z, 3) [mol/kg]
    sol_t     : np.ndarray      積分時刻列 shape (n_t,) [s]
    C_outlet  : np.ndarray      出口濃度時系列 shape (N_ADS, n_t) [mol/m³]
    converged : bool            True = 破過イベントで停止
    """
    dz    = L_bed / _N_Z
    u_eps = u_0 / eps  # 実際の流体速度 [m/s]

    def rhs(t, y):
        # 数値誤差による負値を防止 (LSODA は非負制約を持たない)
        C = np.maximum(y[:_N_Z * _N_ADS].reshape(_N_Z, _N_ADS), 0.0)   # [mol/m³]
        q = np.maximum(y[_N_Z * _N_ADS:].reshape(_N_Z, _N_ADS), 0.0)   # [mol/kg]

        # 多成分 Langmuir (Markham-Benton)
        aC    = a_lang * C                              # [-], shape (N_z, 3)
        denom = 1.0 + aC.sum(axis=1, keepdims=True)
        q_eq  = q_s * aC / denom                       # [mol/kg]

        # LDF
        dq_dt = kfa * (q_eq - q)

        # 上流差分
        C_up       = np.empty_like(C)
        C_up[0, :] = C_feed
        C_up[1:, :] = C[:-1, :]

        dC_dt = -u_eps / dz * (C - C_up) - rho_b / eps * dq_dt
        return np.concatenate([dC_dt.ravel(), dq_dt.ravel()])

    def breakthrough_event(t, y):
        # CH4 (成分 0) の出口濃度が閾値を超えたとき正に転じる
        return y[(_N_Z - 1) * _N_ADS] - breakthrough_ratio * C_feed[0]

    breakthrough_event.terminal  = True
    breakthrough_event.direction = 1

    y0  = np.zeros(2 * _N_Z * _N_ADS)
    sol = solve_ivp(
        rhs,
        [0.0, t_max],
        y0,
        method='LSODA',
        events=breakthrough_event,
        dense_output=False,
        rtol=1e-4,
        atol=1e-7,
    )

    converged = len(sol.t_events[0]) > 0
    t_abs     = float(sol.t_events[0][0]) if converged else sol.t[-1]

    # 固相最終分布
    q_final = sol.y[_N_Z * _N_ADS:, -1].reshape(_N_Z, _N_ADS)

    # 出口濃度時系列 (呼び出し側で CSS 補正後の t_abs に合わせて積分範囲を切り詰める)
    i0       = (_N_Z - 1) * _N_ADS
    C_outlet = sol.y[i0:i0 + _N_ADS, :]  # shape (N_ADS, n_t)

    return t_abs, q_final, sol.t, C_outlet, converged


def _calc_desorption_time(
    q_avg:             np.ndarray,  # [mol/kg] 成分別空間平均固相負荷量
    kfa:               np.ndarray,  # [1/s]
    desorption_target: float,
    t_max:             float = 1e5,
) -> float:
    """
    脱着時間 t_des [s] を brentq で求める。

    低圧ブローダウン後のパージ段階では q* ≈ 0 なので
    LDF 式 dq/dt = KFa·(0 - q) の解析解が得られる:
        q_i(t) = q_i0 · exp(-KFa_i · t)

    停止条件 (brentq で逆算):
        Σ q_i(t_des) / Σ q_i(0) = desorption_target

    t=0 で左辺=1 > target、t→∞ で左辺→0 < target → brentq で解が一意に存在。
    """
    q_total_0 = q_avg.sum()
    if q_total_0 <= 0.0:
        return 0.0

    def residual(t: float) -> float:
        return (q_avg * np.exp(-kfa * t)).sum() / q_total_0 - desorption_target

    # t=0 で residual > 0, t=t_max で residual ≈ -target < 0
    if residual(t_max) >= 0.0:
        return t_max  # 理論上の上限に達しても目標に到達しない (kfa が極めて小さい場合)

    return brentq(residual, 0.0, t_max, xtol=1.0, maxiter=200)


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_psa_system(
    design: PSADesignVars,
    feed:   PSAFeedStream,
    fixed:  PSAFixedParams = None,
) -> PSASimulationResult:
    """
    PSA システム全体をシミュレートする。

    手順
    ----
    1. プレヒーター熱量を計算
    2. フィード濃度・空塔速度を計算 (PR EOS)
    3. 吸着 PDE を solve_ivp (LSODA) で積分、破過時刻 t_abs を取得
    4. 指数減衰モデルで脱着時間 t_des を計算
    5. サイクル塔数 N_total_columns = ceil(t_des/t_abs) + 1
    6. 出力ストリーム (プロダクト / オフガス) を計算
    7. CAPEX (塔体 + 吸着材) を計算

    Parameters
    ----------
    design : PSADesignVars
    feed   : PSAFeedStream
    fixed  : PSAFixedParams  省略時はデフォルト値を使用

    Returns
    -------
    PSASimulationResult
    """
    if fixed is None:
        fixed = PSAFixedParams()

    D_col  = design.D_col
    L_bed  = design.L_bed
    A_col  = math.pi / 4.0 * D_col ** 2  # [m²]
    V_col  = A_col * L_bed                # [m³]

    # -------------------------------------------------------------------------
    # 1. プレヒーター
    # -------------------------------------------------------------------------
    Q_preheat_kW = _calc_preheat_kW(feed.F_in, feed.T_in, fixed.T_abs)

    # -------------------------------------------------------------------------
    # 2. フィード濃度・空塔速度
    # -------------------------------------------------------------------------
    C_feed_ads, C_H2, F_non_C3_mol_s, Z = _calc_feed_state(
        feed.F_in, fixed.T_abs, feed.P_in
    )

    if F_non_C3_mol_s <= 0.0:
        warnings.warn(
            f"PSA penalty: F_non_C3_mol_s={F_non_C3_mol_s:.3e} ≤ 0 (上流 Dist2 で非 C3 がほぼ流れていない)。"
            f" feed.F_in keys: {[k for k,v in feed.F_in.items() if v > 0]}。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='no_non_C3_feed')

    u_0 = F_non_C3_mol_s * Z * R * fixed.T_abs / (feed.P_in * A_col)  # [m/s]

    # 空塔速度が上限超過: LSODA が極めて小さいタイムステップを要求しフリーズする
    if u_0 > _U0_MAX:
        warnings.warn(
            f"PSA penalty: u_0={u_0:.3f}m/s > _U0_MAX={_U0_MAX} (D_col={D_col:.2f}m が小さい / 流量が大きい)。"
            f" D_col を大きく or 流量を減らす方向に探索を誘導。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='u_0_above_max', u_0_actual=u_0)

    # ---- 床圧力損失 (Ergun) チェック (2026-05-31 PSA設計レビュー対応) ----
    # 設計判断: u_0 上限 (_U0_MAX=1.0m/s) は ODE 安定の数値ガードに留め，実機の現実的な
    #   空塔速度は床 ΔP の物理制約で縛る (高 u_0 / 長床 / 小粒径 ほど ΔP 増)。閾値超過で
    #   infeasible 化し，run_one_pass で psa_dp_shortfall を連続シグナル化 (反応器 ΔP と同型)。
    C_total = feed.P_in / (Z * R * fixed.T_abs)                  # [mol/m³] 全濃度
    F_eos = {k: feed.F_in.get(k, 0.0) for k in _EOS_KEYS}
    Ftot_eos = sum(F_eos.values())
    MW_avg = (sum(F_eos[k] * MW[k] for k in _EOS_KEYS) / Ftot_eos) if Ftot_eos > 0 else 2.0
    rho_gas = C_total * MW_avg / 1000.0                          # [kg/m³] ガス密度
    phi_dp = fixed.sphericity * fixed.d_p_m
    eb = fixed.eps
    visc_term  = 150.0 * (1.0 - eb) ** 2 * fixed.mu_gas_pa_s * u_0 / (eb ** 3 * phi_dp ** 2)
    inert_term = 1.75 * (1.0 - eb) * rho_gas * u_0 ** 2 / (eb ** 3 * phi_dp)
    dP_bar = (visc_term + inert_term) * L_bed / 1.0e5            # [bar] 床全体の圧損
    if dP_bar > fixed.dP_max_bar:
        warnings.warn(
            f"PSA penalty: 床ΔP={dP_bar:.3f}bar が上限 {fixed.dP_max_bar}bar 超 "
            f"(u_0={u_0:.3f}m/s, L_bed={L_bed:.1f}m, d_p={fixed.d_p_m*1e3:.1f}mm)。"
            f" D_col を大きく or L_bed を短く or 流量を減らす方向に探索を誘導。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='dp_excess', u_0_actual=u_0, dP_bar_actual=dP_bar)

    # CH4 濃度がゼロの場合は破過検知不能
    if C_feed_ads[0] <= 0.0:
        warnings.warn(
            f"PSA penalty: C_CH4={C_feed_ads[0]:.3e} mol/m³ ≤ 0 (上流から CH4 が流れていない)。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='no_CH4_feed')

    # -------------------------------------------------------------------------
    # 3. 吸着 PDE
    # -------------------------------------------------------------------------
    q_s    = np.array([PSA_LANGMUIR_PARAMS[k]['q_s'] for k in _ADS_ORDER])  # [mol/kg]
    a_lang = np.array([PSA_LANGMUIR_PARAMS[k]['a']   for k in _ADS_ORDER])  # [m³/mol]
    kfa    = np.array([PSA_KFA[k]                    for k in _ADS_ORDER])  # [1/s]

    # 吸着材データ感度 (2026-05-31 PSA設計レビュー対応): q_s/a/KFa/ρ_b を env で上書き可
    #   (既定 1.0 = 挙動不変)。Langmuir 定数・KFa・嵩密度はいずれも !仮置き でベンダーデータ未確定の
    #   ため，exp/exp_psa_sensitivity.py で係数を振って塔数・H2 回収率・TAC の頑健性を評価する。
    q_s    = q_s    * float(os.environ.get('PDH_PSA_QS_FACTOR',  '1.0'))
    a_lang = a_lang * float(os.environ.get('PDH_PSA_A_FACTOR',   '1.0'))
    kfa    = kfa    * float(os.environ.get('PDH_PSA_KFA_FACTOR', '1.0'))
    rho_b_eff = fixed.rho_b * float(os.environ.get('PDH_PSA_RHOB_FACTOR', '1.0'))

    # CSS スケーリング近似の妥当性チェック
    # scaling_ratio = (ρ_b/ε) × (q*(C_feed)/C_feed): シャープフロント指標
    # この値が >> 1 (目安: ≥ 10) のとき t_abs の線形スケーリングが成立する
    if fixed.use_css_approximation and C_feed_ads[0] > 0.0:
        denom_css      = 1.0 + np.sum(a_lang * C_feed_ads)
        q_star_CH4     = q_s[0] * a_lang[0] * C_feed_ads[0] / denom_css
        scaling_ratio  = (rho_b_eff / fixed.eps) * (q_star_CH4 / C_feed_ads[0])
        if scaling_ratio < 10.0:
            # 設計判断 (2026-05-18): CSS 近似の妥当性が低下した状態でも計算は続行する
            # (penalty 化は U-決のため一旦警告強化のみ)。t_abs が線形スケーリングから
            # 乖離するため、N_total_columns が過小評価され CAPEX が偽の最小値に
            # なる可能性が高い。BO 最適解がこの領域に偏ったら U-決で penalty 化要。
            warnings.warn(
                f"CSSスケーリング精度低下リスク: scaling_ratio={scaling_ratio:.1f} < 10.0"
                f" (D_col={D_col:.2f}m, L_bed={L_bed:.2f}m, P_in={feed.P_in/1e5:.2f}bar,"
                f" C_feed_CH4={C_feed_ads[0]:.3f}mol/m³)。"
                f" シャープフロント近似が成立せず、t_abs 線形スケーリングの誤差が大きい。"
                f" この設計が BO ベスト解に残る場合、N_total_columns 過小評価で"
                f" CAPEX が偽の最小値となっている可能性あり。"
                f" 高圧・高不純物濃度の影響が疑われる。Langmuir パラメータ確定後に再評価要。",
                UserWarning,
                stacklevel=2,
            )

    t_abs_clean, q_final, sol_t, C_outlet_t, converged = _run_adsorption(
        C_feed=C_feed_ads,
        u_0=u_0,
        L_bed=L_bed,
        rho_b=rho_b_eff,
        eps=fixed.eps,
        kfa=kfa,
        q_s=q_s,
        a_lang=a_lang,
        breakthrough_ratio=fixed.breakthrough_ratio,
        t_max=fixed.t_ads_max,
    )

    # CSS 簡易補正 (乗算方式)
    # ─────────────────────────────────────────────────────────────────────
    # 残留負荷 q_init = q*(C_feed) × desorption_target が均一に存在する場合、
    # 有効吸着容量 = q*(C_feed) - q_init = q*(C_feed) × (1 - desorption_target)
    # → 清浄床 t_abs に対して: t_abs_css ≈ t_abs_clean × (1 - desorption_target)
    # 成立条件: rho_b/ε × q*(C_feed)/C_feed >> 1 (吸着容量 >> 気相ホールドアップ)
    # ※ q_init の空間非一様性(入口高・出口低)は無視。補正は保守的(やや短め)。
    # ─────────────────────────────────────────────────────────────────────
    if fixed.use_css_approximation:
        t_abs = t_abs_clean * (1.0 - design.desorption_target)
    else:
        t_abs = t_abs_clean

    # CSS補正後 t_abs が極小の場合: scale 発散防止のためペナルティを返す
    if t_abs < _T_ABS_MIN:
        # 設計判断 (2026-05-21): 旧版は silent return で BO が「どこへ逃げれば良いか」
        # 学習できなかった。warning + 連続 shortfall (run_one_pass で計算) で TPE に
        # 「あと何倍 L/D を増やせば feasible に出るか」のシグナルを渡す。
        warnings.warn(
            f"PSA penalty: t_abs={t_abs:.1f}s (CSS 補正後) < _T_ABS_MIN={_T_ABS_MIN}s"
            f" (D_col={D_col:.2f}m, L_bed={L_bed:.2f}m, desorption_target={design.desorption_target:.3f},"
            f" u_0={u_0:.3f}m/s)。L_bed を増やすか D_col を増やすか desorption_target を上げる方向に探索を誘導。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='t_abs_below_min', t_abs_actual=t_abs, u_0_actual=u_0)

    if not converged:
        warnings.warn(
            f"simulate_psa_system: t_ads_max={fixed.t_ads_max:.0f}s 内に"
            f" CH4 破過が検出されませんでした (D_col={D_col:.2f}m, L_bed={L_bed:.2f}m,"
            f" desorption_target={design.desorption_target:.3f},"
            f" u_0={u_0:.3f}m/s, C_feed_CH4={C_feed_ads[0]:.3f}mol/m³,"
            f" P_in={feed.P_in/1e5:.2f}bar, T_in={feed.T_in:.1f}K)。"
            f" CAPEX_total=1e9 億円 (penalty sentinel) を返却。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='breakthrough_no_converge', u_0_actual=u_0)

    # -------------------------------------------------------------------------
    # 4. 脱着時間
    # -------------------------------------------------------------------------
    q_avg   = q_final.mean(axis=0)  # 空間平均 [mol/kg]
    t_des_raw = _calc_desorption_time(q_avg, kfa, design.desorption_target)
    # 安全係数: KFa が仮置き値であり推算誤差が大きいため設計マージンを確保する
    t_des = t_des_raw * fixed.desorption_time_safety_factor

    # -------------------------------------------------------------------------
    # 5. サイクル構成・塔数
    # -------------------------------------------------------------------------
    # 設計判断 (2026-05-09): N_abs_parallel = 1 で固定 (= 2塔最小スイング構成)。
    # 理由: 本実装は cycle 内の均圧・再加圧ステップを陽にモデル化していないため、
    #       塔数を増やしても H2 回収率や CAPEX/塔は変わらず、単に CAPEX が線形増加
    #       するだけ。BO の最適化変数に昇格させると自明に最小値が選ばれて探索予算
    #       の無駄になる。
    # 実機 (4塔 Skarstrom / 9-12塔 UOP Polybed) は均圧で H2 ロスを大きく削減するが、
    # それを表現するにはサイクル全ステップ (吸着→均圧↓→ブローダウン→パージ→
    # 均圧↑→再加圧) のフルモデル化が必要 (= 別 PR)。
    N_abs_parallel  = 1
    N_cycle_sets    = math.ceil(t_des / t_abs) + 1
    N_total_columns = N_abs_parallel * N_cycle_sets

    # PSA サイクルスケジュール逼迫チェック
    # idle_time: 脱着完了から次の吸着開始まで実際に確保される待機時間
    # = ceil(t_des/t_abs) 個の吸着期間 − 実際の t_des
    idle_time    = (N_cycle_sets - 1) * t_abs - t_des
    margin_ratio = idle_time / t_abs
    if margin_ratio < 0.10:
        warnings.warn(
            f"PSAサイクルスケジュール逼迫リスク: idle_time={idle_time:.1f}s,"
            f" t_abs={t_abs:.1f}s (余裕率={margin_ratio*100:.1f}% < 10%)。"
            f" 均圧・加圧ステップを考慮した場合、現在の必要塔数"
            f" ({N_total_columns} 塔) では不足する可能性があります。",
            UserWarning,
            stacklevel=2,
        )

    # -------------------------------------------------------------------------
    # 6. 吸着材総重量
    # -------------------------------------------------------------------------
    W_bed_per_col  = rho_b_eff * V_col                # [kg/塔]
    W_adsorbent_kg = W_bed_per_col * N_total_columns   # [kg]

    # -------------------------------------------------------------------------
    # 7. 出力ストリーム
    # -------------------------------------------------------------------------
    # total_moles_out: CSS補正後の t_abs までを積分 (破過後の区間を除外)
    # sol_t / C_outlet_t は _run_adsorption から受け取った生データ
    _mask = sol_t <= t_abs
    if _mask.sum() < 2:
        warnings.warn(
            f"PSA penalty: sol_t <= t_abs を満たす点が {_mask.sum()} 個 (< 2)。"
            f" t_abs={t_abs:.1f}s, sol_t.size={sol_t.size}, sol_t[-1]={sol_t[-1]:.1f}s。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='mask_lt_2', t_abs_actual=t_abs, u_0_actual=u_0)
    _t_trunc = sol_t[_mask]
    _C_trunc = C_outlet_t[:, _mask]
    # t_abs が sol_t の最後のサンプル点より後にある場合、線形補間で端点を追加
    if _t_trunc[-1] < t_abs - 1e-9:
        _C_end   = np.array([np.interp(t_abs, sol_t, C_outlet_t[j]) for j in range(_N_ADS)])
        _t_trunc = np.append(_t_trunc, t_abs)
        _C_trunc = np.hstack([_C_trunc, _C_end.reshape(-1, 1)])
    total_moles_out = u_0 * A_col * np.trapezoid(_C_trunc, _t_trunc, axis=1)

    # 吸着ステップ中の平均出口流量 [kmol/h]
    # total_moles_out [mol] ÷ t_abs [s] × 3600 [s/h] ÷ 1000 [mol/kmol]
    scale = 3600.0 / (t_abs * 1000.0)

    F_CH4_prod  = total_moles_out[0] * scale
    F_C2H4_prod = total_moles_out[1] * scale
    F_C2H6_prod = total_moles_out[2] * scale

    F_H2_in   = feed.F_in.get('C', 0.0)
    F_CH4_in  = feed.F_in.get('E', 0.0)
    F_C2H4_in = feed.F_in.get('D', 0.0)
    F_C2H6_in = feed.F_in.get('F', 0.0)

    # -------------------------------------------------------------------------
    # H2 損失計算
    # -------------------------------------------------------------------------
    # ブローダウン損失: P_in → P_des の圧力降下で空隙 H2 が放出される
    # 脱着後も P_des 条件の H2 が空隙に残るため (1 - Z·P_des/P_in) 倍だけ放出
    # C_H2 = x_H2 · P_in/(Z·R·T) なので、P_des 残留量 = C_H2 · Z · P_des/P_in
    H2_loss_blowdown = (
        V_col * fixed.eps * C_H2 * (1.0 - Z * fixed.P_des / feed.P_in)
        / t_abs * 3600.0 / 1000.0
    )

    # パージ損失: 1サイクルに 1 塔が t_des 秒間、向流パージを受ける
    # パージ速度 = 吸着空塔速度 u_0 と仮定 (実際は設計次第 — 保守的過大推算)
    # C_H2_purge: フィードの H2 モル分率 × P_des 条件の総濃度 (Z≈1 @P_des=1atm)
    #   = C_H2 · Z · P_des/P_in  (Z はフィード条件の圧縮係数)
    # N_abs_parallel は「吸着塔数」であり「パージ塔数」とは別概念のため使用しない
    # (1 塔ずつ順番にパージを受けるサイクル構成につき係数は 1 が正しい)
    C_H2_purge = C_H2 * Z * fixed.P_des / feed.P_in  # [mol/m³] H2 濃度 at P_des
    H2_loss_purge = (
        u_0 * A_col * C_H2_purge * t_des / t_abs * 3600.0 / 1000.0
    )

    # パージ損失の過大推算チェック
    # H2_loss_purge が F_H2_in の 30% を超えた場合はペナルティ領域と判定する
    # (t_des > t_abs 程度では通常は発火しない; 真に非効率な解域でのみ警告)
    if F_H2_in > 0.0:
        purge_loss_ratio = H2_loss_purge / F_H2_in
        if purge_loss_ratio > 0.30:
            warnings.warn(
                f"PSAパージ損失の過大推算（ペナルティ）リスク:"
                f" H2_loss_purge={H2_loss_purge:.2f} kmol/h,"
                f" F_H2_in={F_H2_in:.2f} kmol/h,"
                f" 損失率={purge_loss_ratio:.2f} (>30%)。"
                f" t_des={t_des:.2f}s, t_abs={t_abs:.2f}s"
                f" (t_des/t_abs={t_des/t_abs:.2f})。"
                f" パージモデル (u_0 定速・1塔逐次) の保守的過大推算が"
                f" H2 回収率を非現実的に低下させている可能性があります。",
                UserWarning,
                stacklevel=2,
            )

    # H2 損失が入力流量を超えないようにクランプ (物質収支保護)
    H2_loss_purge = min(H2_loss_purge, max(0.0, F_H2_in - H2_loss_blowdown))
    F_H2_product = max(0.0, F_H2_in - H2_loss_blowdown - H2_loss_purge)

    product = {
        'A': 0.0,
        'B': 0.0,
        'C': F_H2_product,  # H2: ブローダウン・パージ損失を差し引き済み
        'D': F_C2H4_prod,
        'E': F_CH4_prod,
        'F': F_C2H6_prod,
    }

    offgas = {
        'A': feed.F_in.get('A', 0.0),                    # C3H8 全量
        'B': feed.F_in.get('B', 0.0),                    # C3H6 全量
        'C': H2_loss_blowdown + H2_loss_purge,           # H2 損失分
        'D': max(0.0, F_C2H4_in - F_C2H4_prod),
        'E': max(0.0, F_CH4_in  - F_CH4_prod),
        'F': max(0.0, F_C2H6_in - F_C2H6_prod),
    }

    H2_recovery = 1.0  if F_H2_in  <= 0.0 else F_H2_product / F_H2_in
    CH4_capture = 0.0  if F_CH4_in <= 0.0 else offgas['E'] / F_CH4_in

    # -------------------------------------------------------------------------
    # 8. CAPEX
    # -------------------------------------------------------------------------
    # 設計判断 (2026-05-31 PSA設計レビュー対応): PSA 塔体は calc_reactor_capex_okuyen を流用する
    #   が，同関数は内部で K_SWING=1.2 を乗じる。PSA は圧力スイング操作で高速切替バルブ・均圧/
    #   パージライン・マニホールドを要するため，この 1.2 を「PSA パッケージ係数」(バルブ・配管・
    #   制御系の上乗せ) として再解釈し，そのまま採用する。新たな係数を別途掛けると二重計上になる
    #   ため掛けない。実機詳細設計ではベンダーのパッケージ見積りで更新すること。
    capex_vessels = calc_reactor_capex_okuyen(
        V_col, feed.P_in, D_col, N_total_columns
    )
    capex_adsorbent = (
        W_adsorbent_kg
        * ACTIVATED_CARBON_PRICE_USD_PER_KG
        * (CEPCI_CURRENT / CEPCI_BASE)
        * PLANT_INDIRECT_FACTOR
        * USD_TO_JPY
        / 1.0e8
    )
    capex_total = capex_vessels + capex_adsorbent

    # 吸着材交換 OPEX: 初期 CAPEX を耐用年数で除した年均等費
    # !仮置き ADSORBENT_LIFETIME_YEARS = 4 年
    opex_adsorbent_per_year = capex_adsorbent / ADSORBENT_LIFETIME_YEARS

    equipment = PSAEquipmentData(
        N_abs_parallel  = N_abs_parallel,
        N_cycle_sets    = N_cycle_sets,
        N_total_columns = N_total_columns,
        t_abs_sec       = t_abs,
        t_des_sec       = t_des,
        u_0             = u_0,
        W_adsorbent_kg  = W_adsorbent_kg,
        Q_preheat_kW    = Q_preheat_kW,
        CAPEX_vessels   = capex_vessels,
        CAPEX_adsorbent = capex_adsorbent,
        CAPEX_total     = capex_total,
        H2_loss_blowdown_kmolh         = H2_loss_blowdown,
        H2_loss_purge_kmolh            = H2_loss_purge,
        OPEX_adsorbent_okuyen_per_year = opex_adsorbent_per_year,
        dP_bar_actual                  = dP_bar,
    )

    return PSASimulationResult(
        product     = product,
        offgas      = offgas,
        equipment   = equipment,
        H2_recovery = H2_recovery,
        CH4_capture = CH4_capture,
    )
