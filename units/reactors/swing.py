"""
PDH スイング反応器システム シミュレーター

触媒失活を伴うプロパン脱水素断熱PFR（スイング操作）を模擬し、
後段分離工程への時間平均ストリームと装置コスト情報を出力する。

========================================================================
触媒モデル: Cr2O3-Al2O3 (Catofin プロセス相当)
========================================================================
コンテスト要項 §3-3 で提供される a (失活係数) データは「架空触媒」のもの
だが、その挙動 (700°C で 30 min で a≈0.04 と分単位で急速失活、400°C で
a≈0.81) は工業 Cr2O3-Al2O3 触媒 (Catofin プロセス, ABB Lummus 1986〜) の
バッチ再生型 PDH の物理と整合する。
Pt-Sn/Al2O3 (UOP Oleflex) は本来 CCR (連続触媒再生) 方式で日〜週オーダー
の緩慢失活であり、本実装のスイング再生方式 (t_regen=30 min) + a データの
分単位失活とは物理整合しない。よって以下の経済パラメータも Cr2O3 系の
文献値を採用 (CATALYST_*, rho_b):
  - 単価 $23/kg ≒ 3,654 円/kg : Mangalindan et al. (2025) ACS Eng. Au
  - 寿命 2.5 年 (2-3 年中央値) : Ni (2022) TUM 博士論文 p.20 Catofin 解説
  - bulk density 900 kg/m³    : Chauruka (2021) Leeds 博士論文 p.86 γ-Al2O3
環境懸念: Cr2O3 触媒は六価 Cr の生成/排出規制が厳しい (大気汚染防止法・
水質汚濁防止法)。実プラント設計では Pt-Sn/Pt-Ga 系への置換を検討すべき
だが、本シミュレータでは a データとの物理整合性を優先。レポートで言及予定。
詳細: KNOWN_PLACEHOLDERS.md §A.4 (2026-05-19 卒業エントリ)

使用方法:
    from units.reactors.swing import (
        DesignVars, FeedStream, FixedParams, simulate_swing_reactor_system
    )
    result = simulate_swing_reactor_system(design, feed, fixed)

========================================================================
⚠️ モデル限界 (レポート記載必須項目、2026-05-17 文書化)
========================================================================
本実装は **「単段断熱 PFR + 並列分割 (V_cat_max=200 m³/基)」** モデルである。
商用 PDH プラント (UOP Oleflex / Lummus Catofin) は実態として:
  - 反応器 1 基あたり L = 4〜6 m, D = 3〜5 m
  - **多段直列構成 (3〜4 段)** で、段間で再加熱 (PDH は強吸熱、ΔH ≈ +124 kJ/mol)
  - 段間ヒーターの追加 CAPEX + 燃料 OPEX が運転コストの数十億円規模を占める

これに対し本モデルは:
  - z_cat (例 BO ベスト 21.4 m) を「実質的に多段を 1 段にまとめた等価長さ」とみなす
  - 段間ヒーター CAPEX/OPEX を計上しない
  - 反応器内の温度プロファイルは断熱で T_in からどんどん下がる
  - 並列分割 (V_cat_max 制約) は計算するが、直列分割 (= 多段化) は未実装

結果として本モデルは **CAPEX/OPEX を過小評価** している可能性が高い:
  1. 反応器本体 CAPEX: 単段大型 < 多段小型 (面積効率の差)
  2. 段間ヒーター CAPEX: 計上なし
  3. 段間ヒーター燃料費: 計上なし
  4. 各段で再加熱 → 各段平均温度上昇 → 平均選択率↓・副反応↑

レポートでは:
  - 「本研究の TAC / Profit は単段断熱 PFR 仮定下の値であり、
     実機の多段+段間 reheat 構成では追加 CAPEX/OPEX が見込まれる」旨を明記すること
  - 改修候補: option C (多段反応器モデル + Ergun 圧損 + 触媒粒径 dp 変数化)

なお本実装には以下の制約は既に組み込み済 (2026-05-17 改修):
  - **空塔速度 SV = 0.5〜3.0 m/s** (触媒設計慣行値、SV 範囲外で infeasible)
  - **触媒最大量 V_cat_max = 200 m³/基** (大型固定床の物理上限、超過時並列分割)
========================================================================
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.config import THERMO_DATA
from src.thermo import PDHThermo
from src.kinetics import PDHKinetics
from src.catalyst_model import calculate_activity_a
from src.cost_calculator import calc_reactor_capex_okuyen

_thermo = PDHThermo()
_kinetics = PDHKinetics()

# 成分順序（状態ベクトルのインデックスと対応）
# A=C3H8, B=C3H6, C=H2, D=C2H4, E=CH4, F=C2H6
_COMPS = ['A', 'B', 'C', 'D', 'E', 'F']

# 化学量論行列 stoich[i, j] = 成分i に対する反応j の量論係数
#        r1   r2   r3
_STOICH = np.array([
    [-1,  -1,   0],  # A (C3H8)
    [+1,   0,   0],  # B (C3H6)
    [+1,   0,  -1],  # C (H2)
    [ 0,  +1,  -1],  # D (C2H4)
    [ 0,  +1,   0],  # E (CH4)
    [ 0,   0,  +1],  # F (C2H6)
], dtype=float)

_T_REF = 298.15  # [K] エンタルピー計算基準温度


# ---------------------------------------------------------------------------
# グローバル関数（ODE 内から呼び出す、パラメータ引き回し不要）
# ---------------------------------------------------------------------------

def calc_a(t: float, T: float, P: float) -> float:
    """触媒活性度 a ∈ [0, 1]（t [min], T [K], P [Pa]）"""
    return calculate_activity_a(T - 273.15, t)


def calc_Cp(T: float) -> dict:
    """各成分のモル熱容量 [J/(mol·K)]。keys: 'A','B','C','D','E','F'"""
    return {k: _thermo.calc_cp(k, T) for k in _COMPS}


def calc_rate_constants(T: float) -> dict:
    """反応速度定数・平衡定数の辞書を返す。k1/k2/k3/K_B/K_eq すべて温度依存。"""
    return {
        'k1':   _kinetics._k1(T),
        'k2':   _kinetics._k2(T),
        'k3':   _kinetics._k3(T),
        'K_B':  _kinetics._K_B(T),
        'K_eq': _thermo.calc_keq(T),
    }


# ---------------------------------------------------------------------------
# データクラス定義
# ---------------------------------------------------------------------------

@dataclass
class DesignVars:
    """最適化アルゴリズムが操作する設計変数"""
    T_in:  float   # 反応器入口温度 [K]
    z_cat: float   # 触媒層長さ [m]
    t_cyc: float   # 1サイクル反応フェーズ運転時間 [min]
    D:     float   # 反応器直径 [m]


@dataclass
class FeedStream:
    """入口流体条件"""
    F_in:   Dict[str, float]  # 各成分入口モル流量 [kmol/h]
                               # keys: 'A'(C3H8),'B'(C3H6),'C'(H2),'D'(C2H4),'E'(CH4),'F'(C2H6)
    T_feed: float              # 加熱炉入口（予熱前）原料温度 [K]
    P_in:   float              # 反応器入口圧力 [Pa]


@dataclass
class FixedParams:
    """固定定数・制約条件"""
    # !仮置き: 触媒再生時間 [min]
    #   - コンテスト未指定、工業典型 (Catofin/Oleflex) 15〜60 min の中央値 30 を採用
    #   - ユーザ判断で変更可。値が変わると N_swing_sets (= ceil(t_regen/t_cyc)+1)
    #     経由で W_cat 総量に効く
    t_regen:              float = 30.0

    # 反応器サイズ制約 [m³]
    #   コンテスト仕様: "1基あたりの最大触媒量は 200 m³、それ以上必要なら複数基"
    #   ─ V_cat (床体積) に対する制約として swing.py で N_parallel 計算に使用
    V_cat_max_per_vessel: float = 200.0

    # 床/容器体積比 [-]
    #   コンテスト仕様: "反応器容量 = 触媒量 × 2 の Process Vessel"
    #   ⇒ V_cat / V_vessel = 0.5 ⇒ (1-eps) = 0.5 ⇒ eps = 0.5
    #   ※「eps」は粒子間空隙率ではなく「床体積 / 容器体積」比の意味で本コード内で使用。
    #     V_cat (床) = A × z_cat × (1-eps)、V_vessel (容器) = A × z_cat。
    #     床自体の粒子間空隙率 (~0.4 typical) は本モデルでは陽に扱わない。
    eps:                  float = 0.5

    # 触媒充填密度 (bulk density) [kg/m³]
    #   触媒モデル: Cr2O3-Al2O3 (Catofin プロセス相当)。コンテスト要項 §3-3
    #     提供の a データ (分単位の急速失活) と本実装のスイング再生方式
    #     (t_regen=30 min) の物理整合性から決定。
    #   出典: Chauruka, S. R. (2021). "The formulation and characterisation of
    #     extruded alumina catalyst supports" (Doctoral thesis, University of
    #     Leeds; White Rose eTheses Online)。γ-アルミナ触媒担体の物性表で
    #     Packed Bulk Density 800-1000 g/L (= kg/m³) と記載。
    #   参照 PDF: report_for_processdesign/references/Chauruka_2021_Leeds_alumina_thesis.pdf
    #   採用値: レンジ 800-1000 の中央値 900 kg/m³。Catofin の Cr2O3-Al2O3
    #     触媒は担体 (γ-Al2O3) が bulk density を支配するため、γ-アルミナ
    #     担体の値を採用 (Cr 担持 12-25 wt% で bulk density 微増想定だが
    #     一次データ未取得のため担体値で代用)。
    #   用途: W_cat = V_cat (床体積) × rho_b で触媒重量算出。
    #     公知の式と一致 (床体積 × 充填密度 = 触媒重量)。
    rho_b:                float = 900.0

    # 空塔速度 (superficial velocity) 制約 [m/s]
    # 設計判断 (2026-05-17): 触媒設計アドバイスで「固定床 PFR の SV は通常 1〜3 m/s」
    # 推奨。本実装では下限を 0.5 にやや緩めて BO 探索空間を確保 (低 SV は滞留時間↑で
    # over-conversion ぎみ、高 SV は圧損↑で実装注意)。範囲外なら infeasible (penalty 化)。
    # SV = Q_vol_inlet / (π/4 × D²)、Q_vol は反応器入口 T_in, P_in での理想気体換算。
    SV_min_m_per_s:       float = 0.5
    SV_max_m_per_s:       float = 3.0

    def __post_init__(self) -> None:
        _checks = {
            "t_regen":              self.t_regen,
            "V_cat_max_per_vessel": self.V_cat_max_per_vessel,
            "eps":                  self.eps,
            "rho_b":                self.rho_b,
            "SV_min_m_per_s":       self.SV_min_m_per_s,
            "SV_max_m_per_s":       self.SV_max_m_per_s,
        }
        for name, val in _checks.items():
            if val <= 0:
                raise ValueError(
                    f"FixedParams.{name}={val} は正値でなければなりません。"
                )
        if self.SV_min_m_per_s >= self.SV_max_m_per_s:
            raise ValueError(
                f"FixedParams: SV_min ({self.SV_min_m_per_s}) >= SV_max ({self.SV_max_m_per_s})"
            )


@dataclass
class EffluentStream:
    """後段分離工程への出口流体情報"""
    F_out_avg: Dict[str, float]  # 各成分出口モル流量の時間平均 [kmol/h]
    T_out_avg: float              # 出口温度の時間平均 [K]
    Q_preheat: float              # T_feed → T_in 予熱熱量 [GJ/h]
    P_out:     float              # 出口圧力 [Pa]


@dataclass
class EquipmentCost:
    """装置・経済性情報"""
    V_vessel_actual:       float  # 1基プロセス容器容積 [m³]
    N_parallel:            int    # 200m³制約による並列基数
    N_swing_sets:          int    # 再生をカバーする切り替えセット数
    N_reactors_total:      int    # 総反応器基数
    Catalyst_Weight_Total: float  # システム全体触媒総量 [kg]
    Reactor_CAPEX:         float  # 全基分建設コスト合計 [億円]


@dataclass
class PerformanceMetrics:
    """プロセス指標"""
    Conversion:  float  # プロパン単通反応率（時間平均）[%]
    Selectivity: float  # プロピレン選択率（時間平均）[%]


@dataclass
class SimulationResult:
    """シミュレーション結果（全出力をまとめるルートオブジェクト）"""
    effluent:    EffluentStream
    equipment:   EquipmentCost
    performance: PerformanceMetrics


# ---------------------------------------------------------------------------
# ペナルティ結果（計算不能条件への早期リターン用）
# ---------------------------------------------------------------------------

_PENALTY_CAPEX: float = 1e9  # [億円] 最適化への無効シグナル


def _penalty_result() -> SimulationResult:
    """計算不能条件のときに返すペナルティ SimulationResult。"""
    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg={c: 0.0 for c in _COMPS},
            T_out_avg=0.0,
            Q_preheat=0.0,
            P_out=0.0,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=0.0,
            N_parallel=0,
            N_swing_sets=0,
            N_reactors_total=0,
            Catalyst_Weight_Total=0.0,
            Reactor_CAPEX=_PENALTY_CAPEX,
        ),
        performance=PerformanceMetrics(
            Conversion=0.0,
            Selectivity=0.0,
        ),
    )


# ---------------------------------------------------------------------------
# 内部ヘルパー
# ---------------------------------------------------------------------------

def _reaction_enthalpies(T: float) -> np.ndarray:
    """温度T [K] での各反応エンタルピー ΔH_rxn [J/mol]。shape (3,)"""
    H = {k: THERMO_DATA[k].dHf_298 + _thermo.calc_enthalpy_change(k, _T_REF, T)
         for k in _COMPS}

    dH1 = H['B'] + H['C'] - H['A']  # C3H8 → C3H6 + H2
    dH2 = H['D'] + H['E'] - H['A']  # C3H8 → C2H4 + CH4
    dH3 = H['F'] - H['D'] - H['C']  # C2H4 + H2 → C2H6
    return np.array([dH1, dH2, dH3])


def _ode_axial(z: float, y: np.ndarray,
               a: float, A_cross: float, eps: float, P_in: float
               ) -> np.ndarray:
    """軸方向(z)の常微分方程式。

    State y = [F_A, F_B, F_C, F_D, F_E, F_F, T]  (A=C3H8, B=C3H6, C=H2, D=C2H4, E=CH4, F=C2H6)
    単位: F [mol/s], T [K]

    a : 触媒活性度 [-]。コーキングは入口温度で一律決定されるため、
        呼び出し元（_simulate_one_time）で事前計算した値を受け取る。
    """
    # 負のモル流量をクリップ（数値積分のアンダーシュート対策）
    # 設計判断 (2026-05-18): clip 量が 1e-3 mol/s を超える場合は警告
    # (通常 ODE のアンダーシュートは 1e-9 程度。-1e-3 mol/s = -3.6 kmol/h)
    F_raw = y[:6]
    F = np.maximum(F_raw, 0.0)
    _clip_amount = float(np.sum(np.maximum(-F_raw, 0.0)))
    if _clip_amount > 1.0e-3:
        warnings.warn(
            f"swing._ode_axial: 負モル流量を clip (合計 {_clip_amount:.3e} mol/s)。"
            f" ODE 安定性低下、結果の信頼度に注意。",
            RuntimeWarning, stacklevel=2,
        )
    # 物理的温度範囲に制限（ODE 発散防止）
    # 設計判断 (2026-05-18): [300K, 1500K] は物理パラメータではなく数値ガード。
    #   下限 300K (27°C): 通常運転温度 (700-950K) より遥かに低く、ODE が中間 step で
    #                      負温度域に外れたとき計算続行できないため。
    #   上限 1500K: 通常運転温度の 1.5× 程度。これを超える点を ODE が探索したら
    #                Radau の step が暴走しているサイン。
    # !仮置き 文献値ではなく数値ガード (citation 不要)。
    T_local = float(np.clip(y[6], 300.0, 1500.0))

    F_total = float(np.sum(F))
    if F_total <= 0:
        return np.zeros(7)

    # 局所分圧 [Pa]（圧力損失なしのため P_total = P_in）
    P_local = P_in
    P = {comp: max(float(F[i]) / F_total * P_local, 0.0)
         for i, comp in enumerate(_COMPS)}

    # 反応速度定数（ゼロ除算・駆動力発散防止のため最小値でクリップ）
    rc = calc_rate_constants(T_local)
    k1, k2, k3 = rc['k1'], rc['k2'], rc['k3']
    # 設計判断 (2026-05-18): K_B と K_eq の単位は [Pa] (P_B との比が無次元になる必要)。
    #   K_B (吸着平衡定数): 低温で 0 に近づくと adsorption = 1 + P_B/K_B が発散
    #                       → r1 → 0 (吸着サイト飽和)。floor 1 Pa は通常運転圧 (~0.5 bar
    #                       = 5e4 Pa) より 5 桁小さく、物理計算には実質影響しない数値ガード。
    #   K_eq (反応平衡定数): 低温で 0 に近づくと P_B·P_C/K_eq が発散 → driving_r1 が
    #                       極端に負 → r1 → 大きな負値 (逆反応暴走)。floor 1 Pa は
    #                       同様の数値ガード。
    # citation: 物性値ではなく ODE 数値安定化のための実装判断。
    K_B  = max(rc['K_B'],  1.0)   # [Pa] 低温での吸着項ゼロ除算防止
    K_eq = max(rc['K_eq'], 1.0)   # [Pa] 低温での駆動力発散防止

    # 反応速度 [mol/m³_cat/s]
    driving_r1 = P['A'] - P['B'] * P['C'] / K_eq
    r1 = a * k1 * driving_r1 / (1.0 + P['B'] / K_B)
    r2 = k2 * P['A']
    r3 = k3 * P['D'] * P['C']

    rates = np.array([r1, r2, r3])  # [mol/m³_cat/s]

    # 物質収支: dF_i/dz = (1-eps) * A * Σ stoich[i,j] * r[j]  [mol/(s·m)]
    dFdz = (1.0 - eps) * A_cross * (_STOICH @ rates)

    # エネルギー収支（断熱）: dT/dz = -((1-eps) * A * Σ ΔH_j * r_j) / Σ F_i * Cp_i
    cp_dict = calc_Cp(T_local)
    sum_FCp = sum(max(float(F[i]), 0.0) * cp_dict[comp] for i, comp in enumerate(_COMPS))
    if sum_FCp <= 0:
        dTdz = 0.0
    else:
        dH = _reaction_enthalpies(T_local)
        Q_rxn = float((1.0 - eps) * A_cross * np.dot(dH, rates))  # [J/(m·s)]
        dTdz = -Q_rxn / sum_FCp                                     # [K/m]

    return np.concatenate([dFdz, [dTdz]])


def _simulate_one_time(design: DesignVars, feed: FeedStream,
                       fixed: FixedParams, t_min: float) -> tuple:
    """時刻 t_min [min] における空間積分を実施し (F_out [mol/s], T_out [K]) を返す。"""
    A_cross = math.pi / 4.0 * design.D ** 2  # [m²]

    # 触媒活性：コーキングは入口温度で一律決定。ODE 内では使用しない。
    a = calc_a(t_min, design.T_in, feed.P_in)

    F0 = np.array([feed.F_in.get(c, 0.0) * 1000.0 / 3600.0 for c in _COMPS])
    y0 = np.concatenate([F0, [design.T_in]])

    # 設計判断 (2026-05-08, profile 結果反映):
    # 元値 rtol=1e-5, atol=1e-8 は scipy デフォルト (1e-3, 1e-6) の 100倍厳しく、
    # profile で全体時間の 72% を占めていた。製品純度スペック 99.5wt% に対して
    # ODE 精度 0.01% (rtol=1e-4) で十分余裕があるため緩和して 2-3倍速を狙う。
    # 緩めすぎ (1e-3) は転化率推定がノイズっぽくなるため中庸を選択。
    try:
        sol = solve_ivp(
            fun=lambda z, y: _ode_axial(z, y, a, A_cross, fixed.eps, feed.P_in),
            t_span=(0.0, design.z_cat),
            y0=y0,
            method='Radau',
            rtol=1e-4,
            atol=1e-7,
        )
    except Exception as e:
        warnings.warn(
            f"swing.solve_ivp 例外: {type(e).__name__}: {e} "
            f"(T_in={design.T_in:.1f}K, z_cat={design.z_cat:.2f}m, "
            f"D={design.D:.2f}m, t={t_min:.1f}min, P_in={feed.P_in/1e5:.2f}bar)",
            RuntimeWarning, stacklevel=2,
        )
        return None, None

    if not sol.success:
        warnings.warn(
            f"swing.solve_ivp 収束失敗 (status={sol.status}, message={sol.message!r}) "
            f"(T_in={design.T_in:.1f}K, z_cat={design.z_cat:.2f}m, "
            f"D={design.D:.2f}m, t={t_min:.1f}min, P_in={feed.P_in/1e5:.2f}bar)",
            RuntimeWarning, stacklevel=2,
        )
        return None, None

    return sol.y[:6, -1], float(sol.y[6, -1])


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_swing_reactor_system(
    design: DesignVars,
    feed: FeedStream,
    fixed: FixedParams,
    n_time_samples: int = 20,
) -> SimulationResult:
    """スイング反応器システムをシミュレーションする。

    Parameters
    ----------
    design         : DesignVars   設計変数
    feed           : FeedStream   入口流体条件
    fixed          : FixedParams  固定パラメータ
    n_time_samples : int          時間軸サンプリング点数（デフォルト20）
    """
    # ---- 入力バリデーション ----
    if design.t_cyc <= 0 or design.z_cat <= 0 or design.D <= 0:
        return _penalty_result()
    if design.T_in <= 0 or feed.T_feed <= 0 or feed.P_in <= 0:
        return _penalty_result()
    if any(v < 0 for v in feed.F_in.values()):
        return _penalty_result()
    if sum(feed.F_in.values()) <= 0:
        return _penalty_result()

    # ---- 時間方向サンプリングと空間積分 ----
    if n_time_samples < 2:
        warnings.warn(
            f"n_time_samples={n_time_samples} < 2: 台形積分を行わず t=0 の1点を時間平均として採用します。",
            UserWarning,
            stacklevel=2,
        )
        F_out, T_out = _simulate_one_time(design, feed, fixed, 0.0)
        if F_out is None:
            return _penalty_result()
        F_out_avg_kmolh = {comp: float(F_out[i]) * 3600.0 / 1000.0
                           for i, comp in enumerate(_COMPS)}
        T_out_avg = T_out
    else:
        t_samples = np.linspace(0.0, design.t_cyc, n_time_samples)
        F_out_list, T_out_list = [], []

        for t in t_samples:
            F_out, T_out = _simulate_one_time(design, feed, fixed, float(t))
            if F_out is None:
                return _penalty_result()
            F_out_list.append(F_out)
            T_out_list.append(T_out)

        F_out_arr = np.array(F_out_list)  # (n_time_samples, 6) [mol/s]
        T_out_arr = np.array(T_out_list)  # (n_time_samples,)   [K]

        # 時間平均（台形則）
        _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
        F_out_avg_mol_s = _trapz(F_out_arr, t_samples, axis=0) / design.t_cyc
        T_out_avg = float(_trapz(T_out_arr, t_samples) / design.t_cyc)

        # mol/s → kmol/h
        F_out_avg_kmolh = {comp: float(F_out_avg_mol_s[i]) * 3600.0 / 1000.0
                           for i, comp in enumerate(_COMPS)}

    # ---- 予熱熱量 Q_preheat [GJ/h] ----
    q_w = 0.0
    for comp in _COMPS:
        F_mol_s = feed.F_in.get(comp, 0.0) * 1000.0 / 3600.0
        q_w += F_mol_s * _thermo.calc_enthalpy_change(comp, feed.T_feed, design.T_in)
    Q_preheat_GJh = q_w * 3600.0 / 1e9

    # ---- 装置計算 ----
    A_cross = math.pi / 4.0 * design.D ** 2

    # 空塔速度 (superficial velocity) チェック
    # 設計判断 (2026-05-17): 触媒設計アドバイス (固定床 PFR、1〜3 m/s) に従い、
    # FixedParams.SV_{min,max}_m_per_s 範囲外なら infeasible 化。
    # SV は反応器入口 (T_in, P_in) での理想気体換算体積流量 / 断面積。
    # 並列分割される場合 (V_cat > V_cat_max) は 1 基あたりの断面積で評価する必要が
    # あるが、まず全体断面で評価し、N_parallel ≥ 1 とした上で SV をスケールする。
    n_inlet_mol_s = sum(feed.F_in.values()) * 1000.0 / 3600.0   # kmol/h → mol/s
    if feed.P_in > 0:
        Q_vol_m3_s = n_inlet_mol_s * 8.314 * design.T_in / feed.P_in
    else:
        Q_vol_m3_s = 0.0
    # N_parallel をここで仮計算 (後段と同じロジック)
    V_cat_total = A_cross * design.z_cat * (1.0 - fixed.eps)
    N_parallel = max(math.ceil(V_cat_total / fixed.V_cat_max_per_vessel), 1)
    # 並列分割される場合、1 基あたりの断面積で SV 評価
    SV_m_per_s = Q_vol_m3_s / (A_cross * N_parallel) if A_cross > 0 else 0.0
    if not (fixed.SV_min_m_per_s <= SV_m_per_s <= fixed.SV_max_m_per_s):
        # 範囲外: infeasible として返す
        warnings.warn(
            f"swing reactor: SV={SV_m_per_s:.2f} m/s が範囲 "
            f"[{fixed.SV_min_m_per_s}, {fixed.SV_max_m_per_s}] m/s 外 "
            f"(D={design.D:.2f}m, N_parallel={N_parallel}, T_in={design.T_in:.1f}K, "
            f"P_in={feed.P_in/1e5:.2f}bar) — infeasible 化",
            UserWarning, stacklevel=2,
        )
        return _penalty_result()

    N_swing_sets = math.ceil(fixed.t_regen / design.t_cyc) + 1
    N_reactors_total = N_parallel * N_swing_sets

    V_vessel_actual = (V_cat_total / N_parallel) / (1.0 - fixed.eps)  # [m³]
    # W_cat = (床体積) × N_swing_sets × ρ_bulk
    #   V_cat_total : 床体積 = π/4×D²×z_cat×(1-eps) [m³]、eps は「床/容器比」
    #   N_swing_sets: 再生中の塔含めた全セット分の触媒が必要
    #   rho_b       : 床全体の充填密度 [kg/m³] (粒子間空隙込み)
    catalyst_weight_total = V_cat_total * N_swing_sets * fixed.rho_b  # [kg]

    if V_vessel_actual <= 0:
        return _penalty_result()

    # CAPEX: Bare Module Cost法による推算（縦型プロセス容器）
    try:
        reactor_capex = calc_reactor_capex_okuyen(
            V_vessel_m3=V_vessel_actual,
            P_abs_pa=feed.P_in,
            D_m=design.D,
            N_reactors_total=N_reactors_total,
        )
    except Exception:
        reactor_capex = _PENALTY_CAPEX

    # ---- パフォーマンス指標 ----
    F_A_in  = feed.F_in.get('A', 0.0)
    F_A_out = F_out_avg_kmolh['A']
    F_B_in  = feed.F_in.get('B', 0.0)
    F_B_out = F_out_avg_kmolh['B']

    conversion  = (F_A_in - F_A_out) / F_A_in * 100.0 if F_A_in > 0 else 0.0
    delta_A     = F_A_in - F_A_out
    selectivity = (F_B_out - F_B_in) / delta_A * 100.0 if delta_A > 0 else 0.0
    conversion  = float(np.clip(conversion,  0.0, 100.0))
    selectivity = float(np.clip(selectivity, 0.0, 100.0))

    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg=F_out_avg_kmolh,
            T_out_avg=T_out_avg,
            Q_preheat=Q_preheat_GJh,
            P_out=feed.P_in,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=V_vessel_actual,
            N_parallel=N_parallel,
            N_swing_sets=N_swing_sets,
            N_reactors_total=N_reactors_total,
            Catalyst_Weight_Total=catalyst_weight_total,
            Reactor_CAPEX=reactor_capex,
        ),
        performance=PerformanceMetrics(
            Conversion=conversion,
            Selectivity=selectivity,
        ),
    )
