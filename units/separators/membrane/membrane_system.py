"""
膜分離システム シミュレーター

5ユニット構成:
  1. 気化器 (Vaporizer)         前段蒸留塔底液 → 全量ガス化
  2. フィード圧縮機 (FeedComp)  P_in → P_H
  3. 膜分離モジュール (Membrane) クロスフロー ODE
  4. 製品圧縮機 (ProductComp)   P_L → P_dist
  5. 製品冷却器 (Condenser)     → 飽和液（後段蒸留塔フィード）

接続:
  [前段蒸留塔底液] --MemFeedStream--> simulate_membrane_system
  --> MemSimulationResult.product   --> [後段蒸留塔フィード（飽和液）]
  --> MemSimulationResult.retentate --> [リサイクルまたは排出]

使用方法:
    from units.separators.membrane_system import (
        MemDesignVars, MemFeedStream, MemFixedParams,
        simulate_membrane_system,
    )
    result = simulate_membrane_system(design, feed, fixed)

--------------------------------------------------------------------
【仮定一覧】

■ 文献・仕様書に根拠がある仮定

  [膜パラメータ] Q_A = 40 GPU, 選択性 α = 90
      根拠: Hua et al. (2024) "Unexpectedly High Propylene/Propane
            Separation Performance..." 実測値（室温・大気圧条件）。

  [膜等温仮定] 透過ガス入口温度 = フィード圧縮機出口温度（膜内温度変化なし）
      根拠: 膜パラメータが等温・大気圧条件で測定されており、モデル条件を
            測定条件に合わせることで整合性を確保する（同論文 Gas Permeation
            Measurements セクション）。

  [ユーティリティ温度] T_hot=160°C, T_cold_in=30°C, T_cold_out=40°C
      根拠: 第17回プロセスデザインコンテスト課題 Ver.2.0 のサイト仕様。
            熱媒は入手可能スチームのうち最安の LP Steam (160°C) を選択。

■ 設計判断（根拠あり）

  [製品冷却器の冷媒] 冷却水のみ使用（Case A）
      判断: P_dist を低くしすぎると製品の泡点が T_cold_out (40°C) を下回り
            冷却水では凝縮できなくなる（温度クロスが発生）。PR EOS 試算では
            C3H6 97% 組成の場合 P_dist ≳ 17 bar で冷却水が使えるが、
            組成によって変わるため simulate 内で動的に確認しペナルティ返却。
            冷媒モデルを入れると目的関数が不連続になり最適化に不向きなため
            Case A（ペナルティ制約）を採用。
      TODO: 後段蒸留塔を含む系統最適化に移行する段階で冷媒使用を再検討する。

  [クロスフロー完全混合モデル]
      判断: 産業用スパイラル型モジュールの標準近似。完全混合（タンク型）より
            現実に近く、向流モデルより実装が安定。

  [PR EOS 二成分相互作用係数 kij = 0]
      判断: C3H6/C3H8 は分子構造が近く kij は文献で 0.01 程度と小さいため
            無視しても相平衡計算の精度に大きな影響はない。

■ 暫定値（設計ヒューリスティクス、要メーカー確認）

  U_vap      = 1.0 kW/(m²·K)  気化器（軽質炭化水素蒸発/水蒸気加熱）
                               化工便覧 改訂六版 表6・18  範囲 0.45〜1.14 の中央〜上限値
  U_cond     = 1.0 kW/(m²·K)  冷却器（軽質炭化水素凝縮/冷却水）
                               化工便覧 改訂六版 表6・18  範囲 0.45〜1.14 の中央〜上限値
  eta_comp   = 0.75            圧縮機断熱効率
                               化工便覧 改訂六版 p.333  ポリトロープ効率 0.7〜0.8 の中央値
                               （断熱効率とポリトロープ効率は厳密に異なるが初期設計では同値で近似）
  T_vap_superheat = 5 K        気化器出口の過熱度（露点 + 5K）

■ !仮置き値 — 根拠文献未確定（要確認・要更新）

  A_per_module = 500 m²        モジュール 1 本あたり有効膜面積
                               確認方法: Evonik SEPURAN・UBE 等のデータシートから
                               中空糸寸法（外径・長さ・本数）を取得して算出する。
                               実装箇所: MemFixedParams.A_per_module
                               → 呼び出し時に UserWarning を発行。

  膜モジュール単価 = 50 USD/m² 根拠未確定。高分子膜の概算値を暫定使用。
                               確認方法: Hua et al. (2024) または ZIF-8 膜 TEA 論文。
                               参考: 高分子膜の文献値 $50〜100/m² (Baker & Lokhandwala 2008)
                               実装箇所: cost_parameters.MEM_UNIT_PRICE_USD_PER_M2
                               → calc_mem_capex_okuyen 呼び出し時に UserWarning を発行。
--------------------------------------------------------------------
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.eos import (
    z_factor, residual_enthalpy,
    bubble_point_T, dew_point_T,
    compress_isentropic,
)
from src.thermo import PDHThermo
from src.cost_calculator import calc_he_capex_okuyen, calc_comp_capex_okuyen, calc_mem_capex_okuyen

# 成分順序: index 0 = C3H6 (EOS キー 'B'), index 1 = C3H8 (EOS キー 'A')
_KEYS    = ['B', 'A']
_T_REF   = 298.15        # [K]  エンタルピー基準温度
_GPU_SI  = 3.35e-10      # [mol/(m²·s·Pa)]  1 GPU の SI 換算値
_ATM_BAR = 1.01325       # [bar] 大気圧（ゲージ圧変換用）

_thermo = PDHThermo()


# ---------------------------------------------------------------------------
# データクラス: 入力
# ---------------------------------------------------------------------------

@dataclass
class MemDesignVars:
    """最適化アルゴリズムが操作する設計変数"""
    P_H:    float  # 膜供給側（高圧）圧力 [Pa]
    P_L:    float  # 膜透過側（低圧）圧力 [Pa]
    A_mem:  float  # 総膜面積 [m²]
    # 後段蒸留塔操作圧力 [Pa]（製品圧縮機の出口圧力）
    # 冷媒不使用（Case A）のため T_bp(P_dist, y_perm) > T_cold_out(40°C) が
    # 成立する圧力に制限される。PR EOS による試算:
    #   C3H6 97% 組成の場合 ≳ 17 bar 程度が目安
    #   （組成により変わるため simulate_membrane_system 内で動的に確認）
    # 最適化ソルバー側で下限境界として設定すること。
    P_dist: float

    def __post_init__(self) -> None:
        if self.P_dist <= 0:
            raise ValueError("P_dist は正値でなければなりません。")
        if self.P_dist <= self.P_L:
            raise ValueError(
                f"P_dist ({self.P_dist/1e5:.2f} bar) は P_L ({self.P_L/1e5:.2f} bar) "
                "より大きくなければなりません（製品圧縮機が圧縮できません）。"
            )


@dataclass
class MemFeedStream:
    """前段蒸留塔底液 → 気化器への入力ストリーム"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_in:   float  # 液温度 [K]（飽和液または過冷却液）
    P_in:   float  # 圧力 [Pa]


@dataclass
class MemFixedParams:
    """プロセス固定パラメータ"""
    # 膜性能 — Hua et al. (2024) 実測値
    Q_A_GPU:         float = 40.0    # C3H6 透過度 [GPU]
    alpha:           float = 90.0    # C3H6/C3H8 選択性 [-]
    # !仮置き — メーカーカタログで要確認（Evonik SEPURAN 等のデータシートから算出）
    A_per_module:    float = 500.0   # モジュール 1 本あたり有効膜面積 [m²]（仮置き）
    # 気化器
    T_vap_superheat: float = 5.0     # 露点超過の過熱度 [K]（設計ヒューリスティクス）
    U_vap:           float = 1.0     # 気化器総括伝熱係数 [kW/(m²·K)]
                                     # 軽質炭化水素の蒸発 - 水蒸気加熱
                                     # 根拠: 化工便覧 改訂六版 表6・18  範囲 0.45〜1.14 の中央〜上限値を採用
    T_hot:           float = 433.15  # LP Steam 供給温度 [K] (160°C) — コンテスト仕様
    # 製品冷却器 — 冷媒不使用（Case A）。冷却水のみ使用
    U_cond:          float = 1.0     # 冷却器総括伝熱係数 [kW/(m²·K)]
                                     # 軽質炭化水素の凝縮 - 冷却水
                                     # 根拠: 化工便覧 改訂六版 表6・18  範囲 0.45〜1.14 の中央〜上限値を採用
    T_cold_in:       float = 303.15  # 冷却水入口温度 [K] (30°C) — コンテスト仕様
    T_cold_out:      float = 313.15  # 冷却水出口温度 [K] (40°C) — コンテスト仕様
    # 圧縮機
    eta_comp:        float = 0.75    # 圧縮機断熱効率 [-]
                                     # 根拠: 化工便覧 改訂六版 p.333「5.6.3 気体圧縮機」
                                     #   ポリトロープ効率は普通 0.7〜0.8 程度 → 中央値 0.75 を採用
                                     # 補足: ポリトロープ効率と断熱効率は厳密には異なる（圧縮比が高いほど
                                     #   断熱効率は下がる）が、初期設計では断熱効率 0.75 を用いるのが一般的。
                                     #   同便覧の例題(p.334)でも 0.8 が使用されている。
    # フィード状態フラグ
    #   False : 液フィード（Vaporizer で気化）— デフォルト
    #   True  : ガスフィード（Vaporizer スキップ、フィード圧縮機が直接受け取る）
    # 用途: Hua et al. (2024) の検証範囲 P_H ≤ 9.5 bar を満たすには、Dist2 を低圧
    #       運転して P_H > feed.P_in を確保する必要がある。低圧では C3H8/C3H6 の
    #       泡点が冷却水温度を下回るため液フィードは不可能。ガスフィードに切替える。
    vapor_feed:      bool  = False

    def __post_init__(self) -> None:
        if not 0 < self.eta_comp <= 1.0:
            raise ValueError("eta_comp は (0, 1] でなければなりません。")
        if self.T_hot <= 273.15:
            raise ValueError("T_hot は 273K 超でなければなりません。")
        # 仮置き値の警告はモジュール初回のみ発火 (BO 数千 trial で log noise になるため)。
        _warn_once_A_per_module(self.A_per_module)


# モジュールレベルの "初回のみ警告" フラグ。同じ値で複数回 MemFixedParams が
# インスタンス化されても警告は 1 度だけに抑える (BO で数千回呼ばれるため)。
_A_PER_MODULE_WARNED_VALUES: set = set()


def _warn_once_A_per_module(value: float) -> None:
    if value in _A_PER_MODULE_WARNED_VALUES:
        return
    _A_PER_MODULE_WARNED_VALUES.add(value)
    warnings.warn(
        f"MemFixedParams: A_per_module = {value} m² は仮置き値です。"
        " Evonik SEPURAN 等のデータシートから中空糸寸法を取得して算出後に更新してください。"
        " (この警告は同一値に対し初回のみ発火します)",
        UserWarning,
        stacklevel=3,
    )


# ---------------------------------------------------------------------------
# データクラス: 出力
# ---------------------------------------------------------------------------

@dataclass
class MemRetentateStream:
    """膜非透過ストリーム（C3H8 富化; リサイクルまたは排出）"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_out:  float  # 温度 [K]（フィード圧縮機出口温度 = 膜等温仮定）
    P_out:  float  # 圧力 [Pa] = P_H


@dataclass
class MemProductStream:
    """製品ストリーム → 後段蒸留塔フィード（飽和液）"""
    F_C3H6: float  # C3H6 モル流量 [kmol/h]
    F_C3H8: float  # C3H8 モル流量 [kmol/h]
    T_out:  float  # 泡点温度 [K] @ P_dist（飽和液）
    P_out:  float  # 圧力 [Pa] = P_dist


@dataclass
class MemEquipmentData:
    """機器サイズ・コスト推算テーブル（授業資料 R08-3 形式）

    CAPEX_mem・CAPEX_total は仮置き値を使用中（UserWarning が発行される）。
    cost_parameters.MEM_UNIT_PRICE_USD_PER_M2 を確定させると自動的に更新される。

    各フィールドの意味:
      A_vap / A_cond   : 熱交換器 伝熱面積 [m²]    → 機器タイプ: 熱交換器
      W_feed / W_prod  : 圧縮機 所要動力 [kW]       → 機器タイプ: 圧縮機
      A_mem / n_modules: 膜モジュール面積・本数      → 機器タイプ: 特殊機器
      Pg_*             : ゲージ圧 [barg]             → 圧力補正係数 Fp 計算用
      Q_vap / Q_cond   : OPEX 用熱量 [kW]
    """
    # 気化器
    A_vap:        float  # 伝熱面積 [m²]
    Pg_vap:       float  # ゲージ圧 [barg]
    Q_vap_kW:     float  # 加熱量 [kW]
    # フィード圧縮機
    W_feed_kW:    float  # 所要動力 [kW]
    Pg_feed:      float  # ゲージ圧 [barg]
    # 膜モジュール
    A_mem:        float  # 総膜面積 [m²]
    n_modules:    int    # モジュール本数
    Pg_mem:       float  # ゲージ圧 [barg]
    # 製品圧縮機
    W_prod_kW:    float  # 所要動力 [kW]
    Pg_prod:      float  # ゲージ圧 [barg]
    # 製品冷却器
    A_cond:       float  # 伝熱面積 [m²]
    Pg_cond:      float  # ゲージ圧 [barg]
    Q_cond_kW:    float  # 冷却量 [kW]
    # CAPEX [億円] — Turton Bare Module Cost 法（R08-3.pdf）
    CAPEX_vap:        float = float('nan')  # 気化器
    CAPEX_comp_feed:  float = float('nan')  # フィード圧縮機
    CAPEX_comp_prod:  float = float('nan')  # 製品圧縮機
    CAPEX_cond:       float = float('nan')  # 製品冷却器
    # !仮置き — cost_parameters.MEM_UNIT_PRICE_USD_PER_M2 が確定次第、自動的に更新される
    CAPEX_mem:        float = float('nan')  # 膜モジュール（単価仮置き中）
    CAPEX_total:      float = float('nan')  # 合計 [億円]（CAPEX_mem が仮置きのため暫定値）
    # 設計判断 (2026-05-08): ヒートインテグレーション用ストリーム温度。
    # 気化器: vapor_feed=True のとき in/out は feed.T_in 同値、Q_vap_kW=0。
    T_vap_in_K:   float = float('nan')  # 気化器入口温度 [K] (= feed.T_in)
    T_vap_out_K:  float = float('nan')  # 気化器出口温度 [K] (= 露点+過熱度、または feed.T_in)
    T_cond_in_K:  float = float('nan')  # 製品冷却器入口温度 [K] (= 製品圧縮機出口)
    T_cond_out_K: float = float('nan')  # 製品冷却器出口温度 [K] (= 泡点 @ P_dist)
    # 設計判断 (2026-05-08): 製品冷却器を顕熱(ガスT_in→T_bp)+潜熱(T_bpで凝縮)に分離。
    # HI 抽出時に温度範囲を持つ顕熱ストリームと潜熱ストリームを別扱いするため。
    # Q_cond_kW = Q_cond_sensible_kW + Q_cond_latent_kW (互換性のため両方保持)。
    Q_cond_sensible_kW: float = 0.0   # ガス顕熱 [kW] (T_cond_in → T_cond_out)
    Q_cond_latent_kW:   float = 0.0   # 凝縮潜熱 [kW] (T_cond_out で凝縮)
    # ---- penalty 診断 (2026-05-22 追加、BO の TPE constraints_func 用) ----
    # 設計判断: PSA / Reactor と同じパターン (psa_system.py:279-287 参照)。
    # silent _penalty_result() 経路が BO に「方向のシグナル」を渡せず 80% の trial が
    # 無方向で死ぬ問題への対処。penalty 発火時に「どの条件で死んだか」「actual 値」を
    # 保持し、run_one_pass._compute_mem_shortfall が log10(actual/limit) 等の連続
    # shortfall を計算して TPE 構成制約に流す。通常完走時は penalty_reason='' のまま。
    # 理由ラベル:
    #   ''                     正常完走 (CAPEX_total < threshold)
    #   'invalid_input'        feed.P_in/T_in/P_dist が非正値 (上流バグ系、本来到達不可)
    #   'ph_le_pl'             P_H ≤ P_L (search bounds で防げる)
    #   'invalid_design'       A_mem/P_H/P_L が非正値
    #   'pdist_le_pl'          P_dist ≤ P_L (search bounds で防げる)
    #   'ph_le_pfeed'          P_H ≤ feed.P_in (Mem feed comp が減圧方向) — BO actionable
    #   'neg_feed'/'zero_feed' feed 流量が負/ゼロ (上流 Dist2 が異常組成)
    #   'dew_nan'              feed 露点計算失敗 (組成異常)
    #   'liquid_vaporized'     vapor_feed=False かつ T_in ≥ T_dew (本 run では vapor_feed=True なので発火しない)
    #   'vapor_condensed'      vapor_feed=True かつ T_in < T_dew — BO actionable (Dist2 P 高すぎ)
    #   'vap_exception'/'vap_nan'  気化器計算失敗 (vapor_feed=False 経路、本 run では発火しない)
    #   'feed_comp_exception'  フィード圧縮機計算失敗 (compress_isentropic が ValueError)
    #   'ode_failure'          膜 ODE 積分失敗
    #   'prod_comp_exception'  製品圧縮機計算失敗
    #   'cond_exception'/'cond_nan'  製品冷却器 LMTD クロス等
    #   'bp_le_cold_out'       透過ガス泡点 ≤ 冷却水出口温度 — BO actionable (P_dist3 低すぎ)
    penalty_reason:    str   = ''
    # actionable 経路の実値 (連続 shortfall 計算用、0.0 = 未計測)
    P_H_actual_Pa:     float = 0.0    # ph_le_pfeed / ph_le_pl で実 P_H [Pa]
    P_feed_actual_Pa:  float = 0.0    # ph_le_pfeed で実 feed.P_in [Pa]
    T_dew_actual_K:    float = 0.0    # vapor_condensed で feed 露点 [K]
    T_feed_actual_K:   float = 0.0    # vapor_condensed で feed.T_in [K]
    T_bp_perm_actual_K:float = 0.0    # bp_le_cold_out で透過ガス泡点 [K]
    T_cold_out_actual_K:float = 0.0   # bp_le_cold_out で冷却水出口温度 [K]


@dataclass
class MemSimulationResult:
    """膜分離システム シミュレーション全出力"""
    retentate:   MemRetentateStream
    product:     MemProductStream
    equipment:   MemEquipmentData
    stage_cut:   float  # θ = F_perm / F_feed [-]
    perm_purity: float  # 透過ガス C3H6 モル分率 [-]
    ret_purity:  float  # 非透過ガス C3H6 モル分率 [-]


# ---------------------------------------------------------------------------
# ペナルティ結果
# ---------------------------------------------------------------------------

_PENALTY = 1e9


def _penalty_result(
    reason:              str   = '',
    P_H_actual:          float = 0.0,
    P_feed_actual:       float = 0.0,
    T_dew_actual:        float = 0.0,
    T_feed_actual:       float = 0.0,
    T_bp_perm_actual:    float = 0.0,
    T_cold_out_actual:   float = 0.0,
) -> MemSimulationResult:
    """計算不能な条件のときに返すペナルティ結果。

    Parameters
    ----------
    reason : str
        発火条件のラベル。MemEquipmentData の penalty_reason に格納し、
        run_one_pass._compute_mem_shortfall が連続 shortfall を計算する際に参照。
    P_H_actual, P_feed_actual : float
        ph_le_pfeed / ph_le_pl 経路で実 P_H・実 feed.P_in を保持。
    T_dew_actual, T_feed_actual : float
        vapor_condensed 経路で feed 露点・feed.T_in を保持。
    T_bp_perm_actual, T_cold_out_actual : float
        bp_le_cold_out 経路で透過ガス泡点・冷却水出口温度を保持。

    設計判断 (2026-05-22): 旧版は引数なしの silent penalty で、BO は「どの方向に
    動かせば feasible に出るか」のシグナルを得られず 80% の trial が無方向で死んでいた。
    PSA / Reactor で先に導入したパターン (psa_system.py:_penalty_result) を Mem にも
    展開する。
    """
    zero_stream = MemRetentateStream(0.0, 0.0, 0.0, 0.0)
    zero_prod   = MemProductStream(0.0, 0.0, 0.0, 0.0)
    eq = MemEquipmentData(
        A_vap=0.0, Pg_vap=0.0, Q_vap_kW=0.0,
        W_feed_kW=0.0, Pg_feed=0.0,
        A_mem=0.0, n_modules=0, Pg_mem=0.0,
        W_prod_kW=0.0, Pg_prod=0.0,
        A_cond=0.0, Pg_cond=0.0, Q_cond_kW=0.0,
        CAPEX_total=_PENALTY,
        penalty_reason=reason,
        P_H_actual_Pa=P_H_actual,
        P_feed_actual_Pa=P_feed_actual,
        T_dew_actual_K=T_dew_actual,
        T_feed_actual_K=T_feed_actual,
        T_bp_perm_actual_K=T_bp_perm_actual,
        T_cold_out_actual_K=T_cold_out_actual,
    )
    return MemSimulationResult(zero_stream, zero_prod, eq, 0.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# 内部ヘルパー
# ---------------------------------------------------------------------------

def _h_mol(T: float, P: float, z_C3H6: float, phase: str) -> float:
    """
    C3H6/C3H8 混合物のモルエンタルピー [J/mol]。
    基準: T_ref=298.15K の理想気体。

    H = H_ig(T_ref→T) + H^r(T, P, x, Z_phase)
    """
    x = [z_C3H6, 1.0 - z_C3H6]
    H_ig = sum(
        x[i] * _thermo.calc_enthalpy_change(_KEYS[i], _T_REF, T)
        for i in range(2)
    )
    Z   = z_factor(T, P, x, _KEYS, phase)
    H_r = residual_enthalpy(T, P, x, _KEYS, Z)
    return H_ig + H_r


def _lmtd(dT1: float, dT2: float) -> float:
    """対数平均温度差 LMTD [K]。dT1 ≠ dT2 の場合のみ厳密計算。"""
    if dT1 <= 0 or dT2 <= 0:
        return float('nan')
    if abs(dT1 - dT2) < 1e-3:
        return (dT1 + dT2) / 2.0
    return (dT1 - dT2) / math.log(dT1 / dT2)


# ---------------------------------------------------------------------------
# ユニット 1: 気化器
# ---------------------------------------------------------------------------

def _vaporizer(F_feed_mols: float, z_C3H6: float,
               T_in: float, P_in: float,
               fixed: MemFixedParams):
    """
    気化器: 液フィードを全量ガス化。

    Returns
    -------
    T_vap_out : 出口ガス温度 [K]（露点 + 過熱度）
    Q_vap_kW  : 必要加熱量 [kW]
    A_vap     : 伝熱面積 [m²]
    """
    T_dew = dew_point_T(P_in, [z_C3H6, 1.0 - z_C3H6], _KEYS)
    if math.isnan(T_dew):
        return float('nan'), float('nan'), float('nan')
    T_vap_out = T_dew + fixed.T_vap_superheat

    if T_vap_out >= fixed.T_hot:
        warnings.warn(
            f"気化器: T_vap_out={T_vap_out:.1f}K が T_hot={fixed.T_hot:.1f}K 以上です。"
            " T_hot を上げてください。",
            UserWarning, stacklevel=3,
        )
        return T_vap_out, float('nan'), float('nan')

    H_liq = _h_mol(T_in,      P_in, z_C3H6, 'liquid')
    H_vap = _h_mol(T_vap_out, P_in, z_C3H6, 'vapor')
    Q_vap = F_feed_mols * (H_vap - H_liq)    # [J/s = W]
    Q_vap_kW = Q_vap / 1e3

    dT1 = fixed.T_hot - T_in
    dT2 = fixed.T_hot - T_vap_out
    lmtd = _lmtd(dT1, dT2)
    A_vap = Q_vap_kW / (fixed.U_vap * lmtd) if lmtd > 0 else float('nan')

    return T_vap_out, Q_vap_kW, A_vap


# ---------------------------------------------------------------------------
# ユニット 3: 膜モジュール（クロスフロー ODE）
# ---------------------------------------------------------------------------

def _y_local(x: float, alpha: float, gamma: float) -> float:
    """
    クロスフローモデルにおける局所透過組成 y_local [-] を返す。

    【仮定】完全クロスフロー: 透過ガスは膜面を離れた瞬間に系外へ出ると仮定し、
    フィード側のみ流れ方向に組成変化する。スパイラル型モジュールの標準近似。

    x × P_H と y_local × P_L の分圧差が推進力。
    (1-alpha)*gamma * y² + [(alpha-1)*(x+gamma)+1] * y - alpha*x = 0
    の物理根（正かつ ≤ 1 の根）を返す。
    """
    # ID-06: alpha=1 は選択性なし → 線形退化を避けて物理的に正しい y=x を直返し
    if abs(alpha - 1.0) < 1e-10:
        return x
    a = (1.0 - alpha) * gamma
    b = (alpha - 1.0) * (x + gamma) + 1.0
    c = -alpha * x
    disc = max(0.0, b**2 - 4.0*a*c)
    denom = -b - math.sqrt(disc)
    if abs(denom) < 1e-30:
        return x
    return max(0.0, min(1.0, (2.0*c) / denom))


def _membrane_ode(F_C3H6_feed: float, F_C3H8_feed: float,
                  P_H: float, P_L: float, A_mem: float,
                  Q_A_SI: float, alpha: float):
    """
    クロスフロー膜モジュールの ODE 積分。

    dF_C3H6/dA = −Q_A × (x × P_H − y_local × P_L)
    dF_C3H8/dA = −Q_B × ((1−x) × P_H − (1−y_local) × P_L)

    Returns
    -------
    F_ret_C3H6, F_ret_C3H8 : 非透過ガス流量 [mol/s]
    None, None              : 積分失敗時
    """
    gamma = P_L / P_H
    Q_B   = Q_A_SI / alpha

    def ode(A, F):
        fc  = max(F[0], 1e-12)
        fa  = max(F[1], 1e-12)
        x   = fc / (fc + fa)
        y   = _y_local(x, alpha, gamma)
        J_c = Q_A_SI * (x * P_H - y * P_L)
        J_a = Q_B    * ((1.0-x) * P_H - (1.0-y) * P_L)
        # 負フラックスをクリップ（逆拡散なし）。[0,0] 返却はやめ event に任せる
        return [-max(J_c, 0.0), -max(J_a, 0.0)]

    # ID-02: フラックス枯渇（min(J_c,J_a) が 0 を下回った瞬間）に積分を打ち切る。
    # [0,0] を返し続けると Radau が超大ステップを踏みつつ終了しない場合があるため
    # terminal event で明示的に終了させる。
    def _event_no_flux(A, F):
        fc  = max(F[0], 1e-12)
        fa  = max(F[1], 1e-12)
        x   = fc / (fc + fa)
        y   = _y_local(x, alpha, gamma)
        J_c = Q_A_SI * (x * P_H - y * P_L)
        J_a = Q_B    * ((1.0-x) * P_H - (1.0-y) * P_L)
        return min(J_c, J_a)  # 0 を下から上から切るとき (direction=-1) に停止

    _event_no_flux.terminal  = True
    _event_no_flux.direction = -1

    # 設計判断 (2026-05-08, profile 結果反映):
    # 元値 rtol=1e-5, atol=1e-8 は scipy デフォルトの 100倍厳しい設定。
    # profile で全体時間の 9% を占めていた。膜の透過量積分は反応器ほど厳密性が
    # 要らないため、reactor と統一して rtol=1e-4, atol=1e-7 に緩める。
    try:
        sol = solve_ivp(
            ode,
            t_span=(0.0, A_mem),
            y0=[F_C3H6_feed, F_C3H8_feed],
            method='Radau',
            rtol=1e-4,
            atol=1e-7,
            # ID-01: max_step を設定してステップ数を ~200 に抑えフリーズを防止。
            # 低駆動力・巨大 A_mem 条件で無限に細かいステップを踏む問題を回避する。
            max_step=max(A_mem / 200.0, 0.1),
            events=_event_no_flux,
        )
    except Exception:
        return None, None

    # status 0: t_span 末端まで完走 / status 1: terminal event 発火（正常終了）
    if sol.status == -1:
        return None, None

    return float(sol.y[0, -1]), float(sol.y[1, -1])


# ---------------------------------------------------------------------------
# ユニット 5: 製品冷却器
# ---------------------------------------------------------------------------

def _condenser(F_perm_mols: float, y_C3H6: float,
               T_in: float, P_dist: float,
               fixed: MemFixedParams):
    """
    製品冷却器: 圧縮後の透過ガスを飽和液まで冷却・凝縮。

    設計判断 (2026-05-08, HI 用途追加): 顕熱（ガス T_in → T_bp）と潜熱
    （T_bp での凝縮）を分離して返す。HI 抽出関数が温度範囲を持つ顕熱
    ストリームと潜熱ストリームを別々に扱えるようにするため。
    Q_cond_kW (合計) は従来互換のため残す。

    Returns
    -------
    T_bp                 : 出口（泡点）温度 [K]
    Q_cond_kW            : 必要冷却量 [kW] = 顕熱 + 潜熱
    A_cond               : 伝熱面積 [m²]
    Q_cond_sensible_kW   : ガス顕熱分 [kW] (T_in → T_bp)
    Q_cond_latent_kW     : 凝縮潜熱分 [kW] (T_bp で凝縮)
    """
    T_bp = bubble_point_T(P_dist, [y_C3H6, 1.0 - y_C3H6], _KEYS)
    if math.isnan(T_bp):
        return float('nan'), float('nan'), float('nan'), float('nan'), float('nan')

    # ID-09: 圧縮機出口が既に泡点以下 → エンタルピー差が負になり Q_cond < 0
    if T_in <= T_bp:
        warnings.warn(
            f"製品冷却器: 圧縮機出口温度 {T_in:.1f} K が泡点 {T_bp:.1f} K 以下です。"
            " 透過ガスが既に凝縮状態のため冷却器モデルが無効です。",
            UserWarning, stacklevel=3,
        )
        return T_bp, float('nan'), float('nan'), float('nan'), float('nan')

    # 顕熱: ガス T_in → ガス T_bp
    H_gas_in    = _h_mol(T_in, P_dist, y_C3H6, 'vapor')
    H_gas_at_bp = _h_mol(T_bp, P_dist, y_C3H6, 'vapor')
    # 潜熱: ガス T_bp → 液 T_bp (温度同一)
    H_liq_out   = _h_mol(T_bp, P_dist, y_C3H6, 'liquid')

    Q_cond_sensible_W = F_perm_mols * (H_gas_in    - H_gas_at_bp)   # 顕熱 [W]
    Q_cond_latent_W   = F_perm_mols * (H_gas_at_bp - H_liq_out  )   # 潜熱 [W]
    Q_cond_sensible_kW = Q_cond_sensible_W / 1e3
    Q_cond_latent_kW   = Q_cond_latent_W   / 1e3
    Q_cond_kW          = Q_cond_sensible_kW + Q_cond_latent_kW

    # 向流熱交換器 LMTD
    # ガス側: T_in → T_bp,  冷却水側: T_cold_in → T_cold_out
    dT1 = T_in - fixed.T_cold_out     # ガス入口端
    dT2 = T_bp - fixed.T_cold_in      # ガス出口端（液出口端）
    lmtd = _lmtd(dT1, dT2)
    A_cond = Q_cond_kW / (fixed.U_cond * lmtd) if lmtd > 0 else float('nan')

    return T_bp, Q_cond_kW, A_cond, Q_cond_sensible_kW, Q_cond_latent_kW


# ---------------------------------------------------------------------------
# メイン関数
# ---------------------------------------------------------------------------

def simulate_membrane_system(
    design: MemDesignVars,
    feed: MemFeedStream,
    fixed: MemFixedParams,
) -> MemSimulationResult:
    """
    膜分離システム全体をシミュレーションする。

    Parameters
    ----------
    design : MemDesignVars   設計変数 (P_H, P_L, A_mem, P_dist)
    feed   : MemFeedStream   入力ストリーム（前段蒸留塔底液）
    fixed  : MemFixedParams  固定パラメータ

    Returns
    -------
    MemSimulationResult
    """
    # ---- 入力バリデーション ----
    # ID-08: 非正値チェック
    if feed.P_in <= 0 or feed.T_in <= 0 or design.P_dist <= 0:
        return _penalty_result(reason='invalid_input')
    if design.P_H <= design.P_L:
        warnings.warn("P_H <= P_L: 膜の駆動力がありません。", UserWarning, stacklevel=2)
        return _penalty_result(
            reason='ph_le_pl', P_H_actual=design.P_H, P_feed_actual=design.P_L,
        )
    if design.A_mem <= 0 or design.P_H <= 0 or design.P_L <= 0:
        return _penalty_result(reason='invalid_design')
    # ID-04: 製品圧縮機が減圧方向になる
    if design.P_dist <= design.P_L:
        warnings.warn(
            f"P_dist={design.P_dist/1e5:.2f} bar <= P_L={design.P_L/1e5:.2f} bar:"
            " 製品圧縮機が減圧方向になります。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(reason='pdist_le_pl')
    # ID-05: フィード圧縮機が減圧方向になる
    if design.P_H <= feed.P_in:
        warnings.warn(
            f"P_H={design.P_H/1e5:.2f} bar <= P_in={feed.P_in/1e5:.2f} bar:"
            " フィード圧縮機が減圧方向になります。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(
            reason='ph_le_pfeed', P_H_actual=design.P_H, P_feed_actual=feed.P_in,
        )
    if feed.F_C3H6 < 0 or feed.F_C3H8 < 0:
        return _penalty_result(reason='neg_feed')

    F_total_feed = feed.F_C3H6 + feed.F_C3H8
    if F_total_feed <= 0:
        return _penalty_result(reason='zero_feed')

    z_C3H6_feed = feed.F_C3H6 / F_total_feed  # 供給液中 C3H6 分率

    # ID-03: 液/ガス相整合性チェック
    # dew_point_T は収束失敗時に nan を返す（E-3 修正後）ため try/except は不要
    T_dew_feed = dew_point_T(feed.P_in, [z_C3H6_feed, 1.0 - z_C3H6_feed], _KEYS)
    if math.isnan(T_dew_feed):
        return _penalty_result(reason='dew_nan')
    if not fixed.vapor_feed and feed.T_in >= T_dew_feed:
        warnings.warn(
            f"feed.T_in={feed.T_in:.1f}K が露点 {T_dew_feed:.1f}K 以上です。"
            " 液相フィードを前提とするモデル(vapor_feed=False)と矛盾します。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(
            reason='liquid_vaporized',
            T_dew_actual=T_dew_feed, T_feed_actual=feed.T_in,
        )
    if fixed.vapor_feed and feed.T_in < T_dew_feed:
        warnings.warn(
            f"feed.T_in={feed.T_in:.1f}K が露点 {T_dew_feed:.1f}K 未満です。"
            " ガス相フィードを前提とするモデル(vapor_feed=True)と矛盾します。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(
            reason='vapor_condensed',
            T_dew_actual=T_dew_feed, T_feed_actual=feed.T_in,
        )

    # mol/s 換算（内部計算用）
    F_feed_mols = F_total_feed * 1000.0 / 3600.0   # [mol/s]
    Q_A_SI = fixed.Q_A_GPU * _GPU_SI                # [mol/(m²·s·Pa)]

    # ---- ユニット 1: 気化器（vapor_feed=True ならスキップ）----
    if fixed.vapor_feed:
        # ガスフィード: Vaporizer 不要。フィード温度をそのまま圧縮機入口とする。
        T_vap_out = feed.T_in
        Q_vap_kW  = 0.0
        A_vap     = 0.0
    else:
        try:
            T_vap_out, Q_vap_kW, A_vap = _vaporizer(
                F_feed_mols, z_C3H6_feed, feed.T_in, feed.P_in, fixed
            )
        except ImportError:
            raise
        except Exception:
            return _penalty_result(reason='vap_exception')
        # ID-10: Q_vap_kW だけでなく A_vap も nan チェック（LMTD 温度クロス時に nan になる）
        if math.isnan(Q_vap_kW) or math.isnan(A_vap):
            return _penalty_result(reason='vap_nan')

    # ---- ユニット 2: フィード圧縮機 ----
    try:
        T_feed_comp_out, W_feed_per_mol = compress_isentropic(
            T_vap_out, feed.P_in, design.P_H,
            [z_C3H6_feed, 1.0 - z_C3H6_feed], _KEYS,
            eta=fixed.eta_comp,
        )
    except ImportError:
        raise
    except Exception:
        return _penalty_result(reason='feed_comp_exception')
    W_feed_kW = F_feed_mols * W_feed_per_mol / 1e3   # [kW]

    # ---- ユニット 3: 膜モジュール ----
    F_C3H6_feed_mols = z_C3H6_feed * F_feed_mols
    F_C3H8_feed_mols = (1.0 - z_C3H6_feed) * F_feed_mols

    F_ret_C3H6_mols, F_ret_C3H8_mols = _membrane_ode(
        F_C3H6_feed_mols, F_C3H8_feed_mols,
        design.P_H, design.P_L, design.A_mem,
        Q_A_SI, fixed.alpha,
    )
    if F_ret_C3H6_mols is None:
        return _penalty_result(reason='ode_failure')

    # 非透過・透過流量（mol/s）
    F_ret_C3H6_mols = max(F_ret_C3H6_mols, 0.0)
    F_ret_C3H8_mols = max(F_ret_C3H8_mols, 0.0)
    F_perm_C3H6_mols = F_C3H6_feed_mols - F_ret_C3H6_mols
    F_perm_C3H8_mols = F_C3H8_feed_mols - F_ret_C3H8_mols
    F_perm_total_mols = max(F_perm_C3H6_mols + F_perm_C3H8_mols, 1e-12)

    # 透過・非透過組成
    y_C3H6 = F_perm_C3H6_mols / F_perm_total_mols
    F_ret_total_mols = F_ret_C3H6_mols + F_ret_C3H8_mols
    x_ret_C3H6 = (F_ret_C3H6_mols / max(F_ret_total_mols, 1e-12))
    stage_cut = (F_perm_C3H6_mols + F_perm_C3H8_mols) / F_feed_mols

    # ---- ユニット 4: 製品圧縮機 ----
    # 膜は等温操作と仮定: 透過ガスは T_feed_comp_out, P_L で出る
    try:
        T_prod_comp_out, W_prod_per_mol = compress_isentropic(
            T_feed_comp_out, design.P_L, design.P_dist,
            [y_C3H6, 1.0 - y_C3H6], _KEYS,
            eta=fixed.eta_comp,
        )
    except ImportError:
        raise
    except Exception:
        return _penalty_result(reason='prod_comp_exception')
    W_prod_kW = F_perm_total_mols * W_prod_per_mol / 1e3   # [kW]

    # ---- ユニット 5: 製品冷却器 ----
    try:
        T_bp_perm, Q_cond_kW, A_cond, Q_cond_sens_kW, Q_cond_lat_kW = _condenser(
            F_perm_total_mols, y_C3H6,
            T_prod_comp_out, design.P_dist,
            fixed,
        )
    except ImportError:
        raise
    except Exception:
        return _penalty_result(reason='cond_exception')
    # ID-09/10: Q_cond_kW または A_cond が nan（圧縮機出口が泡点以下・LMTD クロス）
    if math.isnan(Q_cond_kW) or math.isnan(A_cond):
        return _penalty_result(reason='cond_nan')

    # Case A 制約: 製品の泡点が冷却水出口温度を下回ると温度クロスが発生し
    # 冷却水では凝縮できない。冷媒は使用しない設計判断のため、この条件は
    # 物理的に不成立として最適化から除外する（ペナルティ返却）。
    if T_bp_perm <= fixed.T_cold_out:
        warnings.warn(
            f"製品冷却器: 泡点 {T_bp_perm - 273.15:.1f}°C が冷却水出口温度 "
            f"{fixed.T_cold_out - 273.15:.1f}°C 以下です（P_dist が低すぎます）。",
            UserWarning, stacklevel=2,
        )
        return _penalty_result(
            reason='bp_le_cold_out',
            T_bp_perm_actual=T_bp_perm, T_cold_out_actual=fixed.T_cold_out,
        )

    # ---- kmol/h 換算（出力用）----
    to_kmolh = 3600.0 / 1000.0

    # ---- ゲージ圧変換 [barg] ----
    def _pg(P_pa: float) -> float:
        return P_pa / 1e5 - _ATM_BAR

    # ---- モジュール本数 ----
    n_modules = math.ceil(design.A_mem / fixed.A_per_module)

    # ---- CAPEX 推算（Bare Module Cost 法, R08-3.pdf） ----
    try:
        capex_vap       = 0.0 if fixed.vapor_feed else calc_he_capex_okuyen(A_vap)
        capex_comp_feed = calc_comp_capex_okuyen(W_feed_kW)
        capex_comp_prod = calc_comp_capex_okuyen(W_prod_kW)
        capex_cond      = calc_he_capex_okuyen(A_cond)
    except Exception:
        capex_vap = capex_comp_feed = capex_comp_prod = capex_cond = float('nan')
    # !仮置き — MEM_UNIT_PRICE_USD_PER_M2 が根拠文献未確定のため暫定値
    #   calc_mem_capex_okuyen 呼び出し時に UserWarning が発行される
    try:
        capex_mem = calc_mem_capex_okuyen(design.A_mem)
    except Exception:
        capex_mem = float('nan')
    _capex_sum = capex_vap + capex_comp_feed + capex_comp_prod + capex_cond + capex_mem
    # ID-07: 個別 CAPEX が nan（計算例外）のとき合計も nan になり最適化器がハングする
    capex_total = _capex_sum if not math.isnan(_capex_sum) else _PENALTY

    return MemSimulationResult(
        retentate=MemRetentateStream(
            F_C3H6 = F_ret_C3H6_mols * to_kmolh,
            F_C3H8 = F_ret_C3H8_mols * to_kmolh,
            T_out  = T_feed_comp_out,
            P_out  = design.P_H,
        ),
        product=MemProductStream(
            F_C3H6 = F_perm_C3H6_mols * to_kmolh,
            F_C3H8 = F_perm_C3H8_mols * to_kmolh,
            T_out  = T_bp_perm,
            P_out  = design.P_dist,
        ),
        equipment=MemEquipmentData(
            A_vap          = A_vap,
            Pg_vap         = _pg(feed.P_in),
            Q_vap_kW       = Q_vap_kW,
            W_feed_kW      = W_feed_kW,
            Pg_feed        = _pg(design.P_H),
            A_mem          = design.A_mem,
            n_modules      = n_modules,
            Pg_mem         = _pg(design.P_H),
            W_prod_kW      = W_prod_kW,
            Pg_prod        = _pg(design.P_dist),
            A_cond         = A_cond,
            Pg_cond        = _pg(design.P_dist),
            Q_cond_kW      = Q_cond_kW,
            CAPEX_vap      = capex_vap,
            CAPEX_comp_feed= capex_comp_feed,
            CAPEX_comp_prod= capex_comp_prod,
            CAPEX_cond     = capex_cond,
            CAPEX_mem      = capex_mem,
            CAPEX_total    = capex_total,
            T_vap_in_K     = feed.T_in,
            T_vap_out_K    = T_vap_out,
            T_cond_in_K    = T_prod_comp_out,
            T_cond_out_K   = T_bp_perm,
            Q_cond_sensible_kW = Q_cond_sens_kW,
            Q_cond_latent_kW   = Q_cond_lat_kW,
        ),
        stage_cut   = float(np.clip(stage_cut,   0.0, 1.0)),
        perm_purity = float(np.clip(y_C3H6,      0.0, 1.0)),
        ret_purity  = float(np.clip(x_ret_C3H6,  0.0, 1.0)),
    )
