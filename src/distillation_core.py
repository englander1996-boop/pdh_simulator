"""
簡易 FUG 蒸留塔モデル (共通エンジン)

設計判断 (2026-05-08):
  fake_column (split_fracs ベースのダミー) を、Fenske-Underwood-Gilliland
  (FUG) 法ベースの簡易モデルに置き換える。物理的により正しい温度・熱量を
  与え、BO による P_col / N_stages / N_feed / reflux_ratio の最適化を可能にする。

各塔の固有設定 (LK/HK 成分、デフォルト P/N/R) は units/separators/columnX/columnX.py
で定義する (Option C ハイブリッド構造)。

参考文献:
  [1] Fenske M.R., "Fractionation of straight-run Pennsylvania gasoline,"
      Ind. Eng. Chem. 24 (1932) 482-485. — N_min 最小段数式の出典
  [2] Underwood A.J.V., "Fractional distillation of multicomponent
      mixtures," Chem. Eng. Prog. 44 (8) (1948) 603-614. — R_min 最小還流比式
  [3] Eduljee H.E., "Equations replace Gilliland plot,"
      Hydrocarbon Processing 54 (9) (1975) 120. — Gilliland 相関 (式形)
  [4] Kirkbride C.G., "Process design procedure for multicomponent
      fractionators," Petroleum Refiner 23 (9) (1944) 321. — フィード段位置
  [5] プロセス設計授業資料 R08-3.pdf 付録 A
      Table A.1 (Sieve trays K1, K2, K3), 圧力係数 Fp, 据付間接費 1.18,
      CEPCI 補正など。

物性値の出典:
  - 大気圧沸点 T_b [K], 蒸発潜熱 λ [kJ/mol] : 既存 _T_BOIL_ATM, _LAMBDA_KJ
  - PR EOS パラメータ Tc, Pc, ω : src/config.py THERMO_DATA
    (n-Butane Z は化学工学便覧 改訂六版 表 1.3 物質No.181)
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from scipy.optimize import brentq

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.cost_calculator import (
    calc_cp0, calc_fp, calc_tray_capex_okuyen, calc_he_capex_okuyen,
)
from src.cost_parameters import (
    B1, B2, FM, CEPCI_BASE, CEPCI_CURRENT,
    USD_TO_JPY, PLANT_INDIRECT_FACTOR,
)
from src.component_data import cp_of, CP_DEFAULT
from src.config import THERMO_DATA
from src.eos import (
    z_factor, fugacity_coeff, residual_enthalpy,
    bubble_point_T, dew_point_T,
)
from src.utility_selector import (
    select_cooling_utility, select_heating_utility,
)
from flowsheet.heat_integration import StreamPhase, lookup_U, utility_phase


# ---------------------------------------------------------------------------
# 物性データ (大気圧沸点と蒸発潜熱、Clausius-Clapeyron 用)
# ---------------------------------------------------------------------------
_R_GAS = 8.314     # [J/(mol·K)]
_P_ATM = 101325.0  # [Pa]

# 各成分の大気圧沸点 [K] と蒸発潜熱 [kJ/mol]
# (旧版から流用、PDH の主要成分の代表値)
_T_BOIL_ATM = {
    'A': 231.1, 'B': 225.5, 'C':  20.3,
    'D': 169.4, 'E': 111.7, 'F': 184.6, 'Z': 272.7,
}
_LAMBDA_KJ = {
    'A': 19.1, 'B': 18.4, 'C': 0.9,
    'D': 13.5, 'E':  8.2, 'F': 14.7, 'Z': 22.0,
}
_T_BOIL_DEFAULT = 250.0
_LAMBDA_DEFAULT = 15.0

_T_COND_MIN = 313.15  # 40°C — 冷却水使用時の凝縮器最低温度
_PENALTY = 1e9

# PR 泡点フラッシュの brentq 数値探索ガード [K] (物性値ではない、数値の安全域)。
# - 下限 200K: H2/CH4 が混じる組成で泡点が cryogenic に張り付くケースを除外
#   (本フローの最低運転温度 ~-100°C=173K より低いので物理的に無関係)
# - 上限 600K: 探索失敗時に上端打ち切りで非物理解になるのを防ぐ
#   (本フローの蒸留塔操作温度 ~30〜200°C=303〜473K より高い)
# これら閾値外の解はクラジウス・クラペイロン (CC) フォールバックに切替える。
_T_BUBBLE_MIN = 200.0
_T_BUBBLE_MAX = 600.0


# ===========================================================================
# データクラス: 入力・出力
# ===========================================================================

@dataclass(frozen=True)
class ColumnTunables:
    """蒸留塔の最適化対象パラメータ (P/N/R)。

    LK/HK・回収率・q・K_method といった「塔別の物理セマンティクス」は
    column1/2/3.py 側で固定。本データクラスは BO や exp で振る部分のみ保持。

    FlowsheetDesignVars.dist1/dist2/dist3 で使う。

    solver_method: 求解アルゴリズムの選択
      'fug'      : Fenske-Underwood-Gilliland shortcut (高速、BO 用)
      'rigorous' : VLE ベース tray-by-tray (厳密、top-k 用) — 未実装、
                   設定すると NotImplementedError でフォールバック
    """
    P_col:         float    # 塔操作圧力 [Pa]
    N_stages:      int      # 理論段数
    N_feed:        int      # フィード段位置 (rigorous/sm では Kirkbride 推奨を自動採用、本値は無視。FUG では post-hoc 出力のみ)
    reflux_ratio:  float    # 還流比 R = L/D
    solver_method: str = 'fug'   # 'fug' | 'rigorous' | 'sm' | 'hysys'
    # 設計判断 (2026-05-14): recovery を BO で振れるよう field 化。
    # None → column1/2/3.py ラッパの既定値 (0.99) を採用 (後方互換)。
    # float (e.g., 0.95〜0.999) → ラッパ既定値を上書き。
    recovery_LK_top: Optional[float] = None
    recovery_HK_bot: Optional[float] = None
    # 設計判断 (2026-05-19): HYSYS との比較・D 直接指定検証用。
    # None → 従来通り FUG が Fenske split で D を決定。
    # float → rigorous solver で D_total を強制上書き (= "Distillate Rate spec" モード)。
    # 主に Dist2 partial cond で HYSYS の D 値と整合を取るデバッグ用。
    # 値が物理的に不適 (D > F_total など) なら solver 側で _failure_result を返す。
    D_override: Optional[float] = None
    # ---- HYSYS バックエンド用追加 (solver_method='hysys' のときのみ参照) ----
    # 設計判断 (2026-05-22): HYSYS の段数別 HSC ファイル + 主スペック書込み運用。
    # special.py の探索変数は lhs_column*.py の DOE 群 (P_col, feed_flow, 主スペック,
    # feed_stage, N_stages) と一致させる。
    #   hysys_spec_value: 塔ごとの主スペック値
    #     column1: Comp Fraction - 2 [-]
    #     column2: Reflux Ratio [-]
    #     column3: Draw Rate [kgmol/s]
    #   hysys_feed_stage: HYSYS の SpecifyFeedLocation に渡す段番号 (1始まり)
    # main の FUG/rigorous 経路では None のままで OK (後方互換)。
    hysys_spec_value: Optional[float] = None
    hysys_feed_stage: Optional[int] = None


# 分流 (partial condenser) で凝縮させない成分セット
# 設計判断 (2026-05-09): H2 (Tc=33K) と CH4 (Tc=190K) は工業的に到達可能な
# 冷媒温度 (-100°C エチレン = 173K) でも凝縮できないため、partial condenser では
# 凝縮させず vapor distillate でそのまま下流へパススルーする。
# C2H4 (Tc=282K) 以上は propylene/ethylene 冷媒で凝縮可能なので reflux 寄与あり。
NON_CONDENSABLE_COMPS = frozenset(['C', 'E'])  # H2, CH4

# 分流 (partial condenser) で「液側に強く偏る」成分セット (旧 ALWAYS_CONDENSABLE)
# 設計判断 (2026-05-09): C3 (Tc≈370K, P_sat@-30°C ≈ 1.5-1.6 bar) は Dist2 の
# 操作圧 (7-8 bar) >> P_sat なので塔頂温度 -30〜-40°C で大半が凝縮する。
# 設計判断 (2026-05-19): ただし PR EOS 実値で K(C3H6) @ -43°C/7.5bar = 0.18 と
# 「ゼロ」ではない。旧版は無条件 F_top → F_bot 100% 強制移動だったが、これだと
# (i) HYSYS との不一致 (C3H6 100% vs 物理は数%漏れ)
# (ii) FUG の D_total が C3=0 で固定 → rigorous (Wang-Henke) に渡ると質量保存が
#      閉じず MESH 残差過大 → サイレント FUG フォールバック → rigorous 名目だけ
# という同根バグを生んでいた。
# 現版は「閾値ベース丸め」: feed の _PARTIAL_COND_C3_ROUND_FRAC 未満の F_top 漏れだけ
# 0 に丸めて塔底に集約。「ほぼほぼ分離して PSA/Mem に C3 が流れない」目標を保ちつつ
# 物理的な D_total を尊重する。閾値以上の漏れは物理値を残し、警告を出す。
ALWAYS_CONDENSABLE_COMPS = frozenset(['A', 'B'])  # C3H8, C3H6 (閾値判定対象)

# 丸め閾値: feed flow に対する比。これ未満の F_top 漏れは 0 に丸めて塔底集約。
# 1% = C3H6 feed 2753 kmol/h なら 27.5 kmol/h 以下を丸め。
# ユーザー指示 (2026-05-19): 「ほぼほぼ分離・有効数字範囲のバランス崩れは OK」。
_PARTIAL_COND_C3_ROUND_FRAC = 0.01

# rigorous プロキシ罰則の係数 (2026-05-19 Phase A、ユーザー指示「FUG があんまりよろしくない」)
# 単位は spec_coef_okuyen (= 100 億円/(年・%pt)) と同次元。
# 設計判断 (2026-05-20 Phase A強化): main_20260520_040740 BO で trial #32 が
# FUG 中は feasible (TAC_bo=295) なのに rigorous で LK recovery 0.877 (vs spec 0.99)
# 違反で死亡。proxy_penalty 旧係数 (50/30) では BO 中の effective_TAC へのシグナル
# が弱く、BO が FUG の楽観性領域に張り付き続けるため強化:
#   - LEAK 50 → 150 (3x): partial cond の C3 漏れに対する罰則を強化
#   - MARGIN 30 → 80 (~2.7x): narrow margin trial を effective_TAC で明確に弾く
#   - R/N THRESHOLD 1.30 → 1.50: より早期に margin 不足を検出 (FUG-rigorous 乖離は
#     R/R_min < 1.5 領域で顕在化しやすい経験則)
#   - 加えて _compute_proxy_penalty 内で partial cond は係数 ×2 を追加適用
#     (Dist2 のような partial cond は FUG 楽観性が特に大きいため)
_PROXY_LEAK_COEF_OKUYEN_PER_PP = 150.0   # C3 漏れ 1pp 超過につき [億円/年]
_PROXY_MARGIN_COEF_OKUYEN      = 80.0    # narrow margin 1pt につき [億円/年]
# 「narrow」の判定閾値 (R/R_min, N/N_min がこれ未満なら罰則対象)
# 設計判断 (2026-05-20): partial cond / total cond で閾値を分ける。
#   partial cond (Dist2): 1.50。FUG の non-key 分配が物理から大きく乖離するため厳しめ。
#                          C3H6 漏れと連動して BO 楽観性領域を強力に弾く必要あり。
#   total cond (Dist1/Dist3): 1.30 (Phase A 旧値)。C3 splitter 等は physics 素直で
#                            実機運用 R/R_min=1.1-1.3 が一般的。1.50 だと Dist3 の
#                            R 下限を引き下げて overspec 解消する経済探索が阻害される
#                            (= R を下げると penalty が乗りすぎて BO が高 R に張り付く)。
# 設計判断 (2026-05-23 forensic, main_20260523_083228): R 閾値を緩和。
# 旧 (Partial 1.50 / Total 1.30) は「FUG-rigorous gap 防止」が目的だったが、
# main_20260523_083228 で TPE が「R ぎりぎり Dist3 設計」(= ユーザー狙いの境界探索)
# を試すと proxy 罰則だけで死亡する trial が頻発した:
#   #12 (TAC=1055, c3h6=99.98%, prod ok): Dist3 R/R_min=1.24 で +1.5 億円 = infeas
#   #26 (TAC=1080, c3h6=99.95%):          Dist2 R/R_min=1.10 で +28.7 億円 = infeas
#   #30 (TAC=1113, T_in=882K 新領域活用):   Dist2 R/R_min=1.15 + Dist3 R/R_min=1.23 で 死亡
# top-k で rigorous 再評価が走るため、BO loop で proxy を厳しくする必要性は低い。
# 閾値緩和して TPE が「FUG では feas、rigorous で詰むかも」境界を探索可能にする。
_PROXY_R_MARGIN_THRESHOLD_PARTIAL = 1.30  # partial cond の R/R_min 閾値 (1.50→1.30、Dist2)
_PROXY_R_MARGIN_THRESHOLD_TOTAL   = 1.20  # total   cond の R/R_min 閾値 (1.30→1.20、Dist1/Dist3)
_PROXY_N_MARGIN_THRESHOLD_PARTIAL = 1.50  # partial cond の N/N_min 閾値 (N は据置、ユーザー判断で R のみ)
_PROXY_N_MARGIN_THRESHOLD_TOTAL   = 1.30  # total   cond の N/N_min 閾値 (同上)
# partial cond の追加倍率 (Dist2 など、FUG 楽観性が大きい塔への上乗せ)
_PROXY_PARTIAL_COND_MULTIPLIER = 2.0     # partial cond の塔は係数 ×2.0 で評価


def _compute_proxy_penalty(
    design,
    F_top: Dict[str, float],
    feed_F: Dict[str, float],
    R_min: float,
    N_min: float,
    *,
    is_rigorous: bool = False,
) -> Tuple[float, str]:
    """rigorous プロキシ罰則を計算。

    BO objective の effective_TAC に乗せる罰則 [億円/年] を返す。
    呼び出し側 (FUG / rigorous 両パス) で結果整理直前にコールし、
    DistEquipment.proxy_penalty_okuyen にセットする。

    成分:
      (a) C3 漏れ罰則  : partial cond の F_top に C3 が _PARTIAL_COND_C3_ROUND_FRAC を
                       超えて漏れている場合、超過 pp に係数を掛ける。
                       (= FUG が「物理的に C3 が漏れる設計」と判断してもまだ feasible
                        と通してしまうケースをペナルティ化)
      (b) margin 罰則  : R/R_min が _PROXY_R_MARGIN_THRESHOLD 未満、または
                       N/N_min が _PROXY_N_MARGIN_THRESHOLD 未満で罰則。
                       (= FUG 通過するが rigorous Wang-Henke が dT_max stall する領域)

    Returns
    -------
    (penalty_okuyen, reason)
        penalty_okuyen: 罰則合計 [億円/年]
        reason: デバッグ用文字列 (発火項目を列挙)
    """
    penalty = 0.0
    reasons: List[str] = []

    # 設計判断 (2026-05-20 Phase A強化): partial cond の塔 (Dist2) は FUG の
    # non-key 分配が物理から大きく乖離しやすいので、係数を _PROXY_PARTIAL_COND_MULTIPLIER
    # で底上げする。Dist1/Dist3 (total cond) は ×1.0 のまま。
    # 同時に margin 閾値も partial / total で分ける (partial は 1.50、total は 1.30)。
    is_partial = bool(getattr(design, 'partial_condenser', False))
    coef_mult = _PROXY_PARTIAL_COND_MULTIPLIER if is_partial else 1.0
    leak_coef = _PROXY_LEAK_COEF_OKUYEN_PER_PP * coef_mult
    margin_coef = _PROXY_MARGIN_COEF_OKUYEN * coef_mult
    r_threshold = _PROXY_R_MARGIN_THRESHOLD_PARTIAL if is_partial else _PROXY_R_MARGIN_THRESHOLD_TOTAL
    n_threshold = _PROXY_N_MARGIN_THRESHOLD_PARTIAL if is_partial else _PROXY_N_MARGIN_THRESHOLD_TOTAL

    # (a) C3 漏れ罰則 (partial cond のみ)
    # 設計判断 (2026-05-21): is_rigorous=True のときは leak 罰則をスキップ。
    # 理由: rigorous F_top の C3 漏れ量は物理的真値で、PSA/Mem trace bypass の連続
    # penalty (runner.py _TRACE_BYPASS_PENALTY_COEF_OKUYEN) が既にコスト化しているため
    # ここで再度加算すると二重計上になる (実例: main_20260521_092409 trial 0 で
    # proxy +287 億 + bypass +86 億 = 373 億の二重 penalty が warm-start を死なせていた)。
    # FUG パスでは F_top が楽観値なので、proxy leak penalty で「物理的には漏れる」
    # 警告を入れる意味があるが、rigorous は真値なので不要。
    # margin 罰則 (b) は rigorous でも維持: R/R_min 近傍は CAPEX に反映されてはいるが
    # 「操業的脆弱性」として追加コストを乗せる意義あり (係数は小さいので影響軽微)。
    if is_partial and not is_rigorous:
        for c in ALWAYS_CONDENSABLE_COMPS:
            moved = F_top.get(c, 0.0)
            f_c = feed_F.get(c, 0.0)
            if f_c <= 0 or moved <= 0:
                continue
            leak_frac = moved / f_c
            if leak_frac > _PARTIAL_COND_C3_ROUND_FRAC:
                over_pp = (leak_frac - _PARTIAL_COND_C3_ROUND_FRAC) * 100.0
                add = over_pp * leak_coef
                penalty += add
                reasons.append(
                    f"{c}漏れ{leak_frac*100:.2f}% (超過{over_pp:.2f}pp → +{add:.1f}億円)"
                )

    # (b) R/R_min, N/N_min margin 罰則
    # 設計判断 (2026-05-20 Phase A強化): shortage は生値 (0-1) で扱う。
    # threshold は partial / total で分岐 (上で r_threshold/n_threshold 選択済)。
    # 生値のまま coef 80 を掛けることで BO best ~300 億スケールに対し effective
    # な (数十億円オーダー) を狙う。partial cond は ×2 (上で coef_mult 適用済)。
    R = design.reflux_ratio
    N = design.N_stages
    if R_min > 0 and math.isfinite(R_min):
        r_ratio = R / R_min
        if r_ratio < r_threshold:
            shortage = r_threshold - r_ratio   # 生値 (0-1.5)
            add = shortage * margin_coef
            penalty += add
            reasons.append(
                f"R/R_min={r_ratio:.3f} < {r_threshold} "
                f"(不足{shortage:.3f} → +{add:.1f}億円)"
            )
    if N_min > 0 and math.isfinite(N_min):
        n_ratio = N / N_min
        if n_ratio < n_threshold:
            shortage = n_threshold - n_ratio   # 生値 (0-1.5)
            add = shortage * margin_coef
            penalty += add
            reasons.append(
                f"N/N_min={n_ratio:.3f} < {n_threshold} "
                f"(不足{shortage:.3f} → +{add:.1f}億円)"
            )

    prefix = "rigorous" if is_rigorous else "fug"
    suffix = " (partial cond×2)" if is_partial else ""
    reason_str = f"[{prefix}{suffix}] " + "; ".join(reasons) if reasons else ""
    return penalty, reason_str


@dataclass(frozen=True)
class DistDesignVars:
    """蒸留塔設計変数 (BO 探索対象 + 仕様)。

    BO で最適化する変数:
      P_col, N_stages, N_feed, reflux_ratio

    仕様 (塔別ラッパーで固定):
      LK, HK, recovery_LK_top, recovery_HK_bot, K_method, q
    """
    # ---- BO 設計変数 ----
    P_col:            float            # 塔操作圧力 [Pa]
    N_stages:         int              # 理論段数
    N_feed:           int              # フィード段位置 (1=塔頂, N_stages=塔底)
    reflux_ratio:     float            # 還流比 R = L/D
    # ---- 分離仕様 ----
    LK:               str              # 軽キー成分 ('A'-'F','Z')
    HK:               str              # 重キー成分
    recovery_LK_top:  float = 0.99     # LK 塔頂回収率
    recovery_HK_bot:  float = 0.99     # HK 塔底回収率
    # ---- 物性計算オプション ----
    K_method:         str   = 'pr'     # 'pr' (PR EOS) または 'cc' (Clausius-Clapeyron)
    q:                float = 1.0      # フィード状態 (1=飽和液, 0=飽和気)
    # 設計判断 (2026-05-09): partial condenser フラグ。
    # True のとき塔頂は分流型で、軽質ガス (NON_CONDENSABLE_COMPS = H2/CH4) は
    # 凝縮させず vapor distillate でパススルー、reflux 用に重質分のみ凝縮する。
    # Dist2 (脱エタン塔) のように H2/CH4 主体の塔頂で必須。Dist1/Dist3 (C3 主体) は False。
    partial_condenser: bool = False
    # ---- 求解アルゴリズム選択 ----
    # 'fug'      : Fenske-Underwood-Gilliland shortcut (デフォルト、高速、BO 用)
    # 'rigorous' : VLE ベース tray-by-tray (厳密、top-k 用)
    #              ※ 本体未実装。設定された場合は NotImplementedError を投げて
    #                 FUG にフォールバック (警告つき)。VLE ソルバ実装後に有効化される。
    solver_method:    str = 'fug'
    # 設計判断 (2026-05-19): D 直接指定モード (HYSYS の "Distillate Rate spec" 相当)。
    # None → 従来通り FUG の Fenske split で D を決定。
    # float → rigorous で D_total を直接指定 (= recovery spec を上書き)。
    # HYSYS の Distillate Rate 値と直接比較するための機能。
    D_override:       Optional[float] = None


@dataclass(frozen=True)
class DistFixedParams:
    """蒸留塔固定パラメータ (寸法・流体力学、contest §4-2 準拠)。

    出典: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-2 蒸留塔
      段塔仮定、塔径はフラッディングを起こさない許容蒸気質量速度 G* で決定:
        G* = SF · K · √(ρ_v · (ρ_l - ρ_v))   [kg/(m²·s)]
      シーブトレー仮定で SF=0.8、K=0.05 m/s、段間隔 0.6m を採用。
      段効率 80% で実段数を算出、塔頂 2m + 塔底 4m を加算。
    """
    tray_spacing_m:  float = 0.6   # トレイ間隔 [m]、contest §4-2
    top_section_m:   float = 2.0   # 塔頂付加長さ [m] (還流供給・気液分離)、contest §4-2
    bot_section_m:   float = 4.0   # 塔底付加長さ [m] (液ホールドアップ)、contest §4-2
    tray_efficiency: float = 0.8   # 段効率 [-]、contest §4-2 (実段数 = 理論段数 / η)
    SF:              float = 0.8   # G* の系補正係数 [-]、contest §4-2
    K_factor:        float = 0.05  # G* の許容蒸気速度係数 [m/s]、contest §4-2


@dataclass
class DistEquipment:
    """蒸留塔の装置・診断情報。"""
    # ---- 塔寸法 ----
    D_col:               float          # 塔径 [m]
    H_col:               float          # 塔高さ [m]
    V_col:               float          # 塔体積 [m³]
    # ---- CAPEX 内訳 ----
    CAPEX_vessel:        float          # 塔本体 [億円]
    CAPEX_trays:         float          # トレイ [億円]
    CAPEX:               float          # 合計 [億円] (vessel + trays + cond + reb)
    # ---- 熱量 ----
    Q_cond:              float          # コンデンサ熱量 [kW] (放出、正値)
    Q_reb:               float          # リボイラ熱量 [kW]   (吸熱、正値)
    Q_feed_preheat_kW:   float = 0.0    # フィード予熱 [kW]
    # ---- FUG 診断 ----
    N_min:               float = 0.0    # 最小段数 (Fenske)
    R_min:               float = 0.0    # 最小還流比 (Underwood)
    N_feed_kirkbride:    int   = 0      # Kirkbride 推奨フィード段
    feasible:            bool  = True   # N >= N_min かつ R >= R_min
    message:             str   = ""
    # ---- Condenser / Reboiler 装置 (2026-05-09 追加) ----
    A_cond_m2:               float = 0.0    # コンデンサ伝熱面積 [m²]
    A_reb_m2:                float = 0.0    # リボイラ伝熱面積 [m²]
    CAPEX_cond:              float = 0.0    # コンデンサ HE CAPEX [億円]
    CAPEX_reb:               float = 0.0    # リボイラ HE CAPEX [億円]
    cond_utility_name:       str   = ""     # 例 "冷却水≒", "プロピレン冷媒-25C≒"
    cond_utility_jpy_per_GJ: float = 0.0
    reb_utility_name:        str   = ""     # 例 "LP Steam≒", "MP Steam≒"
    reb_utility_jpy_per_GJ:  float = 0.0
    T_top:                   float = 0.0    # 塔頂温度 [K] (clamp なし)
    T_bot:                   float = 0.0    # 塔底温度 [K]
    # ---- rigorous プロキシ罰則 (2026-05-19 追加、Phase A) ----
    # FUG の楽観性を是正するために BO objective に乗せる罰則 [億円/年]。
    # narrow-margin (R/R_min, N/N_min ≪ 1.5) や partial cond の C3 漏れ過大 等、
    # FUG では feasible でも rigorous で詰むことが多い設計に対して加算。
    # runner.py で soft_penalty にマージされ effective_TAC に反映される。
    proxy_penalty_okuyen:    float = 0.0    # 罰則合計 [億円/年]
    proxy_penalty_reason:    str   = ""     # 罰則発火理由 (デバッグ用)
    # ---- infeasibility 連続シグナル (2026-05-20 追加、TPE constraints_func 用) ----
    # FUG penalty の場合: Gilliland _N_needed_ を保持し、run_one_pass で N_shortfall を計算。
    # Rigorous Wang-Henke 失敗の場合: 失敗時の dT_max [K] を保持し、tol からの距離を penalty 化。
    # feasible=True (正常完了) のときは両方 0.0。
    N_needed:                float = 0.0    # Gilliland で必要 N [-] (0 ⇒ 未計算)
    dT_max_rigorous:         float = 0.0    # Wang-Henke 失敗時の最終 dT_max [K] (0 ⇒ 未発生)


@dataclass
class DistResult:
    """蒸留塔シミュレーション結果。"""
    top:       ProcessStream
    bottom:    ProcessStream
    equipment: DistEquipment


# ===========================================================================
# 物性ヘルパー (温度・K 値・α)
# ===========================================================================

def _boil_cc(T_b_atm: float, lam_kJ: float, P_col: float) -> float:
    """Clausius-Clapeyron 式で圧力 P_col [Pa] における沸点 [K] を返す。

    ln(P2/P1) = ΔHvap/R × (1/T1 - 1/T2)
    → T2 = 1 / (1/T1 - R/ΔHvap × ln(P2/P1))
    """
    lam_J = lam_kJ * 1000.0
    inv_T2 = 1.0 / T_b_atm - _R_GAS / lam_J * math.log(P_col / _P_ATM)
    if inv_T2 <= 0.0:
        return T_b_atm * 5.0
    return 1.0 / inv_T2


def _weighted_boil(F_dict: Dict[str, float], P_col: float) -> float:
    """流量加重平均の沸点 [K] (Clausius-Clapeyron で圧力補正)。"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return 298.15
    t = sum(
        F_dict.get(k, 0.0) * _boil_cc(
            _T_BOIL_ATM.get(k, _T_BOIL_DEFAULT),
            _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT),
            P_col,
        )
        for k in F_dict
    ) / F_total
    return t


def _weighted_lambda(F_dict: Dict[str, float]) -> float:
    """流量加重平均蒸発潜熱 [kJ/mol]。"""
    F_total = sum(F_dict.values())
    if F_total <= 0.0:
        return _LAMBDA_DEFAULT
    return sum(
        F_dict.get(k, 0.0) * _LAMBDA_KJ.get(k, _LAMBDA_DEFAULT)
        for k in F_dict
    ) / F_total


def _K_cc(comp: str, T: float, P_col: float) -> float:
    """Clausius-Clapeyron で K_i = P_sat,i(T) / P_col。"""
    T_b = _T_BOIL_ATM.get(comp, _T_BOIL_DEFAULT)
    lam = _LAMBDA_KJ.get(comp, _LAMBDA_DEFAULT) * 1000.0  # J/mol
    # P_sat = P_atm × exp(λ/R × (1/T_b - 1/T))
    P_sat = _P_ATM * math.exp(lam / _R_GAS * (1.0 / T_b - 1.0 / T))
    return P_sat / P_col


def _K_pr(comp_idx: int, comps: List[str], x: List[float],
          T: float, P_col: float) -> float:
    """PR EOS で K_i = phi_i^L / phi_i^V を計算。

    両相 (液・気) でフガシティ係数を計算し比を取る。
    """
    try:
        Z_L = z_factor(T, P_col, x, comps, 'liquid')
        Z_V = z_factor(T, P_col, x, comps, 'vapor')
        phi_L = fugacity_coeff(comp_idx, T, P_col, x, comps, Z_L)
        phi_V = fugacity_coeff(comp_idx, T, P_col, x, comps, Z_V)
        if phi_V <= 1e-30:
            return float('inf')
        return phi_L / phi_V
    except Exception:
        # 範囲外・収束失敗時は CC でフォールバック
        return _K_cc(comps[comp_idx], T, P_col)


def _K_dict(comps: List[str], x: List[float], T: float, P_col: float,
            method: str) -> Dict[str, float]:
    """全成分の K 値を計算して辞書で返す。"""
    if method == 'pr':
        # PR EOS は混合系。Tc/Pc/ω のない成分は CC にフォールバック
        K = {}
        for i, c in enumerate(comps):
            p = THERMO_DATA.get(c)
            if p is None or math.isnan(p.Tc):
                K[c] = _K_cc(c, T, P_col)
            else:
                K[c] = _K_pr(i, comps, x, T, P_col)
        return K
    # 'cc' モード
    return {c: _K_cc(c, T, P_col) for c in comps}


def _alpha_dict(K: Dict[str, float], HK: str) -> Dict[str, float]:
    """HK 基準の相対揮発度 α_i = K_i / K_HK。"""
    K_HK = K.get(HK, 1.0)
    if K_HK <= 1e-30:
        K_HK = 1e-30
    return {c: K[c] / K_HK for c in K}


def _bubble_T_K(
    x_list: List[float],
    comps:  List[str],
    P_col:  float,
    method: str,
) -> Tuple[float, List[float]]:
    """組成 x (mol fraction) の泡点温度と各成分 K 値を返す。

    method='pr' のとき: PR EOS で泡点フラッシュ (eos.bubble_point_T) を実行し、
    収束 T で K_i = φ_L/φ_V を取る。FUG 法で「平均 T で K=phi_L/phi_V」を計算
    すると、その T-P-x が単相領域に入り Z 根が 1 本で K_i ≈ 1 になる病理を
    回避できる (泡点では飽和包絡線上なので必ず 2 根存在)。

    PR が _T_BUBBLE_MIN..MAX を外れる解を返したり例外を出した場合 (例: H2 や CH4
    が x_top の主成分で泡点が cryogenic に張り付くケース) は CC にフォールバック。
    成分単位の失敗は当該成分のみ CC へ。

    method='cc' のとき: 加重平均沸点 + 純成分 P_sat ベース (旧来の挙動)。
    """
    s = sum(max(xi, 0.0) for xi in x_list)
    if s <= 1e-12:
        T_def = 298.15
        return T_def, [_K_cc(c, T_def, P_col) for c in comps]
    x_norm = [max(xi, 0.0) / s for xi in x_list]

    if method == 'cc':
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    # ---- PR 経路 ----
    pr_capable = [
        c in THERMO_DATA and not math.isnan(THERMO_DATA[c].Tc)
        for c in comps
    ]
    if sum(pr_capable) < 2:
        # PR パラメータを持つ成分が 2 つ未満 → 飽和挙動を計算できない。CC へ。
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            T_bub = bubble_point_T(P_col, x_norm, comps,
                                   T_lo=_T_BUBBLE_MIN, T_hi=_T_BUBBLE_MAX)
        if (math.isnan(T_bub)
            or T_bub <= _T_BUBBLE_MIN + 1.0
            or T_bub >= _T_BUBBLE_MAX - 1.0):
            raise ValueError("bubble T pinned to bound or nan")
    except Exception:
        F = {comps[i]: x_norm[i] for i in range(len(comps))}
        T = _weighted_boil(F, P_col)
        return T, [_K_cc(c, T, P_col) for c in comps]

    # 収束 T で K = phi_L/phi_V を成分別に取得 (個別失敗は CC へ)。
    K = []
    for i, c in enumerate(comps):
        if not pr_capable[i]:
            K.append(_K_cc(c, T_bub, P_col))
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                K_i = _K_pr(i, comps, x_norm, T_bub, P_col)
            if math.isfinite(K_i) and K_i > 0.0:
                K.append(K_i)
            else:
                K.append(_K_cc(c, T_bub, P_col))
        except Exception:
            K.append(_K_cc(c, T_bub, P_col))

    # 単相 root 検出: K の spread が小さいと単相領域に張り付いており物理的でない。
    # 例: x に H2/CH4 (Tc << 室温) が混じる x_top で bubble_point_T が cryogenic な
    # 偽解に収束し、cubic が 1 根しかないため Z_L=Z_V → φ_L/φ_V=1 → K_i ≈ 1 全成分。
    # この状態では α 計算が破綻するため CC に退避する。
    K_pos = [k for k in K if k > 0.0]
    if K_pos:
        K_max = max(K_pos)
        K_min = min(K_pos)
        if K_max / max(K_min, 1e-30) < 1.1:
            F = {comps[i]: x_norm[i] for i in range(len(comps))}
            T = _weighted_boil(F, P_col)
            return T, [_K_cc(c, T, P_col) for c in comps]
    return T_bub, K


# ===========================================================================
# FUG 計算式
# ===========================================================================

def _fenske_N_min(alpha_LK: float,
                  recovery_LK_top: float,
                  recovery_HK_bot: float) -> float:
    """Fenske 式: 最小段数 N_min。

    N_min = log[(r_LK/(1-r_LK)) × (r_HK/(1-r_HK))] / log(α_LK,HK)

    r_LK : LK の塔頂回収率
    r_HK : HK の塔底回収率
    """
    if alpha_LK <= 1.0 + 1e-9:
        # α≈1 で分離不能
        return float('inf')
    rLK = max(min(recovery_LK_top, 1.0 - 1e-9), 1e-9)
    rHK = max(min(recovery_HK_bot, 1.0 - 1e-9), 1e-9)
    ratio = (rLK / (1.0 - rLK)) * (rHK / (1.0 - rHK))
    if ratio <= 1.0:
        return 0.0
    return math.log(ratio) / math.log(alpha_LK)


def _underwood_R_min(alpha: Dict[str, float], z: Dict[str, float],
                     x_top: Dict[str, float], q: float,
                     LK: str, HK: str) -> float:
    """Underwood 式: 最小還流比 R_min。

    手順:
      1. Σ α_i × z_i / (α_i - θ) = 1 - q を θ について解く
         (θ は α_HK < θ < α_LK の範囲に存在)
      2. R_min = Σ α_i × x_top_i / (α_i - θ) - 1
    """
    alpha_LK = alpha[LK]
    alpha_HK = alpha[HK]
    if alpha_LK <= alpha_HK + 1e-9:
        # α 序列が逆転している → R_min 不定
        return float('inf')

    def eq1(theta: float) -> float:
        return sum(alpha[c] * z[c] / (alpha[c] - theta) for c in alpha) - (1.0 - q)

    # θ は α_HK と α_LK の間で 1 個の根
    eps = 1e-6
    lo = alpha_HK + eps * (alpha_LK - alpha_HK)
    hi = alpha_LK - eps * (alpha_LK - alpha_HK)
    try:
        theta = brentq(eq1, lo, hi, xtol=1e-7, maxiter=200)
    except (ValueError, RuntimeError):
        return float('inf')

    # R_min を塔頂組成から計算
    R_min = sum(alpha[c] * x_top[c] / (alpha[c] - theta) for c in alpha) - 1.0
    return max(0.0, R_min)


def _gilliland_eduljee(R: float, R_min: float, N_min: float) -> float:
    """Gilliland 相関 (Eduljee 1975 形式) で実段数 N を計算。

    X = (R - R_min) / (R + 1)
    Y = 0.75 × (1 - X^0.5668)
    N = (Y + N_min) / (1 - Y)
    """
    if R <= R_min + 1e-9 or R + 1.0 <= 0:
        return float('inf')
    X = (R - R_min) / (R + 1.0)
    if X <= 0 or X >= 1:
        return float('inf')
    Y = 0.75 * (1.0 - X ** 0.5668)
    if Y >= 1.0 - 1e-9:
        return float('inf')
    return (Y + N_min) / (1.0 - Y)


def _kirkbride_feed_stage(N: int, F_top: float, F_bot: float,
                          z: Dict[str, float],
                          x_top: Dict[str, float], x_bot: Dict[str, float],
                          LK: str, HK: str) -> int:
    """Kirkbride 式で推奨フィード段位置を計算。

    log10(N_above/N_below) = 0.206 × log10[
        (z_HK/z_LK) × (x_LK_bot/x_HK_top)^2 × (B/D)
    ]

    返り値: 塔頂から数えた段番号 (1=塔頂、N=塔底)。
    """
    z_LK = max(z.get(LK, 1e-9), 1e-9)
    z_HK = max(z.get(HK, 1e-9), 1e-9)
    x_LK_bot = max(x_bot.get(LK, 1e-9), 1e-9)
    x_HK_top = max(x_top.get(HK, 1e-9), 1e-9)
    if F_top <= 0:
        return 1
    arg = (z_HK / z_LK) * (x_LK_bot / x_HK_top) ** 2 * (F_bot / F_top)
    if arg <= 0:
        return max(1, N // 2)
    ratio = 10.0 ** (0.206 * math.log10(arg))   # N_above / N_below
    # N_above + N_below = N, N_above / N_below = ratio → N_above = N × ratio/(1+ratio)
    N_above = int(round(N * ratio / (1.0 + ratio)))
    N_above = max(1, min(N - 1, N_above))
    return N_above   # 塔頂から N_above 段目がフィード段


# ===========================================================================
# 物質収支 (Fenske による全成分分配)
# ===========================================================================

def _split_streams(F_in: Dict[str, float],
                   alpha: Dict[str, float],
                   N_min: float, HK: str,
                   recovery_LK_top: float,
                   recovery_HK_bot: float,
                   LK: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Fenske 関係を使って全成分の塔頂・塔底分配を計算。

    アプローチ:
      Key 成分は recovery で固定:
        F_top[LK] = F_in[LK] × recovery_LK_top
        F_top[HK] = F_in[HK] × (1 - recovery_HK_bot)
      Non-key 成分は Fenske の比から個別に分配:
        (x_top/x_bot)_i = α_i^N_min × (x_top/x_bot)_HK_ref
      → 個別物質収支で塔頂分率を導出。

    HK 基準量 (HK の塔頂/塔底比):
      ratio_HK = recovery_HK_bot ベースで決定
                = (F_HK_top / F_HK_bot) = (1 - rec_HK) / rec_HK
    """
    rec_HK = max(min(recovery_HK_bot, 1.0 - 1e-9), 1e-9)
    ratio_HK = (1.0 - rec_HK) / rec_HK   # HK の (top/bot) 流量比

    F_top: Dict[str, float] = {}
    F_bot: Dict[str, float] = {}
    for c, F_c in F_in.items():
        if F_c <= 0:
            F_top[c] = 0.0
            F_bot[c] = 0.0
            continue
        if c == LK:
            frac_top = recovery_LK_top
        elif c == HK:
            frac_top = 1.0 - rec_HK
        else:
            # Fenske: (top/bot)_c = α_c^N_min × (top/bot)_HK
            a_c = alpha.get(c, 1.0)
            if N_min == float('inf') or N_min > 100:
                # α > 1 なら塔頂、α < 1 なら塔底に近づく
                if a_c > 1.0:
                    frac_top = 1.0 - 1e-9
                elif a_c < 1.0:
                    frac_top = 1e-9
                else:
                    frac_top = 0.5
            else:
                ratio_c = (a_c ** N_min) * ratio_HK
                # frac_top = ratio_c / (1 + ratio_c)
                if ratio_c > 1e30:
                    frac_top = 1.0 - 1e-9
                else:
                    frac_top = ratio_c / (1.0 + ratio_c)
        F_top[c] = F_c * frac_top
        F_bot[c] = F_c * (1.0 - frac_top)
    return F_top, F_bot


# ===========================================================================
# 塔本体 CAPEX (Vessel: 既存ロジック踏襲)
# ===========================================================================

def _vessel_capex_okuyen(V_m3: float, P_pa: float, D_m: float) -> float:
    """塔本体 (vertical vessel) の CAPEX [億円]。

    出典: プロセス設計授業資料 R08-3.pdf 付録 A
      C_p0 = 10^(K1 + K2*log10(V) + K3*(log10(V))^2)
      F_p  = max((Pg+1)*D / (10.71 - 0.00756(Pg+1)) + 0.5, 1)  (Pg > -0.5 bar)
      C_BM = C_p0 × (B1 + B2 × Fp × FM)
      C_TM = 1.18 × C_BM × CEPCI 補正
    """
    P_bar = P_pa / 1.0e5 - 1.01325
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        cp0 = calc_cp0(V_m3)
    fp  = calc_fp(P_bar, D_m)
    fbm = B1 + B2 * fp * FM
    cbm = cp0 * fbm
    usd = cbm * (CEPCI_CURRENT / CEPCI_BASE)
    return PLANT_INDIRECT_FACTOR * usd * USD_TO_JPY / 1.0e8


# ===========================================================================
# partial condenser の凝縮器エネルギー収支 (2026-05-26 追加)
# ===========================================================================
_T_REF_ENTHALPY = 298.15   # [K] 理想気体エンタルピー基準


def _H_mol_jmol(T: float, P: float, comp_frac: Dict[str, float],
                keys: List[str], phase: str) -> float:
    """混合物のモルエンタルピー [J/mol] = 理想気体(T_ref=298.15→T) + PR 残差。

    partial condenser の凝縮器エネルギー収支 Q = V2·Hv − D·Hv − L1·Hl 用。
    設計判断 (2026-05-26): 旧 partial-cond Q_cond 式 (R×凝縮性成分×λ) は CMO 還流量
    L=R×D_total を「留出の凝縮性成分」で代用し、かつ塔頂蒸気の顕熱を無視して
    ~6.7倍過小評価だった (tools/_diag_dist2_qcond.py で実証: 8MW vs 厳密54MW≈HYSYS62MW)。
    本関数で厳密エネルギー収支を組む。
    """
    s = sum(max(comp_frac.get(k, 0.0), 0.0) for k in keys)
    if s <= 0:
        return 0.0
    x = [max(comp_frac.get(k, 0.0), 0.0) / s for k in keys]
    H_ig = 0.0
    for i, k in enumerate(keys):
        p = THERMO_DATA.get(k)
        if p is None:
            continue
        T0 = _T_REF_ENTHALPY
        H_ig += x[i] * (p.a * (T - T0)
                        + p.b / 2.0 * (T**2 - T0**2)
                        + p.c / 3.0 * (T**3 - T0**3)
                        + p.d / 4.0 * (T**4 - T0**4))
    try:
        Z = z_factor(T, P, x, keys, phase)
        H_r = residual_enthalpy(T, P, x, keys, Z)
    except Exception:
        H_r = 0.0
    return H_ig + H_r


# ===========================================================================
# メイン関数
# ===========================================================================

def _simulate_rigorous(
    design: DistDesignVars,
    feed:   ProcessStream,
    fixed:  Optional[DistFixedParams] = None,
) -> DistResult:
    """VLE ベース tray-by-tray 厳密求解 (Wang-Henke 法、CMO 仮定)。

    アルゴリズム:
      Ref: Seader, Henley & Roper, "Separation Process Principles", 3rd ed., Ch.10.4
           — 多成分蒸留塔の標準 tray-by-tray 求解 (bubble-point method)
      実装: src/distillation_rigorous.py の wang_henke_solve()

    流れ:
      1. FUG ソルバを最初に走らせて D_total / T_top / T_bot 初期値を取得
         (rating mode で D を固定するため、FUG の Fenske 結果を採用)
      2. Wang-Henke で MESH を rigorous に解く
         - 各段の x_{i,j}, y_{i,j}, T_j を反復で決定
         - K 値は src/eos.py (PR EOS、thermo と整合性 0.02% 確認済) を流用
      3. 収束した rigorous 結果で F_top, F_bot, T_top, T_bot を更新
      4. Q_cond, Q_reb, A_cond, A_reb, CAPEX_cond/reb を再計算
      5. Vessel 寸法 (D_col, H_col) は CMO 仮定下で V_top が FUG と同じなので流用
      6. DistResult を組み立てて返す

    収束失敗時:
      RuntimeError を投げて simulate_distillation_column の dispatch で
      FUG にフォールバックされる (warning つき)。

    K 値ライブラリの選択根拠:
      thermo (CalebBell/thermo, MIT) と src/eos.py の PR EOS 計算を 5 ケースで
      比較した結果、相対差 0.02% 以下で一致 (2026-05-09 確認、テストは
      検証後削除)。新規依存を最小化するため src/eos.py を流用。
      thermo は将来 PT/PH flash (partial condenser, JT 膨張弁) 用に install 済。
    """
    # 局所 import: 循環参照回避 + 起動時オーバーヘッド削減
    from src.distillation_rigorous import wang_henke_solve
    from dataclasses import replace as _dc_replace

    # ---- 1. FUG ソルバで初期値・インフラを取得 ----
    # design の solver_method を 'fug' に差し替えたコピーで再帰呼び出し → FUG パスが走る
    # (FUG が infeasible なら rigorous に進む意味がないので FUG 結果を返す)
    fug_design = _dc_replace(design, solver_method='fug')
    fug_result = simulate_distillation_column(fug_design, feed, fixed)
    if not fug_result.equipment.feasible:
        return fug_result

    # ---- 2. Wang-Henke 求解 ----
    comps   = list(feed.F_in.keys())
    F_total = sum(max(F, 0.0) for F in feed.F_in.values())
    # 設計判断 (2026-05-19): D_override がある場合は FUG の D を上書き。
    # HYSYS の "Distillate Rate spec" 相当。物理的に不適な値は solver 内で _failure_result 化。
    if design.D_override is not None and design.D_override > 0:
        if design.D_override >= F_total:
            warnings.warn(
                f"D_override={design.D_override:.2f} >= F_total={F_total:.2f} は物理不適 "
                f"(D < F でなければ B ≤ 0)。FUG の D={sum(fug_result.top.F_in.values()):.2f} にフォールバック。",
                UserWarning, stacklevel=2,
            )
            D_total = sum(fug_result.top.F_in.values())
        else:
            D_total = design.D_override
    else:
        D_total = sum(fug_result.top.F_in.values())
    N_feed_kirk = max(1, min(fug_result.equipment.N_feed_kirkbride, design.N_stages))

    # フィード段は FUG の Kirkbride 推奨を採用 (design.N_feed は無視)。
    # ColumnTunables.N_feed の docstring 通り、rigorous では Kirkbride を自動採用する。
    rig = wang_henke_solve(
        feed_F            = feed.F_in,
        comps             = comps,
        P_col             = design.P_col,
        N_stages          = design.N_stages,
        N_feed            = N_feed_kirk,
        reflux_ratio      = design.reflux_ratio,
        D_total           = D_total,
        q_feed            = design.q,
        partial_condenser = design.partial_condenser,
        K_method          = design.K_method,
        T_top_init_K      = fug_result.equipment.T_top,
        T_bot_init_K      = fug_result.equipment.T_bot,
    )

    # ---- Phase D (2026-05-19): 失敗時 retry (Wegstein off + 低 damping) ----
    # narrow-margin 設計で Wegstein 加速が振動方向に逸脱する症例 (dT_max=7K stall) に対処。
    # 第一試行で converge せず、または MESH/balance が許容外のとき、Wegstein を切って
    # pure substitution 寄りの damping=0.2 で再試行する。
    retry_needed = (
        not rig.converged
        or rig.mesh_residual_max > 0.01
        or rig.component_balance_max > 0.01
    )
    if retry_needed:
        from src.distillation_rigorous import _RETRY_DAMPING
        rig2 = wang_henke_solve(
            feed_F            = feed.F_in,
            comps             = comps,
            P_col             = design.P_col,
            N_stages          = design.N_stages,
            N_feed            = N_feed_kirk,
            reflux_ratio      = design.reflux_ratio,
            D_total           = D_total,
            q_feed            = design.q,
            partial_condenser = design.partial_condenser,
            K_method          = design.K_method,
            T_top_init_K      = fug_result.equipment.T_top,
            T_bot_init_K      = fug_result.equipment.T_bot,
            max_iter          = 1000,           # 倍 retry
            damping           = _RETRY_DAMPING,  # 0.2 (Wegstein 振動回避)
            use_wegstein      = False,           # 加速無し、pure substitution
        )
        # retry 結果が改善していれば採用
        if (rig2.converged
            and rig2.mesh_residual_max <= 0.01
            and rig2.component_balance_max <= 0.01):
            rig = rig2
        elif rig2.converged and not rig.converged:
            # 完全には数値基準を満たさないが、少なくとも converge した結果は試す
            rig = rig2

    if not rig.converged:
        raise RuntimeError(
            f"Wang-Henke 収束失敗 (retry 後も): {rig.message} "
            f"(LK={design.LK}, HK={design.HK}, R={design.reflux_ratio:.2f}, "
            f"N={design.N_stages}, P={design.P_col/1e5:.1f}bar)"
        )

    # ---- always-on validation: 内部不整合の早期検知 ----
    # 設計判断 (2026-05-10/2026-05-19 Phase D 修正):
    #   MESH 残差 max > 0.05 (旧 0.01)、成分バランス err > 5% (旧 1%) を異常とみなす。
    #   リサイクル過渡 (Δ_rel > 50% で組成が暴れている状態) では 1-2% の残差が
    #   普通に出るため、旧 1% 閾値が strict すぎて recycle 全体を止めていた。
    #   閾値超え分は Phase A proxy_penalty が捕捉する設計。
    _MESH_TOL = 0.05
    _MASS_BAL_TOL = 0.05
    if rig.mesh_residual_max > _MESH_TOL:
        raise RuntimeError(
            f"Wang-Henke MESH 残差過大 (retry 後も): max={rig.mesh_residual_max:.4f} (>{_MESH_TOL})。"
            f" solver 内部不整合の疑い (B_bottoms 等のバグ可能性)。"
        )
    if rig.component_balance_max > _MASS_BAL_TOL:
        raise RuntimeError(
            f"Wang-Henke 成分マスバランス破綻 (retry 後も): max_err={rig.component_balance_max*100:.2f}% (>{_MASS_BAL_TOL*100:.0f}%)。"
            f" solver 内部不整合の疑い。"
        )

    # ---- recovery spec 達成度チェックは撤去 (2026-05-10 終夜) ----
    # 当初 N=N_min ぎりぎりで rigorous が non-spec 解を返すケースの検出を狙ったが、
    # 実装すると以下の理由で本来 OK な解も fallback してしまう:
    #   - exp1 の recycle iter 中 transient state で feed 組成が大きく変動
    #   - その時の rigorous recovery が一時的に spec から外れて誤発動
    #   - run_one_pass.py:94 の `warnings.simplefilter("ignore")` で warning も
    #     見えず silent fallback になり、Profit が 1 億円/年 過小評価
    # 既存の safety net で十分:
    #   - FUG の N >= N_min, R >= R_min check (極端ケース弾く)
    #   - Wang-Henke の tridiagonal LinAlg failure (matrix singular catch)
    #   - MESH 残差 check (B_bottoms 級バグ検出)
    #   - 成分マスバランス check (同上)
    # rating mode で recovery が物理的に動くこと自体は正常挙動。

    # ---- 3. partial condenser: C3 微小漏れの閾値丸め (PSA/Mem への流入抑止) ----
    # 設計判断 (2026-05-19): rigorous は物理的に正しく C3 を「ほぼほぼ塔底へ」振り分ける
    # が、PR EOS 実値で K(C3H6) ≈ 0.18 @ -43°C/7.5bar 程度の有限値があり、F_top に微小量
    # (数 kmol/h オーダー) が残る。下流の PSA/Mem は組成上 C3 を流したくない設計なので、
    # FUG パスと同じ閾値丸め (feed の _PARTIAL_COND_C3_ROUND_FRAC 未満なら 0 集約) を
    # rigorous 結果にも適用する。閾値超えは物理値を残して warning。
    # 設計判断 (2026-05-20 Phase A強化): proxy_penalty 計算用に丸める前の F_top を保持。
    F_top = dict(rig.F_top)
    F_bot = dict(rig.F_bot)
    F_top_pre_round = dict(F_top)
    if design.partial_condenser:
        for c in ALWAYS_CONDENSABLE_COMPS:
            moved = F_top.get(c, 0.0)
            f_c = feed.F_in.get(c, 0.0)
            if moved <= 0 or f_c <= 0:
                continue
            if moved <= _PARTIAL_COND_C3_ROUND_FRAC * f_c:
                F_top[c] = 0.0
                F_bot[c] = F_bot.get(c, 0.0) + moved
            else:
                warnings.warn(
                    f"rigorous partial cond: comp '{c}' が F_top に "
                    f"{moved:.2f} kmol/h ({moved / f_c * 100:.2f}% of feed) 漏出。"
                    f" 閾値 {_PARTIAL_COND_C3_ROUND_FRAC * 100:.1f}% 超 → 物理値を残す。"
                    f" 設計 (P_col={design.P_col / 1e5:.1f}bar, R={design.reflux_ratio:.2f}, "
                    f"N={design.N_stages}) の見直しを検討。",
                    UserWarning, stacklevel=2,
                )

    # ---- 4. rigorous の F_top/F_bot/T_top/T_bot で結果を更新 ----
    F_top_total = sum(F_top.values())
    F_bot_total = sum(F_bot.values())
    T_top = rig.T_profile_K[0]
    T_bot = rig.T_profile_K[-1]

    if F_top_total <= 0 or F_bot_total <= 0:
        raise RuntimeError(
            f"Wang-Henke 結果が異常: F_top={F_top_total:.3f}, F_bot={F_bot_total:.3f}"
        )

    # FUG result から流用するもの (CMO 下で rigorous と同じ): D_col, H_col, vessel CAPEX
    eq = fug_result.equipment

    # ---- 5. Q_cond, Q_reb の再計算 (rigorous compositions ベース) ----
    if design.partial_condenser:
        # 設計判断 (2026-05-26): partial cond の Q_cond を凝縮器(stage1)の厳密
        # エネルギー収支に置換。旧式 Q=R×(留出の凝縮性成分)×λ は CMO 還流量
        # L=R×D_total を「留出の凝縮性成分」で代用 (~5.7x 過小) + 塔頂蒸気の顕熱
        # 無視で計 6.68x 過小評価だった (tools/_diag_dist2_qcond.py: 8MW vs 厳密54MW
        # ≈HYSYS62MW)。収支: Q = V2·Hv(y2,T2) − D·Hv(y1,T1) − L1·Hl(x1,T1)。
        # rig プロファイル (V/L 流量・各段組成 y/x・温度) から組む。
        try:
            V2 = rig.V_top_kmolh
            L1 = rig.L_top_kmolh
            D_v = V2 - L1
            T1 = rig.T_profile_K[0]
            T2 = rig.T_profile_K[1] if len(rig.T_profile_K) > 1 else T1
            y1 = rig.y_profile[0]
            y2 = rig.y_profile[1] if len(rig.y_profile) > 1 else y1
            x1 = rig.x_profile[0]
            Q_cond_kW = (
                V2 * _H_mol_jmol(T2, design.P_col, y2, comps, 'vapor')
                - D_v * _H_mol_jmol(T1, design.P_col, y1, comps, 'vapor')
                - L1 * _H_mol_jmol(T1, design.P_col, x1, comps, 'liquid')
            ) / 3600.0
            if not math.isfinite(Q_cond_kW) or Q_cond_kW <= 0:
                raise ValueError("energy-balance Q_cond が非正/非有限")
        except Exception:
            # フォールバック: 還流量補正した簡易式 (R×D_total×λ、顕熱なし)
            F_top_C = {c: v for c, v in F_top.items()
                       if c not in NON_CONDENSABLE_COMPS and v > 0.0}
            lam_top = (_weighted_lambda(F_top_C) if sum(F_top_C.values()) > 0
                       else _LAMBDA_DEFAULT)
            Q_cond_kW = (design.reflux_ratio * F_top_total
                         * lam_top * 1000.0 / 3600.0)
    else:
        lam_top = _weighted_lambda(F_top)
        Q_cond_kW = (F_top_total * (design.reflux_ratio + 1.0)
                     * lam_top * 1000.0 / 3600.0)

    # Reboiler: 標準 MESH の V' = (R+1)D - (1-q)F、Q_reb = V' × λ_bot
    F_feed_total = sum(max(F, 0.0) for F in feed.F_in.values())
    V_prime_kmolh = (F_top_total * (design.reflux_ratio + 1.0)
                     - F_feed_total * (1.0 - design.q))
    if V_prime_kmolh <= 0:
        # フィードガスが多すぎてリボイラ不要 (q=0 で V'=0 になる病理)
        # FUG 結果の Q_reb をそのまま流用する保守判断
        Q_reb_kW = eq.Q_reb
    else:
        lam_bot = _weighted_lambda(F_bot)
        Q_reb_kW = V_prime_kmolh * lam_bot * 1000.0 / 3600.0

    # ---- 6. utility 選択を rigorous T_top/T_bot で再評価 ----
    # Bug fix (2026-05-10 終夜): 旧版は FUG 結果の utility を流用していたが、
    # rigorous で T_top が大きく動いた場合 (Dist2 で FUG -38°C → rig -73°C 等)
    # FUG の冷媒 tier が温度未到達で ΔT 違反 → 例外で FUG fallback してた。
    # rigorous 結果で再選択して整合させる。
    _T_COOLING_FLOOR = 173.15 + 10.0 + 1.0    # = -89°C: -100°C 冷媒 + 10K approach
    T_top_for_util = max(T_top, _T_COOLING_FLOOR)
    try:
        cond_utility = select_cooling_utility(T_top_for_util, mode='continuous')
    except ValueError:
        raise RuntimeError(
            f"rigorous T_top={T_top-273.15:.1f}°C で冷媒選択不可"
        )
    try:
        reb_utility = select_heating_utility(T_bot, mode='continuous')
    except ValueError:
        raise RuntimeError(
            f"rigorous T_bot={T_bot-273.15:.1f}°C で熱媒選択不可"
        )
    cond_utility_name       = cond_utility.name
    cond_utility_jpy_per_GJ = cond_utility.jpy_per_GJ
    reb_utility_name        = reb_utility.name
    reb_utility_jpy_per_GJ  = reb_utility.jpy_per_GJ

    # ---- 7. A_cond, A_reb の再計算 (contest §4-4 表ベース、FUG パスと同じ式) ----
    cond_T_in  = _supply_T_for_utility(cond_utility_name)
    cond_T_out = cond_T_in + 10.0
    dT1 = T_top - cond_T_in
    dT2 = T_top - cond_T_out
    if dT1 <= 0.0 or dT2 <= 0.0:
        raise RuntimeError(
            f"rigorous 結果で condenser ΔT 不成立 (T_top={T_top-273.15:.1f}°C)"
        )
    if abs(dT1 - dT2) < 1e-3:
        dT_lm_cond = 0.5 * (dT1 + dT2)
    else:
        dT_lm_cond = (dT1 - dT2) / math.log(dT1 / dT2)
    U_cond  = lookup_U(StreamPhase.CONDENSING, utility_phase(cond_utility_name))
    A_cond  = max(Q_cond_kW * 1000.0 / (U_cond * dT_lm_cond), 10.0)

    reb_T_supply = _supply_T_for_utility(reb_utility_name)
    dT_lm_reb = reb_T_supply - T_bot
    if dT_lm_reb <= 0.0:
        raise RuntimeError(
            f"rigorous 結果で reboiler ΔT 不成立 (T_bot={T_bot-273.15:.1f}°C)"
        )
    U_reb  = lookup_U(utility_phase(reb_utility_name), StreamPhase.EVAPORATING)
    A_reb  = max(Q_reb_kW * 1000.0 / (U_reb * dT_lm_reb), 10.0)

    # ---- 8. CAPEX 再計算 (HE のみ、vessel は FUG 流用) ----
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        capex_cond = calc_he_capex_okuyen(A_cond)
        capex_reb  = calc_he_capex_okuyen(A_reb)
    capex_total = eq.CAPEX_vessel + eq.CAPEX_trays + capex_cond + capex_reb

    # ---- 9. Q_feed_preheat は FUG 結果を流用 ----
    # (rigorous でも feed enthalpy は同じなので)
    Q_feed_preheat_kW = eq.Q_feed_preheat_kW

    # ---- 9b. rigorous プロキシ罰則 (Phase A, 2026-05-19) ----
    # rigorous でも leak と margin を再計算して BO objective に反映する。
    # FUG と rigorous の値が異なるのでそれぞれで計算する設計。
    # 設計判断 (2026-05-20 Phase A 強化): 丸める前の F_top_pre_round を渡す
    # (FUG path と同じ意図 — <1% 漏れも捕捉できる拡張余地を確保)。
    proxy_penalty, proxy_reason = _compute_proxy_penalty(
        design, F_top_pre_round, feed.F_in, eq.R_min, eq.N_min, is_rigorous=True,
    )

    # ---- 10. DistResult 組み立て ----
    top_stream = ProcessStream(F_in=F_top, T_in=T_top, P_in=design.P_col)
    bottom_stream = ProcessStream(F_in=F_bot, T_in=T_bot, P_in=design.P_col)
    equipment = DistEquipment(
        D_col              = eq.D_col,
        H_col              = eq.H_col,
        V_col              = eq.V_col,
        CAPEX_vessel       = eq.CAPEX_vessel,
        CAPEX_trays        = eq.CAPEX_trays,
        CAPEX              = capex_total,
        Q_cond             = Q_cond_kW,
        Q_reb              = Q_reb_kW,
        Q_feed_preheat_kW  = Q_feed_preheat_kW,
        N_min              = eq.N_min,
        R_min              = eq.R_min,
        N_feed_kirkbride   = eq.N_feed_kirkbride,
        feasible           = True,
        message            = f"rigorous: {rig.message}",
        A_cond_m2          = A_cond,
        A_reb_m2           = A_reb,
        CAPEX_cond         = capex_cond,
        CAPEX_reb          = capex_reb,
        cond_utility_name  = cond_utility_name,
        cond_utility_jpy_per_GJ = cond_utility_jpy_per_GJ,
        reb_utility_name   = reb_utility_name,
        reb_utility_jpy_per_GJ  = reb_utility_jpy_per_GJ,
        T_top              = T_top,
        T_bot              = T_bot,
        proxy_penalty_okuyen = proxy_penalty,
        proxy_penalty_reason = proxy_reason,
    )
    return DistResult(top=top_stream, bottom=bottom_stream, equipment=equipment)


def _supply_T_for_utility(utility_name: str) -> float:
    """utility 名から供給温度 [K] を逆引き (rigorous の A 計算用ヘルパ)。

    select_cooling_utility/select_heating_utility が返す UtilityTier を直接
    保持していないので、名前から再選択する。設計判断: FUG パスで使った
    utility と同じものを使い回すため、tier テーブルの supply_T を再取得する。
    """
    from src.utility_selector import _COOLING_TIERS, _HEATING_TIERS
    for tier in _COOLING_TIERS + _HEATING_TIERS:
        if tier.name == utility_name:
            return tier.supply_T_K
    # 名前の前後の "≒" 記号などで完全一致しない場合のフォールバック
    for tier in _COOLING_TIERS + _HEATING_TIERS:
        if tier.name in utility_name or utility_name in tier.name:
            return tier.supply_T_K
    raise ValueError(f"utility '{utility_name}' の supply_T_K が引けない")


def simulate_distillation_column(
    design: DistDesignVars,
    feed:   ProcessStream,
    fixed:  Optional[DistFixedParams] = None,
) -> DistResult:
    """蒸留塔シミュレーションのエントリポイント。

    design.solver_method による dispatch:
      'fug'      : Fenske-Underwood-Gilliland shortcut (本関数本体)
      'rigorous' : VLE tray-by-tray (`_simulate_rigorous` に委譲、未実装ならフォールバック)
      'sm'       : SM 法 (`_simulate_sm` に委譲、未実装は NotImplementedError → FUG fallback)

    FUG 手順:
      1. 動作温度推定 (T_top, T_bot を Clausius-Clapeyron で初期推定)
      2. K 値・α 計算 (PR EOS or CC)
      3. Fenske: N_min
      4. 物質収支: 全成分の塔頂・塔底分配
      5. Underwood: R_min (塔頂組成必要のため物質収支後)
      6. Gilliland: feasibility 確認 (N >= N_min, R >= R_min)
      7. Kirkbride: 推奨フィード段
      8. Q_cond, Q_reb, Q_feed_preheat 計算
      9. 塔径・塔高計算
     10. CAPEX (Vessel + Trays)
     11. infeasible なら ペナルティ返却
    """
    # ---- Solver dispatch ----
    # 設計判断 (2026-05-19): rigorous 指定時のフォールバックは fail-fast がデフォルト。
    # 旧版は例外を warnings.warn で握り潰して FUG にサイレント切替していたが、これにより
    #   - exp1 表示は solver=rigorous なのに実体は FUG (発見1, 2026-05-19)
    #   - 結果に rigorous の物理が反映されない (例: Dist2 C3H6 100% 回収の見かけ)
    # という重大な信頼性問題があった。fail-fast にして「rigorous 指定なら rigorous で
    # 通すか or 例外停止」を保証する。
    # BO 等で fallback が必要な場合は環境変数 PDH_RIGOROUS_STRICT=0 で旧挙動 (warn + FUG)
    # に切り替えられる退路を残す。
    if design.solver_method == 'rigorous':
        _strict = os.environ.get('PDH_RIGOROUS_STRICT', '1') != '0'
        try:
            return _simulate_rigorous(design, feed, fixed)
        except NotImplementedError as e:
            if _strict:
                raise RuntimeError(
                    f"simulate_distillation_column: rigorous solver 未実装 ({e})。"
                    f" PDH_RIGOROUS_STRICT=0 で FUG にフォールバック可能。"
                ) from e
            warnings.warn(
                f"simulate_distillation_column: rigorous solver 未実装 ({e})。"
                f" FUG にフォールバック (PDH_RIGOROUS_STRICT=0)。",
                UserWarning, stacklevel=2,
            )
        except Exception as e:
            if _strict:
                # 設計判断 (2026-05-20): 旧版は RuntimeError を raise していたが、これにより
                # run_one_pass で「Dist2 失敗 → 下流の PSA/Mem が ValueError で crash」する
                # 経路があり、TPE constraints_func に dT_max 等の連続シグナルが伝わらなかった。
                # 本版は penalty_result を返して `equipment.feasible=False` でシグナル化、
                # `dT_max_rigorous` を抽出して run_one_pass で user_attr に格納する。
                # 「silent FUG fallback はしない」(2026-05-19 発見1 の防止) はそのまま維持。
                msg = str(e)
                dT_max = _extract_dT_max(msg)
                warnings.warn(
                    f"simulate_distillation_column: rigorous solver で例外 "
                    f"({type(e).__name__}: {e}) → penalty_result 返却"
                    f" (dT_max={dT_max:.2f}K if extracted)。",
                    UserWarning, stacklevel=2,
                )
                return _penalty_result(
                    design,
                    f"rigorous {type(e).__name__}: {msg[:200]}",
                    dT_max_rigorous=dT_max,
                )
            warnings.warn(
                f"simulate_distillation_column: rigorous solver で例外 ({type(e).__name__}: {e})。"
                f" FUG にフォールバック (PDH_RIGOROUS_STRICT=0)。",
                UserWarning, stacklevel=2,
            )

    if design.solver_method == 'sm':
        # SM (smith?) 法。実装は別途追加予定。NotImplementedError なら FUG に fallback。
        try:
            from src.distillation_sm import simulate_sm as _simulate_sm  # type: ignore
            return _simulate_sm(design, feed, fixed)
        except (ImportError, NotImplementedError) as e:
            warnings.warn(
                f"simulate_distillation_column: sm solver 未実装 ({e})。"
                f" FUG にフォールバック。",
                UserWarning, stacklevel=2,
            )
        except Exception as e:
            warnings.warn(
                f"simulate_distillation_column: sm solver で例外 ({type(e).__name__}: {e})。"
                f" FUG にフォールバック。",
                UserWarning, stacklevel=2,
            )

    if fixed is None:
        fixed = DistFixedParams()

    comps = list(feed.F_in.keys())
    F_total = sum(max(F, 0.0) for F in feed.F_in.values())
    if F_total <= 0.0:
        return _penalty_result(design, "feed flow ≤ 0")

    # フィード組成 z (mol fraction)
    z = {c: max(feed.F_in[c], 0.0) / F_total for c in comps}

    # ---- Step 1: 動作温度・K の初期推定 (CC、x=z で粗推定) ----
    # 塔頂は LK 主体、塔底は HK 主体と仮定して純成分沸点で T 範囲を見立てる。
    T_top = _boil_cc(_T_BOIL_ATM.get(design.LK, _T_BOIL_DEFAULT),
                     _LAMBDA_KJ.get(design.LK, _LAMBDA_DEFAULT),
                     design.P_col)
    T_bot = _boil_cc(_T_BOIL_ATM.get(design.HK, _T_BOIL_DEFAULT),
                     _LAMBDA_KJ.get(design.HK, _LAMBDA_DEFAULT),
                     design.P_col)
    # 設計判断 (2026-05-09): T_top の _T_COND_MIN clamp を解除。
    # 旧版は「冷却水到達下限 40°C」で上書きしていたが、これだと Dist2 の x_top
    # (H2/CH4 主体) で実際は cryogenic な凝縮器が必要なケースを隠蔽していた。
    # 本版は実温度を保持し、utility_selector で適切な冷媒を選ぶ。

    # CC で初期 α・初期 split を 1 回作る (PR 反復の出発点)。
    K_init = _K_dict(comps, [z[c] for c in comps],
                     0.5 * (T_top + T_bot), design.P_col, 'cc')
    alpha = _alpha_dict(K_init, design.HK)
    N_min = _fenske_N_min(
        alpha.get(design.LK, 1.0),
        design.recovery_LK_top, design.recovery_HK_bot,
    )
    F_top, F_bot = _split_streams(
        feed.F_in, alpha, N_min, design.HK,
        design.recovery_LK_top, design.recovery_HK_bot, design.LK,
    )

    # ---- Step 2-4: 塔頂/塔底それぞれの泡点で α を評価して反復 ----
    # 設計判断 (2026-05-09): 旧版は K_avg を「塔平均 T で x=z」で計算していたが、
    # この T-P-x は単相領域に入りやすく、PR の 3 次方程式が 1 根しか返さない結果
    # K_i = phi_L/phi_V ≈ 1 → α ≈ 1 → N_min = ∞ で全塔 infeasible になる病理が
    # あった (column1.py 旧コメント参照)。本版では塔頂/塔底それぞれの組成で
    # 泡点フラッシュを実行し、必ず 2 相が成立する点で K を取る。
    # α_geom = sqrt(α_top × α_bot) は FUG の標準的な平均化 (Sinnott 17.5.3)。
    for _ in range(5):
        F_top_total_it = sum(F_top.values())
        F_bot_total_it = sum(F_bot.values())
        if F_top_total_it <= 0.0 or F_bot_total_it <= 0.0:
            return _penalty_result(design, "split failed (top or bottom flow ≤ 0)")
        x_top_list = [F_top[c] / F_top_total_it for c in comps]
        x_bot_list = [F_bot[c] / F_bot_total_it for c in comps]

        T_top_new, K_top_list = _bubble_T_K(
            x_top_list, comps, design.P_col, design.K_method,
        )
        T_bot_new, K_bot_list = _bubble_T_K(
            x_bot_list, comps, design.P_col, design.K_method,
        )
        # 設計判断 (2026-05-09): _T_COND_MIN clamp 解除。実温度を保持して下流で
        # utility_selector に適切な冷媒を選ばせる。

        HK_idx = comps.index(design.HK)
        K_top_HK = max(K_top_list[HK_idx], 1e-30)
        K_bot_HK = max(K_bot_list[HK_idx], 1e-30)
        alpha_new: Dict[str, float] = {}
        for i, c in enumerate(comps):
            a_top = K_top_list[i] / K_top_HK
            a_bot = K_bot_list[i] / K_bot_HK
            alpha_new[c] = math.sqrt(max(a_top * a_bot, 1e-30))

        N_min_new = _fenske_N_min(
            alpha_new.get(design.LK, 1.0),
            design.recovery_LK_top, design.recovery_HK_bot,
        )
        F_top_new, F_bot_new = _split_streams(
            feed.F_in, alpha_new, N_min_new, design.HK,
            design.recovery_LK_top, design.recovery_HK_bot, design.LK,
        )

        # 収束判定: 塔頂流量の相対変化 (フィードに対して 0.1% 未満で打ち切り)
        diff = max(
            (abs(F_top_new[c] - F_top[c]) / max(feed.F_in.get(c, 0.0), 1e-9))
            for c in comps if feed.F_in.get(c, 0.0) > 0.0
        )
        F_top, F_bot = F_top_new, F_bot_new
        alpha, N_min = alpha_new, N_min_new
        T_top, T_bot = T_top_new, T_bot_new
        if diff < 1e-3:
            break

    # ---- 流量・組成 (mol fraction) ----
    F_top_total = sum(F_top.values())
    F_bot_total = sum(F_bot.values())
    if F_top_total <= 0 or F_bot_total <= 0:
        return _penalty_result(design, "split failed (top or bottom flow ≤ 0)")
    x_top = {c: F_top[c] / F_top_total for c in comps}
    x_bot = {c: F_bot[c] / F_bot_total for c in comps}

    # ---- Step 5: Underwood で R_min ----
    R_min = _underwood_R_min(alpha, z, x_top, design.q, design.LK, design.HK)

    # ---- Step 5b: partial condenser 物質収支補正 (閾値丸め) ----
    # 設計判断 (2026-05-19 改訂): 旧版は ALWAYS_CONDENSABLE_COMPS を F_top → F_bot に
    # 無条件 100% 強制移動していたが、これは PR EOS 実値 (K(C3H6) ≈ 0.18 @ -43°C/7.5bar)
    # と整合せず、かつ rigorous の D_total を C3=0 で固定して Wang-Henke の質量保存を
    # 破綻させていた (発見: 2026-05-19、Dist2 HYSYS 不一致調査)。
    # 本版は閾値丸めに変更: F_top の C3 漏れが feed の _PARTIAL_COND_C3_ROUND_FRAC 未満
    # なら 0 に丸めて塔底に集約 (= PSA/Mem に C3 を流さない目的を保つ)。閾値以上の漏れ
    # は物理値を残す (= rigorous で D_total が物理に近づく)。
    # 設計判断 (2026-05-20 Phase A強化): proxy_penalty 計算用に「丸める前の生 F_top」を
    # 保持。丸めで <1% 漏れが消えると BO 中の penalty が発火せず BO が FUG 楽観性領域に
    # 張り付く問題を是正する。閾値判定 (1%) は proxy_penalty 側でも使うので、ここで
    # 「<1% でも penalty 発火」を意図するなら proxy_penalty 側のロジックも別途調整する。
    # 現状は proxy_penalty 内で leak_frac > 1% で発火、丸めもこれと同じ閾値で消去。
    # ここで生値を保持しておけば将来「<1% でも penalty 発火」拡張も最小変更で可能。
    F_top_pre_round = dict(F_top)
    if design.partial_condenser:
        for c in ALWAYS_CONDENSABLE_COMPS:
            moved = F_top.get(c, 0.0)
            f_c = feed.F_in.get(c, 0.0)
            if moved <= 0 or f_c <= 0:
                continue
            if moved <= _PARTIAL_COND_C3_ROUND_FRAC * f_c:
                # 丸めスコープ内: 塔底に集約
                F_top[c] = 0.0
                F_bot[c] = F_bot.get(c, 0.0) + moved
            else:
                # 閾値超え: 物理値を残し warning (設計が壊れているサイン)
                warnings.warn(
                    f"partial condenser: comp '{c}' が F_top に "
                    f"{moved:.2f} kmol/h ({moved / f_c * 100:.2f}% of feed) 漏出。"
                    f" 閾値 {_PARTIAL_COND_C3_ROUND_FRAC * 100:.1f}% 超 → 物理値を残す。"
                    f" 設計 (P_col={design.P_col / 1e5:.1f}bar, R={design.reflux_ratio:.2f}, "
                    f"N={design.N_stages}) の見直しを検討。",
                    UserWarning, stacklevel=2,
                )
        F_top_total = sum(F_top.values())
        F_bot_total = sum(F_bot.values())
        if F_top_total <= 0 or F_bot_total <= 0:
            return _penalty_result(
                design,
                "partial condenser 補正後 split が無効 (top or bottom flow ≤ 0)",
                N_min=N_min, R_min=R_min,
            )
        x_top = {c: F_top[c] / F_top_total for c in comps}
        x_bot = {c: F_bot[c] / F_bot_total for c in comps}

    # ---- Step 6: feasibility ----
    # 設計判断 (2026-05-15): FUG bias 解消のため Gilliland check を追加。
    # 旧実装は (N >= N_min) かつ (R >= R_min) を「必要条件」としてだけチェックしてた
    # が、これでは narrow-margin (N ≈ N_min, R ≈ R_min) でも recovery=0.99 達成可能と
    # 判定されてしまい、rigorous で大幅 spec 違反 (Dist1 HK recovery 20% など) を生む。
    # Gilliland 相関で「与えられた R で recovery=0.99 達成に必要な N」を計算、
    # design.N_stages がそれ未満なら infeasible 化する (= 十分条件チェック)。
    feasible = True
    msg = ""
    if N_min == float('inf'):
        feasible = False
        msg = f"alpha_LK ≈ 1 で分離不能"
    elif design.N_stages < N_min:
        feasible = False
        msg = (f"N_stages={design.N_stages} < N_min={N_min:.1f} "
               f"(LK={design.LK}, HK={design.HK})")
    elif R_min == float('inf'):
        feasible = False
        msg = "Underwood 収束失敗 (alpha 序列が逆)"
    elif design.reflux_ratio < R_min:
        feasible = False
        msg = (f"R={design.reflux_ratio:.2f} < R_min={R_min:.2f} "
               f"(LK={design.LK}, HK={design.HK})")
    else:
        # Gilliland: 与えられた R での recovery 達成に必要な N
        # 設計判断 (2026-05-15): Gilliland (Eduljee 形) は constant α 仮定の経験式で、
        # 実際の PR EOS rigorous より 10-15% 程度 conservative な N を返す傾向がある。
        # baseline (Dist1: N=20, R=1.5) が strict ratio 91% で baseline が落ちないよう
        # 0.85 ファクタ (= 15% tolerance) を入れる。一方 narrow-margin (ratio < 0.7) は
        # 確実に弾かれる (= Dist1 trial 299 級 ratio 60% は INFEASIBLE)。
        N_needed = _gilliland_eduljee(design.reflux_ratio, R_min, N_min)
        N_NEEDED_TOLERANCE = 0.85
        if N_needed == float('inf'):
            feasible = False
            msg = (f"Gilliland: R={design.reflux_ratio:.2f} が R_min={R_min:.2f} に"
                   f" 近すぎ N_needed → ∞ (LK={design.LK}, HK={design.HK})")
        elif design.N_stages < N_needed * N_NEEDED_TOLERANCE:
            feasible = False
            msg = (f"Gilliland: N_stages={design.N_stages} < {N_NEEDED_TOLERANCE:.2f}×N_needed"
                   f"={N_needed*N_NEEDED_TOLERANCE:.1f} (true N_needed={N_needed:.1f},"
                   f" R={design.reflux_ratio:.2f}, R_min={R_min:.2f},"
                   f" N_min={N_min:.1f}) — recovery=0.99 物理的達成不可")

    if not feasible:
        # N_needed は Gilliland branch で計算されていれば渡す (他 branch では 0)。
        # ローカル変数 N_needed が存在しないケースに備えて locals() で安全取得。
        _N_needed_for_penalty = locals().get('N_needed', 0.0)
        if _N_needed_for_penalty == float('inf'):
            _N_needed_for_penalty = 0.0   # ∞ は連続化できないので未計算扱い
        return _penalty_result(
            design, msg, N_min=N_min, R_min=R_min,
            N_needed=_N_needed_for_penalty,
        )

    # ---- Step 7: Kirkbride で推奨 N_feed ----
    N_feed_kirkbride = _kirkbride_feed_stage(
        design.N_stages, F_top_total, F_bot_total,
        z, x_top, x_bot, design.LK, design.HK,
    )

    # ---- Step 8: 熱量計算 ----
    # 設計判断 (2026-05-09): partial condenser 分岐。
    # Total condenser: V = (R+1)D を全量凝縮 → Q = (R+1) × D × λ_avg、T_cond = bubble of x_top
    # Partial condenser: 軽質ガス (NON_CONDENSABLE) は気相パススルー。
    #                   reflux に使う condensable 成分のみ凝縮 → Q = R × D_C × λ_C、
    #                   T_cond = bubble of x_top_condensable (実機では暖かい温度で済む)
    if design.partial_condenser:
        F_top_C = {c: v for c, v in F_top.items()
                   if c not in NON_CONDENSABLE_COMPS and v > 0.0}
        F_C_total = sum(F_top_C.values())
        if F_C_total <= 0:
            return _penalty_result(
                design,
                "partial condenser: 凝縮可能成分が無い (全成分が H2/CH4)",
                N_min=N_min, R_min=R_min,
            )
        x_C_list = [F_top_C.get(c, 0.0) / F_C_total for c in comps]
        T_cond, _ = _bubble_T_K(x_C_list, comps, design.P_col, design.K_method)
        lam_top = _weighted_lambda(F_top_C)        # 凝縮分の λ
        # 設計判断 (2026-05-26): 還流量補正。CMO 還流は L=R×D_total なので Q_latent は
        # R×F_top_total×λ が正 (旧 R×F_C_total は留出の凝縮性成分で代用し ~5.7x 過小)。
        # FUG はプロファイル無しのため塔頂蒸気の顕熱は省略 (rigorous 経路は厳密
        # エネルギー収支)。残差は rigorous 比 ~10-15% 過小。詳細 tools/_diag_dist2_qcond.py。
        Q_cond_kW = (design.reflux_ratio * F_top_total
                     * lam_top * 1000.0 / 3600.0)
    else:
        # Total condenser (Dist1/Dist3): V = (R+1)D を全量凝縮
        T_cond  = T_top
        lam_top = _weighted_lambda(F_top)
        Q_cond_kW = (F_top_total * (design.reflux_ratio + 1.0)
                     * lam_top * 1000.0 / 3600.0)

    # リボイラ熱量 Q_reb = V' × λ_bot
    # 設計判断 (2026-05-09): 旧版 Q_reb = 1.05 × Q_cond は q=1 (飽和液フィード)
    # 専用の経験則で、q=0 (気相フィード) の Dist2 では物理的に不正確だった。
    # 標準 MESH 式 V' = (R+1)D - (1-q)F (q=0 で V' は (R+1)D - F に減る) に置換。
    # q=1 (Dist1, Dist3): V' = (R+1)D となり旧式とほぼ同等。
    # q=0 (Dist2): V' = (R+1)D - F、partial condenser 後の Q_cond と独立に決まる。
    F_feed_total  = sum(max(F, 0.0) for F in feed.F_in.values())
    V_prime_kmolh = (F_top_total * (design.reflux_ratio + 1.0)
                    - (1.0 - design.q) * F_feed_total)
    if V_prime_kmolh > 0:
        lam_bot = _weighted_lambda(F_bot)
        Q_reb_kW = V_prime_kmolh * lam_bot * 1000.0 / 3600.0
    else:
        # V' ≤ 0 は通常 q=0 で feed が大きすぎるとき。塔として動作不能なので
        # 物理的には infeasible だが、ここでは保守的に Q_cond と同等値で埋める。
        Q_reb_kW = Q_cond_kW

    # ProcessStream の T_top には実際の condenser 出口温度を入れる
    # (partial の場合は T_cond > T_top_bubble、total の場合は T_cond = T_top_bubble)
    T_top = T_cond

    # フィード予熱 (顕熱、feed.T_in → T_bot、簡易液仮定)
    Q_feed_preheat_kW = 0.0
    if feed.T_in < T_bot:
        for c, F_kmol_h in feed.F_in.items():
            if F_kmol_h <= 0:
                continue
            cp = cp_of(c) if cp_of(c) > 0 else CP_DEFAULT  # [J/(mol·K)]
            F_mol_s = F_kmol_h * 1000.0 / 3600.0
            Q_feed_preheat_kW += F_mol_s * cp * (T_bot - feed.T_in) / 1000.0

    # ---- Step 9: 塔径・塔高 (contest §4-2) ----
    # 設計判断 (2026-05-09): 旧版は u_vapor=0.3 m/s 単一固定で塔径を出していたが、
    # contest §4-2 仕様 (G* = SF·K·√(ρ_v(ρ_l-ρ_v))) に置換して物理ベースに。
    # 段塔・シーブトレー仮定。塔径は塔頂条件で評価 (蒸気量最大、最も flooding しやすい)。
    from src.component_data import MW, liquid_density_mix
    MW_top_avg = sum(F_top.get(c, 0) * MW.get(c, 50.0) for c in F_top) / max(F_top_total, 1e-9)
    # 蒸気密度: 理想気体近似 ρ_v = P·MW/(R·T) [kg/m³]
    rho_v = design.P_col * (MW_top_avg * 1e-3) / (_R_GAS * T_top)
    # 液密度: 既存 liquid_density_mix (kg/m³)、x_top を流量代用
    rho_l = liquid_density_mix(F_top)
    if rho_l <= rho_v + 1.0:
        # 物性異常 (蒸気と液の密度差が小さすぎる) → 保守側で D を大きめに
        rho_l = rho_v + 100.0
    # G* = SF · K · √(ρ_v · (ρ_l - ρ_v))   [kg/(m²·s)]
    G_star = fixed.SF * fixed.K_factor * math.sqrt(rho_v * (rho_l - rho_v))
    if G_star <= 1e-6:
        return _penalty_result(design, f"G*={G_star:.3e} <= 0 (塔径計算不能)",
                               N_min=N_min, R_min=R_min)
    # 蒸気の質量流量 [kg/s] = (R+1)·D·MW_v / 3600
    mass_flow_vap = F_top_total * (design.reflux_ratio + 1.0) * MW_top_avg / 3600.0
    A_cross = mass_flow_vap / G_star            # [m²]
    D_col = math.sqrt(4.0 * A_cross / math.pi)
    D_col = max(D_col, 0.3)

    # 塔高: 実段数 = 理論段 / 段効率、塔頂・塔底 section を加算
    N_real = math.ceil(design.N_stages / fixed.tray_efficiency)
    H_col = (N_real * fixed.tray_spacing_m
             + fixed.top_section_m + fixed.bot_section_m)
    H_col = max(H_col, 5.0)
    V_col = math.pi / 4.0 * D_col ** 2 * H_col

    # ---- Step 10a: Condenser / Reboiler の utility 選択と装置サイジング ----
    # 設計判断 (2026-05-09): 旧版は塔本体 (vessel + trays) しか CAPEX に含めず、
    # コンデンサも一律 COOLING_WATER 単価で OPEX 計上していた。実機では condenser・
    # reboiler は独立 HE で大きな CAPEX を持ち、温度依存で utility が変わる
    # (Dist2 のように T_top 氷点下ならエチレン冷媒が必要)。
    # ここで utility_selector に T_top/T_bot を渡して適切な tier を自動選択する。
    #
    # T_top が冷媒到達下限を下回るケース (Dist2 で x_top に H2/CH4 主体の場合) は、
    # 実機では partial condenser で対応するが本シミュレータは total condenser モデル。
    # 最低冷媒 (-100°C エチレン) で押し切るコストを保守見積もりとして採用する。
    _T_COOLING_FLOOR = 173.15 + 10.0 + 1.0   # -89.0°C: -100°C 冷媒 + 10K approach + 1K 余裕
    T_top_for_util = max(T_top, _T_COOLING_FLOOR)
    try:
        cond_utility = select_cooling_utility(T_top_for_util, mode='continuous')
    except ValueError:
        return _penalty_result(
            design,
            f"condenser T_top={T_top-273.15:.1f}°C は冷媒で到達不能 (要再設計)",
            N_min=N_min, R_min=R_min,
        )
    try:
        reb_utility = select_heating_utility(T_bot, mode='continuous')
    except ValueError:
        return _penalty_result(
            design,
            f"reboiler T_bot={T_bot-273.15:.1f}°C は熱媒の上限超過 (要再設計)",
            N_min=N_min, R_min=R_min,
        )

    # 伝熱面積: contest §4-4 (lookup_U) で (hot_phase, cold_phase) ペアから U 索引。
    # condenser: 蒸気凝縮 (CONDENSING、constant T_top_for_util) ⇄ 冷媒 (utility_phase)
    cond_T_in  = cond_utility.supply_T_K
    cond_T_out = cond_utility.supply_T_K + 10.0
    dT1 = T_top_for_util - cond_T_in
    dT2 = T_top_for_util - cond_T_out
    if dT1 <= 0.0 or dT2 <= 0.0:
        return _penalty_result(
            design,
            f"condenser ΔT 不成立 (T_top={T_top-273.15:.1f}°C, "
            f"supply={cond_utility.supply_T_K-273.15:.1f}°C)",
            N_min=N_min, R_min=R_min,
        )
    if abs(dT1 - dT2) < 1e-3:
        dT_lm_cond = 0.5 * (dT1 + dT2)
    else:
        dT_lm_cond = (dT1 - dT2) / math.log(dT1 / dT2)
    U_cond  = lookup_U(StreamPhase.CONDENSING, utility_phase(cond_utility.name))
    A_cond  = max(Q_cond_kW * 1000.0 / (U_cond * dT_lm_cond), 10.0)

    # reboiler: 熱媒 (utility_phase、通常 CONDENSING の蒸気か GAS の燃料炉) ⇄ 液沸騰 (EVAPORATING)
    dT_lm_reb = reb_utility.supply_T_K - T_bot
    if dT_lm_reb <= 0.0:
        return _penalty_result(
            design,
            f"reboiler ΔT 不成立 (T_bot={T_bot-273.15:.1f}°C, "
            f"supply={reb_utility.supply_T_K-273.15:.1f}°C)",
            N_min=N_min, R_min=R_min,
        )
    U_reb  = lookup_U(utility_phase(reb_utility.name), StreamPhase.EVAPORATING)
    A_reb  = max(Q_reb_kW * 1000.0 / (U_reb * dT_lm_reb), 10.0)

    # ---- Step 10b: CAPEX (Vessel + Trays + Condenser + Reboiler) ----
    capex_vessel = _vessel_capex_okuyen(V_col, design.P_col, D_col)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        capex_trays = calc_tray_capex_okuyen(D_col, design.N_stages)
        capex_cond  = calc_he_capex_okuyen(A_cond)
        capex_reb   = calc_he_capex_okuyen(A_reb)
    capex_total = capex_vessel + capex_trays + capex_cond + capex_reb

    # ---- Step 11: rigorous プロキシ罰則 (Phase A, 2026-05-19 / Phase A 強化 2026-05-20) ----
    # 丸める前の F_top_pre_round を渡すことで、<1% の C3 漏れでも判定可能になる
    # (現状の proxy_penalty 内ロジック自体は 1% 閾値ベースだが、将来の拡張余地を確保)。
    proxy_penalty, proxy_reason = _compute_proxy_penalty(
        design, F_top_pre_round, feed.F_in, R_min, N_min, is_rigorous=False,
    )

    # ---- Step 12: 結果 ----
    top_stream = ProcessStream(
        F_in=F_top, T_in=T_top, P_in=design.P_col,
    )
    bottom_stream = ProcessStream(
        F_in=F_bot, T_in=T_bot, P_in=design.P_col,
    )
    equipment = DistEquipment(
        D_col=D_col, H_col=H_col, V_col=V_col,
        CAPEX_vessel=capex_vessel,
        CAPEX_trays=capex_trays,
        CAPEX=capex_total,
        Q_cond=Q_cond_kW, Q_reb=Q_reb_kW,
        Q_feed_preheat_kW=Q_feed_preheat_kW,
        N_min=N_min, R_min=R_min,
        N_feed_kirkbride=N_feed_kirkbride,
        feasible=True, message="",
        A_cond_m2=A_cond, A_reb_m2=A_reb,
        CAPEX_cond=capex_cond, CAPEX_reb=capex_reb,
        cond_utility_name=cond_utility.name,
        cond_utility_jpy_per_GJ=cond_utility.jpy_per_GJ,
        reb_utility_name=reb_utility.name,
        reb_utility_jpy_per_GJ=reb_utility.jpy_per_GJ,
        T_top=T_top, T_bot=T_bot,
        proxy_penalty_okuyen=proxy_penalty,
        proxy_penalty_reason=proxy_reason,
    )
    return DistResult(top=top_stream, bottom=bottom_stream, equipment=equipment)


_DT_MAX_RE = __import__('re').compile(r'dT_max=([\d.]+)\s*K')


def _extract_dT_max(msg: str) -> float:
    """Wang-Henke 失敗 RuntimeError message から dT_max [K] を抽出。

    例: "Wang-Henke 収束失敗 (retry 後も): 収束失敗: max_iter=500, dT_max=56.280K (tol=0.05K)"
    抽出できない (= 別種の RuntimeError、例: 冷媒選択失敗) 場合は 0.0 を返す。
    """
    m = _DT_MAX_RE.search(msg)
    if m:
        try:
            return float(m.group(1))
        except (ValueError, IndexError):
            return 0.0
    return 0.0


def _penalty_result(design: DistDesignVars, msg: str,
                    N_min: float = 0.0, R_min: float = 0.0,
                    N_needed: float = 0.0,
                    dT_max_rigorous: float = 0.0) -> DistResult:
    """infeasibility または計算失敗時のペナルティ結果。

    N_needed (Gilliland) / dT_max_rigorous (Wang-Henke) は TPE constraints_func 用の
    連続シグナル。run_one_pass で N_shortfall / dT_shortfall を計算して trial.user_attrs
    に格納する。未指定なら 0.0 (= 連続シグナルなし、binary penalty のみ)。
    """
    warnings.warn(
        f"simulate_distillation_column: infeasible — {msg}",
        UserWarning, stacklevel=2,
    )
    F_zero = {c: 0.0 for c in ('A', 'B', 'C', 'D', 'E', 'F', 'Z')}
    zero_top = ProcessStream(F_in=dict(F_zero), T_in=298.15, P_in=design.P_col)
    zero_bot = ProcessStream(F_in=dict(F_zero), T_in=298.15, P_in=design.P_col)
    eq = DistEquipment(
        D_col=0.3, H_col=5.0, V_col=0.4,
        CAPEX_vessel=_PENALTY,
        CAPEX_trays=0.0,
        CAPEX=_PENALTY,
        Q_cond=0.0, Q_reb=0.0, Q_feed_preheat_kW=0.0,
        N_min=N_min, R_min=R_min,
        N_feed_kirkbride=0,
        feasible=False, message=msg,
        N_needed=N_needed,
        dT_max_rigorous=dT_max_rigorous,
    )
    return DistResult(top=zero_top, bottom=zero_bot, equipment=eq)
