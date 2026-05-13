"""
HEN (Heat Exchanger Network) Synthesis — Stage 2 of heat integration

設計判断 (2026-05-09): Linnhoff の Pinch Design Method (PDM) を自動化。
  Stage 1 (flowsheet/heat_integration.py:pinch_analysis) は targeting のみ
  (Q_H_min, Q_C_min, A_total を理論値として返す)。
  本モジュール Stage 2 では、ストリーム同士の実マッチングを決定して
  「実 HEN 構成・各 HE のサイズ・追加 CAPEX」を計算する。

実装段階 (2026-05-09 時点):
  Phase 2A (本実装、~simple):
    - ピンチ温度で stream を上下分割
    - 各サイドで greedy + tick-off で matching
    - 残り duty は utility 処理 (Stage 1 と同じ tier 配分)
    - LMTD 法で各 HE の A・CAPEX 計算
    - CP rule (FCp_hot ≤ FCp_cold above pinch 等) は厳密チェックせず、
      LMTD > 0 と ΔT > dT_min/2 で feasibility のみ確認
  Phase 2B (将来):
    - CP rule 厳守 + stream split で N_HE_min 厳密達成
    - 教科書解レベル

参考文献:
  [1] 長谷部 伸治, 外輪 健一郎『プロセスシステム工学 (No.4) — 熱交換器
      ネットワークの最適合成』京都大学講義資料 (2025) §4.7
  [2] Linnhoff B., Hindmarsh E., Chem. Eng. Sci. 38 (1983) 745-763
      "The pinch design method for heat exchanger networks"

新規仮定:
  - **既存 HE は撤去せず残す** (物理機器として CAPEX 計上済み)
  - **Stage 2 は追加 HE のみ計算** (recovery exchanger 等、新規 process-process
    matching 用)。既存 HE は utility 接続のままと見なす。
  - これにより top-k TAC = 既存 CAPEX + 追加 HE CAPEX + 実 OPEX となり、
    Stage 1 TAC (= 既存 CAPEX + 理論 OPEX) より大きい構図 (= ユーザー想定通り)
  - LMTD 計算で逆流仮定 (counter-current)
  - U = h_W_m2K (各 HIStream に格納された値) の単純加算なし、両側の
    調和平均的な処理は省略 → 入口側の h を採用 (簡略化)
  - 潜熱 stream (T_top == T_bot) は constant-T 側として扱う

依存:
  - flowsheet/heat_integration.py: HIStream, HIResult, UtilityTier
  - src/cost_calculator.py: calc_he_capex_okuyen
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import math
import warnings

from flowsheet.heat_integration import (
    HIStream, HIResult, UtilityTier, StreamPhase, lookup_U,
    classify_heat_opex_key, calc_hi_opex_okuyen,
)
from src.cost_calculator import calc_he_capex_okuyen


# 設計判断: 数値計算上の許容
_EPS = 1e-9
_DT_MIN_FLOOR = 0.5    # K, ΔT 計算の数値ノイズ下限


# ===========================================================================
# データクラス
# ===========================================================================

@dataclass
class HEMatch:
    """1 個の熱交換器 (process-process or process-utility) の合成結果。"""
    name:          str          # 識別子 (e.g., 'HE_recovery_above_1')
    hot_label:     str          # hot stream 名 (e.g., 'H1_rx_cooler')
    cold_label:    str          # cold stream 名 (e.g., 'C2_rx_preheat')
    Q_kW:          float        # 熱交換量
    T_h_in_K:      float        # hot 入口
    T_h_out_K:     float        # hot 出口
    T_c_in_K:      float        # cold 入口
    T_c_out_K:     float        # cold 出口
    LMTD_K:        float        # 対数平均温度差
    U_W_m2K:       float        # 総括伝熱係数 (両 stream の h の小さい方)
    A_m2:          float        # 伝熱面積
    CAPEX_okuyen:  float        # この 1 機の CAPEX [億円]
    side:          str          # 'above' | 'below' | 'utility_hot' | 'utility_cold'


@dataclass
class HENResult:
    """Stage 2 synthesis 結果。"""
    matches:               List[HEMatch]            # 全 HE
    n_process_HE:          int                      # process-process HE 数
    n_utility_HE:          int                      # process-utility HE 数
    Q_recovered_kW:        float                    # 内部熱回収量 (= sum of process-process Q)
    Q_hot_utility_kW:      float                    # 実 hot utility 必要量 (≥ Q_H_min)
    Q_cold_utility_kW:     float                    # 実 cold utility 必要量
    CAPEX_added_okuyen:    float                    # 追加 HE 分の合計 CAPEX
    OPEX_utility_okuyen:   Dict[str, float]         # utility tier 別 OPEX [億円/年]
    feasible:              bool
    message:               str = ""


# 内部用: stream を pinch で分割した segment
@dataclass
class _Segment:
    name:       str
    is_hot:     bool
    T_top_K:    float       # 高温端
    T_bot_K:    float       # 低温端 (latent なら == T_top)
    F_Cp_kW_K:  float       # 0 なら latent
    Q_kW:       float       # この segment の総 Q
    Q_remain:   float       # tick-off 用残量
    h_W_m2K:    float
    is_latent:  bool
    phase:      str = StreamPhase.LIQUID    # contest §4-4 索引用


# ===========================================================================
# Stream 分割: pinch で上下に切る
# ===========================================================================

def _split_at_pinch(
    streams:        List[HIStream],
    T_pinch_hot_K:  float,        # 実温度の pinch (hot 側)
    T_pinch_cold_K: float,        # 実温度の pinch (cold 側)
) -> Tuple[List[_Segment], List[_Segment], List[_Segment], List[_Segment]]:
    """各 stream を pinch 温度で上下に分割。

    Returns
    -------
    above_hot, above_cold, below_hot, below_cold
        ピンチ上下それぞれの hot/cold segment リスト。
    """
    above_hot:  List[_Segment] = []
    above_cold: List[_Segment] = []
    below_hot:  List[_Segment] = []
    below_cold: List[_Segment] = []

    for s in streams:
        T_top = max(s.T_in_K, s.T_out_K)
        T_bot = min(s.T_in_K, s.T_out_K)
        is_hot = s.is_hot

        # Pinch を hot/cold 側で読み替え
        T_pinch = T_pinch_hot_K if is_hot else T_pinch_cold_K

        # ---- 顕熱 segment ----
        if s.F_Cp_kW_per_K > 0 and abs(T_top - T_bot) > _EPS:
            # 上下にまたがるなら分割、片側だけなら 1 つ
            if T_top > T_pinch + _EPS and T_bot < T_pinch - _EPS:
                # 跨ぐ → 上下に分割
                Q_above = s.F_Cp_kW_per_K * (T_top - T_pinch)
                Q_below = s.F_Cp_kW_per_K * (T_pinch - T_bot)
                seg_above = _Segment(
                    name=s.name + '_above', is_hot=is_hot,
                    T_top_K=T_top, T_bot_K=T_pinch,
                    F_Cp_kW_K=s.F_Cp_kW_per_K,
                    Q_kW=Q_above, Q_remain=Q_above,
                    h_W_m2K=s.h_W_m2K, is_latent=False,
                    phase=s.phase,
                )
                seg_below = _Segment(
                    name=s.name + '_below', is_hot=is_hot,
                    T_top_K=T_pinch, T_bot_K=T_bot,
                    F_Cp_kW_K=s.F_Cp_kW_per_K,
                    Q_kW=Q_below, Q_remain=Q_below,
                    h_W_m2K=s.h_W_m2K, is_latent=False,
                    phase=s.phase,
                )
                (above_hot if is_hot else above_cold).append(seg_above)
                (below_hot if is_hot else below_cold).append(seg_below)
            else:
                # 片側だけ
                Q = s.F_Cp_kW_per_K * (T_top - T_bot)
                seg = _Segment(
                    name=s.name, is_hot=is_hot,
                    T_top_K=T_top, T_bot_K=T_bot,
                    F_Cp_kW_K=s.F_Cp_kW_per_K,
                    Q_kW=Q, Q_remain=Q,
                    h_W_m2K=s.h_W_m2K, is_latent=False,
                    phase=s.phase,
                )
                if T_bot >= T_pinch - _EPS:
                    (above_hot if is_hot else above_cold).append(seg)
                else:
                    (below_hot if is_hot else below_cold).append(seg)

        # ---- 潜熱 segment (constant T) ----
        if s.Q_latent_kW > 0:
            T_phase = s.T_phase_K if s.T_phase_K is not None else 0.5*(T_top+T_bot)
            seg = _Segment(
                name=s.name + ('_latent' if s.F_Cp_kW_per_K > 0 else ''),
                is_hot=is_hot,
                T_top_K=T_phase, T_bot_K=T_phase,
                F_Cp_kW_K=0.0,
                Q_kW=s.Q_latent_kW, Q_remain=s.Q_latent_kW,
                h_W_m2K=s.h_W_m2K, is_latent=True,
                phase=(StreamPhase.CONDENSING if is_hot else StreamPhase.EVAPORATING),
            )
            if T_phase >= T_pinch - _EPS:
                (above_hot if is_hot else above_cold).append(seg)
            else:
                (below_hot if is_hot else below_cold).append(seg)

    return above_hot, above_cold, below_hot, below_cold


# ===========================================================================
# 1 マッチの計算: ΔTmin チェック・LMTD・面積・CAPEX
# ===========================================================================

def _try_match(
    h: _Segment,
    c: _Segment,
    dT_min_K: float,
    side:     str,
    name:     str,
) -> Optional[HEMatch]:
    """hot segment h と cold segment c の matching を試行。

    可能なら HEMatch を作って返す。不可能なら None。

    マッチ量 Q = min(Q_h_remain, Q_c_remain)。
    Tick-off rule: 1 stream を完全消化させる。
    """
    Q = min(h.Q_remain, c.Q_remain)
    if Q < _EPS:
        return None

    # ΔT_in / ΔT_out の計算 (counter-current 仮定)
    # hot 側で Q を提供する範囲: T_top → T_top - Q/FCp (latent なら T_phase)
    # cold 側で Q を吸収する範囲: T_bot → T_bot + Q/FCp (latent なら T_phase)
    if h.is_latent:
        T_h_in  = h.T_top_K
        T_h_out = h.T_top_K
    else:
        T_h_in  = h.T_top_K
        T_h_out = h.T_top_K - Q / h.F_Cp_kW_K

    if c.is_latent:
        T_c_in  = c.T_bot_K
        T_c_out = c.T_bot_K
    else:
        T_c_in  = c.T_bot_K
        T_c_out = c.T_bot_K + Q / c.F_Cp_kW_K

    # 端点 ΔT (counter: hot 入口 vs cold 出口、hot 出口 vs cold 入口)
    dT_hot_end  = T_h_in  - T_c_out
    dT_cold_end = T_h_out - T_c_in

    if dT_hot_end < dT_min_K - _EPS or dT_cold_end < dT_min_K - _EPS:
        return None    # ΔTmin 違反

    # LMTD (どちらも正、近接なら平均)
    if abs(dT_hot_end - dT_cold_end) < _DT_MIN_FLOOR:
        LMTD = 0.5 * (dT_hot_end + dT_cold_end)
    else:
        LMTD = (dT_hot_end - dT_cold_end) / math.log(dT_hot_end / dT_cold_end)
    if LMTD <= 0:
        return None

    # 総括 U: contest §4-4 表 (lookup_U) で両 segment の phase 組合せから決定。
    # Ref: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 (流速によらず固定 U)。
    U_W_m2K = lookup_U(h.phase, c.phase)

    # 面積 [m²]
    A = max(Q * 1000.0 / (U_W_m2K * LMTD), 10.0)   # 下限 10 m² (HE CAPEX 適用範囲)

    # CAPEX [億円]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        CAPEX = calc_he_capex_okuyen(A)

    return HEMatch(
        name=name,
        hot_label=h.name, cold_label=c.name,
        Q_kW=Q,
        T_h_in_K=T_h_in,   T_h_out_K=T_h_out,
        T_c_in_K=T_c_in,   T_c_out_K=T_c_out,
        LMTD_K=LMTD, U_W_m2K=U_W_m2K, A_m2=A,
        CAPEX_okuyen=CAPEX, side=side,
    )


# ===========================================================================
# Greedy + Tick-off matching (ピンチ片側専用)
# ===========================================================================

def _greedy_match_one_side(
    hot_segs:  List[_Segment],
    cold_segs: List[_Segment],
    dT_min_K:  float,
    side:      str,
    name_prefix: str = "HE",
) -> Tuple[List[HEMatch], List[_Segment], List[_Segment]]:
    """ピンチ片側で hot/cold をマッチング。

    Greedy: 残量最大の hot と残量最大の cold を試行。feasible なら確定、
    そうでなければ次の組合せ。tick-off で 1 segment が使い切られたら除外。

    Returns
    -------
    matches, hot_unmatched, cold_unmatched
        unmatched は Q_remain > 0 で残った segment (utility 処理対象)。
    """
    matches: List[HEMatch] = []
    counter = 0

    # コピーで作業 (元 list を変更しない)
    hot  = [s for s in hot_segs  if s.Q_remain > _EPS]
    cold = [s for s in cold_segs if s.Q_remain > _EPS]

    while hot and cold:
        # 残量大きい順に sort
        hot.sort(key=lambda s: -s.Q_remain)
        cold.sort(key=lambda s: -s.Q_remain)

        matched = False
        # 全候補 pair から最初に feasible を採用 (greedy)
        for h in hot:
            for c in cold:
                counter += 1
                m = _try_match(h, c, dT_min_K,
                               side=side, name=f"{name_prefix}_{side}_{counter}")
                if m is None:
                    continue
                # 確定: tick-off
                matches.append(m)
                h.Q_remain -= m.Q_kW
                c.Q_remain -= m.Q_kW
                matched = True
                break
            if matched:
                break

        if not matched:
            # どの pair も feasible ない → 残りは utility へ
            break

        # 消化された segment を除外
        hot  = [s for s in hot  if s.Q_remain > _EPS]
        cold = [s for s in cold if s.Q_remain > _EPS]

    return matches, hot, cold


# ===========================================================================
# 公開関数: HEN 合成
# ===========================================================================

def synthesize_hen(
    streams:        List[HIStream],
    hi_result:      HIResult,
    dT_min_K:       float = 10.0,
    heating_tiers:  Optional[List[UtilityTier]] = None,
    cooling_tiers:  Optional[List[UtilityTier]] = None,
    operating_hours: float = 8000.0,
) -> HENResult:
    """HEN を Pinch Design Method で合成する (Phase 2A)。

    Parameters
    ----------
    streams : List[HIStream]
        flowsheet.heat_integration.extract_streams() の戻り値
    hi_result : HIResult
        Stage 1 (pinch_analysis) の結果。pinch 温度を流用
    dT_min_K : float
        最小接近温度差 (Stage 1 と同値)
    heating_tiers, cooling_tiers : List[UtilityTier]
        utility tier 定義 (None なら OPEX_utility は空)
    operating_hours : float
        年間稼働時間 [h/年]

    Returns
    -------
    HENResult
        全 HE 構成・追加 CAPEX・実 utility OPEX
    """
    # ---- 1. ピンチで分割 ----
    above_hot, above_cold, below_hot, below_cold = _split_at_pinch(
        streams, hi_result.T_pinch_hot_K, hi_result.T_pinch_cold_K,
    )

    # ---- 2. 上下それぞれ greedy match ----
    matches_above, hot_un_above, cold_un_above = _greedy_match_one_side(
        above_hot, above_cold, dT_min_K, side='above',
    )
    matches_below, hot_un_below, cold_un_below = _greedy_match_one_side(
        below_hot, below_cold, dT_min_K, side='below',
    )

    process_matches = matches_above + matches_below

    # ---- 3. 残り duty を utility 接続として記録 ----
    # 上側の hot 残り → cold utility (空冷など) で冷却
    # 上側の cold 残り → hot utility (LP Steam など) で加熱
    # 下側も同様
    Q_hot_unmatched_above  = sum(s.Q_remain for s in hot_un_above)
    Q_cold_unmatched_above = sum(s.Q_remain for s in cold_un_above)
    Q_hot_unmatched_below  = sum(s.Q_remain for s in hot_un_below)
    Q_cold_unmatched_below = sum(s.Q_remain for s in cold_un_below)

    Q_cold_utility = Q_hot_unmatched_above + Q_hot_unmatched_below   # 冷却が必要な hot の残り
    Q_hot_utility  = Q_cold_unmatched_above + Q_cold_unmatched_below # 加熱が必要な cold の残り

    # ---- 4. utility 用 OPEX (Stage 1 hi_result.utility_breakdown を使い回し) ----
    # 設計判断: synthesis 後の utility 必要量は targeting (Q_H_min/Q_C_min) と
    # 一致するのが理論的だが、greedy で完全達成できないと差分が出る。
    # 簡略化: 配分は Stage 1 の tier breakdown をそのまま流用、合計 Q だけ
    # 補正してスケール。
    OPEX_utility: Dict[str, float] = {}
    if heating_tiers is not None and cooling_tiers is not None:
        from src.cost_parameters import FUEL_JPY_PER_GJ, LP_STEAM_JPY_PER_GJ  # noqa
        # tier 配分: hi_result.utility_breakdown (Stage 1 targeting) を採用、
        # ただし合計を実 Q (=Q_hot_utility, Q_cold_utility) にスケール
        target_hot_total  = hi_result.Q_H_min_kW
        target_cold_total = hi_result.Q_C_min_kW
        scale_hot  = (Q_hot_utility  / target_hot_total)  if target_hot_total  > _EPS else 1.0
        scale_cold = (Q_cold_utility / target_cold_total) if target_cold_total > _EPS else 1.0

        all_tiers = {t.name: t for t in heating_tiers + cooling_tiers}
        for tier_name, Q_target in hi_result.utility_breakdown.items():
            tier = all_tiers.get(tier_name)
            if tier is None:
                continue
            scale = scale_hot if tier.is_heating else scale_cold
            Q_actual = Q_target * scale
            cost = Q_actual * 3.6e-3 * operating_hours * tier.jpy_per_GJ / 1.0e8
            OPEX_utility[tier_name] = cost

    # ---- 5. 集計 ----
    Q_recovered = sum(m.Q_kW for m in process_matches)
    CAPEX_added = sum(m.CAPEX_okuyen for m in process_matches)

    # feasibility: greedy で達成できなかった分を判定
    feasible = True
    msg = ""
    Q_recovery_target = (
        sum(s.Q_total_kW for s in streams if s.is_hot)
        - hi_result.Q_C_min_kW
    )
    if Q_recovered < 0.95 * Q_recovery_target:
        feasible = False
        msg = (f"greedy で目標 {Q_recovery_target:.0f} kW のうち "
               f"{Q_recovered:.0f} kW しか回収できず ({Q_recovered/Q_recovery_target*100:.1f}%)")

    return HENResult(
        matches=process_matches,
        n_process_HE=len(process_matches),
        n_utility_HE=0,    # 簡略化のため utility 接続を 1 機ずつ数えない
        Q_recovered_kW=Q_recovered,
        Q_hot_utility_kW=Q_hot_utility,
        Q_cold_utility_kW=Q_cold_utility,
        CAPEX_added_okuyen=CAPEX_added,
        OPEX_utility_okuyen=OPEX_utility,
        feasible=feasible,
        message=msg,
    )


# ===========================================================================
# Economics への適用
# ===========================================================================

def apply_synthesis_to_economics(
    economics,                          # flowsheet.economics.Economics (Stage 1 適用前)
    hen_result:           HENResult,
    operating_hours:      float = 8000.0,
    depreciation_years:   int   = 8,
):
    """既存 Economics に Stage 2 (synthesis) 結果を適用した新 Economics を返す。

    扱い:
      - 既存 HE CAPEX (Cooler/Intercool/Desuper/MemPrecool/Mem冷却器/Dist cond+reb 等)
        は **そのまま据え置き** (物理機器として残置)
      - HEN synthesis で出た **追加 HE CAPEX を加算**
      - 熱系 OPEX を synthesis 後の実 utility breakdown で置換
        (Stage 1 の targeting OPEX より 5-10% 悪化するのが通常)
      - 非熱系 OPEX (触媒・原料費等) はそのまま継承

    Returns
    -------
    Economics
        Stage 2 適用後の新インスタンス
    """
    from copy import deepcopy
    from flowsheet.economics import (
        Economics, apply_hasebe_aggregation, _count_main_equipment,
    )

    # ---- CAPEX: 既存 + HEN 追加 ----
    new_capex = deepcopy(economics.capex)
    if hen_result.CAPEX_added_okuyen > 0:
        new_capex['HEN追加 (process-process HE)'] = hen_result.CAPEX_added_okuyen

    # ---- OPEX: 熱系を synthesis 後 tier 別 OPEX で置換、非熱系は継承 ----
    # Hasebe 集計項は apply_hasebe_aggregation 内部で剥がして再計算する。
    new_opex: Dict[str, float] = {}
    for k, v in economics.opex.items():
        if classify_heat_opex_key(k) is None:
            new_opex[k] = v
    for tier_name, cost in hen_result.OPEX_utility_okuyen.items():
        new_opex[f'Stage2: {tier_name}'] = cost

    # ---- 集計 ----
    new_total_capex = sum(v for v in new_capex.values() if v < 1e6)

    # Hasebe 集計項を再計算 (CAPEX が HEN 追加で増えれば 0.180·C_TM が増、
    # C_UT も Stage2 置換で変化するため両方追従)
    N_eq     = _count_main_equipment(new_capex)
    new_opex = apply_hasebe_aggregation(new_opex, new_total_capex, N_eq)

    new_total_opex  = sum(new_opex.values())
    new_TAC         = new_total_capex / depreciation_years + new_total_opex
    new_profit      = economics.total_revenue - new_TAC
    new_unit_jpy_per_t = (
        new_TAC * 1.0e8 / (economics.annual_kg_C3H6 / 1000.0)
        if economics.annual_kg_C3H6 > 0 else float('inf')
    )

    return Economics(
        capex          =new_capex,
        opex           =new_opex,
        revenue        =deepcopy(economics.revenue),
        total_capex    =new_total_capex,
        total_opex     =new_total_opex,
        total_revenue  =economics.total_revenue,
        TAC            =new_TAC,
        profit         =new_profit,
        annual_kg_C3H6 =economics.annual_kg_C3H6,
        unit_jpy_per_t =new_unit_jpy_per_t,
    )
