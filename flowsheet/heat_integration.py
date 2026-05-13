"""
ヒートインテグレーション (HI) ターゲティング・モジュール

設計判断 (2026-05-08):
  最終的に BO で全体最適化を行う想定で、BO ループ内では「流体ペアリング」を
  行わずに、熱力学的な理論限界値（targeting）だけを返す。
  ペアリングは離散的・組合せ問題なので BO の GP には不向き。
  最適点 x* に対してだけ Stage 2 (post-processing) で実ネットワーク合成を行う。

本モジュール (Stage 1) で計算するもの:
  1. 最少加熱量 Q_H_min, 最少冷却量 Q_C_min, ピンチ温度
       → Problem Table Algorithm（ΔTmin が決まれば一意）
  2. 与熱・受熱複合線, グランドコンポジットカーブ (GCC)
       → 可視化と tier 別 Q 配分の元データ
  3. ユーティリティ tier 別 Q 配分
       → GCC を tier 温度で水平に切って多 tier 配分 (Linnhoff 標準手法)
       → economics 側で OPEX = Σ Q_t × cost_t を計算可能にする
  4. 理想総伝熱面積 A_total （Bath 式）
  5. 最少熱交換器数 N_HE_min （Linnhoff 式）
       → economics 側で CAPEX_HE = N_HE × 固定費 + A_total × 単価 を計算可能にする

参考文献:
  [1] 長谷部 伸治, 外輪 健一郎『プロセスシステム工学 (No.4) — 熱交換器ネット
      ワークの最適合成』京都大学 講義資料、令和7年度 (2025).
      §4.2 エクセルギー、§4.3 T-Q 線図、§4.4 最小接近温度差と最少必要加熱・
      冷却量、§4.5 熱複合線図の分解と最適システム合成、§4.6 エネルギー解析
      (グランドコンポジットカーブ)、§4.7 詳細な検討。
      例題 4.2/4.3 は tests/test_heat_integration.py で再現検証済み。
  [2] Linnhoff B., Hindmarsh E., "The pinch design method for heat exchanger
      networks," Chem. Eng. Sci. 38 (1983) 745-763.
      多 tier 配分・最少 HE 数式 N_HE = N_streams + N_utilities − 1 の出典。

出力は決定論的・滑らかで、ms オーダーの計算時間。BO 評価関数として直接使える。

潜熱の扱い:
  単一の相変化温度 T_phase_K に Q_latent_kW が集中する簡易モデル。
  教科書 [1] 例題 4.2 の蒸気（潜熱 2.0 kJ/g）も表現可能。
  VLE フル対応は将来課題。
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import math


# ===========================================================================
# 総括熱伝達係数 U [W/(m²·K)]: contest §4-4 表
# ===========================================================================
# 出典: 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 熱交換器
#   Q = U · A · ΔT_ln、流速によらず固定 U を使用、9 区分の表で参照。
#   沸点・露点を通過する場合は前後で分割計算。
#
# 設計判断 (2026-05-09): 本プロジェクト用にコンテスト値を採用 (引用)。
#   旧版は化工便覧の per-stream h ベースだったが、contest 仕様の方が
#   各 HE の「相のペア」で直接 U を引ける明確な体系のため切替え。
# ---------------------------------------------------------------------------

class StreamPhase:
    """HE stream の相分類 (contest §4-4 表索引キー)。"""
    GAS         = 'gas'         # 顕熱、ガス相
    LIQUID      = 'liquid'      # 顕熱、液相
    CONDENSING  = 'condensing'  # 潜熱、ガス凝縮
    EVAPORATING = 'evaporating' # 潜熱、液蒸発


# Contest §4-4 表 [W/(m²·K)]、indexed by (hot_phase, cold_phase)
U_TABLE_CONTEST: Dict[Tuple[str, str], float] = {
    (StreamPhase.GAS,         StreamPhase.GAS):           150.0,  # ガス-ガス
    (StreamPhase.LIQUID,      StreamPhase.GAS):           200.0,  # 液-ガス
    (StreamPhase.LIQUID,      StreamPhase.LIQUID):        300.0,  # 液-液
    (StreamPhase.CONDENSING,  StreamPhase.EVAPORATING): 1500.0,   # ガス凝縮-液蒸発
    (StreamPhase.GAS,         StreamPhase.LIQUID):        200.0,  # ガス-液
    (StreamPhase.CONDENSING,  StreamPhase.GAS):           500.0,  # ガス凝縮-ガス
    (StreamPhase.CONDENSING,  StreamPhase.LIQUID):       1000.0,  # ガス凝縮-液
    (StreamPhase.GAS,         StreamPhase.EVAPORATING):  500.0,   # ガス-液蒸発
    (StreamPhase.LIQUID,      StreamPhase.EVAPORATING): 1000.0,   # 液-液蒸発
}

# Contest 表外組合せ用フォールバック (保守側 = U 小さめで A 大きめに)
_U_FALLBACK_W_M2K = 150.0


def lookup_U(hot_phase: str, cold_phase: str) -> float:
    """Contest §4-4 U 表から hot/cold の組合せに対応する U を引く。

    表外 (例: 両方凝縮、両方蒸発、cold 側凝縮 等) はフォールバック値。
    """
    return U_TABLE_CONTEST.get((hot_phase, cold_phase), _U_FALLBACK_W_M2K)


def utility_phase(utility_name: str) -> str:
    """utility tier 名から相分類を推定する。

    - LP/MP/HP Steam: 凝縮潜熱 → CONDENSING
    - 燃料燃焼: 高温燃焼ガス → GAS
    - 空冷: 大気ガス → GAS
    - 冷却水・プロピレン+5C/+15C: 加熱される液 → LIQUID
    - プロピレン-25C/-40C・エチレン-60C/-75C/-100C: 蒸発冷媒 → EVAPORATING
    """
    name = utility_name.lower()
    if 'steam' in name:
        return StreamPhase.CONDENSING
    if '燃料' in utility_name or 'fuel' in name:
        return StreamPhase.GAS
    if '空冷' in utility_name or 'air' in name:
        return StreamPhase.GAS
    if '冷却水' in utility_name or 'cooling water' in name:
        return StreamPhase.LIQUID
    # 冷凍冷媒: 高温側 (+5°C, +15°C) は液想定、低温側 (-25°C 以下) は蒸発潜熱
    if '冷媒' in utility_name:
        if '+' in utility_name:
            return StreamPhase.LIQUID
        return StreamPhase.EVAPORATING
    # フォールバック: 液
    return StreamPhase.LIQUID


# 後方互換性のため h_W_m2K 用の代表値も残置 (Bath 式 _calc_A_total で使用)
# 設計判断 (2026-05-09): h は U 表とは別概念だが、簡略化のため U/2 を仮置き
# (1/U = 1/h_h + 1/h_c、両側等価仮定で h ≈ 2U)。targeting 推算の精度範囲。
H_LIQUID_LIQUID_W_M2K   = 600.0    # ≈ 2 × 300 (液-液)
H_HIGH_TEMP_GAS_W_M2K   = 300.0    # ≈ 2 × 150 (ガス-ガス)
H_VAPORIZER_W_M2K       = 1000.0   # ≈ 2 × 500 (蒸発、保守側)
H_CONDENSER_W_M2K       = 1000.0   # ≈ 2 × 500 (凝縮、保守側)


# ===========================================================================
# データ構造
# ===========================================================================

@dataclass
class HIStream:
    """ヒートインテグレーション対象の単一ストリーム。

    is_hot は T_in > T_out から自動判定する（与熱流体 = 冷却される側 = hot）。

    Attributes
    ----------
    name : str
        ストリーム識別子（"H1_rx_cooler" など）
    T_in_K : float
        入口温度 [K]
    T_out_K : float
        出口温度 [K]
    F_Cp_kW_per_K : float
        顕熱の流量×比熱容量 [kW/K]。Q_sensible = F_Cp × |T_out - T_in|
    Q_latent_kW : float
        潜熱量 [kW]（相変化があるとき。0 で潜熱なし）
    T_phase_K : float | None
        相変化温度 [K]。None の場合は T_in と T_out の中点に置く（簡易）
    h_W_m2K : float
        膜伝達係数 [W/m²/K]。Bath 式での A_total 算出に使用
    """
    name:           str
    T_in_K:         float
    T_out_K:        float
    F_Cp_kW_per_K:  float
    Q_latent_kW:    float = 0.0
    T_phase_K:      Optional[float] = None
    h_W_m2K:        float = H_LIQUID_LIQUID_W_M2K
    # contest §4-4 U 表索引用の相分類 (2026-05-09 追加)
    phase:          str   = StreamPhase.LIQUID

    @property
    def is_hot(self) -> bool:
        return self.T_in_K > self.T_out_K

    @property
    def Q_sensible_kW(self) -> float:
        return self.F_Cp_kW_per_K * abs(self.T_out_K - self.T_in_K)

    @property
    def Q_total_kW(self) -> float:
        return self.Q_sensible_kW + self.Q_latent_kW


@dataclass
class UtilityTier:
    """ユーティリティ tier。

    src.utility_selector.UtilityTier と同じ構造だが、循環 import を避けるため
    本モジュール内に再定義。呼び出し側で適宜変換する。
    """
    name:       str
    supply_T_K: float
    jpy_per_GJ: float
    is_heating: bool


@dataclass
class HIResult:
    """ターゲティング結果。

    economics 側はこの dataclass を受け取って:
      - OPEX_HI = Σ utility_breakdown[t] × cost_t × time
      - CAPEX_HE = N_HE_min × 固定費 + A_total_m2 × 単価
    を計算する。

    可視化用に composite curves と GCC も含める。
    """
    Q_H_min_kW:          float
    Q_C_min_kW:          float
    T_pinch_hot_K:       float       # ピンチ点（実温度、与熱流体側）
    T_pinch_cold_K:      float       # ピンチ点（実温度、受熱流体側）
    A_total_m2:          float
    N_HE_min:            int

    utility_breakdown:   Dict[str, float]                      # {tier_name: Q_kW}
    composite_hot:       List[Tuple[float, float]]             # [(H_kW, T_K), ...]
    composite_cold:      List[Tuple[float, float]]
    GCC:                 List[Tuple[float, float]]             # [(T_shifted_K, H_residual_kW), ...]

    feasible:            bool
    message:             str = ""


# ===========================================================================
# ストリーム展開: 顕熱区間と潜熱区間に分解
# ===========================================================================

@dataclass(frozen=True)
class _Segment:
    """内部用: ストリームを1区間に分解したもの（顕熱 or 潜熱）。

    潜熱区間は T_top_K == T_bot_K（水平線）。
    F_Cp は実 F·Cp（顕熱）または「潜熱を Δε=0.001K に集中させた等価 F·Cp」
    のいずれかだが、本実装では潜熱は別の総熱量で扱い、F_Cp_kW_per_K を 0 にする。
    """
    name:            str
    is_hot:          bool
    T_top_K:         float     # 区間の高温端 [K]
    T_bot_K:         float     # 区間の低温端 [K]
    F_Cp_kW_per_K:   float     # 顕熱区間なら正、潜熱区間なら 0
    Q_kW:            float     # 区間の総熱量 [kW] (絶対値)
    h_W_m2K:         float


def _expand_stream(s: HIStream) -> List[_Segment]:
    """ストリームを顕熱+潜熱の区間に展開する。

    潜熱がない場合は1区間（顕熱のみ）。
    潜熱がある場合は、相変化温度 T_phase で 2 区間（顕熱）+ 1 区間（潜熱）
    に分割する。T_phase が T_in/T_out の範囲外なら片方が消える。
    """
    T_top = max(s.T_in_K, s.T_out_K)
    T_bot = min(s.T_in_K, s.T_out_K)

    if s.Q_latent_kW <= 0.0:
        # 顕熱のみ
        return [_Segment(
            name=s.name, is_hot=s.is_hot,
            T_top_K=T_top, T_bot_K=T_bot,
            F_Cp_kW_per_K=s.F_Cp_kW_per_K,
            Q_kW=s.F_Cp_kW_per_K * (T_top - T_bot),
            h_W_m2K=s.h_W_m2K,
        )]

    T_phase = s.T_phase_K if s.T_phase_K is not None else 0.5 * (T_top + T_bot)
    # 範囲内に clamp
    T_phase = max(T_bot, min(T_top, T_phase))

    segs: List[_Segment] = []

    # 顕熱区間 (T_phase より上)
    if T_top > T_phase + 1e-9:
        segs.append(_Segment(
            name=s.name + "_sens_above", is_hot=s.is_hot,
            T_top_K=T_top, T_bot_K=T_phase,
            F_Cp_kW_per_K=s.F_Cp_kW_per_K,
            Q_kW=s.F_Cp_kW_per_K * (T_top - T_phase),
            h_W_m2K=s.h_W_m2K,
        ))

    # 潜熱区間（水平）
    segs.append(_Segment(
        name=s.name + "_latent", is_hot=s.is_hot,
        T_top_K=T_phase, T_bot_K=T_phase,
        F_Cp_kW_per_K=0.0,
        Q_kW=s.Q_latent_kW,
        h_W_m2K=s.h_W_m2K,
    ))

    # 顕熱区間 (T_phase より下)
    if T_phase > T_bot + 1e-9:
        segs.append(_Segment(
            name=s.name + "_sens_below", is_hot=s.is_hot,
            T_top_K=T_phase, T_bot_K=T_bot,
            F_Cp_kW_per_K=s.F_Cp_kW_per_K,
            Q_kW=s.F_Cp_kW_per_K * (T_phase - T_bot),
            h_W_m2K=s.h_W_m2K,
        ))

    return segs


# ===========================================================================
# 複合線（Composite Curve）
# ===========================================================================

def _build_composite(segments: List[_Segment], hot: bool) -> List[Tuple[float, float]]:
    """与熱 or 受熱の複合線を構築する。

    返り値は (H_kW, T_K) の昇順折れ線。H は累積エンタルピー（0 から開始）。

    潜熱区間（T_top==T_bot）は同じ H で温度ジャンプではなく、エンタルピー側で
    水平に伸びる（つまり (H1, T_phase) → (H1+Q_latent, T_phase)）。
    複合線では「同じ温度で熱量が増える」=「水平区間」になる。
    """
    # 区間に存在する温度をすべて集めてソート（昇順）
    Ts = set()
    for seg in segments:
        Ts.add(seg.T_top_K)
        Ts.add(seg.T_bot_K)
    T_list = sorted(Ts)
    if len(T_list) < 2:
        return []

    H_kW = 0.0
    points: List[Tuple[float, float]] = [(0.0, T_list[0])]

    for i in range(len(T_list) - 1):
        T_low  = T_list[i]
        T_high = T_list[i + 1]
        T_mid  = 0.5 * (T_low + T_high)

        # この区間 (T_low, T_high) に存在するストリームの F·Cp を合計
        # 顕熱寄与のみ（潜熱区間は T_top==T_bot なので含まれない）
        FCp_sum = 0.0
        for seg in segments:
            if seg.T_bot_K <= T_low + 1e-9 and seg.T_top_K >= T_high - 1e-9:
                FCp_sum += seg.F_Cp_kW_per_K
        Q_sens = FCp_sum * (T_high - T_low)
        H_kW += Q_sens
        points.append((H_kW, T_high))

        # T_high で潜熱区間が落ちるなら、その瞬間に水平に伸ばす
        # （T_high が、ある潜熱セグメントの T_phase に一致する場合）
        for seg in segments:
            if seg.F_Cp_kW_per_K == 0.0 and abs(seg.T_top_K - T_high) < 1e-9:
                # 潜熱を加算（水平区間: 同じ T で H だけ伸びる）
                H_kW += seg.Q_kW
                points.append((H_kW, T_high))

    return points


# ===========================================================================
# Problem Table Algorithm
# ===========================================================================

def _problem_table(
    hot_segments:  List[_Segment],
    cold_segments: List[_Segment],
    dT_min_K:      float,
) -> Tuple[float, float, float, List[Tuple[float, float]]]:
    """Problem Table Algorithm を実行する。

    Returns
    -------
    Q_H_min_kW, Q_C_min_kW, T_pinch_shifted_K, GCC
        GCC は (T_shifted_K, H_residual_kW) の昇順折れ線。
        T_pinch_shifted は GCC の H_residual=0 となる温度（複数あれば最低温）。

    アルゴリズム:
      1. 全 segment の温度をシフト（hot は -dTmin/2, cold は +dTmin/2）
      2. シフト温度のユニーク値を昇順ソート → 区間境界
      3. 各区間で Σ(F·Cp_cold) - Σ(F·Cp_hot) = ΔCP_net
         区間で潜熱がある場合は Q を直接加算
      4. 各区間の不足熱量 = ΔCP_net × ΔT + Q_latent_cold - Q_latent_hot
      5. 累積（高温から低温へ）
      6. 累積の最大値（最も「不足」が大きい点）が Q_H_min
         そこから 0 にリセットして再積算 → 末尾 = Q_C_min
         0 となる点 = ピンチ温度（シフト）

    潜熱の扱い:
      潜熱は「同じシフト温度で水平に Q が落ちる/上がる」。
      区間境界が潜熱温度に一致する場合に Q を加減算する。
    """
    # シフト温度
    def T_shift(T_K: float, hot: bool) -> float:
        return T_K - 0.5 * dT_min_K if hot else T_K + 0.5 * dT_min_K

    # 全シフト温度を集める
    T_shifts = set()
    for seg in hot_segments:
        T_shifts.add(T_shift(seg.T_top_K, hot=True))
        T_shifts.add(T_shift(seg.T_bot_K, hot=True))
    for seg in cold_segments:
        T_shifts.add(T_shift(seg.T_top_K, hot=False))
        T_shifts.add(T_shift(seg.T_bot_K, hot=False))

    T_list = sorted(T_shifts, reverse=True)   # 高温から低温へ
    if len(T_list) < 2:
        return 0.0, 0.0, T_list[0] if T_list else 0.0, []

    # 各区間 (T_high, T_low) の熱不足量を計算
    # 熱不足量 = (Σ F_Cp_cold - Σ F_Cp_hot) × ΔT + Σ Q_latent_cold - Σ Q_latent_hot
    # 正なら「外部から熱を補う必要がある」、負なら「過剰」
    deficits: List[float] = []      # 各区間の熱不足量 [kW]

    for i in range(len(T_list) - 1):
        T_high = T_list[i]
        T_low  = T_list[i + 1]
        T_mid  = 0.5 * (T_high + T_low)

        FCp_hot_sum = 0.0
        for seg in hot_segments:
            T_top_sh = T_shift(seg.T_top_K, hot=True)
            T_bot_sh = T_shift(seg.T_bot_K, hot=True)
            if T_bot_sh <= T_low + 1e-9 and T_top_sh >= T_high - 1e-9:
                FCp_hot_sum += seg.F_Cp_kW_per_K

        FCp_cold_sum = 0.0
        for seg in cold_segments:
            T_top_sh = T_shift(seg.T_top_K, hot=False)
            T_bot_sh = T_shift(seg.T_bot_K, hot=False)
            if T_bot_sh <= T_low + 1e-9 and T_top_sh >= T_high - 1e-9:
                FCp_cold_sum += seg.F_Cp_kW_per_K

        # 顕熱の不足量（cold が hot より大きいと不足）
        deficit = (FCp_cold_sum - FCp_hot_sum) * (T_high - T_low)

        # 潜熱: T_high にちょうど含まれる潜熱区間を加減算する
        # （T_high で水平区間が起こる）
        for seg in hot_segments:
            if seg.F_Cp_kW_per_K == 0.0 and seg.Q_kW > 0.0:
                T_lat_sh = T_shift(seg.T_top_K, hot=True)
                # この区間の上端 T_high と一致する潜熱を、この区間で消化
                if abs(T_lat_sh - T_high) < 1e-9:
                    deficit -= seg.Q_kW   # hot 潜熱は供給 → 不足が減る
        for seg in cold_segments:
            if seg.F_Cp_kW_per_K == 0.0 and seg.Q_kW > 0.0:
                T_lat_sh = T_shift(seg.T_top_K, hot=False)
                if abs(T_lat_sh - T_high) < 1e-9:
                    deficit += seg.Q_kW   # cold 潜熱は需要 → 不足が増える

        deficits.append(deficit)

    # 累積熱フロー（外部加熱 0 で開始）
    # GCC は「外部加熱 = Q_H_min 投入後の各温度での累積過剰熱量」
    cumsum_zero: List[float] = [0.0]
    for d in deficits:
        cumsum_zero.append(cumsum_zero[-1] + d)

    # 最も大きい正の値 = 「外部加熱なしでは最大このくらい不足する」 = Q_H_min
    Q_H_min = max(cumsum_zero) if cumsum_zero else 0.0
    if Q_H_min < 0:
        Q_H_min = 0.0

    # GCC: 外部加熱 Q_H_min を最高温で投入した後の累積
    # (T_shifted, H_residual)
    # H_residual = Q_H_min - Σ deficits (高温から積算)
    gcc_points: List[Tuple[float, float]] = []
    H = Q_H_min
    gcc_points.append((T_list[0], H))
    for i, d in enumerate(deficits):
        H -= d
        gcc_points.append((T_list[i + 1], H))

    # Q_C_min = GCC の末尾値（最低温で残った余剰熱）
    Q_C_min = gcc_points[-1][1] if gcc_points else 0.0
    if Q_C_min < 0:
        Q_C_min = 0.0

    # ピンチ温度（シフト）: GCC の H_residual が最小 (~0) となる温度
    # ほぼ 0 に近い点を取る
    pinch_idx = min(range(len(gcc_points)), key=lambda k: abs(gcc_points[k][1]))
    T_pinch_shifted = gcc_points[pinch_idx][0]

    return Q_H_min, Q_C_min, T_pinch_shifted, gcc_points


# ===========================================================================
# 総伝熱面積 A_total（Bath 式）
# ===========================================================================

def _calc_A_total(
    composite_hot:  List[Tuple[float, float]],
    composite_cold: List[Tuple[float, float]],
    Q_H_min_kW:     float,
    hot_segments:   List[_Segment],
    cold_segments:  List[_Segment],
    dT_min_K:       float,
) -> float:
    """Bath 式で総伝熱面積を targeting する。

    A_total = Σ_k (1/ΔT_LM,k) × Σ_i (q_i,k / h_i)

    各温度区間 k で、その区間に存在する hot/cold ストリームを集めて計算する。
    複合線を共通の Q 軸で区間分割するのが教科書通りだが、実装簡略化のため
    シフト温度ベースの区間で評価する。

    返り値の単位は m²（kW・K → m²）。
    """
    # シフト温度区間で割る
    def T_shift(T_K: float, hot: bool) -> float:
        return T_K - 0.5 * dT_min_K if hot else T_K + 0.5 * dT_min_K

    T_shifts = set()
    for seg in hot_segments:
        T_shifts.add(T_shift(seg.T_top_K, hot=True))
        T_shifts.add(T_shift(seg.T_bot_K, hot=True))
    for seg in cold_segments:
        T_shifts.add(T_shift(seg.T_top_K, hot=False))
        T_shifts.add(T_shift(seg.T_bot_K, hot=False))

    T_list = sorted(T_shifts, reverse=True)
    if len(T_list) < 2:
        return 0.0

    A_total = 0.0
    for i in range(len(T_list) - 1):
        T_high_sh = T_list[i]
        T_low_sh  = T_list[i + 1]
        # この区間の hot 流体の高低温（実温度）
        T_hot_high  = T_high_sh + 0.5 * dT_min_K
        T_hot_low   = T_low_sh  + 0.5 * dT_min_K
        T_cold_high = T_high_sh - 0.5 * dT_min_K
        T_cold_low  = T_low_sh  - 0.5 * dT_min_K

        # この区間に存在する hot/cold セグメント
        Q_hot_sum_per_h  = 0.0   # Σ q_i / h_i  [kW·m²·K/W]
        Q_cold_sum_per_h = 0.0
        Q_hot_total      = 0.0
        Q_cold_total     = 0.0

        for seg in hot_segments:
            T_top_sh = T_shift(seg.T_top_K, hot=True)
            T_bot_sh = T_shift(seg.T_bot_K, hot=True)
            if T_bot_sh <= T_low_sh + 1e-9 and T_top_sh >= T_high_sh - 1e-9:
                if seg.F_Cp_kW_per_K > 0:
                    q_i = seg.F_Cp_kW_per_K * (T_high_sh - T_low_sh)
                else:
                    # 潜熱は区間 width=0 なので、上端ぴったりの場合のみ全量
                    q_i = seg.Q_kW if abs(T_top_sh - T_high_sh) < 1e-9 else 0.0
                Q_hot_total += q_i
                # h_i は W/m²/K → kW/m²/K に変換 (1/1000)
                h_kW = seg.h_W_m2K / 1000.0
                if h_kW > 0:
                    Q_hot_sum_per_h += q_i / h_kW

        for seg in cold_segments:
            T_top_sh = T_shift(seg.T_top_K, hot=False)
            T_bot_sh = T_shift(seg.T_bot_K, hot=False)
            if T_bot_sh <= T_low_sh + 1e-9 and T_top_sh >= T_high_sh - 1e-9:
                if seg.F_Cp_kW_per_K > 0:
                    q_i = seg.F_Cp_kW_per_K * (T_high_sh - T_low_sh)
                else:
                    q_i = seg.Q_kW if abs(T_top_sh - T_high_sh) < 1e-9 else 0.0
                Q_cold_total += q_i
                h_kW = seg.h_W_m2K / 1000.0
                if h_kW > 0:
                    Q_cold_sum_per_h += q_i / h_kW

        # この区間で熱交換可能な量（小さい方）
        Q_exch = min(Q_hot_total, Q_cold_total)
        if Q_exch <= 1e-9:
            continue

        # 区間内の温度差 ΔT_LM
        # hot 側: T_hot_high → T_hot_low
        # cold 側: T_cold_high → T_cold_low
        # 向流: ΔT1 = T_hot_high - T_cold_high, ΔT2 = T_hot_low - T_cold_low
        dT1 = T_hot_high - T_cold_high
        dT2 = T_hot_low  - T_cold_low
        if dT1 <= 0 or dT2 <= 0:
            continue
        if abs(dT1 - dT2) < 1e-9:
            dT_LM = dT1
        else:
            dT_LM = (dT1 - dT2) / math.log(dT1 / dT2)

        # 区間内の比例配分: q_exch_i = q_i × (Q_exch / Q_hot_total) など
        # ただし、上の Q_hot_sum_per_h は Q_hot_total ベースなので、Q_exch にスケール
        if Q_hot_total > 0 and Q_cold_total > 0:
            scale_hot  = Q_exch / Q_hot_total
            scale_cold = Q_exch / Q_cold_total
            sum_per_h_scaled = Q_hot_sum_per_h * scale_hot + Q_cold_sum_per_h * scale_cold
            A_total += sum_per_h_scaled / dT_LM

    return A_total


# ===========================================================================
# Linnhoff 式: 最少熱交換器数
# ===========================================================================

def _calc_N_HE_min(
    n_hot: int, n_cold: int,
    n_heating_utilities: int, n_cooling_utilities: int,
    has_pinch: bool,
) -> int:
    """Linnhoff 式: N_HE,min = N_streams + N_utilities - 1
    ピンチで分割する場合は (上半分の式) + (下半分の式) で計算する。

    has_pinch=False (閾値問題): 全体で 1 つの式
    has_pinch=True : ピンチ上下で分割
    """
    if not has_pinch:
        # 全ストリーム + ユーティリティ - 1
        return max(0, n_hot + n_cold + n_heating_utilities + n_cooling_utilities - 1)
    # 上半分: 全 hot + 加熱utility + ピンチを跨ぐcold（簡易: 全cold 含めて近似）
    # 教科書 4-7 式: 厳密には上下それぞれの「アクティブな」ストリーム数を数える。
    # 簡易実装として n_hot + n_cold + utilities - 2 (上下分で 2 度引く)
    return max(0, n_hot + n_cold + n_heating_utilities + n_cooling_utilities - 2)


# ===========================================================================
# tier 別 Q 配分
# ===========================================================================

def _interp_gcc(gcc: List[Tuple[float, float]], T_shifted: float) -> float:
    """GCC を T_shifted で線形補間して H_residual を返す。

    GCC は (T_shifted, H_residual) のタプルリスト。任意の順序を許容。
    範囲外は端の値で clamp。
    """
    if not gcc:
        return 0.0
    pts = sorted(gcc, key=lambda p: p[0])
    if T_shifted <= pts[0][0]:
        return pts[0][1]
    if T_shifted >= pts[-1][0]:
        return pts[-1][1]
    for i in range(len(pts) - 1):
        T1, H1 = pts[i]
        T2, H2 = pts[i + 1]
        if T1 <= T_shifted <= T2:
            if abs(T2 - T1) < 1e-12:
                return H1
            return H1 + (H2 - H1) * (T_shifted - T1) / (T2 - T1)
    return pts[-1][1]


def _assign_utilities_to_tiers(
    gcc:           List[Tuple[float, float]],
    Q_H_min_kW:    float,
    Q_C_min_kW:    float,
    T_pinch_shifted_K: float,
    heating_tiers: List[UtilityTier],
    cooling_tiers: List[UtilityTier],
    dT_min_K:      float,
) -> Dict[str, float]:
    """GCC を tier 温度で水平に切って多 tier 配分する。

    Linnhoff 標準手法 (教科書 図 4.27 の応用):
      加熱 (ピンチより上):
        各 heating tier の effective_T_shifted = supply_T_K - dT_min/2
        (このシフト温度より下の領域なら、その tier で加熱可能)
        低品位 (低 supply_T = 安価) な tier から順に「使えるだけ使う」
        各区間で消化される Q = GCC の H 差分

      冷却 (ピンチより下):
        対称的に高品位 (高 supply_T = 安価) な tier から順に充填

    Returns
    -------
    {tier_name: Q_kW} の dict。1 つの tier だけが使われる場合もあれば、
    複数 tier が組み合わせられる場合もある。
    """
    breakdown: Dict[str, float] = {}

    if not gcc:
        return breakdown

    def _normalize(bd: Dict[str, float], target: float,
                   tiers_in_order: List[UtilityTier]):
        """配分合計が target になるよう調整する。

        bd は加熱専用 or 冷却専用の dict (混在させない)。
        tiers_in_order は **削減 priority 順** (= 安価順) のリスト。
        - 合計 > target: tiers_in_order の先頭から削る
                       (pocket 部分は内部熱交換で消化される解釈)
        - 合計 < target: tiers_in_order の末尾 (= 最後のフォールバック) に追加
        """
        total = sum(bd.values())
        over = total - target
        if over > 1e-6:
            for tier in tiers_in_order:
                if over <= 1e-9:
                    break
                if tier.name not in bd:
                    continue
                reduce = min(over, bd[tier.name])
                bd[tier.name] -= reduce
                over -= reduce
            for k in list(bd.keys()):
                if bd[k] < 1e-9:
                    del bd[k]
        elif over < -1e-6:
            fallback = tiers_in_order[-1]
            bd[fallback.name] = bd.get(fallback.name, 0.0) + (-over)

    # ===== 加熱側 =====
    # 各 tier の T_eff_shifted (= supply_T - dT_min/2) で GCC を切り、
    # 低 supply_T (=安価) 順に「累積差分」を加算する。
    # GCC pocket (非単調) で過剰計上が起きるが、最後に _normalize で
    # 合計を Q_H_min に揃える (pocket 部分は内部熱交換で消化扱い)。
    if Q_H_min_kW > 1e-9 and heating_tiers:
        heating_bd: Dict[str, float] = {}
        # 安価順 = supply_T 昇順 (LP < MP < HP < 燃料)
        tiers_sorted = sorted(heating_tiers, key=lambda t: t.supply_T_K)
        total_assigned = 0.0
        for tier in tiers_sorted:
            T_eff_shifted = tier.supply_T_K - 0.5 * dT_min_K
            if T_eff_shifted <= T_pinch_shifted_K + 1e-9:
                continue
            H_at_cut = _interp_gcc(gcc, T_eff_shifted)
            Q_assigned = max(0.0, H_at_cut - total_assigned)
            if Q_assigned > 1e-9:
                heating_bd[tier.name] = heating_bd.get(tier.name, 0.0) + Q_assigned
                total_assigned = max(total_assigned, H_at_cut)
        _normalize(heating_bd, Q_H_min_kW, tiers_sorted)
        breakdown.update(heating_bd)

    # ===== 冷却側 =====
    # 対称的: 高 supply_T (=安価) 順 (空冷 → 冷却水 → プロピレン → エチレン)
    if Q_C_min_kW > 1e-9 and cooling_tiers:
        cooling_bd: Dict[str, float] = {}
        # 安価順 = supply_T 降順 (空冷 308.15 K → エチレン 173.15 K)
        tiers_sorted = sorted(cooling_tiers, key=lambda t: -t.supply_T_K)
        total_assigned = 0.0
        for tier in tiers_sorted:
            T_eff_shifted = tier.supply_T_K + 0.5 * dT_min_K
            if T_eff_shifted >= T_pinch_shifted_K - 1e-9:
                continue
            H_at_cut = _interp_gcc(gcc, T_eff_shifted)
            Q_assigned = max(0.0, H_at_cut - total_assigned)
            if Q_assigned > 1e-9:
                cooling_bd[tier.name] = cooling_bd.get(tier.name, 0.0) + Q_assigned
                total_assigned = max(total_assigned, H_at_cut)
        _normalize(cooling_bd, Q_C_min_kW, tiers_sorted)
        breakdown.update(cooling_bd)

    return breakdown


# ===========================================================================
# 公開 API
# ===========================================================================

def pinch_analysis(
    streams:        List[HIStream],
    dT_min_K:       float = 10.0,
    heating_tiers:  Optional[List[UtilityTier]] = None,
    cooling_tiers:  Optional[List[UtilityTier]] = None,
) -> HIResult:
    """ピンチ解析を実行し、targeting 結果を返す。

    Parameters
    ----------
    streams : List[HIStream]
        対象ストリームのリスト（hot/cold は T_in > T_out で自動判定）
    dT_min_K : float
        最小接近温度差 [K]（教科書標準は 10）
    heating_tiers, cooling_tiers : List[UtilityTier] | None
        ユーティリティ tier。None なら utility_breakdown は空 dict になる。
        BO 内では src.utility_selector の tier を変換して渡す想定。

    Returns
    -------
    HIResult
        Q_H_min, Q_C_min, ピンチ温度, A_total, N_HE_min, utility_breakdown,
        composite curves, GCC を含む。
    """
    if not streams:
        return HIResult(
            Q_H_min_kW=0.0, Q_C_min_kW=0.0,
            T_pinch_hot_K=0.0, T_pinch_cold_K=0.0,
            A_total_m2=0.0, N_HE_min=0,
            utility_breakdown={},
            composite_hot=[], composite_cold=[], GCC=[],
            feasible=True, message="no streams",
        )

    # 1. ストリームを区間に展開
    all_segments = []
    for s in streams:
        all_segments.extend(_expand_stream(s))

    hot_segments  = [seg for seg in all_segments if seg.is_hot]
    cold_segments = [seg for seg in all_segments if not seg.is_hot]

    # 2. 複合線
    composite_hot  = _build_composite(hot_segments,  hot=True)
    composite_cold = _build_composite(cold_segments, hot=False)

    # 3. Problem Table Algorithm
    Q_H_min, Q_C_min, T_pinch_shifted, gcc = _problem_table(
        hot_segments, cold_segments, dT_min_K,
    )

    # 4. ピンチ温度（実温度）
    T_pinch_hot_K  = T_pinch_shifted + 0.5 * dT_min_K
    T_pinch_cold_K = T_pinch_shifted - 0.5 * dT_min_K

    # 5. Bath 式
    A_total = _calc_A_total(
        composite_hot, composite_cold,
        Q_H_min, hot_segments, cold_segments, dT_min_K,
    )

    # 6. tier 別 Q 配分
    breakdown: Dict[str, float] = {}
    if heating_tiers is not None and cooling_tiers is not None:
        breakdown = _assign_utilities_to_tiers(
            gcc, Q_H_min, Q_C_min, T_pinch_shifted,
            heating_tiers, cooling_tiers, dT_min_K,
        )

    # 7. Linnhoff 式: ピンチが存在するか判定
    has_pinch = (Q_H_min > 1e-9 and Q_C_min > 1e-9)
    n_hot  = len(set(s.name for s in streams if s.is_hot))
    n_cold = len(set(s.name for s in streams if not s.is_hot))
    n_hu   = 1 if Q_H_min > 1e-9 else 0
    n_cu   = 1 if Q_C_min > 1e-9 else 0
    N_HE   = _calc_N_HE_min(n_hot, n_cold, n_hu, n_cu, has_pinch)

    # 8. feasible 判定（簡易）
    feasible = True
    msg = ""
    if heating_tiers is not None and Q_H_min > 1e-9:
        # 最高温 tier でも届かない場合は infeasible
        T_cold_max = max(p[0] for p in gcc) + 0.5 * dT_min_K
        if not any(t.supply_T_K - dT_min_K >= T_cold_max for t in heating_tiers):
            feasible = False
            msg = f"加熱: 最高温 tier でも T_cold_max={T_cold_max-273.15:.1f}°C に届かない"
    if cooling_tiers is not None and Q_C_min > 1e-9:
        # 最低温 tier でも届かない場合は infeasible
        # (= 全 hot 最低温度より、最低 supply_T_K + dT_min が高い)
        T_hot_min = min(p[0] for p in gcc) - 0.5 * dT_min_K
        if not any(t.supply_T_K + dT_min_K <= T_hot_min for t in cooling_tiers):
            feasible = False
            extra = (f"冷却: 最低温 tier でも T_hot_min="
                     f"{T_hot_min-273.15:.1f}°C に届かない")
            msg = msg + " | " + extra if msg else extra

    return HIResult(
        Q_H_min_kW=Q_H_min,
        Q_C_min_kW=Q_C_min,
        T_pinch_hot_K=T_pinch_hot_K,
        T_pinch_cold_K=T_pinch_cold_K,
        A_total_m2=A_total,
        N_HE_min=N_HE,
        utility_breakdown=breakdown,
        composite_hot=composite_hot,
        composite_cold=composite_cold,
        GCC=gcc,
        feasible=feasible,
        message=msg,
    )


# ===========================================================================
# 可視化補助（matplotlib 依存。呼び出し側でインポートを確認）
# ===========================================================================

def plot_TQ(result: HIResult, ax=None, T_in_celsius: bool = True):
    """T-Q 線図を描画する。

    Parameters
    ----------
    result : HIResult
    ax : matplotlib.axes.Axes | None
        None なら現在の軸に描画
    T_in_celsius : bool
        True なら縦軸を °C で表示
    """
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()

    def Tdisp(T_K: float) -> float:
        return T_K - 273.15 if T_in_celsius else T_K

    if result.composite_hot:
        H_hot = [p[0] for p in result.composite_hot]
        T_hot = [Tdisp(p[1]) for p in result.composite_hot]
        ax.plot(H_hot, T_hot, 'r-o', label='与熱複合線 (Hot)', linewidth=2)

    if result.composite_cold:
        # 受熱複合線は Q_C_min だけ右にシフトして描く（教科書の作図法）
        offset = result.Q_C_min_kW
        H_cold = [p[0] + offset for p in result.composite_cold]
        T_cold = [Tdisp(p[1]) for p in result.composite_cold]
        ax.plot(H_cold, T_cold, 'b-o', label='受熱複合線 (Cold)', linewidth=2)

    ax.set_xlabel('Q [kW]')
    ax.set_ylabel('T [°C]' if T_in_celsius else 'T [K]')
    ax.set_title(
        f'T-Q 線図  (Q_H_min={result.Q_H_min_kW:.1f}, '
        f'Q_C_min={result.Q_C_min_kW:.1f} kW, '
        f'pinch={Tdisp(result.T_pinch_hot_K):.1f} °C)'
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


# ===========================================================================
# PDH プロセス用 ストリーム抽出
# ===========================================================================

def _make_sensible_stream(
    name: str, T_in_K: float, T_out_K: float, Q_kW: float,
    h_W_m2K: float = H_LIQUID_LIQUID_W_M2K,
    phase:   str   = StreamPhase.LIQUID,
) -> Optional[HIStream]:
    """温度差を持つ顕熱ストリームを HIStream に変換する。

    Q_kW が小さい (< 1e-6) または T_in == T_out の場合は None を返す。
    phase は contest §4-4 U 表索引用 (StreamPhase.GAS or LIQUID)。
    """
    if abs(Q_kW) < 1e-6:
        return None
    dT = abs(T_in_K - T_out_K)
    if dT < 1e-6:
        return None
    return HIStream(
        name=name, T_in_K=T_in_K, T_out_K=T_out_K,
        F_Cp_kW_per_K=abs(Q_kW) / dT,
        h_W_m2K=h_W_m2K,
        phase=phase,
    )


def _make_latent_stream(
    name: str, T_K: float, Q_kW: float, is_hot: bool,
    h_W_m2K: float = H_CONDENSER_W_M2K,
) -> Optional[HIStream]:
    """潜熱のみのストリーム (T_in == T_out) を構築。

    is_hot=True (= 冷却される hot) は CONDENSING (ガス凝縮)、
    is_hot=False (= 加熱される cold) は EVAPORATING (液蒸発)。
    """
    if Q_kW < 1e-6:
        return None
    eps = 1e-3
    if is_hot:
        T_in_K, T_out_K = T_K + eps, T_K - eps
        phase = StreamPhase.CONDENSING
    else:
        T_in_K, T_out_K = T_K - eps, T_K + eps
        phase = StreamPhase.EVAPORATING
    return HIStream(
        name=name, T_in_K=T_in_K, T_out_K=T_out_K,
        F_Cp_kW_per_K=0.0,
        Q_latent_kW=Q_kW,
        T_phase_K=T_K,
        h_W_m2K=h_W_m2K,
        phase=phase,
    )


def _make_sens_plus_latent_stream(
    name: str, T_in_K: float, T_out_K: float,
    Q_sensible_kW: float, Q_latent_kW: float,
    T_phase_K: Optional[float] = None,
    h_W_m2K: float = H_VAPORIZER_W_M2K,
    phase:   Optional[str] = None,
) -> Optional[HIStream]:
    """顕熱+潜熱を持つストリーム (例: 液→気の気化器) を構築。

    phase: 顕熱区間の相分類 (None なら自動: hot=GAS, cold=LIQUID)。
    潜熱区間は別途 _make_latent_stream で出すべきだが、簡略化のため
    本関数は「顕熱+潜熱合算」した HIStream として返す。
    targeting Bath 式では潜熱と顕熱が分けて扱われるが、phase は顕熱側に当てる。
    """
    Q_sens_abs = abs(Q_sensible_kW)
    Q_lat_abs  = abs(Q_latent_kW)
    if Q_sens_abs < 1e-6 and Q_lat_abs < 1e-6:
        return None
    dT = abs(T_in_K - T_out_K)
    F_Cp = Q_sens_abs / dT if dT > 1e-6 and Q_sens_abs > 1e-6 else 0.0
    if phase is None:
        is_hot = T_in_K > T_out_K
        phase  = StreamPhase.GAS if is_hot else StreamPhase.LIQUID
    return HIStream(
        name=name, T_in_K=T_in_K, T_out_K=T_out_K,
        F_Cp_kW_per_K=F_Cp,
        Q_latent_kW=Q_lat_abs,
        T_phase_K=T_phase_K,
        h_W_m2K=h_W_m2K,
        phase=phase,
    )


def extract_streams(
    one_pass:    Dict,
    swing_T_in:  float,
) -> List[HIStream]:
    """run_one_pass の戻り値から HI 対象ストリームを抽出する。

    設計判断 (2026-05-08):
      各ユニットの T_in/T_out/Q は run_one_pass 結果と equipment dataclass から
      直接読み取り、新たな仮定は導入しない。Q が小さいストリームは自動的に除外。
      PSA preheat は符号で hot/cold を判定。

    Parameters
    ----------
    one_pass : Dict
        run_one_pass の戻り値 (dict)
    swing_T_in : float
        反応器入口温度 [K] (= design.swing.T_in、Q_preheat の出口側温度)

    Returns
    -------
    List[HIStream]
        ストリーム名は H1/H2/...(hot) または C1/C2/...(cold) 慣例。
    """
    R = one_pass
    streams: List[HIStream] = []

    # ---- H1: Cooler (反応器出口冷却、ガス顕熱) ----
    cooled_eq  = R['cooled'].equipment
    s = _make_sensible_stream(
        'H1_rx_cooler',
        T_in_K=R['rx_out'].T_in,
        T_out_K=R['cooled'].outlet.T_in,
        Q_kW=cooled_eq.Q_sensible_kW,
        h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
        phase=StreamPhase.GAS,
    )
    if s: streams.append(s)
    # 潜熱がある場合 (現状の cooler は phase_change=False で 0)
    if cooled_eq.Q_latent_kW > 1e-6:
        s = _make_latent_stream(
            'H1_rx_cooler_latent',
            T_K=0.5*(R['rx_out'].T_in + R['cooled'].outlet.T_in),
            Q_kW=cooled_eq.Q_latent_kW, is_hot=True,
        )
        if s: streams.append(s)

    # ---- H2: Intercool (圧縮ガス顕熱) ----
    ic_eq = R['intercool'].equipment
    s = _make_sensible_stream(
        'H2_intercool',
        T_in_K=R['comp2a'].outlet.T_in,
        T_out_K=R['intercool'].outlet.T_in,
        Q_kW=ic_eq.Q_sensible_kW, h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
        phase=StreamPhase.GAS,
    )
    if s: streams.append(s)

    # ---- H3: Desuper (Comp2b → Dist2、ガス顕熱) ----
    ds_eq = R['desuper'].equipment
    s = _make_sensible_stream(
        'H3_desuper',
        T_in_K=R['comp2b'].outlet.T_in,
        T_out_K=R['desuper'].outlet.T_in,
        Q_kW=ds_eq.Q_sensible_kW, h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
        phase=StreamPhase.GAS,
    )
    if s: streams.append(s)

    # ---- C1: mem_precool (液→気、顕熱+潜熱) ----
    # 設計判断 (2026-05-09): 顕熱区間は LIQUID、潜熱区間は EVAPORATING で別 stream に分離。
    mp_eq = R['mem_precool'].equipment
    if mp_eq.Q_sensible_kW > 1e-6:
        s = _make_sensible_stream(
            'C1_mem_precool_sens',
            T_in_K=R['r2'].bottom.T_in,
            T_out_K=R['mem_precool'].outlet.T_in,
            Q_kW=mp_eq.Q_sensible_kW,
            h_W_m2K=H_LIQUID_LIQUID_W_M2K,
            phase=StreamPhase.LIQUID,
        )
        if s: streams.append(s)
    if mp_eq.Q_latent_kW > 1e-6:
        s = _make_latent_stream(
            'C1_mem_precool_lat',
            T_K=R['mem_precool'].outlet.T_in,
            Q_kW=mp_eq.Q_latent_kW, is_hot=False,    # → EVAPORATING
            h_W_m2K=H_VAPORIZER_W_M2K,
        )
        if s: streams.append(s)

    # ---- C2: Reactor preheat (ガス顕熱、38°C → 950K の最大 cold) ----
    Q_preheat_GJh = R['r_rx'].effluent.Q_preheat
    Q_preheat_kW  = Q_preheat_GJh * 1e9 / 3600.0 / 1000.0
    T_feed = R['reactor_inlet'].T_in
    s = _make_sensible_stream(
        'C2_rx_preheat',
        T_in_K=T_feed, T_out_K=swing_T_in, Q_kW=Q_preheat_kW,
        h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
        phase=StreamPhase.GAS,
    )
    if s: streams.append(s)

    # ---- Dist1/2/3: condenser (潜熱), reboiler (潜熱), feed preheat (顕熱) ----
    for label, key in [('Dist1', 'r1'), ('Dist2', 'r2'), ('Dist3', 'r3')]:
        eq = R[key].equipment
        T_top = R[key].top.T_in
        T_bot = R[key].bottom.T_in
        # condenser (与熱、潜熱)
        s = _make_latent_stream(f'H_{label}_cond', T_K=T_top, Q_kW=eq.Q_cond, is_hot=True)
        if s: streams.append(s)
        # reboiler (受熱、潜熱)
        s = _make_latent_stream(f'C_{label}_reb', T_K=T_bot, Q_kW=eq.Q_reb, is_hot=False)
        if s: streams.append(s)
        # feed preheat (受熱、顕熱: feed→T_bot)
        # フィード入口温度・相は塔ごとに異なる:
        #   Dist1: pump1 outlet → 液 (高圧液)
        #   Dist2: desuper outlet (50°C) → ガス (圧縮蒸気)
        #   Dist3: r_mem product → 液 (膜保留側、加圧液)
        if label == 'Dist1':
            T_feed_in    = R['pump1'].outlet.T_in
            feed_phase   = StreamPhase.LIQUID
            feed_h       = H_LIQUID_LIQUID_W_M2K
        elif label == 'Dist2':
            T_feed_in    = R['desuper'].outlet.T_in
            feed_phase   = StreamPhase.GAS
            feed_h       = H_HIGH_TEMP_GAS_W_M2K
        else:  # Dist3
            T_feed_in    = R['r_mem'].product.T_out
            feed_phase   = StreamPhase.LIQUID
            feed_h       = H_LIQUID_LIQUID_W_M2K
        s = _make_sensible_stream(
            f'C_{label}_feed_preheat',
            T_in_K=T_feed_in, T_out_K=T_bot, Q_kW=eq.Q_feed_preheat_kW,
            h_W_m2K=feed_h,
            phase=feed_phase,
        )
        if s: streams.append(s)

    # ---- Mem vaporizer (受熱) と cooler (与熱) ----
    meq = R['r_mem'].equipment
    # Mem 気化器 (vapor_feed=True なら Q_vap_kW=0、自動的にスキップされる)
    if meq.Q_vap_kW > 1e-6 and not _is_nan(meq.T_vap_in_K):
        s = _make_sens_plus_latent_stream(
            'C_mem_vaporizer',
            T_in_K=meq.T_vap_in_K, T_out_K=meq.T_vap_out_K,
            Q_sensible_kW=0.0, Q_latent_kW=meq.Q_vap_kW,  # 全て潜熱として扱う (簡易)
            T_phase_K=0.5*(meq.T_vap_in_K + meq.T_vap_out_K),
            h_W_m2K=H_VAPORIZER_W_M2K,
        )
        if s: streams.append(s)
    # Mem 製品冷却器 (与熱): 顕熱 (T_cond_in → T_cond_out のガス冷却) と
    # 潜熱 (T_cond_out で凝縮) を別ストリームに分割
    if not _is_nan(meq.T_cond_in_K) and not _is_nan(meq.T_cond_out_K):
        # 顕熱区間 (ガス冷却)
        if meq.Q_cond_sensible_kW > 1e-6:
            s = _make_sensible_stream(
                'H_mem_cooler_gas',
                T_in_K=meq.T_cond_in_K, T_out_K=meq.T_cond_out_K,
                Q_kW=meq.Q_cond_sensible_kW, h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
                phase=StreamPhase.GAS,
            )
            if s: streams.append(s)
        # 潜熱区間 (T_cond_out での凝縮)
        if meq.Q_cond_latent_kW > 1e-6:
            s = _make_latent_stream(
                'H_mem_cooler_cond',
                T_K=meq.T_cond_out_K,
                Q_kW=meq.Q_cond_latent_kW, is_hot=True,
                h_W_m2K=H_CONDENSER_W_M2K,
            )
            if s: streams.append(s)

    # ---- PSA preheat (符号で判定) ----
    peq = R['r_psa'].equipment
    Q_psa = peq.Q_preheat_kW
    if abs(Q_psa) > 1e-3:
        T_psa_in  = R['r2'].top.T_in
        T_psa_out = 298.15        # PSAFixedParams.T_abs デフォルト
        if Q_psa > 0:
            # 加熱 (受熱、cold)
            s = _make_sensible_stream(
                'C_psa_preheat', T_in_K=T_psa_in, T_out_K=T_psa_out,
                Q_kW=Q_psa, h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
                phase=StreamPhase.GAS,
            )
        else:
            # 冷却 (与熱、hot): T_in > T_out
            s = _make_sensible_stream(
                'H_psa_precool', T_in_K=T_psa_in, T_out_K=T_psa_out,
                Q_kW=abs(Q_psa), h_W_m2K=H_HIGH_TEMP_GAS_W_M2K,
                phase=StreamPhase.GAS,
            )
        if s: streams.append(s)

    return streams


def _is_nan(x: float) -> bool:
    return x != x   # NaN != NaN


# ===========================================================================
# プロジェクト標準 tier の取得 (src.utility_selector から変換)
# ===========================================================================

# ===========================================================================
# OPEX 計算ヘルパ (BO 評価関数・比較レポートで再利用)
# ===========================================================================

# economics.py の opex キー → 熱の種別 マッピング
#
# 設計判断 (2026-05-09): economics.py が utility 名を動的に埋め込むスタイル
# (例: 'Dist1リボイラ (LP Steam≒)') に変わったので、固定 key 表は使えない。
# prefix ベースで判定する (より長く・特異な prefix を先に置く)。
#
# economics.py の opex キー命名が変わったら本テーブルを更新すること。
#
# 規約: 順序が意味を持つ (前から走査して先に match した prefix を採用)。
# 例: 'PSA予熱冷却 (...)' は 'PSA予熱' より先に書いて cold と判定させる。
_OPEX_PREFIX_HEAT_KIND_ORDERED: List[Tuple[str, str]] = [
    # PSA予熱: 熱方向に応じて分岐 (具体度の高い順に並べる)
    ('PSA予熱冷却',         'cold'),   # T_in > 25°C のとき (現状はほぼ発生しない)
    ('PSA予熱',             'hot'),    # T_in < 25°C のとき (Dist2 partial cond で発生)
    # 蒸留塔 リボイラ・フィード予熱・コンデンサ
    ('Dist1リボイラ',       'hot'),
    ('Dist2リボイラ',       'hot'),
    ('Dist3リボイラ',       'hot'),
    ('Dist1フィード予熱',   'hot'),
    ('Dist2フィード予熱',   'hot'),
    ('Dist3フィード予熱',   'hot'),
    ('Dist1コンデンサ',     'cold'),
    ('Dist2コンデンサ',     'cold'),
    ('Dist3コンデンサ',     'cold'),
    # 個別熱交換器 (utility_selector で tier 名動的)
    ('Cooler',              'cold'),   # 反応器出口冷却
    ('Intercool',           'cold'),   # 圧縮段間冷却
    ('Desuper',             'cold'),   # Comp2b → Dist2 desuperheater
    ('MemPrecool',          'hot'),    # 液→気化 (加熱)
]

# 完全一致が必要な fixed key (utility 名を含まない単純 key)
_OPEX_FIXED_HEAT_KEY_KIND: Dict[str, str] = {
    'Mem気化器蒸気':       'hot',
    'Reactor予熱燃料':     'hot',
    'Mem冷却器冷水':       'cold',
}


def classify_heat_opex_key(key: str) -> Optional[str]:
    """economics.opex のキーから熱の種別 ('hot'/'cold') を返す。

    熱関連でない (電力・触媒・吸着剤・原料費等) は None。
    """
    # prefix を順に試す (順序が意味を持つ)
    for prefix, kind in _OPEX_PREFIX_HEAT_KIND_ORDERED:
        if key.startswith(prefix):
            return kind
    return _OPEX_FIXED_HEAT_KEY_KIND.get(key)


def summarize_current_heat_opex(opex: Dict[str, float]) -> Dict[str, object]:
    """economics.opex から熱関連項目を集計する。

    Parameters
    ----------
    opex : Dict[str, float]
        Economics.opex (各項目: 億円/年)

    Returns
    -------
    Dict with keys:
        'hot_total':       加熱系 OPEX 合計 [億円/年]
        'cold_total':      冷却系 OPEX 合計 [億円/年]
        'total':           hot + cold [億円/年]
        'hot_breakdown':   {key: 億円/年} hot 内訳
        'cold_breakdown':  {key: 億円/年} cold 内訳
    """
    hot_total = 0.0
    cold_total = 0.0
    hot_breakdown: Dict[str, float] = {}
    cold_breakdown: Dict[str, float] = {}
    for k, v in opex.items():
        kind = classify_heat_opex_key(k)
        if kind == 'hot':
            hot_total += v
            hot_breakdown[k] = v
        elif kind == 'cold':
            cold_total += v
            cold_breakdown[k] = v
    return {
        'hot_total':      hot_total,
        'cold_total':     cold_total,
        'total':          hot_total + cold_total,
        'hot_breakdown':  hot_breakdown,
        'cold_breakdown': cold_breakdown,
    }


def calc_hi_opex_okuyen(
    hi_result:               HIResult,
    heating_tiers:           List[UtilityTier],
    cooling_tiers:           List[UtilityTier],
    operating_hours_per_year: float = 8000.0,
) -> Dict[str, object]:
    """HIResult.utility_breakdown から tier 別 OPEX [億円/年] を計算する。

    BO 評価関数で「HI 後の熱関連 OPEX」を即時計算するための関数。
    変換式: Q [kW] × 3.6e-3 [GJ/(kW·h)] × hours × 円/GJ / 1e8 → 億円/年

    Parameters
    ----------
    hi_result : HIResult
        pinch_analysis の戻り値
    heating_tiers, cooling_tiers : List[UtilityTier]
        tier 単価情報 (get_default_utility_tiers() で取得可能)
    operating_hours_per_year : float
        年間稼働時間 [h]、デフォルト 8000 (= プロジェクト標準)

    Returns
    -------
    Dict with keys:
        'hot_total':      加熱 OPEX 合計 [億円/年]
        'cold_total':     冷却 OPEX 合計 [億円/年]
        'total':          総熱関連 OPEX [億円/年]
        'breakdown':      {tier_name: 億円/年}
        'unmatched':      {tier_name: Q_kW}  (tier リストに該当なし)
    """
    all_tiers = {t.name: t for t in heating_tiers + cooling_tiers}
    hot_total  = 0.0
    cold_total = 0.0
    breakdown: Dict[str, float] = {}
    unmatched: Dict[str, float] = {}

    for tier_name, Q_kW in hi_result.utility_breakdown.items():
        tier = all_tiers.get(tier_name)
        if tier is None:
            unmatched[tier_name] = Q_kW
            continue
        # Q [kW] × 3.6e-3 GJ/kWh × hours × 円/GJ / 1e8 → 億円/年
        cost = Q_kW * 3.6e-3 * operating_hours_per_year * tier.jpy_per_GJ / 1.0e8
        breakdown[tier_name] = cost
        if tier.is_heating:
            hot_total += cost
        else:
            cold_total += cost

    return {
        'hot_total':  hot_total,
        'cold_total': cold_total,
        'total':      hot_total + cold_total,
        'breakdown':  breakdown,
        'unmatched':  unmatched,
    }


def apply_hi_to_economics(
    economics,                                  # flowsheet.economics.Economics
    hi_result:                  HIResult,
    heating_tiers:              List[UtilityTier],
    cooling_tiers:              List[UtilityTier],
    operating_hours_per_year:   float = 8000.0,
    depreciation_years:         int   = 8,
):
    """既存 Economics に HI (pinch targeting) を適用した新 Economics を返す。

    熱系 OPEX (utility による加熱・冷却) を HI tier 別 OPEX で置換する。
    非熱系 OPEX (触媒交換、吸着剤交換、原料費、電力) と Revenue・CAPEX は
    そのまま継承する。

    設計判断 (2026-05-09):
      Phase 1+2 では「HI targeting only」を実装。物理 HE network は合成せず、
      理論限界 OPEX (Q_H_min, Q_C_min を tier 配分したもの) を採用する。
      CAPEX は据え置き (実機の HE はそのまま、HI 適用で安くなった OPEX のみ反映)。
      Phase 3 で HEN CAPEX 置換を実装予定。

    Parameters
    ----------
    economics : Economics
        HI なし元データ (calculate_economics の戻り値)
    hi_result : HIResult
        pinch_analysis の戻り値
    heating_tiers, cooling_tiers : List[UtilityTier]
    operating_hours_per_year : float
    depreciation_years : int

    Returns
    -------
    Economics
        HI 適用後の新インスタンス (元のは破壊しない)
    """
    # 循環 import 回避のため遅延 import
    from copy import deepcopy
    from flowsheet.economics import (
        Economics, apply_hasebe_aggregation, _count_main_equipment,
    )

    # ---- HI 後 OPEX を構築 ----
    new_opex: Dict[str, float] = {}

    # 非熱系かつ Hasebe 集計項でない生エントリのみコピー (例: 電力、触媒、吸着剤、原料費)
    # Hasebe 集計項は apply_hasebe_aggregation 内部で剥がして再計算する。
    for k, v in economics.opex.items():
        if classify_heat_opex_key(k) is None:
            new_opex[k] = v

    # HI tier 別 OPEX を「HI:」プレフィクス付きで追加 (用役費 = Hasebe C_UT に該当)
    hi_opex_calc = calc_hi_opex_okuyen(
        hi_result, heating_tiers, cooling_tiers,
        operating_hours_per_year=operating_hours_per_year,
    )
    for tier_name, cost in hi_opex_calc['breakdown'].items():
        new_opex[f'HI: {tier_name}'] = cost

    # ---- Hasebe 集計項を再計算 (C_UT が HI で変化したため delta が変わる) ----
    N_eq     = _count_main_equipment(economics.capex)
    new_opex = apply_hasebe_aggregation(new_opex, economics.total_capex, N_eq)

    # ---- 集計値の再計算 ----
    new_total_opex    = sum(new_opex.values())
    new_TAC           = economics.total_capex / depreciation_years + new_total_opex
    new_profit        = economics.total_revenue - new_TAC
    new_unit_jpy_per_t = (
        new_TAC * 1.0e8 / (economics.annual_kg_C3H6 / 1000.0)
        if economics.annual_kg_C3H6 > 0 else float('inf')
    )

    return Economics(
        capex          =deepcopy(economics.capex),
        opex           =new_opex,
        revenue        =deepcopy(economics.revenue),
        total_capex    =economics.total_capex,
        total_opex     =new_total_opex,
        total_revenue  =economics.total_revenue,
        TAC            =new_TAC,
        profit         =new_profit,
        annual_kg_C3H6 =economics.annual_kg_C3H6,
        unit_jpy_per_t =new_unit_jpy_per_t,
    )


def get_default_utility_tiers() -> Tuple[List[UtilityTier], List[UtilityTier]]:
    """プロジェクト標準のユーティリティ tier 定義を取得する。

    src.utility_selector の `_HEATING_TIERS` / `_COOLING_TIERS` を
    本モジュールの `UtilityTier` 型に変換して返す。

    BO ループや exp スクリプト、テストで毎回 tier を手書きする代わりに
    この関数を呼ぶ。tier 定義の変更は src.utility_selector 一箇所で完結する。

    Returns
    -------
    (heating_tiers, cooling_tiers) : Tuple[List[UtilityTier], List[UtilityTier]]
        - heating_tiers: LP/MP/HP Steam, 燃料燃焼 (低 supply_T → 高 supply_T)
        - cooling_tiers: 空冷, 冷却水, プロピレン冷媒, エチレン冷媒
    """
    from src.utility_selector import _HEATING_TIERS, _COOLING_TIERS

    heating = [
        UtilityTier(
            name=t.name,
            supply_T_K=t.supply_T_K,
            jpy_per_GJ=t.jpy_per_GJ,
            is_heating=True,
        )
        for t in _HEATING_TIERS
    ]
    cooling = [
        UtilityTier(
            name=t.name,
            supply_T_K=t.supply_T_K,
            jpy_per_GJ=t.jpy_per_GJ,
            is_heating=False,
        )
        for t in _COOLING_TIERS
    ]
    return heating, cooling


def plot_GCC(result: HIResult, ax=None, T_in_celsius: bool = True):
    """グランドコンポジットカーブを描画する。"""
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()

    def Tdisp(T_K: float) -> float:
        return T_K - 273.15 if T_in_celsius else T_K

    if result.GCC:
        H = [p[1] for p in result.GCC]
        T = [Tdisp(p[0]) for p in result.GCC]
        ax.plot(H, T, 'g-o', label='GCC', linewidth=2)
        # ピンチ点に印
        T_pinch_disp = Tdisp(0.5 * (result.T_pinch_hot_K + result.T_pinch_cold_K))
        ax.axhline(T_pinch_disp, color='gray', linestyle='--', alpha=0.5,
                   label=f'pinch (shifted)')

    ax.set_xlabel('H_residual [kW]')
    ax.set_ylabel('T_shifted [°C]' if T_in_celsius else 'T_shifted [K]')
    ax.set_title(
        f'グランドコンポジットカーブ  '
        f'(Q_H_min={result.Q_H_min_kW:.1f}, Q_C_min={result.Q_C_min_kW:.1f} kW)'
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax
