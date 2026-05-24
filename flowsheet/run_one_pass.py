"""
PDH プロセスの 1 パス計算 (リサイクル収束の 1 反復に相当)。

入力:
  - tear streams (Mem retentate / Dist3 bottom それぞれの A=C3H8, B=C3H6)
  - Fresh LPG 流量 (C3H8, C4H10)
  - 設計変数バンドル (FlowsheetDesignVars)
  - 運転条件 (OperatingConfig)

戻り値: dict (各ユニットの結果と tear stream 更新値・温度)

外側のソルバ (flowsheet/solver.py) がこの関数を反復呼び出しして
リサイクルを収束させる。
"""

import math
import os
import warnings
from typing import Dict, Tuple

from stream.stream import ProcessStream
from units.utils.mixer import mix_streams
from units.utils.cooler import simulate_cooler
from flowsheet.heat_integration import StreamPhase
from units.utils.compressor import simulate_compressor
from units.utils.pump import simulate_pump
from units.utils.expansion_valve import simulate_jt_expansion
from units.separators.column1.column1 import simulate_column1
from units.separators.column2.column2 import simulate_column2
from units.separators.column3.column3 import simulate_column3
from units.reactors.swing import (
    FeedStream as SwingFeed, FixedParams as SwingFixed,
    simulate_swing_reactor_system,
)
from units.separators.psa.psa_system import (
    PSAFeedStream, PSAFixedParams, simulate_psa_system,
)
from units.separators.membrane.membrane_system import (
    MemFeedStream, MemFixedParams, simulate_membrane_system,
)

from flowsheet.design import FlowsheetDesignVars
from config.load import OperatingConfig
from src.cost_parameters import PENALTY_CAPEX_THRESHOLD_OKUYEN

_ZERO = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'D': 0.0, 'E': 0.0, 'F': 0.0, 'Z': 0.0}


# ---------------------------------------------------------------------------
# trace bypass (2026-05-19 追加、モデル簡略化に対する数値補正)
# ---------------------------------------------------------------------------
#
# 設計判断 (2026-05-19): PSA / Mem の design モデルは入口組成に「主要成分のみ」
# の前提を置いている (PSA は CH4 破過計算、Mem は C3H6/C3H8 二成分透過)。
# 上流 (Dist2 partial cond) で微量の不純物 (例: C3H6 漏れ → PSA 入口に C3,
# C2H6 漏れ → Mem 入口に C2) が流れ込むと簡略モデルが破綻し _penalty_result
# を返してしまう。
#
# 本来は PSA / Mem を多成分対応に拡張すべき (TODO) だが、当面は orchestration
# 層で「閾値未満の微量成分を design 計算から除き、マスバランスは保ったまま下流
# (PSA: offgas, Mem: retentate=recycle) に直接ルーティング」する近似処置を入れる。
# これは「物理装置の追加」ではなく「シミュレータの数値処理」であり、物理嘘の
# force-move とは区別される (= マスは保たれる、設計計算の適用範囲を明示)。
#
# 閾値超過時は warning を出して「もう微量とは言えない量」を明示。
_TRACE_BYPASS_FRAC = 0.01    # 入口総モル流量に対する微量判定閾値 (1%)
# 設計判断 (2026-05-19 確定): 1% で固定。
# 当初 5% / 15% への拡張を試したが、これは「Dist2 設計の不味さ (10%+ 漏れ) を
# シミュレータの bypass で隠す」物理嘘になっていた。本来の意図は:
#   - Dist2 を設計でしっかり詰めて、物理的に C3 漏れ <1% を達成する
#   - 残った微量 (≤ 1%) を bypass で吸収 (= 数値処理として透過的)
# つまり 1% 閾値は「設計が正しいことの保証」、超えたら設計が悪い (=本物の警告)。
# BO は探索範囲を制約することで <1% 領域を探すよう誘導 (search_space.py を見直し)。
_PSA_TRACE_COMPS = ('A', 'B')                # PSA で許容しない: C3H8, C3H6
_MEM_TRACE_COMPS = ('C', 'D', 'E', 'F')      # Mem で許容しない: H2, CH4, C2H4, C2H6


# 設計判断 (2026-05-22): iter ごとにどのユニットで penalty が出たかを stderr/stdout に
# 出すロガー。環境変数 PDH_PER_UNIT_LOG=1 で有効化。BO ループで巨大ログにならないよう
# デフォルト OFF。
_PER_UNIT_LOG = os.environ.get('PDH_PER_UNIT_LOG', '0') == '1'

def _log_unit_failure(unit_name: str, equipment) -> None:
    """ユニットの penalty 発火を stderr に1行で出す (PDH_PER_UNIT_LOG=1 時のみ)。"""
    if not _PER_UNIT_LOG:
        return
    import sys as _sys
    feasible = getattr(equipment, 'feasible', None)
    capex = getattr(equipment, 'CAPEX', None)
    capex_total = getattr(equipment, 'CAPEX_total', None)
    reactor_capex = getattr(equipment, 'Reactor_CAPEX', None)
    msg = (getattr(equipment, 'message', '') or
           getattr(equipment, 'penalty_reason', '') or '')
    _sys.stderr.write(
        f"[PENALTY] {unit_name}: feasible={feasible} "
        f"CAPEX={capex} CAPEX_total={capex_total} R_CAPEX={reactor_capex}"
        f" msg={str(msg)[:150]}\n"
    )
    _sys.stderr.flush()


def _apply_trace_bypass(
    F_in:           Dict[str, float],
    trace_comps:    tuple,
    threshold_frac: float,
    label:          str,
) -> Tuple[Dict[str, float], Dict[str, float], float]:
    """入口流量から指定成分の微量 (= total の threshold_frac 未満) を抽出。

    Parameters
    ----------
    F_in : dict
        入口の成分別モル流量 [kmol/h]。
    trace_comps : tuple
        微量判定対象の成分キー ('A'〜'F','Z')。
    threshold_frac : float
        F_in 全量に対する微量判定の閾値 (= 0.01 で 1%)。
    label : str
        warning 用のユニット識別子 ("PSA" / "Mem")。

    Returns
    -------
    (cleaned_F, bypass_F, max_excess_frac)
        cleaned_F      : 微量分を除いた入口流量 (ユニットの design 計算へ)。
        bypass_F       : 除かれた微量分 (ユニットの出口に合算してマスバランス保持)。
        max_excess_frac: trace_comps の中で frac > threshold だった成分の
                         (frac - threshold) の最大値 [-]。超過なしなら 0。
                         BO objective の連続 penalty 化に runner.py で利用。
    """
    # 設計判断 (2026-05-20、ユーザー方針): 閾値超過時の warning は出さない。
    # 過去版は「PSA/Mem の多成分対応化を検討」と促していたが、多成分化は実施しない
    # 方針確定 (理由: 1-2% 漏れは有効数字範囲で許容、多成分化はコスト対効果悪い)。
    # max_excess_frac は将来 BO penalty 化する余地として戻り値に残す (現状未使用)。
    F_total = sum(max(F, 0.0) for F in F_in.values())
    cleaned: Dict[str, float] = dict(F_in)
    bypass:  Dict[str, float] = {c: 0.0 for c in F_in}
    max_excess: float = 0.0
    if F_total <= 0:
        return cleaned, bypass, max_excess
    for c in trace_comps:
        v = F_in.get(c, 0.0)
        if v <= 0:
            continue
        frac = v / F_total
        if frac <= threshold_frac:
            # 微量 → bypass
            cleaned[c] = 0.0
            bypass[c]  = v
        else:
            # 閾値超え: 物理値を残して通過させる (warning 出さない)
            excess = frac - threshold_frac
            if excess > max_excess:
                max_excess = excess
    return cleaned, bypass, max_excess


# ---------------------------------------------------------------------------
# Warning 集約ヘルパ
# ---------------------------------------------------------------------------
# 旧版 (〜2026-05-17) は warnings.simplefilter("ignore") で全 warning を抑制し、
# fallback の発生 (PR EOS Z=1 fallback、brentq 偽根、Wang-Henke MESH 残差超過に
# よる FUG fallback 等) が静かに進行していた。これが原因で「数値結果は正常だが
# 実は silent fallback が走っていて TAC が 1 億円/年単位で過小評価」という事案が
# 2026-05-10 に発覚した (STATUS_2026-05-10.md 参照)。
#
# 本実装 (2026-05-18) では catch_warnings(record=True) + simplefilter("always") で
# warning を抑制せず捕捉し、ラベル付きで結果 dict に格納する。runner.py が
# failure_reason に集約することで、BO log で fallback 発火が追跡可能になる。
class _CapturedWarning:
    """source ラベル付きの warning エントリ。"""
    __slots__ = ('source', 'category', 'message', 'filename', 'lineno')

    def __init__(self, source: str, w):
        self.source   = source
        self.category = w.category.__name__
        self.message  = str(w.message)
        self.filename = w.filename
        self.lineno   = w.lineno

    def __repr__(self) -> str:
        return f"[{self.source}] {self.category}: {self.message}"


def _capture_warnings(source: str, captured_list: list):
    """`with _capture_warnings("Dist1", warnings_log):` の形で使うヘルパ。

    内部の simplefilter は "always" を設定するため、warning は抑制されずに
    captured_list に append される。元の warning ストリームには出ない (record=True)。
    """
    import contextlib

    @contextlib.contextmanager
    def _ctx():
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            yield
            for w in ws:
                captured_list.append(_CapturedWarning(source, w))

    return _ctx()


# ---------------------------------------------------------------------------
# 早期 penalty 終了用ヘルパ (reactor で SV check 等が NG だった場合)
# ---------------------------------------------------------------------------

class _PenaltyEquipment:
    """downstream units の equipment.CAPEX_total / CAPEX_total / 等を捏造するスタブ。"""
    CAPEX_total: float = 1e9   # solver の penalty_hit 検査が拾えるよう sentinel
    CAPEX:       float = 1e9
    Reactor_CAPEX: float = 1e9


class _PenaltyResult:
    """downstream units の擬似結果。equipment 属性経由で sentinel CAPEX を返す。"""
    equipment = _PenaltyEquipment()
    top    = ProcessStream(F_in=dict(_ZERO), T_in=298.15, P_in=1e5)
    bottom = ProcessStream(F_in=dict(_ZERO), T_in=298.15, P_in=1e5)
    outlet = ProcessStream(F_in=dict(_ZERO), T_in=298.15, P_in=1e5)


def _build_penalty_one_pass_result(
    r_rx, reactor_inlet, dist1_top_rx, recycle_dist3, recycle_mem,
    r1, fresh, pump1,
):
    """reactor が penalty 返却したとき、下流装置を全部スタブにして solver に渡す形を構築。

    solver は `results['r_rx'].equipment.Reactor_CAPEX >= PENALTY_CAPEX_THRESHOLD_OKUYEN`
    を penalty_hit と判定するため、本関数を経由しても適切に penalty 化される。下流装置の
    P_in=0 例外を回避することが主目的。
    """
    stub = _PenaltyResult()
    result = dict(
        pump1=pump1, r1=r1, dist1_top_rx=dist1_top_rx,
        fresh=fresh,
        reactor_inlet=reactor_inlet,
        r_rx=r_rx, rx_out=ProcessStream(F_in=dict(_ZERO), T_in=298.15, P_in=1e5),
        cooled=stub,
        comp2a=stub, intercool=stub, comp2b=stub,
        desuper=stub,
        r2=stub,
        r_psa=stub, mem_precool=stub, r_mem=stub, r3=stub,
        tear_dist3_new={'A': 0.0, 'B': 0.0},
        tear_mem_new  ={'A': 0.0, 'B': 0.0},
        T_d3_new=298.15, T_mem_new=298.15,
        warnings_captured=[],
        trace_bypass_psa_excess=0.0,
        trace_bypass_mem_excess=0.0,
        dist1_N_shortfall=0.0, dist2_N_shortfall=0.0, dist3_N_shortfall=0.0,
        dist1_dT_shortfall=0.0, dist2_dT_shortfall=0.0, dist3_dT_shortfall=0.0,
        # Mem shortfall (Reactor 失敗で Mem 未到達 = 0 で伝播)
        mem_ph_shortfall=0.0, mem_bp_shortfall=0.0,
        mem_phase_shortfall=0.0, mem_other_shortfall=0.0,
        # 観測ラベル (2026-05-22)
        first_failed_unit='r_rx',
    )
    # 上流の r1 + reactor の penalty_reason / SV_actual を抽出。下流は stub なので無視。
    result.update(_EMPTY_UNIT_DIAG)
    result.update(_extract_unit_diagnostics(r1=r1, r_rx=r_rx))
    return result


# ---------------------------------------------------------------------------
# 蒸留塔 penalty 早期検出 (2026-05-20 追加)
# ---------------------------------------------------------------------------
# 設計判断: 旧版は r1/r2/r3 が penalty (CAPEX=1e9 sentinel, feasible=False) で
# 全成分流量ゼロの top/bottom を返したあと下流の simulate_jt_expansion / PSA / Mem
# に流れ、`ValueError: expansion_valve: 全成分流量がゼロです` で例外 crash していた
# (BO #001 trial 等で多発)。runner.py は例外を catch するが、TPE constraints_func に
# 渡る情報が「solver_failure 固定」になり連続シグナルが失われる。
# 本ヘルパは塔 penalty を早期検出して solver に penalty_hit で抜けてもらう経路に
# 変換し、同時に N_shortfall / dT_shortfall を計算して runner.py 経由で
# trial.user_attrs に格納する (TPE が「あとどれだけ N/dT が足りない領域か」を学習)。
_RIG_TOL_K = 0.05   # Wang-Henke 収束 tol (src/distillation_rigorous.py と同期)
# PSA penalty 連続シグナルの基準 (psa_system.py の _T_ABS_MIN, _U0_MAX と同期)
_PSA_T_ABS_MIN_S = 60.0
_PSA_U0_MAX_MS   = 1.0
# Reactor SV penalty 連続シグナルの基準 (swing.py の FixedParams.SV_{min,max}_m_per_s と同期)
_REACTOR_SV_MIN_MS = 0.5
_REACTOR_SV_MAX_MS = 3.0
# Mem penalty 連続シグナル用の安全マージン (2026-05-22)
#   ph_le_pfeed:  P_H が feed.P_in に何倍近いか → log10 比
#   bp_le_cold:   bp が cold_out に何 K 不足か → 対数比 (T 単位)
#   vapor_cond:   T_in が T_dew に何 K 不足か → 対数比
_MEM_T_MARGIN_K = 5.0   # bp - cold_out / T_in - dew で「健全」とみなす最小マージン


def _compute_reactor_shortfall(r_rx) -> Dict[str, float]:
    """Reactor penalty 経路から TPE constraints_func 用の連続 shortfall を計算。

    Returns
    -------
    dict
        - 'reactor_sv_shortfall' [-] : SV 範囲外への log10 距離。
          SV < min → log10(SV_MIN / SV)、SV > max → log10(SV / SV_MAX)、範囲内 → 0。
        - 'reactor_other_shortfall' [-] : SV 以外の penalty (input_invalid, sim_failure,
          volume_zero, capex_exception) で 1.0、正常完走で 0。

    設計判断 (2026-05-21): PSA shortfall と同じパターン。SV < min なら D を小さく or
    並列数を増やす方向、SV > max なら D を大きく or 並列数を減らす方向に TPE が
    学習できるよう正負を log10 比で連続化。
    """
    out = {
        'reactor_sv_shortfall':    0.0,
        'reactor_other_shortfall': 0.0,
    }
    eq = getattr(r_rx, 'equipment', None)
    if eq is None:
        return out
    reason = getattr(eq, 'penalty_reason', '') or ''
    if reason == '':
        return out  # 正常完走
    if reason == 'sv_out_of_range':
        sv = getattr(eq, 'SV_actual', 0.0) or 0.0
        if sv > 0 and sv < _REACTOR_SV_MIN_MS:
            out['reactor_sv_shortfall'] = math.log10(max(_REACTOR_SV_MIN_MS / sv, 1.0))
        elif sv > _REACTOR_SV_MAX_MS:
            out['reactor_sv_shortfall'] = math.log10(max(sv / _REACTOR_SV_MAX_MS, 1.0))
        else:
            out['reactor_sv_shortfall'] = 1.0  # 想定外 (SV=0 等)
    else:
        # input_invalid / sim_failure / volume_zero / capex_exception
        out['reactor_other_shortfall'] = 1.0
    return out


def _compute_psa_shortfall(r_psa) -> Dict[str, float]:
    """PSA penalty 経路から TPE constraints_func 用の連続 shortfall を計算。

    Returns
    -------
    dict
        - 'psa_t_abs_shortfall' [-] : log10(_T_ABS_MIN / t_abs_actual)。
          t_abs_below_min 経路でのみ > 0 (例: t_abs=37s → log10(60/37)≈0.21)。
          mask_lt_2 でも t_abs_actual が伝わるので 0 になりうる。
        - 'psa_u_0_shortfall'  [-] : log10(u_0_actual / _U0_MAX)。
          u_0_above_max 経路でのみ > 0 (例: u_0=1.5m/s → log10(1.5/1.0)≈0.18)。
        - 'psa_feed_shortfall' [-] : feed 異常 (no_non_C3_feed / no_CH4_feed) で 1.0。
          これらは上流 (Dist2 rigorous) が異常組成を返した場合に発火する稀ケース。

    設計判断 (2026-05-21): silent _penalty_result() の代替シグナル。BO の TPE
    constraints_func で feasible 境界への接近度を学習させる目的。正常完走時は
    すべて 0.0 を返す (penalty_reason='')。
    """
    out = {
        'psa_t_abs_shortfall': 0.0,
        'psa_u_0_shortfall':   0.0,
        'psa_feed_shortfall':  0.0,
    }
    eq = getattr(r_psa, 'equipment', None)
    if eq is None:
        return out
    reason = getattr(eq, 'penalty_reason', '') or ''
    if reason == '':
        return out  # 正常完走

    t_abs = getattr(eq, 't_abs_actual_s', 0.0) or 0.0
    u_0   = getattr(eq, 'u_0_actual',     0.0) or 0.0

    if reason == 't_abs_below_min' and t_abs > 0:
        out['psa_t_abs_shortfall'] = math.log10(max(_PSA_T_ABS_MIN_S / t_abs, 1.0))
    elif reason == 'u_0_above_max' and u_0 > 0:
        out['psa_u_0_shortfall'] = math.log10(max(u_0 / _PSA_U0_MAX_MS, 1.0))
    elif reason in ('no_non_C3_feed', 'no_CH4_feed'):
        out['psa_feed_shortfall'] = 1.0
    elif reason == 'breakthrough_no_converge':
        # 物理的には t_abs > t_ads_max (= 7200s) → 「破過遅すぎ」= 設計過大
        # BO 視点では desorption_target を下げる or L を下げる方向。
        # 連続シグナルは作りにくいので t_abs_shortfall に弱い負方向 (= -0.5) を入れる
        # ... よりは別エントリで 1.0 を入れて TPE に「ここは避ける」のみ伝える。
        out['psa_t_abs_shortfall'] = 1.0
    elif reason == 'mask_lt_2':
        # 通常は t_abs_below_min と同じ症状の派生
        if t_abs > 0:
            out['psa_t_abs_shortfall'] = math.log10(max(_PSA_T_ABS_MIN_S / t_abs, 1.0))
        else:
            out['psa_t_abs_shortfall'] = 1.0
    return out


def _compute_mem_shortfall(r_mem) -> Dict[str, float]:
    """Mem penalty 経路から TPE constraints_func 用の連続 shortfall を計算。

    Returns
    -------
    dict
        - 'mem_ph_shortfall'    [-] : P_H が feed.P_in に対し不足している log10 比。
              ph_le_pfeed 経路で `log10(feed.P_in / P_H)` > 0 → BO に「P_H を上げる
              or Dist2 圧力を下げる」シグナル。ph_le_pl は P_L=1atm 固定運用上ほぼ
              発火しないが、念のため 1.0 を入れる (search bounds で防げる領域)。
        - 'mem_bp_shortfall'    [-] : 透過ガス泡点が冷却水出口温度を下回る不足比。
              bp_le_cold_out 経路で `log10((T_cold_out + margin) / T_bp_perm)` > 0
              → BO に「P_dist3 を上げる or 透過組成を C3H6 純度寄りにする」シグナル。
        - 'mem_phase_shortfall' [-] : ガスフィード前提に対し T_in が露点を下回る不足比。
              vapor_condensed 経路で `log10(T_dew / T_in)` > 0 → BO に「Dist2 圧力を
              下げる (= 露点を下げる)」シグナル。T_in は config 固定なので露点側を動かす。
        - 'mem_other_shortfall' [-] : 上記以外の numerical 失敗 (vap/ode/comp/cond
              exception, nan, invalid_input, neg_feed 等) で 1.0、正常完走で 0。

    設計判断 (2026-05-22): PSA / Reactor で先行導入したパターンに揃える。
    silent _penalty_result() 経路 (membrane_system.py の 16 箇所) が BO に方向を
    渡せず、main_20260522_005631 で 80% の trial が PSA/Mem CAPEX sentinel hit で
    無方向に死んでいた問題への対処。Mem feed は Dist2 塔底 → mem_precool で
    config.temperature.mem_feed_K (= 323.15K 固定) まで加熱されてから入るので、
    露点・bp 関連の shortfall は実質的に Dist2 圧力 / Dist3 圧力を動かすシグナル
    として効く。
    """
    out = {
        'mem_ph_shortfall':    0.0,
        'mem_bp_shortfall':    0.0,
        'mem_phase_shortfall': 0.0,
        'mem_other_shortfall': 0.0,
    }
    eq = getattr(r_mem, 'equipment', None)
    if eq is None:
        return out
    reason = getattr(eq, 'penalty_reason', '') or ''
    if reason == '':
        return out  # 正常完走

    P_H        = getattr(eq, 'P_H_actual_Pa',      0.0) or 0.0
    P_feed     = getattr(eq, 'P_feed_actual_Pa',   0.0) or 0.0
    T_dew      = getattr(eq, 'T_dew_actual_K',     0.0) or 0.0
    T_feed     = getattr(eq, 'T_feed_actual_K',    0.0) or 0.0
    T_bp       = getattr(eq, 'T_bp_perm_actual_K', 0.0) or 0.0
    T_cold_out = getattr(eq, 'T_cold_out_actual_K',0.0) or 0.0

    if reason == 'ph_le_pfeed' and P_H > 0 and P_feed > 0:
        out['mem_ph_shortfall'] = math.log10(max(P_feed / P_H, 1.0))
    elif reason == 'ph_le_pl':
        # P_L=1atm 固定運用ではほぼ発火しないが、search bounds が壊れた場合の保険
        out['mem_ph_shortfall'] = 1.0
    elif reason == 'bp_le_cold_out' and T_bp > 0 and T_cold_out > 0:
        # 「マージン込みで bp が cold_out + margin を下回る分」を log10 で連続化
        target = T_cold_out + _MEM_T_MARGIN_K
        out['mem_bp_shortfall'] = math.log10(max(target / T_bp, 1.0))
    elif reason == 'vapor_condensed' and T_dew > 0 and T_feed > 0:
        # T_in < T_dew → log10(T_dew / T_in)。Dist2 P を下げる方向のシグナル
        out['mem_phase_shortfall'] = math.log10(max(T_dew / T_feed, 1.0))
    elif reason == 'liquid_vaporized' and T_dew > 0 and T_feed > 0:
        # 本 run では vapor_feed=True なので発火しないが、対称性のため入れる
        out['mem_phase_shortfall'] = math.log10(max(T_feed / T_dew, 1.0))
    else:
        # invalid_input / invalid_design / pdist_le_pl / neg_feed / zero_feed /
        # dew_nan / vap_exception / vap_nan / feed_comp_exception / ode_failure /
        # prod_comp_exception / cond_exception / cond_nan
        out['mem_other_shortfall'] = 1.0
    return out


def _compute_dist_shortfalls(col_key: str, col_result, col_design) -> Dict[str, float]:
    """塔 penalty 結果から N_shortfall / dT_shortfall を計算。

    col_key   : 'r1' | 'r2' | 'r3'
    col_result: DistResult (equipment.N_needed / dT_max_rigorous を保持)
    col_design: ColumnTunables (N_stages を持つ)
    """
    idx = col_key[-1]
    out = {f'dist{idx}_N_shortfall': 0.0, f'dist{idx}_dT_shortfall': 0.0}
    eq = getattr(col_result, 'equipment', None)
    if eq is None:
        return out
    # FUG (Gilliland infeasible) → N_needed > 0
    N_needed = getattr(eq, 'N_needed', 0.0) or 0.0
    N_stages = max(int(getattr(col_design, 'N_stages', 0) or 0), 1)
    if N_needed > 0:
        out[f'dist{idx}_N_shortfall'] = max(0.0, (N_needed - N_stages) / N_stages)
    # Rigorous (Wang-Henke 収束失敗) → dT_max_rigorous > 0
    dT_max = getattr(eq, 'dT_max_rigorous', 0.0) or 0.0
    if dT_max > 0:
        # tol からの「比」で正規化 (tol=0.05K に対して 56K → shortfall=1120)。
        # スケールが過大なので log 圧縮で 0-10 程度に抑える。
        # 例: dT_max=56K → log10(56/0.05) ≈ log10(1120) ≈ 3.05
        out[f'dist{idx}_dT_shortfall'] = math.log10(max(dT_max / _RIG_TOL_K, 1.0))
    return out


# ---------------------------------------------------------------------------
# 観測ラベル抽出 (2026-05-22 L1 観測強化、ユーザー要望)
# ---------------------------------------------------------------------------
# 設計判断: shortfall (連続値) は既に user_attrs に流れているが、各装置が持つ
# 「どの penalty_reason ラベルで死んだか」「実 actual 値はいくつだったか」は
# trial に届いていない。BO 走行中に「Mem.bp_le_cold_out (T_bp=305K < T_cold=313K)」
# のような具体ラベル + 数値を live 表示できるよう、result dict に構造化フィールド
# を埋める。表示は callbacks.py 側。

# 装置別診断フィールドの default 値 (run_one_pass の戻り dict キーを揃える)。
# str フィールドは '' で空、float は 0.0 で「未計測/未該当」を示す。
_EMPTY_UNIT_DIAG: Dict[str, object] = {
    'reactor_penalty_reason':   '',
    'reactor_SV_actual_m_s':    0.0,
    'psa_penalty_reason':       '',
    'psa_t_abs_actual_s':       0.0,
    'psa_u_0_actual_m_s':       0.0,
    'mem_penalty_reason':       '',
    'mem_P_H_actual_Pa':        0.0,
    'mem_P_feed_actual_Pa':     0.0,
    'mem_T_dew_actual_K':       0.0,
    'mem_T_feed_actual_K':      0.0,
    'mem_T_bp_perm_actual_K':   0.0,
    'mem_T_cold_out_actual_K':  0.0,
    'r1_penalty_msg':           '',
    'r1_N_needed':              0.0,
    'r1_dT_max_K':              0.0,
    'r2_penalty_msg':           '',
    'r2_N_needed':              0.0,
    'r2_dT_max_K':              0.0,
    'r3_penalty_msg':           '',
    'r3_N_needed':              0.0,
    'r3_dT_max_K':              0.0,
}


def _extract_unit_diagnostics(**unit_results) -> Dict[str, object]:
    """各装置の equipment から penalty_reason + key actual 値を 1 dict に抽出。

    Parameters
    ----------
    unit_results : kwargs
        r1, r_rx, r2, r_psa, r_mem, r3 のうち計算済みのものを渡す。
        None や stub (penalty_reason='') は無視され default 値のまま。

    Returns
    -------
    dict
        _EMPTY_UNIT_DIAG と同じキーセット。BO の trial.user_attrs にそのまま展開可能。
    """
    d: Dict[str, object] = dict(_EMPTY_UNIT_DIAG)

    # Reactor (swing.py で penalty_reason ラベル + SV_actual を埋め込み済み)
    r_rx = unit_results.get('r_rx')
    eq = getattr(r_rx, 'equipment', None) if r_rx is not None else None
    if eq is not None:
        reason = getattr(eq, 'penalty_reason', '') or ''
        if reason:
            d['reactor_penalty_reason'] = reason
            sv = getattr(eq, 'SV_actual', 0.0) or 0.0
            if sv > 0:
                d['reactor_SV_actual_m_s'] = float(sv)

    # PSA (psa_system.py で penalty_reason + t_abs_actual_s + u_0_actual を埋め込み)
    r_psa = unit_results.get('r_psa')
    eq = getattr(r_psa, 'equipment', None) if r_psa is not None else None
    if eq is not None:
        reason = getattr(eq, 'penalty_reason', '') or ''
        if reason:
            d['psa_penalty_reason'] = reason
            t_abs = getattr(eq, 't_abs_actual_s', 0.0) or 0.0
            u_0   = getattr(eq, 'u_0_actual',     0.0) or 0.0
            if t_abs > 0:
                d['psa_t_abs_actual_s'] = float(t_abs)
            if u_0 > 0:
                d['psa_u_0_actual_m_s'] = float(u_0)

    # Mem (membrane_system.py で penalty_reason + 6 つの T/P actual を埋め込み)
    r_mem = unit_results.get('r_mem')
    eq = getattr(r_mem, 'equipment', None) if r_mem is not None else None
    if eq is not None:
        reason = getattr(eq, 'penalty_reason', '') or ''
        if reason:
            d['mem_penalty_reason'] = reason
            for src_attr, dst_key in (
                ('P_H_actual_Pa',       'mem_P_H_actual_Pa'),
                ('P_feed_actual_Pa',    'mem_P_feed_actual_Pa'),
                ('T_dew_actual_K',      'mem_T_dew_actual_K'),
                ('T_feed_actual_K',     'mem_T_feed_actual_K'),
                ('T_bp_perm_actual_K',  'mem_T_bp_perm_actual_K'),
                ('T_cold_out_actual_K', 'mem_T_cold_out_actual_K'),
            ):
                v = getattr(eq, src_attr, 0.0) or 0.0
                if v > 0:
                    d[dst_key] = float(v)

    # 蒸留塔 r1/r2/r3 (DistEquipment.feasible=False で penalty、message + N_needed + dT_max)
    for col_key in ('r1', 'r2', 'r3'):
        col = unit_results.get(col_key)
        eq = getattr(col, 'equipment', None) if col is not None else None
        if eq is None:
            continue
        if getattr(eq, 'feasible', True):
            continue  # 正常塔は無視
        msg = getattr(eq, 'message', '') or ''
        if msg:
            d[f'{col_key}_penalty_msg'] = str(msg)[:80]
        n_needed = getattr(eq, 'N_needed', 0.0) or 0.0
        if n_needed > 0:
            d[f'{col_key}_N_needed'] = float(n_needed)
        dt_max = getattr(eq, 'dT_max_rigorous', 0.0) or 0.0
        if dt_max > 0:
            d[f'{col_key}_dT_max_K'] = float(dt_max)
    return d


def _determine_first_failed_unit(**unit_results) -> str:
    """パイプライン順 (r1→r_rx→r2→r_psa→r_mem→r3) で最初に penalty を起こしたユニットを返す。

    PSA/Mem は早期 return しないため、success path でも r_psa/r_mem に penalty_reason
    が刺さっている可能性がある (solver は CAPEX sentinel で penalty_hit を後から判定する)。
    本関数はそれを取りこぼさず、Dist3 まで素通りした trial でも「実は PSA で死んでた」
    を正しくラベル付けする。

    Returns
    -------
    str : 'r1' | 'r_rx' | 'r2' | 'r_psa' | 'r_mem' | 'r3' | ''
          '' は run_one_pass 内 penalty なし (solver/spec レイヤーで判定する)
    """
    r1 = unit_results.get('r1')
    eq = getattr(r1, 'equipment', None) if r1 is not None else None
    if eq is not None and not getattr(eq, 'feasible', True):
        return 'r1'

    r_rx = unit_results.get('r_rx')
    eq = getattr(r_rx, 'equipment', None) if r_rx is not None else None
    if eq is not None and (
        (getattr(eq, 'penalty_reason', '') or '') != '' or
        (getattr(eq, 'Reactor_CAPEX', 0.0) or 0.0) >= PENALTY_CAPEX_THRESHOLD_OKUYEN
    ):
        return 'r_rx'

    r2 = unit_results.get('r2')
    eq = getattr(r2, 'equipment', None) if r2 is not None else None
    if eq is not None and not getattr(eq, 'feasible', True):
        return 'r2'

    r_psa = unit_results.get('r_psa')
    eq = getattr(r_psa, 'equipment', None) if r_psa is not None else None
    if eq is not None and (getattr(eq, 'penalty_reason', '') or '') != '':
        return 'r_psa'

    r_mem = unit_results.get('r_mem')
    eq = getattr(r_mem, 'equipment', None) if r_mem is not None else None
    if eq is not None and (getattr(eq, 'penalty_reason', '') or '') != '':
        return 'r_mem'

    r3 = unit_results.get('r3')
    eq = getattr(r3, 'equipment', None) if r3 is not None else None
    if eq is not None and not getattr(eq, 'feasible', True):
        return 'r3'

    return ''


def _build_penalty_after_column(
    failed_col_key: str,
    design,
    **upstream,
) -> Dict:
    """蒸留塔 penalty 時に「上流の有効ピース + 下流の stub」で one_pass dict を構築。

    Parameters
    ----------
    failed_col_key : 'r1' | 'r2' | 'r3'
    design         : FlowsheetDesignVars (N_stages 取得用)
    upstream       : 既に計算済みの有効ピース (fresh, pump1, r1, ...) を kwargs で

    Returns
    -------
    dict : run_one_pass の返却形と同じ key 構成。dist{1,2,3}_{N,dT}_shortfall を含む。
    """
    stub = _PenaltyResult()
    zero_stream = ProcessStream(F_in=dict(_ZERO), T_in=298.15, P_in=1e5)
    base = dict(
        pump1=stub, r1=stub, dist1_top_rx=zero_stream,
        fresh=stub,
        reactor_inlet=zero_stream,
        r_rx=stub, rx_out=zero_stream,
        cooled=stub,
        comp2a=stub, intercool=stub, comp2b=stub,
        desuper=stub,
        r2=stub,
        r_psa=stub, mem_precool=stub, r_mem=stub, r3=stub,
        tear_dist3_new={'A': 0.0, 'B': 0.0},
        tear_mem_new  ={'A': 0.0, 'B': 0.0},
        T_d3_new=298.15, T_mem_new=298.15,
        warnings_captured=[],
        trace_bypass_psa_excess=0.0,
        trace_bypass_mem_excess=0.0,
        dist1_N_shortfall=0.0, dist2_N_shortfall=0.0, dist3_N_shortfall=0.0,
        dist1_dT_shortfall=0.0, dist2_dT_shortfall=0.0, dist3_dT_shortfall=0.0,
        # Mem shortfall (上流失敗時の既定値、Mem に到達できなかった = 0 で伝播)
        mem_ph_shortfall=0.0, mem_bp_shortfall=0.0,
        mem_phase_shortfall=0.0, mem_other_shortfall=0.0,
    )
    base.update(upstream)
    # solver の penalty_hit 検出経路 (r_psa/r_mem/r_rx) に乗せるため、failed が
    # r1/r2/r3 でも r_psa の CAPEX を sentinel に設定。r_psa は stub (CAPEX_total=1e9)
    # でデフォルト sentinel になっているため追加処理不要。

    # shortfall 計算
    failed_col_design = {
        'r1': design.dist1, 'r2': design.dist2, 'r3': design.dist3,
    }.get(failed_col_key)
    if failed_col_design is not None:
        col_result = upstream.get(failed_col_key)
        if col_result is not None:
            shortfalls = _compute_dist_shortfalls(failed_col_key, col_result, failed_col_design)
            base.update(shortfalls)
    # 観測ラベル (2026-05-22): どの装置で詰まったかを明示。
    # 上流に passed kwargs のうち unit result だけ抜き出して diag 抽出。
    base['first_failed_unit'] = failed_col_key
    base.update(_EMPTY_UNIT_DIAG)
    unit_kwargs = {k: v for k, v in upstream.items()
                   if k in ('r1', 'r_rx', 'r2', 'r_psa', 'r_mem', 'r3')}
    base.update(_extract_unit_diagnostics(**unit_kwargs))
    return base


def run_one_pass(
    tear_dist3:   dict,
    tear_mem:     dict,
    T_d3:         float,
    T_mem:        float,
    F_C3H8_feed:  float,
    F_C4H10_feed: float,
    design:       FlowsheetDesignVars,
    config:       OperatingConfig,
) -> dict:
    """1 パスのプロセスシミュレーションを実行。

    Parameters
    ----------
    tear_dist3 : dict {'A','B'} [kmol/h]
        Dist3 塔底由来のリサイクル (前反復値)
    tear_mem : dict {'A','B'} [kmol/h]
        Mem 保留側由来のリサイクル (前反復値)
    T_d3, T_mem : float [K]
        各リサイクルの温度 (前反復値)
    F_C3H8_feed, F_C4H10_feed : float [kmol/h]
        Fresh LPG 中 C3H8, C4H10 流量 (外側で調整)
    design : FlowsheetDesignVars
        設計変数 (SWING, PSA, MEM)
    config : OperatingConfig
        運転条件 (圧力, 温度, 製品仕様, 原料状態)

    Returns
    -------
    dict
        全ユニット結果 + tear_*_new + T_*_new
    """
    P_rx = config.pressure.reactor_inlet_Pa  # 反応器入口・膨張弁後の圧力 (contest 規定 0.5 bar)

    # 1 パス全体の warning 捕捉用 (旧 simplefilter("ignore") の代替)
    warnings_captured: list = []

    # ---- Fresh LPG 原料 (30°C 飽和液、C3H8:C4H10 = 9:1, ~9.97 bar) ----
    fresh = ProcessStream(
        F_in={'A': F_C3H8_feed, 'Z': F_C4H10_feed,
              'B': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=config.feed.T_K, P_in=config.feed.P_Pa,
    )

    # ---- Step 1: Pump1 → Dist1 (Fresh のみ) ----
    # 原料は液 (30°C 飽和液) のため、Dist1 (17 bar) への昇圧は液送ポンプで行う。
    # contest §3-3-3 「加圧すべき箇所には、ポンプ(液) を入れること」に従う。
    # 設計判断 (2026-05-09): 出口圧力は design.dist1.P_col に同期する
    # (operating.toml の pump1_out_Pa は backward compat のため残置だが未使用)。
    pump1 = simulate_pump(fresh, P_out_target=design.dist1.P_col)
    with _capture_warnings("Dist1", warnings_captured):
        r1 = simulate_column1(pump1.outlet, tunables=design.dist1)

    # 設計判断 (2026-05-20): Dist1 FUG が _penalty_result (Gilliland infeasible 等) を
    # 返したとき、r1.top.F_in は全成分ゼロ → 下流 simulate_jt_expansion が
    # `ValueError: 全成分流量がゼロ` で crash する経路があった (BO trial #0,1,9,12,...で多発)。
    # 早期に penalty 経路に分岐し、shortfall を TPE 用 user_attr に格納する。
    if not getattr(r1.equipment, 'feasible', True):
        _log_unit_failure('Dist1 (r1)', r1.equipment)
        result = _build_penalty_after_column(
            'r1', design,
            fresh=fresh, pump1=pump1, r1=r1,
        )
        result['warnings_captured'] = warnings_captured
        return result

    # 塔頂を反応器圧力 (0.5 bar) に膨張 (C4 除去済みの C3 主成分)
    # 設計判断 (2026-05-08): 旧版は P を書き換えるだけで T を維持していたが
    # 物理的には等エンタルピー (JT) 膨張で温度低下する。反応器入口プレヒート
    # Q_preheat の見積もり精度を上げるため simulate_jt_expansion を経由する。
    # 膨張弁本体はコストフリー (装置 CAPEX/OPEX に計上しない、配管中の絞り弁)。
    dist1_top_rx = simulate_jt_expansion(
        ProcessStream(F_in=dict(r1.top.F_in), T_in=r1.top.T_in, P_in=r1.top.P_in),
        P_out=P_rx,
    )

    # ---- リサイクルストリーム (膨張弁経由、コストなし) ----
    # tear stream の元圧力:
    #   recycle_dist3 : Dist3 塔底 = design.dist3.P_col
    #   recycle_mem   : Mem 保留   = design.mem.P_H
    recycle_dist3 = simulate_jt_expansion(
        ProcessStream(
            F_in={**_ZERO, 'A': tear_dist3['A'], 'B': tear_dist3['B']},
            T_in=T_d3, P_in=design.dist3.P_col,
        ),
        P_out=P_rx,
    )
    recycle_mem = simulate_jt_expansion(
        ProcessStream(
            F_in={**_ZERO, 'A': tear_mem['A'], 'B': tear_mem['B']},
            T_in=T_mem, P_in=design.mem.P_H,
        ),
        P_out=P_rx,
    )

    # ---- Reactor 入口で合流 ----
    reactor_inlet = mix_streams([dist1_top_rx, recycle_dist3, recycle_mem])

    # ---- Step 2: Swing Reactor ----
    swing_feed = SwingFeed(
        F_in=reactor_inlet.F_in,
        T_feed=reactor_inlet.T_in,
        P_in=reactor_inlet.P_in,
    )
    r_rx = simulate_swing_reactor_system(design.swing, swing_feed, SwingFixed())

    # 設計判断 (2026-05-17): reactor が penalty 返却 (例: SV 範囲外、V_cat 異常等で
    # _penalty_result()) のとき、effluent は F=0, T=0, P=0。このまま下流の cooler/
    # compressor に流すと「P_in=0」で ValueError 発生 → solver.py の penalty_hit
    # 検査まで到達できない。早期に「penalty 状態で zero 流」のダミーを返して下流の
    # 全装置を penalty 結果でスキップする。solver 側で Reactor_CAPEX >=
    # PENALTY_CAPEX_THRESHOLD_OKUYEN をもって penalty_hit と判定する。
    # 設計判断 (2026-05-21): Reactor SV penalty 経路の連続シグナルを抽出。
    # 旧版は silent _penalty_result() で BO が「D を上下どちらに動かせば良いか」
    # 分からなかった。reactor_sv_shortfall を user_attr → constraints_func に渡す。
    reactor_shortfalls = _compute_reactor_shortfall(r_rx)
    if r_rx.equipment.Reactor_CAPEX >= PENALTY_CAPEX_THRESHOLD_OKUYEN:
        _log_unit_failure('Reactor (r_rx)', r_rx.equipment)
        result = _build_penalty_one_pass_result(
            r_rx, reactor_inlet, dist1_top_rx, recycle_dist3, recycle_mem,
            r1=r1, fresh=fresh, pump1=pump1,
        )
        result.update(reactor_shortfalls)
        return result

    rx_out = ProcessStream(
        F_in=r_rx.effluent.F_out_avg,
        T_in=r_rx.effluent.T_out_avg,
        P_in=r_rx.effluent.P_out,
    )

    # ---- Step 3: Cooler → Comp2 (2 段+段間冷却) → Dist2 ----
    # 反応器圧力 0.5 bar → Dist2 圧力 8.5 bar = 圧縮比 17:1。
    # 単段では断熱温度上昇が過大になる (~T_in×4) ため等圧縮比 √17≈4.12 の 2 段
    # に分割し、段間で T_intercool まで冷却して動力と機械的負荷を低減する。
    # 設計判断 (2026-05-09): Comp2 最終出口圧力は design.dist2.P_col に同期する
    # (operating.toml の comp2_out_Pa は backward compat のため残置だが未使用)。
    cooled = simulate_cooler(
        rx_out,
        T_out_target=config.temperature.cooler_after_reactor_K,
        process_phase=StreamPhase.GAS,
    )
    P_in_comp2  = cooled.outlet.P_in
    P_out_final = design.dist2.P_col
    P_mid       = math.sqrt(P_in_comp2 * P_out_final)
    T_intercool = config.temperature.cooler_after_reactor_K
    # 設計判断 (2026-05-09): Comp2b 出口 (~151°C) は Dist2 dew point (~50°C @ 8.5 bar)
    # を遥かに上回る超加熱蒸気。Dist2 partial condenser に直接入れると顕熱を冷凍冷媒で
    # 処理することになるため、工業実機の desuperheater (冷却水 HE) で dew 直上まで冷却。
    # Q ~12 MW を冷却水 (60 円/GJ) で除去 → 冷凍冷媒 (~1820 円/GJ) より遥かに安い。
    T_dist2_feed_K = 323.15   # 50°C: 8.5bar dew (~40-50°C) より少し上、5K margin
    with _capture_warnings("Comp2/Dist2", warnings_captured):
        comp2a    = simulate_compressor(cooled.outlet, P_out_target=P_mid)
        intercool = simulate_cooler(comp2a.outlet, T_out_target=T_intercool,
                                    process_phase=StreamPhase.GAS)
        comp2b    = simulate_compressor(intercool.outlet, P_out_target=P_out_final)
        desuper   = simulate_cooler(comp2b.outlet, T_out_target=T_dist2_feed_K,
                                    process_phase=StreamPhase.GAS)
        r2        = simulate_column2(desuper.outlet, tunables=design.dist2)

    # 設計判断 (2026-05-20): Dist2 rigorous (Wang-Henke) 収束失敗時 → penalty_result。
    # 下流 PSA / Mem は r2.top / r2.bottom がゼロ流量で組成計算が破綻するため早期 return。
    # dT_max_rigorous が equipment に格納されているので shortfall を TPE に伝える。
    if not getattr(r2.equipment, 'feasible', True):
        _log_unit_failure('Dist2 (r2)', r2.equipment)
        result = _build_penalty_after_column(
            'r2', design,
            fresh=fresh, pump1=pump1, r1=r1, dist1_top_rx=dist1_top_rx,
            reactor_inlet=reactor_inlet, r_rx=r_rx, rx_out=rx_out,
            cooled=cooled,
            comp2a=comp2a, intercool=intercool, comp2b=comp2b,
            desuper=desuper, r2=r2,
        )
        result['warnings_captured'] = warnings_captured
        # Reactor が成功していても shortfall フィールドは正常値 0 で伝播 (key 必須)
        result.update(reactor_shortfalls)
        return result

    # ---- Step 4: PSA ----
    # trace bypass (2026-05-19): PSA design モデルは C3 を扱えないため、入口の
    # C3 微量分 (≤ 1% of total) を design 計算から除き、後で offgas に合算する。
    # 2026-05-20: 閾値超過分 (= max_excess_frac) を BO penalty 用に runner.py へ伝播。
    psa_in_cleaned, psa_bypass, psa_trace_excess = _apply_trace_bypass(
        r2.top.F_in, _PSA_TRACE_COMPS, _TRACE_BYPASS_FRAC, label='PSA',
    )
    psa_feed = PSAFeedStream(
        F_in=psa_in_cleaned, T_in=r2.top.T_in, P_in=r2.top.P_in,
    )
    with _capture_warnings("PSA", warnings_captured):
        r_psa = simulate_psa_system(design.psa, psa_feed, PSAFixedParams())
    # PSA penalty を log (PDH_PER_UNIT_LOG=1 時のみ stderr)
    if getattr(r_psa.equipment, 'CAPEX_total', 0) >= PENALTY_CAPEX_THRESHOLD_OKUYEN:
        _log_unit_failure('PSA (r_psa)', r_psa.equipment)
    # 設計判断 (2026-05-21): PSA silent penalty 経路に対する連続 shortfall を抽出。
    # solver.py:191 が CAPEX sentinel で penalty_hit を判定するが、その情報のみだと
    # BO は「どう逃げれば良いか」分からない。psa_t_abs_shortfall 等を計算して
    # one_pass dict に積み、objective.py 経由で TPE constraints_func に届ける。
    psa_shortfalls = _compute_psa_shortfall(r_psa)
    # bypass 分を offgas に合算 (マスバランス保持)。r_psa.offgas は Dict[str,float]。
    if any(v > 0 for v in psa_bypass.values()):
        for c, v in psa_bypass.items():
            if v > 0:
                r_psa.offgas[c] = r_psa.offgas.get(c, 0.0) + v

    # ---- Step 5: Membrane ----
    # Dist2 を 8.5 bar 運転にした影響で塔底液の泡点が ~20°C まで下がり、冷却水
    # では液状態を保てない。膜の P_H <= 9.5 bar (Hua et al. 2024) を守るため、
    # 塔底液を mem_feed_K まで気化・過熱してガスフィードで膜へ送る。
    # 設計判断 (2026-05-08): 旧版は感熱のみで潜熱無視 → Mem 気化器 OPEX が 0
    # になっていた既知バグ。phase_change=True で潜熱を加算する。
    mem_precool = simulate_cooler(
        r2.bottom,
        T_out_target=config.temperature.mem_feed_K,
        phase_change=True,
        process_phase=StreamPhase.LIQUID,    # 顕熱区間: 液相加熱、潜熱区間は EVAPORATING に自動切替
    )
    # trace bypass (2026-05-19): Mem design モデルは C3H6/C3H8 二成分のみを扱う。
    # 上流 (Dist2 bot) に微量の non-C3 (H2/CH4/C2H4/C2H6) が混在する場合、
    # それを抽出して retentate (= recycle) に合算する。閾値超え時は warning。
    # 2026-05-20: 閾値超過分 (= max_excess_frac) を BO penalty 用に runner.py へ伝播。
    mem_in_cleaned, mem_bypass, mem_trace_excess = _apply_trace_bypass(
        mem_precool.outlet.F_in, _MEM_TRACE_COMPS, _TRACE_BYPASS_FRAC, label='Mem',
    )
    mem_feed = MemFeedStream(
        F_C3H6=mem_in_cleaned.get('B', 0.),
        F_C3H8=mem_in_cleaned.get('A', 0.),
        T_in=mem_precool.outlet.T_in,
        P_in=mem_precool.outlet.P_in,
    )
    with _capture_warnings("Mem", warnings_captured):
        r_mem = simulate_membrane_system(design.mem, mem_feed, MemFixedParams(vapor_feed=True))
    # Mem penalty を log (PDH_PER_UNIT_LOG=1 時のみ stderr)
    if getattr(r_mem.equipment, 'CAPEX_total', 0) >= PENALTY_CAPEX_THRESHOLD_OKUYEN:
        _log_unit_failure('Mem (r_mem)', r_mem.equipment)
    # 設計判断 (2026-05-22): Mem silent penalty 経路に対する連続 shortfall を抽出。
    # solver.py:191 が CAPEX sentinel で penalty_hit を判定するが、その情報のみだと
    # BO は「どう逃げれば良いか」分からない。mem_ph_shortfall / mem_bp_shortfall /
    # mem_phase_shortfall / mem_other_shortfall を one_pass dict に積み、
    # objective.py 経由で TPE constraints_func に届ける。
    mem_shortfalls = _compute_mem_shortfall(r_mem)
    # mem_bypass は recycle に合算 (= reactor 入口 mixer で扱われる)。tear_mem 構造は
    # 'A', 'B' しか持たないので、bypass 分は別チャネルで管理し reactor_inlet 直前で合流。
    # mem_bypass を後段で reactor_inlet に注入するため一旦保持。

    # ---- Step 6: Dist3 ----
    mem_to_dist3 = ProcessStream(
        F_in={'A': r_mem.product.F_C3H8, 'B': r_mem.product.F_C3H6,
              'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=r_mem.product.T_out, P_in=r_mem.product.P_out,
    )
    with _capture_warnings("Dist3", warnings_captured):
        r3 = simulate_column3(mem_to_dist3, tunables=design.dist3)

    # 設計判断 (2026-05-20): Dist3 penalty 早期検出。Dist3 失敗時は tear_dist3 が
    # ゼロ確定するので次反復で recycle_dist3 expansion が ValueError を再発する。
    if not getattr(r3.equipment, 'feasible', True):
        _log_unit_failure('Dist3 (r3)', r3.equipment)
        result = _build_penalty_after_column(
            'r3', design,
            fresh=fresh, pump1=pump1, r1=r1, dist1_top_rx=dist1_top_rx,
            reactor_inlet=reactor_inlet, r_rx=r_rx, rx_out=rx_out,
            cooled=cooled,
            comp2a=comp2a, intercool=intercool, comp2b=comp2b,
            desuper=desuper, r2=r2,
            r_psa=r_psa, mem_precool=mem_precool, r_mem=r_mem, r3=r3,
        )
        result['warnings_captured'] = warnings_captured
        result['trace_bypass_psa_excess'] = psa_trace_excess
        result['trace_bypass_mem_excess'] = mem_trace_excess
        # PSA / Reactor / Mem shortfall も伝播 (Dist3 失敗経路でも上流 penalty 情報を残す)
        result.update(psa_shortfalls)
        result.update(reactor_shortfalls)
        result.update(mem_shortfalls)
        return result

    # ---- tear stream の更新値 ----
    tear_dist3_new = {
        'A': r3.bottom.F_in.get('A', 0.0),
        'B': r3.bottom.F_in.get('B', 0.0),
    }
    tear_mem_new = {
        'A': r_mem.retentate.F_C3H8,
        'B': r_mem.retentate.F_C3H6,
    }
    T_d3_new  = r3.bottom.T_in
    T_mem_new = r_mem.retentate.T_out

    # 観測ラベル (2026-05-22): success path でも PSA/Mem は早期 return しないため、
    # 「Dist3 まで素通りしたが実は r_psa/r_mem に penalty_reason が刺さってる」trial が
    # ある (solver は CAPEX sentinel で後から penalty_hit を判定する経路)。
    # ここで第一失敗ユニットを特定し、診断ラベル群を user_attrs 用に出しておく。
    first_failed = _determine_first_failed_unit(
        r1=r1, r_rx=r_rx, r2=r2, r_psa=r_psa, r_mem=r_mem, r3=r3,
    )
    unit_diag = _extract_unit_diagnostics(
        r1=r1, r_rx=r_rx, r2=r2, r_psa=r_psa, r_mem=r_mem, r3=r3,
    )
    result = dict(
        pump1=pump1, r1=r1, dist1_top_rx=dist1_top_rx,
        reactor_inlet=reactor_inlet,
        r_rx=r_rx, rx_out=rx_out,
        cooled=cooled,
        comp2a=comp2a, intercool=intercool, comp2b=comp2b,
        desuper=desuper,
        r2=r2,
        r_psa=r_psa, mem_precool=mem_precool, r_mem=r_mem, r3=r3,
        tear_dist3_new=tear_dist3_new, tear_mem_new=tear_mem_new,
        T_d3_new=T_d3_new, T_mem_new=T_mem_new,
        warnings_captured=warnings_captured,
        trace_bypass_psa_excess=psa_trace_excess,
        trace_bypass_mem_excess=mem_trace_excess,
        # 正常完走時は全 0 (= TPE constraints_func で feasible 領域シグナル)
        dist1_N_shortfall=0.0, dist2_N_shortfall=0.0, dist3_N_shortfall=0.0,
        dist1_dT_shortfall=0.0, dist2_dT_shortfall=0.0, dist3_dT_shortfall=0.0,
        # PSA shortfall (正常完走時は全 0、penalty 経路でのみ > 0)
        **psa_shortfalls,
        # Reactor shortfall (正常完走時は全 0、SV 範囲外等の penalty 経路でのみ > 0)
        **reactor_shortfalls,
        # Mem shortfall (正常完走時は全 0、ph/bp/phase/other の penalty 経路で > 0)
        **mem_shortfalls,
        # 観測ラベル (2026-05-22)
        first_failed_unit=first_failed,
    )
    result.update(unit_diag)
    return result
