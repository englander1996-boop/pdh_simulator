"""アダプタ共通ヘルパー。

3塔アダプタで重複する低レベル処理:
  - スプレッドシート (SPR-1) の A1/A2/A3 書込み (圧力 + 流量)
  - feed ストリームへの組成・温度・圧力の書込み
  - 還流比のスペック取得 (英/日/InternalRefluxRatio フォールバック)
  - 塔頂・塔底の基本出力 (温度・圧力・流量・組成) の回収
  - エンタルピー流量 [MW]
"""
from __future__ import annotations

from typing import Dict, Optional

from units.vle.hysys.adapters.components import (
    PDH_TO_HYSYS_KEYWORDS, build_pdh_to_hysys_index,
)
from units.vle.hysys.adapters.types import ColumnResult


# フィード圧力 = 塔圧 + 50 kPa という 3塔共通の規約。
FEED_PRESSURE_MARGIN_KPA = 50.0

# HYSYS の "Empty" sentinel value。Solver が未走 / 未収束のとき返る。
HYSYS_EMPTY_VALUE = -32767.0
_EMPTY_TOL = 1.0


def is_hysys_empty(value) -> bool:
    """HYSYS の empty sentinel を判定。Solver 結果が無効なら True。"""
    try:
        return abs(float(value) - HYSYS_EMPTY_VALUE) < _EMPTY_TOL
    except (TypeError, ValueError):
        return True


def write_pressure_spr1(ss, col_p_kpa: float) -> None:
    """SPR-1 の A1/A2 (= 凝縮器圧力 / リボイラ圧力) に塔圧を書込み。"""
    ss.Cell("A1").CellValue = float(col_p_kpa)
    ss.Cell("A2").CellValue = float(col_p_kpa)


def write_feed_flow_kgmolh(ss, feed_stream, flow_kgmolh: float) -> None:
    """SPR-1 A3 に kgmole/h、feed.MolarFlow に kgmole/s (= /3600) で書込み。

    column1, column2 用。column3 は kgmol/s で書く別関数を使う。
    """
    ss.Cell("A3").CellValue = float(flow_kgmolh)
    feed_stream.MolarFlow.Value = float(flow_kgmolh) / 3600.0


def write_feed_flow_kgmols(ss, feed_stream, flow_kgmols: float) -> None:
    """SPR-1 A3 にも kgmol/s 直書き、feed.MolarFlow にも同値。

    column3 用。
    """
    ss.Cell("A3").CellValue = float(flow_kgmols)
    feed_stream.MolarFlow.Value = float(flow_kgmols)


def write_feed_pressure(feed_stream, p_kpa: float) -> None:
    feed_stream.Pressure.Value = float(p_kpa)


def write_feed_temperature(feed_stream, t_c: float) -> None:
    feed_stream.Temperature.Value = float(t_c)


def write_feed_composition_pdh(feed_stream, comps, pdh_flows: Dict[str, float]) -> None:
    """PDH 成分流量 (kmol/h) を mol fraction に正規化して HYSYS feed に書込む。

    HYSYS Components 集合を index で走査し、PDH キーごとに対応する index を
    特定して current_fracs[idx] を上書きする。マッチしないキーは無視。
    """
    total = sum(max(v, 0.0) for v in pdh_flows.values())
    if total <= 0.0:
        return
    fractions = {k: max(v, 0.0) / total for k, v in pdh_flows.items()}

    current = list(feed_stream.ComponentMolarFractionValue)
    n = comps.Count
    for j in range(n):
        name = comps.Item(j).Name
        # PDH キーのどれかにマッチするか判定
        for pdh_key, kws in PDH_TO_HYSYS_KEYWORDS.items():
            if any(kw in name for kw in kws):
                current[j] = fractions.get(pdh_key, 0.0)
                break
        else:
            current[j] = 0.0
    # 数値誤差で合計が 1 から微妙にズレるので正規化
    s = sum(current)
    if s > 0:
        current = [x / s for x in current]
    feed_stream.ComponentMolarFractionValue = current


def get_reflux_ratio(col) -> Optional[float]:
    """還流比を取得 (英名 → 日本語名 → InternalRefluxRatio の順でフォールバック)。"""
    cf_specs = col.ColumnFlowsheet.Specifications
    for name in ("Reflux Ratio", "還流比"):
        try:
            return float(cf_specs.Item(name).CurrentValue)
        except Exception:
            pass
    try:
        return float(col.ColumnFlowsheet.InternalRefluxRatio(1))
    except Exception:
        return None


def enthalpy_flow_mw(stream) -> Optional[float]:
    """ストリームのエンタルピー流量 [MW]。

    MolarFlow [kgmol/s] × MolarEnthalpy [kJ/kgmol] / 1000 = MW
    取得に失敗したら None。
    """
    try:
        return float(stream.MolarFlow.Value) * float(stream.MolarEnthalpy.Value) / 1000.0
    except Exception:
        return None


def read_pdh_composition(stream, comps) -> Dict[str, float]:
    """HYSYS ストリームの mol fraction を PDH キー → 値 の辞書で返す。

    マッチしない成分があれば PDH キーから抜ける (キーが落ちる)。
    """
    fractions = list(stream.ComponentMolarFractionValue)
    idx_map = build_pdh_to_hysys_index(comps)
    return {k: float(fractions[i]) for k, i in idx_map.items()}


class HysysEmptyOutputError(RuntimeError):
    """HYSYS Solver が empty value を返した (= 未収束 / 未走)。"""


def fill_outputs(
    result: ColumnResult,
    feed_stream, top_stream, bot_stream,
    qc_stream, qr_stream,
    col, comps,
) -> None:
    """ColumnResult に基本出力 (Q, T, P, F, 組成, R, エンタルピー) を埋める。

    成功時は result.success=True にする。
    Solver が走らずに empty value (-32767) を返している場合は HysysEmptyOutputError を投げる。
    例外は呼び出し側で処理する想定 (raise する)。
    """
    # まず top_stream.Temperature が empty でないかを確認。
    # HYSYS は Solver 未走でも例外を出さず、参照すると EMPTY_VALUE (-32767) を返す。
    # これを feasible として扱うと下流の utility 選択で T=-32767°C で死ぬ。
    top_T = top_stream.Temperature.Value
    if is_hysys_empty(top_T):
        raise HysysEmptyOutputError(
            "top.Temperature が empty (Solver 未収束 / 未走)"
        )
    bot_T = bot_stream.Temperature.Value
    if is_hysys_empty(bot_T):
        raise HysysEmptyOutputError(
            "bottom.Temperature が empty (Solver 未収束 / 未走)"
        )
    top_F = top_stream.MolarFlow.Value
    bot_F = bot_stream.MolarFlow.Value
    if is_hysys_empty(top_F) or is_hysys_empty(bot_F):
        raise HysysEmptyOutputError(
            "top/bottom.MolarFlow が empty (Solver 未収束)"
        )
    qc_v = qc_stream.HeatFlow.Value
    qr_v = qr_stream.HeatFlow.Value
    if is_hysys_empty(qc_v) or is_hysys_empty(qr_v):
        raise HysysEmptyOutputError(
            "Qc / Qr が empty (Solver 未収束)"
        )

    result.qc_mw = float(qc_v) / 1000.0
    result.qr_mw = float(qr_v) / 1000.0

    result.top_flow_kmolh = float(top_F) * 3600.0
    result.top_temp_c = float(top_T)
    result.top_pressure_kpa = float(top_stream.Pressure.Value)
    result.top_composition = read_pdh_composition(top_stream, comps)

    result.bot_flow_kmolh = float(bot_F) * 3600.0
    result.bot_temp_c = float(bot_T)
    result.bot_pressure_kpa = float(bot_stream.Pressure.Value)
    result.bot_composition = read_pdh_composition(bot_stream, comps)

    result.reflux_ratio = get_reflux_ratio(col)
    result.feed_enthalpy_mw = enthalpy_flow_mw(feed_stream)
    result.top_enthalpy_mw = enthalpy_flow_mw(top_stream)
    result.bot_enthalpy_mw = enthalpy_flow_mw(bot_stream)

    result.success = True
