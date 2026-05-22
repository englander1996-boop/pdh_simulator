"""Dist3 (C3 スプリッタ) の HYSYS アダプタ。

入力契約: Column3Input
主スペック: Draw Rate (フォールバック: Distillate Rate / 留出液流量)
流量単位 : kgmol/s で SPR-1 にも feed.MolarFlow にも書込み (column1/2 と違う)
書込順序: 圧力 → 流量・組成 → Draw Rate → フィード段 (最後)
"""
from __future__ import annotations

from typing import Optional

from units.vle.hysys.session import (
    HysysSession, HysysConvergenceError,
    set_feed_stage, run_column_and_wait,
)
from units.vle.hysys.adapters.base import (
    FEED_PRESSURE_MARGIN_KPA,
    write_pressure_spr1, write_feed_flow_kgmols,
    write_feed_pressure, write_feed_temperature,
    write_feed_composition_pdh, fill_outputs,
)
from units.vle.hysys.adapters.components import hysys_index_of
from units.vle.hysys.adapters.types import Column3Input, ColumnResult


class Column3Adapter:
    COLUMN_NAME = "T-100"
    SPREADSHEET_NAME = "SPR-1"
    STREAMS = {
        "feed":   "Feed",
        "top":    "Top",
        "bottom": "Bottom",
        "qc":     "Qc",
        "qr":     "Qr",
    }
    # lhs_column3.py / run_single_column3.py のフォールバック順 (旧運用、現在は使わない)
    DRAW_RATE_SPEC_NAMES = ("Draw Rate", "Distillate Rate", "留出液流量")
    # 設計判断 (2026-05-22): Draw Rate (絶対量 spec) は iter 変動に弱いので、
    # Comp Fraction (= 製品純度比率 spec) に切替え。HSC inspection で確認済み (Active=True)。
    PURITY_SPEC_NAMES = ("Comp Fraction",)

    @classmethod
    def run(
        cls,
        sess: HysysSession,
        spec: Column3Input,
        timeout_sec: float = 120.0,
        run_id: Optional[int] = None,
    ) -> ColumnResult:
        result = ColumnResult()
        if sess.popup_monitor is not None and run_id is not None:
            sess.popup_monitor.set_current_run(run_id)

        try:
            fs = sess.flowsheet
            ss = fs.Operations.Item(cls.SPREADSHEET_NAME)
            col = fs.Operations.Item(cls.COLUMN_NAME)
            feed = fs.MaterialStreams.Item(cls.STREAMS["feed"])
            top  = fs.MaterialStreams.Item(cls.STREAMS["top"])
            bot  = fs.MaterialStreams.Item(cls.STREAMS["bottom"])
            qc   = fs.EnergyStreams.Item(cls.STREAMS["qc"])
            qr   = fs.EnergyStreams.Item(cls.STREAMS["qr"])
            comps = sess.components()

            sess.solver_stop()

            # 1. 圧力 (lhs_column3.py 規約: feed.Pressure = 塔圧 + 50 kPa を必ず書込み)
            write_pressure_spr1(ss, spec.col_p_kpa)
            feed_p = (
                spec.feed_pressure_kpa
                if spec.feed_pressure_kpa is not None
                else spec.col_p_kpa + FEED_PRESSURE_MARGIN_KPA
            )
            write_feed_pressure(feed, feed_p)

            # 2. 流量 (kgmol/s 直書き)・組成
            write_feed_flow_kgmols(ss, feed, spec.feed_flow_kgmols)
            cls._write_composition(feed, comps, spec)
            if spec.feed_temperature_c is not None:
                write_feed_temperature(feed, spec.feed_temperature_c)

            # 3. スペック: exp1 rigorous と同じく "recovery × feed_LK" から Draw Rate を
            # 動的計算 (iter 変動に追従)。
            # 値が 0-1 範囲なら recovery_LK_top として解釈し、feed の C3H6 (LK='B') 流量
            # から Draw Rate [kgmol/s] を毎回計算。1.0 超なら旧来の絶対量 spec。
            spec_val = float(spec.draw_rate_kgmols)
            if 0.0 < spec_val <= 1.0:
                # feed の C3H6 (Propene) 流量 [kgmol/s] を取得
                j_propene = hysys_index_of(comps, "B")
                feed_flow_kgmols = float(feed.MolarFlow.Value)
                if j_propene is not None:
                    feed_fracs = list(feed.ComponentMolarFractionValue)
                    c3h6_kgmols = feed_flow_kgmols * float(feed_fracs[j_propene])
                else:
                    c3h6_kgmols = feed_flow_kgmols
                draw_rate = min(spec_val * c3h6_kgmols, feed_flow_kgmols * 0.99)
                draw_rate = max(draw_rate, 1e-6)
                import sys as _s
                _s.stderr.write(f"[Column3] feed_kgmols={feed_flow_kgmols:.4f}, "
                              f"c3h6_kgmols={c3h6_kgmols:.4f}, "
                              f"target_recovery={spec_val:.3f}, "
                              f"draw_rate_set={draw_rate:.4f} kgmol/s\n")
                _s.stderr.flush()
                _set_draw_rate(col, draw_rate, cls.DRAW_RATE_SPEC_NAMES)
            else:
                _set_draw_rate(col, spec_val, cls.DRAW_RATE_SPEC_NAMES)

            # 4. フィード段 (必ず最後)
            ok_feed = set_feed_stage(col, feed, int(spec.feed_stage))
            if not ok_feed:
                result.message = f"FeedStage {spec.feed_stage} 設定失敗 (続行)"

            # 5. ソルバ起動 + 収束待機
            sess.solver_start()
            converged = run_column_and_wait(col, sess.case, cls.COLUMN_NAME, timeout=timeout_sec)
            if not converged:
                raise HysysConvergenceError(
                    f"Column3 タイムアウト (timeout={timeout_sec}s)"
                )

            # 6. 結果回収
            fill_outputs(result, feed, top, bot, qc, qr, col, comps)

        except HysysConvergenceError as e:
            result.success = False
            result.message = str(e)
        except Exception as e:
            import traceback
            result.success = False
            result.message = f"Column3 COM エラー: {type(e).__name__}: {e}\n{traceback.format_exc()}"
        finally:
            if sess.popup_monitor is not None and run_id is not None:
                count, details = sess.popup_monitor.consume_run_events(run_id)
                result.warning_count = count
                result.warning_details = details

        return result

    # ----- 内部 -----
    @staticmethod
    def _write_composition(feed_stream, comps, spec: Column3Input) -> None:
        """組成書込み。propane_frac が指定されていれば Propane/Propene の 2 成分のみ
        (lhs_column3.py の運用)、そうでなければ全成分上書き。
        """
        if spec.propane_frac is not None:
            # Propane と Propene/Propylene のみ操作 (旧 lhs_column3.py 互換)
            propene_frac = 1.0 - float(spec.propane_frac)
            current = list(feed_stream.ComponentMolarFractionValue)
            j_propane = hysys_index_of(comps, "A")
            j_propene = hysys_index_of(comps, "B")
            if j_propane is not None:
                current[j_propane] = float(spec.propane_frac)
            if j_propene is not None:
                current[j_propene] = propene_frac
            s = sum(current)
            if s > 0:
                current = [x / s for x in current]
            feed_stream.ComponentMolarFractionValue = current
        elif spec.feed_composition:
            write_feed_composition_pdh(feed_stream, comps, spec.feed_composition)


def _set_draw_rate(col, value_kgmols: float, candidates) -> bool:
    """Draw Rate を順次フォールバックで設定。"""
    cf_specs = col.ColumnFlowsheet.Specifications
    for name in candidates:
        try:
            cf_specs.Item(name).GoalValue = float(value_kgmols)
            return True
        except Exception:
            pass
    return False


def _set_purity(col, value_frac: float, candidates) -> bool:
    """Comp Fraction (純度) spec を設定。"""
    cf_specs = col.ColumnFlowsheet.Specifications
    for name in candidates:
        try:
            cf_specs.Item(name).GoalValue = float(value_frac)
            return True
        except Exception:
            pass
    return False
