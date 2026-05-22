"""Dist1 (脱ブタン塔) の HYSYS アダプタ。

入力契約: Column1Input (col_p_kpa, feed_flow_kmolh, comp_frac_2, feed_stage,
                       feed_composition[PDH キー → kmol/h])
主スペック: Comp Fraction - 2 (lhs_column1.py の運用に準拠)
書込順序 : 圧力 → 流量・組成 → Comp Fraction → フィード段 (最後)
"""
from __future__ import annotations

from typing import Optional

from units.vle.hysys.session import (
    HysysSession, HysysConvergenceError,
    set_feed_stage, run_column_and_wait,
)
from units.vle.hysys.adapters.base import (
    FEED_PRESSURE_MARGIN_KPA,
    write_pressure_spr1, write_feed_flow_kgmolh,
    write_feed_pressure, write_feed_temperature,
    write_feed_composition_pdh, fill_outputs,
)
from units.vle.hysys.adapters.types import Column1Input, ColumnResult


class Column1Adapter:
    COLUMN_NAME = "T-100"
    SPREADSHEET_NAME = "SPR-1"
    STREAMS = {
        "feed":   "Feed",
        "top":    "Top",
        "bottom": "Bottom",
        "qc":     "Qc",
        "qr":     "Qr",
    }
    SPEC_NAME = "Comp Fraction - 2"

    @classmethod
    def run(
        cls,
        sess: HysysSession,
        spec: Column1Input,
        timeout_sec: float = 120.0,
        run_id: Optional[int] = None,
    ) -> ColumnResult:
        """1ケース実行。HysysSession は open 済み前提。

        例外 (HYSYS COM エラー) はキャッチして ColumnResult.success=False で返す
        (上位 BO ループの停止を防ぐ)。
        """
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

            # ソルバ停止 → 入力 → 起動 のシーケンス
            sess.solver_stop()

            # 1. 圧力
            #  - SPR-1 A1/A2 に塔圧を書込み
            #  - feed.Pressure = 塔圧 + 50 kPa で必ず書込み (lhs_column1.py 規約)
            #    feed を書かないと HSC 現状値より塔圧が高い場合に逆流扱いで
            #    Solver が起動拒否する (column2 で実証)
            write_pressure_spr1(ss, spec.col_p_kpa)
            feed_p = (
                spec.feed_pressure_kpa
                if spec.feed_pressure_kpa is not None
                else spec.col_p_kpa + FEED_PRESSURE_MARGIN_KPA
            )
            write_feed_pressure(feed, feed_p)

            # 2. 流量・組成
            write_feed_flow_kgmolh(ss, feed, spec.feed_flow_kmolh)
            if spec.feed_composition:
                write_feed_composition_pdh(feed, comps, spec.feed_composition)
            if spec.feed_temperature_c is not None:
                write_feed_temperature(feed, spec.feed_temperature_c)

            # 3. スペック (Comp Fraction - 2)
            col.ColumnFlowsheet.Specifications.Item(cls.SPEC_NAME).GoalValue = float(spec.comp_frac_2)

            # 4. フィード段 (必ず最後)
            ok_feed = set_feed_stage(col, feed, int(spec.feed_stage))
            if not ok_feed:
                result.message = f"FeedStage {spec.feed_stage} 設定失敗 (続行)"

            # 5. ソルバ起動 + 収束待機
            sess.solver_start()
            converged = run_column_and_wait(col, sess.case, cls.COLUMN_NAME, timeout=timeout_sec)
            if not converged:
                raise HysysConvergenceError(
                    f"Column1 タイムアウト (timeout={timeout_sec}s)"
                )

            # 6. 結果回収
            fill_outputs(result, feed, top, bot, qc, qr, col, comps)

        except HysysConvergenceError as e:
            result.success = False
            result.message = str(e)
        except Exception as e:
            import traceback
            result.success = False
            result.message = f"Column1 COM エラー: {type(e).__name__}: {e}\n{traceback.format_exc()}"
        finally:
            if sess.popup_monitor is not None and run_id is not None:
                count, details = sess.popup_monitor.consume_run_events(run_id)
                result.warning_count = count
                result.warning_details = details

        return result
