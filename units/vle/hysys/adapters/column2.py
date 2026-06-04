"""Dist2 (脱エタン塔) の HYSYS アダプタ。

入力契約: Column2Input
主スペック: Reflux Ratio
ストリーム差分: top = "Topv" (column1/3 は "Top") ← partial cond で vapor distillate
書込順序: 圧力 → 流量・組成 → Reflux Ratio → フィード段 (最後)
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
    write_feed_composition_pdh, fill_outputs, is_hysys_empty,
)
from units.vle.hysys.adapters.types import Column2Input, ColumnResult


class Column2Adapter:
    COLUMN_NAME = "T-100"
    SPREADSHEET_NAME = "SPR-1"
    STREAMS = {
        "feed":   "Feed",
        "top":    "Topv",      # column2 だけ vapor distillate
        "bottom": "Bottom",
        "qc":     "Qc",
        "qr":     "Qr",
    }
    SPEC_NAME = "Reflux Ratio"

    # ---- cold分岐 安定化 (継続法 warm-start) ----
    # partial-condenser 脱エタン塔 (Reflux Ratio スペック) は多重定常解を持ち、HYSYS が
    # 初期値依存で warm 分岐 (塔頂 ~-85℃, 物理解) と cold 分岐 (塔頂 -128〜-150℃, 非物理) の
    # どちらにも収束しうる。同一入力でも run 毎に success↔fail が揺れる主因。対策として、
    # 目標 R で解いた塔頂が COLD_BRANCH_C より冷たい (= cold分岐) or 未収束なら、warm 解が
    # 一意な低 R (WARMSTART_REFLUX_MIN) から目標 R へ刻み (WARMSTART_RAMP_STEP) で継続収束させ、
    # warm 分岐に乗せ直す。warm 解が真に存在しない設計はそのまま cold/未収束で返る (= 真の infeasible)。
    COLD_BRANCH_C        = -105.0   # warm解は~-85℃。これより冷たい塔頂は cold分岐とみなす
    WARMSTART_REFLUX_MIN = 6.0      # 継続法の起点 (warm解が一意になりやすい低還流比)
    WARMSTART_RAMP_STEP  = 4.0      # 起点→目標へ上げる刻み

    @classmethod
    def run(
        cls,
        sess: HysysSession,
        spec: Column2Input,
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

            # 1. 圧力 (lhs_column2.py 規約: feed.Pressure = 塔圧 + 50 kPa を必ず書込み)
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

            # 3. スペック (Reflux Ratio)
            col.ColumnFlowsheet.Specifications.Item(cls.SPEC_NAME).GoalValue = float(spec.reflux_ratio)

            # 4. フィード段 (必ず最後)
            ok_feed = set_feed_stage(col, feed, int(spec.feed_stage))
            if not ok_feed:
                result.message = f"FeedStage {spec.feed_stage} 設定失敗 (続行)"

            # 5. ソルバ起動 + 収束待機
            sess.solver_start()
            converged = run_column_and_wait(col, sess.case, cls.COLUMN_NAME, timeout=timeout_sec)

            # 5b. cold分岐検出 + 継続法 warm-start リカバリ
            #     塔頂が cold分岐 (< COLD_BRANCH_C) か未収束なら、低 R から目標 R へ継続して
            #     warm 分岐に乗せ直す。warm 解が無い設計はそのまま未収束/cold で返る。
            def _t_top_c():
                try:
                    v = top.Temperature.Value
                    return None if is_hysys_empty(v) else float(v)
                except Exception:
                    return None

            tt = _t_top_c()
            R_target = float(spec.reflux_ratio)
            if ((not converged) or (tt is not None and tt < cls.COLD_BRANCH_C)) \
                    and cls.WARMSTART_REFLUX_MIN < R_target:
                spec_obj = col.ColumnFlowsheet.Specifications.Item(cls.SPEC_NAME)
                # 低 R → 目標 R の継続ラダー (warm 分岐を保ったまま昇順収束)
                ladder, r = [], cls.WARMSTART_REFLUX_MIN
                while r < R_target:
                    ladder.append(r)
                    r += cls.WARMSTART_RAMP_STEP
                ladder.append(R_target)
                for r_step in ladder:
                    sess.solver_stop()
                    spec_obj.GoalValue = float(r_step)
                    sess.solver_start()
                    converged = run_column_and_wait(
                        col, sess.case, cls.COLUMN_NAME, timeout=timeout_sec)
                tt2 = _t_top_c()
                note = (f"cold分岐検出(T_top={tt}°C)→継続法warm-start "
                        f"(R:{cls.WARMSTART_REFLUX_MIN}→{R_target}, 復帰後T_top={tt2}°C)")
                result.message = (result.message + " | " + note) if result.message else note

            if not converged:
                raise HysysConvergenceError(
                    f"Column2 タイムアウト/未収束 (timeout={timeout_sec}s)"
                )

            # 6. 結果回収
            fill_outputs(result, feed, top, bot, qc, qr, col, comps)

        except HysysConvergenceError as e:
            result.success = False
            result.message = str(e)
        except Exception as e:
            import traceback
            result.success = False
            result.message = f"Column2 COM エラー: {type(e).__name__}: {e}\n{traceback.format_exc()}"
        finally:
            if sess.popup_monitor is not None and run_id is not None:
                count, details = sess.popup_monitor.consume_run_events(run_id)
                result.warning_count = count
                result.warning_details = details

        return result
