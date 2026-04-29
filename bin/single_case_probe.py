"""
1ケースだけHYSYSに流して結果を返す最小プローブ。
状態汚染を避けるため毎回新規ケースで実行する。
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import win32com.client as win32

if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

_HERE = Path(__file__).parent
HSC_PATH = str(_HERE / "C3C4_Splitter_Base.hsc")
HYSYS_EMPTY_VALUE = -32767.0


def _is_empty(v) -> bool:
    try:
        return np.isnan(float(v)) or abs(float(v) - HYSYS_EMPTY_VALUE) < 1.0
    except Exception:
        return True


def set_spec(specs, name, value):
    spec = specs.Item(name)
    for attr in ("GoalValue", "Value", "TargetValue"):
        try:
            setattr(spec, attr, value)
            return attr
        except Exception:
            pass
    raise RuntimeError(f"{name} に値を設定できません")


def run_one(feed_value, rr_value, draw_value, timeout_sec=30):
    hysys = win32.Dispatch("HYSYS.Application")
    hysys.Visible = True
    case = hysys.SimulationCases.Open(HSC_PATH)
    try:
        feed = case.Flowsheet.MaterialStreams.Item("Feed")
        col = case.Flowsheet.Operations.Item("T-100")
        specs = col.ColumnFlowsheet.Specifications
        dist = case.Flowsheet.MaterialStreams.Item("Distillate")
        bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")

        case.Solver.CanSolve = False
        time.sleep(0.5)

        feed.MolarFlow.Value = feed_value
        rr_attr = set_spec(specs, "Reflux Ratio", rr_value)
        draw_attr = set_spec(specs, "Draw Rate", draw_value)

        case.Solver.CanSolve = True

        deadline = time.monotonic() + timeout_sec
        first_valid = None
        while time.monotonic() < deadline:
            try:
                temp = dist.Temperature.Value
                if not _is_empty(temp):
                    if first_valid is None:
                        first_valid = temp
                    if abs(float(temp) - float(first_valid)) < 0.1:
                        fracs_dist = dist.ComponentMolarFractionValue
                        fracs_bottom = bottoms.ComponentMolarFractionValue
                        return {
                            "status": "success",
                            "rr_attr": rr_attr,
                            "draw_attr": draw_attr,
                            "top_temp": round(float(dist.Temperature.Value), 4),
                            "bottom_temp": round(float(bottoms.Temperature.Value), 4),
                            "propane_top": round(float(fracs_dist[0]), 8),
                            "butane_bottom": round(float(fracs_bottom[1]), 8),
                        }
                    first_valid = temp
            except Exception as e:
                return {"status": f"error: {e}"}
            time.sleep(1)
        return {"status": "timeout", "rr_attr": rr_attr, "draw_attr": draw_attr}
    finally:
        try:
            case.Close(False)
        except Exception:
            pass
        try:
            hysys.Quit()
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feed", type=float, required=True)
    parser.add_argument("--rr", type=float, required=True)
    parser.add_argument("--draw", type=float, required=True)
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    print("=" * 60)
    print("single case probe")
    print("=" * 60)
    print(f"feed={args.feed}, rr={args.rr}, draw={args.draw}, timeout={args.timeout}")
    result = run_one(args.feed, args.rr, args.draw, args.timeout)
    print(result)


if __name__ == "__main__":
    main()
