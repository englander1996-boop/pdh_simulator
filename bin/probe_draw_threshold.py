"""
draw の下限だけをざっくり確認する小テスト。
"""
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

CASES = [
    {"name": "draw_90", "feed": 100.0, "rr": 2.0, "draw": 90.0},
    {"name": "draw_88", "feed": 100.0, "rr": 2.0, "draw": 88.0},
    {"name": "draw_86", "feed": 100.0, "rr": 2.0, "draw": 86.0},
]


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


def run_case(case, params):
    feed = case.Flowsheet.MaterialStreams.Item("Feed")
    col = case.Flowsheet.Operations.Item("T-100")
    specs = col.ColumnFlowsheet.Specifications
    dist = case.Flowsheet.MaterialStreams.Item("Distillate")
    bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")

    case.Solver.CanSolve = False
    time.sleep(0.5)
    feed.MolarFlow.Value = params["feed"]
    rr_attr = set_spec(specs, "Reflux Ratio", params["rr"])
    draw_attr = set_spec(specs, "Draw Rate", params["draw"])
    case.Solver.CanSolve = True

    deadline = time.monotonic() + 30
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
            return {"status": f"error: {e}", "rr_attr": rr_attr, "draw_attr": draw_attr}
        time.sleep(1)
    return {"status": "timeout", "rr_attr": rr_attr, "draw_attr": draw_attr}


def main():
    print("=" * 60)
    print("draw下限プローブ")
    print("=" * 60)
    hysys = win32.Dispatch("HYSYS.Application")
    hysys.Visible = True
    case = hysys.SimulationCases.Open(HSC_PATH)
    print(f"opened: {HSC_PATH}")

    results = []
    try:
        for i, params in enumerate(CASES, 1):
            print(f"\n[{i}/{len(CASES)}] {params['name']} -> feed={params['feed']}, rr={params['rr']}, draw={params['draw']}")
            try:
                res = run_case(case, params)
            except Exception as e:
                res = {"status": f"error: {e}"}
            results.append((params["name"], res))
            print(f"    result: {res}")
    finally:
        try:
            case.Close(False)
        except Exception:
            pass
        try:
            hysys.Quit()
        except Exception:
            pass

    print("\nSUMMARY")
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
