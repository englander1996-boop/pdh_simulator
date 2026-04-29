"""
1ケースのみを実行する小テスト。
ハング切り分けのため、外側プロセスからタイムアウト管理される想定。
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import win32com.client as win32

import run_sampling as rs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--feed", type=float, required=True)
    parser.add_argument("--rr", type=float, required=True)
    parser.add_argument("--draw", type=float, required=True)
    parser.add_argument("--solver-timeout", type=int, default=120)
    parser.add_argument("--stall-timeout", type=int, default=25)
    args = parser.parse_args()

    rs.SOLVER_TIMEOUT_SEC = args.solver_timeout
    rs.SOLVER_STALL_ABORT_SEC = args.stall_timeout

    here = Path(__file__).parent
    hsc_path = str(here / "C3C4_Splitter_Base.hsc")

    row = pd.Series(
        {
            "LHS_Feed_Flow": args.feed,
            "HYSYS_Reflux_Ratio": args.rr,
            "HYSYS_Draw_Rate": args.draw,
        }
    )

    hysys_app = None
    case = None

    try:
        hysys_app = win32.Dispatch("HYSYS.Application")
        hysys_app.Visible = False
        case = hysys_app.SimulationCases.Open(hsc_path)
        comp_indices = rs.find_comp_indices(case)
        result = rs.run_hysys_simulation(hysys_app, case, row, comp_indices)

        payload = {
            "feed": args.feed,
            "rr": args.rr,
            "draw": args.draw,
            **result,
        }
        print("JSON_RESULT=" + json.dumps(payload, ensure_ascii=False))

    finally:
        if case is not None:
            try:
                case.Close(False)
            except Exception:
                pass
        if hysys_app is not None:
            try:
                hysys_app.Quit()
            except Exception:
                pass


if __name__ == "__main__":
    main()
