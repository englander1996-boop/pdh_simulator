"""
段階的デバッグ用ミニテスト。
1) 固定条件を3回連続実行
2) 近傍9点(3x3)を実行

目的: 小さなテストを高速反復し、収束する安全領域を把握する。
"""

import sys
from pathlib import Path

import pandas as pd
import win32com.client as win32

import run_sampling as rs

if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).parent
HSC_PATH = str(HERE / "C3C4_Splitter_Base.hsc")


def build_row(feed_flow: float, reflux_ratio: float, draw_rate: float) -> pd.Series:
    return pd.Series(
        {
            "LHS_Feed_Flow": feed_flow,
            "HYSYS_Reflux_Ratio": reflux_ratio,
            "HYSYS_Draw_Rate": draw_rate,
        }
    )


def main() -> None:
    # 小テスト用に待機時間を短縮
    rs.SOLVER_TIMEOUT_SEC = 120
    rs.SOLVER_STALL_ABORT_SEC = 25

    hysys_app = None
    case = None

    fixed_plan = [
        (100.0, 2.0, 90.0),
        (100.0, 2.0, 90.0),
        (100.0, 2.0, 90.0),
    ]

    neighborhood_plan = [
        (95.0, 1.8, 85.0),
        (95.0, 2.0, 85.0),
        (95.0, 2.2, 85.0),
        (100.0, 1.8, 90.0),
        (100.0, 2.0, 90.0),
        (100.0, 2.2, 90.0),
        (105.0, 1.8, 95.0),
        (105.0, 2.0, 95.0),
        (105.0, 2.2, 95.0),
    ]

    all_cases = [("fixed", i + 1, c) for i, c in enumerate(fixed_plan)]
    all_cases += [("near", i + 1, c) for i, c in enumerate(neighborhood_plan)]

    records = []

    try:
        print("HYSYS接続中...", flush=True)
        hysys_app = win32.Dispatch("HYSYS.Application")
        hysys_app.Visible = True
        case = hysys_app.SimulationCases.Open(HSC_PATH)
        print(f"HYSYSファイルを開きました: {HSC_PATH}", flush=True)

        print("コンポーネントインデックス特定中...", flush=True)
        comp_indices = rs.find_comp_indices(case)

        print("\n--- ミニテスト開始 ---", flush=True)
        print(f"固定3件 + 近傍9件 = {len(all_cases)} 件", flush=True)

        for mode, idx, (feed, rr, dr) in all_cases:
            print("\n" + "=" * 64, flush=True)
            print(f"[{mode} #{idx}] Feed={feed:.2f}, RR={rr:.2f}, Draw={dr:.2f}", flush=True)
            row = build_row(feed, rr, dr)
            result = rs.run_hysys_simulation(hysys_app, case, row, comp_indices)
            status = result.get("Status", "Unknown")

            records.append(
                {
                    "mode": mode,
                    "index": idx,
                    "feed": feed,
                    "rr": rr,
                    "draw": dr,
                    **result,
                }
            )
            print(f"=> {status}", flush=True)

        df = pd.DataFrame(records)
        out_csv = HERE / "mini_test_results.csv"
        df.to_csv(out_csv, index=False)

        n_total = len(df)
        n_success = int((df["Status"] == "Success").sum())
        n_timeout = int((df["Status"] == "Timeout").sum())
        n_other = n_total - n_success - n_timeout

        print("\n--- 集計 ---", flush=True)
        print(f"総件数   : {n_total}", flush=True)
        print(f"Success  : {n_success}", flush=True)
        print(f"Timeout  : {n_timeout}", flush=True)
        print(f"Other    : {n_other}", flush=True)
        print(f"保存先   : {out_csv}", flush=True)

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
