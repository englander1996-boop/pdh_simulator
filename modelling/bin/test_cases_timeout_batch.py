"""
ケースごとに別プロセスで実行し、タイムアウトで切り捨てるバッチ診断。
"""

import json
import subprocess
from pathlib import Path

import pandas as pd


def run_case(feed: float, rr: float, draw: float, case_timeout_sec: int = 45) -> dict:
    script = Path(__file__).parent / "test_single_case.py"
    cmd = [
        "python",
        str(script),
        "--feed",
        str(feed),
        "--rr",
        str(rr),
        "--draw",
        str(draw),
        "--solver-timeout",
        "120",
        "--stall-timeout",
        "25",
    ]

    try:
        cp = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=case_timeout_sec,
        )
        out = cp.stdout or ""
        marker = "JSON_RESULT="
        for line in out.splitlines()[::-1]:
            if line.startswith(marker):
                data = json.loads(line[len(marker) :])
                data["ProcessStatus"] = "Completed"
                data["ReturnCode"] = cp.returncode
                return data

        return {
            "feed": feed,
            "rr": rr,
            "draw": draw,
            "Status": "NoJsonResult",
            "ProcessStatus": "Completed",
            "ReturnCode": cp.returncode,
            "StdoutTail": "\n".join(out.splitlines()[-10:]),
            "StderrTail": "\n".join((cp.stderr or "").splitlines()[-10:]),
        }

    except subprocess.TimeoutExpired:
        return {
            "feed": feed,
            "rr": rr,
            "draw": draw,
            "Status": "ProcessTimeout",
            "ProcessStatus": "Timeout",
            "ReturnCode": -1,
        }


def main() -> None:
    cases = [
        (100.0, 2.0, 90.0),
        (95.0, 1.8, 85.0),
        (95.0, 2.0, 85.0),
        (95.0, 2.2, 85.0),
        (100.0, 1.8, 90.0),
        (100.0, 2.2, 90.0),
        (105.0, 1.8, 95.0),
        (105.0, 2.0, 95.0),
        (105.0, 2.2, 95.0),
    ]

    rows = []
    for i, (feed, rr, draw) in enumerate(cases, start=1):
        print(f"[{i}/{len(cases)}] feed={feed}, rr={rr}, draw={draw}", flush=True)
        row = run_case(feed, rr, draw, case_timeout_sec=45)
        rows.append(row)
        print(f"    -> Status={row.get('Status')} Process={row.get('ProcessStatus')}", flush=True)

    df = pd.DataFrame(rows)
    out_csv = Path(__file__).parent / "mini_timeout_batch_results.csv"
    df.to_csv(out_csv, index=False)

    print("\n--- 集計 ---")
    print(df[["feed", "rr", "draw", "Status", "ProcessStatus"]].to_string(index=False))
    print(f"\n保存先: {out_csv}")


if __name__ == "__main__":
    main()
