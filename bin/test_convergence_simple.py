"""
シンプルな収束テスト: 固定条件1件だけ実行し、結果表示。
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

def _is_hysys_empty(v) -> bool:
    try:
        return np.isnan(float(v)) or abs(float(v) - HYSYS_EMPTY_VALUE) < 1.0
    except Exception:
        return True

def main():
    print("=" * 60)
    print("シンプル収束テスト")
    print("=" * 60)

    # [1] HYSYS接続
    print("\n[1] HYSYS接続...")
    hysys_app = win32.Dispatch("HYSYS.Application")
    hysys_app.Visible = True
    case = hysys_app.SimulationCases.Open(HSC_PATH)
    print("    OK")

    # [2] オブジェクト取得
    print("\n[2] オブジェクト取得...")
    feed = case.Flowsheet.MaterialStreams.Item("Feed")
    col = case.Flowsheet.Operations.Item("T-100")
    dist = case.Flowsheet.MaterialStreams.Item("Distillate")
    bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")
    specs = col.ColumnFlowsheet.Specifications
    print("    OK")

    # [3] ソルバー停止
    print("\n[3] ソルバー停止...")
    case.Solver.CanSolve = False
    time.sleep(1)
    print("    OK")

    # [4] 固定条件を設定
    print("\n[4] 入力設定...")
    feed.MolarFlow.Value = 100.0
    
    # スペック値設定を複数の属性で試す
    spec_rr = specs.Item("Reflux Ratio")
    spec_dr = specs.Item("Draw Rate")
    
    for attr in ("SpecValue", "Value", "GoalValue", "TargetValue"):
        try:
            setattr(spec_rr, attr, 2.0)
            print(f"    Reflux Ratio: {attr} で設定成功")
            break
        except Exception:
            pass
    
    for attr in ("SpecValue", "Value", "GoalValue", "TargetValue"):
        try:
            setattr(spec_dr, attr, 90.0)
            print(f"    Draw Rate: {attr} で設定成功")
            break
        except Exception:
            pass
    
    print("    Feed=100, RR=2.0, DrawRate=90")

    # [5] ソルバー起動
    print("\n[5] ソルバー起動...")
    case.Solver.CanSolve = True
    print("    CanSolve=True（非ブロッキング）")

    # [6] Temperature 有効値到達を待つ
    print("\n[6] Temperature 変化を監視...")
    deadline = time.time() + 180
    first_valid = None
    stable_count = 0

    while time.time() < deadline:
        try:
            temp = dist.Temperature.Value
            if not _is_hysys_empty(temp):
                if first_valid is None:
                    first_valid = temp
                    print(f"    ✓ 初回有効値: {temp:.2f} °C")
                if abs(temp - first_valid) < 0.1:
                    stable_count += 1
                    if stable_count >= 3:
                        print(f"    ✓ 安定収束: {temp:.2f} °C")
                        break
                else:
                    stable_count = 0
                    first_valid = temp
        except Exception as e:
            pass
        
        time.sleep(1)
        print(".", end="", flush=True)

    print("\n")

    # [7] 結果読み取り
    print("[7] 結果読み取り...")
    try:
        fracs_dist = dist.ComponentMolarFractionValue
        fracs_bottom = bottoms.ComponentMolarFractionValue
        print(f"    Propane_Purity_Top     = {fracs_dist[0]:.6f}")
        print(f"    Butane_Purity_Bottom   = {fracs_bottom[1]:.6f}")
        print(f"    Top_Temperature        = {dist.Temperature.Value:.2f}")
        print(f"    Bottom_Temperature     = {bottoms.Temperature.Value:.2f}")
        print("    ✓ Success!")
    except Exception as e:
        print(f"    ✗ FAILED: {e}")

    # [8] 終了
    print("\n[8] 終了...")
    try:
        case.Close(False)
        hysys_app.Quit()
    except Exception:
        pass
    print("    OK")

if __name__ == "__main__":
    main()
