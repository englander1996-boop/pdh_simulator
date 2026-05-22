r"""HYSYS COM 接続の最小スモークテスト。

目的: PDH シミュレーター本体を動かす前に、HYSYS が起動して HSC を開けるか
だけを切り分けて確認する。失敗した場合の原因を 4 段階に分けて報告する:

  [1] win32com.client.Dispatch("HYSYS.Application") が成功するか
      → HYSYS 自体がインストールされているか / COM 登録されているか
  [2] HYSYS.SimulationCases.Open() で HSC が開けるか
      → ライセンスが利用可能か / HSC ファイルが破損していないか
  [3] 蒸留塔 (T-100) と Spreadsheet (SPR-1) を取得できるか
      → HSC 内のオブジェクト名が想定通りか
  [4] Solver.IsSolving / Solver.CanSolve が読めるか
      → ソルバ制御 API が動くか

使い方: .\.venv\Scripts\python.exe tools\_smoke_test_hysys.py
"""
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from units.vle.hysys.registry import HsysRegistry


def _h(title: str) -> None:
    print(f"\n--- {title} ---")


def main() -> int:
    # ---- Step 0: registry 確認 ----
    _h("Step 0: HsysRegistry")
    try:
        reg = HsysRegistry()
        print(reg.summary())
    except Exception as e:
        print(f"NG: {e}")
        return 1

    # column1/30.hsc を対象に最小確認
    try:
        hsc_path = reg.get_path("column1", 30)
    except Exception as e:
        print(f"NG: get_path: {e}")
        return 2
    print(f"対象 HSC: {hsc_path}")

    # ---- Step 1: Dispatch ----
    _h("Step 1: win32com.client.Dispatch('HYSYS.Application')")
    try:
        import pythoncom
        import win32com.client as win32
        pythoncom.CoInitialize()
        t0 = time.time()
        app = win32.Dispatch("HYSYS.Application")
        print(f"OK: Dispatch 成功 ({time.time()-t0:.2f}s)")
        try:
            ver = app.Version
            print(f"   HYSYS Version = {ver}")
        except Exception:
            pass
    except Exception as e:
        print(f"NG: {type(e).__name__}: {e}")
        print(traceback.format_exc())
        return 3

    # ---- Step 2: SimulationCases.Open ----
    _h("Step 2: HSC オープン")
    case = None
    try:
        app.Visible = False
        t0 = time.time()
        case = app.SimulationCases.Open(str(hsc_path))
        print(f"OK: HSC オープン成功 ({time.time()-t0:.2f}s)")
        time.sleep(2.0)
    except Exception as e:
        print(f"NG: {type(e).__name__}: {e}")
        try:
            app.Quit()
        except Exception:
            pass
        return 4

    # ---- Step 3: 塔 / Spreadsheet 取得 ----
    _h("Step 3: T-100 / SPR-1 取得")
    try:
        fs = case.Flowsheet
        ops_count = fs.Operations.Count
        print(f"   Operations.Count = {ops_count}")
        col = fs.Operations.Item("T-100")
        ss  = fs.Operations.Item("SPR-1")
        feed = fs.MaterialStreams.Item("Feed")
        print(f"OK: T-100 / SPR-1 / Feed 取得成功")
        print(f"   feed.Pressure.Value = {feed.Pressure.Value}")
        print(f"   ss.Cell('A1').CellValue = {ss.Cell('A1').CellValue}")
    except Exception as e:
        print(f"NG: {type(e).__name__}: {e}")

    # ---- Step 4: Solver API ----
    _h("Step 4: Solver API")
    try:
        is_solving = case.Solver.IsSolving
        can_solve  = case.Solver.CanSolve
        print(f"OK: Solver.IsSolving={is_solving}, Solver.CanSolve={can_solve}")
    except Exception as e:
        print(f"NG: {type(e).__name__}: {e}")

    # ---- 終了処理 ----
    _h("終了処理")
    try:
        case.Close(False)
        app.Quit()
        print("OK: ケースクローズ + HYSYS Quit 成功")
    except Exception as e:
        print(f"クローズ時警告: {e}")

    print("\n==== スモークテスト完了 ====")
    return 0


if __name__ == "__main__":
    sys.exit(main())
