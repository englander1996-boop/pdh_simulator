"""
HYSYS COM API 診断スクリプト。
ベースケース値のまま1回だけ実行し、
- CanSolve=True の挙動（ブロッキングか否か）
- 収束検出方法
- ComponentMolarFraction のアクセス方法
を確認する。LHSループなし・入力変更なし。
"""
import sys
import time
import pythoncom
from pathlib import Path

import win32com.client as win32

if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

_HERE    = Path(__file__).parent
HSC_PATH = str(_HERE / "C3C4_Splitter_Base.hsc")


def main():
    print("=" * 60)
    print("HYSYS COM API 診断スクリプト")
    print("=" * 60)

    # ── 接続 ────────────────────────────────────────────────────
    print("\n[A] HYSYS 接続...")
    hysys_app = win32.Dispatch("HYSYS.Application")
    hysys_app.Visible = True
    case = hysys_app.SimulationCases.Open(HSC_PATH)
    print(f"    OK: {HSC_PATH}")

    # ── オブジェクト取得 ────────────────────────────────────────
    print("\n[B] オブジェクト取得...")
    feed    = case.Flowsheet.MaterialStreams.Item("Feed")
    col     = case.Flowsheet.Operations.Item("T-100")
    dist    = case.Flowsheet.MaterialStreams.Item("Distillate")
    bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")
    qr      = case.Flowsheet.EnergyStreams.Item("Qr")
    qc      = case.Flowsheet.EnergyStreams.Item("Qc")
    print("    OK")

    # ── CanSolve=False ──────────────────────────────────────────
    print("\n[C] CanSolve=False...")
    case.Solver.CanSolve = False
    time.sleep(1)
    print("    OK")

    # ── CanSolve=True（ブロッキング計測）───────────────────────
    print("\n[D] CanSolve=True（何秒かかるか計測）...")
    t0 = time.time()
    case.Solver.CanSolve = True
    elapsed = time.time() - t0
    print(f"    CanSolve=True から返るまで: {elapsed:.2f} 秒")
    if elapsed < 2:
        print("    → 非ブロッキング（HYSYSが非同期で計算中の可能性）")
    else:
        print("    → ブロッキング（HYSYS が計算してから返ってきた）")

    # ── 少し待って Temperature を確認 ──────────────────────────
    print("\n[E] 5秒待機後に Temperature を確認...")
    time.sleep(5)
    try:
        t = dist.Temperature.Value
        print(f"    dist.Temperature.Value = {t}")
        if abs(t - (-32767)) < 1:
            print("    → -32767（未計算）")
        else:
            print("    → 有効値！（収束済み）")
    except Exception as e:
        print(f"    FAILED: {e}")

    # さらに待つ
    for wait in [10, 30, 60]:
        print(f"\n[E2] さらに{wait}秒待機後に Temperature を確認...")
        time.sleep(wait)
        try:
            t = dist.Temperature.Value
            b = bottoms.Temperature.Value
            qr_val = qr.HeatFlow.Value
            print(f"    dist.Temperature   = {t}")
            print(f"    bottoms.Temperature = {b}")
            print(f"    Reboiler Duty      = {qr_val}")
            if abs(t - (-32767)) > 1:
                print("    → 収束を確認！以降の診断に進みます")
                break
        except Exception as e:
            print(f"    FAILED: {e}")
    else:
        print("\n    → 105秒待っても -32767。収束していないか、COM経由で取得できない。")

    # ── ComponentMolarFraction の診断 ───────────────────────────
    print("\n[F] ComponentMolarFraction の構造診断...")
    try:
        col_obj = dist.ComponentMolarFraction
        print(f"    type: {type(col_obj)}")

        # Count の確認
        try:
            n = col_obj.Count
            print(f"    Count = {n}")

            for i in range(n):
                try:
                    item = col_obj.Item(i)
                    try:
                        name = item.Name
                    except Exception:
                        name = "(名前なし)"
                    try:
                        val = item.Value
                    except Exception:
                        val = "(値取得失敗)"
                    print(f"    Item({i}): name={name!r}, value={val}")
                except Exception as ie:
                    print(f"    Item({i}) FAILED: {ie}")

        except Exception as ce:
            print(f"    .Count FAILED: {ce}")

            # SafeArray やタプルとして試す
            try:
                arr = list(col_obj)
                print(f"    list 変換成功: {arr[:5]}")
            except Exception as le:
                print(f"    list 変換 FAILED: {le}")

            # 直接インデックス
            for i in range(5):
                try:
                    val = col_obj[i]
                    print(f"    col_obj[{i}] = {val}")
                except Exception as ie:
                    print(f"    col_obj[{i}] FAILED: {ie}")
                    break

    except Exception as e:
        print(f"    ComponentMolarFraction FAILED: {e}")

    # ComponentMolarFractionValue の診断
    print("\n[G] ComponentMolarFractionValue の診断...")
    try:
        val_obj = dist.ComponentMolarFractionValue
        print(f"    type: {type(val_obj)}")
        try:
            arr = list(val_obj)
            print(f"    list 変換成功: {arr[:5]}")
        except Exception:
            pass
        for i in range(5):
            try:
                v = val_obj[i]
                print(f"    val_obj[{i}] = {v}")
            except Exception as ie:
                print(f"    val_obj[{i}] FAILED: {ie}")
                break
    except Exception as e:
        print(f"    ComponentMolarFractionValue FAILED: {e}")

    # ── ソルバー状態プロパティの全列挙 ─────────────────────────
    print("\n[H] case.Solver の利用可能プロパティ...")
    try:
        solver_attrs = sorted(a for a in dir(case.Solver) if not a.startswith("_"))
        print(f"    {solver_attrs}")
    except Exception as e:
        print(f"    FAILED: {e}")

    # ── 列ソルバーの確認 ──────────────────────────────────────
    print("\n[I] col.ColumnFlowsheet.Solver の確認...")
    try:
        col_solver = col.ColumnFlowsheet.Solver
        print(f"    type: {type(col_solver)}")
        print(f"    dir: {sorted([a for a in dir(col_solver) if not a.startswith('_')])}")
        if hasattr(col_solver, "Converged"):
            print(f"    Converged = {col_solver.Converged}")
        else:
            print(f"    Converged プロパティなし")
    except Exception as e:
        print(f"    FAILED: {e}")

    # ── スペック一覧の全列挙 ──────────────────────────────────
    print("\n[J] Column Specs の全列挙...")
    try:
        specs_dict = {
            "Specs": col.ColumnFlowsheet.Specs,
            "Specifications": col.ColumnFlowsheet.Specifications
        }
    except Exception as e:
        specs_dict = {}
        for name in ["Specs", "Specifications"]:
            try:
                specs_dict[name] = getattr(col.ColumnFlowsheet, name)
            except Exception:
                pass

    for location, specs in specs_dict.items():
        print(f"\n    [{location}]")
        if specs is not None:
            try:
                n = specs.Count
                print(f"    Count = {n}")
                for i in range(n):
                    try:
                        spec = specs.Item(i)
                        try:
                            spec_name = spec.Name
                        except Exception:
                            spec_name = "(名前なし)"
                        try:
                            spec_value = spec.SpecValue
                        except Exception:
                            spec_value = "(値取得失敗)"
                        print(f"    [{i}] {spec_name!r} = {spec_value}")
                    except Exception as e:
                        print(f"    [{i}] FAILED: {e}")
            except Exception as e:
                print(f"    FAILED to iterate: {e}")
        else:
            print(f"    None")

    # ── 終了 ────────────────────────────────────────────────────
    print("\n[Z] 後片付け（ファイルを閉じてHYSYSを終了）...")
    try:
        case.Close(False)
    except Exception:
        pass
    try:
        hysys_app.Quit()
    except Exception:
        pass
    print("    完了")


if __name__ == "__main__":
    main()
