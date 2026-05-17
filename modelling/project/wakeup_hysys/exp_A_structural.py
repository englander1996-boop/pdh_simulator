"""
実験A: 段数・フィード段のCOM書き換え検証
目的:
  1. NumTrays / NumberOfStages の読み書きができるか
  2. FeedStage を複数のCOMパスで書き換えられるか
  3. Reset() 後にソルバーが収束するか
"""
import time
import win32com.client as win32
from pathlib import Path

HSC_PATH = str(Path(__file__).parent / "C3C4_Splitter_Base.hsc")
HYSYS_EMPTY = -32767.0

TARGET_STAGES    = 45   # 変更先の段数（元は40段）
TARGET_FEED_STAGE = 22  # 変更先のフィード段


def is_empty(v):
    try:
        return abs(float(v) - HYSYS_EMPTY) < 1.0
    except Exception:
        return True


def wait_convergence(dist, timeout=240):
    """dist.Temperature.Value が3回連続して安定したら収束とみなす。"""
    deadline = time.monotonic() + timeout
    ref = None
    stable = 0
    print("    収束待機中", end="", flush=True)
    while time.monotonic() < deadline:
        try:
            t = float(dist.Temperature.Value)
            if not is_empty(t):
                if ref is None:
                    ref = t
                if abs(t - ref) < 0.1:
                    stable += 1
                    if stable >= 3:
                        print(f"\n    → {t:.2f}°C で収束確認")
                        return True, t
                else:
                    stable = 0
                    ref = t
        except Exception:
            pass
        time.sleep(2)
        print(".", end="", flush=True)
    print("\n    → タイムアウト")
    return False, None


def try_read(obj, *props):
    """プロパティ名を順に試して最初に成功した値を返す。"""
    for p in props:
        try:
            return p, getattr(obj, p)
        except Exception:
            pass
    return None, None


def try_write(obj, value, *props):
    """プロパティ名を順に試して書き込む。成功したプロパティ名を返す。"""
    for p in props:
        try:
            setattr(obj, p, value)
            return p
        except Exception:
            pass
    return None


def main():
    sep = "=" * 60
    print(sep)
    print("実験A: 段数・フィード段 COM書き換え検証")
    print(sep)

    hysys = win32.Dispatch("HYSYS.Application")
    hysys.Visible = True
    case = hysys.SimulationCases.Open(HSC_PATH)
    print(f"ファイル: {HSC_PATH}\n")

    col    = case.Flowsheet.Operations.Item("T-100")
    col_fs = col.ColumnFlowsheet
    dist   = case.Flowsheet.MaterialStreams.Item("Distillate")

    try:
        main_tower = col_fs.Operations.Item("Main Tower")
    except Exception as e:
        print(f"[ERROR] 'Main Tower' が見つかりません: {e}")
        print("col_fs 内の機器一覧:")
        for i in range(col_fs.Operations.Count):
            op = col_fs.Operations.Item(i)
            print(f"  [{i}] {op.Name}")
        return

    # ─── [1] 現在の段数を読み取る ───────────────────────────────────
    print("[1] 段数読み取り")
    prop_stages, val_stages = try_read(main_tower, "NumberOfStages", "NumTrays")
    if prop_stages:
        print(f"    Main Tower.{prop_stages} = {val_stages}  ← 読み取りOK")
    else:
        print("    段数読み取り失敗（NumberOfStages / NumTrays 両方NG）")

    # ─── [2] 現在のフィード段を複数パスで読み取る ────────────────────
    print("\n[2] フィード段読み取り（複数パスを試す）")
    feed_stage_readers = [
        ("Main Tower.FeedStreams.Item('Feed').FeedStage",
         lambda: main_tower.FeedStreams.Item("Feed").FeedStage),
        ("Main Tower.FeedStreams.Item(0).FeedStage",
         lambda: main_tower.FeedStreams.Item(0).FeedStage),
        ("col.FeedConnections.Item(0).StepNumber",
         lambda: col.FeedConnections.Item(0).StepNumber),
        ("col.FeedConnections.Item(0).FeedStage",
         lambda: col.FeedConnections.Item(0).FeedStage),
    ]
    working_reader = None
    for path_name, reader in feed_stage_readers:
        try:
            val = reader()
            print(f"    OK  {path_name} = {val}")
            if working_reader is None:
                working_reader = (path_name, reader)
        except Exception as e:
            print(f"    NG  {path_name} → {e}")

    # ─── [3] ベースケース収束確認 ────────────────────────────────────
    print("\n[3] ベースケース収束確認（変更なし）")
    case.Solver.CanSolve = True
    ok, temp = wait_convergence(dist, timeout=120)
    print(f"    結果: {'Success' if ok else 'Timeout'}, Top_Temp={temp}")

    # ─── [4] 段数書き換えテスト ───────────────────────────────────────
    print(f"\n[4] 段数書き換えテスト: {val_stages} → {TARGET_STAGES}")
    case.Solver.CanSolve = False
    time.sleep(1)

    written_prop = try_write(main_tower, TARGET_STAGES, "NumberOfStages", "NumTrays")
    if written_prop:
        _, readback = try_read(main_tower, written_prop)
        changed = readback == TARGET_STAGES
        print(f"    書き込みプロパティ: {written_prop}")
        print(f"    読み戻し値: {readback}  → {'変更成功!' if changed else '変更失敗（値が変わっていない）'}")
    else:
        print("    段数書き込み失敗（全プロパティNG）")

    # ─── [5] フィード段書き換えテスト ────────────────────────────────
    print(f"\n[5] フィード段書き換えテスト: → {TARGET_FEED_STAGE}")
    feed_stage_writers = [
        ("Main Tower.FeedStreams.Item('Feed').FeedStage",
         lambda v: setattr(main_tower.FeedStreams.Item("Feed"), "FeedStage", v)),
        ("Main Tower.FeedStreams.Item(0).FeedStage",
         lambda v: setattr(main_tower.FeedStreams.Item(0), "FeedStage", v)),
        ("col.FeedConnections.Item(0).StepNumber",
         lambda v: setattr(col.FeedConnections.Item(0), "StepNumber", v)),
        ("col.FeedConnections.Item(0).FeedStage",
         lambda v: setattr(col.FeedConnections.Item(0), "FeedStage", v)),
    ]
    for path_name, writer in feed_stage_writers:
        try:
            writer(TARGET_FEED_STAGE)
            # 読み戻し
            readback = None
            if working_reader:
                try:
                    readback = working_reader[1]()
                except Exception:
                    pass
            print(f"    OK  {path_name}  読み戻し={readback}")
            break
        except Exception as e:
            print(f"    NG  {path_name} → {e}")

    # ─── [6] Reset + 収束確認 ─────────────────────────────────────────
    print(f"\n[6] Reset() → CanSolve=True → 収束待機（最大240秒）")
    print("    ※ コールドスタートのため時間がかかります")
    try:
        col.Reset()
        time.sleep(2)
    except Exception as e:
        print(f"    Reset() エラー: {e}")

    case.Solver.CanSolve = True
    ok, temp = wait_convergence(dist, timeout=240)
    print(f"    結果: {'Success' if ok else 'Timeout/NonConverged'}, Top_Temp={temp}")

    print(f"\n{sep}")
    print("実験A 終了")
    print(sep)
    try:
        case.Close(False)
        hysys.Quit()
    except Exception:
        pass


if __name__ == "__main__":
    main()
