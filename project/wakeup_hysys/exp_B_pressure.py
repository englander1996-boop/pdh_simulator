"""
実験B: フィード圧力と塔圧力プロファイルの関係検証
目的:
  1. 現在の塔圧力プロファイルを読み取る
  2. フィード圧力をフィード段推定圧力より低く設定するとHYSYSが止まるか確認
  3. 収束する圧力の範囲を把握する
"""
import time
import win32com.client as win32
from pathlib import Path

HSC_PATH = str(Path(__file__).parent / "C3C4_Splitter_Base.hsc")
HYSYS_EMPTY = -32767.0


def is_empty(v):
    try:
        return abs(float(v) - HYSYS_EMPTY) < 1.0
    except Exception:
        return True


def wait_convergence(dist, timeout=90):
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
    print("\n    → タイムアウト/未収束")
    return False, None


def reset_to_base(case, feed, base_p):
    """ベースケースの圧力に戻してからソルバーを起動し、収束させる。"""
    case.Solver.CanSolve = False
    time.sleep(0.5)
    feed.Pressure.Value = base_p
    case.Solver.CanSolve = True
    # 次のケースのため収束を待つ（簡易）
    time.sleep(15)
    case.Solver.CanSolve = False
    time.sleep(0.5)


def main():
    sep = "=" * 60
    print(sep)
    print("実験B: フィード圧力と塔圧力プロファイルの関係検証")
    print(sep)

    hysys = win32.Dispatch("HYSYS.Application")
    hysys.Visible = True
    case = hysys.SimulationCases.Open(HSC_PATH)
    print(f"ファイル: {HSC_PATH}\n")

    col    = case.Flowsheet.Operations.Item("T-100")
    col_fs = col.ColumnFlowsheet
    feed   = case.Flowsheet.MaterialStreams.Item("Feed")
    dist   = case.Flowsheet.MaterialStreams.Item("Distillate")

    # ─── [1] 現在の圧力プロファイルを読み取る ────────────────────────
    print("[1] 圧力プロファイルの読み取り")

    base_feed_p = float(feed.Pressure.Value)
    print(f"    Feed.Pressure           = {base_feed_p:.2f} kPa")

    cond_p = reb_p = None
    for name in ("Condenser", "Reboiler"):
        try:
            p = float(col_fs.Operations.Item(name).Pressure.Value)
            print(f"    {name}.Pressure      = {p:.2f} kPa")
            if name == "Condenser":
                cond_p = p
            else:
                reb_p = p
        except Exception as e:
            print(f"    {name}.Pressure      → エラー: {e}")

    # フィード段位置を取得して、フィード段の推定圧力を計算
    total_stages = feed_stage = None
    main_tower = None
    try:
        main_tower = col_fs.Operations.Item("Main Tower")
        for prop in ("NumberOfStages", "NumTrays"):
            try:
                total_stages = int(getattr(main_tower, prop))
                print(f"    Main Tower.{prop}    = {total_stages}")
                break
            except Exception:
                pass
    except Exception as e:
        print(f"    Main Tower 取得エラー: {e}")

    for getter_name, getter in [
        ("Main Tower.FeedStreams.Item('Feed').FeedStage",
         lambda: int(main_tower.FeedStreams.Item("Feed").FeedStage) if main_tower else None),
        ("col.FeedConnections.Item(0).StepNumber",
         lambda: int(col.FeedConnections.Item(0).StepNumber)),
    ]:
        try:
            feed_stage = getter()
            if feed_stage is not None:
                print(f"    フィード段              = {feed_stage}  (via {getter_name})")
                break
        except Exception:
            pass

    # フィード段推定圧力
    estimated_feed_tray_p = None
    if cond_p and reb_p and total_stages and feed_stage:
        estimated_feed_tray_p = cond_p + (feed_stage - 1) / (total_stages - 1) * (reb_p - cond_p)
        print(f"\n    ─── 推定フィード段圧力 = {estimated_feed_tray_p:.2f} kPa ───")
        print(f"    フィード圧との差       = {base_feed_p - estimated_feed_tray_p:+.2f} kPa")
    else:
        print("    フィード段推定圧力: 計算できません（段数またはフィード段が不明）")
        estimated_feed_tray_p = base_feed_p  # フォールバック

    # ─── [2] ベースケース確認 ─────────────────────────────────────────
    print(f"\n[2] ベースケース（Feed_P = {base_feed_p:.1f} kPa）")
    case.Solver.CanSolve = True
    ok, temp = wait_convergence(dist, timeout=90)
    print(f"    結果: {'Success' if ok else 'Timeout'}, Top_Temp={temp}")

    # ─── [3] フィード圧力を変化させて実験 ───────────────────────────
    # 推定フィード段圧力を基準に、下から上へと段階的にテスト
    test_cases = [
        ("推定フィード段圧 -200 kPa（明らかに低い）",  estimated_feed_tray_p - 200),
        ("推定フィード段圧 -100 kPa",                  estimated_feed_tray_p - 100),
        ("推定フィード段圧 -50 kPa",                   estimated_feed_tray_p -  50),
        ("推定フィード段圧 と同じ",                     estimated_feed_tray_p),
        ("推定フィード段圧 +50 kPa",                   estimated_feed_tray_p +  50),
        ("推定フィード段圧 +100 kPa（明らかに高い）",  estimated_feed_tray_p + 100),
        ("ベース圧力に戻す",                            base_feed_p),
    ]

    results = []
    for label, test_p in test_cases:
        print(f"\n[3] Feed_P = {test_p:.1f} kPa  ({label})")

        # ベース圧力に戻してから収束させる（前のケースの状態をリセット）
        reset_to_base(case, feed, base_feed_p)

        # テスト圧力を設定
        try:
            feed.Pressure.Value = test_p
            readback = float(feed.Pressure.Value)
            print(f"    設定値: {test_p:.1f} → 読み戻し: {readback:.1f} kPa")
        except Exception as e:
            print(f"    圧力設定エラー: {e}")
            results.append((label, test_p, "WriteError", None))
            continue

        case.Solver.CanSolve = True
        ok, temp = wait_convergence(dist, timeout=90)
        status = "Success" if ok else "Timeout/NonConverged"
        print(f"    結果: {status}, Top_Temp={temp}")
        results.append((label, test_p, status, temp))

    # ─── サマリー ─────────────────────────────────────────────────────
    print(f"\n{sep}")
    print("結果サマリー")
    print(sep)
    if estimated_feed_tray_p:
        print(f"  推定フィード段圧力: {estimated_feed_tray_p:.1f} kPa")
    print(f"  {'Feed_P [kPa]':>15}  {'Status':<20}  Top_Temp")
    print(f"  {'-'*55}")
    for label, p, status, temp in results:
        temp_str = f"{temp:.2f}°C" if temp is not None else "---"
        print(f"  {p:>15.1f}  {status:<20}  {temp_str}   ({label})")

    try:
        case.Close(False)
        hysys.Quit()
    except Exception:
        pass
    print(f"\n{sep}")
    print("実験B 終了")


if __name__ == "__main__":
    main()
