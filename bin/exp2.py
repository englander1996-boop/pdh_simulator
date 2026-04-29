import win32com.client
import os
import time

# 実験1で成功したパスを使用
HYSYS_FILE_PATH = r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc" 
FEED_STREAM_NAME = "Feed"

def test_hysys_write_standalone():
    print("1. HYSYSアプリケーションを起動中...")
    try:
        hysys_app = win32com.client.Dispatch("HYSYS.Application")
        hysys_app.Visible = True 
        print(" -> 起動成功")
    except Exception as e:
        print(f" -> [エラー] 起動失敗: {e}")
        return

    print("2. ケースファイルを開いています...")
    try:
        case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
        print(f" -> ケースが開かれました: {case.Title}")
    except Exception as e:
        print(f" -> [エラー] ケースを開くのに失敗: {e}")
        return
    
    try:
        flowsheet = case.Flowsheet
        feed_stream = flowsheet.MaterialStreams.Item(FEED_STREAM_NAME)
        
        # --- 実験ポイント: ソルバーの停止 ---
        print("3. ソルバーを一時停止します...")
        case.Solver.CanSolve = False
        
        # --- 実験ポイント: 書き込み ---
        new_pressure = 1600.0
        print(f"4. 圧力を {new_pressure} に書き込み中...")
        feed_stream.Pressure.Value = new_pressure
        
        # --- 実験ポイント: ソルバー再開と収束待ち ---
        print("5. ソルバーを再開し、収束を待ちます...")
        case.Solver.CanSolve = True
        
        # 収束待ち
        start_time = time.time()
        converged = False
        while time.time() - start_time < 10:
            if not case.Solver.IsSolving:
                converged = True
                break
            time.sleep(0.5)
            
        if converged:
            final_p = feed_stream.Pressure.Value
            print(f" -> 成功！ 書き込み後の圧力: {final_p}")
        else:
            print(" -> タイムアウト：計算が終わりませんでした。")

    except Exception as e:
        print(f" -> [エラー] 書き込み処理中に失敗: {e}")
        
    finally:
        print("6. 検証のため、HYSYSはそのままにして終了します。")

if __name__ == "__main__":
    test_hysys_write_standalone()