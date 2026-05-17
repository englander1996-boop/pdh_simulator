import win32com.client
import os

# テスト用設定：実際のパスに書き換えてください
HYSYS_FILE_PATH = r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc" 
FEED_STREAM_NAME = "Feed" # 実際のフィードストリーム名に合わせてください

def test_hysys_connection():
    print("1. HYSYSアプリケーションを起動中...")
    try:
        # 起動と可視化（裏で動かさず、必ず画面に出す）
        hysys_app = win32com.client.Dispatch("HYSYS.Application")
        hysys_app.Visible = True 
        print(" -> アプリケーションの起動に成功しました。")
    except Exception as e:
        print(f" -> [エラー] HYSYSの起動に失敗: {e}")
        return

    print(f"2. ケースファイルを開いています: {HYSYS_FILE_PATH}")
    try:
        # ファイルの存在確認
        if not os.path.exists(HYSYS_FILE_PATH):
            raise FileNotFoundError(f"ファイルが見つかりません: {HYSYS_FILE_PATH}")
            
        case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
        print(f" -> ケースが開かれました: {case.Title}")
    except Exception as e:
        print(f" -> [エラー] ケースを開くのに失敗: {e}")
        return

    print(f"3. ストリーム '{FEED_STREAM_NAME}' の圧力を取得中...")
    try:
        # フローシートオブジェクトの取得
        flowsheet = case.Flowsheet
        
        # ストリームへのアクセス
        feed_stream = flowsheet.MaterialStreams.Item(FEED_STREAM_NAME)
        
        # 圧力の取得（HYSYSの内部単位系に注意）
        # .Value で数値を取得します
        pressure_value = feed_stream.Pressure.Value
        
        print(f" -> 取得成功！ フィード圧力: {pressure_value}")
        
    except Exception as e:
        print(f" -> [エラー] 変数の読み取りに失敗: {e}")
    
    finally:
        print("4. クリーンアップを行わずに終了します（HYSYSの画面を確認するため）。")
        # 今回はわざと閉じません。HYSYSがどんな状態か目視するためです。

if __name__ == "__main__":
    test_hysys_connection()