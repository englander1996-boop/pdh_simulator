import win32com.client as win32
import os

# 保存したファイルの絶対パスを指定してください
file_path = os.path.abspath(r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc")

try:
    # HYSYS アプリケーションを起動
    hysys = win32.Dispatch("HYSYS.Application")
    
    # ファイルを開く
    case = hysys.SimulationCases.Open(file_path)
    case.Visible = True  # HYSYSの画面を表示する
    
    print(f"無事に読み込めました: {file_path}")
    
    # 例：蒸留塔（T-100）にアクセスする場合
    # column = case.Flowsheet.Operations.Item("T-100")
    
except Exception as e:
    print(f"エラーが発生しました: {e}")