import win32com.client as win32
import time
import os

def test_change_trays():
    # 開きたいHYSYSファイルの絶対パス（rをつけることでエスケープシーケンスを無効化）
    file_path = r"Z:\pdh_simulator\project\wakeup_hysys\C3C4_Splitter_Base.hsc"
    
    # ファイルが存在するか念のためチェック
    if not os.path.exists(file_path):
        print(f"エラー: 指定されたファイルが見つかりません。\nパス: {file_path}")
        return

    print("HYSYSアプリケーションを起動・接続中...")
    try:
        # HYSYSアプリケーションを捕まえる（起動していない場合は新しく起動する）
        hyApp = win32.Dispatch("HYSYS.Application")
        hyApp.Visible = True  # HYSYSの画面を表示させる
        
        # 指定したファイルを直接開く
        print(f"ファイルを開いています: {file_path}")
        # ActiveDocumentではなく、パスを指定してOpenメソッドを使う
        hyCase = hyApp.SimulationCases.Open(file_path)
        
        # ファイルが開ききるまで少し待機
        time.sleep(3)
        
        # メインフローシートから蒸留塔「T-100」を取得
        col = hyCase.Flowsheet.Operations.Item("T-100")
        
        # 蒸留塔内部のサブフローシート環境にアクセス
        col_fs = col.ColumnFlowsheet
        
        # ★修正箇所：サブフローシート内の「塔本体（Main Tower）」を取得
        main_tower = col_fs.Operations.Item("Main Tower")
        
        # col_fsではなく、main_towerに対して段数を問い合わせる
        current_trays = main_tower.NumTrays
        print(f"接続成功！ 現在のT-100の段数: {current_trays}段")
        
        # 構造変更前にソルバーを停止
        print("ソルバーを一時停止します...")
        hyCase.Solver.CanSolve = False
        time.sleep(1)
        
        # ==========================================
        # テスト1：段数を変えてみる（40段 → 45段）
        # ==========================================
        new_trays = 45
        print(f"\nテスト開始: 段数を {new_trays} 段に変更します...")
        
        # ★修正箇所：main_towerに対して段数を書き換える
        main_tower.NumTrays = new_trays
        print("Pythonからの書き換え命令を送信しました。")
        time.sleep(1)
        
        # ==========================================
        # 発見ポイント：再計算前に「リセット」をかける
        # ==========================================
        print("塔のプロファイルをリセット（コールドスタート）します...")
        col.Reset()  # 蒸留塔のリセットボタンを押すのと同じ処理
        time.sleep(1)
        
        # ==========================================
        # テスト2：ソルバーを再起動してどうなるか観察
        # ==========================================
        print("ソルバーを再起動して計算を回します...")
        hyCase.Solver.CanSolve = True
        
        print("\n【観察ポイント】")
        print("HYSYSの画面を見てください。")
        print("1. T-100の段数は45段になっていますか？")
        print("2. 自動でリセットがかかり、無事に収束（緑色）しましたか？")
        
    except Exception as e:
        print(f"\nエラーが発生しました:\n{e}")
    finally:
        # 万が一エラーが起きてもソルバーはONに戻しておく
        try:
            hyCase.Solver.CanSolve = True
        except:
            pass

if __name__ == "__main__":
    test_change_trays()