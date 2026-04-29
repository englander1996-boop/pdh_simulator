import win32com.client
import os
import time

HYSYS_FILE_PATH = r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc" 
COLUMN_NAME = "T-100" # ※もし実際の塔の名前が違えば修正してください

def test_hysys_column_access():
    print("1. HYSYSアプリケーションを起動中...")
    hysys_app = win32com.client.Dispatch("HYSYS.Application")
    hysys_app.Visible = True 
    case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
    
    try:
        flowsheet = case.Flowsheet
        
        print(f"2. 蒸留塔 '{COLUMN_NAME}' へのアクセスを試みます...")
        column = flowsheet.Operations.Item(COLUMN_NAME)
        # 蒸留塔内部のサブフローシートを取得
        col_fs = column.ColumnFlowsheet 
        
        print("3. 塔に設定されているスペック（仕様）一覧を取得します...")
        specs = col_fs.Specifications
        spec_names = []
        for i in range(specs.Count):
            spec_names.append(specs.Item(i).Name)
        print(f" -> 取得されたスペック名: {spec_names}")
        
        # --- 実験ポイント: ソルバーの停止 ---
        case.Solver.CanSolve = False
        
        # 「還流比」のスペック名を探す（デフォルトは "Reflux Ratio" のことが多い）
        target_spec = "Reflux Ratio" 
        
        if target_spec in spec_names:
            rr_spec = specs.Item(target_spec)
            current_rr = rr_spec.Goal.Value
            print(f"4. 現在の '{target_spec}' の値: {current_rr}")
            
            new_rr = current_rr + 0.1 # 少しだけ値を変更してみる
            print(f" -> 値を {new_rr} に書き込みます...")
            rr_spec.Goal.Value = new_rr
            
            print("5. ソルバーを再開し、収束を待ちます...")
            case.Solver.CanSolve = True
            
            start_time = time.time()
            converged = False
            while time.time() - start_time < 15:
                # 塔自体の収束状況も確認できるが、まずはケース全体のソルバーで判定
                if not case.Solver.IsSolving:
                    converged = True
                    break
                time.sleep(0.5)
                
            if converged:
                final_rr = specs.Item(target_spec).Goal.Value
                print(f" -> 成功！ 書き込み後の '{target_spec}': {final_rr}")
            else:
                print(" -> タイムアウト：計算が終わりませんでした。")
        else:
            print(f" -> [警告] '{target_spec}' という名前のスペックが見つかりません。")
            print(" ※HYSYS上の実際のスペック名を確認し、スクリプト内の target_spec を修正する必要があります。")

    except Exception as e:
        print(f" -> [エラー] 蒸留塔の操作中に失敗: {e}")
        
    finally:
        print("6. 検証終了。HYSYSの画面を確認してください。")

if __name__ == "__main__":
    test_hysys_column_access()