import win32com.client

HYSYS_FILE_PATH = r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc"
COLUMN_NAME = "T-100"

def explore_column_internals():
    print("1. HYSYSを起動してケースを開きます...")
    hysys_app = win32com.client.Dispatch("HYSYS.Application")
    hysys_app.Visible = True
    
    try:
        case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
        column = case.Flowsheet.Operations.Item(COLUMN_NAME)
        col_fs = column.ColumnFlowsheet
        
        print(f"2. 蒸留塔 '{COLUMN_NAME}' 内部の機器一覧を取得します...")
        op_count = col_fs.Operations.Count
        print(f" -> 内部機器の数: {op_count} 個")
        
        for i in range(op_count):
            op = col_fs.Operations.Item(i)
            # TypeNameが取得できないオブジェクトもあるため、try-exceptで囲む
            try:
                op_type = op.TypeName
            except:
                op_type = "Unknown"
                
            print(f"  [{i}] 名前: '{op.Name}' | 種類: {op_type}")
            
    except Exception as e:
        print(f" -> [エラー]: {e}")
    finally:
        print("3. 探索終了。")

if __name__ == "__main__":
    explore_column_internals()