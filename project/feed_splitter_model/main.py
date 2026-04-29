import win32com.client
import pandas as pd
import numpy as np
from scipy.stats import qmc
import time
import os

# ==========================================
# 設定値
# ==========================================
HYSYS_FILE_PATH = r"Z:\pdh_simulator\project\feed_splitter_model\C3C4_Splitter_Base.hsc"
COLUMN_NAME = "T-100"
FEED_STREAM_NAME = "Feed"
CSV_FILE_PATH = "lhs_results_fixed_stages.csv"
TARGET_SUCCESS = 80
MAX_SAMPLES = 300

TOP_STREAM_NAME = "Distillate"
BTM_STREAM_NAME = "Bottoms"
PROPANE_IDX = 0 
BUTANE_IDX = 1

# ==========================================
# ステップ1: LHS生成（連続変数のみ）
# ==========================================
def generate_lhs_conditions(num_samples=300):
    bounds = {
        'LHS_Feed_Pressure': (1200, 2500), # kPa
        'LHS_Feed_Flow': (70, 130),        # kgmole/h
        'LHS_Reflux_Ratio': (1.0, 5.0),
        'LHS_D_F_Ratio': (0.85, 0.95)
    }
    
    keys = list(bounds.keys())
    lower_bounds = [bounds[k][0] for k in keys]
    upper_bounds = [bounds[k][1] for k in keys]
    
    sampler = qmc.LatinHypercube(d=len(keys))
    sample = sampler.random(n=num_samples)
    scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)
    
    df = pd.DataFrame(scaled_sample, columns=keys)
    
    # 実入力値の計算（段数関連は削除）
    df['HYSYS_Feed_Pressure'] = df['LHS_Feed_Pressure']
    df['HYSYS_Reboiler_Pressure'] = df['HYSYS_Feed_Pressure']
    df['HYSYS_Condenser_Pressure'] = df['HYSYS_Feed_Pressure'] - 20.0
    df['HYSYS_Reflux_Ratio'] = df['LHS_Reflux_Ratio']
    df['HYSYS_Draw_Rate'] = df['LHS_Feed_Flow'] * df['LHS_D_F_Ratio']
    
    return df

# ==========================================
# ステップ2: HYSYS実行（トポロジー変更なし）
# ==========================================
def run_hysys_simulation(case, row_data):
    flowsheet = case.Flowsheet
    column = flowsheet.Operations.Item(COLUMN_NAME)
    col_fs = column.ColumnFlowsheet
    
    result = {
        'Status': 'NonConverged',
        'Top_Propane_Purity': np.nan,
        'Btm_Butane_Purity': np.nan,
        'Reboiler_Duty': np.nan,
        'Condenser_Duty': np.nan,
        'Top_Temp': np.nan,
        'Btm_Temp': np.nan
    }
    
    try:
        case.Solver.CanSolve = False
        
        # 1. フィード条件入力
        feed_stream = flowsheet.MaterialStreams.Item(FEED_STREAM_NAME)
        feed_stream.Pressure.Value = row_data['HYSYS_Feed_Pressure']
        feed_stream.MolarFlow.Value = row_data['LHS_Feed_Flow']
        
        # 2. 圧力プロファイル入力
        try:
            col_fs.Operations.Item("Condenser").Pressure.Value = row_data['HYSYS_Condenser_Pressure']
            col_fs.Operations.Item("Reboiler").Pressure.Value = row_data['HYSYS_Reboiler_Pressure']
        except:
            pass 
            
        # 3. スペック入力
        specs = col_fs.Specifications
        specs.Item("Reflux Ratio").Goal.Value = row_data['HYSYS_Reflux_Ratio']
        specs.Item("Draw Rate").Goal.Value = row_data['HYSYS_Draw_Rate']
        
        # 4. ソルバー待機 (120秒)
        case.Solver.CanSolve = True
        start_time = time.time()
        while time.time() - start_time < 120:
            if not case.Solver.IsSolving:
                break
            time.sleep(1)
            
        # 5. 結果取得
        if col_fs.Converged:
            result['Status'] = 'Success'
            top_stream = flowsheet.MaterialStreams.Item(TOP_STREAM_NAME)
            btm_stream = flowsheet.MaterialStreams.Item(BTM_STREAM_NAME)
            
            result['Top_Propane_Purity'] = top_stream.ComponentMolarFractionValue(PROPANE_IDX)
            result['Btm_Butane_Purity'] = btm_stream.ComponentMolarFractionValue(BUTANE_IDX)
            result['Top_Temp'] = top_stream.Temperature.Value
            result['Btm_Temp'] = btm_stream.Temperature.Value
            
            try:
                result['Condenser_Duty'] = col_fs.Operations.Item("Condenser").HeatFlow.Value
                result['Reboiler_Duty'] = col_fs.Operations.Item("Reboiler").HeatFlow.Value
            except:
                pass
        else:
            result['Status'] = 'NonConverged'
            case.Solver.CanSolve = False 
            
    except Exception as e:
        print(f"  [Error] {e}")
        result['Status'] = 'Error'
        case.Solver.CanSolve = False
        
    return pd.Series(result)

# ==========================================
# ステップ3: メインループ
# ==========================================
def main_sampling_loop():
    print("1. LHSサンプリング生成...")
    df_lhs = generate_lhs_conditions(MAX_SAMPLES)
    
    print("2. HYSYS起動...")
    hysys_app = win32com.client.Dispatch("HYSYS.Application")
    hysys_app.Visible = True 
    case = None
    success_count = 0
    
    try:
        case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
        print(f" -> ケースオープン: {case.Title}")
        
        for idx, row in df_lhs.iterrows():
            print(f"\n--- 試行 {idx + 1}/{MAX_SAMPLES} ---")
            print(f"入力: P={row['HYSYS_Feed_Pressure']:.1f}, F={row['LHS_Feed_Flow']:.1f}, RR={row['HYSYS_Reflux_Ratio']:.2f}")
            
            result_series = run_hysys_simulation(case, row)
            print(f"結果: {result_series['Status']}")
            
            combined_row = pd.concat([row, result_series]).to_frame().T
            write_header = not os.path.exists(CSV_FILE_PATH)
            combined_row.to_csv(CSV_FILE_PATH, mode='a', header=write_header, index=False)
            
            if result_series['Status'] == 'Success':
                success_count += 1
                print(f"★ 成功: {success_count} / {TARGET_SUCCESS}")
                
            if success_count >= TARGET_SUCCESS:
                print("\n目標到達。終了します。")
                break
                
    except Exception as e:
        print(f"\n[致命的エラー] {e}")
    finally:
        if case is not None:
            case.Close() 
        hysys_app.Quit()
        print("HYSYS終了。")

if __name__ == "__main__":
    main_sampling_loop()