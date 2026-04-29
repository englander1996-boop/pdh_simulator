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
CSV_FILE_PATH = "lhs_results.csv"
TARGET_SUCCESS = 80
MAX_SAMPLES = 300

# ※抽出対象のストリーム名は実際のHYSYSの設定に合わせて適宜変更してください
TOP_STREAM_NAME = "Distillate"
BTM_STREAM_NAME = "Bottoms"
# ※成分インデックスも実際のHYSYSのComponent Listに合わせてください（例: Propaneが0, Butaneが1など）
PROPANE_IDX = 0 
BUTANE_IDX = 1

# ==========================================
# ステップ1: LHSによるデータフレーム生成関数
# ==========================================
def generate_lhs_conditions(num_samples=300):
    # 探索変数の境界値 (Lower bounds, Upper bounds)
    bounds = {
        'LHS_Total_Stages': (20, 80),
        'LHS_Feed_Stage_Ratio': (0.3, 0.7),
        'LHS_Feed_Pressure': (1200, 2500), # kPa
        'LHS_Feed_Flow': (70, 130),        # kgmole/h
        'LHS_Reflux_Ratio': (1.0, 5.0),
        'LHS_D_F_Ratio': (0.85, 0.95)
    }
    
    keys = list(bounds.keys())
    lower_bounds = [bounds[k][0] for k in keys]
    upper_bounds = [bounds[k][1] for k in keys]
    
    # LHSサンプリング (scipy.stats.qmc)
    sampler = qmc.LatinHypercube(d=len(keys))
    sample = sampler.random(n=num_samples)
    scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)
    
    df = pd.DataFrame(scaled_sample, columns=keys)
    
    # 実入力値の計算と追加
    df['HYSYS_Total_Stages'] = df['LHS_Total_Stages'].round().astype(int)
    df['HYSYS_Feed_Stage'] = (df['HYSYS_Total_Stages'] * df['LHS_Feed_Stage_Ratio']).round().astype(int)
    df['HYSYS_Feed_Pressure'] = df['LHS_Feed_Pressure']
    df['HYSYS_Reboiler_Pressure'] = df['HYSYS_Feed_Pressure']
    df['HYSYS_Condenser_Pressure'] = df['HYSYS_Feed_Pressure'] - 20.0
    df['HYSYS_Reflux_Ratio'] = df['LHS_Reflux_Ratio']
    df['HYSYS_Draw_Rate'] = df['LHS_Feed_Flow'] * df['LHS_D_F_Ratio']
    
    return df

# ==========================================
# ステップ2: HYSYSシミュレーション実行関数
# ==========================================
def run_hysys_simulation(case, row_data):
    flowsheet = case.Flowsheet
    column = flowsheet.Operations.Item(COLUMN_NAME)
    col_fs = column.ColumnFlowsheet
    
    # 戻り値の初期化
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
        # 1. ソルバーを一時停止
        case.Solver.CanSolve = False
        
        # 2. フィード条件の入力
        feed_stream = flowsheet.MaterialStreams.Item(FEED_STREAM_NAME)
        feed_stream.Pressure.Value = row_data['HYSYS_Feed_Pressure']
        feed_stream.MolarFlow.Value = row_data['LHS_Feed_Flow']
        
        # 3. 塔の構造的変数の入力（トポロジー制約を考慮した安全な書き換え）
        main_tower = col_fs.Operations.Item("Main Tower")
        feed_conn = column.FeedConnections.Item(0)
        
        target_stages = int(row_data['HYSYS_Total_Stages'])
        target_feed_stage = int(row_data['HYSYS_Feed_Stage'])
        
        # バージョンごとのプロパティ名吸収
        try:
            current_stages = main_tower.NumTrays
            prop_name = "NumTrays"
        except:
            current_stages = main_tower.NumberOfStages
            prop_name = "NumberOfStages"
            
        # ★極めて重要：クラッシュ防止のための順序制御★
        if target_stages > current_stages:
            # 段数を増やす場合は、先に全体段数を増やしてからフィード段を下げる
            setattr(main_tower, prop_name, target_stages)
            feed_conn.StepNumber = target_feed_stage
        else:
            # 段数を減らす場合は、先にフィード段を上げてから全体段数を減らす
            feed_conn.StepNumber = target_feed_stage
            setattr(main_tower, prop_name, target_stages)
            
        # 4. 圧力プロファイルの入力
        condenser = col_fs.Operations.Item("Condenser")
        reboiler = col_fs.Operations.Item("Reboiler")
        # 機器自体の圧力設定（機器がサポートしていない場合はStream側に設定が必要）
        try:
            condenser.Pressure.Value = row_data['HYSYS_Condenser_Pressure']
            reboiler.Pressure.Value = row_data['HYSYS_Reboiler_Pressure']
        except:
            pass # すでに設定されているか、プロパティが異なる場合はスキップ
            
        # 5. スペック（設計仕様）の入力
        specs = col_fs.Specifications
        specs.Item("Reflux Ratio").Goal.Value = row_data['HYSYS_Reflux_Ratio']
        specs.Item("Draw Rate").Goal.Value = row_data['HYSYS_Draw_Rate']
        
        # 6. ソルバー再開と待機 (最大120秒)
        case.Solver.CanSolve = True
        
        start_time = time.time()
        while time.time() - start_time < 120:
            if not case.Solver.IsSolving:
                break
            time.sleep(1)
            
        # 7. 結果の判定と取得
        if column.ColumnFlowsheet.Converged:
            result['Status'] = 'Success'
            
            top_stream = flowsheet.MaterialStreams.Item(TOP_STREAM_NAME)
            btm_stream = flowsheet.MaterialStreams.Item(BTM_STREAM_NAME)
            
            # 純度取得 (ComponentMolarFractionValue のインデックスは要確認)
            result['Top_Propane_Purity'] = top_stream.ComponentMolarFractionValue(PROPANE_IDX)
            result['Btm_Butane_Purity'] = btm_stream.ComponentMolarFractionValue(BUTANE_IDX)
            
            # 温度取得
            result['Top_Temp'] = top_stream.Temperature.Value
            result['Btm_Temp'] = btm_stream.Temperature.Value
            
            # 熱負荷取得（コンデンサーとリボイラーのEnergy Stream、または機器のHeatFlowから）
            try:
                result['Condenser_Duty'] = condenser.HeatFlow.Value
                result['Reboiler_Duty'] = reboiler.HeatFlow.Value
            except:
                pass
                
        else:
            # 計算は終わったが収束していない場合
            result['Status'] = 'NonConverged'
            case.Solver.CanSolve = False # 次のイテレーションのためにリセット
            
    except Exception as e:
        print(f"  [Error] シミュレーション実行中に例外発生: {e}")
        result['Status'] = 'Error'
        case.Solver.CanSolve = False
        
    return pd.Series(result)

# ==========================================
# ステップ3: メインループ
# ==========================================
def main_sampling_loop():
    print("1. LHSサンプリングデータを生成中...")
    df_lhs = generate_lhs_conditions(MAX_SAMPLES)
    print(" -> 生成完了。先頭5行:")
    print(df_lhs.head())
    print("-" * 50)
    
    print("2. HYSYSを起動しています...")
    hysys_app = win32com.client.Dispatch("HYSYS.Application")
    # バックグラウンドで回す場合は False にしますが、デバッグ中は True が安全です
    hysys_app.Visible = True 
    
    case = None
    success_count = 0
    
    try:
        case = hysys_app.SimulationCases.Open(HYSYS_FILE_PATH)
        print(f" -> ケースオープン完了: {case.Title}")
        
        for idx, row in df_lhs.iterrows():
            print(f"--- 試行 {idx + 1}/{MAX_SAMPLES} ---")
            print(f"入力: Stages={row['HYSYS_Total_Stages']}, FeedStage={row['HYSYS_Feed_Stage']}, RR={row['HYSYS_Reflux_Ratio']:.2f}")
            
            # HYSYS関数の実行
            result_series = run_hysys_simulation(case, row)
            print(f"結果: {result_series['Status']}")
            
            # 行データと結果を結合
            combined_row = pd.concat([row, result_series]).to_frame().T
            
            # CSVへ追記保存 (ヘッダーは初回のみ)
            write_header = not os.path.exists(CSV_FILE_PATH)
            combined_row.to_csv(CSV_FILE_PATH, mode='a', header=write_header, index=False)
            
            if result_series['Status'] == 'Success':
                success_count += 1
                print(f"★ 成功カウント: {success_count} / {TARGET_SUCCESS}")
                
            if success_count >= TARGET_SUCCESS:
                print(f"\n目標成功回数 ({TARGET_SUCCESS}回) に到達しました。ループを終了します。")
                break
                
    except Exception as e:
        print(f"\n[致命的エラー] メインループ処理中にエラー発生: {e}")
        
    finally:
        print("\n3. クリーンアップ処理を実行中...")
        if case is not None:
            # 保存せずにケースを閉じる
            case.Close() 
        hysys_app.Quit()
        print(" -> HYSYSを終了しました。データ収集完了です！")

if __name__ == "__main__":
    main_sampling_loop()