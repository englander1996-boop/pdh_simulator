"""
階層的サロゲートモデル構築用 HYSYSデータ収集スクリプト (Step 1〜3 統合版)
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import win32com.client as win32
from scipy.stats.qmc import LatinHypercube, scale

# Windowsコンソールでの日本語文字化け対策
if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

# ── 設定・定数 ─────────────────────────────────────────────────────────────────

_HERE = Path(__file__).parent
HSC_PATH = str(_HERE / "C3C4_Splitter_Base.hsc")
CSV_PATH = str(_HERE / "lhs_results.csv")

TARGET_SUCCESS = 80
SOLVER_TIMEOUT_SEC = 120
SOLVER_POLL_INTERVAL = 1.0

# 探索変数名（LHS生成順）
LHS_VARIABLE_NAMES = [
    "LHS_Total_Stages",
    "LHS_Feed_Stage_Ratio",
    "LHS_Feed_Pressure",
    "LHS_Feed_Flow",
    "LHS_Reflux_Ratio",
    "LHS_D_F_Ratio",
]

# 各変数の下限・上限
LOWER_BOUNDS = [20,   0.3, 1200,  70, 1.0, 0.85]
UPPER_BOUNDS = [80,   0.7, 2500, 130, 5.0, 0.95]

_RESULT_KEYS = [
    "Propane_Purity_Top",
    "Butane_Purity_Bottom",
    "Reboiler_Duty",
    "Condenser_Duty",
    "Top_Temperature",
    "Bottom_Temperature",
]

# ── Step 1: LHSサンプリング ───────────────────────────────────────────────────

def generate_lhs_conditions(n_samples: int = 300, seed: int = 42) -> pd.DataFrame:
    """
    LHSによるHYSYSシミュレーション条件一覧を生成する。
    """
    sampler = LatinHypercube(d=len(LHS_VARIABLE_NAMES), seed=seed)
    unit_samples = sampler.random(n=n_samples)

    # [0,1) → 実スケールへ変換
    scaled = scale(unit_samples, l_bounds=LOWER_BOUNDS, u_bounds=UPPER_BOUNDS)

    df = pd.DataFrame(scaled, columns=LHS_VARIABLE_NAMES)

    # HYSYS入力列の計算
    df["HYSYS_Total_Stages"] = df["LHS_Total_Stages"].round().astype(int)

    df["HYSYS_Feed_Stage"] = (
        df["HYSYS_Total_Stages"] * df["LHS_Feed_Stage_Ratio"]
    ).round().astype(int)

    df["HYSYS_Feed_Pressure"]      = df["LHS_Feed_Pressure"]
    df["HYSYS_Reboiler_Pressure"]  = df["LHS_Feed_Pressure"]
    df["HYSYS_Condenser_Pressure"] = df["LHS_Feed_Pressure"] - 20

    df["HYSYS_Reflux_Ratio"] = df["LHS_Reflux_Ratio"]
    df["HYSYS_Draw_Rate"]    = df["LHS_Feed_Flow"] * df["LHS_D_F_Ratio"]

    return df

# ── Step 2: HYSYS シミュレーション実行 ───────────────────────────────────────

def _get_specs_collection(column_flowsheet):
    """HYSYS のバージョン差を吸収して、スペック集合を返す。"""
    for attr in ("Specs", "Specifications"):
        try:
            specs = getattr(column_flowsheet, attr)
            if specs is not None:
                return specs
        except Exception:
            continue
    raise RuntimeError("ColumnFlowsheet から Specs / Specifications を取得できません")

def run_hysys_simulation(hysys_app, case, row_data: pd.Series) -> dict:
    """
    HYSYS に LHSサンプル1行分の条件を入力し、収束後に結果を取得して返す。
    """
    result = {key: np.nan for key in _RESULT_KEYS}
    result["Status"] = "NonConverged"

    try:
        # COMオブジェクト取得
        feed    = case.Flowsheet.MaterialStreams.Item("Feed")
        col     = case.Flowsheet.Operations.Item("T-100")
        main_ts = col.ColumnFlowsheet.Operations.Item("Main TS")
        cond    = col.ColumnFlowsheet.Operations.Item("Condenser")
        reboi   = col.ColumnFlowsheet.Operations.Item("Reboiler")
        specs   = _get_specs_collection(col.ColumnFlowsheet)

        dist    = case.Flowsheet.MaterialStreams.Item("Distillate")
        bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")
        qr      = case.Flowsheet.EnergyStreams.Item("Qr")
        qc      = case.Flowsheet.EnergyStreams.Item("Qc")

        # 変数入力前にメインソルバーを停止
        case.Solver.CanSolve = False

        # 1. Feed ストリームへの入力
        feed.Pressure.Value  = float(row_data["HYSYS_Feed_Pressure"])   # kPa
        feed.MolarFlow.Value = float(row_data["LHS_Feed_Flow"])         # kgmole/h (LHS元値)

        # 2. 蒸留塔 T-100 の構造パラメータ
        main_ts.NumberOfStages = int(row_data["HYSYS_Total_Stages"])
        feed_stage = int(row_data["HYSYS_Feed_Stage"])
        col.ColumnFlowsheet.Operations.Item("Main TS").FeedStreams.Item("Feed").FeedStage = feed_stage

        reboi.Pressure.Value = float(row_data["HYSYS_Reboiler_Pressure"])   # kPa
        cond.Pressure.Value  = float(row_data["HYSYS_Condenser_Pressure"])  # kPa

        # 3. スペック（Monitor画面）への入力
        specs.Item("Reflux Ratio").SpecValue = float(row_data["HYSYS_Reflux_Ratio"])
        specs.Item("Draw Rate").SpecValue    = float(row_data["HYSYS_Draw_Rate"])    # kgmole/h

        # ソルバー起動
        case.Solver.CanSolve = True

        # 収束待機（最大 SOLVER_TIMEOUT_SEC 秒）
        t_start = time.time()
        converged = False
        while (time.time() - t_start) < SOLVER_TIMEOUT_SEC:
            if col.ColumnFlowsheet.Solver.Converged:
                converged = True
                break
            time.sleep(SOLVER_POLL_INTERVAL)

        if not converged:
            case.Solver.CanSolve = False
            result["Status"] = "Timeout"
            return result

        # 結果取得
        try:
            propane_fractions = dist.ComponentMolarFractionValue
            result["Propane_Purity_Top"]   = propane_fractions[0]
        except Exception:
            result["Propane_Purity_Top"]   = dist.ComponentMolarFractions.Item("Propane").Value
        
        try:
            butane_fractions = bottoms.ComponentMolarFractionValue
            result["Butane_Purity_Bottom"] = butane_fractions[1]
        except Exception:
            result["Butane_Purity_Bottom"] = bottoms.ComponentMolarFractions.Item("n-Butane").Value
        
        result["Reboiler_Duty"]        = qr.HeatFlow.Value
        result["Condenser_Duty"]       = qc.HeatFlow.Value
        result["Top_Temperature"]      = dist.Temperature.Value
        result["Bottom_Temperature"]   = bottoms.Temperature.Value
        result["Status"] = "Success"

    except Exception as exc:
        result["Status"] = "NonConverged"
        print(f"  [HYSYS ERROR] {exc}", flush=True)

    return result

# ── Step 3: メインループ ───────────────────────────────────────────────────────

def main_sampling_loop():
    print("=" * 60)
    print("HYSYS データ収集メインループ開始")
    print("=" * 60)

    # 1. LHSデータフレーム生成
    print("\n[1] LHS条件の生成中...")
    df_conditions = generate_lhs_conditions(n_samples=300, seed=42)
    print(f"    生成完了: {len(df_conditions)} 件")
    print("    (先頭5行プレビュー)")
    print(df_conditions.head())

    if not os.path.exists(HSC_PATH):
        raise FileNotFoundError(f"HYSYSファイルが見つかりません: {HSC_PATH}")

    hysys_app = None
    case = None
    success_count = 0

    # CSVヘッダー書き込み要否
    write_header = not os.path.exists(CSV_PATH)

    print("\n[2] HYSYS接続...")
    try:
        hysys_app = win32.Dispatch("HYSYS.Application")
        hysys_app.Visible = True
        case = hysys_app.SimulationCases.Open(HSC_PATH)
        print(f"    ファイルを開きました: {HSC_PATH}")

        print("\n[3] ループ処理開始...")
        for idx, row in df_conditions.iterrows():
            print(f"\n--- 試行 {idx + 1}/{len(df_conditions)} ---")
            
            # シミュレーション実行
            result_dict = run_hysys_simulation(hysys_app, case, row)
            
            status = result_dict.get("Status")
            print(f"    結果 Status: {status}")

            if status == "Success":
                success_count += 1
                print(f"    成功回数: {success_count} / {TARGET_SUCCESS}")
                print(f"    塔頂温度: {result_dict.get('Top_Temperature'):.2f} °C, 塔底温度: {result_dict.get('Bottom_Temperature'):.2f} °C")

            # データの結合と保存
            row_result = pd.Series(result_dict)
            combined_row = pd.concat([row, row_result])
            combined_df = combined_row.to_frame().T
            
            combined_df.to_csv(CSV_PATH, mode='a', header=write_header, index=False)
            write_header = False

            # 目標回数に達したら終了
            if success_count >= TARGET_SUCCESS:
                print(f"\n目標成功回数 ({TARGET_SUCCESS}回) に到達しました。ループを終了します。")
                break

    except Exception as e:
        print(f"\n予期せぬエラーが発生しました: {e}")
    finally:
        print("\n[4] クリーンアップ...")
        if case is not None:
            try:
                case.Close(False)  # 保存せずに閉じる
                print("    ケースを閉じました。")
            except Exception as e:
                print(f"    ケースクローズ失敗: {e}")
        if hysys_app is not None:
            try:
                hysys_app.Quit()
                print("    HYSYSを終了しました。")
            except Exception as e:
                print(f"    HYSYS終了失敗: {e}")

if __name__ == "__main__":
    main_sampling_loop()