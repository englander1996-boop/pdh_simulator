"""
C3/C4蒸留塔サロゲートモデル構築用データ収集スクリプト。
LHSサンプリング → HYSYS実行 → CSV保存 を一括実行する。

実行方法:
    python run_sampling.py
"""

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


# ════════════════════════════════════════════════════════════════════════════════
# ステップ1: LHSサンプル生成
# ════════════════════════════════════════════════════════════════════════════════

LHS_VARIABLE_NAMES = [
    "LHS_Total_Stages",
    "LHS_Feed_Stage_Ratio",
    "LHS_Feed_Pressure",
    "LHS_Feed_Flow",
    "LHS_Reflux_Ratio",
    "LHS_D_F_Ratio",
]

LOWER_BOUNDS = [20,   0.3, 1200,  70, 1.0, 0.85]
UPPER_BOUNDS = [80,   0.7, 2500, 130, 5.0, 0.95]


def generate_lhs_conditions(n_samples: int = 300, seed: int = 42) -> pd.DataFrame:
    """LHSによるHYSYSシミュレーション条件一覧を生成する。"""
    sampler = LatinHypercube(d=len(LHS_VARIABLE_NAMES), seed=seed)
    unit_samples = sampler.random(n=n_samples)
    scaled = scale(unit_samples, l_bounds=LOWER_BOUNDS, u_bounds=UPPER_BOUNDS)
    df = pd.DataFrame(scaled, columns=LHS_VARIABLE_NAMES)

    df["HYSYS_Total_Stages"]       = df["LHS_Total_Stages"].round().astype(int)
    df["HYSYS_Feed_Stage"]         = (df["HYSYS_Total_Stages"] * df["LHS_Feed_Stage_Ratio"]).round().astype(int)
    df["HYSYS_Feed_Pressure"]      = df["LHS_Feed_Pressure"]
    df["HYSYS_Reboiler_Pressure"]  = df["LHS_Feed_Pressure"]
    df["HYSYS_Condenser_Pressure"] = df["LHS_Feed_Pressure"] - 20
    df["HYSYS_Reflux_Ratio"]       = df["LHS_Reflux_Ratio"]
    df["HYSYS_Draw_Rate"]          = df["LHS_Feed_Flow"] * df["LHS_D_F_Ratio"]

    return df


# ════════════════════════════════════════════════════════════════════════════════
# ステップ2: HYSYS 1ケース実行
# ════════════════════════════════════════════════════════════════════════════════

SOLVER_TIMEOUT_SEC   = 120
SOLVER_POLL_INTERVAL = 1.0

_RESULT_KEYS = [
    "Propane_Purity_Top",
    "Butane_Purity_Bottom",
    "Reboiler_Duty",
    "Condenser_Duty",
    "Top_Temperature",
    "Bottom_Temperature",
]


def run_hysys_simulation(hysys_app, case, row_data: pd.Series) -> dict:
    """HYSYSにLHSサンプル1行分の条件を入力し、収束後に結果を取得して返す。"""
    result = {key: np.nan for key in _RESULT_KEYS}
    result["Status"] = "NonConverged"

    try:
        feed    = case.Flowsheet.MaterialStreams.Item("Feed")
        col     = case.Flowsheet.Operations.Item("T-100")
        main_ts = col.ColumnFlowsheet.Operations.Item("Main TS")
        cond    = col.ColumnFlowsheet.Operations.Item("Condenser")
        reboi   = col.ColumnFlowsheet.Operations.Item("Reboiler")
        specs   = col.ColumnFlowsheet.Specs

        dist    = case.Flowsheet.MaterialStreams.Item("Distillate")
        bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")
        qr      = case.Flowsheet.EnergyStreams.Item("Qr")
        qc      = case.Flowsheet.EnergyStreams.Item("Qc")

        case.Solver.CanSolve = False

        # 1. Feed ストリームへの入力
        feed.Pressure.Value  = float(row_data["HYSYS_Feed_Pressure"])
        feed.MolarFlow.Value = float(row_data["LHS_Feed_Flow"])

        # 2. 蒸留塔 T-100 の構造パラメータ（塔段数 → フィード段 → 圧力の順）
        main_ts.NumberOfStages = int(row_data["HYSYS_Total_Stages"])
        feed_stage = int(row_data["HYSYS_Feed_Stage"])
        col.ColumnFlowsheet.Operations.Item("Main TS").FeedStreams.Item("Feed").FeedStage = feed_stage
        reboi.Pressure.Value = float(row_data["HYSYS_Reboiler_Pressure"])
        cond.Pressure.Value  = float(row_data["HYSYS_Condenser_Pressure"])

        # 3. スペック（Monitor画面）への入力
        specs.Item("Reflux Ratio").SpecValue = float(row_data["HYSYS_Reflux_Ratio"])
        specs.Item("Draw Rate").SpecValue    = float(row_data["HYSYS_Draw_Rate"])

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
        result["Propane_Purity_Top"]   = dist.ComponentMolarFractions.Item("Propane").Value
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


# ════════════════════════════════════════════════════════════════════════════════
# ステップ3: メインループ
# ════════════════════════════════════════════════════════════════════════════════

_HERE          = Path(__file__).parent
HSC_PATH       = str(_HERE / "C3C4_Splitter_Base.hsc")
OUTPUT_CSV     = _HERE / "lhs_results.csv"
TARGET_SUCCESS = 80
N_SAMPLES      = 300
SEED           = 42


def main_sampling_loop(
    hsc_path:       str  = HSC_PATH,
    output_csv:     Path = OUTPUT_CSV,
    n_samples:      int  = N_SAMPLES,
    target_success: int  = TARGET_SUCCESS,
    seed:           int  = SEED,
) -> None:
    """
    LHSサンプルを1行ずつHYSYSに投入し、収束・非収束を問わずCSVに追記する。
    target_success 件の "Success" が得られた時点でループを終了する。
    """
    df = generate_lhs_conditions(n_samples=n_samples, seed=seed)
    print(f"LHSサンプル生成完了: {len(df)} 件 (目標収束数: {target_success} 件)\n")

    hysys_app = None
    case      = None
    n_success = 0
    n_trials  = 0

    try:
        hysys_app = win32.Dispatch("HYSYS.Application")
        hysys_app.Visible = True
        case = hysys_app.SimulationCases.Open(hsc_path)
        print(f"HYSYSファイルを開きました: {hsc_path}\n")

        for _, row_data in df.iterrows():
            n_trials += 1
            print(f"[試行 {n_trials:>3}/{len(df)}] 実行中... ", end="", flush=True)

            sim_result = run_hysys_simulation(hysys_app, case, row_data)

            # 入力列と出力列を結合して1行に
            row_out = row_data.to_dict()
            row_out.update(sim_result)

            # アペンドモードでCSV保存（ファイル新規作成時のみヘッダーを書く）
            write_header = not output_csv.exists()
            pd.DataFrame([row_out]).to_csv(
                output_csv, mode="a", header=write_header, index=False
            )

            status = sim_result["Status"]
            if status == "Success":
                n_success += 1
                print(f"Success  ({n_success}/{target_success})", flush=True)
            else:
                print(status, flush=True)

            if n_success >= target_success:
                print(f"\n目標 {target_success} 件の収束データ取得完了。ループを終了します。")
                break

    finally:
        # ケースを保存せずに閉じ、HYSYSを終了
        if case is not None:
            try:
                case.Close(False)
            except Exception:
                pass
        if hysys_app is not None:
            try:
                hysys_app.Quit()
            except Exception:
                pass

    print(f"\n--- 完了 ---")
    print(f"  試行数   : {n_trials} 件")
    print(f"  収束     : {n_success} 件")
    print(f"  非収束等 : {n_trials - n_success} 件")
    print(f"  保存先   : {output_csv}")


if __name__ == "__main__":
    main_sampling_loop()
