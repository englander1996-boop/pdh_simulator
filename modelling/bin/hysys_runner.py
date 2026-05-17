"""
HYSYS COMインターフェース経由でシミュレーションを1ケース実行する関数モジュール。
lhs_sampling.py で生成した pd.Series (1行) を受け取り、結果を dict で返す。
"""

import time
import numpy as np
import pandas as pd

# ── 定数 ───────────────────────────────────────────────────────────────────────

SOLVER_TIMEOUT_SEC    = 120   # 収束待機の上限秒数
SOLVER_POLL_INTERVAL  = 1.0   # 収束確認ポーリング間隔 (秒)

# 出力キー一覧（Status を除く数値列）
_RESULT_KEYS = [
    "Propane_Purity_Top",
    "Butane_Purity_Bottom",
    "Reboiler_Duty",
    "Condenser_Duty",
    "Top_Temperature",
    "Bottom_Temperature",
]


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


# ── メイン関数 ──────────────────────────────────────────────────────────────────

def run_hysys_simulation(hysys_app, case, row_data: pd.Series) -> dict:
    """
    HYSYS に LHSサンプル1行分の条件を入力し、収束後に結果を取得して返す。

    Parameters
    ----------
    hysys_app : HYSYS.Application COM object
        メインループ側で Dispatch した HYSYS アプリケーションオブジェクト
    case : HYSYS SimulationCase COM object
        hysys_app.SimulationCases.Open() で取得したケースオブジェクト
    row_data : pd.Series
        generate_lhs_conditions() が返すデータフレームの1行
        (LHS_* 列と HYSYS_* 列の両方を含む)

    Returns
    -------
    dict
        キー: _RESULT_KEYS の各出力値 + "Status"
        Status: "Success" | "NonConverged" | "Timeout"
        未収束・エラー時は数値列がすべて np.nan
    """

    # 初期値は全 NaN (後で上書き)
    result = {key: np.nan for key in _RESULT_KEYS}
    result["Status"] = "NonConverged"

    try:
        # ── COMオブジェクト取得 ────────────────────────────────────────────────
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

        # ── 変数入力前にメインソルバーを停止 ──────────────────────────────────
        case.Solver.CanSolve = False

        # ── 1. Feed ストリームへの入力 ─────────────────────────────────────────
        feed.Pressure.Value  = float(row_data["HYSYS_Feed_Pressure"])   # kPa
        feed.MolarFlow.Value = float(row_data["LHS_Feed_Flow"])         # kgmole/h (LHS元値)

        # ── 2. 蒸留塔 T-100 の構造パラメータ ──────────────────────────────────
        # 塔段数を先に変更してからフィード段を設定する（順序依存あり）
        main_ts.NumberOfStages = int(row_data["HYSYS_Total_Stages"])
        # フィード段（入力段）の変更
        # 注: HYSYSの段番号は1から始まるため、HYSYS_Feed_Stage は 1 以上 NumberOfStages 以下である前提
        feed_stage = int(row_data["HYSYS_Feed_Stage"])
        col.ColumnFlowsheet.Operations.Item("Main TS").FeedStreams.Item("Feed").FeedStage = feed_stage

        # 圧力はリボイラ → コンデンサの順に設定（逆流エラー防止）
        reboi.Pressure.Value = float(row_data["HYSYS_Reboiler_Pressure"])   # kPa
        cond.Pressure.Value  = float(row_data["HYSYS_Condenser_Pressure"])  # kPa

        # ── 3. スペック（Monitor画面）への入力 ────────────────────────────────
        specs.Item("Reflux Ratio").SpecValue = float(row_data["HYSYS_Reflux_Ratio"])
        specs.Item("Draw Rate").SpecValue    = float(row_data["HYSYS_Draw_Rate"])    # kgmole/h

        # ── ソルバー起動 ───────────────────────────────────────────────────────
        case.Solver.CanSolve = True

        # ── 収束待機（最大 SOLVER_TIMEOUT_SEC 秒）──────────────────────────────
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

        # ── 結果取得 ───────────────────────────────────────────────────────────
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
        # COMエラー・型エラーなど予期せぬ例外
        result["Status"] = "NonConverged"
        print(f"  [HYSYS ERROR] {exc}", flush=True)

    return result
