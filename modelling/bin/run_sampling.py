"""
C3/C4蒸留塔サロゲートモデル構築用データ収集スクリプト。
LHSサンプリング → HYSYS実行 → CSV保存 を一括実行する。

実行方法:
    python run_sampling.py
"""

import argparse
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
    "LHS_Feed_Flow",
    "LHS_Reflux_Ratio",
    "LHS_D_F_Ratio",
]

LOWER_BOUNDS = [ 90, 1.5, 0.85]
UPPER_BOUNDS = [110, 2.5, 0.95]

DEBUG_FIXED_FLOWS = {
    "LHS_Feed_Flow": 100.0,
    "HYSYS_Reflux_Ratio": 2.0,
    "HYSYS_Draw_Rate": 90.0,
}


def generate_lhs_conditions(
    n_samples: int = 300,
    seed: int = 42,
    base_total_stages: int = 40,
) -> pd.DataFrame:
    """LHSによるHYSYSシミュレーション条件一覧を生成する。"""
    sampler = LatinHypercube(d=len(LHS_VARIABLE_NAMES), seed=seed)
    unit_samples = sampler.random(n=n_samples)
    scaled = scale(unit_samples, l_bounds=LOWER_BOUNDS, u_bounds=UPPER_BOUNDS)
    df = pd.DataFrame(scaled, columns=LHS_VARIABLE_NAMES)

    df["HYSYS_Total_Stages"] = base_total_stages
    df["HYSYS_Reflux_Ratio"] = df["LHS_Reflux_Ratio"]
    df["HYSYS_Draw_Rate"]    = df["LHS_Feed_Flow"] * df["LHS_D_F_Ratio"]

    return df


def generate_fixed_condition(base_total_stages: int) -> pd.DataFrame:
    """既知の収束確認用に、1行だけ固定条件を生成する。"""
    df = pd.DataFrame([DEBUG_FIXED_FLOWS])
    df["HYSYS_Total_Stages"] = base_total_stages
    df["LHS_D_F_Ratio"] = df["HYSYS_Draw_Rate"] / df["LHS_Feed_Flow"]
    return df


# ════════════════════════════════════════════════════════════════════════════════
# ステップ2: HYSYS 1ケース実行
# ════════════════════════════════════════════════════════════════════════════════

SOLVER_TIMEOUT_SEC = 180   # CanSolve=True 後に収束を待つ最大秒数
SOLVER_POLL_INTERVAL = 1.0
SOLVER_STALL_ABORT_SEC = 30  # 変化が止まったら早期終了する秒数
HYSYS_EMPTY_VALUE  = -32767.0  # HYSYSが「未計算」を示す番兵値

_RESULT_KEYS = [
    "Propane_Purity_Top",
    "Butane_Purity_Bottom",
    "Reboiler_Duty",
    "Condenser_Duty",
    "Top_Temperature",
    "Bottom_Temperature",
]


def _is_hysys_empty(v) -> bool:
    """HYSYSの未計算番兵値か NaN かを判定する。"""
    try:
        return np.isnan(float(v)) or abs(float(v) - HYSYS_EMPTY_VALUE) < 1.0
    except Exception:
        return True


def find_comp_indices(case) -> dict:
    """FluidPackage からコンポーネント名→インデックスの辞書を返す。
    HYSYS COM では stream.ComponentMolarFraction(index) で取得するため、
    Propane と n-Butane のインデックスを起動時に一度だけ特定する。
    """
    pkg   = case.Flowsheet.FluidPackage
    n     = pkg.Components.Count
    names = [pkg.Components.Item(i).Name for i in range(n)]
    print(f"    FluidPackage コンポーネント: {names}", flush=True)

    indices = {}
    for i, name in enumerate(names):
        low = name.lower()
        if "propane" in low or low in ("c3", "c3h8"):
            indices.setdefault("Propane", i)
        if "butane" in low and "iso" not in low:
            indices.setdefault("n-Butane", i)

    missing = [k for k in ("Propane", "n-Butane") if k not in indices]
    if missing:
        raise RuntimeError(f"コンポーネントが見つかりません: {missing}。全コンポーネント: {names}")

    print(f"    Propane index = {indices['Propane']}, n-Butane index = {indices['n-Butane']}", flush=True)
    return indices


def _get_component_molar_fraction(stream, component_index: int):
    """HYSYS の成分モル分率をタプルのインデックスで取得する。"""
    try:
        fractions = stream.ComponentMolarFractionValue
        if isinstance(fractions, (tuple, list)):
            return fractions[component_index]
    except Exception as e:
        pass

    try:
        collection = stream.ComponentMolarFractions
        return collection.Item(component_index).Value
    except Exception as e:
        pass

    raise RuntimeError(f"成分モル分率を取得できません: index {component_index}")


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


def _wait_for_column_convergence(case, col, timeout_sec: int = SOLVER_TIMEOUT_SEC) -> tuple[bool, str]:
    """列ソルバーの収束を同一スレッドで待つ。
    返されない場合は Temperature の -32767 → 有効値変化で判定する。"""
    deadline = time.monotonic() + timeout_sec
    stall_deadline = time.monotonic() + SOLVER_STALL_ABORT_SEC
    last_valid_value = None
    found_valid = False

    print(f"    [6]  Temperature 有効値到達を待機中", end="", flush=True)
    while time.monotonic() < deadline:
        try:
            dist = case.Flowsheet.MaterialStreams.Item("Distillate")
            temp_val = dist.Temperature.Value

            is_empty = _is_hysys_empty(temp_val)

            if not is_empty:
                if not found_valid:
                    found_valid = True
                    stall_deadline = time.monotonic() + SOLVER_STALL_ABORT_SEC
                if last_valid_value is not None and abs(float(temp_val) - float(last_valid_value)) < 0.1:
                    return True, f"Top_Temp={temp_val:.2f} (stable)"
                last_valid_value = temp_val
                stall_deadline = time.monotonic() + SOLVER_STALL_ABORT_SEC
            elif found_valid and time.monotonic() >= stall_deadline:
                return False, f"temperature lost for {SOLVER_STALL_ABORT_SEC}s"

        except Exception as exc:
            pass

        time.sleep(SOLVER_POLL_INTERVAL)
        print(".", end="", flush=True)

    if found_valid:
        return True, f"Top_Temp={last_valid_value:.2f} (timeout but valid)"
    return False, f"timeout after {timeout_sec}s (no valid temp)"


def run_hysys_simulation(hysys_app, case, row_data: pd.Series, comp_indices: dict) -> dict:
    """HYSYSにLHSサンプル1行分の条件を入力し、収束後に結果を取得して返す。"""
    result = {key: np.nan for key in _RESULT_KEYS}
    result["Status"] = "NonConverged"

    # ── [1] オブジェクト取得 ─────────────────────────────────────────────
    print(f"    [1] オブジェクト取得中...", flush=True)
    try:
        feed    = case.Flowsheet.MaterialStreams.Item("Feed")
        col     = case.Flowsheet.Operations.Item("T-100")
        specs   = _get_specs_collection(col.ColumnFlowsheet)
        dist    = case.Flowsheet.MaterialStreams.Item("Distillate")
        bottoms = case.Flowsheet.MaterialStreams.Item("Bottoms")
        qr      = case.Flowsheet.EnergyStreams.Item("Qr")
        qc      = case.Flowsheet.EnergyStreams.Item("Qc")
        print(f"    [1] OK", flush=True)
    except Exception as e:
        print(f"    [1] FAILED: {e}", flush=True)
        return result

    # ── [2] ソルバー停止 ──────────────────────────────────────────────────
    print(f"    [2] ソルバー停止...", flush=True)
    try:
        case.Solver.CanSolve = False
        time.sleep(1)
        print(f"    [2] OK", flush=True)
    except Exception as e:
        print(f"    [2] FAILED: {e}", flush=True)
        return result

    # ── [3] 入力設定 ──────────────────────────────────────────────────────
    feed_flow = float(row_data["LHS_Feed_Flow"])
    rr        = float(row_data["HYSYS_Reflux_Ratio"])
    dr        = float(row_data["HYSYS_Draw_Rate"])
    print(f"    [3] 入力値: Feed={feed_flow:.2f} kmol/h, RR={rr:.3f}, DrawRate={dr:.2f} kmol/h", flush=True)
    try:
        feed.MolarFlow.Value = feed_flow
        print(f"    [3] FeedFlow OK", flush=True)
    except Exception as e:
        print(f"    [3] FeedFlow FAILED: {e}", flush=True)
        return result

    def _set_spec(name: str, value: float) -> None:
        spec_obj = specs.Item(name)
        for attr in ("GoalValue", "Value", "TargetValue", "SpecValue"):
            try:
                setattr(spec_obj, attr, value)
                return
            except Exception:
                pass
        raise RuntimeError(f"'{name}' スペック値書き込み失敗")

    print(f"    [4] Reflux Ratio スペック設定...", flush=True)
    try:
        _set_spec("Reflux Ratio", rr)
        print(f"    [4] OK", flush=True)
    except Exception as e:
        print(f"    [4] FAILED: {e}", flush=True)
        return result

    print(f"    [5] Draw Rate スペック設定...", flush=True)
    try:
        _set_spec("Draw Rate", dr)
        print(f"    [5] OK", flush=True)
    except Exception as e:
        print(f"    [5] FAILED: {e}", flush=True)
        return result

    # ── [6] ソルバー起動 → 収束待機 ──────────────────────────────────────
    print(f"    [6] ソルバー起動 → 収束待ち（最大 {SOLVER_TIMEOUT_SEC}s）...", flush=True)
    try:
        case.Solver.CanSolve = True
    except Exception as e:
        print(f"    [6] CanSolve=True FAILED: {e}", flush=True)
        return result

    converged, detail = _wait_for_column_convergence(case, col, timeout_sec=SOLVER_TIMEOUT_SEC)
    if not converged:
        print(f"\n    [6] Timeout ({detail})", flush=True)
        result["Status"] = "Timeout"
        return result

    print(f"\n    [6] ソルバー収束完了 ({detail})", flush=True)

    # ── [7] 結果読み取り（-32767 チェック付き）────────────────────────────
    print(f"    [7] 結果読み取り...", flush=True)
    time.sleep(1)  # 念のため少し待つ

    props_ok = True

    # ComponentMolarFraction(index) ← HYSYS V14 の正式アクセス方法（整数インデックス）
    try:
        val = _get_component_molar_fraction(dist, comp_indices["Propane"])
        if _is_hysys_empty(val):
            raise ValueError(f"未計算値: {val}")
        result["Propane_Purity_Top"] = val
        print(f"    [7a] Propane_Purity_Top = {val:.4f}", flush=True)
    except Exception as e:
        print(f"    [7a] FAILED: {e}", flush=True)
        props_ok = False

    try:
        val = _get_component_molar_fraction(bottoms, comp_indices["n-Butane"])
        if _is_hysys_empty(val):
            raise ValueError(f"未計算値: {val}")
        result["Butane_Purity_Bottom"] = val
        print(f"    [7b] Butane_Purity_Bottom = {val:.4f}", flush=True)
    except Exception as e:
        print(f"    [7b] FAILED: {e}", flush=True)
        props_ok = False

    try:
        result["Reboiler_Duty"] = qr.HeatFlow.Value
        if _is_hysys_empty(result["Reboiler_Duty"]):
            raise ValueError(f"未計算値: {result['Reboiler_Duty']}")
        print(f"    [7c] Reboiler_Duty = {result['Reboiler_Duty']:.2f}", flush=True)
    except Exception as e:
        print(f"    [7c] FAILED: {e}", flush=True)
        props_ok = False

    try:
        result["Condenser_Duty"] = qc.HeatFlow.Value
        if _is_hysys_empty(result["Condenser_Duty"]):
            raise ValueError(f"未計算値: {result['Condenser_Duty']}")
        print(f"    [7d] Condenser_Duty = {result['Condenser_Duty']:.2f}", flush=True)
    except Exception as e:
        print(f"    [7d] FAILED: {e}", flush=True)
        props_ok = False

    try:
        result["Top_Temperature"] = dist.Temperature.Value
        if _is_hysys_empty(result["Top_Temperature"]):
            raise ValueError(f"未計算値: {result['Top_Temperature']}")
        print(f"    [7e] Top_Temperature = {result['Top_Temperature']:.2f}", flush=True)
    except Exception as e:
        print(f"    [7e] FAILED: {e}", flush=True)
        props_ok = False

    try:
        result["Bottom_Temperature"] = bottoms.Temperature.Value
        if _is_hysys_empty(result["Bottom_Temperature"]):
            raise ValueError(f"未計算値: {result['Bottom_Temperature']}")
        print(f"    [7f] Bottom_Temperature = {result['Bottom_Temperature']:.2f}", flush=True)
    except Exception as e:
        print(f"    [7f] FAILED: {e}", flush=True)
        props_ok = False

    result["Status"] = "Success" if props_ok else "PartialError"
    return result


# ════════════════════════════════════════════════════════════════════════════════
# ステップ3: メインループ
# ════════════════════════════════════════════════════════════════════════════════

_HERE          = Path(__file__).parent
HSC_PATH       = str(_HERE / "C3C4_Splitter_Base.hsc")
OUTPUT_CSV     = _HERE / "lhs_results.csv"
TARGET_SUCCESS = 5   # 5件の収束データを取る
N_SAMPLES      = 30  # 30サンプルで試す
SEED           = 42


def main_sampling_loop(
    hsc_path:       str  = HSC_PATH,
    output_csv:     Path = OUTPUT_CSV,
    n_samples:      int  = N_SAMPLES,
    target_success: int  = TARGET_SUCCESS,
    seed:           int  = SEED,
    use_fixed_case:  bool = False,
) -> None:
    """LHSサンプルまたは固定1ケースをHYSYSに投入し、結果をCSVに追記する。"""
    hysys_app   = None
    case        = None
    n_success   = 0
    n_trials    = 0
    first_write = True

    try:
        print("HYSYS接続中...", flush=True)
        hysys_app = win32.Dispatch("HYSYS.Application")
        hysys_app.Visible = True
        case = hysys_app.SimulationCases.Open(hsc_path)
        print(f"HYSYSファイルを開きました: {hsc_path}")

        _col_tmp = case.Flowsheet.Operations.Item("T-100")
        base_total_stages = _col_tmp.ColumnFlowsheet.Operations.Item("Main Tower").NumberOfStages
        print(f"ベースケース段数: {base_total_stages} 段（固定）")

        print("コンポーネントインデックス特定中...")
        comp_indices = find_comp_indices(case)
        print()

        if use_fixed_case:
            df = generate_fixed_condition(base_total_stages=base_total_stages)
            print("固定条件サンプル生成完了: 1 件 (まず接続確認用)\n")
        else:
            df = generate_lhs_conditions(
                n_samples=n_samples, seed=seed, base_total_stages=base_total_stages
            )
            print(f"LHSサンプル生成完了: {len(df)} 件 (目標収束数: {target_success} 件)\n")

        for _, row_data in df.iterrows():
            n_trials += 1
            print(f"\n{'='*60}", flush=True)
            print(f"[試行 {n_trials}/{len(df)}]", flush=True)

            sim_result = run_hysys_simulation(hysys_app, case, row_data, comp_indices)

            row_out = row_data.to_dict()
            row_out.update(sim_result)

            pd.DataFrame([row_out]).to_csv(
                output_csv, mode="w" if first_write else "a", header=first_write, index=False
            )
            first_write = False

            status = sim_result["Status"]
            if status == "Success":
                n_success += 1
                success_goal = 1 if use_fixed_case else target_success
                print(f"→ Success ({n_success}/{success_goal})", flush=True)
            else:
                print(f"→ {status}", flush=True)

            if (not use_fixed_case) and n_success >= target_success:
                print(f"\n目標 {target_success} 件の収束データ取得完了。ループを終了します。")
                break

    finally:
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
    parser = argparse.ArgumentParser(description="HYSYS 向け LHS サンプリング / 固定ケース検証")
    parser.add_argument("--fixed", action="store_true", help="固定の1ケースだけを実行して接続確認する")
    parser.add_argument("--n-samples", type=int, default=N_SAMPLES, help="LHS サンプル数")
    parser.add_argument("--target-success", type=int, default=TARGET_SUCCESS, help="LHS で目標とする収束数")
    parser.add_argument("--seed", type=int, default=SEED, help="乱数シード")
    parser.add_argument("--output-csv", type=Path, default=OUTPUT_CSV, help="出力CSVパス")
    args = parser.parse_args()

    main_sampling_loop(
        output_csv=args.output_csv,
        n_samples=args.n_samples,
        target_success=args.target_success,
        seed=args.seed,
        use_fixed_case=args.fixed,
    )
