"""
LHSサンプリングによるHYSYS蒸留塔シミュレーション条件生成スクリプト。
HYSYSとの通信部分は含まない。条件データフレームの生成のみを行う。
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats.qmc import LatinHypercube, scale

# Windowsコンソールでの日本語文字化け対策
if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")


# ── 変数定義 ──────────────────────────────────────────────────────────────────

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


# ── LHSサンプル生成関数 ───────────────────────────────────────────────────────

def generate_lhs_conditions(n_samples: int = 300, seed: int = 42) -> pd.DataFrame:
    """
    LHSによるHYSYSシミュレーション条件一覧を生成する。

    Parameters
    ----------
    n_samples : int
        生成するサンプル数（デフォルト 300）
    seed : int
        再現性のための乱数シード

    Returns
    -------
    pd.DataFrame
        LHSサンプル列（LHS_*）とHYSYS入力列（HYSYS_*）を含むデータフレーム
    """
    sampler = LatinHypercube(d=len(LHS_VARIABLE_NAMES), seed=seed)
    unit_samples = sampler.random(n=n_samples)  # shape: (n_samples, 6)、値域 [0, 1)

    # [0,1) → 実スケールへ変換
    scaled = scale(unit_samples, l_bounds=LOWER_BOUNDS, u_bounds=UPPER_BOUNDS)

    df = pd.DataFrame(scaled, columns=LHS_VARIABLE_NAMES)

    # ── HYSYS入力列の計算 ─────────────────────────────────────────────────────

    df["HYSYS_Total_Stages"] = df["LHS_Total_Stages"].round().astype(int)

    df["HYSYS_Feed_Stage"] = (
        df["HYSYS_Total_Stages"] * df["LHS_Feed_Stage_Ratio"]
    ).round().astype(int)

    df["HYSYS_Feed_Pressure"]      = df["LHS_Feed_Pressure"]
    df["HYSYS_Reboiler_Pressure"]  = df["LHS_Feed_Pressure"]
    df["HYSYS_Condenser_Pressure"] = df["LHS_Feed_Pressure"] - 20

    df["HYSYS_Reflux_Ratio"] = df["LHS_Reflux_Ratio"]

    df["HYSYS_Draw_Rate"] = df["LHS_Feed_Flow"] * df["LHS_D_F_Ratio"]

    return df


# ── 実行テスト ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    df = generate_lhs_conditions(n_samples=300, seed=42)

    print(f"生成サンプル数: {len(df)} 件")
    print(f"カラム一覧: {list(df.columns)}\n")
    print("先頭5行:")
    print(df.head())
