"""
flowsheet.evaluate() を cProfile でプロファイリングし、関数別累積時間を出力する。

設計判断 (2026-05-08): 収束高速化施策 (Wegstein, run_one_pass の分割キャッシュ等) を
実装する前に「真のボトルネックがどこにあるか」を確認するためのツール。
ボトルネックが tear 反復 (= solver) でなく PSA/Mem の内部反復にあるなら、
tear 法をどう変えても効かない。

使い方:
    python tools/profile_flowsheet.py
    python tools/profile_flowsheet.py --top 30 --output profile.txt

設計変数は exp/exp1.py と同じ値を使う (リファクタ後の挙動を測定するため)。
"""

import argparse
import cProfile
import io
import os
import pstats
import sys

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from config.load import load_operating_config
from flowsheet import evaluate, FlowsheetDesignVars
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars


def _exp1_design() -> FlowsheetDesignVars:
    """exp/exp1.py の設計変数と同じ値。リファクタ前後の比較用に固定。"""
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=900.0, z_cat=15.0, t_cyc=15.0, D=5.0),
        psa  =PSADesignVars(D_col=3.0, L_bed=20.0, desorption_target=0.35),
        mem  =MemDesignVars(P_H=9.5e5, P_L=1.0e5, A_mem=100000.0, P_dist=20.0e5),
        dist1=ColumnTunables(P_col=17.0e5, N_stages=20, N_feed=10, reflux_ratio=1.5),
        dist2=ColumnTunables(P_col= 8.5e5, N_stages=20, N_feed=10, reflux_ratio=6.0),
        dist3=ColumnTunables(P_col=20.0e5, N_stages=200, N_feed=100, reflux_ratio=12.0),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--top', type=int, default=30,
                    help='累積時間順に上位 N 関数を表示 (default: 30)')
    ap.add_argument('--sort', default='cumulative',
                    choices=['cumulative', 'tottime', 'ncalls'],
                    help='ソート基準 (default: cumulative)')
    ap.add_argument('--output', default=None,
                    help='詳細レポートをファイルに保存 (省略時は stdout のみ)')
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

    config = load_operating_config()
    design = _exp1_design()

    # ---- プロファイル実行 ----
    print("プロファイリング開始 (verbose=False で実行)...")
    profiler = cProfile.Profile()
    profiler.enable()
    result = evaluate(design, config, verbose=False)
    profiler.disable()
    print("完了。")
    print()

    # ---- サマリ (まず結果が妥当か確認) ----
    print("=" * 64)
    print("  evaluate() の結果サマリ (sanity check)")
    print("=" * 64)
    print(f"  is_feasible      : {result.is_feasible}")
    print(f"  effective_TAC    : {result.effective_TAC:.4f} 億円/年")
    if result.economics:
        print(f"  TAC (実コスト)   : {result.economics.TAC:.4f} 億円/年")
    if result.specs:
        print(f"  C3H6 純度        : {result.specs.c3h6_purity_wtfrac*100:.3f} wt%"
              f" {'✓' if result.specs.c3h6_pass else '✗'}")
        print(f"  H2 純度          : {result.specs.h2_purity_molfrac*100:.3f} mol%"
              f" {'✓' if result.specs.h2_pass else '✗'}")
        print(f"  生産量           : {result.specs.production_kmol_h:.2f} kmol/h"
              f" {'✓' if result.specs.production_pass else '✗'}")
    if not result.is_feasible:
        print(f"  failure_reason   : {result.failure_reason}")
    print()

    # ---- プロファイル統計 ----
    print("=" * 64)
    print(f"  cProfile 上位 {args.top} 関数 (sort={args.sort})")
    print("=" * 64)

    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).sort_stats(args.sort)
    stats.print_stats(args.top)
    report = stream.getvalue()
    print(report)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"\n詳細レポートを {args.output} に保存しました。")


if __name__ == '__main__':
    main()
