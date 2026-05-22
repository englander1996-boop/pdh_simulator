r"""Dist3 (column3) を HYSYS バックエンドで単独実行するスモークテスト。

feed は exp/exp2.py の Dist3 feed と同じ (Mem 透過後 49°C 飽和液):
  C3H8=32.3, C3H6=1191.7 kmol/h
  T=322.15K (49°C), P=20 bar
"""
import os
import sys
import time

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

# 循環 import 回避
import flowsheet  # noqa: F401

from stream.stream import ProcessStream
from src.distillation_core import ColumnTunables, DistFixedParams
from units.separators.column3.column3 import simulate_column3
from units.vle.hysys.provider import shutdown_default_provider


def _h(title: str) -> None:
    print(f"\n{'='*60}\n  {title}\n{'='*60}")


def main() -> int:
    feed = ProcessStream(
        F_in={
            'A': 32.3,      # C3H8
            'B': 1191.7,    # C3H6
            'C': 0.0,
            'D': 0.0,
            'E': 0.0,
            'F': 0.0,
            'Z': 0.0,
        },
        T_in=322.15,
        P_in=20.0e5,
    )

    # Draw Rate: 留出液流量 (kgmol/s 単位、HYSYS 内部単位)
    # exp1 best 設計の C3H6 製品 1188 kmol/h ≒ 0.33 kgmol/s を狙う
    draw_rate_kmolh = 1188.0
    tunables = ColumnTunables(
        P_col=1800.0 * 1000.0,
        N_stages=150,               # HSC column3/150.hsc
        N_feed=1,
        reflux_ratio=12.0,
        solver_method='hysys',
        hysys_spec_value=draw_rate_kmolh / 3600.0,   # kgmol/s
        hysys_feed_stage=100,
    )

    _h("入力")
    print(f"  feed total: {sum(feed.F_in.values()):.1f} kmol/h")
    for k, v in feed.F_in.items():
        if v > 0:
            print(f"    {k}: {v:.1f} kmol/h")
    print(f"  T={feed.T_in-273.15:.1f}C, P={feed.P_in/1e5:.1f} bar")
    print(f"  ColumnTunables: P={tunables.P_col/1e5:.1f} bar, N={tunables.N_stages}")
    print(f"                  Draw Rate={draw_rate_kmolh} kmol/h ({tunables.hysys_spec_value:.5f} kgmol/s)")
    print(f"                  FeedStage={tunables.hysys_feed_stage}")

    _h("HYSYS 実行")
    t0 = time.time()
    try:
        result = simulate_column3(feed, tunables, DistFixedParams())
    except Exception as e:
        import traceback
        print(f"NG: {type(e).__name__}: {e}")
        print(traceback.format_exc())
        shutdown_default_provider()
        return 1
    elapsed = time.time() - t0

    _h("結果")
    eq = result.equipment
    print(f"  実行時間          : {elapsed:.1f} 秒")
    print(f"  feasible          : {eq.feasible}")
    print(f"  message           : {eq.message[:200] if eq.message else ''}")
    print(f"  T_top / T_bot     : {eq.T_top-273.15:.2f} / {eq.T_bot-273.15:.2f} C")
    print(f"  Q_cond / Q_reb    : {eq.Q_cond:.0f} / {eq.Q_reb:.0f} kW")
    print(f"  A_cond / A_reb    : {eq.A_cond_m2:.0f} / {eq.A_reb_m2:.0f} m2")
    print(f"  Cool utility      : {eq.cond_utility_name} ({eq.cond_utility_jpy_per_GJ:.1f} 円/GJ)")
    print(f"  Heat utility      : {eq.reb_utility_name} ({eq.reb_utility_jpy_per_GJ:.1f} 円/GJ)")
    print(f"  CAPEX (HE only)   : {eq.CAPEX:.4f} 億円")
    print()
    print(f"  Top   流量合計    : {sum(result.top.F_in.values()):.2f} kmol/h")
    for k, v in result.top.F_in.items():
        if v > 0.01:
            print(f"      {k}: {v:.2f} kmol/h")
    print(f"  Bot   流量合計    : {sum(result.bottom.F_in.values()):.2f} kmol/h")
    for k, v in result.bottom.F_in.items():
        if v > 0.01:
            print(f"      {k}: {v:.2f} kmol/h")

    feed_total = sum(feed.F_in.values())
    out_total  = sum(result.top.F_in.values()) + sum(result.bottom.F_in.values())
    rel_err    = abs(feed_total - out_total) / feed_total * 100 if feed_total > 0 else 0
    print(f"\n  物質収支 (feed vs top+bot): {feed_total:.2f} vs {out_total:.2f}  "
          f"(rel_err={rel_err:.3f} %)")

    if not eq.feasible:
        print("\n[NG] feasible=False")
        shutdown_default_provider()
        return 2
    if eq.Q_cond <= 0 or eq.Q_reb <= 0:
        print("\n[NG] Q が 0 以下")
        shutdown_default_provider()
        return 3
    if rel_err > 1.0:
        print(f"\n[WARN] 物質収支誤差 {rel_err:.2f}% > 1%")

    print("\n[OK] Dist3 HYSYS バックエンド単独実行 成功")
    shutdown_default_provider()
    return 0


if __name__ == "__main__":
    sys.exit(main())
