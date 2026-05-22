r"""Dist1 (column1) を HYSYS バックエンドで単独実行するスモークテスト。

フローシート全体評価とは切り離して、simulate_column1 が HYSYS 経路で
DistResult を正しく返すかだけを確認する。

feed は exp/exp2.py の Dist1 feed と同じ (Pump1 後 Fresh LPG 想定):
  C3H8 = 1763.8, C4H10 = 196.0 kmol/h, 30°C, 17 bar

期待動作:
  - HYSYS V14 起動
  - column1/40.hsc オープン
  - Comp Fraction-2=0.99, FeedStage=20 で計算
  - DistResult が feasible=True で返る
  - 塔頂・塔底の流量・温度・組成、Q_cond, Q_reb が妥当

使い方: .\.venv\Scripts\python.exe tools\_smoke_test_column1_hysys.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

# 循環 import 回避: distillation_core が flowsheet.heat_integration を、
# flowsheet/__init__ が design 経由で distillation_core を読むため、
# flowsheet を先に初期化する必要がある (exp2.py と同じパターン)。
import flowsheet  # noqa: F401

from stream.stream import ProcessStream
from src.distillation_core import ColumnTunables, DistFixedParams
from units.separators.column1.column1 import simulate_column1
from units.vle.hysys.provider import shutdown_default_provider


def _h(title: str) -> None:
    print(f"\n{'='*60}\n  {title}\n{'='*60}")


def main() -> int:
    # ---- feed: exp2 の Dist1 ベース ----
    feed = ProcessStream(
        F_in={
            'A': 1763.8,    # C3H8
            'B': 0.0,       # C3H6
            'C': 0.0,       # H2
            'D': 0.0,       # C2H4
            'E': 0.0,       # CH4
            'F': 0.0,       # C2H6
            'Z': 196.0,     # C4H10
        },
        T_in=303.15,        # 30°C
        P_in=17.0e5,        # 17 bar
    )

    # ---- 設計変数 ----
    tunables = ColumnTunables(
        P_col=1700.0 * 1000.0,    # 1700 kPa = 17 bar
        N_stages=40,              # HSC column1/40.hsc
        N_feed=1,                 # HYSYS 経路では未使用
        reflux_ratio=2.0,         # HYSYS 経路では未使用
        solver_method='hysys',
        hysys_spec_value=0.99,    # Comp Fraction - 2
        hysys_feed_stage=20,
    )

    _h("入力")
    print(f"  feed: C3H8={feed.F_in['A']:.1f}, C4H10={feed.F_in['Z']:.1f} kmol/h")
    print(f"        T={feed.T_in-273.15:.1f}°C, P={feed.P_in/1e5:.1f} bar")
    print(f"  ColumnTunables: P={tunables.P_col/1e5:.1f} bar, N={tunables.N_stages}")
    print(f"                  Comp Fraction-2={tunables.hysys_spec_value}, FeedStage={tunables.hysys_feed_stage}")

    _h("HYSYS 実行")
    t0 = time.time()
    try:
        result = simulate_column1(feed, tunables, DistFixedParams())
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
    print(f"  message           : {eq.message}")
    print(f"  T_top / T_bot     : {eq.T_top-273.15:.2f} / {eq.T_bot-273.15:.2f} °C")
    print(f"  Q_cond / Q_reb    : {eq.Q_cond:.0f} / {eq.Q_reb:.0f} kW")
    print(f"  A_cond / A_reb    : {eq.A_cond_m2:.0f} / {eq.A_reb_m2:.0f} m2")
    print(f"  Cool utility      : {eq.cond_utility_name} ({eq.cond_utility_jpy_per_GJ:.1f} 円/GJ)")
    print(f"  Heat utility      : {eq.reb_utility_name} ({eq.reb_utility_jpy_per_GJ:.1f} 円/GJ)")
    print(f"  CAPEX (HE only)   : {eq.CAPEX:.4f} 億円")
    print()
    print(f"  Top   流量合計    : {sum(result.top.F_in.values()):.2f} kmol/h")
    print(f"  Top   組成 (key):  C3H8={result.top.F_in.get('A',0):.2f}, "
          f"C4H10={result.top.F_in.get('Z',0):.2f} kmol/h")
    print(f"  Bot   流量合計    : {sum(result.bottom.F_in.values()):.2f} kmol/h")
    print(f"  Bot   組成 (key):  C3H8={result.bottom.F_in.get('A',0):.2f}, "
          f"C4H10={result.bottom.F_in.get('Z',0):.2f} kmol/h")
    print()

    # 物質収支検査
    feed_total = sum(feed.F_in.values())
    out_total  = sum(result.top.F_in.values()) + sum(result.bottom.F_in.values())
    rel_err    = abs(feed_total - out_total) / feed_total * 100 if feed_total > 0 else 0
    print(f"  物質収支 (feed vs top+bot): {feed_total:.2f} vs {out_total:.2f}  "
          f"(rel_err={rel_err:.3f} %)")

    # 妥当性 check
    if not eq.feasible:
        print("\n[NG] feasible=False、HYSYS 計算が失敗している")
        shutdown_default_provider()
        return 2
    if eq.Q_cond <= 0 or eq.Q_reb <= 0:
        print("\n[NG] Q が 0 以下、計算結果が異常")
        shutdown_default_provider()
        return 3
    if rel_err > 1.0:
        print(f"\n[WARN] 物質収支誤差 {rel_err:.2f}% > 1%、組成書込みの問題の可能性")

    print("\n[OK] Dist1 HYSYS バックエンド単独実行 成功")
    shutdown_default_provider()
    return 0


if __name__ == "__main__":
    sys.exit(main())
