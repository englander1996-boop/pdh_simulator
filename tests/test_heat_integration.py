"""
flowsheet/heat_integration.py の検証テスト

題材は授業資料の例題 4.2 / 4.3。
- 例題 4.2: 受熱 a,b と与熱 A,B、ΔTmin=10°C → Q_H_min = 30 kJ/s が解答
- 例題 4.3: 同データを 7 区間に分割、ピンチ位置などの内部チェック

参考文献:
  長谷部 伸治, 外輪 健一郎『プロセスシステム工学 (No.4) — 熱交換器ネット
  ワークの最適合成』京都大学 講義資料、令和7年度 (2025).
  §4.4 例題 4.2 (p.4-7), §4.5 例題 4.3 (p.4-10).
"""

import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from flowsheet.heat_integration import (
    HIStream, UtilityTier, pinch_analysis,
)


def C(t_celsius: float) -> float:
    """°C → K"""
    return t_celsius + 273.15


# ===========================================================================
# 例題 4.2 のストリーム定義
# ===========================================================================
#
# 教科書 p.4-7:
#   受熱 a: 40 → 90 °C,  Cp=5.0 kJ/kg/K, F=0.6 kg/s → F·Cp = 3.0 kW/K, Q=150 kW
#   受熱 b: 80 → 110°C,  Cp=5.0,         F=1.2     → F·Cp = 6.0 kW/K, Q=180 kW
#   与熱 A: 125→ 80 °C,  Cp=4.0,         F=1.0     → F·Cp = 4.0 kW/K, Q=180 kW
#   与熱 B: 100→ 60 °C,  Cp=4.0,         F=1.0     → F·Cp = 4.0 kW/K, Q=160 kW
#
# ΔTmin = 10°C のとき:
#   解答: Q_H_min = 30 kJ/s = 30 kW
# ===========================================================================

def example_4_2_streams():
    """例題 4.2 の 4 ストリーム"""
    return [
        HIStream(name='a', T_in_K=C(40),  T_out_K=C(90),  F_Cp_kW_per_K=3.0),
        HIStream(name='b', T_in_K=C(80),  T_out_K=C(110), F_Cp_kW_per_K=6.0),
        HIStream(name='A', T_in_K=C(125), T_out_K=C(80),  F_Cp_kW_per_K=4.0),
        HIStream(name='B', T_in_K=C(100), T_out_K=C(60),  F_Cp_kW_per_K=4.0),
    ]


# ===========================================================================
# 例題 4.2: Q_H_min = 30 kW
# ===========================================================================

def test_example_4_2_Q_H_min():
    """例題 4.2: ΔTmin=10°C のとき Q_H_min = 30 kW"""
    streams = example_4_2_streams()
    result = pinch_analysis(streams, dT_min_K=10.0)

    assert result.feasible, f"infeasible: {result.message}"
    assert abs(result.Q_H_min_kW - 30.0) < 1e-6, (
        f"Q_H_min={result.Q_H_min_kW:.3f} kW, 期待 30.0 kW"
    )


def test_example_4_2_Q_C_min():
    """例題 4.2: 熱バランスから Q_C_min を逆算

    総受熱量 Σ Q_cold = 150 + 180 = 330 kW
    総与熱量 Σ Q_hot  = 180 + 160 = 340 kW
    熱バランス: Q_H_min + Σ Q_hot = Q_C_min + Σ Q_cold
        → Q_C_min = Q_H_min + 340 - 330 = 30 + 10 = 40 kW
    """
    streams = example_4_2_streams()
    result = pinch_analysis(streams, dT_min_K=10.0)
    assert abs(result.Q_C_min_kW - 40.0) < 1e-6, (
        f"Q_C_min={result.Q_C_min_kW:.3f} kW, 期待 40.0 kW"
    )


def test_example_4_2_pinch_temperature():
    """例題 4.2: ピンチ温度

    教科書例題 4.3 のグランドコンポジットカーブ 図 4.27 から、
    ピンチは温度 80°C (cold) / 90°C (hot) に位置する。
    """
    streams = example_4_2_streams()
    result = pinch_analysis(streams, dT_min_K=10.0)
    # cold 側ピンチが 80°C 付近のはず
    assert abs(result.T_pinch_cold_K - C(80)) < 1.0, (
        f"T_pinch_cold={result.T_pinch_cold_K - 273.15:.1f}°C, 期待 80°C 付近"
    )
    assert abs(result.T_pinch_hot_K - C(90)) < 1.0, (
        f"T_pinch_hot={result.T_pinch_hot_K - 273.15:.1f}°C, 期待 90°C 付近"
    )


def test_heat_balance():
    """熱バランス: Q_H_min + Σ Q_hot = Q_C_min + Σ Q_cold が常に成立"""
    streams = example_4_2_streams()
    result = pinch_analysis(streams, dT_min_K=10.0)
    Q_hot_total  = sum(s.Q_total_kW for s in streams if s.is_hot)
    Q_cold_total = sum(s.Q_total_kW for s in streams if not s.is_hot)
    lhs = result.Q_H_min_kW + Q_hot_total
    rhs = result.Q_C_min_kW + Q_cold_total
    assert abs(lhs - rhs) < 1e-6, (
        f"熱バランス違反: Q_H_min({result.Q_H_min_kW}) + Q_hot({Q_hot_total}) "
        f"= {lhs} ≠ Q_C_min({result.Q_C_min_kW}) + Q_cold({Q_cold_total}) = {rhs}"
    )


# ===========================================================================
# ΔTmin = 0 (理想限界)
# ===========================================================================

def test_dT_min_zero_threshold_problem():
    """ΔTmin → 0 のとき Q_H_min は単純な熱バランスに近づく

    ΣQ_cold = 330, ΣQ_hot = 340 → ΔTmin=0 でもピンチがあれば Q_H_min ≠ 0 だが、
    最低限 Q_H_min ≤ 30 (例題 4.2 の値) になる。
    """
    streams = example_4_2_streams()
    result = pinch_analysis(streams, dT_min_K=0.0)
    assert result.Q_H_min_kW <= 30.0 + 1e-6


# ===========================================================================
# ΔTmin を増やすと Q_H_min が単調増加
# ===========================================================================

def test_dT_min_monotonic():
    """ΔTmin を増やすと Q_H_min は単調増加する（一般的性質）"""
    streams = example_4_2_streams()
    Q_prev = -1.0
    for dT in [0.0, 5.0, 10.0, 15.0, 20.0]:
        result = pinch_analysis(streams, dT_min_K=dT)
        assert result.Q_H_min_kW >= Q_prev - 1e-9, (
            f"ΔTmin={dT}: Q_H_min={result.Q_H_min_kW} が単調増加違反 "
            f"(前値 {Q_prev})"
        )
        Q_prev = result.Q_H_min_kW


# ===========================================================================
# 単一ストリームのみ（hot のみ or cold のみ）
# ===========================================================================

def test_only_hot_streams():
    """与熱のみ: Q_H_min=0, Q_C_min = Σ Q_hot"""
    streams = [
        HIStream(name='A', T_in_K=C(125), T_out_K=C(80), F_Cp_kW_per_K=4.0),
        HIStream(name='B', T_in_K=C(100), T_out_K=C(60), F_Cp_kW_per_K=4.0),
    ]
    result = pinch_analysis(streams, dT_min_K=10.0)
    assert abs(result.Q_H_min_kW) < 1e-6
    Q_hot_total = sum(s.Q_total_kW for s in streams)
    assert abs(result.Q_C_min_kW - Q_hot_total) < 1e-6


def test_only_cold_streams():
    """受熱のみ: Q_C_min=0, Q_H_min = Σ Q_cold"""
    streams = [
        HIStream(name='a', T_in_K=C(40), T_out_K=C(90),  F_Cp_kW_per_K=3.0),
        HIStream(name='b', T_in_K=C(80), T_out_K=C(110), F_Cp_kW_per_K=6.0),
    ]
    result = pinch_analysis(streams, dT_min_K=10.0)
    assert abs(result.Q_C_min_kW) < 1e-6
    Q_cold_total = sum(s.Q_total_kW for s in streams)
    assert abs(result.Q_H_min_kW - Q_cold_total) < 1e-6


# ===========================================================================
# tier 配分
# ===========================================================================

def test_utility_tier_assignment():
    """tier 別配分: 例題 4.2 で蒸気 (140°C) と冷却水 (30°C) を渡すと
    Q_H_min がすべて蒸気に、Q_C_min がすべて冷却水に振られる
    """
    streams = example_4_2_streams()
    heating_tiers = [
        UtilityTier(name='LP Steam', supply_T_K=C(140), jpy_per_GJ=1800, is_heating=True),
    ]
    cooling_tiers = [
        UtilityTier(name='冷却水', supply_T_K=C(30), jpy_per_GJ=60, is_heating=False),
    ]
    result = pinch_analysis(
        streams, dT_min_K=10.0,
        heating_tiers=heating_tiers, cooling_tiers=cooling_tiers,
    )
    assert 'LP Steam' in result.utility_breakdown
    assert abs(result.utility_breakdown['LP Steam'] - 30.0) < 1e-6
    assert '冷却水' in result.utility_breakdown
    assert abs(result.utility_breakdown['冷却水'] - 40.0) < 1e-6


# ===========================================================================
# 潜熱: 蒸気凝縮による加熱
# ===========================================================================

def test_latent_heat_steam():
    """蒸気凝縮で受熱を満たす: 110°C 受熱 30 kW を 140°C 蒸気で加熱

    冷たい流体 (cold): 100→110°C, F·Cp=3.0 → Q=30 kW
    熱い流体 (hot, 蒸気): 140°C で潜熱 30 kW (T_in=T_out=140°C は数値上扱えないので
                          141→139°C の極小温度差 + 潜熱で表現するか、F_Cp=0 + 潜熱のみで表現)
    """
    streams = [
        HIStream(name='cold', T_in_K=C(100), T_out_K=C(110), F_Cp_kW_per_K=3.0),
        HIStream(name='steam', T_in_K=C(141), T_out_K=C(139),
                 F_Cp_kW_per_K=0.0, Q_latent_kW=30.0, T_phase_K=C(140)),
    ]
    result = pinch_analysis(streams, dT_min_K=10.0)
    # ΔTmin=10°C: 蒸気 140°C → cold 130°C 以下にしか使えないが、cold の最高は 110°C なのでOK
    # Q_H_min は cold の不足分 = 30 - 30(蒸気) = 0
    assert abs(result.Q_H_min_kW) < 1e-6
    # 蒸気の余剰は 0（全量使用）
    assert abs(result.Q_C_min_kW) < 1e-6


if __name__ == '__main__':
    # 直接実行用: pytest -v 相当
    import traceback

    tests = [
        ('Q_H_min',                  test_example_4_2_Q_H_min),
        ('Q_C_min',                  test_example_4_2_Q_C_min),
        ('pinch temperature',        test_example_4_2_pinch_temperature),
        ('heat balance',             test_heat_balance),
        ('dT_min=0',                 test_dT_min_zero_threshold_problem),
        ('dT_min monotonic',         test_dT_min_monotonic),
        ('only hot streams',         test_only_hot_streams),
        ('only cold streams',        test_only_cold_streams),
        ('utility tier assignment',  test_utility_tier_assignment),
        ('latent heat steam',        test_latent_heat_steam),
    ]
    n_pass = 0
    n_fail = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  [OK]   {name}")
            n_pass += 1
        except AssertionError as e:
            print(f"  [FAIL] {name}: {e}")
            n_fail += 1
        except Exception as e:
            print(f"  [ERR]  {name}: {e}")
            traceback.print_exc()
            n_fail += 1
    print(f"\n結果: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
