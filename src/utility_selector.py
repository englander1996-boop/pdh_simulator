"""
ターゲット温度から冷媒/熱媒を自動選択し、その単価を返す。

設計判断 (2026-05-08):
  contest.md §2-3 のユーティリティ階層に従って、ターゲット温度を引数に
  最も適した tier を選ぶ。
  「離散選択」モードと「連続補間」モードの 2 つを提供する。

  (1) 離散モード ("discrete"):
      実プラントの段階的な refrigerant system に即した選択。
      ある温度を境にコストが跳ねる (階段関数) → BO の Gaussian Process で扱いにくい。
      診断・最終設計検証に使う。

  (2) 連続モード ("continuous", デフォルト):
      離散 tier の (supply_T+approach, cost) 点列を区分線形補間。
      連続・単調・C0 の関数になり、最適化器に優しい。
      最適化中はこちらを使い、最終設計時に discrete で round するのが標準パターン。

冷却 tier (供給温度の高い順 = 安い順、target ≥ supply + approach で利用可):
    空冷                : 35°C       30 円/GJ
    冷却水              : 30°C       60 円/GJ
    プロピレン冷媒+15°C  :         300 円/GJ
    プロピレン冷媒+5°C   :         500 円/GJ
    プロピレン冷媒-25°C  :        1200 円/GJ
    プロピレン冷媒-40°C  :        2000 円/GJ
    エチレン冷媒-60°C    :        3500 円/GJ
    エチレン冷媒-75°C    :        5500 円/GJ
    エチレン冷媒-100°C   :        9000 円/GJ

加熱 tier (target ≤ supply - approach で利用可):
    LP Steam : 160°C    1800 円/GJ
    MP Steam : 186°C    2200 円/GJ
    HP Steam : 230°C    2800 円/GJ
    燃料燃焼 : 制約なし  1500 円/GJ  (反応器予熱炉用)
"""

from dataclasses import dataclass
from typing import List, Tuple

from src.cost_parameters import (
    AIR_COOLING_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ,
    PROPYLENE_REFRIG_15C_JPY_PER_GJ,
    PROPYLENE_REFRIG_5C_JPY_PER_GJ,
    PROPYLENE_REFRIG_M25C_JPY_PER_GJ,
    PROPYLENE_REFRIG_M40C_JPY_PER_GJ,
    ETHYLENE_REFRIG_M60C_JPY_PER_GJ,
    ETHYLENE_REFRIG_M75C_JPY_PER_GJ,
    ETHYLENE_REFRIG_M100C_JPY_PER_GJ,
    LP_STEAM_JPY_PER_GJ,
    MP_STEAM_JPY_PER_GJ,
    HP_STEAM_JPY_PER_GJ,
    FUEL_JPY_PER_GJ,
    FURNACE_EFFICIENCY,
)


# 設計判断: ターゲット温度と冷媒/熱媒供給温度の最低伝熱差 [K]。
# 経験則として 10K (logmean ベースで現実的に到達できる下限)。
_MIN_APPROACH_TEMP_K = 10.0


@dataclass(frozen=True)
class UtilityTier:
    """ユーティリティ tier の定義。

    name : 人間可読名 (OPEX 内訳表示や OPEX 計算で識別子になる)
    supply_T_K : 供給温度 [K]  (冷媒なら供給温度、熱媒なら飽和温度)
    jpy_per_GJ : 単価 [円/GJ]  (連続モードでは補間値が入る)
    is_heating : True=加熱, False=冷却
    """
    name:        str
    supply_T_K:  float
    jpy_per_GJ:  float
    is_heating:  bool


# ---------------------------------------------------------------------------
# 冷却 tier 一覧 (供給温度の高い順 = コストの安い順)
# 設計判断: 単調順序を保つことが重要。これにより
#   (a) discrete モード: 「最初に approach を満たす tier」= 「最も安い適切な tier」
#   (b) continuous モード: 補間が単調になり、温度が下がるほどコストが上がる
# ---------------------------------------------------------------------------
_COOLING_TIERS: List[UtilityTier] = [
    UtilityTier('空冷',                308.15, AIR_COOLING_JPY_PER_GJ,           False),
    UtilityTier('冷却水',              303.15, COOLING_WATER_JPY_PER_GJ,         False),
    UtilityTier('プロピレン冷媒+15C',  288.15, PROPYLENE_REFRIG_15C_JPY_PER_GJ,  False),
    UtilityTier('プロピレン冷媒+5C',   278.15, PROPYLENE_REFRIG_5C_JPY_PER_GJ,   False),
    UtilityTier('プロピレン冷媒-25C',  248.15, PROPYLENE_REFRIG_M25C_JPY_PER_GJ, False),
    UtilityTier('プロピレン冷媒-40C',  233.15, PROPYLENE_REFRIG_M40C_JPY_PER_GJ, False),
    UtilityTier('エチレン冷媒-60C',    213.15, ETHYLENE_REFRIG_M60C_JPY_PER_GJ,  False),
    UtilityTier('エチレン冷媒-75C',    198.15, ETHYLENE_REFRIG_M75C_JPY_PER_GJ,  False),
    UtilityTier('エチレン冷媒-100C',   173.15, ETHYLENE_REFRIG_M100C_JPY_PER_GJ, False),
]

# ---------------------------------------------------------------------------
# 加熱 tier 一覧 (供給温度の低い順 = コストの安い順 + 燃料は別枠)
# 燃料は供給温度上限なし扱い (直接燃焼炉で 1000°C 級まで届く)。
# 反応器予熱炉のような高温用途では LP/MP/HP Steam を超えて燃料を選ぶ。
# ---------------------------------------------------------------------------
_HEATING_TIERS: List[UtilityTier] = [
    UtilityTier('LP Steam',  433.15, LP_STEAM_JPY_PER_GJ, True),     # 160°C
    UtilityTier('MP Steam',  459.15, MP_STEAM_JPY_PER_GJ, True),     # 186°C
    UtilityTier('HP Steam',  503.15, HP_STEAM_JPY_PER_GJ, True),     # 230°C
    # 燃料単価は η_furnace=0.85 (化工便覧 18·4·3 表 18·11) で実効単価に補正
    UtilityTier('燃料燃焼', 1273.15, FUEL_JPY_PER_GJ / FURNACE_EFFICIENCY, True),     # ~1000°C 上限
]


# ===========================================================================
# 内部ヘルパ: 連続補間
# ===========================================================================

def _interp_linear(x: float, points: List[Tuple[float, float]]) -> Tuple[float, int]:
    """ソート済み (x, y) 列を線形補間。x が範囲外なら端の値で clamp。

    Returns
    -------
    (y_interp, segment_index)
        segment_index は補間に使った区間 (どの 2 tier 間か) を返す。デバッグ用。
    """
    if x <= points[0][0]:
        return points[0][1], 0
    if x >= points[-1][0]:
        return points[-1][1], len(points) - 1

    for i in range(len(points) - 1):
        x1, y1 = points[i]
        x2, y2 = points[i + 1]
        if x1 <= x <= x2:
            frac = (x - x1) / (x2 - x1) if x2 > x1 else 0.0
            return y1 + frac * (y2 - y1), i
    return points[-1][1], len(points) - 1   # フォールバック (到達しない想定)


def _closest_tier(target_T_K: float, tiers: List[UtilityTier]) -> UtilityTier:
    """target_T_K に対応する tier を返す (連続モードでの supply_T 提供用)。

    設計判断 (2026-05-09): 「最近傍 tier」だと target が tier 境界を僅かに下回る
    場合に approach 不足で ΔT_lm 不成立になる。下流で ΔT を計算する用途では、
    target を物理的に到達可能な「最も安い (=最も暖かい) feasible tier」を選ぶべき。

    Cooling: supply_T + approach ≤ target を満たす中で supply_T 最大の tier。
    Heating: supply_T - approach ≥ target を満たす中で supply_T 最小の tier。
    feasible なものが無ければ最も冷たい (cooling) / 熱い (heating) tier を返す。
    """
    if not tiers:
        return None  # type: ignore
    is_heating = tiers[0].is_heating

    if is_heating:
        feasible = [t for t in tiers
                    if t.supply_T_K - _MIN_APPROACH_TEMP_K >= target_T_K]
        if feasible:
            return min(feasible, key=lambda t: t.supply_T_K)
        # 最も熱い tier (= 燃料燃焼) を返す
        return max(tiers, key=lambda t: t.supply_T_K)

    # cooling
    feasible = [t for t in tiers
                if t.supply_T_K + _MIN_APPROACH_TEMP_K <= target_T_K]
    if feasible:
        return max(feasible, key=lambda t: t.supply_T_K)
    # 最も冷たい tier (= エチレン冷媒-100°C) を返す
    return min(tiers, key=lambda t: t.supply_T_K)


# ===========================================================================
# 公開 API
# ===========================================================================

def select_cooling_utility(target_T_K: float, mode: str = "continuous") -> UtilityTier:
    """target に冷却するために必要なユーティリティを選ぶ。

    Parameters
    ----------
    target_T_K : float
        冷却後の目標温度 [K]
    mode : "continuous" | "discrete"
        "continuous" (デフォルト): 区分線形補間で連続コスト
        "discrete": 階段関数で離散選択 (診断用)
    """
    if mode == "discrete":
        # 元来の離散選択ロジック: 供給温度高い順 (=安い順) に試して、
        # 最初に approach を満たす tier を採用
        for tier in _COOLING_TIERS:
            if target_T_K - _MIN_APPROACH_TEMP_K > tier.supply_T_K:
                return tier
        raise ValueError(
            f"target_T_K = {target_T_K:.2f} K ({target_T_K-273.15:.1f}°C) は"
            f" 提供される最低冷媒温度 (-100°C エチレン冷媒) より低く冷却不能"
        )

    if mode == "continuous":
        # 連続補間: (effective_T, cost) を effective_T 昇順に並べて線形補間
        # effective_T = supply_T + approach (= 「target がこの値以下では使えない」境界)
        points = sorted(
            [(t.supply_T_K + _MIN_APPROACH_TEMP_K, t.jpy_per_GJ) for t in _COOLING_TIERS],
            key=lambda p: p[0],
        )
        if target_T_K < points[0][0]:
            # 最低 tier 以下 → エチレン-100°C でも届かない
            raise ValueError(
                f"target_T_K = {target_T_K:.2f} K ({target_T_K-273.15:.1f}°C) は"
                f" 提供される最低冷媒温度 (-100°C エチレン冷媒) より低く冷却不能"
            )
        cost, _ = _interp_linear(target_T_K, points)
        closest = _closest_tier(target_T_K, _COOLING_TIERS)
        return UtilityTier(
            name      =f"{closest.name}≒",   # 「≒」で連続補間値であることを示す
            supply_T_K=closest.supply_T_K,
            jpy_per_GJ=cost,
            is_heating=False,
        )

    raise ValueError(f"Unknown mode: {mode!r}. Expected 'continuous' or 'discrete'.")


def select_heating_utility(target_T_K: float, mode: str = "continuous") -> UtilityTier:
    """target に加熱するために必要なユーティリティを選ぶ。"""
    if mode == "discrete":
        for tier in _HEATING_TIERS:
            if target_T_K + _MIN_APPROACH_TEMP_K < tier.supply_T_K:
                return tier
        raise ValueError(
            f"target_T_K = {target_T_K:.2f} K ({target_T_K-273.15:.1f}°C) は"
            f" 全熱媒の供給温度を超えている。設計を見直すこと"
        )

    if mode == "continuous":
        # 加熱は target が下がるほど安価 (LP < MP < HP < 燃料)
        # effective_T = supply_T - approach (= target がこの値以上では使えない上限)
        points = sorted(
            [(t.supply_T_K - _MIN_APPROACH_TEMP_K, t.jpy_per_GJ) for t in _HEATING_TIERS],
            key=lambda p: p[0],
        )
        if target_T_K > points[-1][0]:
            raise ValueError(
                f"target_T_K = {target_T_K:.2f} K ({target_T_K-273.15:.1f}°C) は"
                f" 全熱媒の供給温度を超えている。設計を見直すこと"
            )
        cost, _ = _interp_linear(target_T_K, points)
        closest = _closest_tier(target_T_K, _HEATING_TIERS)
        return UtilityTier(
            name      =f"{closest.name}≒",
            supply_T_K=closest.supply_T_K,
            jpy_per_GJ=cost,
            is_heating=True,
        )

    raise ValueError(f"Unknown mode: {mode!r}. Expected 'continuous' or 'discrete'.")


def select_utility(
    target_T_K: float,
    inlet_T_K:  float,
    mode:       str = "continuous",
) -> UtilityTier:
    """入口温度と出口温度から冷却 or 加熱を判定し、最適な tier を返す。

    inlet > target なら冷却、 inlet < target なら加熱。
    """
    if target_T_K < inlet_T_K:
        return select_cooling_utility(target_T_K, mode=mode)
    if target_T_K > inlet_T_K:
        return select_heating_utility(target_T_K, mode=mode)
    # 同温 — 名目的に冷却水を選ぶが Q≈0 なので OPEX も≈0
    return _COOLING_TIERS[1]
