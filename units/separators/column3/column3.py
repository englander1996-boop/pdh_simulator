"""
Dist3: C3 スプリッター (FUG ベース)

設計判断 (2026-05-08):
  Membrane 透過後 (20 bar, 飽和液) からポリマーグレード C3H6 を回収。
  PDH プロセスで最もエネルギー集約的な塔 (C3H6/C3H8 揮発度差が小さい)。
  塔頂: C3H6 製品 (≥99.5 wt%)
  塔底: C3H8 + 微量 C3H6 → リサイクル

key 成分:
  LK = 'B' (C3H6): 塔頂回収率 99%
  HK = 'A' (C3H8): 塔底回収率 99%
  α (20 bar): C3H6/C3H8 ≈ 1.05〜1.10 (極小)
  → N_min ≈ 80〜100、R_min ≈ 5〜8、N=200, R=12 で実機相当

q = 1 (飽和液、Mem 透過後の冷却液)
K_method = 'pr' (PR EOS、α が極小なので精度重視)
"""

import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.distillation_core import (
    ColumnTunables, DistDesignVars, DistFixedParams, DistResult,
    simulate_distillation_column,
)
from stream.stream import ProcessStream


# 既定 tunables: BO/exp で上書きされない場合に使う。
# 設計判断 (2026-05-09):
#   PR で α(C3H6/C3H8 @20bar) ≈ 1.07 (CC は 1.10 と過大推定) のため R_min が
#   7.22 → 10.08 に上昇。R = 12.0 (= R_min × 1.19) は経済最適レンジの下限近傍。
#   Dist3 は OPEX 支配的 (Q_reb ~80MW) で R を下げる効果が最も大きく、
#   BO で振るときは下限 11.0 程度 (margin 1.09) を限度に。
_DEFAULT_TUNABLES = ColumnTunables(
    P_col        = 20.0e5,
    N_stages     = 200,        # C3H6/C3H8 分離は 150〜250 段が典型
    N_feed       = 100,
    reflux_ratio = 12.0,
)
_DEFAULT_FIXED = DistFixedParams()


# 製品純度 spec (C3H6 ≥ 99.5 wt% — contest spec)。
# C3H6 (MW=42.08) と C3H8 (MW=44.10) は MW 差 < 5% なので mol 分率と wt 分率は
# ほぼ等価とみなし、ここでは mol 分率基準で recovery を逆算する。
_DIST3_PURITY_SPEC_MOL = 0.995

# rigorous との乖離余裕。Gilliland (Eduljee) は constant α 仮定の経験式で、
# rigorous PR EOS よりやや楽観 (= 同 N で実 recovery がやや低い) になる傾向あり。
# 動的 recovery 計算で得た値に 1/SAFETY_FACTOR を掛けて純度 spec の余裕を持たせる。
# 1.2 は経験的 magic constant (=実 spec margin 20%)、要 rigorous 比較で調整。
_DIST3_PURITY_SAFETY_FACTOR = 1.2


def _dist3_dynamic_recovery_HK_bot(
    F_LK_feed: float,
    F_HK_feed: float,
    rec_LK_top: float,
) -> float:
    """製品純度 spec から recovery_HK_bot を動的計算 (Dist3 特有のロジック)。

    設計判断 (2026-05-19):
      旧版は recovery_HK_bot=0.99 hardcode で、Dist3 のフィード C3H8 量が極少
      (~1%) のとき製品純度 100% の overspec を生んでいた。Gilliland feasibility
      check (src/distillation_core.py:1060-1101) がこの高 recovery 達成のための
      N_stages を強制し、BO が N=170+ を選ばざるを得ない構造になっていた。
      本関数はフィード組成から「純度 99.5 wt% spec を満たす最低 recovery_HK_bot」
      を逆算し、Gilliland check の基準を真の spec 制約に揃える。

    物質収支:
      F_LK_top = rec_LK_top × F_LK_feed
      F_HK_top = (1 - rec_HK_bot) × F_HK_feed
      purity   = F_LK_top / (F_LK_top + F_HK_top) ≥ purity_spec
      → rec_HK_bot ≥ 1 - rec_LK_top × (F_LK_feed/F_HK_feed) × ((1-p)/p)
    """
    target = 1.0 - (1.0 - _DIST3_PURITY_SPEC_MOL) / _DIST3_PURITY_SAFETY_FACTOR
    F_LK = max(F_LK_feed, 1e-9)
    F_HK = max(F_HK_feed, 1e-9)
    rec_HK_bot = 1.0 - rec_LK_top * (F_LK / F_HK) * ((1.0 - target) / target)
    # 物理的に有効な範囲に clip:
    #   下限 0.5: フィード HK が極少のとき rec_HK_bot が負になり得る → 0.5 で頭打ち
    #            (recovery < 0.5 は LK/HK が逆になっており物理的に意味がない)
    #   上限 0.999: 旧 hardcode 0.99 より厳しい指定が来た場合の保守的 cap
    return max(0.5, min(0.999, rec_HK_bot))


def simulate_column3(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,
    fixed:    DistFixedParams | None = None,
) -> DistResult:
    """Dist3 (C3 スプリッター) をシミュレーションする。

    LK/HK/K_method/q は本ラッパで固定: LK='B' (C3H6), HK='A' (C3H8), K_method='pr', q=1.0
    recovery_LK_top は None なら 0.99 既定 (BO 等で上書き可)。
    recovery_HK_bot は None なら **製品純度 99.5 wt% から動的計算** (2026-05-19〜)。
    float 明示なら尊重 (BO で振りたいときの後方互換、要注意: 過剰な高値は overspec を再発させる)。
    """
    t = tunables if tunables is not None else _DEFAULT_TUNABLES
    rec_LK_top = t.recovery_LK_top if t.recovery_LK_top is not None else 0.99
    if t.recovery_HK_bot is not None:
        rec_HK_bot = t.recovery_HK_bot
    else:
        rec_HK_bot = _dist3_dynamic_recovery_HK_bot(
            F_LK_feed = feed.F_in.get('B', 0.0),   # C3H6
            F_HK_feed = feed.F_in.get('A', 0.0),   # C3H8
            rec_LK_top = rec_LK_top,
        )
    design = DistDesignVars(
        P_col           = t.P_col,
        N_stages        = t.N_stages,
        N_feed          = t.N_feed,
        reflux_ratio    = t.reflux_ratio,
        LK              = 'B',       # C3H6 (製品)
        HK              = 'A',       # C3H8 (リサイクル)
        recovery_LK_top = rec_LK_top,
        recovery_HK_bot = rec_HK_bot,
        K_method        = 'pr',
        q               = 1.0,
        solver_method   = t.solver_method,
    )
    return simulate_distillation_column(
        design, feed,
        fixed if fixed is not None else _DEFAULT_FIXED,
    )
