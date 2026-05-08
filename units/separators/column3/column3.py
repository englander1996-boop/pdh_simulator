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


def simulate_column3(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,
    fixed:    DistFixedParams | None = None,
) -> DistResult:
    """Dist3 (C3 スプリッター) をシミュレーションする。

    LK/HK/回収率/K_method/q は本ラッパで固定:
      LK='B' (C3H6), HK='A' (C3H8), recovery 0.99/0.99, K_method='pr', q=1.0
    """
    t = tunables if tunables is not None else _DEFAULT_TUNABLES
    design = DistDesignVars(
        P_col           = t.P_col,
        N_stages        = t.N_stages,
        N_feed          = t.N_feed,
        reflux_ratio    = t.reflux_ratio,
        LK              = 'B',       # C3H6 (製品)
        HK              = 'A',       # C3H8 (リサイクル)
        recovery_LK_top = 0.99,
        recovery_HK_bot = 0.99,
        K_method        = 'pr',
        q               = 1.0,
    )
    return simulate_distillation_column(
        design, feed,
        fixed if fixed is not None else _DEFAULT_FIXED,
    )
