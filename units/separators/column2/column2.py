"""
Dist2: 脱エタン塔 (FUG ベース)

設計判断 (2026-05-08):
  反応器出口の冷却・圧縮ガス (8.5 bar, 47°C ガス相) から軽質ガスを分離。
  塔頂: H2 + CH4 + C2H4 → PSA 原料
  塔底: C3H8 + C3H6 + C2H6 → Membrane 系へ

key 成分:
  LK = 'D' (C2H4): 塔頂回収率 95%
  HK = 'A' (C3H8): 塔底回収率 98%
  設計判断: C2H6 は微量 (~0.1 mol%) で key にすると Underwood の R_min が
  過大評価される。一方 C2H4 は塔頂主要成分でフィード組成も多く、
  「軽質 vs C3」境界の代表として LK にすると R_min が現実的な値 (1〜2) になる。
  α (8.5 bar): C2H4/C3H8 ≈ 15-16 (大、分離容易)
  軽質ガス (H2, CH4) は K 値がさらに大 → 自動的に塔頂へ
  C2H6 は α=10.5 で C2H4 より塔底寄り → 仕様的には fake_column2 と同じ挙動

q = 0 (気フィード、Comp2b 出口の高圧ガス)
  → 水素含む軽質ガスが多いため、塔頂部分はガスのまま下降せず処理される
K_method = 'pr' (PR EOS、軽質ガスの K 値計算が大事)
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
#   PR で α(C2H4/C3H8 @8.5bar) は CC より大きく出るが、Dist2 入口組成は
#   フローシート上の運転状態 (single pass か recycle ありか) で z_LK=C2H4 の
#   比率が 0.26〜数% まで変動する。z_LK が小さいフィードでは Underwood の
#   R_min が大きくなり、運転状態をまたいで feasible にするには余裕が必要。
#   R = 6.0 は single-pass (R_min ~4.8) で margin 1.25、recycle 想定で margin 大。
_DEFAULT_TUNABLES = ColumnTunables(
    P_col        = 8.5e5,
    N_stages     = 20,
    N_feed       = 10,
    reflux_ratio = 6.0,
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column2(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,
    fixed:    DistFixedParams | None = None,
) -> DistResult:
    """Dist2 (脱エタン塔) をシミュレーションする。

    LK/HK/回収率/K_method/q は本ラッパで固定:
      LK='D' (C2H4), HK='A' (C3H8), recovery 0.95/0.98, K_method='pr', q=0.0
    x_top に H2/CH4 が混じる場合、PR 泡点フラッシュが cryogenic に張り付く
    ため distillation_core._bubble_T_K が CC へ自動フォールバックする。
    """
    t = tunables if tunables is not None else _DEFAULT_TUNABLES
    design = DistDesignVars(
        P_col           = t.P_col,
        N_stages        = t.N_stages,
        N_feed          = t.N_feed,
        reflux_ratio    = t.reflux_ratio,
        LK              = 'D',       # C2H4
        HK              = 'A',       # C3H8
        recovery_LK_top = 0.95,
        recovery_HK_bot = 0.98,
        K_method        = 'pr',
        q               = 0.0,       # 気フィード (Comp2b 出口)
    )
    return simulate_distillation_column(
        design, feed,
        fixed if fixed is not None else _DEFAULT_FIXED,
    )
