"""
Dist2: 脱エタン塔 (ダミーモデル)

スイング反応器出口の冷却・昇圧ガスから H2・軽質ガスを分離する。
  塔頂: H2(C) + CH4(E) + C2H4(D) → PSA 原料
  塔底: C3H8(A) + C3H6(B) + C2H6(F) → Membrane

デフォルト設計条件:
  P_col = 8.5 bar, N_stages = 20, reflux = 2.0
  A→0.02, B→0.01, C→0.99, D→0.95, E→0.98, F→0.05
"""

import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.distillation_core import (
    DistDesignVars, DistFixedParams, DistResult, simulate_distillation_column,
)
from stream.stream import ProcessStream

_DEFAULT_DESIGN = DistDesignVars(
    P_col        = 8.5e5,
    N_stages     = 20,
    reflux_ratio = 2.0,
    split_fracs  = {
        'A': 0.02,   # C3H8  → 塔底
        'B': 0.01,   # C3H6  → 塔底
        'C': 0.99,   # H2    → 塔頂
        'D': 0.95,   # C2H4  → 塔頂
        'E': 0.98,   # CH4   → 塔頂
        'F': 0.05,   # C2H6  → 塔底
        'Z': 0.01,   # C4H10 → 塔底
    },
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column2(
    feed:   ProcessStream,
    design: DistDesignVars | None = None,
    fixed:  DistFixedParams | None = None,
) -> DistResult:
    """Dist2 (脱エタン塔) をシミュレーションする。"""
    return simulate_distillation_column(
        design if design is not None else _DEFAULT_DESIGN,
        feed,
        fixed  if fixed  is not None else _DEFAULT_FIXED,
    )
