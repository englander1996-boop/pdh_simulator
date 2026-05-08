"""
Dist1: LPG 脱ブタン塔 (ダミーモデル)

LPG 原料（常温常圧）から C4+ 重質分を除去し、C3 以下を反応器へ送る。
  塔頂: C3H8(A) 主体 → Mixer → Reactor
  塔底: C4H10(Z) → 製品 / 廃棄

デフォルト設計条件:
  P_col = 17 bar, N_stages = 25, reflux = 3.0
  A→0.98, B→0.99, C→1.0, D→1.0, E→1.0, F→0.99, Z→0.02
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
    P_col        = 17.0e5,
    N_stages     = 25,
    reflux_ratio = 3.0,
    split_fracs  = {
        'A': 0.98,   # C3H8  → 塔頂
        'B': 0.99,   # C3H6  → 塔頂
        'C': 1.00,   # H2    → 塔頂
        'D': 1.00,   # C2H4  → 塔頂
        'E': 1.00,   # CH4   → 塔頂
        'F': 0.99,   # C2H6  → 塔頂
        'Z': 0.02,   # C4H10 → 塔底
    },
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column1(
    feed:   ProcessStream,
    design: DistDesignVars | None = None,
    fixed:  DistFixedParams | None = None,
) -> DistResult:
    """Dist1 (LPG 脱ブタン塔) をシミュレーションする。"""
    return simulate_distillation_column(
        design if design is not None else _DEFAULT_DESIGN,
        feed,
        fixed  if fixed  is not None else _DEFAULT_FIXED,
    )
