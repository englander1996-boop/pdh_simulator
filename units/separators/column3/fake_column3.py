"""
Dist3: C3 スプリッター (ダミーモデル)

膜分離透過液からポリマーグレードのプロピレン(C3H6)を回収する。
PDH プロセスで最もエネルギー集約的なカラム（C3H6/C3H8 揮発度差が小さい）。
  塔頂: C3H6(B) 製品 (99.5%+ 純度)
  塔底: C3H8(A) → Mixer → Reactor へリサイクル

デフォルト設計条件:
  P_col = 20 bar, N_stages = 200, reflux = 12.0
  B→0.99 (製品回収), A→0.02 (C3H8 は塔底へ)
"""

import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.distillation_core import (
    DistDesignVars, DistFixedParams, DistResult, simulate_distillation_column,
)
from units.utils.process_stream import ProcessStream

_DEFAULT_DESIGN = DistDesignVars(
    P_col        = 20.0e5,
    N_stages     = 200,    # C3H6/C3H8 分離は 150-250 段が典型的
    reflux_ratio = 12.0,   # 高還流比（C3 スプリッターは通常 10-15）
    split_fracs  = {
        'A': 0.02,   # C3H8  → 塔底（リサイクル）
        'B': 0.99,   # C3H6  → 塔頂（製品）
        'C': 1.00,   # H2    → 塔頂
        'D': 0.99,   # C2H4  → 塔頂
        'E': 1.00,   # CH4   → 塔頂
        'F': 0.05,   # C2H6  → 塔底
        'Z': 0.01,   # C4H10 → 塔底
    },
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column3(
    feed:   ProcessStream,
    design: DistDesignVars | None = None,
    fixed:  DistFixedParams | None = None,
) -> DistResult:
    """Dist3 (C3 スプリッター) をシミュレーションする。"""
    return simulate_distillation_column(
        design if design is not None else _DEFAULT_DESIGN,
        feed,
        fixed  if fixed  is not None else _DEFAULT_FIXED,
    )
