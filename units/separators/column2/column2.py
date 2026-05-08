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
    DistDesignVars, DistFixedParams, DistResult, simulate_distillation_column,
)
from stream.stream import ProcessStream


_DEFAULT_DESIGN = DistDesignVars(
    P_col            = 8.5e5,
    N_stages         = 20,
    N_feed           = 10,
    # 設計判断 (2026-05-08): FUG R_min ≈ 5.7 (CC) のため余裕を持って 8.0。
    # 旧 fake_column2 の R=2.0 は split_fracs ベースの簡略計算で、Underwood
    # 法の物理から見ると分離仕様 (D 塔頂 95%, A 塔底 98%) を満たすには
    # R≧6 が必要 (q=0 ガスフィード、軽質ガスが大量)。
    reflux_ratio     = 8.0,
    LK               = 'D',       # C2H4 (塔頂主要軽質成分の代表)
    HK               = 'A',       # C3H8 (塔底主要重質成分の代表)
    recovery_LK_top  = 0.95,
    recovery_HK_bot  = 0.98,
    # 設計判断 (2026-05-08): K_method='cc' 暫定採用 (column1 と同じ理由)
    K_method         = 'cc',
    q                = 0.0,        # 気フィード
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
