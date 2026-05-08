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
    DistDesignVars, DistFixedParams, DistResult, simulate_distillation_column,
)
from stream.stream import ProcessStream


_DEFAULT_DESIGN = DistDesignVars(
    P_col            = 20.0e5,
    N_stages         = 200,        # C3H6/C3H8 分離は 150〜250 段が典型
    N_feed           = 100,
    # 設計判断 (2026-05-08): reflux スイープで TAC 最低点を採用。
    # FUG R_min ≈ 7.22 (CC, LK=B HK=A) に対し R = 7.7 (= R_min × 1.07)。
    # Dist3 は OPEX 支配的 (Q_reb ~80MW) なので R を下げる効果が最も大きい。
    # ぎりぎりまで攻めた値、BO で R を振るときは下限 7.5 程度に留めること。
    reflux_ratio     = 7.7,        # 標準は 10〜15 だが TAC 最適化のため攻めの値
    LK               = 'B',         # C3H6 (製品)
    HK               = 'A',         # C3H8 (リサイクル)
    recovery_LK_top  = 0.99,        # C3H6 製品回収率 99%
    recovery_HK_bot  = 0.99,
    # 設計判断 (2026-05-08): K_method='cc' 暫定採用 (column1 と同じ理由)。
    # Dist3 は α が極小なので本来 PR が望ましいが、現状 PR の単相判定
    # 問題が解決するまで CC で運用 (alpha_LK ≈ 1.07-1.10 が CC で出る)。
    K_method         = 'cc',
    q                = 1.0,
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
