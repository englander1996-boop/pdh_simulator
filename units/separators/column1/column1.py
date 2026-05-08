"""
Dist1: 脱ブタン塔 (FUG ベース)

設計判断 (2026-05-08):
  Pump1 後の Fresh LPG (30°C, 17 bar 液) から C4H10 を分離。
  塔頂: C3H8 + 軽質ガス → 反応器系へ (膨張弁経由 0.5 bar)
  塔底: C4H10 + 微量 C3 → 廃棄

key 成分:
  LK = 'A' (C3H8): 塔頂回収率 99%
  HK = 'Z' (C4H10): 塔底回収率 99%
  α (17 bar): C3H8/C4H10 ≈ 3〜4 (大、分離容易)
  → N_min ≈ 7、R_min ≈ 0.3〜0.5、N=20, R=1 で十分余裕

q = 1 (液フィード、Pump1 後の過冷却液)
K_method = 'pr' (PR EOS デフォルト)

注意: 旧 fake_column1 は split_fracs ベースだったが、本実装で FUG に置換。
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
    P_col            = 17.0e5,
    N_stages         = 20,
    N_feed           = 10,        # Kirkbride 推奨値はランタイムで計算、ここはデフォ
    # 設計判断 (2026-05-08): reflux スイープで TAC 最低点を採用。
    # FUG R_min ≈ 0.44 (CC) に対し R = 0.6 (= R_min × 1.36)。
    # Dist1 は α 大 (C3/C4 ~3.4) で R_min 小、攻めても余裕あり。
    reflux_ratio     = 0.6,
    LK               = 'A',       # C3H8
    HK               = 'Z',       # C4H10
    recovery_LK_top  = 0.99,
    recovery_HK_bot  = 0.99,
    # 設計判断 (2026-05-08): K_method='cc' を暫定採用。
    # 'pr' (PR EOS) は 17 bar の C3/C4 系で z_factor が単相 root を返し
    # K_i ≈ 1 になる問題があり別途調査要。CC は fake_column と同じ
    # 物性値ベースで動作確認済み (alpha_LK ≈ 3.4)。
    K_method         = 'cc',
    q                = 1.0,        # 飽和液 (Pump1 後の Fresh)
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column1(
    feed:   ProcessStream,
    design: DistDesignVars | None = None,
    fixed:  DistFixedParams | None = None,
) -> DistResult:
    """Dist1 (脱ブタン塔) をシミュレーションする。"""
    return simulate_distillation_column(
        design if design is not None else _DEFAULT_DESIGN,
        feed,
        fixed  if fixed  is not None else _DEFAULT_FIXED,
    )
