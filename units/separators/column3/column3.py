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
    # 設計判断 (2026-05-09): K_method='pr' 切替に伴い再チューニング。
    # PR で α(C3H6/C3H8 @20bar) ≈ 1.07 (CC は 1.10 と過大推定) のため R_min が
    # 7.22 → 10.08 に上昇。R = 12.0 (= R_min × 1.19) は経済最適レンジの下限近傍。
    # 旧版の R=7.7 は CC 基準では feasible だったが PR では infeasible。
    # Dist3 は OPEX 支配的 (Q_reb ~80MW) で R を下げる効果が最も大きく、
    # BO で振るときは下限 11.0 程度 (margin 1.09) を限度に。
    reflux_ratio     = 12.0,
    LK               = 'B',         # C3H6 (製品)
    HK               = 'A',         # C3H8 (リサイクル)
    recovery_LK_top  = 0.99,        # C3H6 製品回収率 99%
    recovery_HK_bot  = 0.99,
    # 設計判断 (2026-05-09): K_method='pr' に復帰。Dist3 は α が極小 (1.05〜1.10)
    # で R_min/N_min が α に超敏感なため、本来 PR が必須。distillation_core 側を
    # bubble-point ベースに改修したことで「単相 root → K=1」病理が解消され、
    # PR が安定動作する。CC との比較では PR の方が α が ~1〜2% 高めに出るため
    # R_min が小さくなり、reflux_ratio の最適化余地が広がる方向。
    K_method         = 'pr',
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
