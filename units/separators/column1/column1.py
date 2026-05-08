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
    # 設計判断 (2026-05-09): K_method='pr' 切替に伴い再チューニング。
    # PR で α(C3/C4 @17bar) ≈ 2.3 (CC は 3.4 と過大推定) のため R_min が
    # 0.44 → 0.95 に上昇。R = 1.5 (= R_min × 1.57) で reflux スイープ後の
    # TAC 最低点近傍。Dist1 は OPEX 影響が小さいので margin を広めに取る。
    reflux_ratio     = 1.5,
    LK               = 'A',       # C3H8
    HK               = 'Z',       # C4H10
    recovery_LK_top  = 0.99,
    recovery_HK_bot  = 0.99,
    # 設計判断 (2026-05-09): K_method='pr' に復帰。
    # 旧版 (2026-05-08) は「塔平均 T で K=phi_L/phi_V」を計算しており、
    # この T-P-x が単相領域に入って Z 根が 1 本 → K_i ≈ 1 になる病理で
    # CC へ退避していた。distillation_core 側を塔頂/塔底それぞれの泡点
    # フラッシュで K を取るように改修したため (alpha_geom = sqrt(top×bot))、
    # PR を本筋に戻す。Dist1 は alpha 大 (C3/C4 ≈ 3〜4) なので CC でも動くが、
    # 高圧下の PR の方が物性精度が高い。
    K_method         = 'pr',
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
