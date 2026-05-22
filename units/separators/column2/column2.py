"""
Dist2: 脱エタン塔 (FUG ベース、真正 deethanizer 構成)

設計判断 (2026-05-09):
  反応器出口の冷却・圧縮ガス (8.5 bar, 47°C ガス相) から軽質成分を分離。
  塔頂: H2 + CH4 + C2H4 + C2H6  → PSA 原料 (含むオフガス → 燃料)
  塔底: C3H8 + C3H6 (clean C3 のみ) → Membrane 系へ

key 成分 (2026-05-09 改訂、真正 deethanizer 構成):
  LK = 'F' (C2H6): 塔頂回収率 99%  ← 「最も重い軽キー」 (= 塔頂に行ってほしい中で最も重い)
  HK = 'A' (C3H8): 塔底回収率 99%
  α (8.5 bar): C2H6/C3H8 ≈ 3.9 (CC ~3.5、PR ~4)

  旧版 (2026-05-08) は LK='D' (C2H4) としていたが、これだと C2H6 が非キー扱いで
  Fenske split が α^N_min × ratio_HK の経験式に依存し、N_min ≈ 2 の Dist2 では
  C2H6 が塔底に 20% 漏洩 (= 85 kmol/h) → Mem 入口に C2 が混入し質量保存が破綻。
  実機の deethanizer は LK = 最も重い軽質成分 (C2H6) を取るのが標準で、
  recovery_LK_top=0.99 で C2H6 が 99% 塔頂に固定される。C2H4・CH4・H2 は α 大で
  自動的に塔頂、C3H6 は HK=A の隣で α≈1.13 のため非キーながら 98% 塔底へ。

非キー成分の自動分配 (Fenske 後):
  C2H4 (D): α/α_LK ≈ 1.6 → 塔頂 99.9% (lighter than LK)
  CH4  (E): α/α_LK ≈ 6   → 塔頂 100%
  H2   (C): α/α_LK ≈ ∞   → 塔頂 100%
  C3H6 (B): α/α_LK ≈ 0.29 → 塔底 98% (heavier than HK=A の隣)
  C4H10(Z): α/α_LK ≈ 0    → 塔底 100% (Dist1 で除去済みで Dist2 入口にはほぼ無)

q = 0 (気フィード、Desuperheater で 50°C まで冷やした飽和蒸気)
K_method = 'pr' (PR EOS、α 計算精度重視)
partial_condenser = True (H2/CH4 凝縮不能のため分流型)
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
# 設計判断 (2026-05-09 改訂、真正 deethanizer 化):
#   PR で α(C2H6/C3H8 @8.5bar) ≈ 3.9。recovery 99/99 で N_min ≈ 7 (Fenske)。
#   R_min は Underwood で feed 組成依存だが、概ね 5 程度 (推算)。
#   margin 1.4× で R = 7.0 を採用。N=20 段で margin 約 3× (N_min=7 に対し)。
_DEFAULT_TUNABLES = ColumnTunables(
    P_col        = 8.5e5,
    N_stages     = 20,
    N_feed       = 10,
    reflux_ratio = 7.0,
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column2(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,
    fixed:    DistFixedParams | None = None,
) -> DistResult:
    """Dist2 (脱エタン塔) をシミュレーションする。

    LK/HK/回収率/K_method/q/partial_condenser は本ラッパで固定:
      LK='F' (C2H6), HK='A' (C3H8), recovery 0.99/0.99,
      K_method='pr', q=0.0, partial_condenser=True

    塔頂は H2/CH4/C2H4/C2H6 主体 → PSA へ。
    塔底は C3H8/C3H6 のみ (clean C3) → Membrane へ。
    質量保存が塔内で自然に閉じるよう、LK = 最も重い軽質成分 (C2H6) を選定。
    """
    t = tunables if tunables is not None else _DEFAULT_TUNABLES

    # HYSYS バックエンド経路 (special.py 用)。
    if t.solver_method == 'hysys':
        from units.vle.hysys.provider import solve_column2_via_hysys
        return solve_column2_via_hysys(feed, t, fixed if fixed is not None else _DEFAULT_FIXED)

    rec_LK_top = t.recovery_LK_top if t.recovery_LK_top is not None else 0.99
    rec_HK_bot = t.recovery_HK_bot if t.recovery_HK_bot is not None else 0.99
    design = DistDesignVars(
        P_col           = t.P_col,
        N_stages        = t.N_stages,
        N_feed          = t.N_feed,
        reflux_ratio    = t.reflux_ratio,
        LK              = 'F',       # C2H6 (最も重い軽質成分)
        HK              = 'A',       # C3H8
        recovery_LK_top = rec_LK_top,  # C2H6 塔頂回収率 (既定 0.99、BO で振れる)
        recovery_HK_bot = rec_HK_bot,
        K_method        = 'pr',
        q               = 0.0,       # 気フィード (Desuper 後 50°C 飽和蒸気)
        # 設計判断 (2026-05-09): Dist2 は partial condenser (分流型)。塔頂に
        # H2/CH4 が大量に来るため total condenser だと T_cond が cryogenic に
        # なる病理がある。実機どおり H2/CH4 はそのまま vapor distillate で抜き、
        # 残りの C2 だけを propylene 冷媒で凝縮して reflux する分流モデル。
        partial_condenser = True,
        solver_method   = t.solver_method,
        D_override      = t.D_override,
    )
    return simulate_distillation_column(
        design, feed,
        fixed if fixed is not None else _DEFAULT_FIXED,
    )
