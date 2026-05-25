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
    ColumnTunables, DistDesignVars, DistFixedParams, DistResult,
    simulate_distillation_column,
)
from stream.stream import ProcessStream


# 既定 tunables: BO/exp で上書きされない場合に使う。
# 設計判断 (2026-05-09):
#   PR で α(C3/C4 @17bar) ≈ 2.3 (CC は 3.4 と過大推定) のため R_min が
#   0.44 → 0.95 に上昇。R = 1.5 (= R_min × 1.57) で reflux スイープ後の
#   TAC 最低点近傍。Dist1 は OPEX 影響が小さいので margin を広めに取る。
_DEFAULT_TUNABLES = ColumnTunables(
    P_col        = 17.0e5,
    N_stages     = 20,
    N_feed       = 10,        # Kirkbride 推奨値はランタイムで計算、ここはデフォ
    reflux_ratio = 1.5,
)
_DEFAULT_FIXED = DistFixedParams()


def simulate_column1(
    feed:     ProcessStream,
    tunables: ColumnTunables | None = None,
    fixed:    DistFixedParams | None = None,
) -> DistResult:
    """Dist1 (脱ブタン塔) をシミュレーションする。

    tunables は P_col/N_stages/N_feed/reflux_ratio のみを保持。
    LK/HK/回収率/K_method/q は本ラッパで固定:
      LK='A' (C3H8), HK='Z' (C4H10), recovery 0.99/0.99, K_method='pr', q=1.0
    """
    t = tunables if tunables is not None else _DEFAULT_TUNABLES

    # HYSYS バックエンド経路 (special.py 用)。tunables.solver_method == 'hysys' のとき
    # provider 側で HSC ファイルを選択し、Column1Adapter で実行 → DistResult 返却。
    # LK/HK/recovery は HYSYS では使われないので以下の組立はスキップする。
    if t.solver_method == 'hysys':
        from units.vle.hysys.provider import solve_column1_via_hysys
        return solve_column1_via_hysys(feed, t, fixed if fixed is not None else _DEFAULT_FIXED)

    # SM (surrogate) 経路。学習済み GPR で HYSYS 解を近似 (2026-05-25)。
    # In_CompFraction2 = t.hysys_spec_value を使うので hysys 系の tunables を流用。
    if t.solver_method == 'sm':
        from src.distillation_sm import solve_column1_via_sm
        return solve_column1_via_sm(feed, t, fixed if fixed is not None else _DEFAULT_FIXED)

    # tunables.recovery_* が None なら 0.99 (既定)、float ならその値を採用 (BO 用)
    rec_LK_top = t.recovery_LK_top if t.recovery_LK_top is not None else 0.99
    rec_HK_bot = t.recovery_HK_bot if t.recovery_HK_bot is not None else 0.99
    design = DistDesignVars(
        P_col           = t.P_col,
        N_stages        = t.N_stages,
        N_feed          = t.N_feed,
        reflux_ratio    = t.reflux_ratio,
        LK              = 'A',       # C3H8
        HK              = 'Z',       # C4H10
        recovery_LK_top = rec_LK_top,
        recovery_HK_bot = rec_HK_bot,
        # 設計判断 (2026-05-09): K_method='pr' に復帰。distillation_core 側を
        # bubble-point ベースに改修したため単相 root 病理は解消済み。
        K_method        = 'pr',
        q               = 1.0,       # 飽和液 (Pump1 後の Fresh)
        solver_method   = t.solver_method,
        D_override      = t.D_override,
    )
    return simulate_distillation_column(
        design, feed,
        fixed if fixed is not None else _DEFAULT_FIXED,
    )
