"""
optimization/penalty_scale.py — Adaptive penalty scaling (2026-05-21)

trial 進行に応じて soft penalty 係数 (trace_bypass / proxy / spec) を動的に強化する仕組み。

戦略:
  - 序盤 (trial 0-30%): scale=0.3  — 探索を広げ、infeasible 領域も評価する余地を残す
  - 中盤 (30-70%):      scale=1.0  — 標準
  - 終盤 (70-100%):     scale=3.0  — feasible 領域に強制収束

実装方針:
  ・本モジュールに module-level 変数 _CURRENT_SCALE を持ち、objective.py から
    set_scale() で更新、各 penalty 計算箇所が get_scale() を multiply する。
  ・Optuna が単一プロセスである前提 (n_jobs=1)。並列化時はスケール共有問題あり
    (process-local global なので、複数プロセスで個別 scale → 影響軽微と判断)。
"""

from typing import Callable


# 現在の penalty scale。trial 開始時に objective.py が更新する。
_CURRENT_SCALE: float = 1.0


def set_scale(scale: float) -> None:
    """current scale を上書き (objective が trial 開始時に呼ぶ)。"""
    global _CURRENT_SCALE
    _CURRENT_SCALE = float(scale)


def get_scale() -> float:
    """current scale を読み出し (penalty 計算箇所が multiply する)。"""
    return _CURRENT_SCALE


# ---------------------------------------------------------------------------
# Default schedule
# ---------------------------------------------------------------------------

def default_schedule(trial_num: int, n_total: int) -> float:
    """既定 schedule: 線形 0.2 → 3.0 (2026-05-23 forensic, 境界探索強化)。

    変遷:
      - 旧旧 (3段 step 0.3/1.0/3.0): 段差で TPE が混乱、feas=0 のときも序盤甘さが無効
      - 旧 (線形 0.5→3.0, 2026-05-22): 序盤も「弱くはない」 → 境界 trial が同等扱いに
      - 現 (線形 0.2→3.0, 2026-05-23): 序盤を更に弱く → spec ぎりぎりの trial も TPE
        には「悪くない TAC」として認識させる。Dist3 recovery 変数解放と組合せ、
        「purity 99.5% ギリギリ」設計を境界 trial として TPE に学ばせる狙い。
        feas 率↓ は許容、best TAC 優先 (ユーザー判断 2026-05-23)。
    """
    if n_total <= 1:
        return 1.0
    pct = trial_num / n_total
    return 0.2 + (3.0 - 0.2) * pct


def linear_schedule(trial_num: int, n_total: int,
                    start: float = 0.3, end: float = 3.0) -> float:
    """線形 schedule: trial 0 で start、trial N で end まで線形補間。"""
    if n_total <= 1:
        return end
    pct = trial_num / n_total
    return start + (end - start) * pct
