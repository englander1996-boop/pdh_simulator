"""
tear stream 加速法のパッケージ。

config/operating.toml の solver.inner.method で選択可能 (env PDH_RECYCLE_ACCEL で上書き可):
  - "successive_substitution": α-relaxation (旧来法、ベースライン)
  - "wegstein": Wegstein 加速 (成分ごとスカラー、推奨デフォルト)

注記 (2026-05-27): Anderson 加速も試したが、exp4 ベンチで真値 (tight tol) に対し +3.3% ズレた点に
着地 / 0.1% では収束不能だったため**却下・削除**。詳細は memory:project-recycle-accel-benchmark。
"""

from flowsheet.solvers.base import TearAccelerator
from flowsheet.solvers.successive_substitution import SuccessiveSubstitution
from flowsheet.solvers.wegstein import Wegstein


def make_accelerator(method: str, config) -> TearAccelerator:
    """config から TearAccelerator インスタンスを生成。

    Parameters
    ----------
    method : str
        "successive_substitution" or "wegstein"
    config : InnerSolverSpec
        solver.inner セクションの dataclass
    """
    if method == "successive_substitution":
        return SuccessiveSubstitution(alpha=config.relax)
    elif method == "wegstein":
        return Wegstein(
            q_min     =config.wegstein_q_min,
            q_max     =config.wegstein_q_max,
            alpha_init=config.relax,
        )
    else:
        raise ValueError(
            f"Unknown solver.inner.method: {method!r}. "
            f"Expected 'successive_substitution' or 'wegstein'."
        )


__all__ = [
    "TearAccelerator",
    "SuccessiveSubstitution",
    "Wegstein",
    "make_accelerator",
]
