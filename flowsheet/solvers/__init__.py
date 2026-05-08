"""
tear stream 加速法のパッケージ。

config/operating.toml の solver.inner.method で選択可能:
  - "successive_substitution": α-relaxation (旧来法、ベースライン)
  - "wegstein": Wegstein 加速 (推奨)
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
