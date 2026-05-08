"""
逐次置換 (Successive Substitution) + アンダーリラックス。

リファクタ前から使われていた最も基本的な収束加速法。
  x_next = α × x_computed + (1 - α) × x_current

α は 0 < α ≤ 1 で、小さいほど安定だが収束が遅い。実測で α=0.5 のとき
リサイクル率の高い系では収束比 ~0.96 (200 反復で 1.5 桁収束) と非常に遅い。

Wegstein 法より遅いが、履歴を持たないため第 1 反復から使え、安定。
比較ベースラインとして残している。
"""

from flowsheet.solvers.base import TearAccelerator


class SuccessiveSubstitution(TearAccelerator):
    def __init__(self, alpha: float = 0.5):
        if not (0.0 < alpha <= 1.0):
            raise ValueError(f"alpha must be in (0, 1], got {alpha}")
        self.alpha = alpha

    def step(self, x_current, x_computed):
        a = self.alpha
        return {k: a * x_computed[k] + (1.0 - a) * x_current[k] for k in x_current}
