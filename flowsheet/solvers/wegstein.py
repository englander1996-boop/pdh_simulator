"""
Wegstein 法 — 履歴 2 点から各変数の最適加速係数を自動計算する SS 拡張。

背景:
  逐次置換 (α=0.5) は内側反復が多くかかる。本系は単調収束 (振動なし) で
  収束比 0.96 → λ ≈ 0.92 と推定される。Wegstein はこの λ を 1 反復ごとに
  実測して最適加速係数 q を自動計算するため、人手チューニング不要。

理論:
  各変数 k について、前回 (x_old, f_old) と今回 (x_curr, f_curr) から
    s_k = (f_curr - f_old) / (x_curr - x_old)     ← 局所的な λ 推定
    q_k = s_k / (s_k - 1)                         ← 最適加速係数
    x_next = q_k × x_curr + (1 - q_k) × f_curr
  s が 1 に近い (= 遅い収束) ほど q が大きく負になり、強い over-relaxation。

安全装置:
  - 第 1 反復は履歴がないため α-relaxation にフォールバック
  - q が極端に負だと発散するので q_min ≥ -5 (標準) でクランプ
  - x_curr - x_old が 0 に近いと数値破綻するので閾値で α-relax にフォールバック

参考: J. F. Boston (1979) "An improved algorithm for solving distillation
columns" あるいは Aspen の DEFAULT 加速法もこれ。
"""

from typing import Dict, Optional, Tuple

from flowsheet.solvers.base import TearAccelerator


class Wegstein(TearAccelerator):
    def __init__(
        self,
        q_min:      float = -5.0,
        q_max:      float =  0.0,
        alpha_init: float =  0.5,
        zero_eps:   float =  1e-12,
    ):
        """
        Parameters
        ----------
        q_min : float
            q の下限 (負ほど強い over-relaxation)。標準値 -5.0。
            -10 以下にすると非線形系で発散リスク。
        q_max : float
            q の上限。0.0 のときは少なくとも完全置換 (α=1) 相当を保証。
            正の値は under-relaxation で安定だが遅い。
        alpha_init : float
            第 1 反復 (履歴なし) で使う α-relax の係数。
        zero_eps : float
            分母が 0 に近いとき α-relax にフォールバックする閾値。
        """
        if q_min > q_max:
            raise ValueError(f"q_min ({q_min}) > q_max ({q_max})")
        self.q_min      = q_min
        self.q_max      = q_max
        self.alpha_init = alpha_init
        self.zero_eps   = zero_eps
        self._prev: Optional[Tuple[Dict[str, float], Dict[str, float]]] = None

    def reset(self) -> None:
        self._prev = None

    def step(self, x_current, x_computed):
        if self._prev is None:
            # 第 1 反復: 履歴がないので α-relax で履歴を仕込む
            self._prev = (dict(x_current), dict(x_computed))
            a = self.alpha_init
            return {k: a * x_computed[k] + (1.0 - a) * x_current[k] for k in x_current}

        x_old, f_old = self._prev
        result = {}
        for k in x_current:
            dx = x_current[k] - x_old[k]
            df = x_computed[k] - f_old[k]

            if abs(dx) < self.zero_eps:
                # 動きがない → α-relax にフォールバック
                a = self.alpha_init
                result[k] = a * x_computed[k] + (1.0 - a) * x_current[k]
                continue

            s = df / dx
            if abs(s - 1.0) < self.zero_eps:
                # s ≈ 1 で発散域、q_max にクランプ
                q = self.q_max
            else:
                q = s / (s - 1.0)
                q = max(self.q_min, min(self.q_max, q))

            # x_next = q × x_curr + (1 - q) × x_computed
            #   q < 0: over-relaxation (computed よりさらに進む)
            #   q = 0: 完全置換 (= SS with α=1)
            #   q > 0: under-relaxation
            result[k] = q * x_current[k] + (1.0 - q) * x_computed[k]

        # 履歴更新
        self._prev = (dict(x_current), dict(x_computed))
        return result
