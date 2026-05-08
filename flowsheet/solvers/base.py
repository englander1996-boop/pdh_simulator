"""
Tear stream 加速法の抽象基底クラス。

リサイクル収束における 1 反復の更新ルールを定義する戦略 (strategy pattern)。
具体実装:
  - SuccessiveSubstitution: x_next = α × computed + (1-α) × current  (旧来法)
  - Wegstein: 履歴 2 点から最適加速係数 q を自動計算

設計判断 (2026-05-08, profile 結果反映):
  内側反復ロジック (Python) は全体時間の <1% で、ODE 積分が 99% 以上を占める。
  ただし反復回数を減らせば ODE 呼び出し数も比例して減るため、Wegstein は
  「反復回数を 70 → 20 に減らす」効果で実質 3-4 倍速をもたらす。
"""

from abc import ABC, abstractmethod
from typing import Dict


class TearAccelerator(ABC):
    """tear stream 加速法のインターフェイス。

    Solver は 1 反復ごとに step() を呼び、新しい tear 推定値を受け取る。
    state は flat dict (例: {'d3_A': ..., 'd3_B': ..., 'mem_A': ..., 'mem_B': ...,
    'T_d3': ..., 'T_mem': ...}) を想定。Solver 側でパック/アンパックする。
    """

    @abstractmethod
    def step(self, x_current: Dict[str, float],
             x_computed: Dict[str, float]) -> Dict[str, float]:
        """次の反復の tear 推定値を返す。

        Parameters
        ----------
        x_current : dict
            現在の tear 推定値 (前反復で混合済み)
        x_computed : dict
            x_current を入力にして run_one_pass を呼んだ結果の生 tear 出力 (混合前)

        Returns
        -------
        dict
            次の反復で使う tear 推定値 (混合後)
        """
        ...

    def reset(self) -> None:
        """加速法の内部履歴をリセット。新しい外側 iter ごとに呼ぶ。"""
        pass
