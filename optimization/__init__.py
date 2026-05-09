"""
optimization モジュール — top-k 候補の re-evaluation 用の重い計算群。

設計判断 (2026-05-09):
  BO loop 内では計算負荷が許される範囲のものを `flowsheet/` に置き、
  top-k 候補だけ走らせる「重い・離散的・組合せ的」な処理を本モジュールに置く。

現時点の含有物:
  - hen_synthesis.py: HEN (Heat Exchanger Network) の Pinch Design Method 合成
                      Stage 2 (= Stage 1 ピンチ targeting に対する詳細設計)

将来の含有物候補:
  - bo_runner.py: Bayesian Optimization 本体 (search space, GP, acquisition)
  - sensitivity.py: top-k 候補に対する設計変数感度解析
  - report.py: top-k 比較レポート生成
"""

from optimization.hen_synthesis import (
    HEMatch,
    HENResult,
    synthesize_hen,
    apply_synthesis_to_economics,
)

__all__ = [
    'HEMatch',
    'HENResult',
    'synthesize_hen',
    'apply_synthesis_to_economics',
]
