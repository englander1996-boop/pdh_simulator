"""
simulation モジュール — exp 用の定型処理 (結果表示など) を集約。

設計変数の値は exp/exp*.py 側に直書きする。本モジュールは「実験ごとに
変えない部分」(表示関数、ヘッダ整形等) を提供して exp 側を実験本体に集中させる。

使用例:
    from simulation import display_full_results, hdr
    ...
    result = evaluate(design, config, verbose=True)
    display_full_results(result, design, config)
"""

from simulation.display import (
    display_full_results, hdr, show_stream,
    show_input_snapshot,
    show_hi_summary, show_stage2_synthesis,
    show_final_summary_box,
)
from simulation.exp_runner import (
    run_with_capture, run_exp, outer_iter_progress,
)

__all__ = [
    'display_full_results',
    'hdr',
    'show_stream',
    'show_input_snapshot',
    'show_hi_summary',
    'show_stage2_synthesis',
    'show_final_summary_box',
    # exp_runner
    'run_with_capture', 'run_exp', 'outer_iter_progress',
]
