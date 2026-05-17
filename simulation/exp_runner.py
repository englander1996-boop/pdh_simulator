"""
simulation/exp_runner.py — exp1 等の「実験スクリプト」共通ランナー

exp1.py / 将来の他 exp スクリプトで共通の処理:
  - stdout を StringIO にキャプチャ → outputs/<label>_<ts>.txt に保存
  - キャプチャ中も stderr に進捗を流し続ける (黙々と待たされる感を回避)
  - 進捗 ticker (バックグラウンドスレッドで N 秒ごとに残り時間ヒントを stderr へ)
  - 例外/失敗時もキャプチャ済みログは残す

設計判断:
  - exp1.py に直書きしていた logging/ticker/save 処理を関数化、別 exp でも再利用可
  - stdout キャプチャ vs 直出力は save_output フラグで切替
"""

import contextlib
import io
import sys
import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional, TypeVar

T = TypeVar('T')


def run_with_capture(
    run_fn:                 Callable[[], T],
    label:                  str,
    output_dir:             Path | str,
    save_output:            bool = True,
    progress_phrase:        Callable[[str], str] | None = None,
    expected_phases:        int = 6,
    tick_interval_sec:      float = 5.0,
) -> T:
    """run_fn() を実行し、stdout キャプチャ + 進捗 ticker + ファイル保存を一括処理。

    Parameters
    ----------
    run_fn : callable
        引数なしで呼べる関数 (キャプチャしたい処理)
    label : str
        出力ファイル名のプレフィックス (例 'exp1' → outputs/exp1_<ts>.txt)
    output_dir : Path | str
        出力ディレクトリ。
    save_output : bool
        True (デフォルト) で stdout を txt に保存、False ならターミナル直出力。
    progress_phrase : callable | None
        ticker 内で「これまでにキャプチャされた stdout」を見て進捗文字列を返す関数。
        例: lambda text: f"外側 iter {text.count('[外側 iter ')} / ~{expected_phases}"
        None なら経過秒数のみ。
    expected_phases : int
        進捗 % 推定の分母 (例 exp1 の外側 iter 6 回想定)
    tick_interval_sec : float
        ticker の更新間隔 [秒]。

    Returns
    -------
    run_fn() の戻り値。
    """
    if not save_output:
        return run_fn()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d%H%M")
    txt_path = output_dir / f'{label}_{timestamp}.txt'

    print(f"[{label}] 計算中... (出力先 {txt_path})", file=sys.stderr, flush=True)

    stop_flag = threading.Event()
    buf = io.StringIO()

    def _ticker():
        elapsed = 0
        while not stop_flag.wait(tick_interval_sec):
            elapsed += int(tick_interval_sec)
            text = buf.getvalue()
            if progress_phrase is not None:
                phase = progress_phrase(text)
            else:
                phase = "実行中"
            print(f"  ... {elapsed}s 経過、{phase}",
                  file=sys.stderr, flush=True)

    ticker_thread = threading.Thread(target=_ticker, daemon=True)
    ticker_thread.start()

    try:
        with contextlib.redirect_stdout(buf):
            result = run_fn()
    finally:
        stop_flag.set()
        ticker_thread.join(timeout=1.0)

    captured = buf.getvalue()
    txt_path.write_text(captured, encoding='utf-8')
    n_lines = len(captured.splitlines())
    size_kb = txt_path.stat().st_size / 1024
    print(f"[{label}] 完了 → {txt_path}  ({n_lines} 行、{size_kb:.1f} KB)",
          file=sys.stderr, flush=True)
    return result


def outer_iter_progress(expected_iters: int = 6) -> Callable[[str], str]:
    """exp1 用 progress_phrase: 「外側 iter N/~6 (推定 X%)」形式。"""
    def _phrase(text: str) -> str:
        n_outer = text.count('[外側 iter ')
        pct = min(95, max(5, int(n_outer * 100 / expected_iters)))
        return (f"外側 iter {n_outer}/~{expected_iters} (推定 ~{pct}%)"
                if n_outer > 0 else "初期化中")
    return _phrase


def run_exp(
    label:                  str,
    eval_callable:          Callable[[], 'FlowsheetResult'],   # type: ignore[name-defined]
    output_dir:             Path | str = 'outputs',
    save_output:            bool = True,
    expected_outer_iters:   int = 6,
):
    """exp1 / exp 系スクリプトのデフォルト実行ラッパ。

    Parameters
    ----------
    label : str
        出力ファイルプレフィックス (例: 'exp1')
    eval_callable : callable
        引数なしで呼べる「シミュレーション本体」関数。FlowsheetResult を返す想定。
    output_dir, save_output : 共通設定
    expected_outer_iters : int
        進捗 % 推定の分母 (= 通常の外側ループ収束回数)

    Returns
    -------
    FlowsheetResult: eval_callable の戻り値
    """
    return run_with_capture(
        run_fn          = eval_callable,
        label           = label,
        output_dir      = output_dir,
        save_output     = save_output,
        progress_phrase = outer_iter_progress(expected_outer_iters),
        expected_phases = expected_outer_iters,
        tick_interval_sec = 5.0,
    )
