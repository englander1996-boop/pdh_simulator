"""HYSYS COM 経由の VLE/蒸留塔バックエンド。

公開 API:
  HysysSession            : HSC ファイルを開閉する context manager
  HysysPopupMonitor       : モーダル警告の自動押下
  HysysConnectionError    : 接続失敗
  HysysConvergenceError   : 収束失敗・タイムアウト
  HsysRegistry            : 段数 → HSC パス解決
  set_feed_stage          : フィード段書込み (SpecifyFeedLocation 方式)
  run_column_and_wait     : ソルバ実行＋収束待機
"""

from units.vle.hysys.session import (
    HysysSession,
    HysysPopupMonitor,
    HysysConnectionError,
    HysysConvergenceError,
    set_feed_stage,
    run_column_and_wait,
)
from units.vle.hysys.registry import HsysRegistry

__all__ = [
    "HysysSession",
    "HysysPopupMonitor",
    "HysysConnectionError",
    "HysysConvergenceError",
    "set_feed_stage",
    "run_column_and_wait",
    "HsysRegistry",
]
