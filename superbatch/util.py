# -*- coding: utf-8 -*-
"""superbatch.util — ログ・時間整形の小物。"""
import sys
import datetime

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass


def log(msg):
    """タイムスタンプ付きで flush 出力 ([[feedback_detailed_logging]] 準拠)。"""
    print(f"[{datetime.datetime.now():%m-%d %H:%M:%S}] {msg}", flush=True)


def fmt_dur(sec):
    """秒 → H:MM:SS / M:SS。"""
    sec = int(max(sec, 0))
    h, r = divmod(sec, 3600)
    m, s = divmod(r, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"
