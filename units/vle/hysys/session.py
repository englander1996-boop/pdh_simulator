r"""HYSYS COM 接続の低レイヤラッパー。

役割:
  - HysysSession: HYSYS アプリ起動 / HSC オープン / 終了の context manager
  - HysysPopupMonitor: モーダル警告 (HYSYS/Aspen) を別スレッドで自動押下
  - set_feed_stage: SpecifyFeedLocation 方式によるフィード段書込み
  - run_column_and_wait: COM の Solver 起動 → GUI 「実行」ボタンクリックの
                          多段フォールバック付き収束待機

設計判断:
  - 構造パラメータ (段数) は COM で変更不可なので、段数別 HSC を切替える方式。
    段数解決は registry.py に委ねる。
  - フィード段は他入力をすべて済ませた「最後」に SpecifyFeedLocation で
    書き込まないと収束しない (G:\サロゲートモデル\main で確立された手順)。
  - ポップアップは別スレッドで poll しないと、ソルバ実行中の警告で計算が止まる。

元コードからの移植元:
  - HysysPopupMonitor: G:\サロゲートモデル\main\column1\lhs_column1.py
  - set_feed_stage / run_column_and_wait / click_run_button_gui:
      G:\サロゲートモデル\main\column{1,2,3}\run_single_column*.py
"""
from __future__ import annotations

import ctypes
import ctypes.wintypes
import datetime
import os
import threading
import time
from pathlib import Path
from typing import Optional, Tuple

import pythoncom
import win32com.client as win32
import win32con
import win32gui
import win32process


# ---------------------------------------------------------------------------
# 例外
# ---------------------------------------------------------------------------

class HysysError(RuntimeError):
    """HYSYS バックエンド共通の基底例外。"""


class HysysConnectionError(HysysError):
    """HYSYS アプリ起動 / HSC オープンに失敗。"""


class HysysConvergenceError(HysysError):
    """ソルバ非収束 / タイムアウト / 結果取得失敗。"""


# ---------------------------------------------------------------------------
# win32 補助 (process name 取得用)
# ---------------------------------------------------------------------------

_KERNEL32 = ctypes.windll.kernel32
_KERNEL32.OpenProcess.argtypes = [
    ctypes.wintypes.DWORD,
    ctypes.wintypes.BOOL,
    ctypes.wintypes.DWORD,
]
_KERNEL32.OpenProcess.restype = ctypes.wintypes.HANDLE
_KERNEL32.QueryFullProcessImageNameW.argtypes = [
    ctypes.wintypes.HANDLE,
    ctypes.wintypes.DWORD,
    ctypes.wintypes.LPWSTR,
    ctypes.POINTER(ctypes.wintypes.DWORD),
]
_KERNEL32.QueryFullProcessImageNameW.restype = ctypes.wintypes.BOOL
_KERNEL32.CloseHandle.argtypes = [ctypes.wintypes.HANDLE]
_KERNEL32.CloseHandle.restype = ctypes.wintypes.BOOL

_BM_CLICK = getattr(win32con, "BM_CLICK", 0x00F5)
_GW_OWNER = getattr(win32con, "GW_OWNER", 4)


# ---------------------------------------------------------------------------
# swap_case 安定化待ち時間 (2026-05-25 短縮)
# ---------------------------------------------------------------------------
# 設計判断 (2026-05-25): swap_case は HSC を Close→Open し直すたびに HYSYS 内部の
# 初期化を待つため固定 sleep を挟む。過去 (〜2026-05-22) は empty-value (Solver 未走)
# バグ回避で保守的に計 7.0s 積まれていたが、1 パスで塔切替 2-3 回 × 内側 ~10 反復の
# ため総時間の支配項になっていた (exp3 290s のうち ~200s が swap sleep)。
# cold 解の結果は待ち時間に依存しない (= 結果中立) ので、安定動作を保てる範囲で短縮。
# 短すぎた場合は fill_outputs (base.py) が empty-value を検出して明示的に失敗扱いに
# するので、無音の精度劣化ではなく可視の penalty として現れる。
# 必要なら下の値を戻すだけで旧挙動に復帰できる。
_SWAP_SLEEP_AFTER_CLOSE_S     = 0.5   # 旧 1.0
_SWAP_SLEEP_AFTER_OPEN_S      = 2.0   # 旧 3.0
_SWAP_SLEEP_AFTER_ACTIVATE_S  = 1.0   # 旧 2.0
_SWAP_SLEEP_AFTER_ACTIVEDOC_S = 0.5   # 旧 1.0


# ---------------------------------------------------------------------------
# HysysPopupMonitor (移植: lhs_column1.HysysPopupMonitor)
# ---------------------------------------------------------------------------

class HysysPopupMonitor:
    """HYSYS/Aspen の自動 dismiss スレッド。

    別スレッドで定期的に全ウィンドウを列挙し、HYSYS/Aspen のモーダルを
    検出すると優先順位に従ってボタンを PostMessage で押す。
    sample 単位の警告ログは set_current_run / consume_run_events で取得。

    DOE ループ中・BO 評価中は常時 start() しておく。
    """

    BUTTON_PRIORITY = (
        "OK", "Yes", "はい", "Continue", "続行",
        "Ignore", "無視", "Close", "閉じる",
    )
    TITLE_KEYWORDS = ("hysys", "aspen")
    PROCESS_KEYWORDS = ("hysys", "aspenhysys")
    PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
    MAX_PROCESS_PATH = 32768
    RETRY_SECONDS = 1.0

    def __init__(self, poll_interval: float = 0.2):
        self.poll_interval = poll_interval
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()
        self._current_run_id = None
        self._events_by_run: dict = {}
        self._last_attempt_by_hwnd: dict = {}

    # -------- ライフサイクル --------
    def start(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._stop_event.clear()
        self._thread = threading.Thread(
            target=self._watch_loop,
            name="HysysPopupMonitor",
            daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        if self._thread and self._thread.is_alive():
            self._thread.join(timeout=max(1.0, self.poll_interval * 3))

    # -------- run スコープのイベント集約 --------
    def set_current_run(self, run_id) -> None:
        with self._lock:
            self._current_run_id = run_id
            if run_id is not None:
                self._events_by_run.setdefault(run_id, [])

    def consume_run_events(self, run_id) -> Tuple[int, str]:
        with self._lock:
            events = self._events_by_run.pop(run_id, [])
        if not events:
            return 0, ""
        details = []
        for ev in events:
            title = self._clean_text(ev.get("title")) or "(no title)"
            message = self._clean_text(ev.get("message"))
            action = self._clean_text(ev.get("action")) or "WM_CLOSE"
            if len(message) > 300:
                message = message[:300] + "...truncated"
            if message:
                details.append(f'{ev["time"]} title="{title}" text="{message}" action="{action}"')
            else:
                details.append(f'{ev["time"]} title="{title}" action="{action}"')
        details_text = " | ".join(details)
        if len(details_text) > 30000:
            details_text = details_text[:30000] + "...truncated"
        return len(events), details_text

    # -------- 内部: 走査ループ --------
    def _watch_loop(self) -> None:
        while not self._stop_event.is_set():
            try:
                self.scan_once()
            except Exception:
                pass
            self._stop_event.wait(self.poll_interval)

    def scan_once(self) -> None:
        def callback(hwnd, _):
            if self._stop_event.is_set():
                return False
            try:
                if self._should_dismiss(hwnd):
                    self._dismiss(hwnd)
            except Exception:
                pass
            return True
        win32gui.EnumWindows(callback, None)

    def _should_dismiss(self, hwnd) -> bool:
        if not win32gui.IsWindowVisible(hwnd):
            return False
        if not self._looks_like_popup(hwnd):
            return False
        title = self._safe_window_text(hwnd).lower()
        proc = self._process_name_for_window(hwnd).lower()
        return (
            any(k in title for k in self.TITLE_KEYWORDS)
            or any(k in proc for k in self.PROCESS_KEYWORDS)
        )

    def _looks_like_popup(self, hwnd) -> bool:
        if self._safe_class_name(hwnd) == "#32770":
            return True
        owner = win32gui.GetWindow(hwnd, _GW_OWNER)
        return bool(owner and self._find_preferred_button(hwnd)[0])

    def _dismiss(self, hwnd) -> None:
        now = time.monotonic()
        if now - self._last_attempt_by_hwnd.get(hwnd, 0) < self.RETRY_SECONDS:
            return
        self._last_attempt_by_hwnd[hwnd] = now

        title = self._safe_window_text(hwnd)
        message = self._message_text(hwnd)
        button_hwnd, button_text = self._find_preferred_button(hwnd)
        action = button_text or "WM_CLOSE"

        if button_hwnd:
            try:
                win32gui.PostMessage(button_hwnd, _BM_CLICK, 0, 0)
            except Exception:
                win32gui.PostMessage(hwnd, win32con.WM_CLOSE, 0, 0)
                action = "WM_CLOSE"
        else:
            win32gui.PostMessage(hwnd, win32con.WM_CLOSE, 0, 0)

        event = {
            "time": datetime.datetime.now().isoformat(timespec="seconds"),
            "title": title,
            "message": message,
            "action": action,
        }
        with self._lock:
            run_id = self._current_run_id
            if run_id is not None:
                self._events_by_run.setdefault(run_id, []).append(event)

    def _find_preferred_button(self, hwnd):
        buttons = []
        for child, cls, text in self._child_controls(hwnd):
            if cls.lower() == "button" and text.strip():
                buttons.append((child, text.strip(), self._normalize_button_text(text)))
        for preferred in self.BUTTON_PRIORITY:
            np_pref = self._normalize_button_text(preferred)
            for child, text, ntext in buttons:
                if np_pref in ntext:
                    return child, text
        return None, ""

    def _message_text(self, hwnd) -> str:
        texts = []
        for _, cls, text in self._child_controls(hwnd):
            if not text.strip() or cls.lower() == "button":
                continue
            texts.append(text.strip())
        return " / ".join(dict.fromkeys(texts))

    def _child_controls(self, hwnd):
        controls = []

        def cb(child, _):
            try:
                if win32gui.IsWindowVisible(child):
                    controls.append((
                        child,
                        self._safe_class_name(child),
                        self._safe_window_text(child),
                    ))
            except Exception:
                pass
            return True
        try:
            win32gui.EnumChildWindows(hwnd, cb, None)
        except Exception:
            pass
        return controls

    def _process_name_for_window(self, hwnd) -> str:
        try:
            _, pid = win32process.GetWindowThreadProcessId(hwnd)
        except Exception:
            return ""
        return self._process_name(pid)

    def _process_name(self, pid) -> str:
        handle = _KERNEL32.OpenProcess(self.PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not handle:
            return ""
        try:
            size = ctypes.wintypes.DWORD(self.MAX_PROCESS_PATH)
            buf = ctypes.create_unicode_buffer(self.MAX_PROCESS_PATH)
            ok = _KERNEL32.QueryFullProcessImageNameW(handle, 0, buf, ctypes.byref(size))
            if not ok:
                return ""
            return os.path.basename(buf.value)
        finally:
            _KERNEL32.CloseHandle(handle)

    def _safe_window_text(self, hwnd) -> str:
        try:
            return win32gui.GetWindowText(hwnd)
        except Exception:
            return ""

    def _safe_class_name(self, hwnd) -> str:
        try:
            return win32gui.GetClassName(hwnd)
        except Exception:
            return ""

    @staticmethod
    def _normalize_button_text(text) -> str:
        return str(text).replace("&", "").replace("...", "").strip().lower()

    @staticmethod
    def _clean_text(text) -> str:
        return " ".join(str(text or "").split())


# ---------------------------------------------------------------------------
# ヘルパー関数: フィード段書込み・ソルバ実行
# ---------------------------------------------------------------------------

def set_feed_stage(col, feed_stream, feed_stage: int) -> bool:
    r"""フィード段を feed_stage に変更する。

    手順 (G:\サロゲートモデル\main\column1\lhs_column1.py の実績手順):
      1. col.ColumnFlowsheet.Operations.Item(0) で Main Tower 取得
      2. tower.SpecifyFeedLocation(int_feed, stage) を試行 (stage は 0始まり/1始まり両方試行)
      3. ダメなら DeleteFeedStream + AddFeedStream フォールバック

    Returns
    -------
    bool : 書込み成功なら True
    """
    col_fs = col.ColumnFlowsheet
    feed_set = False
    try:
        tower = col_fs.Operations.Item(0)
        try:
            int_feed = col_fs.MaterialStreams.Item(feed_stream.Name)
        except Exception:
            int_feed = col_fs.FeedStreams.Item(feed_stream.Name)

        for stage_arg in (feed_stage - 1, feed_stage):
            if feed_set:
                break
            try:
                tower.SpecifyFeedLocation(int_feed, stage_arg)
                feed_set = True
            except Exception:
                pass

        if not feed_set:
            for stage_arg in (feed_stage - 1, feed_stage):
                if feed_set:
                    break
                try:
                    tower.DeleteFeedStream(int_feed)
                    tower.AddFeedStream(int_feed, stage_arg)
                    feed_set = True
                except Exception:
                    pass
    except Exception:
        pass
    return feed_set


def click_run_button_gui(column_name: str) -> None:
    """win32gui で HYSYS カラムダイアログ内の「実行」ボタンを PostMessage する。

    COM の col.Run() が無いビルド向けのフォールバック。
    """
    RUN_LABELS = ("実行", "実行(U)", "Run", "&Run")

    def _find_in_window(hwnd, result):
        if not win32gui.IsWindowVisible(hwnd):
            return True
        try:
            title = win32gui.GetWindowText(hwnd)
        except Exception:
            return True
        if column_name not in title:
            return True

        def _find_button(child, _):
            try:
                cls = win32gui.GetClassName(child)
                text = win32gui.GetWindowText(child)
                if cls.lower() == "button" and any(lbl in text for lbl in RUN_LABELS):
                    result.append(child)
            except Exception:
                pass
            return True

        try:
            win32gui.EnumChildWindows(hwnd, _find_button, None)
        except Exception:
            pass
        return True

    buttons = []
    try:
        win32gui.EnumWindows(_find_in_window, buttons)
    except Exception:
        pass
    if buttons:
        try:
            win32gui.PostMessage(buttons[0], _BM_CLICK, 0, 0)
        except Exception:
            pass


def run_column_and_wait(col, case, column_name: str, timeout: float = 120.0,
                        startup_grace: float = 10.0) -> bool:
    """フィード段変更後にソルバ起動し、収束まで待つ。

    COM API → GUI ボタンクリックの順でフォールバック。

    設計判断 (2026-05-22 改訂):
      旧版は `case.Solver.IsSolving` (= メイン flowsheet の Solver) を polling していたが、
      これは塔の sub-flowsheet とは別物で、塔だけ動かしている時は常に False のままだった。
      → column2 で「収束済み」 GUI 表示なのに私のコードは「Solver 未起動」と誤判定し
        empty value で失敗していた。
      本版は **塔の sub-flowsheet Solver** (`col.ColumnFlowsheet.Solver`) を polling 対象に
      切替え、IsSolving + Converged の両方で判定する。

    Returns
    -------
    bool : timeout 内に塔が収束したら True
    """
    import os as _os
    _timing = _os.environ.get('PDH_HYSYS_TIMING', '0') == '1'
    _t_start_all = time.monotonic()
    _observed_solving = False

    ran = False
    for runner in (
        lambda: col.Run(),
        lambda: col.ColumnFlowsheet.Solver.Run(),
        lambda: col.ColumnFlowsheet.Run(),
    ):
        try:
            runner()
            ran = True
            break
        except Exception:
            pass

    if not ran:
        click_run_button_gui(column_name)

    def _timing_report(result_label: str) -> None:
        if _timing:
            import sys as _sys
            _sys.stderr.write(
                f"[TIMING]   run_column_and_wait({column_name}): "
                f"{time.monotonic() - _t_start_all:5.2f}s  "
                f"IsSolving観測={_observed_solving}  結果={result_label}\n"
            )
            _sys.stderr.flush()

    # 塔の sub-flowsheet Solver を取得
    try:
        col_solver = col.ColumnFlowsheet.Solver
    except Exception:
        col_solver = None

    if col_solver is None:
        # fallback: case.Solver で見る (旧挙動、最低限のサポート)
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                if not case.Solver.IsSolving:
                    _timing_report("fallback-done")
                    return True
            except Exception:
                _timing_report("fallback-exc")
                return True
            time.sleep(0.5)
        _timing_report("fallback-timeout")
        return False

    # Phase 1: IsSolving=True を観測 (startup_grace 内)。観測できなくても先へ進む
    # (HYSYS は CanSolve=True にした瞬間に計算が完了するケースもあるため fail-fast にしない)。
    start_at = time.monotonic()
    grace_deadline = start_at + startup_grace
    while time.monotonic() < grace_deadline:
        try:
            if col_solver.IsSolving:
                _observed_solving = True
                break
        except Exception:
            break
        time.sleep(0.1)

    # Phase 2: IsSolving=False に戻るのを待つ + Converged 確認
    deadline = start_at + timeout
    while time.monotonic() < deadline:
        try:
            solving = col_solver.IsSolving
        except Exception:
            solving = False
        if not solving:
            # Converged プロパティで最終判定
            try:
                _timing_report("converged-check")
                return bool(col_solver.Converged)
            except Exception:
                # Converged 取れない場合は IsSolving=False をもって完了とみなす
                _timing_report("issolving-false")
                return True
        time.sleep(0.5)
    _timing_report("timeout")
    return False


# ---------------------------------------------------------------------------
# HysysSession (context manager)
# ---------------------------------------------------------------------------

class HysysSession:
    """HYSYS アプリ + 1 つの HSC ケースを扱う context manager。

    Usage
    -----
        with HysysSession(hsc_path, visible=False) as sess:
            case = sess.case
            fs   = case.Flowsheet
            col  = fs.Operations.Item("T-100")
            ...

    終了時 (__exit__) に case.Close(False) → app.Quit() を行う。例外発生時も
    確実にクローズするため try/finally を内蔵。PopupMonitor も自動 start/stop。
    """

    def __init__(
        self,
        hsc_path: os.PathLike | str,
        visible: bool = False,
        popup_monitor: bool = True,
        keep_open: bool = False,
    ):
        """
        Parameters
        ----------
        keep_open : bool
            True のとき __exit__ でも HYSYS / HSC を閉じない。デバッグ用 (ユーザーが
            GUI で状態を確認したい場合)。プロセスが残るので使用後に手動で HYSYS を
            閉じる必要あり。
        """
        self.hsc_path = Path(hsc_path).resolve()
        self.visible = visible
        self._use_popup_monitor = popup_monitor
        self._keep_open = keep_open
        self._app = None
        self._case = None
        self.popup_monitor: Optional[HysysPopupMonitor] = None
        self._co_initialized = False

    # -------- context manager --------
    def __enter__(self) -> "HysysSession":
        # スレッド単位で必須。並列化時に重要。
        try:
            pythoncom.CoInitialize()
            self._co_initialized = True
        except Exception:
            self._co_initialized = False

        try:
            self._app = win32.Dispatch("HYSYS.Application")
            self._app.Visible = bool(self.visible)
        except Exception as e:
            self._safe_cleanup()
            raise HysysConnectionError(f"HYSYS 起動失敗: {e}") from e

        if not self.hsc_path.exists():
            self._safe_cleanup()
            raise HysysConnectionError(f"HSC が見つからない: {self.hsc_path}")

        try:
            self._case = self._app.SimulationCases.Open(str(self.hsc_path))
        except Exception as e:
            self._safe_cleanup()
            raise HysysConnectionError(
                f"HSC オープン失敗: {self.hsc_path}: {e}"
            ) from e

        # HSC オープン直後の安定化待ち (run_single_column3.py より、time.sleep(2.0))。
        # GUI 描画完了とコンポーネント初期化を保証するため。
        time.sleep(2.0)

        # 設計判断 (2026-05-22): Visible=True でも SimulationCases.Open() で開いた
        # case は HYSYS V14 GUI に自動表示されない仕様。case.Activate() で明示的に
        # アクティブ化する (失敗してもエラーにせず、計算には影響なし)。
        if self.visible:
            for method_name in ("Activate", "Show"):
                try:
                    method = getattr(self._case, method_name)
                    method()
                    break
                except Exception:
                    pass

        if self._use_popup_monitor:
            self.popup_monitor = HysysPopupMonitor()
            self.popup_monitor.start()

        return self

    def __exit__(self, exc_type, exc, tb):
        if self._keep_open:
            # GUI 確認用: ポップアップ監視だけ止めて、case / app は残す
            if self.popup_monitor is not None:
                try:
                    self.popup_monitor.stop()
                except Exception:
                    pass
                self.popup_monitor = None
            return False
        self._safe_cleanup()
        return False  # 例外を握り潰さない

    def _safe_cleanup(self) -> None:
        if self.popup_monitor is not None:
            try:
                self.popup_monitor.stop()
            except Exception:
                pass
            self.popup_monitor = None

        if self._case is not None:
            try:
                self._case.Solver.CanSolve = False
            except Exception:
                pass
            try:
                self._case.Close(False)
            except Exception:
                pass
            self._case = None

        if self._app is not None:
            try:
                self._app.Quit()
            except Exception:
                pass
            self._app = None

        if self._co_initialized:
            try:
                pythoncom.CoUninitialize()
            except Exception:
                pass
            self._co_initialized = False

    # -------- ショートカット --------
    @property
    def app(self):
        if self._app is None:
            raise HysysConnectionError("未接続 (with ブロック外で参照)")
        return self._app

    @property
    def case(self):
        if self._case is None:
            raise HysysConnectionError("未接続 (with ブロック外で参照)")
        return self._case

    @property
    def flowsheet(self):
        return self.case.Flowsheet

    def solver_stop(self) -> None:
        self.case.Solver.CanSolve = False

    def solver_start(self) -> None:
        self.case.Solver.CanSolve = True

    def components(self):
        return self.case.BasisManager.FluidPackages.Item(0).Components

    # ---- アプリ使い回しで case だけ差替え (連続起動時の COM 安定化用) ----
    def swap_case(self, new_hsc_path) -> None:
        """app を保持したまま case を Close → 新規 Open する。

        フローシート評価で Dist1 → Dist2 → Dist3 と HSC を切替える時、app ごと
        再起動すると COM が不安定になるので、case のみ差替え。

        設計判断 (2026-05-22 改訂): swap 後の sleep を 5秒に延長。HYSYS は HSC を
        Open した直後にバックグラウンドで内部初期化を行うため、即座に書込み・Solver 起動
        すると empty value (Solver 未走) が返るケースが column3 → swap で再現した。
        """
        new_path = Path(new_hsc_path).resolve()
        if not new_path.exists():
            raise HysysConnectionError(f"HSC が見つからない: {new_path}")
        # 既存 case を閉じる
        if self._case is not None:
            try:
                self._case.Solver.CanSolve = False
            except Exception:
                pass
            try:
                self._case.Close(False)
            except Exception:
                pass
            self._case = None
        # close 後の安定化
        time.sleep(_SWAP_SLEEP_AFTER_CLOSE_S)
        # 新規 case を開く (リトライ 3 回)
        last_exc = None
        for attempt in range(3):
            try:
                self._case = self._app.SimulationCases.Open(str(new_path))
                last_exc = None
                break
            except Exception as e:
                last_exc = e
                time.sleep(2.0)   # リトライ前に少し待つ
        if last_exc is not None:
            raise HysysConnectionError(
                f"swap_case で HSC オープン失敗 (3回試行): {new_path}: {last_exc}"
            ) from last_exc
        self.hsc_path = new_path
        # 安定化待ち
        time.sleep(_SWAP_SLEEP_AFTER_OPEN_S)
        # 設計判断 (2026-05-22): swap_case 後は **常に** Activate を呼ぶ。
        # visible=False でも Activate しないと ActiveDocument が前の case のまま残って
        # Solver 起動が前のセッションの状態を引きずる (Dist3 で empty value 発生の原因)。
        for method_name in ("Activate", "Show"):
            try:
                method = getattr(self._case, method_name)
                method()
                break
            except Exception:
                pass
        # Activate 後の追加安定化待ち (HSC 大きいほど必要)
        time.sleep(_SWAP_SLEEP_AFTER_ACTIVATE_S)
        # 念のため ActiveDocument を明示設定 (HYSYS API では set 可能なケースあり)
        try:
            self._app.ActiveDocument = self._case
        except Exception:
            pass
        time.sleep(_SWAP_SLEEP_AFTER_ACTIVEDOC_S)
