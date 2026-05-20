"""BO ループ用の表示 callback。

study.optimize の callbacks= に渡して、各 trial 完了時に compact な状態表示を行う。

設計判断 (2026-05-20、ユーザー要望):
  Optuna デフォルトの logger は "Trial N finished with value: X and parameters: {…}"
  を全 param 込みで 1 行に出すため可読性が低い。
  独自 callback で 1 trial = 4-5 行の構造化された表示に置き換える。
    - 1 行目: status (★ BEST / ✓ feas / ✗ infeas) + TAC + 経過秒
    - 2-4 行目: 全 design vars をユニット別にグループ化
    - 5 行目: progress (完了/全体, feasibility 率, elapsed, ETA, pace)

  ETA は直近 N=20 trial の duration 中央値で計算 (累積平均は trial 間ばらつきで暴れる)。
"""

from __future__ import annotations

import time
from collections import deque
from typing import Callable

import optuna


_BAR_W = 30  # progress bar の幅 (文字)


def _fmt_dur(sec: float) -> str:
    """秒を HH:MM:SS / MM:SS / SSs 形式に。"""
    if sec is None or sec < 0:
        return '--:--'
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    if m > 0:
        return f"{m:d}:{s:02d}"
    return f"{s:d}s"


def _fmt_failure(reason: str) -> str:
    """failure_reason から短い理由ラベルを抽出。"""
    if not reason:
        return ""
    r = str(reason)
    if '生産量' in r:
        return '生産量未達'
    if 'C3H6 純度' in r or '純度' in r:
        return '純度不足'
    if 'H2 純度' in r:
        return 'H2 純度不足'
    if 'プロキシ' in r:
        return 'proxy 罰則'
    if 'trace bypass' in r:
        return 'trace bypass'
    if 'strict recovery' in r:
        return 'strict rec NG'
    if 'penalty' in r.lower() or '罰則' in r or 'CAPEX' in r:
        return 'solver/penalty'
    # 先頭 30 文字に省略
    short = r.split('|')[0].strip()
    return short[:30] + ('…' if len(short) > 30 else '')


def _fmt_vars(params: dict) -> list[str]:
    """全 design vars をユニット別 3 行にまとめる。"""
    g = params.get

    rx = (
        f"Reactor: T={g('T_in_K', 0):.0f}K z={g('z_cat_m', 0):.1f}m "
        f"t={g('t_cyc_min', 0):.1f}min D={g('D_reactor_m', 0):.2f}m"
    )
    psa = (
        f"PSA: D={g('D_psa_col_m', 0):.2f}m L={g('L_psa_bed_m', 0):.1f}m "
        f"des={g('desorption_target', 0):.3f}"
    )
    mem = (
        f"Mem: P_H={g('P_H_Pa', 0)/1e5:.2f}bar A={g('A_mem_m2', 0):.2e}m²"
    )
    line1 = f"       {rx}  |  {psa}  |  {mem}"

    d1 = f"Dist1: P={g('P_dist1_Pa', 0)/1e5:.1f}bar N={g('N_dist1', 0)} R={g('reflux_dist1', 0):.2f}"
    d2 = f"Dist2: P={g('P_dist2_Pa', 0)/1e5:.2f}bar N={g('N_dist2', 0)} R={g('reflux_dist2', 0):.2f}"
    d3 = f"Dist3: P={g('P_dist3_Pa', 0)/1e5:.2f}bar N={g('N_dist3', 0)} R={g('reflux_dist3', 0):.2f}"
    line2 = f"       {d1}   {d2}   {d3}"

    rec_lk = g('rec_LK_top_dist2', None)
    rec_hk = g('rec_HK_bot_dist2', None)
    rec_lk_s = f"{rec_lk:.4f}" if rec_lk is not None else '—'
    rec_hk_s = f"{rec_hk:.4f}" if rec_hk is not None else '—'
    line3 = (
        f"       F_fresh={g('F_C3H8_fresh_kmol_h', 0):.1f}  "
        f"rec_LK_d2={rec_lk_s}  rec_HK_d2={rec_hk_s}"
    )

    return [line1, line2, line3]


def make_compact_callback(n_trials_total: int) -> Callable[[optuna.Study, optuna.trial.FrozenTrial], None]:
    """compact 表示用の callback を生成。

    Parameters
    ----------
    n_trials_total : int
        全体の試行数 (進捗バーと ETA 計算に使用)。

    Returns
    -------
    callback : callable
        study.optimize(callbacks=[callback]) に渡せる関数。

    使用例:
      from optimization.callbacks import make_compact_callback
      cb = make_compact_callback(n_trials_total=cfg.n_trials)
      study.optimize(obj, n_trials=cfg.n_trials, callbacks=[cb])

    state は closure で保持する (グローバル汚染しない):
      - start_time : study 開始時刻 (time.monotonic)
      - prev_best  : 直前のベスト値 (改善幅算出)
      - n_feas     : これまでの feasible trial 数
      - recent_dur : 直近 N trial の所要時間 deque (ETA 中央値計算)
    """
    state = {
        'start_time': None,
        'prev_best':  float('inf'),
        'n_feas':     0,
        'n_done':     0,
        'recent_dur': deque(maxlen=20),
    }

    def _callback(study: optuna.Study, trial: optuna.trial.FrozenTrial) -> None:
        # 初回呼び出しで開始時刻を記録
        if state['start_time'] is None:
            state['start_time'] = time.monotonic()

        # trial duration (Optuna が datetime オブジェクトを持っている)
        dur_sec = 0.0
        if trial.datetime_start is not None and trial.datetime_complete is not None:
            dur_sec = (trial.datetime_complete - trial.datetime_start).total_seconds()
        state['recent_dur'].append(dur_sec)
        state['n_done'] += 1

        # feasibility / best 判定
        is_feas = bool(trial.user_attrs.get('is_feasible', False))
        if is_feas:
            state['n_feas'] += 1

        val = trial.value if trial.value is not None else float('inf')
        is_new_best = is_feas and val < state['prev_best']
        delta = state['prev_best'] - val if is_new_best else None
        if is_new_best:
            state['prev_best'] = val

        # status badge
        if is_new_best:
            badge = '★ BEST  '
        elif is_feas:
            badge = '✓ feas  '
        else:
            badge = '✗ infeas'

        # trial 値表示 (infeas で 10000 は省略表示)
        if val >= 9999.0:
            val_s = '   ----'
        else:
            val_s = f"{val:9.2f}"

        # 改善幅
        if delta is not None and delta < 1e9:
            delta_s = f" (Δ-{delta:.1f})"
        else:
            delta_s = ""

        # 失敗理由 (infeas のみ)
        if not is_feas:
            reason = _fmt_failure(trial.user_attrs.get('failure_reason', ''))
            reason_s = f"  reason={reason}" if reason else ""
        else:
            reason_s = ""

        # 1 行目
        line0 = (
            f"[#{trial.number:03d}] {badge}  TAC={val_s}{delta_s}"
            f"{reason_s}  {dur_sec:5.1f}s"
        )

        # vars 3 行
        vars_lines = _fmt_vars(trial.params)

        # progress 行 (ETA は直近 trial duration の中央値で算出 = ばらつきに頑健)
        elapsed = time.monotonic() - state['start_time']
        n = state['n_done']
        if state['recent_dur']:
            sorted_d = sorted(state['recent_dur'])
            median_dur = sorted_d[len(sorted_d) // 2]
        else:
            median_dur = 0.0
        remaining = max(n_trials_total - n, 0)
        eta = remaining * median_dur if median_dur > 0 else 0.0

        pct = 100.0 * n / max(n_trials_total, 1)
        filled = int(_BAR_W * n / max(n_trials_total, 1))
        bar = '█' * filled + '░' * (_BAR_W - filled)
        feas_pct = 100.0 * state['n_feas'] / max(n, 1)
        best_s = f"{state['prev_best']:.2f}" if state['prev_best'] < 1e9 else '----'

        progress = (
            f"       [{bar}] {n}/{n_trials_total} ({pct:.0f}%)  "
            f"feas {state['n_feas']}/{n} ({feas_pct:.0f}%)  "
            f"elapsed {_fmt_dur(elapsed)}  ETA {_fmt_dur(eta)}  "
            f"pace {median_dur:.1f}s/trial  best {best_s}"
        )

        # 一括出力 (flush で tqdm との競合を避ける)
        print(line0, flush=True)
        for ln in vars_lines:
            print(ln, flush=True)
        print(progress, flush=True)
        print('', flush=True)  # trial 間の空行で区切り

    return _callback
