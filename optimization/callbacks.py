"""BO ループ用の表示 callback。

study.optimize の callbacks= に渡して、各 trial 完了時に compact な状態表示を行う。

Optuna デフォルトの logger は "Trial N finished with value: X and parameters: {…}" を
全 param 込みで 1 行に出すため可読性が低い。独自 callback で 1 trial = 4-5 行の構造化
された表示に置き換える。
    - 1 行目: status (★ BEST / ✓ feas / ✗ infeas) + TAC + reason + 経過秒
    - 2-4 行目: 全 design vars をユニット別にグループ化
    - 5 行目: progress (完了/全体, feasibility 率, elapsed, ETA, pace, top fails tally)

ETA は直近 N=20 trial の duration 中央値で計算 (累積平均は trial 間ばらつきで暴れる)。

失敗理由は trial.user_attrs の failure_unit + 各装置の penalty_reason + key actual 値を
読んで「Mem.bp_le_cold_out (T_bp=305<313)」形式に構造化表示する。さらに closure state に
Counter を持ち、top 5 failure mode を progress 行に追加表示することで「100 trial で
prod_under=58 件 → F_fresh 高すぎ」のようなパターンを走行中に把握できる。
"""

from __future__ import annotations

import time
from collections import Counter, deque
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


def _fmt_reason_from_trial(trial: optuna.trial.FrozenTrial) -> str:
    """failure_unit + 各装置の penalty_reason / actual 値から「どこで詰まったか」を構造化。

    優先順位:
      1. failure_unit が 'r_xxx' / 'spec_xxx' / 'solver_xxx' / 'timeout' / 'exception:...' の
         いずれかなら、その装置/段階の詳細情報を user_attrs から拾って「Unit.reason (key=val)」
         形式に整形。
      2. failure_unit が空 (legacy) なら failure_reason 文字列を 30 文字以内に圧縮。
      3. feasible (failure_unit='success') なら ''。

    例:
      'r_mem' + mem_penalty_reason='bp_le_cold_out' + T_bp=305 + T_cold=313
        → "Mem.bp_le_cold_out (T_bp=305K<313K)"
      'r_psa' + psa_penalty_reason='t_abs_below_min' + t_abs=37
        → "PSA.t_abs_below_min (t_abs=37s<60s)"
      'r_rx' + reactor_penalty_reason='sv_out_of_range' + SV=4.03
        → "Rx.sv_out_of_range (SV=4.03>3.0)"
      'spec_production_under' + production_kmol_h=1077 + target=1186 + under_pp=9.2
        → "spec.prod_under (1077/1186, -9.2pp)"
      'solver_inner_diverge'  → "solver_inner_diverge"
      'timeout'               → "timeout (>120s)"
      'exception:Mem'         → "exception:Mem"
    """
    a = trial.user_attrs
    fu = a.get('failure_unit', '') or ''

    if fu == '' or fu == 'success':
        return ''

    # ---- 蒸留塔 ----
    if fu in ('r1', 'r2', 'r3'):
        idx = fu[-1]
        msg = a.get(f'r{idx}_penalty_msg', '')
        n_needed = a.get(f'r{idx}_N_needed', 0.0) or 0.0
        dt_max   = a.get(f'r{idx}_dT_max_K', 0.0) or 0.0
        parts = []
        if n_needed > 0:
            parts.append(f"N_needed={n_needed:.0f}")
        if dt_max > 0:
            parts.append(f"dT={dt_max:.1f}K")
        if msg and not parts:
            # 数値出ない場合は msg を 40 文字に切る
            parts.append(str(msg)[:40])
        suffix = f" ({', '.join(parts)})" if parts else ""
        return f"Dist{idx}{suffix}"

    # ---- Reactor ----
    if fu == 'r_rx':
        reason = a.get('reactor_penalty_reason', '') or 'penalty'
        sv = a.get('reactor_SV_actual_m_s', 0.0) or 0.0
        if reason == 'sv_out_of_range' and sv > 0:
            # 範囲は run_one_pass._REACTOR_SV_MIN_MS=0.5 / _MAX_MS=3.0
            if sv < 0.5:
                return f"Rx.sv_out_of_range (SV={sv:.2f}<0.5)"
            if sv > 3.0:
                return f"Rx.sv_out_of_range (SV={sv:.2f}>3.0)"
        return f"Rx.{reason}"

    # ---- PSA ----
    if fu == 'r_psa':
        reason = a.get('psa_penalty_reason', '') or 'penalty'
        t_abs = a.get('psa_t_abs_actual_s', 0.0) or 0.0
        u_0   = a.get('psa_u_0_actual_m_s', 0.0) or 0.0
        if reason == 't_abs_below_min' and t_abs > 0:
            return f"PSA.t_abs_below_min (t_abs={t_abs:.0f}s<60s)"
        if reason == 'u_0_above_max' and u_0 > 0:
            return f"PSA.u_0_above_max (u_0={u_0:.2f}m/s>1.0)"
        if reason == 'mask_lt_2' and t_abs > 0:
            return f"PSA.mask_lt_2 (t_abs={t_abs:.0f}s)"
        return f"PSA.{reason}"

    # ---- Mem ----
    if fu == 'r_mem':
        reason = a.get('mem_penalty_reason', '') or 'penalty'
        if reason == 'bp_le_cold_out':
            t_bp = a.get('mem_T_bp_perm_actual_K', 0.0) or 0.0
            t_co = a.get('mem_T_cold_out_actual_K', 0.0) or 0.0
            if t_bp > 0 and t_co > 0:
                return f"Mem.bp_le_cold_out (T_bp={t_bp:.0f}K<{t_co:.0f}K+5)"
        elif reason == 'ph_le_pfeed':
            p_h    = a.get('mem_P_H_actual_Pa', 0.0) or 0.0
            p_feed = a.get('mem_P_feed_actual_Pa', 0.0) or 0.0
            if p_h > 0 and p_feed > 0:
                return f"Mem.ph_le_pfeed (P_H={p_h/1e5:.1f}<{p_feed/1e5:.1f}bar)"
        elif reason == 'vapor_condensed':
            t_dew  = a.get('mem_T_dew_actual_K', 0.0) or 0.0
            t_feed = a.get('mem_T_feed_actual_K', 0.0) or 0.0
            if t_dew > 0 and t_feed > 0:
                return f"Mem.vapor_condensed (T_in={t_feed:.0f}<T_dew={t_dew:.0f}K)"
        return f"Mem.{reason}"

    # ---- solver / timeout / guard ----
    if fu == 'timeout':
        import os as _os
        budget = _os.environ.get('PDH_TRIAL_TIME_BUDGET_SEC', '120')
        return f"timeout (>{budget}s)"
    if fu.startswith('solver_'):
        return fu
    if fu == 'solver_penalty':
        return "solver_penalty (unit unknown)"

    # ---- strict_recovery ----
    if fu.startswith('strict_recovery_'):
        return fu  # 例: strict_recovery_r2_lk

    # ---- spec ----
    if fu == 'spec_production_under':
        prod   = a.get('production_kmol_h', 0.0) or 0.0
        target = a.get('target_kmol_h',     0.0) or 0.0
        pp     = a.get('production_under_pp', 0.0) or 0.0
        if prod > 0 and target > 0:
            return f"spec.prod_under ({prod:.0f}/{target:.0f}kmol/h, -{pp:.1f}pp)"
        return "spec.prod_under"
    if fu == 'spec_production_over':
        prod   = a.get('production_kmol_h', 0.0) or 0.0
        target = a.get('target_kmol_h',     0.0) or 0.0
        pp     = a.get('production_over_pp', 0.0) or 0.0
        if prod > 0 and target > 0:
            return f"spec.prod_over ({prod:.0f}/{target:.0f}kmol/h, +{pp:.1f}pp)"
        return "spec.prod_over"
    if fu == 'spec_c3h6_purity':
        c3 = a.get('c3h6_purity_wtfrac', 0.0) or 0.0
        if c3 > 0:
            # 閾値は config 依存。ここは user_attrs しか持たず threshold を知らないので、
            # ハードコード比較を書かず実測値だけ示す (h2 行と同様)。
            return f"spec.c3h6_purity ({c3*100:.2f}wt%)"
        return "spec.c3h6_purity"
    if fu == 'spec_h2_purity':
        h2 = a.get('h2_purity_molfrac', 0.0) or 0.0
        if h2 > 0:
            return f"spec.h2_purity ({h2*100:.2f}mol%)"
        return "spec.h2_purity"

    # ---- exception ----
    if fu.startswith('exception:'):
        return fu  # 例: exception:Mem

    # ---- fallback (legacy: failure_unit が空 / 未知ラベル) ----
    reason = a.get('failure_reason', '')
    if not reason:
        return fu or ''
    r = str(reason)
    if '生産量' in r:
        return '生産量未達'
    if 'C3H6 純度' in r or '純度' in r:
        return '純度不足'
    if 'プロキシ' in r:
        return 'proxy 罰則'
    if 'trace bypass' in r:
        return 'trace bypass'
    if 'strict recovery' in r:
        return 'strict rec NG'
    if 'penalty' in r.lower() or '罰則' in r or 'CAPEX' in r:
        return 'solver/penalty'
    short = r.split('|')[0].strip()
    return short[:30] + ('…' if len(short) > 30 else '')


def _fmt_vars(params: dict, user_attrs: dict | None = None) -> list[str]:
    """全 design vars をユニット別 3 行にまとめる。

    user_attrs が渡されると 3 行目末尾に「yield=XX.X% (prod/fresh)」を追加する。
    yield = production_kmol_h / F_C3H8_fresh_used_kmol_h × 100%
    PDH 反応 C3H8 → C3H6 + H2 は 1:1 mol、production は C3H6 [kmol/h]、fresh は C3H8 [kmol/h]。
    工業的目安: per-pass conversion 30-50%、recycle 込み overall yield 85-90%。
    production または F_fresh が 0/未計算なら yield 表示は省略 (上流 fail で値が出ない場合)。
    """
    g = params.get

    # 反応器: 軸流 (z_cat_m/D_reactor_m) / Catofin (N_online/L_bed_m) / 径方向流
    # (D_inner_m/bed_thickness_m/H_m) を params のキーで判定して表示を切替。既定は Catofin。
    if 'z_cat_m' in params:
        rx = (
            f"Reactor(axial): T={g('T_in_K', 0):.0f}K z={g('z_cat_m', 0):.1f}m "
            f"t={g('t_cyc_min', 0):.1f}min D={g('D_reactor_m', 0):.2f}m"
        )
    elif 'N_online' in params:
        rx = (
            f"Reactor(catofin): T={g('T_in_K', 0):.0f}K t={g('t_cyc_min', 0):.1f}min "
            f"D={g('D_reactor_m', 0):.2f}m L_bed={g('L_bed_m', 0):.2f}m "
            f"N={g('N_online', 0)} dp={g('d_p_mm', 0):.1f}mm"
        )
    else:
        rx = (
            f"Reactor(radial): T={g('T_in_K', 0):.0f}K t={g('t_cyc_min', 0):.1f}min "
            f"Di={g('D_inner_m', 0):.1f}m dr={g('bed_thickness_m', 0):.2f}m H={g('H_m', 0):.1f}m"
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

    # 収率表示: user_attrs から production/F_fresh_used を取って表示。
    if user_attrs is not None:
        prod = user_attrs.get('production_kmol_h', 0.0) or 0.0
        f_used = user_attrs.get('F_C3H8_fresh_used_kmol_h', 0.0) or 0.0
        if prod > 0 and f_used > 0:
            yield_pct = prod / f_used * 100.0
            line3 += f"  yield={yield_pct:.1f}% ({prod:.0f}/{f_used:.0f})"

    return [line1, line2, line3]


def _fmt_tally(tally: Counter, top_k: int = 5) -> str:
    """failure_unit Counter の top-K を 'key=count key=count ...' 形式に。"""
    if not tally:
        return ''
    top = tally.most_common(top_k)
    return ' '.join(f"{k}={v}" for k, v in top)


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

    state は closure で保持する (グローバル汚染しない):
      - start_time : study 開始時刻 (time.monotonic)
      - prev_best  : 直前のベスト値 (改善幅算出)
      - n_feas     : これまでの feasible trial 数
      - recent_dur : 直近 N trial の所要時間 deque (ETA 中央値計算)
      - tally      : failure_unit 別の累計 Counter (top fails 表示用)
    """
    state = {
        'start_time': None,
        'prev_best':  float('inf'),
        'n_feas':     0,
        'n_done':     0,
        'recent_dur': deque(maxlen=20),
        'tally':      Counter(),
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

        # failure_unit を tally に積む。feasible でも failure_unit='success' で計上
        # (top fails には success も入るが _fmt_tally で見やすく)。
        fu = trial.user_attrs.get('failure_unit', '') or ''
        if fu:
            state['tally'][fu] += 1
        elif not is_feas:
            # legacy: failure_unit が無い (新コード以前の trial) → 'legacy' バケットに
            state['tally']['legacy'] += 1

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

        # 失敗理由 (infeas のみ、新フォーマット)
        if not is_feas:
            reason = _fmt_reason_from_trial(trial)
            reason_s = f"  reason={reason}" if reason else ""
        else:
            reason_s = ""

        # 1 行目
        line0 = (
            f"[#{trial.number:03d}] {badge}  TAC={val_s}{delta_s}"
            f"{reason_s}  {dur_sec:5.1f}s"
        )

        # vars 3 行 (3 行目末尾に yield 表示、user_attrs から計算)
        vars_lines = _fmt_vars(trial.params, trial.user_attrs)

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

        # top fails 累計を progress の次行に。
        # 例: top fails: r_mem=12  r_psa=8  spec_production_under=58  r_rx=5  timeout=3
        tally_s = _fmt_tally(state['tally'], top_k=5)
        tally_line = f"       top fails: {tally_s}" if tally_s else ""

        # 一括出力 (flush で tqdm との競合を避ける)
        print(line0, flush=True)
        for ln in vars_lines:
            print(ln, flush=True)
        print(progress, flush=True)
        if tally_line:
            print(tally_line, flush=True)
        print('', flush=True)  # trial 間の空行で区切り

    return _callback
