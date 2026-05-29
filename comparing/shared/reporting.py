r"""
comparing.shared.reporting — main.py (旧 special.py) の I/O 機構の忠実移植。

main.py (旧 special.py) 内のモジュール関数 (_make_main_callback / _save_trials_csv /
_save_best_reports / _write_readme) を、ファイル改名 (2026-05-29: 旧 special.py→main.py) に巻き込まれないよう
ここへ忠実に移植する。ライブログの体裁・CSV 列・top-N 詳細レポート・README ガイドの
出力形式は special と同一に保ち、エラー調査・分析のしやすさを引き継ぐ。

依存 (流用):
  optimization.callbacks : _fmt_dur, _fmt_reason_from_trial, _fmt_tally
  simulation             : display_full_results, show_input_snapshot
  comparing.shared.space : build_design (best 再評価用)
  comparing.shared.simulator : raw_evaluate, CONFIG, EVAL_KWARGS_DEFAULT
"""

import os
import io
import csv
import json
import contextlib
from collections import Counter as _Counter, deque as _deque
import time as _time

import optuna

from optimization.callbacks import _fmt_dur, _fmt_reason_from_trial, _fmt_tally
from simulation import display_full_results, show_input_snapshot

from comparing.shared import space
from comparing.shared import simulator

_BAR_W = 30


# ===========================================================================
# ライブログ callback (special._make_special_callback の忠実移植)
# ===========================================================================
def make_callback(n_total: int):
    """trial 毎の compact ライブログ callback を生成 (special._make_special_callback と同形式)。

    1 trial = status 行 + 変数 4 行 (Rx/PSA/Mem/F, Dist1 SM, Dist2 HY, Dist3 SM) +
    purity/prod/yield 行 + progress バー (feas率/elapsed/ETA/pace/best) + top fails tally。
    flush=True なのでリダイレクトでも live に書き出される。
    """
    state = {'start': None, 'prev_best': float('inf'), 'n_feas': 0,
             'n_done': 0, 'recent': _deque(maxlen=20), 'tally': _Counter()}

    def cb(study, trial):
        if state['start'] is None:
            state['start'] = _time.monotonic()
        dur = 0.0
        if trial.datetime_start is not None and trial.datetime_complete is not None:
            dur = (trial.datetime_complete - trial.datetime_start).total_seconds()
        state['recent'].append(dur); state['n_done'] += 1
        a = trial.user_attrs
        is_feas = bool(a.get('is_feasible', False))
        if is_feas:
            state['n_feas'] += 1
        fu = a.get('failure_unit', '') or ('legacy' if not is_feas else '')
        if fu:
            state['tally'][fu] += 1
        val = trial.value if trial.value is not None else float('inf')
        new_best = is_feas and val < state['prev_best']
        delta = state['prev_best'] - val if new_best else None
        if new_best:
            state['prev_best'] = val
        badge = '★ BEST  ' if new_best else ('✓ feas  ' if is_feas else '✗ infeas')
        val_s = '   ----' if val >= 9999.0 else f"{val:9.2f}"
        delta_s = f" (Δ-{delta:.1f})" if (delta is not None and delta < 1e9) else ""
        reason_s = ("  reason=" + _fmt_reason_from_trial(trial)) if not is_feas else ""
        line0 = f"[#{trial.number:03d}] {badge}  TAC={val_s}{delta_s}{reason_s}  {dur:5.1f}s"

        p = trial.params
        v0 = (f"Rx: T={p.get('T_in_K',0):.0f}K z={p.get('z_cat_m',0):.1f} t={p.get('t_cyc_min',0):.1f} "
              f"D={p.get('D_reactor_m',0):.2f} | PSA: D={p.get('D_psa_col_m',0):.2f} L={p.get('L_psa_bed_m',0):.1f} "
              f"des={p.get('desorption_target',0):.3f} | Mem: P_H={p.get('P_H_Pa',0)/1e5:.2f}bar "
              f"A={p.get('A_mem_m2',0):.2e} | F={p.get('F_C3H8_fresh_kmol_h',0):.0f}")
        v1 = (f"Dist1(SM): P={p.get('col1_p_kpa',0):.0f}kPa N={p.get('col1_n_stages',0)} "
              f"feed={p.get('col1_feed_stage',0)} cf={p.get('col1_comp_frac_2',0):.3f}")
        v2 = (f"Dist2(HY): P={p.get('col2_p_kpa',0):.0f}kPa N={p.get('col2_n_stages',0)} "
              f"fr={p.get('col2_feed_ratio',0):.2f} R={p.get('col2_reflux_ratio',0):.1f}")
        v3 = (f"Dist3(SM): P={p.get('col3_p_kpa',0):.0f}kPa N={p.get('col3_n_stages',0)} "
              f"fr={p.get('col3_feed_ratio',0):.2f}")
        pur = a.get('c3h6_purity_wtfrac'); prod = a.get('production_kmol_h')
        fused = a.get('F_C3H8_fresh_used_kmol_h')
        outs = ""
        if pur and prod:
            outs = f"       -> purity={float(pur)*100:.2f}wt% prod={float(prod):.0f}kmol/h"
            if fused:
                outs += f" yield={float(prod)/float(fused)*100:.1f}%"

        elapsed = _time.monotonic() - state['start']; n = state['n_done']
        med = (sorted(state['recent'])[len(state['recent']) // 2] if state['recent'] else 0.0)
        eta = max(n_total - n, 0) * med
        pct = 100.0 * n / max(n_total, 1)
        filled = int(_BAR_W * n / max(n_total, 1))
        bar = '█' * filled + '░' * (_BAR_W - filled)
        feas_pct = 100.0 * state['n_feas'] / max(n, 1)
        best_s = f"{state['prev_best']:.2f}" if state['prev_best'] < 1e9 else '----'
        prog = (f"       [{bar}] {n}/{n_total} ({pct:.0f}%)  feas {state['n_feas']}/{n} "
                f"({feas_pct:.0f}%)  elapsed {_fmt_dur(elapsed)} ETA {_fmt_dur(eta)} "
                f"pace {med:.1f}s best {best_s}")
        tally_s = _fmt_tally(state['tally'], 5)

        print(line0, flush=True)
        print("       " + v0, flush=True)
        print("       " + v1, flush=True)
        print("       " + v2, flush=True)
        print("       " + v3, flush=True)
        if outs:
            print(outs, flush=True)
        print(prog, flush=True)
        if tally_s:
            print(f"       top fails: {tally_s}", flush=True)
        print('', flush=True)

    return cb


# ===========================================================================
# 結果サマリ
# ===========================================================================
def summarize(study: optuna.Study):
    """(complete, feasible, best) を返す。best は feasible 優先、無ければ TAC 最小。"""
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
                and t.value is not None]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    best = None
    if feasible:
        best = min(feasible, key=lambda t: t.value)
    elif complete:
        best = min(complete, key=lambda t: t.value)
    return complete, feasible, best


# ===========================================================================
# trials.csv (special._save_trials_csv の忠実移植)
# ===========================================================================
def save_trials_csv(study: optuna.Study, path: str) -> None:
    param_keys: list = []
    attr_keys: list = []
    for t in study.trials:
        for k in t.params:
            if k not in param_keys:
                param_keys.append(k)
        for k in t.user_attrs:
            if k not in attr_keys:
                attr_keys.append(k)
    header = ['number', 'value', 'state'] + param_keys + [f'attr.{k}' for k in attr_keys]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        for t in study.trials:
            row = [t.number, (t.value if t.value is not None else ''), t.state.name]
            row += [t.params.get(k, '') for k in param_keys]
            row += [t.user_attrs.get(k, '') for k in attr_keys]
            w.writerow(row)


# ===========================================================================
# top-N 詳細レポート (special._save_best_reports の忠実移植)
# ===========================================================================
def save_best_reports(study: optuna.Study, out_dir: str, top_n: int,
                      *, eval_kwargs: dict = None, build_design=None,
                      console_top1: bool = True) -> list:
    """上位 top_n の feasible 候補を再評価し CAPEX/OPEX/spec/HI 内訳を top{rank}_trial{N}.txt に保存。

    special._save_best_reports と同形式。feasible が無ければスキップ (penalty result でのクラッシュ回避)。
    eval_kwargs 既定 = simulator.EVAL_KWARGS_DEFAULT。build_design 既定 = space.build_design。
    """
    eval_kwargs = eval_kwargs or dict(simulator.EVAL_KWARGS_DEFAULT)
    build_design = build_design or space.build_design

    comp = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
    if not feas:
        print("  feasible trial が無いため詳細レポートはスキップ (CSV/JSON は保存済み)", flush=True)
        return []
    cand = sorted(feas, key=lambda t: t.value)[:top_n]
    saved = []
    for rank, t in enumerate(cand, 1):
        try:
            design = build_design(t.params)
            F_fresh = float(t.params.get('F_C3H8_fresh_kmol_h'))
            res = simulator.raw_evaluate(design, F_fresh=F_fresh, **eval_kwargs)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                print(f"# comparing top{rank}  trial #{t.number}  "
                      f"effective_TAC={t.value:.2f} 億円/年  "
                      f"feasible={t.user_attrs.get('is_feasible')}")
                print("#" + "=" * 70)
                show_input_snapshot(design, simulator.CONFIG, eval_kwargs)
                display_full_results(res, design, simulator.CONFIG)
            report_text = buf.getvalue()
            path = os.path.join(out_dir, f"top{rank}_trial{t.number}.txt")
            with open(path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            saved.append(path)
            if rank == 1 and console_top1:
                print("\n" + "=" * 72, flush=True)
                print(f"  ★ ベスト候補 詳細レポート (top1 / trial #{t.number}) ─ コンソール表示", flush=True)
                print("=" * 72, flush=True)
                print(report_text, flush=True)
            print(f"  top{rank} 詳細レポート(CAPEX/OPEX/spec内訳): {path}", flush=True)
        except Exception as e:
            print(f"  top{rank} レポート生成失敗 (trial #{t.number}): {type(e).__name__}: {e}", flush=True)
    return saved


# ===========================================================================
# README ガイド (special._write_readme の忠実移植 + 問題 P の記述)
# ===========================================================================
def write_readme(out_dir: str, *, method: str, p_codes: str, description: str,
                 study: optuna.Study, best, saved_reports: list,
                 settings: dict, extra_lines: list = None) -> None:
    """run subdir に README.md (結果の見方ガイド) を出力。special._write_readme と同形式 + 再現した問題の記述。"""
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    top1_name = os.path.basename(saved_reports[0]) if saved_reports else None

    L = []
    L.append(f"# {method} — 再現した問題: {p_codes}")
    L.append("")
    L.append(description.strip())
    L.append("")
    L.append("## まず見るべきファイル (推奨順)")
    L.append("")
    if top1_name:
        L.append(f"1. **`{top1_name}`** ─ ★この手法の best 候補の詳細 "
                 "(CAPEX/OPEX/spec/HI 内訳 + 入力スナップショット)。")
        L.append(f"   - `top2_*` / `top3_*` は次点候補。同形式で比較できる。")
    else:
        L.append("1. **`top*_trial*.txt`** ─ ★best 詳細 (今回は feasible 無しで未生成)。")
    L.append("2. **`best.json`** ─ best trial の params + 診断値。BO (outputs/special_*) との突合用。")
    L.append("3. **`trials.csv`** ─ 全評価履歴。手法の探索軌跡を Excel/pandas で確認。")
    L.append("")
    L.append("## この手法の設定")
    L.append("")
    for k, v in settings.items():
        L.append(f"- {k} = {v}")
    L.append("")
    L.append("## ベスト要約")
    L.append("")
    L.append(f"- 完了評価 = {len(complete)} / feasible = {len(feasible)}")
    if best is not None:
        tag = "feasible ✓" if best.user_attrs.get('is_feasible', False) else "infeasible ✗ (feasible 無し、TAC 最小)"
        L.append(f"- ベスト: **trial #{best.number}** ({tag})")
        L.append(f"- effective_TAC = **{best.value:.3f}** 億円/年")
        try:
            _pur  = float(best.user_attrs.get('c3h6_purity_wtfrac'))
            _prod = float(best.user_attrs.get('production_kmol_h'))
            _ff   = float(best.params.get('F_C3H8_fresh_kmol_h'))
            L.append(f"- purity = {_pur*100:.2f} wt%, 生産量 = {_prod:.1f} kmol/h, "
                     f"F_fresh = {_ff:.1f} kmol/h, 収率 = {_prod/_ff*100:.1f}%")
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    else:
        L.append("- 完了評価なし")
    if extra_lines:
        L.append("")
        L.extend(extra_lines)
    L.append("")
    L.append("## BO との比較 (ユーザ側で実施)")
    L.append("")
    L.append("本手法の `best.json` と BO の `outputs/special_*/best.json` を突合し、")
    L.append("ΔTAC = TAC(本手法) − TAC(BO) [億円/年, %] を「この欠陥の損失」として算出する。")
    with open(os.path.join(out_dir, 'README.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
