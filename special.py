r"""
special.py — HYSYS + SM バックエンドでの PDH プロセス全変数最適化 (制約付き Optuna BO)

exp3.py を評価関数として、**全 21 設計変数**(反応器・PSA・膜・原料 + 蒸留塔 3 つ)を
Bayesian Optimization で探索する。main.py (FUG/rigorous, 全フローシート最適化) の HYSYS+SM 版。

設計判断 (2026-05-25): main.py の BO 成功インフラを移植 + HYSYS/SM 特性に適応。
  借用 (main.py): QMC→TPE 2相サンプラ・constraints_func(連続制約)・penalty_scale・
    _store_diagnostics・compact callback(flush 付きライブログ, callbacks.py)・
    詳細レポート保存(display_full_results)。
  適応 (HYSYS/SM): Dist1/Dist3=SM(学習済み GPR, ~瞬時)、Dist2=HYSYS。
    探索 bounds は main.py の forensic 値 ∩ SM 分類器 feasible 領域。
    純度は SM Dist3 が 99.5 mol%=99.497 wt% 固定 → spec を 99.45 wt% に緩和(決定A)。
    塔本体 CAPEX は provider 側で FUG と同式で計算済み(N/還流が CAPEX に効く)。

21 変数:
  反応器(4): T_in_K, z_cat_m, t_cyc_min, D_reactor_m
  PSA(3)   : D_psa_col_m, L_psa_bed_m, desorption_target
  膜(2)    : P_H_Pa, A_mem_m2   (P_L=1atm 固定、mem.P_dist=Dist3圧 同期)
  原料(1)  : F_C3H8_fresh_kmol_h
  Dist1(4) : col1_p_kpa, col1_n_stages, col1_feed_stage, col1_comp_frac_2  (SM)
  Dist2(4) : col2_p_kpa, col2_n_stages, col2_feed_ratio, col2_reflux_ratio (HYSYS)
  Dist3(3) : col3_p_kpa, col3_n_stages, col3_feed_ratio                    (SM, spec なし)

出力:
  outputs/special_<ts>_trials.csv     : 全 trial の params + 診断
  outputs/special_<ts>_best.json      : best trial 要約
  outputs/special_<ts>_top{1..N}_*.txt: 上位候補の詳細レポート (CAPEX/OPEX/spec/HI 内訳)
  stdout(リダイレクト推奨): compact callback による trial 毎ライブログ

使い方:  .\.venv\Scripts\python.exe special.py > outputs\special_run.log 2>&1
  (Python は 3.13。flush 付きログなのでリダイレクトでも live に書き出される。)
"""

import os
import sys
import csv
import json
import time
import datetime
import contextlib
from typing import Optional

os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '300')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import optuna

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate, FlowsheetResult
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

from optimization.study import make_sampler, run_optimization, _default_constraints_func
from optimization.objective import _store_diagnostics
from optimization.penalty_scale import set_scale, default_schedule
from simulation import display_full_results, show_input_snapshot


# ===========================================================================
# § 1. BO 設定
# ===========================================================================
N_TRIALS    = 300            # main.py 準拠。全21変数なので 300 推奨 (~1-1.5h 目安)
N_STARTUP   = 50             # QMC Sobol 広域カバレッジ (以降 TPE)
SEED        = 42
N_JOBS      = 1              # HYSYS COM + penalty_scale global のため 1
N_TOPK      = 3              # 詳細レポートを出す上位候補数

USE_SQLITE_STORAGE = False
STUDY_NAME         = "pdh_hysys_sm_special"


# ===========================================================================
# § 2. 固定値 / 評価オプション
# ===========================================================================
P_L_Pa       = 1.0e5         # 膜透過側圧力 (大気圧固定、main.py と同じ)
APPLY_HI     = True
APPLY_STAGE2 = True
HI_DT_MIN_K  = 10.0


# ===========================================================================
# § 3. 探索範囲 — 全 21 変数  形式: (low, high, scale, type)
#   bounds は main.py の forensic 調整値 ∩ SM 分類器 feasible 領域 (2026-05-25 プローブ)。
#   SM 崖 (予測不能域) は除外、解の縁 (収束ぎりぎり) は含める。
# ===========================================================================
SEARCH_SPACE = {
    # ----- 反応器 (Swing) — main.py 準拠 -----
    "T_in_K":              (880.0,  940.0,  'linear', 'float'),
    "z_cat_m":             (15.0,   40.0,   'linear', 'float'),
    "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
    "D_reactor_m":         (7.0,    10.0,   'linear', 'float'),

    # ----- PSA — main.py 準拠 -----
    "D_psa_col_m":         (2.9,    5.0,    'linear', 'float'),
    "L_psa_bed_m":         (22.0,   30.0,   'linear', 'float'),
    "desorption_target":   (0.22,   0.40,   'linear', 'float'),

    # ----- 膜 (P_L 固定、P_dist は Dist3 圧と同期) — main.py 準拠 -----
    "P_H_Pa":              (7.5e5,  9.5e5,  'linear', 'float'),
    "A_mem_m2":            (5.0e4,  3.0e5,  'log',    'float'),

    # ----- 原料 (外側ループ skip、override) — main.py 準拠 -----
    "F_C3H8_fresh_kmol_h": (1380.0, 1500.0, 'linear', 'float'),

    # ----- Dist1 (SM model1: N30-60/P1600-2200/feed_stage10-39/CF0.9-0.999) -----
    # feed_stage は SM feas ≥22 (プローブ: <21 で 0%)。範囲は (22,28) 固定: N 下限 30 でも
    # N-2=28 以下で常に有効になり、動的 search space を避けて TPE(multivariate) をフル活用。
    # (SM feas 域 ≤39 のうち 22-28 を採用。Dist1 は cost driver でないため探索損失は許容。)
    "col1_p_kpa":          (1600.0, 2000.0, 'linear', 'float'),
    "col1_n_stages":       (30,     60,     'linear', 'int'  ),
    "col1_feed_stage":     (22,     28,     'linear', 'int'  ),
    "col1_comp_frac_2":    (0.90,   0.999,  'linear', 'float'),

    # ----- Dist2 (HYSYS 脱エタン塔)。収束 envelope 狭、縁を含む (N≈45 安/N=75-80 頑健) -----
    "col2_p_kpa":          (500.0,  620.0,  'linear', 'float'),  # 膜 P_H 未満
    "col2_n_stages":       (44,     80,     'linear', 'int'  ),
    "col2_feed_ratio":     (0.40,   0.60,   'linear', 'float'),
    "col2_reflux_ratio":   (8.0,    13.0,   'linear', 'float'),

    # ----- Dist3 (SM model3: N69-200/P1600-2200, spec なし)。feas: N≥115, P≤1900 -----
    "col3_p_kpa":          (1600.0, 1900.0, 'linear', 'float'),
    "col3_n_stages":       (115,    160,    'linear', 'int'  ),
    "col3_feed_ratio":     (0.60,   0.90,   'linear', 'float'),
}

# SM/HYSYS の feed_stage 絶対 bounds (ratio→絶対段変換後に clamp)。
_FEED_STAGE_ABS = {
    "col2": (2, 9999),    # HYSYS は N-2 のみ
    "col3": (70, 180),    # model3 feas 域 (feed_stage ≥70)
}


# ===========================================================================
# § 4. 出力
# ===========================================================================
OUTPUT_DIR = 'outputs'


# ===========================================================================
# ↓↓↓ 以下はパイプライン (通常触らない) ↓↓↓
# ===========================================================================

import dataclasses as _dc
_CONFIG = load_operating_config()
# 決定A (2026-05-25): SM Dist3 は 99.5 mol%=99.497 wt% 固定。mol↔wt 差 0.003pp を吸収するため
# purity 閾値を 99.45 wt% に緩和 (SM の実力を尊重)。詳細は project_sm_integration メモ。
_CONFIG = _dc.replace(_CONFIG, spec=_dc.replace(_CONFIG.spec, c3h6_min_wtfrac=0.9945))


def _feed_stage_from_ratio(ratio: float, n: int, lo: int, hi: int) -> int:
    fs = int(round(ratio * n))
    hi_eff = min(hi, n - 2)
    lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))


def _suggest_params(trial: optuna.trial.Trial) -> dict:
    """SEARCH_SPACE 全変数を suggest して params dict を返す。

    全変数が固定範囲なので multivariate TPE をフル活用できる (動的 search space なし)。
    col2/col3 の feed は ratio で suggest し、_build_design で絶対段に変換 (N 依存 clamp)。
    """
    p: dict = {}
    for name, (low, high, scale, typ) in SEARCH_SPACE.items():
        if typ == 'int':
            p[name] = trial.suggest_int(name, int(low), int(high), log=(scale == 'log'))
        else:
            p[name] = trial.suggest_float(name, float(low), float(high), log=(scale == 'log'))
    return p


def _build_design(p: dict) -> FlowsheetDesignVars:
    """params dict (21 変数) から FlowsheetDesignVars を構築。trial 非依存 (best 再評価でも使用)。"""
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages']); fs2 = _feed_stage_from_ratio(p['col2_feed_ratio'], n2, *_FEED_STAGE_ABS['col2'])
    n3 = int(p['col3_n_stages']); fs3 = _feed_stage_from_ratio(p['col3_feed_ratio'], n3, *_FEED_STAGE_ABS['col3'])
    p3_kpa = float(p['col3_p_kpa'])
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=p['T_in_K'], z_cat=p['z_cat_m'],
                          t_cyc=p['t_cyc_min'], D=p['D_reactor_m']),
        psa=PSADesignVars(D_col=p['D_psa_col_m'], L_bed=p['L_psa_bed_m'],
                          desorption_target=p['desorption_target']),
        mem=MemDesignVars(P_H=p['P_H_Pa'], P_L=P_L_Pa, A_mem=p['A_mem_m2'],
                          P_dist=p3_kpa * 1000.0),
        dist1=ColumnTunables(
            P_col=float(p['col1_p_kpa']) * 1000.0, N_stages=n1, N_feed=1, reflux_ratio=2.0,
            solver_method='sm', hysys_spec_value=float(p['col1_comp_frac_2']), hysys_feed_stage=fs1),
        dist2=ColumnTunables(
            P_col=float(p['col2_p_kpa']) * 1000.0, N_stages=n2, N_feed=1,
            reflux_ratio=float(p['col2_reflux_ratio']),
            solver_method='hysys', hysys_spec_value=float(p['col2_reflux_ratio']), hysys_feed_stage=fs2),
        dist3=ColumnTunables(
            P_col=p3_kpa * 1000.0, N_stages=n3, N_feed=1, reflux_ratio=12.0,
            solver_method='sm', hysys_spec_value=0.99, hysys_feed_stage=fs3),  # spec は SM 未使用
    )


def objective(trial: optuna.trial.Trial) -> float:
    scale = default_schedule(trial.number, N_TRIALS)
    set_scale(scale)
    trial.set_user_attr('penalty_scale', scale)

    params = _suggest_params(trial)
    design = _build_design(params)
    F_fresh = float(params['F_C3H8_fresh_kmol_h'])

    t0 = time.perf_counter()
    result: FlowsheetResult = evaluate(
        design, _CONFIG, verbose=False,
        apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K, apply_stage2=APPLY_STAGE2,
        F_C3H8_override=F_fresh,
    )
    trial.set_user_attr('wallclock_sec', time.perf_counter() - t0)
    _store_diagnostics(trial, result)
    trial.set_user_attr('F_C3H8_fresh_used_kmol_h', F_fresh)
    return result.effective_TAC


def constraints_func(trial: optuna.trial.FrozenTrial):
    """main.py 既定の連続制約 (feas / 生産量方向 / 反応器SV / PSA / 膜 / 塔 shortfall)。

    全変数最適化なので反応器/PSA/膜の shortfall 信号も活性化し、TPE が上流の方向も学習する。
    純度は SM で不変のため制約化しない (定数=無意味)。
    """
    return _default_constraints_func(trial)


# ===========================================================================
# ライブログ用 compact callback (main.py の make_compact_callback 相当、flush 付き)
# ===========================================================================
from collections import Counter as _Counter, deque as _deque
import time as _time
from optimization.callbacks import _fmt_dur, _fmt_reason_from_trial, _fmt_tally

_BAR_W = 30


def _make_special_callback(n_total: int):
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
# レポート / 保存
# ===========================================================================

def _save_trials_csv(study: optuna.Study, path: str) -> None:
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


def _save_best_reports(study: optuna.Study, base_path: str, top_n: int) -> list:
    """上位 top_n 候補を再評価して exp3 形式の詳細レポート (CAPEX/OPEX/spec/HI 内訳) を保存。

    main.py の display_best_full / top-k レポート相当。feasible 優先、無ければ TAC 最小。
    """
    comp = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
    # 設計判断: display_full_results は economics 前提なので feasible のみ対象。
    # feasible が無ければスキップ (penalty result でのクラッシュ回避)。
    if not feas:
        print("  feasible trial が無いため詳細レポートはスキップ (CSV/JSON は保存済み)", flush=True)
        return []
    cand = sorted(feas, key=lambda t: t.value)[:top_n]
    eval_kwargs = dict(apply_hi=APPLY_HI, apply_stage2=APPLY_STAGE2, hi_dT_min_K=HI_DT_MIN_K)
    saved = []
    for rank, t in enumerate(cand, 1):
        try:
            design = _build_design(t.params)
            F_fresh = float(t.params.get('F_C3H8_fresh_kmol_h'))
            res = evaluate(design, _CONFIG, verbose=False,
                           apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
                           apply_stage2=APPLY_STAGE2, F_C3H8_override=F_fresh)
            path = f"{base_path}_top{rank}_trial{t.number}.txt"
            with open(path, 'w', encoding='utf-8') as f, contextlib.redirect_stdout(f):
                print(f"# special.py top{rank}  trial #{t.number}  "
                      f"effective_TAC(BO)={t.value:.2f} 億円/年  "
                      f"feasible={t.user_attrs.get('is_feasible')}")
                print("#" + "=" * 70)
                show_input_snapshot(design, _CONFIG, eval_kwargs)
                display_full_results(res, design, _CONFIG)
            saved.append(path)
            print(f"  top{rank} 詳細レポート(CAPEX/OPEX/spec内訳): {path}", flush=True)
        except Exception as e:
            print(f"  top{rank} レポート生成失敗 (trial #{t.number}): {type(e).__name__}: {e}", flush=True)
    return saved


def main():
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    print(f"==== special.py: PDH HYSYS+SM 全21変数 制約付き BO ====", flush=True)
    print(f"  N_TRIALS={N_TRIALS}, N_STARTUP(QMC)={N_STARTUP}, seed={SEED}, top-k report={N_TOPK}", flush=True)
    print(f"  Dist1/Dist3 = SM, Dist2 = HYSYS / 上流(反応器・PSA・膜)・F_fresh も変数化", flush=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    sampler = make_sampler('tpe', SEED, N_STARTUP, constraints_func=constraints_func)
    storage = None
    if USE_SQLITE_STORAGE:
        storage = f"sqlite:///{os.path.join(OUTPUT_DIR, f'special_{ts}.db')}"
    study = optuna.create_study(
        study_name=STUDY_NAME if not USE_SQLITE_STORAGE else f"{STUDY_NAME}_{ts}",
        sampler=sampler, direction='minimize',
        storage=storage, load_if_exists=bool(storage),
    )

    run_optimization(
        study, objective, n_trials=N_TRIALS,
        show_progress_bar=False, n_jobs=N_JOBS,
        callbacks=[_make_special_callback(N_TRIALS)],
    )

    # ---- 結果サマリ ----
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    print("\n==== 結果 ====", flush=True)
    print(f"  完了 trial: {len(complete)} / feasible: {len(feasible)}", flush=True)
    best = None
    if feasible:
        best = min(feasible, key=lambda t: t.value); tag = "feasible best"
    elif complete:
        best = min(complete, key=lambda t: t.value); tag = "best (feasible 無し)"
    if best is not None:
        print(f"  {tag}: trial #{best.number}  effective_TAC={best.value:.2f} 億円/年", flush=True)
        print(f"    purity={best.user_attrs.get('c3h6_purity_wtfrac','?')} "
              f"prod={best.user_attrs.get('production_kmol_h','?')} "
              f"F_fresh={best.params.get('F_C3H8_fresh_kmol_h','?')}", flush=True)
        for k, v in best.params.items():
            print(f"    {k} = {v}", flush=True)

    # ---- 保存 ----
    base = os.path.join(OUTPUT_DIR, f'special_{ts}')
    _save_trials_csv(study, base + '_trials.csv')
    print(f"\n  trial 履歴 CSV: {base}_trials.csv", flush=True)
    if best is not None:
        with open(base + '_best.json', 'w', encoding='utf-8') as f:
            json.dump({'number': best.number, 'effective_TAC': best.value,
                       'params': best.params,
                       'user_attrs': {k: v for k, v in best.user_attrs.items()}},
                      f, ensure_ascii=False, indent=2, default=str)
        print(f"  best JSON: {base}_best.json", flush=True)
        # 上位候補の詳細レポート (CAPEX/OPEX/spec/HI 内訳)
        print(f"\n  上位 {N_TOPK} 候補の詳細レポートを生成中...", flush=True)
        _save_best_reports(study, base, N_TOPK)

    try:
        from units.vle.hysys.provider import shutdown_default_provider
        shutdown_default_provider()
    except Exception:
        pass


if __name__ == "__main__":
    main()
