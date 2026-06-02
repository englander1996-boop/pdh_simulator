r"""
sub/sub2.py — (SM, rigorous, SM) バックエンドでの PDH 全変数 制約付き BO (マルチコア並列)

並列対応アーカイブ (parallel kind='sub2')。現行 BO の本丸は ../main.py (HYSYS+SM)。
exp4 / exp_tin_sweep / exp_dist2_pressure / exp_recost_359 は本ファイル(sub2)を import して
感度解析に使う。

main.py (sm, hysys, sm) のフォーク。**Dist2 を HYSYS → rigorous に置換**することで全塔が
pure Python になり、HYSYS の単一 COM 制約から解放されて **マルチプロセス並列** (N_WORKERS>1)
が可能になる。Dist1/Dist3 は学習済み SM(GPR)、Dist2 は in-house rigorous (Wang-Henke)。
Stage2(HEN 合成)を **全 trial で実行** (apply_stage2=True) するので、BO の目的関数が
「本物のプラント TAC」= top-k 再評価とほぼ一致する。

構成の比較:
  - main (sm,hysys,sm): 精度◎だが HYSYS COM で並列不可・遅い。
  - sub1 (fug,rig,fug→top-k rig): 並列可だが BO ループの Dist1/Dist3 が fug(低精度)、
    Stage2 は top-k のみ。BO↔再評価に乖離。
  - 本ファイル sub2 (sm,rig,sm)+Stage2 全 trial: SM が Dist1/Dist3 を HYSYS 精度・fug 速度で置換、
    Dist2 は精度保護のため rigorous 据え置き、全部 pure Python ゆえ 6 worker 並列可。
    → main の精度構造 × sub1 の並列性 × Stage2 を毎回。

Dist1 SM レンジ事情:
  SM column1 は N=30-60 で学習。本ファイルは col1_n_stages=(30,60) を採用し SM 学習域に合わせる。

21 変数 (Dist2 ブロックは rigorous 用に組み替え):
  反応器(4): T_in_K, z_cat_m, t_cyc_min, D_reactor_m
  PSA(3)   : D_psa_col_m, L_psa_bed_m, desorption_target
  膜(2)    : P_H_Pa, A_mem_m2
  原料(1)  : F_C3H8_fresh_kmol_h
  Dist1(4) : col1_p_kpa, col1_n_stages, col1_feed_stage, col1_comp_frac_2          (SM)
  Dist2(5) : col2_p_kpa, col2_n_stages, col2_reflux_ratio, col2_rec_LK, col2_rec_HK (rigorous)
  Dist3(3) : col3_p_kpa, col3_n_stages, col3_feed_ratio                            (SM)

出力 (outputs/final_<ts>/ にまとめる):
  trials.csv / best.json / top{1..N}_*.txt(CAPEX/OPEX/spec/HI 内訳) / README.md /
  optuna.db / _worker{0..N}.log(並列時の worker 別ライブログ)
  進捗の俯瞰は別ターミナルで:  python tools\monitor_main.py   /   tools\top_designs.py

使い方:  .\.venv\Scripts\python.exe sub\sub2.py > outputs\sub2_run.log 2>&1
"""

import os
import sys
import csv
import json
import time
import datetime
import contextlib
import io

# per-trial 時間予算 200s。Stage2(HEN) を全 trial 実行すると median~115s/max~170s かかり、
# 短いと過半が timeout する。200s で max~170s をカバーし timeout 率を <~20% に抑える。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '200')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

# repo root は 1 階層上。
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
N_TRIALS    = 360            # 6 worker × 60。全21変数。
N_STARTUP   = 50             # QMC Sobol 広域カバレッジ (共有 study 総完了数で TPE 切替)
SEED        = 42
N_JOBS      = 1              # スレッド並列は penalty_scale global のため不可。並列は N_WORKERS。
# 全塔 pure Python (SM + rigorous、HYSYS 不使用) なのでマルチプロセス並列が可能。各 worker
# 単スレッドで penalty_scale/GIL 問題なし、constant_liar=True で冗長サンプリング抑制 (parallel.py)。
# QMC startup フェーズは Sobol が結果非依存ゆえ並列で品質ゼロ劣化。TPE フェーズの同時数 6 が
# async バッチ=staleness の安全圏上端。1 で単一プロセス。
N_WORKERS   = 6
N_TOPK      = 3              # 詳細レポートを出す上位候補数

STUDY_PREFIX = "pdh"         # study 名 = pdh_<ts> (monitor_main.py / top_designs.py が認識)


# ===========================================================================
# § 2. 固定値 / 評価オプション
# ===========================================================================
P_L_Pa       = 1.0e5         # 膜透過側圧力 (大気圧固定)
APPLY_HI     = True
APPLY_STAGE2 = True          # ★ 全 trial で Stage2(HEN 合成) = 本物のプラント TAC を最適化
HI_DT_MIN_K  = 10.0


# ===========================================================================
# § 3. 探索範囲 — 21 変数  形式: (low, high, scale, type)
# ===========================================================================
SEARCH_SPACE = {
    # ----- 反応器 (Swing) -----
    "T_in_K":              (880.0,  940.0,  'linear', 'float'),
    "z_cat_m":             (15.0,   40.0,   'linear', 'float'),
    "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
    "D_reactor_m":         (7.0,    10.0,   'linear', 'float'),

    # ----- PSA -----
    "D_psa_col_m":         (2.9,    5.0,    'linear', 'float'),
    "L_psa_bed_m":         (22.0,   30.0,   'linear', 'float'),
    "desorption_target":   (0.22,   0.40,   'linear', 'float'),

    # ----- 膜 -----
    "P_H_Pa":              (7.5e5,  9.5e5,  'linear', 'float'),
    "A_mem_m2":            (5.0e4,  3.0e5,  'log',    'float'),

    # ----- 原料 (外側ループ skip、override) -----
    # Dist2 が rigorous で収率挙動が変わりうるため、生産量バンド (1129-1248) 達成状況を見て調整する。
    "F_C3H8_fresh_kmol_h": (1500.0, 1750.0, 'linear', 'float'),

    # ----- Dist1 (SM model1: N30-60/P1600-2200/feed_stage10-39/CF0.9-0.999) -----
    "col1_p_kpa":          (1600.0, 2000.0, 'linear', 'float'),
    "col1_n_stages":       (30,     60,     'linear', 'int'  ),   # SM 学習域
    "col1_feed_stage":     (22,     28,     'linear', 'int'  ),
    "col1_comp_frac_2":    (0.90,   0.999,  'linear', 'float'),

    # ----- Dist2 (rigorous 脱エタン塔) -----
    # rigorous 用域 (rigorous は N20-40 で収束、recovery spec で分離を直接制御)。
    # feed 段は Kirkbride 自動 (N_feed=1 placeholder)。
    "col2_p_kpa":          (500.0,  800.0,  'linear', 'float'),  # 5-8 bar (膜 P_H 未満を維持)
    "col2_n_stages":       (20,     40,     'linear', 'int'  ),
    "col2_reflux_ratio":   (6.0,    10.0,   'linear', 'float'),
    "col2_rec_LK_top":     (0.95,   0.999,  'linear', 'float'),  # C2H6 → top
    "col2_rec_HK_bot":     (0.998,  0.9995, 'linear', 'float'),  # C3H8 → bot

    # ----- Dist3 (SM model3: N69-200/P1600-2200, spec なし) -----
    "col3_p_kpa":          (1600.0, 1900.0, 'linear', 'float'),
    "col3_n_stages":       (115,    160,    'linear', 'int'  ),
    "col3_feed_ratio":     (0.60,   0.90,   'linear', 'float'),
}

# SM feed_stage 絶対 bounds (ratio→絶対段変換後に clamp)
_FEED_STAGE_ABS = {
    "col3": (70, 180),    # model3 feas 域 (feed_stage ≥70)
}

OUTPUT_DIR = 'outputs'


# ===========================================================================
# ↓↓↓ 以下はパイプライン (通常触らない) ↓↓↓
# ===========================================================================

import dataclasses as _dc
_CONFIG = load_operating_config()
# SM Dist3 は 99.5 mol%=99.497 wt% 固定。mol↔wt 差 0.003pp を吸収するため
# purity 閾値を 99.45 wt% に緩和。
_CONFIG = _dc.replace(_CONFIG, spec=_dc.replace(_CONFIG.spec, c3h6_min_wtfrac=0.9945))

_N_FEED_PLACEHOLDER = 1   # rigorous/sm は core 側で Kirkbride 推奨を自動採用、本値は無視


def _feed_stage_from_ratio(ratio: float, n: int, lo: int, hi: int) -> int:
    fs = int(round(ratio * n))
    hi_eff = min(hi, n - 2)
    lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))


def _suggest_params(trial: optuna.trial.Trial) -> dict:
    """SEARCH_SPACE 全変数を suggest して params dict を返す (固定範囲 = multivariate TPE 活用)。"""
    p: dict = {}
    for name, (low, high, scale, typ) in SEARCH_SPACE.items():
        if typ == 'int':
            p[name] = trial.suggest_int(name, int(low), int(high), log=(scale == 'log'))
        else:
            p[name] = trial.suggest_float(name, float(low), float(high), log=(scale == 'log'))
    return p


def _build_design(p: dict) -> FlowsheetDesignVars:
    """params dict (21 変数) から FlowsheetDesignVars を構築。trial 非依存 (best 再評価でも使用)。

    Dist1/Dist3 = SM、Dist2 = rigorous (recovery spec + Kirkbride feed)。
    """
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages'])
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
        # Dist2 = rigorous: reflux + recovery spec で分離を制御、feed は Kirkbride 自動。
        dist2=ColumnTunables(
            P_col=float(p['col2_p_kpa']) * 1000.0, N_stages=n2, N_feed=_N_FEED_PLACEHOLDER,
            reflux_ratio=float(p['col2_reflux_ratio']), solver_method='rigorous',
            recovery_LK_top=float(p['col2_rec_LK_top']),
            recovery_HK_bot=float(p['col2_rec_HK_bot'])),
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
    """main.py 既定の連続制約 (feas / 生産量方向 / 反応器SV / PSA / 膜 / 塔 shortfall)。"""
    return _default_constraints_func(trial)


# ===========================================================================
# ライブログ用 compact callback (flush 付き、5 行 + 進捗 + ETA + 失敗 tally)
# ===========================================================================
from collections import Counter as _Counter, deque as _deque
import time as _time
from optimization.callbacks import _fmt_dur, _fmt_reason_from_trial, _fmt_tally

_BAR_W = 30


def _make_final_callback(n_total: int):
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
        v2 = (f"Dist2(rig): P={p.get('col2_p_kpa',0):.0f}kPa N={p.get('col2_n_stages',0)} "
              f"R={p.get('col2_reflux_ratio',0):.2f} recLK={p.get('col2_rec_LK_top',0):.4f} "
              f"recHK={p.get('col2_rec_HK_bot',0):.4f}")
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
    """上位 top_n 候補を再評価して exp3 形式の詳細レポート (CAPEX/OPEX/spec/HI 内訳) を保存。"""
    comp = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
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
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                print(f"# final.py top{rank}  trial #{t.number}  "
                      f"effective_TAC(BO)={t.value:.2f} 億円/年  "
                      f"feasible={t.user_attrs.get('is_feasible')}")
                print("#" + "=" * 70)
                show_input_snapshot(design, _CONFIG, eval_kwargs)
                display_full_results(res, design, _CONFIG)
            report_text = buf.getvalue()
            path = f"{base_path}_top{rank}_trial{t.number}.txt"
            with open(path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            saved.append(path)
            if rank == 1:
                print("\n" + "=" * 72, flush=True)
                print(f"  ★ ベスト候補 詳細レポート (top1 / trial #{t.number}) ─ コンソール表示",
                      flush=True)
                print("=" * 72, flush=True)
                print(report_text, flush=True)
            print(f"  top{rank} 詳細レポート(CAPEX/OPEX/spec内訳): {path}", flush=True)
        except Exception as e:
            print(f"  top{rank} レポート生成失敗 (trial #{t.number}): {type(e).__name__}: {e}", flush=True)
    return saved


def _write_readme(out_dir: str, ts: str, best, n_complete: int, n_feas: int) -> None:
    """run subdir に README.md (結果の見方ガイド) を出力。"""
    L = []
    L.append(f"# PDH final.py run — {ts}  (SM, rigorous, SM) + Stage2 全trial / {N_WORKERS} worker 並列")
    L.append("")
    L.append("## まず見るべきファイル")
    L.append("1. **`top1_trial*.txt`** ─ ★最終ベストの全詳細 (CAPEX/OPEX/spec/HI/HEN 内訳)。")
    L.append("2. **`best.json`** ─ ベスト trial の params + 経済値 (再現用)。")
    L.append("3. **`trials.csv`** ─ 全 trial 履歴 (収率/N_d3/利益の分析は tools/top_designs.py が便利)。")
    L.append("")
    L.append("## 進捗の俯瞰 (run 中・別ターミナル)")
    L.append("```")
    L.append("python tools\\monitor_main.py            # 全worker合算の進捗・ETA")
    L.append("python tools\\top_designs.py --k 10      # feasible を 収率/N_d3/利益 で一覧")
    L.append("```")
    L.append("")
    L.append("## この run の構成")
    L.append(f"- 塔: Dist1=SM(N30-60), Dist2=rigorous, Dist3=SM  /  Stage2(HEN)=全trial実行")
    L.append(f"- N_TRIALS={N_TRIALS}, N_STARTUP={N_STARTUP}, N_WORKERS={N_WORKERS}, seed={SEED}")
    L.append(f"- 完了 trial={n_complete}, feasible={n_feas}")
    if best is not None:
        L.append(f"- ベスト: trial #{best.number}, effective_TAC={best.value:.3f} 億円/年")
    L.append("")
    with open(os.path.join(out_dir, 'README.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


def main():
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = os.path.join(OUTPUT_DIR, f'final_{ts}')   # run 専用サブディレクトリ
    os.makedirs(out_dir, exist_ok=True)
    print("=" * 72, flush=True)
    print(f"==== final.py: PDH (SM, rigorous, SM) 全21変数 制約付き BO ====", flush=True)
    print(f"  N_TRIALS={N_TRIALS}, N_STARTUP(QMC)={N_STARTUP}, seed={SEED}, "
          f"N_WORKERS={N_WORKERS}, top-k report={N_TOPK}", flush=True)
    print(f"  Dist1/Dist3 = SM, Dist2 = rigorous, Stage2(HEN) = 全trial / 上流も変数化", flush=True)
    print(f"  出力先: {out_dir}/", flush=True)
    print(f"  詳細ログ: {out_dir}/_worker*.log  (coordinator は 60 秒毎に集約進捗を表示)", flush=True)
    print(f"  進捗(別ターミナル): .\\.venv\\Scripts\\python.exe tools\\monitor_main.py", flush=True)
    print(f"                      .\\.venv\\Scripts\\python.exe tools\\top_designs.py --k 10", flush=True)
    print("=" * 72, flush=True)

    sampler = make_sampler('tpe', SEED, N_STARTUP, constraints_func=constraints_func)
    # 並列(N_WORKERS>1)時は worker 間で study 共有のため SQLite 必須。
    _para = N_WORKERS > 1
    study_name = f"{STUDY_PREFIX}_{ts}"
    _db_path = os.path.join(out_dir, 'optuna.db')
    storage = f"sqlite:///{_db_path}" if _para else None
    study = optuna.create_study(
        study_name=study_name, sampler=sampler, direction='minimize',
        storage=storage, load_if_exists=bool(storage),
    )

    t_start = datetime.datetime.now()
    if _para:
        print(f"  並列実行: {N_WORKERS} worker (各々 pure Python、HYSYS なし)。"
              f"worker ログ: {out_dir}/_worker*.log", flush=True)
        from optimization.parallel import spawn_workers
        spawn_workers(
            kind='sub2', study_name=study_name, storage_url=storage,
            db_path=_db_path, n_workers=N_WORKERS, n_trials_total=N_TRIALS,
            n_startup=N_STARTUP, base_seed=SEED, out_dir=out_dir,
        )
    else:
        run_optimization(
            study, objective, n_trials=N_TRIALS,
            show_progress_bar=False, n_jobs=N_JOBS,
            callbacks=[_make_final_callback(N_TRIALS)],
        )
    t_bo = (datetime.datetime.now() - t_start).total_seconds()

    # ---- 結果サマリ ----
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    print(f"\n==== 結果 ({t_bo:.0f}秒) ====", flush=True)
    print(f"  完了 trial: {len(complete)} / feasible: {len(feasible)}", flush=True)
    best = None
    if feasible:
        best = min(feasible, key=lambda t: t.value); tag = "feasible best"
    elif complete:
        best = min(complete, key=lambda t: t.value); tag = "best (feasible 無し)"
    if best is not None:
        print(f"  {tag}: trial #{best.number}  effective_TAC={best.value:.2f} 億円/年", flush=True)
        try:
            _pur  = float(best.user_attrs.get('c3h6_purity_wtfrac'))
            _prod = float(best.user_attrs.get('production_kmol_h'))
            _ff   = float(best.params.get('F_C3H8_fresh_kmol_h'))
            print(f"    purity={_pur*100:.2f}wt%  prod={_prod:.1f}kmol/h  "
                  f"F_fresh={_ff:.1f}kmol/h  yield={_prod/_ff*100:.1f}%", flush=True)
        except Exception:
            pass

    # ---- 保存 ----
    base = os.path.join(out_dir, 'final')
    _save_trials_csv(study, os.path.join(out_dir, 'trials.csv'))
    print(f"\n  trial 履歴 CSV: {out_dir}/trials.csv", flush=True)
    if best is not None:
        with open(os.path.join(out_dir, 'best.json'), 'w', encoding='utf-8') as f:
            json.dump({'number': best.number, 'effective_TAC': best.value,
                       'params': best.params,
                       'user_attrs': {k: v for k, v in best.user_attrs.items()}},
                      f, ensure_ascii=False, indent=2, default=str)
        print(f"  best JSON: {out_dir}/best.json", flush=True)
        print(f"\n  上位 {N_TOPK} 候補の詳細レポートを生成中...", flush=True)
        _save_best_reports(study, base, N_TOPK)

    _write_readme(out_dir, ts, best, len(complete), len(feasible))
    print(f"\n  README: {out_dir}/README.md", flush=True)
    print("=" * 72, flush=True)
    print(f"完了。成果物: {os.path.abspath(out_dir)}/", flush=True)
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
