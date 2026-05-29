r"""
comparing.shared.harness — 素朴手法を「enqueue 駆動の Optuna study」として走らせる共通エンジン。

問題のあるやり方 (逐次1次元・グリッド・部分最適化・後置ピンチ) は「探索点をどう生成して
enqueue するか」で表現する。実行・記録は special.py と同じ Optuna 経路に乗せることで、
ライブログ・trials.csv・top-N 詳細レポート・README を忠実に流用する。

学習はしない (BO ではない) ので sampler は RandomSampler。各手法は **全 21 変数を enqueue** する
(掃引する変数は格子点、それ以外は固定ベースライン値) ため sampler は実際には引かれず決定的。

提供物:
  suggest_all(trial, control_space)         : 21 設計変数 (+制御変数) を suggest
  make_objective(...)                       : objective(trial) を生成 (_store_diagnostics 込み)
  new_study(seed)                           : RandomSampler の study
  run_batch(study, objective, points, cb)   : points を enqueue→評価し、その batch の trial を返す
  best_of(trials)                           : batch 内 best (feasible 優先, TAC 最小)
  finalize(study, ...)                      : trials.csv/best.json/top-N/README を保存 (special.main 流)
  new_run_dir(method)                       : comparing/results/<method>_<ts>/ を作る
  save_table_csv(rows, path, fieldnames)    : 手法固有の表 (コスト曲線・dTmin 掃引) を CSV 化
"""

import os
import json
import time
import datetime

import optuna

from optimization.objective import _store_diagnostics
from optimization.penalty_scale import set_scale
from optimization.study import make_sampler, run_optimization

from comparing.shared import space
from comparing.shared import simulator
from comparing.shared import reporting

optuna.logging.set_verbosity(optuna.logging.WARNING)

_RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results')

# 素朴手法の penalty scale は固定 1.0 (BO のような curriculum は使わない。
# default_schedule の中盤値に相当し、feasible best の TAC にはほぼ影響しない)。
PENALTY_SCALE = 1.0


def new_run_dir(method: str) -> str:
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    d = os.path.join(_RESULTS_DIR, f'{method}_{ts}')
    os.makedirs(d, exist_ok=True)
    return d


def suggest_all(trial: optuna.trial.Trial, control_space: dict = None) -> dict:
    """SEARCH_SPACE 全 21 変数 (+ control_space の float 制御変数) を suggest して dict で返す。"""
    control_space = control_space or {}
    p: dict = {}
    for name, (low, high, scale, typ) in space.SEARCH_SPACE.items():
        if typ == 'int':
            p[name] = trial.suggest_int(name, int(low), int(high), log=(scale == 'log'))
        else:
            p[name] = trial.suggest_float(name, float(low), float(high), log=(scale == 'log'))
    for name, (low, high) in control_space.items():
        p[name] = trial.suggest_float(name, float(low), float(high))
    return p


def make_objective(*, backend: dict = None, penalty_scale: float = PENALTY_SCALE,
                   control_space: dict = None, eval_opts=None):
    """objective(trial) -> effective_TAC を生成。

    eval_opts: callable(params) -> dict(apply_hi, hi_dT_min_K, apply_stage2)。
      None なら simulator.EVAL_KWARGS_DEFAULT 固定 (special.py と同じ apply_hi/stage2/dTmin=10)。
      後置ピンチ手法は params['hi_dT_min_K'] / apply_hi を読む eval_opts を渡す。
    """
    control_space = control_space or {}

    def _opts(params):
        if eval_opts is not None:
            return eval_opts(params)
        return dict(simulator.EVAL_KWARGS_DEFAULT)

    def objective(trial: optuna.trial.Trial) -> float:
        set_scale(penalty_scale)
        trial.set_user_attr('penalty_scale', penalty_scale)
        params = suggest_all(trial, control_space)
        design_params = {k: params[k] for k in space.PARAM_NAMES}
        design = space.build_design(design_params, backend)
        F = float(params['F_C3H8_fresh_kmol_h'])
        opts = _opts(params)

        t0 = time.perf_counter()
        result = simulator.raw_evaluate(design, F_fresh=F, **opts)
        trial.set_user_attr('wallclock_sec', time.perf_counter() - t0)
        _store_diagnostics(trial, result)
        trial.set_user_attr('F_C3H8_fresh_used_kmol_h', F)
        trial.set_user_attr('apply_hi', opts.get('apply_hi'))
        trial.set_user_attr('hi_dT_min_K', opts.get('hi_dT_min_K'))
        trial.set_user_attr('apply_stage2', opts.get('apply_stage2'))
        return result.effective_TAC

    return objective


def new_study(seed: int = 42) -> optuna.Study:
    """学習なしの study (RandomSampler)。探索点は enqueue で外から与える。"""
    sampler = make_sampler('random', seed, 0)
    return optuna.create_study(sampler=sampler, direction='minimize')


def run_batch(study: optuna.Study, objective, points: list, callback) -> list:
    """points (params dict のリスト) を enqueue して評価し、この batch で走った trial を返す。

    各 point は全 21 変数を含むこと (掃引変数=格子点、他=固定値)。callback は手法全体で
    1 個を使い回す (state が累積し progress/ETA が手法全体で連続する)。
    """
    before = len(study.trials)
    for p in points:
        study.enqueue_trial(p, skip_if_exists=False)
    run_optimization(study, objective, n_trials=len(points), show_progress_bar=False,
                     callbacks=[callback], n_jobs=1)
    return study.trials[before:]


def best_of(trials: list):
    """batch 内 best trial (feasible 優先、TAC 最小)。無ければ None。"""
    if not trials:
        return None
    feas = [t for t in trials if t.value is not None and t.user_attrs.get('is_feasible', False)]
    pool = feas or [t for t in trials if t.value is not None]
    if not pool:
        return None
    return min(pool, key=lambda t: t.value)


def save_table_csv(rows: list, path: str, fieldnames: list) -> None:
    """手法固有の表 (逐次のコスト曲線、dTmin 掃引表など) を CSV に。"""
    import csv
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fieldnames})


def finalize(study: optuna.Study, *, method: str, p_codes: str, description: str,
             settings: dict, extra_lines: list = None, eval_kwargs: dict = None,
             build_design=None, top_n: int = 3) -> tuple:
    """trials.csv / best.json / top-N 詳細 / README を 1 run dir に保存 (special.main 流)。(out_dir, best) を返す。"""
    out_dir = new_run_dir(method)
    complete, feasible, best = reporting.summarize(study)

    print("\n==== 結果 ====", flush=True)
    print(f"  完了評価: {len(complete)} / feasible: {len(feasible)}", flush=True)
    if best is not None:
        tag = "feasible best" if best.user_attrs.get('is_feasible') else "best (feasible 無し)"
        print(f"  {tag}: trial #{best.number}  effective_TAC={best.value:.2f} 億円/年", flush=True)
        try:
            _pur = float(best.user_attrs.get('c3h6_purity_wtfrac'))
            _prod = float(best.user_attrs.get('production_kmol_h'))
            _ff = float(best.params.get('F_C3H8_fresh_kmol_h'))
            print(f"    purity={_pur*100:.2f}wt% prod={_prod:.1f}kmol/h F_fresh={_ff:.1f} "
                  f"yield={_prod/_ff*100:.1f}%", flush=True)
        except (TypeError, ValueError, ZeroDivisionError):
            pass

    reporting.save_trials_csv(study, os.path.join(out_dir, 'trials.csv'))
    if best is not None:
        with open(os.path.join(out_dir, 'best.json'), 'w', encoding='utf-8') as f:
            json.dump({'method': method, 'p_codes': p_codes, 'number': best.number,
                       'effective_TAC': best.value, 'params': best.params,
                       'user_attrs': {k: v for k, v in best.user_attrs.items()}},
                      f, ensure_ascii=False, indent=2, default=str)
        print(f"  best JSON: {os.path.join(out_dir, 'best.json')}", flush=True)

    saved = reporting.save_best_reports(study, out_dir, top_n,
                                        eval_kwargs=eval_kwargs, build_design=build_design)
    reporting.write_readme(out_dir, method=method, p_codes=p_codes, description=description,
                           study=study, best=best, saved_reports=saved,
                           settings=settings, extra_lines=extra_lines)
    print(f"\n成果物: {os.path.abspath(out_dir)}/", flush=True)
    return out_dir, best
