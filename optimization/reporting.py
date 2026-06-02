"""
optimization/reporting.py — 最適化結果の CSV / JSON / 比較表出力

pandas 不依存:
  - 標準ライブラリ csv / json で十分。numpy/scipy は環境にあるが、依存を最小化。
  - CSV は Excel / Optuna Dashboard / Jupyter の何処でも読める汎用形式。
  - JSON は best 1 件のみ (人間可読のサマリ)。
  - top-k 比較表は txt (BO vs 再評価の差分が一目でわかる形式)。
"""

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import List, Iterable, Optional
import optuna

from optimization.topk import TopKEntry


def save_trials_csv(
    study:     optuna.Study,
    path:      Path | str,
    param_keys: Optional[List[str]] = None,
) -> int:
    """全 trial を CSV に書き出す。

    Parameters
    ----------
    study : optuna.Study
    path : Path | str
        出力ファイルパス。
    param_keys : list[str] | None
        params 列の順序を固定したい場合に指定。None なら最初の COMPLETE trial の
        params キーから推定 (順序は dict 挿入順)。

    Returns
    -------
    int : 書き出した行数 (header 除く)。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    trials = study.trials
    if not trials:
        path.write_text('', encoding='utf-8')
        return 0

    if param_keys is None:
        param_keys = _infer_param_keys(trials)
    attr_keys = _infer_attr_keys(trials)

    fixed_cols = ['trial_number', 'state', 'value', 'datetime_start', 'datetime_complete']
    columns = fixed_cols + [f'param_{k}' for k in param_keys] + [f'attr_{k}' for k in attr_keys]

    n_rows = 0
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for t in trials:
            row = [
                t.number,
                t.state.name,
                t.value if t.value is not None else '',
                t.datetime_start.isoformat() if t.datetime_start else '',
                t.datetime_complete.isoformat() if t.datetime_complete else '',
            ]
            row += [t.params.get(k, '') for k in param_keys]
            row += [t.user_attrs.get(k, '') for k in attr_keys]
            writer.writerow(row)
            n_rows += 1
    return n_rows


def save_best_json(study: optuna.Study, path: Path | str) -> None:
    """ベスト trial を JSON に書き出す。

    完了済み trial が無い場合は空 JSON を出力。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        best = study.best_trial
        payload = {
            'trial_number':       best.number,
            'effective_TAC':      best.value,
            'params':             dict(best.params),
            'user_attrs':         dict(best.user_attrs),
            'datetime_complete':  best.datetime_complete.isoformat() if best.datetime_complete else None,
            'study_name':         study.study_name,
            'n_trials_total':     len(study.trials),
        }
    except ValueError:
        payload = {'study_name': study.study_name, 'n_trials_total': len(study.trials),
                   'note': '完了済み trial がありません'}
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')


def save_topk_report(
    entries:    List[TopKEntry],
    path:       Path | str,
    param_keys: Optional[List[str]] = None,
) -> None:
    """top-k 比較レポート (txt) を書き出す。

    BO ループの effective_TAC と、rigorous + Stage 2 再評価後の effective_TAC を
    並べて表示。差分 (Δ) が大きい候補は FUG bias の影響が大きかった事を示唆する。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not entries:
        path.write_text('top-k entries が空です。\n', encoding='utf-8')
        return

    if param_keys is None:
        param_keys = list(entries[0].params.keys())

    lines = []
    lines.append("=" * 100)
    lines.append("top-k 再評価レポート (BO 結果 vs rigorous + Stage 2 再評価)")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"{'rank':>5} {'trial':>6} {'TAC_bo':>12} {'TAC_re':>12} {'Δ':>10} "
                 f"{'feas_bo':>8} {'feas_re':>8}  failure_reason_re")
    lines.append("-" * 100)
    for e in entries:
        delta = e.effective_TAC_re - e.effective_TAC_bo
        lines.append(
            f"{e.rank:>5} {e.trial_number:>6} "
            f"{e.effective_TAC_bo:>12.3f} {e.effective_TAC_re:>12.3f} {delta:>+10.3f} "
            f"{str(e.is_feasible_bo):>8} {str(e.is_feasible_re):>8}  "
            f"{e.failure_reason_re}"
        )

    lines.append("")
    lines.append("=" * 100)
    lines.append("各候補の設計変数:")
    lines.append("=" * 100)
    for e in entries:
        lines.append("")
        lines.append(f"--- rank {e.rank} (trial #{e.trial_number}) ---")
        for k in param_keys:
            v = e.params.get(k)
            if isinstance(v, float):
                lines.append(f"  {k:<22} = {v:.6g}")
            else:
                lines.append(f"  {k:<22} = {v}")
        lines.append(f"  TAC_bo (BO)            = {e.effective_TAC_bo:.4f} 億円/年")
        lines.append(f"  TAC_re (再評価)         = {e.effective_TAC_re:.4f} 億円/年")
        if e.result.economics is not None:
            lines.append(f"  TAC (再評価, raw)       = {e.result.economics.TAC:.4f} 億円/年")
            lines.append(f"  Revenue (再評価, raw)   = {e.result.economics.total_revenue:.4f} 億円/年")
        if e.result.economics_hi is not None:
            lines.append(f"  TAC (HI 後)            = {e.result.economics_hi.TAC:.4f} 億円/年")
        if e.result.economics_synth is not None:
            lines.append(f"  TAC (Stage 2 後)        = {e.result.economics_synth.TAC:.4f} 億円/年")
            lines.append(f"  Profit (Stage 2 後)     = "
                         f"{e.result.economics_synth.total_revenue - e.result.economics_synth.TAC:.4f} 億円/年")
        if e.result.specs is not None:
            lines.append(f"  C3H6 純度               = {e.result.specs.c3h6_purity_wtfrac*100:.4f} wt%")
            lines.append(f"  H2 純度                 = {e.result.specs.h2_purity_molfrac*100:.4f} mol%")
            lines.append(f"  生産量                  = {e.result.specs.production_kmol_h:.2f} kmol/h "
                         f"(target: {e.result.specs.target_kmol_h:.2f})")
        if e.failure_reason_re:
            lines.append(f"  failure_reason          = {e.failure_reason_re}")

    path.write_text("\n".join(lines) + "\n", encoding='utf-8')


# ---------------------------------------------------------------------------
# 内部ヘルパ
# ---------------------------------------------------------------------------

def _infer_param_keys(trials: Iterable[optuna.trial.FrozenTrial]) -> List[str]:
    """完了済み trial から params のキー集合を抽出 (順序維持)。"""
    seen: dict[str, None] = {}
    for t in trials:
        if t.state == optuna.trial.TrialState.COMPLETE:
            for k in t.params.keys():
                if k not in seen:
                    seen[k] = None
    return list(seen.keys())


def _infer_attr_keys(trials: Iterable[optuna.trial.FrozenTrial]) -> List[str]:
    """全 trial から user_attrs のキー集合を抽出 (順序維持)。"""
    seen: dict[str, None] = {}
    for t in trials:
        for k in t.user_attrs.keys():
            if k not in seen:
                seen[k] = None
    return list(seen.keys())
