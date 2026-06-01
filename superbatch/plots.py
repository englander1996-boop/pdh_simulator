# -*- coding: utf-8 -*-
"""superbatch.plots — matplotlib による可視化 (PNG をバッチ dir の plots/ に出力)。

文字化け回避のためプロット内ラベルは英語。各 run の trials.csv は runs/ 配下を読む。

生成図:
  1. best_tac_by_sampler.png : sampler 別 best-TAC 分布 (箱ひげ+各 run 点)。BO が効くか一目で。
  2. convergence.png         : 各 run の feasible best-so-far vs trial 重ね描き。収束の一致を見る。
  3. best_tac_hist.png       : 全 feasible run の best-TAC ヒストグラム。
  4. param_convergence.png   : 主要変数で TPE best がどれだけ 1 点に集まるか (strip)。大域性の傍証。
"""
import os
import csv
import json

from . import config, manifest


def _read_trials_curve(out_dir_basename):
    """run の trials.csv (runs/ 配下) から feasible best-so-far の (trial, TAC) 系列。"""
    d = manifest.run_dir(out_dir_basename)
    if not d:
        return [], []
    path = os.path.join(d, 'trials.csv')
    if not os.path.exists(path):
        return [], []
    xs, ys, best = [], [], float('inf')
    for row in csv.DictReader(open(path, encoding='utf-8')):
        if row.get('attr.is_feasible') != 'True':
            continue
        try:
            n = int(row['number'])
            v = float(row['value'])
        except (ValueError, KeyError):
            continue
        best = min(best, v)
        xs.append(n)
        ys.append(best)
    return xs, ys


def _fig_best_tac_by_sampler(plt, recs):
    feas = manifest.feasible_recs(recs)
    data, positions, labels, pts = [], [], [], []
    for i, samp in enumerate(config.SAMPLERS):
        vals = [r['effective_TAC'] for r in feas if r['sampler'] == samp]
        if not vals:
            continue
        data.append(vals)
        positions.append(i)
        labels.append(f"{samp}\n(n={len(vals)})")
        pts.append((i, samp, vals))
    if not data:
        return None
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(data, positions=positions, widths=0.5, showfliers=False)
    for i, samp, vals in pts:
        jitter = [i + (((j % 7) - 3) * 0.03) for j in range(len(vals))]
        ax.scatter(jitter, vals, color=config.COLORS.get(samp, 'k'),
                   alpha=0.7, zorder=3, s=30)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel('best effective_TAC [oku-yen/yr]  (lower = better)')
    ax.set_title('best-TAC by sampler\n(TPE/CMA-ES below random => BO works)')
    ax.grid(axis='y', alpha=0.3)
    return fig


def _fig_convergence(plt, recs):
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = {}
    for r in recs:
        xs, ys = _read_trials_curve(r.get('out_dir'))
        if not xs:
            continue
        samp = r.get('sampler', '?')
        ax.plot(xs, ys, color=config.COLORS.get(samp, 'k'), alpha=0.5, lw=1.2)
        plotted[samp] = plotted.get(samp, 0) + 1
    if not plotted:
        return None
    handles = [plt.Line2D([0], [0], color=config.COLORS[s], lw=2,
                          label=f"{s} (n={plotted.get(s, 0)})")
               for s in config.SAMPLERS if plotted.get(s)]
    if handles:
        ax.legend(handles=handles)
    ax.set_xlabel('trial number')
    ax.set_ylabel('feasible best-so-far TAC [oku-yen/yr]')
    ax.set_title('convergence per run (overlaid)\nsame settling level => global hint')
    ax.grid(alpha=0.3)
    return fig


def _fig_hist(plt, recs):
    vals = [r['effective_TAC'] for r in manifest.feasible_recs(recs)]
    if not vals:
        return None
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(vals, bins=min(20, max(5, len(vals))), color='tab:blue', alpha=0.8,
            edgecolor='white')
    ax.axvline(min(vals), color='red', ls='--', lw=1.5, label=f'min={min(vals):.1f}')
    ax.set_xlabel('best effective_TAC [oku-yen/yr]')
    ax.set_ylabel('number of runs')
    ax.set_title(f'best-TAC histogram (all feasible runs, n={len(vals)})')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    return fig


def _fig_param_convergence(plt, recs):
    tpe_feas = [r for r in manifest.feasible_recs(recs) if r['sampler'] == 'tpe']
    sets = [(r, manifest.read_params(r['out_dir'])) for r in tpe_feas]
    sets = [(r, p) for r, p in sets if p]
    if len(sets) < 2:
        return None
    keys = [k for k in config.KEY_VARS
            if all(isinstance(p.get(k), (int, float)) for _, p in sets)]
    if not keys:
        return None
    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(2.2 * n, 4.2))
    if n == 1:
        axes = [axes]
    for ax, k in zip(axes, keys):
        vals = [float(p[k]) for _, p in sets]
        jitter = [((j % 5) - 2) * 0.04 for j in range(len(vals))]
        ax.scatter(jitter, vals, color='tab:blue', alpha=0.7, s=28)
        ax.set_xticks([])
        ax.set_title(k, fontsize=8)
        ax.grid(axis='y', alpha=0.3)
    fig.suptitle(f'TPE best params spread (n={len(sets)} runs)\n'
                 'tight vertical band => converged to same basin')
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def generate_all():
    """全図を生成して plots/ に保存する。matplotlib 不在やデータ不足は握りつぶす。"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  matplotlib 不在のため可視化スキップ ({type(e).__name__})。", flush=True)
        return

    recs = list(manifest.load_done().values())
    if not recs:
        print("  可視化: 完了 run が無いためスキップ。", flush=True)
        return
    os.makedirs(config.PLOTS_DIR, exist_ok=True)
    jobs = [
        ('best_tac_by_sampler.png', _fig_best_tac_by_sampler),
        ('convergence.png',         _fig_convergence),
        ('best_tac_hist.png',       _fig_hist),
        ('param_convergence.png',   _fig_param_convergence),
    ]
    for fname, fn in jobs:
        try:
            fig = fn(plt, recs)
        except Exception as e:
            print(f"  {fname}: スキップ ({type(e).__name__}: {e})", flush=True)
            continue
        if fig is None:
            print(f"  {fname}: データ不足でスキップ", flush=True)
            continue
        path = os.path.join(config.PLOTS_DIR, fname)
        fig.savefig(path, dpi=120, bbox_inches='tight')
        plt.close(fig)
        print(f"  保存: {path}", flush=True)
