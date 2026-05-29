r"""analyze_seed_robustness.py — run_seed_robustness.py のバッチ(TPE seed散らし + random対照群)を分析。

入力:
  verification/seedrobust_<ts>/manifest.json  (run_seed_robustness.py が生成。meta.arms / meta.seeds / runs[] を持つ)
  各 run の生成物 outputs/main_<ts>/{best.json, trials.csv}

主役の2つの問い (前バッチ analyze_main_batch.py の「HYSYS ノイズ床」とは別物):
  (1) BO の正当性: ペア化シードで TPE が学習なし random を安定して上回るか。
      → win rate(TPE が勝ったシード割合)・改善幅・feasible 発見率の差。
  (2) 再現性/初期化ロバスト性: TPE が初期引きに依らず同じ best 設計に収束するか。
      → best TAC の散らばり・設計変数の CV(全シードで同値に張り付くレバー)・PCA(別 basin か)。

  注意: シードを散らしたので startup相(Sobol)の点列は run ごとに変わる。よって前バッチの
  「同一入力での HYSYS ノイズ床」は本バッチでは測れない(=この分析はそれを扱わない)。

出力 (verification/seedrobust_<ts>/analysis/):
  report.md               : 総合レポート (日本語、最初に開く)
  paired_comparison.csv   : シードごと TPE vs random の best TAC ペア表
  summary_<arm>.csv       : アームごとの per-run best 要約
  variable_stability_tpe.csv : TPE best 設計の変数収束安定性 (CV)
  best_designs_tpe.csv    : TPE best 設計の全21変数
  *.png                   : ペア散布・収束帯比較・best TAC 分布・変数CV・PCA

使い方:
  .\.venv\Scripts\python.exe verification\analyze_seed_robustness.py --batch verification\seedrobust_<ts>
"""

import sys
import json
import argparse
import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

# 同ディレクトリの analyze_main_batch.py のヘルパを再利用 (テスト済みの読込/集計/プロット)。
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_main_batch as amb   # noqa: E402

REPO = Path(__file__).resolve().parents[1]


# ===========================================================================
# § 1. データ読み込み (manifest からアーム別に run を展開)
# ===========================================================================
def _resolve_subdir(rec: dict) -> Path:
    """manifest の run レコードから生成物 subdir の実パスを得る。"""
    if rec.get('subdir_abspath'):
        p = Path(rec['subdir_abspath'])
        if p.exists():
            return p
    if rec.get('subdir'):
        return REPO / 'outputs' / rec['subdir']
    return Path('___missing___')


def _load_by_arm(batch_dir: Path) -> tuple:
    """(runs_by_arm, arms, seeds, manifest) を返す。runs は amb._load_run dict に seed/sampler を付与。"""
    mf = batch_dir / 'manifest.json'
    if not mf.exists():
        print(f"  [エラー] manifest が無い: {mf}", flush=True)
        sys.exit(1)
    with open(mf, encoding='utf-8') as f:
        manifest = json.load(f)
    meta = manifest.get('meta', {})
    arms = list(meta.get('arms', ['tpe', 'random']))
    seeds = list(meta.get('seeds', []))

    runs_by_arm = {a: [] for a in arms}
    for rec in manifest.get('runs', []):
        arm = rec.get('sampler')
        if arm not in runs_by_arm:
            continue
        subdir = _resolve_subdir(rec)
        if not subdir.exists():
            print(f"  [skip] 成果物なし: job#{rec.get('job')} {arm} seed={rec.get('seed')}", flush=True)
            continue
        run = amb._load_run(subdir)
        if run['best'] is None and run['trials'] is None:
            print(f"  [skip] best/trials なし: {subdir.name}", flush=True)
            continue
        run['seed'] = rec.get('seed')
        run['sampler'] = arm
        runs_by_arm[arm].append(run)
    return runs_by_arm, arms, seeds, manifest


def _feasible_best_tac(run: dict):
    """run の feasible な best TAC (best.json は feasible 優先選定なので is_feasible で判定)。"""
    b = run.get('best') or {}
    ua = b.get('user_attrs', {}) or {}
    if ua.get('is_feasible') is True:
        return amb._to_float(b.get('effective_TAC'))
    return None


def _global_best(runs_by_arm: dict) -> dict:
    """全 run(全アーム・全シード)横断の feasible 最良設計 = マルチスタート最適化の成果。

    本バッチは「ロバスト性検証」と同時に「BO で一番いい結果を探す」多スタート探索でもある。
    全 run の best.json から feasible かつ最小 TAC のものを選ぶ。
    """
    best = None
    for arm, runs in runs_by_arm.items():
        for r in runs:
            tac = _feasible_best_tac(r)
            if tac is None:
                continue
            if best is None or tac < best['best_TAC']:
                b = r['best'] or {}
                best = {
                    'sampler': arm,
                    'seed': r.get('seed'),
                    'subdir': r['subdir'],
                    'best_trial': b.get('number'),
                    'best_TAC': tac,
                    'params': b.get('params', {}),
                    'user_attrs': b.get('user_attrs', {}),
                }
    return best or {}


# ===========================================================================
# § 2. ペア化 BO vs random 比較
# ===========================================================================
def _paired_table(runs_by_arm: dict) -> pd.DataFrame:
    """シードごとに各アームの feasible best TAC を並べる。TPE と random の paired 比較用。"""
    seed_best = {}   # seed -> {arm: feasible best TAC or None}
    for arm, runs in runs_by_arm.items():
        for r in runs:
            seed_best.setdefault(r['seed'], {})[arm] = _feasible_best_tac(r)
    rows = []
    for seed in sorted(seed_best):
        d = seed_best[seed]
        tpe = d.get('tpe')
        rnd = d.get('random')
        row = {
            'seed': seed,
            'tpe_best_TAC': tpe,
            'random_best_TAC': rnd,
            'tpe_feasible': tpe is not None,
            'random_feasible': rnd is not None,
        }
        if tpe is not None and rnd is not None:
            row['diff_random_minus_tpe'] = rnd - tpe        # 正 = TPE の方が安い(良い)
            row['tpe_improve_%'] = (rnd - tpe) / rnd * 100 if rnd else np.nan
            row['tpe_wins'] = tpe < rnd
        else:
            row['diff_random_minus_tpe'] = np.nan
            row['tpe_improve_%'] = np.nan
            row['tpe_wins'] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _paired_stats(paired: pd.DataFrame) -> dict:
    """ペア表から head-to-head 統計を要約。"""
    both = paired.dropna(subset=['diff_random_minus_tpe'])
    out = {
        'n_seeds': int(len(paired)),
        'n_paired_feasible': int(len(both)),
        'tpe_feasible_rate': float(paired['tpe_feasible'].mean()) if len(paired) else np.nan,
        'random_feasible_rate': float(paired['random_feasible'].mean()) if len(paired) else np.nan,
    }
    if len(both):
        out['tpe_win_rate'] = float(both['tpe_wins'].mean())
        out['median_improve_%'] = float(both['tpe_improve_%'].median())
        out['mean_improve_%'] = float(both['tpe_improve_%'].mean())
        out['median_diff'] = float(both['diff_random_minus_tpe'].median())
    return out


def _plot_paired_scatter(paired: pd.DataFrame, out: Path):
    both = paired.dropna(subset=['diff_random_minus_tpe'])
    if both.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 6))
    x = both['random_best_TAC'].values
    y = both['tpe_best_TAC'].values
    ax.scatter(x, y, c='#4c78a8', s=55, edgecolor='k', zorder=3)
    lo = float(min(x.min(), y.min())); hi = float(max(x.max(), y.max()))
    pad = (hi - lo) * 0.05 or 1.0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'r--', lw=1.3, label='y = x (tie)')
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel('random best TAC [oku-yen/yr]')
    ax.set_ylabel('TPE best TAC [oku-yen/yr]')
    ax.set_title('Paired BO vs random (same seed)\nbelow diagonal = TPE better')
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def _plot_convergence_compare(curves_by_arm: dict, out: Path):
    colors = {'tpe': '#1f4e79', 'random': '#e45756'}
    fig, ax = plt.subplots(figsize=(8, 5))
    for arm, curves in curves_by_arm.items():
        if not curves:
            continue
        maxlen = max((len(x) for x, _ in curves.values()), default=0)
        grid = np.arange(maxlen)
        stacked = np.full((len(curves), maxlen), np.nan)
        for k, (_, (x, y)) in enumerate(curves.items()):
            stacked[k, :len(y)] = y
        if not maxlen:
            continue
        med = np.nanmedian(stacked, axis=0)
        q1 = np.nanpercentile(stacked, 25, axis=0)
        q3 = np.nanpercentile(stacked, 75, axis=0)
        c = colors.get(arm, '0.4')
        ax.fill_between(grid, q1, q3, color=c, alpha=0.18)
        ax.plot(grid, med, color=c, lw=2.0, label=f'{arm} median (n={len(curves)})')
    ax.set_xlabel('trial number')
    ax.set_ylabel('best-so-far feasible TAC [oku-yen/yr]')
    ax.set_title('Convergence bands: TPE vs random (across seeds)')
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


# ===========================================================================
# § 3. レポート
# ===========================================================================
def _write_report(out_dir, batch_dir, arms, seeds, runs_by_arm, paired, pstats,
                  per_arm_summary, stab_tpe, loc_tpe, global_best):
    L = []
    A = L.append
    A(f"# seed散らし + random対照群 バッチ 分析レポート")
    A("")
    A(f"- 生成日時: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    A(f"- batch: `{batch_dir.name}`")
    A(f"- 設計: **ペア化シード**({len(seeds)} シードを各アームに使用) / arms = {arms}")
    A(f"- 解析対象 run: " + " / ".join(f"{a}={len(runs_by_arm.get(a, []))}" for a in arms))
    A("")
    A("> 2つの問いを検証する: **(1) BO の正当性** = 同一シードで TPE が学習なし random を安定して"
      "上回るか / **(2) 再現性** = TPE が初期引きに依らず同じ best 設計へ収束するか。")
    A("")
    A("## 0. 結論サマリ")
    A("")
    A(f"- feasible 発見率: **TPE {pstats.get('tpe_feasible_rate', float('nan'))*100:.0f}%** "
      f"vs random {pstats.get('random_feasible_rate', float('nan'))*100:.0f}% "
      f"(seed {pstats['n_seeds']} 中)")
    if pstats.get('n_paired_feasible', 0) > 0:
        A(f"- **BO 勝率 (両者 feasible な {pstats['n_paired_feasible']} ペア)**: "
          f"TPE が random に勝ったシード = **{pstats['tpe_win_rate']*100:.0f}%**")
        A(f"- **改善幅**: TPE は random 比で中央 **{pstats['median_improve_%']:.2f}%** "
          f"(平均 {pstats['mean_improve_%']:.2f}%) 安い "
          f"= 中央 {pstats['median_diff']:.2f} 億円/年")
    else:
        A("- 両アームとも feasible なペアが無く、head-to-head 比較不可。")
    feas_tpe = [_feasible_best_tac(r) for r in runs_by_arm.get('tpe', [])]
    feas_tpe = [v for v in feas_tpe if v is not None]
    if feas_tpe:
        arr = np.array(feas_tpe)
        spread_pp = (arr.max() - arr.min()) / np.median(arr) * 100 if np.median(arr) else float('nan')
        A(f"- **TPE 収束 (再現性)**: best TAC = median {np.median(arr):.2f} / std {arr.std(ddof=1) if len(arr)>1 else 0:.2f} "
          f"/ range [{arr.min():.2f}, {arr.max():.2f}] → 散らばり幅 **{spread_pp:.1f}%**")
        A(f"  (この幅が小さいほど『どの初期引きでも同じ水準に収束』= 結果は引き運の産物でない)")
    if global_best:
        A(f"- **best-of-N (全 run 横断の最良設計)** = **{global_best['best_TAC']:.2f} 億円/年** "
          f"(arm={global_best['sampler']} / seed={global_best['seed']} / trial #{global_best['best_trial']})")
        gua = global_best.get('user_attrs', {}) or {}
        if gua.get('c3h6_purity_wtfrac'):
            A(f"  purity {float(gua['c3h6_purity_wtfrac'])*100:.2f} wt% / "
              f"production {amb._to_float(gua.get('production_kmol_h')):.1f} kmol/h "
              f"→ 詳細: `outputs/{global_best['subdir']}/top1_*.txt`")
    A("")
    A("## 1. BO vs random (ペア化 head-to-head) — BO の正当性")
    A("")
    A("![paired scatter](paired_scatter.png)")
    A("")
    A("同一シードでの TPE(縦) と random(横) の best TAC。対角線(y=x)より**下**にある点は TPE の勝ち。"
      "点が一貫して下側なら『BO の学習が単純探索を上回っている』直接証拠。")
    A("")
    A("![convergence compare](convergence_compare.png)")
    A("")
    A("best-so-far(feasible) の中央/IQR帯をアーム別に重ねた図。TPE の帯が random より低く"
      "(=安く)、かつ早く下がるほど BO の探索効率が高い。")
    A("")
    A("| seed | TPE best | random best | 差(rnd−tpe) | TPE改善% | TPE勝ち |")
    A("|---|---|---|---|---|---|")
    for _, r in paired.iterrows():
        tpe = f"{r['tpe_best_TAC']:.2f}" if pd.notna(r['tpe_best_TAC']) else '—(infeas)'
        rnd = f"{r['random_best_TAC']:.2f}" if pd.notna(r['random_best_TAC']) else '—(infeas)'
        diff = f"{r['diff_random_minus_tpe']:+.2f}" if pd.notna(r['diff_random_minus_tpe']) else '—'
        imp = f"{r['tpe_improve_%']:+.2f}" if pd.notna(r['tpe_improve_%']) else '—'
        win = {True: '✓', False: '✗'}.get(r['tpe_wins'], '—') if pd.notna(r['tpe_wins']) else '—'
        A(f"| {int(r['seed'])} | {tpe} | {rnd} | {diff} | {imp} | {win} |")
    A("")
    A("## 2. 収束性 / 再現性 (TPE アーム)")
    A("")
    A("![best TAC hist](best_tac_hist_tpe.png)")
    A("")
    A("TPE の各シード best TAC の分布。狭い1峰なら初期化ロバスト、広い/多峰なら最終解が引き運次第。")
    A("")
    A("### 設計変数の収束安定性 (CV = std/|mean|, 昇順=安定して同値に収束する効くレバー)")
    A("")
    A("| 変数 | mean | std | CV | min | max |")
    A("|---|---|---|---|---|---|")
    for _, row in stab_tpe.sort_values('CV', na_position='last').iterrows():
        A(f"| {row['variable']} | {row['mean']:.4g} | {row['std']:.4g} | "
          f"{row['CV']:.4f} | {row['min']:.4g} | {row['max']:.4g} |")
    A("")
    A("![variable CV](variable_cv_tpe.png)")
    A("")
    A("CV が小さい変数 = どのシードでも同じ値へ収束する『BO が確信して当てているレバー』。"
      "大きい変数 = 初期化/ノイズで動く緩い変数(目的値に効きにくいか、等価最適が複数)。")
    A("")
    A("![best designs PCA](best_designs_pca_tpe.png)")
    A("")
    A(f"- best 設計のユニーク basin 目安(標準化ペアワイズ距離): "
      f"平均 {loc_tpe.get('pairwise_mean', float('nan')):.2f} / 最大 {loc_tpe.get('pairwise_max', float('nan')):.2f}")
    A("PCA で1塊なら全シードが同一 basin にロバスト収束(理想)。複数塊なら初期化次第で別局所解。")
    A("")
    A("## 3. アーム別 per-run 要約")
    A("")
    for arm in arms:
        s = per_arm_summary.get(arm)
        if s is None or s.empty:
            continue
        feas = s[s['feasible'] == True]['best_TAC'].dropna()
        A(f"### {arm} (n={len(s)})")
        if len(feas):
            A(f"- feasible {len(feas)}/{len(s)} / best-of-N {feas.min():.2f} / "
              f"median {feas.median():.2f} / mean {feas.mean():.2f} / std {feas.std() if len(feas)>1 else 0:.2f}")
        else:
            A(f"- feasible 0/{len(s)}")
        A("")
    A("## 4. best-of-N 設計 (全 run 横断の最良) — BO で一番いい結果")
    A("")
    if global_best:
        A(f"- arm={global_best['sampler']} / seed={global_best['seed']} / trial #{global_best['best_trial']} / "
          f"effective_TAC = **{global_best['best_TAC']:.2f} 億円/年**")
        A(f"- 詳細レポート: `outputs/{global_best['subdir']}/top1_*.txt`  (CAPEX/OPEX/spec/HI 内訳)")
        A("")
        A("```json")
        A(json.dumps(global_best.get('params', {}), ensure_ascii=False, indent=2))
        A("```")
    else:
        A("- feasible な設計を持つ run が無い。")
    A("")
    A("## 付属ファイル")
    A("")
    A("- `global_best.json` … 全 run 横断の最良設計 (params + 診断)")
    A("- `paired_comparison.csv` … シードごと TPE vs random ペア表")
    A("- `summary_<arm>.csv` … アームごと per-run best 要約")
    A("- `best_designs_tpe.csv` / `variable_stability_tpe.csv` … TPE 収束分析")
    A("")
    with open(out_dir / 'report.md', 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


# ===========================================================================
# § 4. main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description='seed散らし+random対照群バッチの分析')
    ap.add_argument('--batch', type=str, required=True, help='verification/seedrobust_<ts> ディレクトリ')
    ap.add_argument('--out', type=str, default=None, help='出力先 (既定: batch_dir/analysis)')
    args = ap.parse_args()

    batch_dir = Path(args.batch)
    if not batch_dir.is_absolute():
        batch_dir = REPO / batch_dir

    print(f"==== analyze_seed_robustness: {batch_dir.name} ====", flush=True)
    runs_by_arm, arms, seeds, manifest = _load_by_arm(batch_dir)
    for a in arms:
        print(f"  {a}: {len(runs_by_arm.get(a, []))} run", flush=True)

    out_dir = Path(args.out) if args.out else (batch_dir / 'analysis')
    if not out_dir.is_absolute():
        out_dir = REPO / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- ペア化比較 ----
    print("  [1/5] ペア化 BO vs random", flush=True)
    paired = _paired_table(runs_by_arm)
    paired.to_csv(out_dir / 'paired_comparison.csv', index=False, encoding='utf-8-sig')
    pstats = _paired_stats(paired)

    # ---- アーム別 per-run 要約 ----
    print("  [2/5] アーム別 per-run 要約", flush=True)
    per_arm_summary = {}
    for arm in arms:
        runs = runs_by_arm.get(arm, [])
        if not runs:
            per_arm_summary[arm] = pd.DataFrame()
            continue
        s = amb._per_run_summary(runs)
        s.insert(1, 'seed', [r.get('seed') for r in runs])
        s.to_csv(out_dir / f'summary_{arm}.csv', index=False, encoding='utf-8-sig')
        per_arm_summary[arm] = s

    # ---- TPE 収束分析 ----
    print("  [3/5] TPE 収束分析 (変数CV / PCA)", flush=True)
    tpe_runs = runs_by_arm.get('tpe', [])
    if tpe_runs:
        best_df_tpe = amb._best_designs_frame(tpe_runs)
        best_df_tpe.to_csv(out_dir / 'best_designs_tpe.csv', index=False, encoding='utf-8-sig')
        stab_tpe = amb._variable_stability(best_df_tpe)
        stab_tpe.to_csv(out_dir / 'variable_stability_tpe.csv', index=False, encoding='utf-8-sig')
        loc_tpe = amb._local_optima_summary(best_df_tpe)
    else:
        best_df_tpe = pd.DataFrame(columns=amb.PARAM_COLS)
        stab_tpe = pd.DataFrame(columns=['variable', 'mean', 'std', 'CV', 'min', 'max'])
        loc_tpe = {}

    # ---- プロット ----
    print("  [4/5] 図の生成 (PNG)", flush=True)
    _plot_paired_scatter(paired, out_dir / 'paired_scatter.png')
    curves_by_arm = {a: amb._convergence_curves(runs_by_arm.get(a, [])) for a in arms}
    _plot_convergence_compare(curves_by_arm, out_dir / 'convergence_compare.png')
    if tpe_runs:
        amb._plot_best_tac_hist(per_arm_summary['tpe'], out_dir / 'best_tac_hist_tpe.png')
        amb._plot_variable_cv(stab_tpe, out_dir / 'variable_cv_tpe.png')
        amb._plot_best_pca(best_df_tpe, out_dir / 'best_designs_pca_tpe.png')

    # ---- best-of-N (マルチスタート最適化の成果) ----
    global_best = _global_best(runs_by_arm)
    with open(out_dir / 'global_best.json', 'w', encoding='utf-8') as f:
        json.dump(global_best, f, ensure_ascii=False, indent=2, default=str)

    # ---- レポート ----
    print("  [5/5] レポート生成", flush=True)
    _write_report(out_dir, batch_dir, arms, seeds, runs_by_arm, paired, pstats,
                  per_arm_summary, stab_tpe, loc_tpe, global_best)

    print(f"\n==== 完了 ====", flush=True)
    print(f"  📌 レポート: {out_dir / 'report.md'}", flush=True)
    if pstats.get('n_paired_feasible', 0):
        print(f"  BO 勝率 {pstats['tpe_win_rate']*100:.0f}% / 中央改善 {pstats['median_improve_%']:.2f}%", flush=True)
    if global_best:
        print(f"  best-of-N = {global_best['best_TAC']:.2f} 億円/年 "
              f"(arm={global_best['sampler']}, seed={global_best['seed']})", flush=True)


if __name__ == '__main__':
    main()
