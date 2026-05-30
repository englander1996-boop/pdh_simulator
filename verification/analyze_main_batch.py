r"""analyze_main_batch.py — run_main_batch.py が回した N 回分の main.py 結果をまとめて分析。

入力:
  verification/batch_<ts>/manifest.json  (run_main_batch.py が生成)
  各 run の生成物 outputs/main_<ts>/{best.json, trials.csv}  (manifest が場所を指す)

出力先の方針 (2026-05-29):
  各 run の生成物は outputs/ のまま。分析結果は verification/batch_<ts>/analysis/ に出す。

本バッチの実験設計:
  全 run が同一シード(SEED=42)。main.py は HYSYS(Dist2 COM)以外完全決定的なので、
  run 間の差は「HYSYS 非決定性が BO 軌道に与える揺らぎ」だけ。これを N=48 でサンプリングし、
  オプティマイザのロバスト性・局所解の分かれ方を測る。特に:
    - startup相 (trial #0〜#49 の Sobol QMC) は全 run でパラメータが完全一致する。
      その目的値ばらつき = 同一入力に対する HYSYS 再現性の直接測定 (ノイズ床)。
    - TPE相 (trial #50〜) の発散 = ノイズが探索軌道/局所解選択をどう変えるか。

出力 (verification/batch_<ts>/analysis/):
  report.md                 : 総合レポート (日本語、最初に開く)
  summary_per_run.csv       : run ごとの best 要約
  best_designs.csv          : N 個の best 設計の全21変数
  variable_stability.csv    : 変数ごとの収束安定性 (mean/std/CV)
  startup_noise.csv         : startup trial の run 間ばらつき (HYSYS ノイズ床)
  global_best.json          : 全 run 横断の最良設計 (best-of-N)
  *.png                     : 分布・収束帯・変数CV・局所解PCA・ノイズ床 の図

使い方:
  .\.venv\Scripts\python.exe verification\analyze_main_batch.py --batch verification\batch_<ts>
  .\.venv\Scripts\python.exe verification\analyze_main_batch.py --glob "outputs/main_*"   # manifest 無しで直接集約
"""

import sys
import json
import glob
import argparse
import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')   # 画面なしで PNG 保存
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

REPO = Path(__file__).resolve().parents[1]

# main.py SEARCH_SPACE の変数 (best 設計の変数分布解析に使う)。
# 反応器は REACTOR_KIND で軸流(z_cat/D_reactor)/径方向流(D_inner/bed_thickness/H)が
# 切り替わるため、両方のキーを列挙する (使われない側は CSV で空欄)。
PARAM_COLS = [
    'T_in_K', 't_cyc_min',
    'z_cat_m', 'D_reactor_m',                    # 軸流
    'D_inner_m', 'bed_thickness_m', 'H_m',       # 径方向流
    'D_psa_col_m', 'L_psa_bed_m', 'desorption_target',
    'P_H_Pa', 'A_mem_m2', 'F_C3H8_fresh_kmol_h',
    'col1_p_kpa', 'col1_n_stages', 'col1_feed_stage', 'col1_comp_frac_2',
    'col2_p_kpa', 'col2_n_stages', 'col2_feed_ratio', 'col2_reflux_ratio',
    'col3_p_kpa', 'col3_n_stages', 'col3_feed_ratio',
]
DEFAULT_N_STARTUP = 50   # main.py の N_STARTUP (Sobol QMC 相の長さ)


# ===========================================================================
# § 1. データ読み込み
# ===========================================================================
def _coerce_bool(series: pd.Series) -> pd.Series:
    """trials.csv の 'attr.is_feasible' は 'True'/'False' 文字列。bool へ変換。"""
    return series.astype(str).str.strip().str.lower().isin(('true', '1', 'yes'))


def _to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def _discover_run_dirs(args) -> tuple:
    """(run_subdir パスのリスト, manifest dict or None, batch_dir) を返す。

    run の生成物は outputs/ にある。manifest の subdir_abspath を最優先で使い、
    無ければ outputs/<subdir> にフォールバックする。
    """
    if args.batch:
        batch_dir = Path(args.batch)
        if not batch_dir.is_absolute():
            batch_dir = REPO / batch_dir
        mf = batch_dir / 'manifest.json'
        if mf.exists():
            with open(mf, encoding='utf-8') as f:
                manifest = json.load(f)
            dirs = []
            for r in manifest.get('runs', []):
                if not r.get('subdir'):
                    continue
                p = Path(r['subdir_abspath']) if r.get('subdir_abspath') else (REPO / 'outputs' / r['subdir'])
                if not p.exists():
                    p = REPO / 'outputs' / r['subdir']
                if p.exists():
                    dirs.append(p)
            return dirs, manifest, batch_dir
        # manifest 無し → outputs/ の main_* を拾う
        dirs = sorted(Path(p) for p in glob.glob(str(REPO / 'outputs' / 'main_*')) if Path(p).is_dir())
        return dirs, None, batch_dir
    # --glob 指定
    pattern = args.glob or str(REPO / 'outputs' / 'main_*')
    dirs = sorted(Path(p) for p in glob.glob(pattern) if Path(p).is_dir())
    out_dir = REPO / 'verification' / f'batch_glob_{datetime.datetime.now():%Y%m%d_%H%M%S}'
    return dirs, None, out_dir


def _load_run(subdir: Path) -> dict:
    """1 run の best.json と trials.csv を読み込み、解析用 dict を返す。"""
    rec = {'subdir': subdir.name, 'path': subdir, 'best': None, 'trials': None}
    bj = subdir / 'best.json'
    if bj.exists():
        try:
            with open(bj, encoding='utf-8') as f:
                rec['best'] = json.load(f)
        except Exception as e:
            print(f"  [警告] best.json 読込失敗 {subdir.name}: {e}", flush=True)
    tc = subdir / 'trials.csv'
    if tc.exists():
        try:
            df = pd.read_csv(tc)
            if 'value' in df.columns:
                df['value'] = df['value'].apply(_to_float)
            if 'attr.is_feasible' in df.columns:
                df['feasible'] = _coerce_bool(df['attr.is_feasible'])
            else:
                df['feasible'] = False
            rec['trials'] = df
        except Exception as e:
            print(f"  [警告] trials.csv 読込失敗 {subdir.name}: {e}", flush=True)
    return rec


# ===========================================================================
# § 2. 各種分析
# ===========================================================================
def _per_run_summary(runs: list) -> pd.DataFrame:
    """run ごとの best 要約 (best.json ベース、main.py の選定ロジックと一致)。"""
    rows = []
    for i, r in enumerate(runs):
        b = r['best'] or {}
        ua = b.get('user_attrs', {}) or {}
        params = b.get('params', {}) or {}
        df = r['trials']
        n_complete = int((df['state'] == 'COMPLETE').sum()) if df is not None and 'state' in df else np.nan
        n_feas = int(df['feasible'].sum()) if df is not None else np.nan
        prod = _to_float(ua.get('production_kmol_h'))
        ff = _to_float(ua.get('F_C3H8_fresh_used_kmol_h') or params.get('F_C3H8_fresh_kmol_h'))
        rows.append({
            'run_index': i,
            'subdir': r['subdir'],
            'best_trial': b.get('number'),
            'best_TAC': _to_float(b.get('effective_TAC')),
            'feasible': ua.get('is_feasible'),
            'purity_wt%': _to_float(ua.get('c3h6_purity_wtfrac')) * 100 if ua.get('c3h6_purity_wtfrac') else np.nan,
            'production_kmol_h': prod,
            'F_fresh_kmol_h': ff,
            'yield_%': (prod / ff * 100) if (prod and ff) else np.nan,
            'n_complete': n_complete,
            'n_feasible': n_feas,
        })
    return pd.DataFrame(rows)


def _best_designs_frame(runs: list) -> pd.DataFrame:
    """N 個の best 設計の全21変数を 1 行/run で。"""
    rows = []
    for i, r in enumerate(runs):
        b = r['best'] or {}
        params = b.get('params', {}) or {}
        row = {'run_index': i, 'best_TAC': _to_float(b.get('effective_TAC')),
               'best_trial': b.get('number')}
        for c in PARAM_COLS:
            row[c] = _to_float(params.get(c))
        rows.append(row)
    return pd.DataFrame(rows)


def _variable_stability(best_df: pd.DataFrame) -> pd.DataFrame:
    """変数ごとの収束安定性 (mean/std/CV/min/max)。CV 降順=シード(HYSYS)依存が強い順。"""
    rows = []
    for c in PARAM_COLS:
        v = best_df[c].dropna().values
        if len(v) == 0:
            continue
        mean = float(np.mean(v)); std = float(np.std(v, ddof=1)) if len(v) > 1 else 0.0
        cv = std / abs(mean) if mean != 0 else np.nan
        rows.append({'variable': c, 'mean': mean, 'std': std,
                     'CV': cv, 'min': float(np.min(v)), 'max': float(np.max(v)),
                     'n': len(v)})
    df = pd.DataFrame(rows).sort_values('CV', ascending=False, na_position='last')
    return df


def _startup_noise(runs: list, n_startup: int) -> tuple:
    """startup相 (trial # < n_startup) の run 間ばらつき = HYSYS ノイズ床。

    全 run で同一パラメータのはずの各 startup trial について、目的値(value)の
    run 間 std/range を集計。さらにパラメータが本当に一致しているか検証する。
    返り値: (startup_df, param_match_ok, max_param_reldiff, divergence_onset)
    """
    value_by_trial = {}     # trial番号 -> [各run の value]
    params_by_trial = {}    # trial番号 -> [各run の param ベクトル]
    for r in runs:
        df = r['trials']
        if df is None or 'number' not in df.columns:
            continue
        for _, row in df[df['number'] < n_startup].iterrows():
            n = int(row['number'])
            value_by_trial.setdefault(n, []).append(_to_float(row.get('value')))
            pv = [_to_float(row.get(c)) for c in PARAM_COLS if c in df.columns]
            params_by_trial.setdefault(n, []).append(pv)

    rows = []
    param_match_ok = True
    max_param_reldiff = 0.0
    for n in sorted(value_by_trial):
        vals = np.array([v for v in value_by_trial[n] if not np.isnan(v)])
        pmat = np.array(params_by_trial[n], dtype=float)
        if pmat.shape[0] > 1:
            col_range = np.nanmax(pmat, axis=0) - np.nanmin(pmat, axis=0)
            col_scale = np.maximum(np.abs(np.nanmean(pmat, axis=0)), 1e-9)
            reldiff = float(np.nanmax(col_range / col_scale))
            max_param_reldiff = max(max_param_reldiff, reldiff)
            if reldiff > 1e-6:
                param_match_ok = False
        if len(vals) >= 1:
            rows.append({
                'trial': n,
                'n_runs': len(vals),
                'value_mean': float(np.mean(vals)),
                'value_std': float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                'value_min': float(np.min(vals)),
                'value_max': float(np.max(vals)),
                'value_range': float(np.max(vals) - np.min(vals)),
            })
    startup_df = pd.DataFrame(rows)

    divergence_onset = None
    if not startup_df.empty:
        diverged = startup_df[startup_df['value_std'] > 1e-6]
        if not diverged.empty:
            divergence_onset = int(diverged['trial'].min())
    return startup_df, param_match_ok, max_param_reldiff, divergence_onset


def _convergence_curves(runs: list) -> dict:
    """run ごとの best-so-far(feasible TAC) 曲線。trial番号→累積最小。"""
    curves = {}
    for i, r in enumerate(runs):
        df = r['trials']
        if df is None or 'number' not in df.columns:
            continue
        d = df.sort_values('number')
        feas_val = d['value'].where(d['feasible'], other=np.inf).values
        best_so_far = np.minimum.accumulate(feas_val)
        best_so_far = np.where(np.isfinite(best_so_far), best_so_far, np.nan)
        curves[i] = (d['number'].values, best_so_far)
    return curves


def _failure_tally(runs: list) -> pd.Series:
    """全 run 合算の failure_unit 集計。"""
    counts = {}
    for r in runs:
        df = r['trials']
        if df is None or 'attr.failure_unit' not in df.columns:
            continue
        for fu in df['attr.failure_unit'].dropna().astype(str):
            fu = fu.strip()
            if fu and fu.lower() != 'nan':
                counts[fu] = counts.get(fu, 0) + 1
    return pd.Series(counts).sort_values(ascending=False)


# ===========================================================================
# § 3. プロット (ラベルは matplotlib 既定フォント互換のため英数字)
# ===========================================================================
def _plot_best_tac_hist(summary: pd.DataFrame, out: Path):
    feas = summary[summary['feasible'] == True]['best_TAC'].dropna()
    allv = summary['best_TAC'].dropna()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if len(allv):
        ax.hist(allv, bins=min(15, max(5, len(allv))), color='#cccccc',
                edgecolor='k', label=f'all runs (n={len(allv)})')
    if len(feas):
        ax.hist(feas, bins=min(15, max(5, len(feas))), color='#4c78a8',
                edgecolor='k', alpha=0.8, label=f'feasible (n={len(feas)})')
        ax.axvline(feas.min(), color='r', ls='--', lw=1.5,
                   label=f'best-of-N = {feas.min():.1f}')
    ax.set_xlabel('best effective_TAC [oku-yen/yr]')
    ax.set_ylabel('# runs')
    ax.set_title('Distribution of per-run best TAC (same seed, HYSYS noise)')
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def _plot_convergence(curves: dict, out: Path):
    fig, ax = plt.subplots(figsize=(8, 5))
    maxlen = max((len(x) for x, _ in curves.values()), default=0)
    grid = np.arange(maxlen)
    stacked = np.full((len(curves), maxlen), np.nan)
    for k, (i, (x, y)) in enumerate(curves.items()):
        ax.plot(x, y, color='0.7', lw=0.7, alpha=0.6)
        stacked[k, :len(y)] = y
    if maxlen:
        med = np.nanmedian(stacked, axis=0)
        q1 = np.nanpercentile(stacked, 25, axis=0)
        q3 = np.nanpercentile(stacked, 75, axis=0)
        ax.fill_between(grid, q1, q3, color='#4c78a8', alpha=0.25, label='IQR')
        ax.plot(grid, med, color='#1f4e79', lw=2.0, label='median')
    ax.set_xlabel('trial number')
    ax.set_ylabel('best-so-far feasible TAC [oku-yen/yr]')
    ax.set_title('Convergence bands across runs (same seed)')
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def _plot_variable_cv(stab: pd.DataFrame, out: Path):
    d = stab.dropna(subset=['CV']).copy()
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.barh(d['variable'][::-1], d['CV'][::-1], color='#e45756')
    ax.set_xlabel('CV (std/|mean|) of best design across runs')
    ax.set_title('Design-variable convergence stability\n(high CV = seed/HYSYS-sensitive lever)')
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def _plot_best_pca(best_df: pd.DataFrame, out: Path):
    """N 個の best 設計を標準化して PCA(2成分) 散布。局所解の basin を可視化。"""
    X = best_df[PARAM_COLS].values.astype(float)
    mask = ~np.isnan(X).any(axis=1)
    X = X[mask]
    tac = best_df['best_TAC'].values[mask]
    if X.shape[0] < 3:
        return
    mu = X.mean(axis=0); sd = X.std(axis=0); sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    pc = U[:, :2] * S[:2]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc = ax.scatter(pc[:, 0], pc[:, 1], c=tac, cmap='viridis', s=60, edgecolor='k')
    var = (S ** 2) / (S ** 2).sum()
    ax.set_xlabel(f'PC1 ({var[0]*100:.0f}% var)')
    ax.set_ylabel(f'PC2 ({var[1]*100:.0f}% var)')
    ax.set_title('Best designs in PCA space (clusters = local optima basins)')
    fig.colorbar(sc, label='best TAC [oku-yen/yr]')
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def _plot_startup_noise(startup_df: pd.DataFrame, out: Path):
    if startup_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(startup_df['trial'], startup_df['value_range'], color='#72b7b2')
    ax.set_xlabel('startup trial number (identical params across runs)')
    ax.set_ylabel('value range across runs [oku-yen/yr]')
    ax.set_title('HYSYS noise floor: objective spread on identical inputs')
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


# ===========================================================================
# § 4. レポート生成
# ===========================================================================
def _local_optima_summary(best_df: pd.DataFrame) -> dict:
    """best 設計の散らばりから局所解の分かれ方を要約。"""
    X = best_df[PARAM_COLS].values.astype(float)
    mask = ~np.isnan(X).any(axis=1)
    X = X[mask]
    out = {'n': int(X.shape[0])}
    if X.shape[0] >= 2:
        mu = X.mean(axis=0); sd = X.std(axis=0); sd[sd == 0] = 1.0
        Z = (X - mu) / sd
        d = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(axis=2))
        iu = np.triu_indices(len(Z), k=1)
        out['pairwise_mean'] = float(d[iu].mean())
        out['pairwise_max'] = float(d[iu].max())
    out['n_unique_best_trial'] = int(best_df['best_trial'].nunique(dropna=True))
    return out


def _write_report(out_dir: Path, batch_dir: Path, runs: list, summary: pd.DataFrame,
                  best_df: pd.DataFrame, stab: pd.DataFrame, startup_df: pd.DataFrame,
                  param_match_ok: bool, max_param_reldiff: float, divergence_onset,
                  fail_tally: pd.Series, global_best: dict, loc: dict, n_startup: int):
    L = []
    A = L.append
    A(f"# main.py {len(runs)}回バッチ 集計分析レポート")
    A("")
    A(f"- 生成日時: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    A(f"- batch: `{batch_dir.name}`  / 解析対象 run 数: **{len(runs)}**")
    A(f"- 各 run の生成物: `outputs/main_<ts>/`  (本分析は manifest 経由で参照)")
    A(f"- 実験設計: **同一シード(SEED=42)で main.py を {len(runs)} 回**。"
      f"HYSYS(Dist2)以外は完全決定的なので、run 間差は HYSYS 非決定性に起因。")
    A("")
    A("## 0. 結論サマリ")
    A("")
    feas = summary[summary['feasible'] == True]['best_TAC'].dropna()
    if len(feas):
        A(f"- **best-of-N (feasible)** = **{feas.min():.2f} 億円/年** "
          f"(run #{int(summary.loc[summary['best_TAC'].idxmin(),'run_index'])})")
        A(f"- per-run best TAC (feasible {len(feas)} run): "
          f"median {feas.median():.2f} / mean {feas.mean():.2f} / std {feas.std():.2f} / "
          f"IQR [{feas.quantile(.25):.2f}, {feas.quantile(.75):.2f}] / "
          f"range [{feas.min():.2f}, {feas.max():.2f}]")
        spread_pp = (feas.max() - feas.min()) / feas.median() * 100 if feas.median() else float('nan')
        A(f"- ばらつき幅 = best の **{spread_pp:.1f}%** (HYSYS ノイズが BO 最終解に与える振れ幅)")
    else:
        A("- feasible な best を持つ run が無い。")
    A(f"- feasible 到達 run: {int((summary['feasible']==True).sum())}/{len(summary)}")
    A("")
    A("## 1. HYSYS ノイズ床 (startup相の再現性)")
    A("")
    A(f"startup相 (trial #0〜#{n_startup-1} の Sobol QMC) は**全 run でパラメータが同一**。"
      f"その目的値の run 間ばらつきが、同一入力に対する HYSYS の非決定性そのもの。")
    A("")
    A(f"- パラメータ一致検証: **{'一致 (期待どおり)' if param_match_ok else '不一致あり'}** "
      f"(最大相対差 {max_param_reldiff:.2e})")
    if not startup_df.empty:
        A(f"- startup trial の目的値 range: 中央 {startup_df['value_range'].median():.3f} / "
          f"最大 {startup_df['value_range'].max():.3f} 億円/年")
        n_identical = int((startup_df['value_std'] <= 1e-6).sum())
        A(f"- 完全一致(std≦1e-6)した startup trial: {n_identical}/{len(startup_df)}")
    if divergence_onset is not None:
        A(f"- **発散開始 trial = #{divergence_onset}**: ここで初めて目的値が run 間でばらつく "
          f"(以降 TPE が異なる履歴を学習し軌道が分岐)。")
    else:
        A(f"- startup相では目的値の発散を検出せず (HYSYS が高度に再現的)。")
    A("")
    A("![startup noise](startup_noise.png)")
    A("")
    A("## 2. オプティマイザのロバスト性")
    A("")
    A("![best TAC hist](best_tac_hist.png)")
    A("")
    A("![convergence](convergence_bands.png)")
    A("")
    A("各 run の best-so-far(feasible)曲線と中央/IQR帯。帯が狭ければ HYSYS ノイズに対して "
      "BO はロバスト、広ければ最終解が運に左右される。")
    A("")
    A("## 3. 設計変数の収束安定性")
    A("")
    A("N 個の best 設計について変数ごとの分布。CV(=std/|mean|)が小さい変数は "
      "どの run でも同じ値に収束する『効くレバー』、大きい変数は HYSYS ノイズで動く『緩い変数』。")
    A("")
    A("| 変数 | mean | std | CV | min | max |")
    A("|---|---|---|---|---|---|")
    for _, row in stab.iterrows():
        A(f"| {row['variable']} | {row['mean']:.4g} | {row['std']:.4g} | "
          f"{row['CV']:.4f} | {row['min']:.4g} | {row['max']:.4g} |")
    A("")
    A("![variable CV](variable_cv.png)")
    A("")
    A("## 4. 局所解 (BO の有用性)")
    A("")
    A(f"- best 設計のユニーク trial番号数: {loc.get('n_unique_best_trial')}")
    if 'pairwise_mean' in loc:
        A(f"- 標準化空間でのペアワイズ距離: 平均 {loc['pairwise_mean']:.2f} / 最大 {loc['pairwise_max']:.2f} "
          f"(大きいほど別々の basin に着地=多峰的)")
    A("")
    A("![best designs PCA](best_designs_pca.png)")
    A("")
    A("PCA 散布で点が複数の塊に分かれていれば、BO は HYSYS ノイズ次第で異なる局所解 basin に "
      "着地している(=単一大域解への信頼性は限定的)。1 塊なら同一 basin にロバストに収束。")
    A("")
    A("## 5. 失敗ユニット集計 (全 run 合算)")
    A("")
    if len(fail_tally):
        A("| failure_unit | count |")
        A("|---|---|")
        for k, v in fail_tally.head(15).items():
            A(f"| {k} | {int(v)} |")
    else:
        A("- failure_unit の記録なし。")
    A("")
    A("## 6. best-of-N 設計 (全 run 横断の最良)")
    A("")
    if global_best:
        A(f"- run #{global_best.get('run_index')} / trial #{global_best.get('best_trial')} / "
          f"effective_TAC = **{global_best.get('best_TAC'):.2f} 億円/年**")
        ua = global_best.get('user_attrs', {}) or {}
        if ua.get('c3h6_purity_wtfrac'):
            A(f"- purity {float(ua['c3h6_purity_wtfrac'])*100:.2f} wt% / "
              f"production {_to_float(ua.get('production_kmol_h')):.1f} kmol/h")
        A(f"- 詳細レポート: `outputs/{global_best.get('subdir')}/top1_*.txt`")
        A("")
        A("```json")
        A(json.dumps(global_best.get('params', {}), ensure_ascii=False, indent=2))
        A("```")
    A("")
    A("## 付属ファイル")
    A("")
    A("- `summary_per_run.csv` … run ごとの best 要約")
    A("- `best_designs.csv` … N 個 best の全21変数")
    A("- `variable_stability.csv` … 変数ごと収束安定性")
    A("- `startup_noise.csv` … startup相の HYSYS ノイズ床")
    A("- `global_best.json` … best-of-N の params + 診断")
    A("")
    with open(out_dir / 'report.md', 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


# ===========================================================================
# § 5. main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description='main.py N回バッチの集計分析')
    ap.add_argument('--batch', type=str, default=None, help='verification/batch_<ts> ディレクトリ')
    ap.add_argument('--glob', type=str, default=None, help='run subdir の glob (manifest 無し時)')
    ap.add_argument('--n-startup', type=int, default=DEFAULT_N_STARTUP,
                    help=f'startup相の長さ (既定 {DEFAULT_N_STARTUP})')
    ap.add_argument('--out', type=str, default=None, help='出力先 (既定: batch_dir/analysis)')
    args = ap.parse_args()

    if not args.batch and not args.glob:
        ap.error('--batch か --glob のいずれかを指定してください。')

    dirs, manifest, batch_dir = _discover_run_dirs(args)
    if not dirs:
        print("  [エラー] 解析対象の run subdir が見つかりません。", flush=True)
        sys.exit(1)

    print(f"==== analyze_main_batch: {len(dirs)} run を集計 ====", flush=True)
    runs = []
    for d in dirs:
        r = _load_run(d)
        if r['best'] is None and r['trials'] is None:
            print(f"  [skip] 成果物なし: {d.name}", flush=True)
            continue
        runs.append(r)
    print(f"  有効 run: {len(runs)}", flush=True)
    if not runs:
        sys.exit(1)

    out_dir = Path(args.out) if args.out else (batch_dir / 'analysis')
    if not out_dir.is_absolute():
        out_dir = REPO / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- 分析 ----
    print("  [1/7] per-run サマリ", flush=True)
    summary = _per_run_summary(runs)
    summary.to_csv(out_dir / 'summary_per_run.csv', index=False, encoding='utf-8-sig')

    print("  [2/7] best 設計の変数表", flush=True)
    best_df = _best_designs_frame(runs)
    best_df.to_csv(out_dir / 'best_designs.csv', index=False, encoding='utf-8-sig')

    print("  [3/7] 変数収束安定性", flush=True)
    stab = _variable_stability(best_df)
    stab.to_csv(out_dir / 'variable_stability.csv', index=False, encoding='utf-8-sig')

    print("  [4/7] HYSYS ノイズ床 (startup相)", flush=True)
    startup_df, pmatch, max_reldiff, div_onset = _startup_noise(runs, args.n_startup)
    startup_df.to_csv(out_dir / 'startup_noise.csv', index=False, encoding='utf-8-sig')

    print("  [5/7] 収束曲線 / 失敗集計", flush=True)
    curves = _convergence_curves(runs)
    fail_tally = _failure_tally(runs)

    print("  [6/7] best-of-N / 局所解", flush=True)
    feas_rows = summary[summary['feasible'] == True].dropna(subset=['best_TAC'])
    if not feas_rows.empty:
        gi = int(feas_rows['best_TAC'].idxmin())
    else:
        gi = int(summary['best_TAC'].dropna().idxmin()) if summary['best_TAC'].notna().any() else 0
    grun = runs[int(summary.loc[gi, 'run_index'])]
    gb = grun['best'] or {}
    global_best = {
        'run_index': int(summary.loc[gi, 'run_index']),
        'subdir': grun['subdir'],
        'best_trial': gb.get('number'),
        'best_TAC': _to_float(gb.get('effective_TAC')),
        'params': gb.get('params', {}),
        'user_attrs': gb.get('user_attrs', {}),
    }
    with open(out_dir / 'global_best.json', 'w', encoding='utf-8') as f:
        json.dump(global_best, f, ensure_ascii=False, indent=2, default=str)
    loc = _local_optima_summary(best_df)

    print("  [7/7] 図の生成 (PNG)", flush=True)
    _plot_best_tac_hist(summary, out_dir / 'best_tac_hist.png')
    _plot_convergence(curves, out_dir / 'convergence_bands.png')
    _plot_variable_cv(stab, out_dir / 'variable_cv.png')
    _plot_best_pca(best_df, out_dir / 'best_designs_pca.png')
    _plot_startup_noise(startup_df, out_dir / 'startup_noise.png')

    _write_report(out_dir, batch_dir, runs, summary, best_df, stab, startup_df,
                  pmatch, max_reldiff, div_onset, fail_tally, global_best, loc, args.n_startup)

    print(f"\n==== 完了 ====", flush=True)
    print(f"  📌 レポート: {out_dir / 'report.md'}", flush=True)
    print(f"  ・CSV/PNG/global_best.json も同ディレクトリに保存", flush=True)


if __name__ == '__main__':
    main()
