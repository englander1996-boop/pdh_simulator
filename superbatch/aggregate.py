# -*- coding: utf-8 -*-
"""superbatch.aggregate — 全 run の集計。

画面出力に加えて、バッチ dir に summary.txt (集計テキスト) と summary.csv (1 run 1 行) を
書き出す。散乱防止のため出力はすべてバッチ親フォルダ配下にまとまる。
"""
import os
import statistics

from . import config, manifest


def _tac_stats(tacs):
    tacs = sorted(tacs)
    n = len(tacs)
    return {'n': n, 'min': tacs[0], 'med': statistics.median(tacs),
            'max': tacs[-1], 'std': (statistics.pstdev(tacs) if n > 1 else 0.0)}


def build_summary_lines(recs):
    """集計テキストを行リストで組み立てる (画面/ファイル共用)。"""
    L = []
    L.append("=" * 72)
    L.append("=== super_main 集計 ===")
    L.append("=" * 72)
    if not recs:
        L.append("  manifest が空。まだ完了した run がありません。")
        return L

    for samp in config.SAMPLERS:
        sub = [r for r in recs if r['sampler'] == samp]
        if not sub:
            continue
        feas = [r['effective_TAC'] for r in sub
                if r.get('is_feasible') and r.get('effective_TAC') is not None]
        if feas:
            st = _tac_stats(feas)
            L.append(f"\n[{samp}] feasible {st['n']}/{len(sub)} run  best-TAC: "
                     f"min={st['min']:.2f}  中央={st['med']:.2f}  max={st['max']:.2f}  "
                     f"std={st['std']:.2f}")
        else:
            L.append(f"\n[{samp}] feasible 0/{len(sub)} run")

    feas_all = manifest.feasible_recs(recs)

    if feas_all:
        champ = min(feas_all, key=lambda r: r['effective_TAC'])
        L.append("\n--- champion (★ 要 exp3 再検証、min はノイズで楽観バイアスあり) ---")
        L.append(f"  TAC={champ['effective_TAC']:.2f} 億円/年  "
                 f"sampler={champ['sampler']} seed={champ['seed']}")
        L.append(f"  best.json: runs/{champ['out_dir']}/best.json")
        L.append(f"  → 再現確認: 上記 params を exp/exp3.py に転記し PDH_HYSYS_FORCE_COLD=1 で "
                 f"3-5 回評価し TAC の mean±std を見る (単一 min を結果にしない)。")

    med_by = {}
    for samp in config.SAMPLERS:
        v = [r['effective_TAC'] for r in feas_all if r['sampler'] == samp]
        if v:
            med_by[samp] = statistics.median(v)
    if 'random' in med_by and ('tpe' in med_by or 'cmaes' in med_by):
        parts = [f"{s} 中央 {med_by[s]:.2f}" for s in config.SAMPLERS if s in med_by]
        L.append(f"\n  BO 正当性: {' / '.join(parts)}  "
                 f"(tpe/cmaes が random より低ければ最適化が効いている)")

    L += _param_convergence_lines(feas_all)
    L.append("=" * 72)
    return L


def _param_convergence_lines(feas_all):
    """TPE feasible の best params が 1 点に集まるか (CV) を行リストで返す。"""
    tpe_feas = [r for r in feas_all if r['sampler'] == 'tpe']
    param_sets = [p for p in (manifest.read_params(r['out_dir']) for r in tpe_feas) if p]
    if len(param_sets) < 3:
        return []
    keys = [k for k in param_sets[0]
            if all(isinstance(p.get(k), (int, float)) for p in param_sets)]
    rows = []
    for k in keys:
        vals = [float(p[k]) for p in param_sets]
        mean = statistics.mean(vals)
        sd = statistics.pstdev(vals)
        cv = (sd / abs(mean)) if mean else float('inf')
        rows.append((cv, k, mean, sd))
    L = [f"\n--- best params の収束一致 (TPE feasible n={len(param_sets)}、CV=std/|mean|) ---",
         "    CV 小=全 seed が同じ値に収束(大域の傍証) / CV 大=ばらつく(多峰 or 無関係)"]
    for cv, k, mean, sd in sorted(rows):
        flag = ' ←収束' if cv < 0.05 else (' ←ばらつき大' if cv > 0.25 else '')
        L.append(f"    {k:22s} mean={mean:12.3f}  std={sd:10.3f}  CV={cv:5.2f}{flag}")
    return L


def print_summary():
    """集計を画面表示し、summary.txt / summary.csv に書き出す。"""
    recs = list(manifest.load_done().values())
    lines = build_summary_lines(recs)
    print("\n" + "\n".join(lines), flush=True)
    try:
        if config.SUMMARY_TXT:
            os.makedirs(os.path.dirname(config.SUMMARY_TXT), exist_ok=True)
            with open(config.SUMMARY_TXT, 'w', encoding='utf-8') as f:
                f.write("\n".join(lines) + "\n")
        manifest.write_summary_csv(recs)
    except Exception as e:
        print(f"  (summary 書き出しに失敗: {type(e).__name__}: {e})", flush=True)
