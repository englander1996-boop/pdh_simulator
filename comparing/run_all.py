# -*- coding: utf-8 -*-
r"""comparing/run_all.py — 全 comparing ケースを一括実行し BO との ΔTAC 比較を出力する単一エントリ。

これ 1 つを HYSYS のある PC で起動すれば、`comparing/case_*/` の全ケースを直列実行し、
各手法のベスト TAC を BO ベスト (outputs/main_*/best.json) と突き合わせて
  ΔTAC = TAC(本手法) − TAC(BO)   [億円/年, %]
を算出、比較表(CSV)・棒グラフ(PNG)・レポート貼付用 Markdown 表を生成する。

  実行:  .\.venv\Scripts\python.exe comparing\run_all.py
  下見:  .\.venv\Scripts\python.exe comparing\run_all.py --dry-run     (HYSYS を使わず計画だけ表示)
  一部:  .\.venv\Scripts\python.exe comparing\run_all.py --only case_rep_styrene2025,case_rep_eo2025
  基準:  .\.venv\Scripts\python.exe comparing\run_all.py --baseline outputs\main_20260601_150117\best.json

注意: 各ケースは Dist2=HYSYS で全系評価を数十回回すため 1 ケース ~1.5h、全 22 ケースで十数時間。
HYSYS は単一 COM インスタンスのため**直列実行のみ**(並列不可)。結果は 1 ケース完了ごとに
逐次 CSV へ追記するので、途中で止めても完了分は残る。
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import importlib.util
import json
import os
import sys
import time
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

_RESULTS_DIR = os.path.join(_HERE, 'results')

# レポート確定の BO ベスト (trial #194, outputs/main_20260601_150117)。outputs/ は
# git-ignore でラボ PC に無いため、レポートと同一基準で比較できるようリポジトリに同梱する。
# これを既定の baseline とし、--baseline で上書き可。同梱が無ければ最新 outputs/main_* に fallback。
_REPORT_BASELINE = os.path.join(_HERE, 'baseline_best.json')

# ---------------------------------------------------------------------------
# ケースのメタ情報 (匿名: テーマ/年・成熟度・手法。著者名は出さない)
#   成熟度 = 12 − 検出された欠陥(P)数。低いほど欠陥が多く BO との ΔTAC が大きい想定。
#   出典・対応は comparing/CASES.md / REPORT_METHODS_ANALYSIS.md を参照。
# ---------------------------------------------------------------------------
CASE_META = {
    # 欠陥部品 (手法そのものの単体再現)
    'case_p01_subsystem':   ('—',           '—',    '部分最適化 (ブロック別座標降下)'),
    'case_p02_pinch':       ('—',           '—',    '後置ピンチ (HIなし設計→後付け)'),
    'case_p04_sequential':  ('—',           '—',    '1次元逐次 (座標降下) + 整数の連続扱い'),
    'case_p05_grid':        ('—',           '—',    '粗グリッド総当たり'),
    'case_p06_multistart':  ('—',           '—',    '多始点 (大域性未保証)'),
    'case_p10_fug':         ('—',           '—',    'FUG 短絡 vs HYSYS 真値の乖離'),
    'case_p12_converge':    ('—',           '—',    '収束未検証 (1巡打切り)'),
    'case_combo_typical':   ('—',           '—',    '典型ワークフロー全体 (逐次→後置ピンチ)'),
    # 実在レポート手法の匿名再現 (テーマ/年, 成熟度, 手法)
    'case_rep_styrene2025':  ('スチレン/2025',     '2/12', '反応器ブロックのみ1次元 (部分最適化+後置ピンチ)'),
    'case_rep_butadiene2019':('ブタジエン/2019',    '2/12', '原料/空気比で反応器のみ1次元'),
    'case_rep_propro2019':   ('プロピオネート/2019','3/12', '反応条件 3×3 粗グリッド'),
    'case_rep_toluene2023':  ('トルエン/2023',     '3/12', '反応器 多変数フルグリッド'),
    'case_rep_pdh2025a':     ('PDH/2025',         '4/12', 'C3 塔を 圧力→段数→feed 逐次'),
    'case_rep_pdh2024a':     ('PDH/2024',         '4/12', '反応器 転化率を 形状/温度 で1次元'),
    'case_rep_pdh2025b':     ('PDH/2025',         '4/12', '1塔の段数のみ1次元'),
    'case_rep_dme2025':      ('DME/2025',         '4/12', '反応器 段数×温度×圧力 グリッド'),
    'case_rep_eo2019':       ('EO/2019',          '5/12', '4変数 座標降下'),
    'case_rep_methanol2025': ('メタノール/2025',   '5/12', '反応器→塔 の逐次'),
    'case_rep_eo2025':       ('EO/2025',          '7/12', '入口温度+塔段数 逐次 (熱統合内・検証あり)'),
    'case_rep_ammonia2025':  ('アンモニア/2025',   '—',    '(case README / docstring 参照)'),
    'case_rep_furfural2022': ('フルフラール/2022', '—',    '(case README / docstring 参照)'),
    'case_rep_mek2019':      ('MEK/2019',         '—',    '(case README / docstring 参照)'),
}


def discover_cases():
    """comparing/case_*/main.py を名前順に列挙して [(case_name, main_path), ...] を返す。"""
    out = []
    for path in sorted(glob.glob(os.path.join(_HERE, 'case_*', 'main.py'))):
        out.append((os.path.basename(os.path.dirname(path)), path))
    return out


def _load_module(case_name: str, main_path: str):
    """case の main.py をファイルパスから import (case ディレクトリは __init__.py を持たないため)。"""
    spec = importlib.util.spec_from_file_location(f'comparing_case_{case_name}', main_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _docstring_head(mod) -> str:
    doc = (mod.__doc__ or '').strip().splitlines()
    return doc[0].strip() if doc else ''


def load_baseline_tac(baseline_path: str | None) -> tuple[float, str]:
    """BO ベストの effective_TAC を返す。baseline 未指定なら最新の outputs/main_*/best.json。"""
    if baseline_path is None:
        if os.path.exists(_REPORT_BASELINE):
            baseline_path = _REPORT_BASELINE          # レポートと同一基準 (trial #194)
        else:
            cands = sorted(glob.glob(os.path.join(_REPO_ROOT, 'outputs', 'main_*', 'best.json')))
            if not cands:
                raise FileNotFoundError(
                    'BO ベースラインが見つかりません。--baseline で best.json を指定してください '
                    '(例 outputs/main_20260601_150117/best.json)')
            baseline_path = cands[-1]
    with open(baseline_path, encoding='utf-8') as f:
        data = json.load(f)
    tac = float(data.get('effective_TAC') or data.get('value'))
    return tac, baseline_path


def _newest_result_dir(before: set) -> str | None:
    """run() 実行前後で comparing/results/ に増えたディレクトリ(最新)を返す。"""
    if not os.path.isdir(_RESULTS_DIR):
        return None
    now = {os.path.join(_RESULTS_DIR, d) for d in os.listdir(_RESULTS_DIR)}
    new = [d for d in (now - before) if os.path.isdir(d)]
    if not new:
        return None
    return max(new, key=os.path.getmtime)


def _extract_tac(best, out_dir):
    """run() の戻り best (optuna trial) か out_dir/best.json から TAC・feasible・診断を取り出す。"""
    if best is not None and getattr(best, 'value', None) is not None:
        ua = getattr(best, 'user_attrs', {}) or {}
        return (float(best.value), ua.get('is_feasible'),
                ua.get('production_kmol_h'), ua.get('c3h6_purity_wtfrac'))
    if out_dir and os.path.exists(os.path.join(out_dir, 'best.json')):
        with open(os.path.join(out_dir, 'best.json'), encoding='utf-8') as f:
            d = json.load(f)
        ua = d.get('user_attrs', {}) or {}
        tac = d.get('effective_TAC') or d.get('value')
        return (float(tac) if tac is not None else None, ua.get('is_feasible'),
                ua.get('production_kmol_h'), ua.get('c3h6_purity_wtfrac'))
    return (None, None, None, None)


def run_case(case_name: str, main_path: str):
    """1 ケースを実行し結果 dict を返す。例外は捕捉して error フィールドに格納。"""
    before = set()
    if os.path.isdir(_RESULTS_DIR):
        before = {os.path.join(_RESULTS_DIR, d) for d in os.listdir(_RESULTS_DIR)}
    t0 = time.perf_counter()
    rec = {'case': case_name, 'TAC': None, 'feasible': None,
           'production_kmol_h': None, 'purity_wt': None, 'out_dir': None,
           'sec': None, 'error': ''}
    try:
        mod = _load_module(case_name, main_path)
        runner = getattr(mod, 'run', None) or getattr(mod, 'main', None)
        res = runner()
        out_dir, best = (res if isinstance(res, tuple) and len(res) == 2
                         else (_newest_result_dir(before), None))
        tac, feas, prod, pur = _extract_tac(best, out_dir)
        rec.update(TAC=tac, feasible=feas, production_kmol_h=prod,
                   purity_wt=pur, out_dir=out_dir)
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {e}'
        traceback.print_exc()
    rec['sec'] = round(time.perf_counter() - t0, 1)
    return rec


def _fmt(v, nd=2):
    return f'{v:.{nd}f}' if isinstance(v, (int, float)) else '—'


def write_outputs(rows, baseline_tac, baseline_path, out_csv, out_md, out_png):
    # ΔTAC 付与
    for r in rows:
        if isinstance(r['TAC'], (int, float)):
            r['dTAC'] = r['TAC'] - baseline_tac
            r['dTAC_pct'] = 100.0 * r['dTAC'] / baseline_tac
        else:
            r['dTAC'] = r['dTAC_pct'] = None
        meta = CASE_META.get(r['case'], ('—', '—', ''))
        r['theme'], r['maturity'], r['method'] = meta

    # CSV
    fields = ['case', 'theme', 'maturity', 'method', 'feasible', 'TAC',
              'dTAC', 'dTAC_pct', 'production_kmol_h', 'purity_wt', 'sec', 'error', 'out_dir']
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})

    # 並べ替え (ΔTAC 降順 = 損失大きい順。None は末尾)
    ranked = sorted(rows, key=lambda r: (r['dTAC'] is None, -(r['dTAC'] or 0)))

    # Markdown (レポート §11.3 貼付用)
    L = [f'# 他手法との ΔTAC 比較 (BO ベスト基準 TAC = {baseline_tac:.1f} 億円/年)',
         f'',
         f'基準 BO: `{os.path.relpath(baseline_path, _REPO_ROOT)}`  /  生成: run_all.py',
         f'ΔTAC = TAC(本手法) − TAC(BO)。正の値ほど BO に対する損失が大きい。',
         f'',
         '| ケース | テーマ/年 | 成熟度 | 手法 | feasible | TAC [億円/年] | ΔTAC [億円/年] | ΔTAC [%] |',
         '|---|---|---|---|:--:|--:|--:|--:|']
    for r in ranked:
        L.append('| {c} | {t} | {m} | {me} | {f} | {tac} | {d} | {p} |'.format(
            c=r['case'], t=r['theme'], m=r['maturity'], me=r['method'],
            f=('✓' if r['feasible'] else ('✗' if r['feasible'] is False else '—')),
            tac=_fmt(r['TAC'], 1), d=_fmt(r['dTAC'], 1), p=_fmt(r['dTAC_pct'], 1)))
    with open(out_md, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')

    # 棒グラフ (グレースケール・グリッドなし・四辺内向き目盛り = レポート図様式)
    plotted = [r for r in ranked if isinstance(r['dTAC'], (int, float))]
    if plotted:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib import font_manager
            for fnt in ['Yu Gothic', 'Meiryo', 'MS Gothic', 'Noto Sans CJK JP']:
                if any(fnt.lower() == e.name.lower() for e in font_manager.fontManager.ttflist):
                    matplotlib.rcParams['font.family'] = fnt
                    break
            matplotlib.rcParams.update({'axes.unicode_minus': False, 'axes.grid': False,
                                        'xtick.direction': 'in', 'ytick.direction': 'in',
                                        'xtick.top': True, 'ytick.right': True})
            labels = [f"{r['case'].replace('case_', '')}" for r in plotted]
            vals = [r['dTAC'] for r in plotted]
            fig, ax = plt.subplots(figsize=(8.5, max(4.0, 0.34 * len(plotted) + 1)))
            ax.barh(range(len(vals)), vals, color='0.45', edgecolor='black', linewidth=0.6)
            ax.set_yticks(range(len(vals)))
            ax.set_yticklabels(labels, fontsize=9)
            ax.invert_yaxis()
            ax.axvline(0, color='black', lw=0.8)
            ax.set_xlabel('ΔTAC = TAC(本手法) − TAC(BO)  [億円/年]', fontsize=11)
            ax.set_title('欠陥的最適化手法の BO に対する損失 ΔTAC', fontsize=12)
            fig.tight_layout()
            fig.savefig(out_png, dpi=200, bbox_inches='tight')
            fig.savefig(out_png.replace('.png', '.pdf'), bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            print(f'  (図の生成をスキップ: {e})', flush=True)

    return ranked


def main():
    ap = argparse.ArgumentParser(description='全 comparing ケースを実行し BO との ΔTAC を比較')
    ap.add_argument('--dry-run', action='store_true', help='HYSYS を使わず計画のみ表示 (import 検証込み)')
    ap.add_argument('--only', default='', help='実行するケースをカンマ区切りで限定')
    ap.add_argument('--baseline', default=None, help='BO best.json のパス (既定=最新 outputs/main_*/best.json)')
    args = ap.parse_args()

    cases = discover_cases()
    if args.only:
        want = {s.strip() for s in args.only.split(',') if s.strip()}
        cases = [(n, p) for n, p in cases if n in want]

    print('=' * 72, flush=True)
    print(f'  comparing/run_all — {len(cases)} ケース', flush=True)
    print('=' * 72, flush=True)

    if args.dry_run:
        try:
            tac, path = load_baseline_tac(args.baseline)
            print(f'  BO ベースライン: {os.path.relpath(path, _REPO_ROOT)}  (TAC={tac:.1f} 億円/年)')
        except Exception as e:
            print(f'  [警告] BO ベースライン未検出: {e}')
        print(f'  目安実行時間: 1 ケース ~1.5h × {len(cases)} = ~{len(cases) * 1.5:.0f}h (HYSYS 直列)')
        print('  --- ケース import 検証 ---')
        ok = 0
        for name, p in cases:
            try:
                mod = _load_module(name, p)
                has = 'run' if hasattr(mod, 'run') else ('main' if hasattr(mod, 'main') else '✗')
                meta = CASE_META.get(name, ('—', '—', _docstring_head(mod)))
                print(f'    [{has:>4}] {name:24s} {meta[0]:14s} 成熟度 {meta[1]:5s} {meta[2]}')
                ok += 1
            except Exception as e:
                print(f'    [ NG ] {name:24s} import 失敗: {type(e).__name__}: {e}')
        print(f'  import OK: {ok}/{len(cases)}。--dry-run 終了 (実走なし)。')
        return

    baseline_tac, baseline_path = load_baseline_tac(args.baseline)
    print(f'  BO ベースライン: {os.path.relpath(baseline_path, _REPO_ROOT)}  (TAC={baseline_tac:.1f} 億円/年)', flush=True)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out_csv = os.path.join(_RESULTS_DIR, f'comparison_{ts}.csv')
    out_md = os.path.join(_RESULTS_DIR, f'comparison_{ts}.md')
    out_png = os.path.join(_RESULTS_DIR, f'comparison_{ts}.png')

    rows = []
    t_start = time.perf_counter()
    for i, (name, p) in enumerate(cases, 1):
        print('\n' + '#' * 72, flush=True)
        print(f'# [{i}/{len(cases)}] {name} 実行中 ...', flush=True)
        print('#' * 72, flush=True)
        rec = run_case(name, p)
        rows.append(rec)
        # 1 ケース完了ごとに途中経過を保存 (中断耐性)
        write_outputs(rows, baseline_tac, baseline_path, out_csv, out_md, out_png)
        elapsed = time.perf_counter() - t_start
        eta = elapsed / i * (len(cases) - i)
        tac_s = _fmt(rec['TAC'], 1)
        d_s = _fmt(rec.get('dTAC'), 1)
        status = rec['error'] or f"TAC={tac_s}  ΔTAC={d_s}  feasible={rec['feasible']}"
        print(f'# 完了 [{i}/{len(cases)}] {name}: {status}  ({rec["sec"]}s)', flush=True)
        print(f'# 経過 {elapsed/3600:.1f}h / 残り目安 {eta/3600:.1f}h', flush=True)

    ranked = write_outputs(rows, baseline_tac, baseline_path, out_csv, out_md, out_png)

    print('\n' + '=' * 72, flush=True)
    print('  ΔTAC 比較 (損失大きい順)', flush=True)
    print('=' * 72, flush=True)
    print(f'  {"case":24s} {"成熟度":>5s} {"TAC":>9s} {"ΔTAC":>9s} {"ΔTAC%":>7s}  feas', flush=True)
    for r in ranked:
        print(f'  {r["case"]:24s} {r["maturity"]:>5s} {_fmt(r["TAC"],1):>9s} '
              f'{_fmt(r.get("dTAC"),1):>9s} {_fmt(r.get("dTAC_pct"),1):>7s}  '
              f'{("✓" if r["feasible"] else ("✗" if r["feasible"] is False else "—"))}', flush=True)
    print('\n  成果物:', flush=True)
    print(f'    CSV : {out_csv}', flush=True)
    print(f'    MD  : {out_md}   (レポート §11.3 貼付用)', flush=True)
    print(f'    図  : {out_png} / .pdf', flush=True)


if __name__ == '__main__':
    main()
