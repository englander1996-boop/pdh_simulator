# -*- coding: utf-8 -*-
"""superbatch.manifest — 完了 run の記録 (resume) と best.json / params の読み出し。

manifest は 1 行 1 run の JSONL。1 run の rec:
  {idx, sampler, seed, returncode, out_dir(basename), effective_TAC, is_feasible, dur_sec, ts}
"""
import os
import csv
import glob
import json

from . import config


def run_dir(out_dir_basename):
    """run の成果物 dir (runs/main_<ts>) のフルパス。"""
    return os.path.join(config.RUNS_DIR, out_dir_basename) if out_dir_basename else None


def load_done():
    """manifest を読み、{(sampler, seed): rec} を返す (resume 用)。"""
    done = {}
    if os.path.exists(config.MANIFEST):
        for line in open(config.MANIFEST, encoding='utf-8'):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
                done[(r['sampler'], r['seed'])] = r
            except Exception:
                pass
    return done


def append(rec):
    """1 run の結果を manifest に追記する。"""
    os.makedirs(os.path.dirname(config.MANIFEST), exist_ok=True)
    with open(config.MANIFEST, 'a', encoding='utf-8') as f:
        f.write(json.dumps(rec, ensure_ascii=False) + '\n')


def snapshot_main_dirs():
    """現在存在する outputs/main_* dir の集合 (run 前後の差分で成果物を特定する)。"""
    return set(glob.glob(os.path.join(config.OUTPUTS, 'main_*')))


def read_best(out_dir_fullpath):
    """best.json (フルパス dir) → (effective_TAC, is_feasible)。無ければ (None, None)。"""
    bp = os.path.join(out_dir_fullpath, 'best.json')
    if not os.path.exists(bp):
        return None, None
    try:
        d = json.load(open(bp, encoding='utf-8'))
    except Exception:
        return None, None
    return d.get('effective_TAC'), bool(d.get('user_attrs', {}).get('is_feasible', False))


def read_params(out_dir_basename):
    """manifest の out_dir (basename) から best.json の params dict を読む (runs/ 配下)。無ければ {}。"""
    d = run_dir(out_dir_basename)
    if not d:
        return {}
    bp = os.path.join(d, 'best.json')
    if not os.path.exists(bp):
        return {}
    try:
        return json.load(open(bp, encoding='utf-8')).get('params', {}) or {}
    except Exception:
        return {}


def feasible_recs(recs):
    """feasible かつ TAC のある rec だけ返す。"""
    return [r for r in recs
            if r.get('is_feasible') and r.get('effective_TAC') is not None]


def write_summary_csv(recs, path=None):
    """1 run 1 行の一覧 CSV を書き出す (散乱防止の集約ファイル)。"""
    path = path or config.SUMMARY_CSV
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ['idx', 'sampler', 'seed', 'is_feasible', 'effective_TAC',
            'returncode', 'dur_sec', 'out_dir', 'ts']
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in sorted(recs, key=lambda x: x.get('idx', 0)):
            w.writerow([r.get(c, '') for c in cols])
