# -*- coding: utf-8 -*-
"""superbatch.runner — main.py を 1 本ずつ直列起動する (HYSYS 競合ガード込み)。

⚠ HYSYS は単一インスタンス共有 ([[project_hysys_no_parallel]])。main.py を並列起動すると
   壊れるため、本モジュールは subprocess を直列に走らせる。手動 main.py との同時実行も不可。
"""
import os
import glob
import time
import shutil
import datetime
import statistics
import subprocess

from . import config, manifest
from .util import log, fmt_dur


def active_main_run():
    """直近 90s 以内に更新された trials_live.jsonl があれば、別の main.py 稼働中とみなしてパスを返す。

    HYSYS 競合防止のガード。super_main 起動時に呼ぶ (自分が回し始めた後は自分が唯一の起動元)。
    """
    now = time.time()
    for p in glob.glob(os.path.join(config.OUTPUTS, 'main_*', 'trials_live.jsonl')):
        try:
            if now - os.path.getmtime(p) < 90:
                return p
        except OSError:
            pass
    return None


def run_one(idx, total, sampler, seed, remaining, durations):
    """main.py を 1 回 subprocess 起動し、結果を manifest に追記して (rec, dur) を返す。

    env で PDH_SEED / PDH_SAMPLER / PDH_N_TRIALS を渡す。出力 (秒精度 ts の main_* dir) は
    起動前後の dir 差分で特定する (直列実行なので 1 run = 1 新 dir で衝突しない)。
    """
    env = dict(os.environ)
    env['PDH_SEED'] = str(seed)
    env['PDH_SAMPLER'] = sampler
    env['PDH_N_TRIALS'] = str(config.N_TRIALS_PER_RUN)
    env['PYTHONIOENCODING'] = 'utf-8'   # 子の 📌/日本語を utf-8 で logfile へ
    env['PYTHONUNBUFFERED'] = '1'

    logfile = os.path.join(config.LOG_DIR, f"run{idx:02d}_{sampler}_seed{seed}.log")
    med = statistics.median(durations) if durations else config.EST_RUN_SEC
    eta = remaining * med
    log(f"=== run {idx}/{total}  sampler={sampler} seed={seed}  "
        f"(残り {remaining} run, ETA ~{fmt_dur(eta)}) ===")
    log(f"    main.py 起動。ライブ BO ログ: {logfile}  (別窓で `tail -f` 可)")

    before = manifest.snapshot_main_dirs()
    t0 = time.time()
    with open(logfile, 'w', encoding='utf-8') as lf:
        proc = subprocess.run([config.PY, '-u', 'main.py'],
                              cwd=config.REPO, env=env,
                              stdout=lf, stderr=subprocess.STDOUT)
    dur = time.time() - t0
    rc = proc.returncode

    new_dirs = sorted(manifest.snapshot_main_dirs() - before)
    out_dir = new_dirs[-1] if new_dirs else None
    tac, feas = (manifest.read_best(out_dir) if out_dir else (None, None))

    # 成果物を outputs/super_main/runs/ へ移動 (outputs 直下に main_* を散乱させない)
    basename = None
    if out_dir:
        basename = os.path.basename(out_dir)
        try:
            os.makedirs(config.RUNS_DIR, exist_ok=True)
            dest = os.path.join(config.RUNS_DIR, basename)
            if os.path.abspath(out_dir) != os.path.abspath(dest):
                if os.path.exists(dest):
                    shutil.rmtree(dest, ignore_errors=True)
                shutil.move(out_dir, dest)
        except Exception as e:
            log(f"    ⚠ 成果物の runs/ 移動に失敗 ({type(e).__name__}: {e})。元の場所に残す。")

    rec = {
        'idx': idx, 'sampler': sampler, 'seed': seed, 'returncode': rc,
        'out_dir': os.path.basename(out_dir) if out_dir else None,
        'effective_TAC': tac, 'is_feasible': feas, 'dur_sec': round(dur, 1),
        'ts': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    manifest.append(rec)

    tac_s = f"{tac:.2f}" if tac is not None else "----"
    log(f"    完了 rc={rc}  TAC={tac_s}  feasible={feas}  所要 {fmt_dur(dur)}  "
        f"dir={rec['out_dir']}")
    if rc != 0:
        log(f"    ⚠ main.py が rc={rc} で異常終了。{logfile} を確認 (続行する)。")
    elif out_dir is None:
        log(f"    ⚠ 新しい main_* dir が見つからない。best.json 未生成かも。{logfile} を確認。")
    return rec, dur
