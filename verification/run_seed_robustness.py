r"""run_seed_robustness.py — main.py を「シード散らし(TPE) + 学習なし対照群(random)」で回すバッチドライバ。

狙い (2026-05-29 ユーザー指示):
  前バッチ run_main_batch.py は SEED=42 固定で、startup相(Sobol QMC)の50点が全 run 同一
  → 実質「HYSYS ノイズだけ」を測っていた。今回はそれと逆に **BO の乱数も散らす**。
    (1) 再現性/初期化ロバスト性: どの初期引きから始めても同じ best 設計に収束するか。
        収束するなら「結果は1回の引き運の産物ではない」= BO の信頼性の正当性になる。
    (2) BO の正当性そのもの: 同じ予算・同じシードで「学習なし random 探索」も回し、
        TPE が安定して random を上回るなら「TPE の学習が効いている」直接証拠になる。

非対称ペア化設計 (本ドライバの核):
  TPE は全シードで回す(=best探しの multistart shot を最大化)。random は冒頭 K(=--control)
  シードだけで回す。random の K シードは TPE と同一シードなので paired 比較が成立する。
    → best探し(TPE 全シード)と BO 正当性検証(K ペアで BO vs random)を1バッチで両立。
  シードは OS エントロピー(secrets)から生成し manifest に記録 → 後で完全再現できる。

main.py への依存 (2026-05-29 追加の env 上書きを利用):
  main.py:79-86 で SEED/SAMPLER を env 上書き可能にした(既定 tpe/42、単体起動は不変)。
  本ドライバは run ごとに PDH_SEED(散らしたシード)・PDH_SAMPLER('tpe'|'random') を設定して起動する。
  それ以外の main.py のチューニングには一切触れない。

HYSYS 制約 (run_main_batch.py と同じ):
  HYSYS.Application は単一インスタンスの COM サーバー → 並列不可。各 run は別プロセスで
  main.py を素のまま逐次起動する(プロセス分離で HYSYS がクリーン起動・前 run 汚染なし・
  1 run クラッシュしても次は無傷)。

成果物 (verification/seedrobust_<ts>/):
  manifest.json   : meta(seeds[], arms) + runs[]{job, seed, sampler, subdir, exit_code, wall_sec, best_TAC, feasible}
                    各 run 完了ごとに上書き保存 = 中断耐性・--resume 判定に使用。
  run_NN.log      : 各 run の main.py 全出力。
  batch.log       : バッチ全体の進捗ログ。
  各 run の生成物 (main_<ts>/) は outputs/ にそのまま (main.py 既定、移動しない)。

使い方 (Python は .venv の 3.13):
  .\.venv\Scripts\python.exe verification\run_seed_robustness.py                  # TPE24 + random8 = 32 run (既定)
  .\.venv\Scripts\python.exe verification\run_seed_robustness.py --seeds 24 --control 8
  .\.venv\Scripts\python.exe verification\run_seed_robustness.py --arms tpe        # 検証なし、best探しのみ(全TPE)
  .\.venv\Scripts\python.exe verification\run_seed_robustness.py --resume verification\seedrobust_<ts>
  .\.venv\Scripts\python.exe verification\run_seed_robustness.py --dry-run

逐次実行。1 run ~30-60分。TPE24+random8=32 run で約 16-32 時間が目安。
別ターミナルで verification\seedrobust_<ts>\run_NN.log を tail すると個別 run の trial 進捗が live で見える。
"""

import os
import sys
import json
import time
import secrets
import argparse
import datetime
import subprocess
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')
except Exception:
    pass


# ===========================================================================
# § 1. 設定
# ===========================================================================
N_SEEDS_DEFAULT   = 24           # TPE を回すシード数 (=best探しの multistart shot 数)。
N_CONTROL_DEFAULT = 8            # うち random 対照群も回すシード数 (paired 比較用の冒頭 K シード)。
ARMS_DEFAULT      = ('tpe', 'random')   # TPE(本命) + 学習なし対照群(random)
SEED_MIN, SEED_MAX = 1, 2**31 - 2     # 生成シードの範囲 (Optuna/numpy が安全に扱う正の int32)

REPO    = Path(__file__).resolve().parents[1]
OUTPUTS = REPO / 'outputs'
VERIFY  = REPO / 'verification'
MAIN_PY = REPO / 'main.py'
VENV_PY = REPO / '.venv' / 'Scripts' / 'python.exe'

_ECHO_MARKERS = ('★', '====', '成果物', 'Traceback', 'Error', 'ERROR', '結果')
_ETA_ECHO_EVERY = 20


# ===========================================================================
# § 2. ユーティリティ
# ===========================================================================
def _now_ts() -> str:
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')


def _fmt_dur(sec: float) -> str:
    sec = int(sec)
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}:{m:02d}:{s:02d}"
    if m:
        return f"{m}:{s:02d}"
    return f"{s}s"


def _resolve_python() -> str:
    if VENV_PY.exists():
        return str(VENV_PY)
    print(f"  [警告] .venv python が無い ({VENV_PY})。sys.executable で代替する。", flush=True)
    return sys.executable


def _snapshot_main_dirs() -> set:
    return {p.name for p in OUTPUTS.glob('main_*') if p.is_dir()}


def _read_best_json(subdir: Path) -> dict:
    bj = subdir / 'best.json'
    if not bj.exists():
        return {}
    try:
        with open(bj, encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return {}


class _Tee:
    """バッチ進捗を batch.log とコンソールの両方へ flush 付きで書く。"""
    def __init__(self, log_path: Path):
        self._f = open(log_path, 'a', encoding='utf-8')

    def __call__(self, msg: str = '') -> None:
        print(msg, flush=True)
        self._f.write(msg + '\n')
        self._f.flush()

    def close(self) -> None:
        try:
            self._f.close()
        except Exception:
            pass


def _build_jobs(seeds: list, arms: tuple, n_control: int) -> list:
    """非対称ジョブ列を作る。TPE は全シード、random は冒頭 n_control シードだけ。

    狙いの両立:
      - best探し: TPE を全シードで回す(multistart shot 最大化)。
      - BO の正当性検証: 冒頭 n_control シードは random も回す(同一シード=paired 比較成立)。
    並び順は「ペア(tpe→random)を先頭に、残りを TPE 単独」。途中で止めても paired 比較が揃う。
    """
    use_tpe = 'tpe' in arms
    use_rnd = 'random' in arms
    jobs = []
    for i, s in enumerate(seeds):
        if use_tpe:
            jobs.append({'seed': int(s), 'sampler': 'tpe'})
        if use_rnd and i < n_control:
            jobs.append({'seed': int(s), 'sampler': 'random'})
    return jobs


# ===========================================================================
# § 3. 1 run の実行
# ===========================================================================
def _run_one(job_index: int, n_jobs: int, seed: int, sampler: str,
             batch_dir: Path, pylog) -> dict:
    """main.py を 1 回サブプロセスで起動 (PDH_SEED/PDH_SAMPLER を env 注入)。生成 subdir を特定して記録。"""
    py = _resolve_python()
    run_log_path = batch_dir / f'run_{job_index:02d}.log'

    before = _snapshot_main_dirs()
    pylog(f"\n{'='*72}")
    pylog(f"  ▶ run #{job_index:02d}/{n_jobs:02d}  sampler={sampler}  seed={seed}  "
          f"start {datetime.datetime.now():%H:%M:%S}")
    pylog(f"    cmd: {py} main.py  (PDH_SAMPLER={sampler}, PDH_SEED={seed}, cwd={REPO})")
    pylog(f"    log: {run_log_path}")
    pylog(f"{'='*72}")

    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    env['PDH_SEED'] = str(seed)          # main.py:79 がこれを読む(seed散らし)
    env['PDH_SAMPLER'] = sampler         # main.py:80 がこれを読む('tpe'|'random')

    t0 = time.perf_counter()
    exit_code = None
    eta_count = 0
    with open(run_log_path, 'w', encoding='utf-8') as logf:
        try:
            proc = subprocess.Popen(
                [py, 'main.py'],
                cwd=str(REPO),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=env,
                encoding='utf-8',
                errors='replace',
                bufsize=1,
            )
            for line in proc.stdout:
                logf.write(line)
                logf.flush()
                s = line.rstrip('\n')
                if 'ETA' in s:
                    eta_count += 1
                    if eta_count % _ETA_ECHO_EVERY == 0:
                        print(f"    [#{job_index:02d} {sampler}] {s.strip()}", flush=True)
                elif any(m in s for m in _ECHO_MARKERS):
                    print(f"    [#{job_index:02d} {sampler}] {s.strip()}", flush=True)
            proc.wait()
            exit_code = proc.returncode
        except Exception as e:
            logf.write(f"\n[run_seed_robustness] サブプロセス起動失敗: {type(e).__name__}: {e}\n")
            pylog(f"    ✗ 起動失敗: {type(e).__name__}: {e}")
            exit_code = -999
    wall = time.perf_counter() - t0

    after = _snapshot_main_dirs()
    new_dirs = sorted(after - before)
    subdir_name = new_dirs[-1] if new_dirs else None
    if len(new_dirs) > 1:
        pylog(f"    [警告] 新規 main_* が複数検出: {new_dirs} → 最新を採用: {subdir_name}")

    best_tac = feasible = best_num = None
    if subdir_name:
        bj = _read_best_json(OUTPUTS / subdir_name)
        if bj:
            best_tac = bj.get('effective_TAC')
            best_num = bj.get('number')
            feasible = (bj.get('user_attrs', {}) or {}).get('is_feasible')

    rec = {
        'job': job_index,
        'seed': int(seed),
        'sampler': sampler,
        'subdir': subdir_name,
        'subdir_abspath': str(OUTPUTS / subdir_name) if subdir_name else None,
        'exit_code': exit_code,
        'wall_sec': round(wall, 1),
        'best_trial_number': best_num,
        'best_TAC': best_tac,
        'feasible': feasible,
        'finished_at': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    tac_s = f"{best_tac:.2f}" if isinstance(best_tac, (int, float)) else '----'
    feas_s = {True: 'feasible✓', False: 'infeasible✗'}.get(feasible, '?')
    pylog(f"  ◀ run #{job_index:02d} 完了  exit={exit_code}  {_fmt_dur(wall)}  "
          f"[{sampler} seed={seed}]  best_TAC={tac_s}({feas_s})  subdir=outputs/{subdir_name}")
    return rec


# ===========================================================================
# § 4. manifest 保存 / 読込
# ===========================================================================
def _save_manifest(batch_dir: Path, records: list, meta: dict) -> None:
    with open(batch_dir / 'manifest.json', 'w', encoding='utf-8') as f:
        json.dump({'meta': meta, 'runs': records}, f, ensure_ascii=False, indent=2)


def _load_manifest(batch_dir: Path) -> dict:
    mf = batch_dir / 'manifest.json'
    if not mf.exists():
        return {}
    with open(mf, encoding='utf-8') as f:
        return json.load(f)


# ===========================================================================
# § 5. main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description='main.py を seed散らし(TPE)+対照群(random) で回すバッチ')
    ap.add_argument('--seeds', type=int, default=N_SEEDS_DEFAULT,
                    help=f'TPE を回すシード数(=best探しの shot 数) (既定 {N_SEEDS_DEFAULT})')
    ap.add_argument('--control', type=int, default=N_CONTROL_DEFAULT,
                    help=f'うち random 対照群も回す冒頭シード数(paired 比較用) (既定 {N_CONTROL_DEFAULT})')
    ap.add_argument('--arms', type=str, default=','.join(ARMS_DEFAULT),
                    help=f"アーム (カンマ区切り)。既定 '{','.join(ARMS_DEFAULT)}'。'tpe' 単独も可")
    ap.add_argument('--resume', type=str, default=None,
                    help='既存 seedrobust_<ts> を指定して未完了分を再開 (記録済みシードを再利用)')
    ap.add_argument('--dry-run', action='store_true', help='ジョブ列を表示するだけで実行しない')
    args = ap.parse_args()

    arms = tuple(a.strip() for a in args.arms.split(',') if a.strip())
    for a in arms:
        if a not in ('tpe', 'random'):
            ap.error(f"未知の arm: {a!r} (許容: 'tpe' | 'random')")

    OUTPUTS.mkdir(exist_ok=True)
    VERIFY.mkdir(exist_ok=True)

    # ---- batch dir + シード決定 (新規 or 再開) ----
    if args.resume:
        batch_dir = Path(args.resume)
        if not batch_dir.is_absolute():
            batch_dir = REPO / batch_dir
        if not batch_dir.exists():
            print(f"  [エラー] 再開先が存在しない: {batch_dir}", flush=True)
            sys.exit(1)
        manifest = _load_manifest(batch_dir)
        meta = manifest.get('meta', {})
        records = manifest.get('runs', [])
        seeds = meta.get('seeds', [])
        arms = tuple(meta.get('arms', arms))         # 再開時は記録済み構成を尊重
        n_control = int(meta.get('n_control', len(seeds)))
        jobs = _build_jobs(seeds, arms, n_control)
        start_index = len(records)
        print(f"  ▶ 再開: {batch_dir}  (記録済み {len(records)}/{len(jobs)} run)", flush=True)
    else:
        ts = _now_ts()
        batch_dir = VERIFY / f'seedrobust_{ts}'
        batch_dir.mkdir(parents=True, exist_ok=True)
        # OS エントロピーから重複なしシードを生成 (再現用に manifest へ保存)。
        n_control = max(0, min(args.control, args.seeds))   # random は TPE シードの部分集合
        rng = secrets.SystemRandom()
        seeds = []
        while len(seeds) < args.seeds:
            s = rng.randint(SEED_MIN, SEED_MAX)
            if s not in seeds:
                seeds.append(s)
        jobs = _build_jobs(seeds, arms, n_control)
        records = []
        start_index = 0
        meta = {
            'created_at': datetime.datetime.now().isoformat(timespec='seconds'),
            'experiment': 'seed_robustness + random control (asymmetric: TPE all seeds, random first K)',
            'n_seeds': args.seeds,
            'n_control': n_control,
            'arms': list(arms),
            'seeds': seeds,
            'n_runs_total': len(jobs),
            'note': ('TPE は全シードで multistart(best探し)、random は冒頭 K シードだけで paired 比較。'
                     'TPE が learning なし random を安定して上回るか(BOの正当性)・初期引きに依らず'
                     '同じ解に収束するか(再現性)・全 run 横断 best-of-N(最良設計) を検証。'),
            'run_outputs_location': 'outputs/main_<ts>/',
            'main_py': str(MAIN_PY),
            'python': _resolve_python(),
        }

    if args.dry_run:
        print(f"  [dry-run] batch_dir = {batch_dir}")
        print(f"  [dry-run] arms = {arms} / TPE seeds = {len(seeds)} / random(control) seeds = {n_control}")
        print(f"  [dry-run] seeds = {seeds}")
        print(f"  [dry-run] 総 run 数 = {len(jobs)} (start_index={start_index})")
        for j, job in enumerate(jobs):
            mark = 'done' if j < start_index else ' '
            print(f"    [{mark}] run #{j:02d}: PDH_SAMPLER={job['sampler']}  PDH_SEED={job['seed']}")
        return

    pylog = _Tee(batch_dir / 'batch.log')
    pylog(f"==== run_seed_robustness: {len(jobs)} run "
          f"(TPE {len(seeds)} seeds + random {meta.get('n_control', 0)} control seeds) ====")
    pylog(f"  batch_dir(管理/分析) = {batch_dir}")
    pylog(f"  run 生成物            = {OUTPUTS}/main_<ts>/")
    pylog(f"  python                = {meta['python']}")
    pylog(f"  シード                = {seeds}")
    pylog(f"  開始 index = {start_index} / 総 {len(jobs)}")
    pylog(f"  ※逐次実行。1 run ~30-60分。中断時は --resume {batch_dir} で再開可(シードは記録済みを再利用)。")

    _save_manifest(batch_dir, records, meta)

    batch_t0 = time.perf_counter()
    durations = [r['wall_sec'] for r in records if isinstance(r.get('wall_sec'), (int, float))]
    for i in range(start_index, len(jobs)):
        job = jobs[i]
        rec = _run_one(i, len(jobs), job['seed'], job['sampler'], batch_dir, pylog)
        records.append(rec)
        _save_manifest(batch_dir, records, meta)
        if isinstance(rec.get('wall_sec'), (int, float)):
            durations.append(rec['wall_sec'])

        # ---- 進捗 + ETA + アーム別 best-of-batch ----
        n_remain = len(jobs) - (i + 1)
        med = sorted(durations)[len(durations) // 2] if durations else 0.0
        eta = n_remain * med
        elapsed = time.perf_counter() - batch_t0
        bob = {}
        for arm in arms:
            tacs = [r['best_TAC'] for r in records
                    if r.get('sampler') == arm and r.get('feasible') is True
                    and isinstance(r.get('best_TAC'), (int, float))]
            bob[arm] = f"{min(tacs):.2f}" if tacs else '----'
        bob_s = '  '.join(f"{a}:{bob[a]}" for a in arms)
        pylog(f"  ── 進捗 {i+1}/{len(jobs)}  elapsed {_fmt_dur(elapsed)}  ETA {_fmt_dur(eta)}  "
              f"pace {_fmt_dur(med)}/run  best-of-batch(feasible) [{bob_s}]\n")

    pylog(f"\n==== バッチ完了 ====")
    ok = [r for r in records if r.get('exit_code') == 0 and r.get('subdir')]
    pylog(f"  成功 run: {len(ok)}/{len(records)}")
    for arm in arms:
        tacs = [(r['seed'], r['best_TAC']) for r in records
                if r.get('sampler') == arm and r.get('feasible') is True
                and isinstance(r.get('best_TAC'), (int, float))]
        if tacs:
            bs, bt = min(tacs, key=lambda x: x[1])
            pylog(f"  [{arm}] best-of-batch (feasible): seed={bs}  effective_TAC={bt:.2f} 億円/年  "
                  f"(feasible {len(tacs)}/{sum(1 for r in records if r.get('sampler')==arm)})")
    pylog(f"  manifest: {batch_dir / 'manifest.json'}")
    pylog(f"\n  次の分析コマンド:")
    pylog(f"    {_resolve_python()} verification\\analyze_seed_robustness.py --batch {batch_dir}")
    pylog.close()


if __name__ == '__main__':
    main()
