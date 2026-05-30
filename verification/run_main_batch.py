r"""run_main_batch.py — main.py(HYSYS+SM 全変数 BO、反応器 REACTOR_KIND で軸流21/径方向流22)を N 回逐次実行するバッチドライバ。

目的 (2026-05-29 ユーザー指示):
  同一シード(SEED=42)のまま main.py を 48 回回し、各 run の成果物を保存、
  後段の analyze_main_batch.py でまとめて分析する。

出力先の方針 (2026-05-29):
  - 各 run の生成物 (main_<ts>/) は **outputs/ にそのまま** (main.py の既定出力先を尊重)。
  - バッチ管理ファイル (manifest/各runログ/batch.log) と分析結果は **verification/** に集約。
    → 本スクリプトは verification/ 配下に置く。

設計の核 (本リポ精読に基づく判断):
  - main.py は SEED=42 固定 + HYSYS 以外は完全決定的 (反応器/PSA/膜/SM/再循環ソルバ/
    コスト計算はいずれも乱数源ゼロ)。よって run 間の差は実質「HYSYS(Dist2 COM)の
    非決定性が BO 軌道に与える揺らぎ」だけ。これを 48 回サンプリングして
    オプティマイザのロバスト性・局所解の分かれ方を測るのが本バッチの狙い。
  - HYSYS.Application はユーザーセッションに単一インスタンスの COM サーバーで、
    複数 worker が Dispatch すると同じ HYSYS を共有して結果が汚染する (main.py:85-88)。
    → 並列は不可。各 run は「別プロセスで main.py を素のまま起動」する。
      プロセスを分けることで HYSYS が毎回クリーンに起動し、前 run の汚染も持ち越さない。
      1 run がクラッシュしても次 run は無傷で続行できる(隔離)。
  - main.py は触らない (手厚くチューニング済み)。各 run は通常どおり
    outputs/main_<ts>/ を新規生成する。本ドライバは実行前後で outputs/ を
    スナップショット比較して「今 run が作った subdir」を特定し、manifest に記録する。

成果物 (verification/batch_<ts>/):
  manifest.json        : run ごとの {index, seed, subdir(outputs内), exit_code, wall_sec, best_TAC, feasible}
                         (各 run 完了ごとに上書き保存 = 中断耐性・再開判定に使用)
  run_NN.log           : 各 run の main.py 全出力 (flush 付きライブログをそのまま保存)
  batch.log            : バッチ全体の進捗ログ

使い方 (Python は .venv の 3.13):
  .\.venv\Scripts\python.exe verification\run_main_batch.py                  # 48 回 (既定)
  .\.venv\Scripts\python.exe verification\run_main_batch.py --runs 10        # 回数変更
  .\.venv\Scripts\python.exe verification\run_main_batch.py --resume verification\batch_<ts>
  .\.venv\Scripts\python.exe verification\run_main_batch.py --dry-run        # 起動コマンド確認のみ

逐次実行のため総時間は (1 run ~30-60分) × N。48 回で約 24-48 時間が目安。
flush 付きライブログなので、別ターミナルで verification\batch_<ts>\run_NN.log を tail すると
個別 run の trial 進捗が live で見える。本ドライバ自身も batch.log に進捗を書き出す。
"""

import os
import sys
import glob
import json
import time
import argparse
import datetime
import subprocess
from pathlib import Path

# Windows コンソール(cp932)で ▶/★ 等の非ASCII進捗記号を print してもクラッシュしないよう、
# main.py:52-55 と同様に stdout/stderr を utf-8 へ再設定する。
try:
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')
except Exception:
    pass


# ===========================================================================
# § 1. 設定
# ===========================================================================
N_RUNS_DEFAULT = 48          # ユーザー指示 (2026-05-29): 同一シードで 48 回
SEED_FIXED     = 42          # main.py の SEED と一致。同一シードで回す(本バッチの実験設計)。

REPO    = Path(__file__).resolve().parents[1]
OUTPUTS = REPO / 'outputs'        # 各 run の main_<ts>/ はここに出る (main.py 既定、変更しない)
VERIFY  = REPO / 'verification'   # バッチ管理ファイル・分析結果はここに集約
MAIN_PY = REPO / 'main.py'
VENV_PY = REPO / '.venv' / 'Scripts' / 'python.exe'

# コンソールへ echo する行のフィルタ (全行はログファイルへ、要点だけ画面へ)。
_ECHO_MARKERS = ('★', '====', '成果物', 'Traceback', 'Error', 'ERROR', '結果')
_ETA_ECHO_EVERY = 20         # 'ETA' を含む進捗バー行は 20 本に 1 本だけ画面へ


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
    """実行に使う Python を解決。.venv の 3.13 を最優先 (リポの前提環境)。"""
    if VENV_PY.exists():
        return str(VENV_PY)
    print(f"  [警告] .venv python が無い ({VENV_PY})。sys.executable で代替する。", flush=True)
    return sys.executable


def _snapshot_main_dirs() -> set:
    """現在 outputs/ 配下に存在する main_* subdir の集合。"""
    return {p.name for p in OUTPUTS.glob('main_*') if p.is_dir()}


def _read_best_json(subdir: Path) -> dict:
    """run subdir の best.json を読む (無ければ空 dict)。"""
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


# ===========================================================================
# § 3. 1 run の実行
# ===========================================================================
def _run_one(run_index: int, batch_dir: Path, pylog) -> dict:
    """main.py を 1 回サブプロセスで起動し、生成された run subdir(outputs内)を特定して記録を返す。

    出力は run_NN.log に全文保存し、要点だけコンソールへ echo する。
    """
    py = _resolve_python()
    run_log_path = batch_dir / f'run_{run_index:02d}.log'

    before = _snapshot_main_dirs()
    pylog(f"\n{'='*72}")
    pylog(f"  ▶ run #{run_index:02d}/{N_runs_global:02d}  seed={SEED_FIXED}  start {datetime.datetime.now():%H:%M:%S}")
    pylog(f"    cmd: {py} main.py   (cwd={REPO})")
    pylog(f"    log: {run_log_path}")
    pylog(f"{'='*72}")

    # main.py は flush 付きで stdout に出すので、行単位でストリームしてログ保存＋間引き echo。
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    # 設計判断: SEED は main.py 内の定数(=42)をそのまま使う(同一シード実験)。env では上書きしない。

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
                        print(f"    [#{run_index:02d}] {s.strip()}", flush=True)
                elif any(m in s for m in _ECHO_MARKERS):
                    print(f"    [#{run_index:02d}] {s.strip()}", flush=True)
            proc.wait()
            exit_code = proc.returncode
        except Exception as e:
            logf.write(f"\n[run_main_batch] サブプロセス起動失敗: {type(e).__name__}: {e}\n")
            pylog(f"    ✗ 起動失敗: {type(e).__name__}: {e}")
            exit_code = -999
    wall = time.perf_counter() - t0

    # 生成された subdir を特定 (実行後に増えた main_* )。main_<ts> は outputs/ に残す。
    after = _snapshot_main_dirs()
    new_dirs = sorted(after - before)
    subdir_name = new_dirs[-1] if new_dirs else None  # 通常ちょうど 1 個
    if len(new_dirs) > 1:
        pylog(f"    [警告] 新規 main_* が複数検出: {new_dirs} → 最新を採用: {subdir_name}")

    best_tac = None
    feasible = None
    best_num = None
    if subdir_name:
        bj = _read_best_json(OUTPUTS / subdir_name)
        if bj:
            best_tac = bj.get('effective_TAC')
            best_num = bj.get('number')
            ua = bj.get('user_attrs', {}) or {}
            feasible = ua.get('is_feasible')

    rec = {
        'index': run_index,
        'seed': SEED_FIXED,
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
    pylog(f"  ◀ run #{run_index:02d} 完了  exit={exit_code}  {_fmt_dur(wall)}  "
          f"best_TAC={tac_s}({feas_s})  subdir=outputs/{subdir_name}")
    return rec


# ===========================================================================
# § 4. manifest 保存 / 再開判定
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
N_runs_global = N_RUNS_DEFAULT   # echo 用 (run_one から参照)


def main():
    global N_runs_global
    ap = argparse.ArgumentParser(description='main.py を N 回逐次実行するバッチドライバ')
    ap.add_argument('--runs', type=int, default=N_RUNS_DEFAULT, help=f'実行回数 (既定 {N_RUNS_DEFAULT})')
    ap.add_argument('--resume', type=str, default=None,
                    help='既存 batch_<ts> ディレクトリを指定して未完了分を再開')
    ap.add_argument('--dry-run', action='store_true', help='起動コマンドを表示するだけで実行しない')
    args = ap.parse_args()

    N_runs_global = args.runs
    OUTPUTS.mkdir(exist_ok=True)   # main.py の出力先
    VERIFY.mkdir(exist_ok=True)    # バッチ管理・分析の出力先

    # ---- batch dir 決定 (新規 or 再開) ----
    if args.resume:
        batch_dir = Path(args.resume)
        if not batch_dir.is_absolute():
            batch_dir = REPO / batch_dir
        if not batch_dir.exists():
            print(f"  [エラー] 再開先が存在しない: {batch_dir}", flush=True)
            sys.exit(1)
        manifest = _load_manifest(batch_dir)
        records = manifest.get('runs', [])
        done_ok = [r for r in records if r.get('exit_code') == 0 and r.get('subdir')]
        start_index = len(records)
        print(f"  ▶ 再開: {batch_dir}  (記録済み {len(records)} run / 成功 {len(done_ok)})", flush=True)
    else:
        ts = _now_ts()
        batch_dir = VERIFY / f'batch_{ts}'
        batch_dir.mkdir(parents=True, exist_ok=True)
        records = []
        start_index = 0

    meta = {
        'created_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'n_runs_target': args.runs,
        'seed_fixed': SEED_FIXED,
        'note': '同一シード42で main.py を N 回逐次実行 (HYSYS 非決定性による BO ロバスト性の測定)',
        'run_outputs_location': 'outputs/main_<ts>/ (各 run の生成物はここ)',
        'main_py': str(MAIN_PY),
        'python': _resolve_python(),
    }

    if args.dry_run:
        print(f"  [dry-run] batch_dir = {batch_dir}")
        print(f"  [dry-run] 各 run コマンド: {_resolve_python()} main.py  (cwd={REPO}, env PYTHONIOENCODING=utf-8)")
        print(f"  [dry-run] 各 run 生成物 = outputs/main_<ts>/  (検出のみ、移動しない)")
        print(f"  [dry-run] 実行回数 = {args.runs} (start_index={start_index})")
        return

    pylog = _Tee(batch_dir / 'batch.log')
    pylog(f"==== run_main_batch: main.py を {args.runs} 回逐次実行 (seed={SEED_FIXED} 固定) ====")
    pylog(f"  batch_dir(管理/分析) = {batch_dir}")
    pylog(f"  run 生成物            = {OUTPUTS}/main_<ts>/")
    pylog(f"  python                = {meta['python']}")
    pylog(f"  開始 index = {start_index} / 目標 {args.runs}")
    pylog(f"  ※逐次実行。1 run ~30-60分。総時間 = N×(1 run)。中断時は --resume {batch_dir} で再開可。")

    _save_manifest(batch_dir, records, meta)

    batch_t0 = time.perf_counter()
    durations = [r['wall_sec'] for r in records if isinstance(r.get('wall_sec'), (int, float))]
    for i in range(start_index, args.runs):
        rec = _run_one(i, batch_dir, pylog)
        records.append(rec)
        _save_manifest(batch_dir, records, meta)
        if isinstance(rec.get('wall_sec'), (int, float)):
            durations.append(rec['wall_sec'])

        # ---- バッチ進捗 + ETA ----
        n_remain = args.runs - (i + 1)
        med = sorted(durations)[len(durations) // 2] if durations else 0.0
        eta = n_remain * med
        elapsed = time.perf_counter() - batch_t0
        feas_tacs = [r['best_TAC'] for r in records
                     if r.get('feasible') is True and isinstance(r.get('best_TAC'), (int, float))]
        bob = f"{min(feas_tacs):.2f}" if feas_tacs else '----'
        pylog(f"  ── 進捗 {i+1}/{args.runs}  elapsed {_fmt_dur(elapsed)}  ETA {_fmt_dur(eta)}  "
              f"pace {_fmt_dur(med)}/run  best-of-batch(feasible) {bob}\n")

    pylog(f"\n==== バッチ完了 ====")
    ok = [r for r in records if r.get('exit_code') == 0 and r.get('subdir')]
    pylog(f"  成功 run: {len(ok)}/{len(records)}")
    feas_tacs = [(r['index'], r['best_TAC']) for r in records
                 if r.get('feasible') is True and isinstance(r.get('best_TAC'), (int, float))]
    if feas_tacs:
        bi, bt = min(feas_tacs, key=lambda x: x[1])
        pylog(f"  best-of-batch (feasible): run #{bi:02d}  effective_TAC={bt:.2f} 億円/年")
    pylog(f"  manifest: {batch_dir / 'manifest.json'}")
    pylog(f"\n  次の分析コマンド:")
    pylog(f"    {_resolve_python()} verification\\analyze_main_batch.py --batch {batch_dir}")
    pylog.close()


if __name__ == '__main__':
    main()
