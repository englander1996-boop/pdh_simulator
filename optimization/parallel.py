"""マルチプロセス並列最適化 (multi-core) のための worker 起動・worker 本体。

Optuna の n_jobs(スレッド) は GIL で CPU 律速の本問題を速くできず、かつ
penalty_scale(module global) がスレッド競合で壊れるため不可。代わりに **複数プロセスが
同じ study(共有 SQLite storage) を load_if_exists で共有** する分散最適化を使う。
各プロセスは単スレッド = GIL 無関係・penalty_scale はプロセス分離で安全。worker 同士は
DB 経由で running trial が見えるので constant_liar=True が効き、冗長サンプリングを抑える。
QMC startup は worker 毎に seed をずらして重複回避。

正しさ(学習方向)について:
  並列 async TPE は「方向が壊れる」のではなく「sample 効率がやや落ちる」だけ
  (各提案が in-flight の N-1 trial を知らない=少し古いモデル)。constant_liar と
  控えめな worker 数(3-4)でほぼ逐次同等の品質に到達し、wall-clock は ~N 倍速。

役割:
  - coordinator (sub1 の run_pipeline / main.main / sub2.main 内から呼ぶ): study を
    SQLite で作成 → spawn_workers() で N worker を起動・待機 → 既存の top-k/レポートをそのまま実行。
  - worker (本ファイルを `python -m optimization.parallel --worker ...` で起動):
    共有 study を load_study し、自前の sampler(TPE constant_liar + QMC seed=base+i) で
    担当 n_trials を最適化。

kind は対象ファイルを表す: 'main'(HYSYS+SM) / 'sub1'(FUG) / 'sub2'(SM/rigorous/SM)。
  main(HYSYS) は単一 COM で並列不可 (N_WORKERS=1) のため kind='main' 経路は実質 sub1/sub2 用。
"""
from __future__ import annotations

import os
import sys
import math
import sqlite3
import argparse
import subprocess
from typing import Callable

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
# worker から sub/ 配下のモジュールを bare import するため path 追加。
_SUB = os.path.join(_ROOT, 'sub')
if _SUB not in sys.path:
    sys.path.insert(0, _SUB)


# ---------------------------------------------------------------------------
# SQLite を WAL mode に (複数プロセス書込みの lock 競合を緩和)
# ---------------------------------------------------------------------------

def set_sqlite_wal(db_path: str) -> None:
    """SQLite db を WAL journal mode に設定 (永続)。

    WAL は「複数リーダー + 単一ライター同時可」で、Optuna の trial 書込み(短時間)が
    並列 worker から来ても lock 競合が起きにくくなる。journal_mode は db の永続属性
    なので 1 度設定すれば worker の接続にも効く。
    """
    try:
        con = sqlite3.connect(db_path, timeout=100.0)
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA synchronous=NORMAL;")
        con.close()
    except Exception as e:  # WAL 失敗は致命ではない (timeout 長めで続行)
        print(f"[parallel] set_sqlite_wal 失敗 (続行): {e}", flush=True)


def _storage_with_timeout(storage_url: str):
    """worker 用に connect timeout を長めに取った RDBStorage を返す。

    WAL でも瞬間的な書込み競合はあるので、busy timeout を 100s 取って待たせる
    (待てば成功する。即エラーで trial を落とさない)。
    """
    import optuna
    return optuna.storages.RDBStorage(
        url=storage_url,
        engine_kwargs={"connect_args": {"timeout": 100.0}},
    )


# ---------------------------------------------------------------------------
# objective 再構築 (worker プロセス内で kind から作る)
# ---------------------------------------------------------------------------

def _build_objective(kind: str) -> Callable:
    """worker プロセスで objective 関数を再構築する。

    kind='sub1' : sub1.SEARCH_SPACE / SOLVER_BO 等から make_objective で構築 (FUG)。
    kind='main' : main.objective をそのまま使う (HYSYS/SM、module global 依存)。
    kind='sub2' : sub2.objective をそのまま使う (SM/rigorous/SM + Stage2 全trial)。
    import 時に pipeline/main() は走らない (どれも if __name__=='__main__' ガード)。
    """
    if kind == 'sub1':
        import sub1 as M
        from config.load import load_operating_config
        from optimization.objective import make_objective
        cfg = load_operating_config()
        return make_objective(
            search_space=M.SEARCH_SPACE,
            solver_assignment=M.SOLVER_BO,
            config=cfg,
            apply_hi=M.APPLY_HI,
            apply_stage2=False,                 # BO ループでは Stage2 を回さない
            hi_dT_min_K=M.HI_DT_MIN_K,
            strict_recovery_check=M.STRICT_RECOVERY_BO,
            recovery_tolerance=M.RECOVERY_TOLERANCE,
            n_trials_total=M.N_TRIALS,
        )
    if kind == 'main':
        import main as S
        return S.objective
    if kind == 'sub2':
        # sub2 = (sm, rigorous, sm) + Stage2 全trial。main と同じく module global
        # (penalty_scale) 依存だが worker はプロセス分離で安全。
        import sub2 as Fm
        return Fm.objective
    raise ValueError(f"未知の kind: {kind!r} ('main' | 'sub1' | 'sub2')")


# ---------------------------------------------------------------------------
# worker 本体
# ---------------------------------------------------------------------------

def _run_worker(kind: str, study_name: str, storage_url: str,
                n_trials: int, seed: int, n_startup: int) -> None:
    """1 worker プロセスの処理: 共有 study を load して担当 n_trials を最適化。"""
    import optuna
    from optimization.study import make_sampler

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    # worker 毎に seed をずらす → QMC(Sobol) の点が重複しない / TPE の乱数も分離。
    # constant_liar=True で pending trial を悲観評価 → worker 同士が群がらない。
    sampler = make_sampler('tpe', seed=seed, n_startup=n_startup,
                           constraints_func=None, constant_liar=True)
    storage = _storage_with_timeout(storage_url)
    study = optuna.load_study(study_name=study_name, storage=storage, sampler=sampler)
    objective = _build_objective(kind)

    # compact callback (任意): live ログを worker 個別ファイル(=stdout)へ。
    callbacks = []
    try:
        if kind == 'sub1':
            from optimization.callbacks import make_compact_callback
            callbacks = [make_compact_callback(n_trials_total=n_trials)]
        elif kind == 'sub2':
            import sub2 as Fm
            callbacks = [Fm._make_final_callback(n_trials)]
        else:  # kind == 'main' (HYSYS+SM)
            import main as S
            callbacks = [S._make_main_callback(n_trials)]
    except Exception:
        callbacks = []

    study.optimize(
        objective, n_trials=n_trials, n_jobs=1,
        catch=(Exception,), callbacks=callbacks, show_progress_bar=False,
    )
    print(f"[worker seed={seed}] 完了: {n_trials} trial", flush=True)


# ---------------------------------------------------------------------------
# coordinator: worker を N 個 spawn して待つ
# ---------------------------------------------------------------------------

def split_trials(n_total: int, n_workers: int) -> list:
    """n_total を n_workers にできるだけ均等配分 (合計 = n_total)。"""
    base = n_total // n_workers
    rem = n_total % n_workers
    return [base + (1 if i < rem else 0) for i in range(n_workers)]


def spawn_workers(
    kind: str,
    study_name: str,
    storage_url: str,
    db_path: str,
    n_workers: int,
    n_trials_total: int,
    n_startup: int,
    base_seed: int,
    out_dir: str,
    low_priority: bool = True,
) -> list:
    """N worker サブプロセスを起動して全完了まで待つ。

    各 worker は `python -m optimization.parallel --worker ...` で起動し、stdout を
    out_dir/_worker{i}.log にリダイレクト。戻り値は worker log パスのリスト。

    low_priority=True (既定): worker を BelowNormal 優先度で起動し、ユーザーの前景
    作業 (IDE/HYSYS GUI 等) に CPU を譲る。8コア飽和してもマシンが軽いまま。
    """
    set_sqlite_wal(db_path)
    # worker を低優先度で起動し対話作業を妨げない (Windows)。
    _flags = 0
    if low_priority and os.name == 'nt':
        _flags = getattr(subprocess, 'BELOW_NORMAL_PRIORITY_CLASS', 0)
    per = split_trials(n_trials_total, n_workers)
    print(f"[parallel] kind={kind} workers={n_workers} trials={per} "
          f"(計 {sum(per)})  storage={storage_url}", flush=True)

    procs = []
    logs = []
    for i in range(n_workers):
        seed = base_seed + i * 1000
        logf = os.path.join(out_dir, f"_worker{i}.log")
        fout = open(logf, "w", encoding="utf-8")
        cmd = [
            sys.executable, "-m", "optimization.parallel", "--worker",
            "--kind", kind, "--study-name", study_name, "--storage", storage_url,
            "--n-trials", str(per[i]), "--seed", str(seed), "--n-startup", str(n_startup),
        ]
        env = dict(os.environ)
        env["PYTHONIOENCODING"] = "utf-8"
        # main の worker は各々 HYSYS インスタンスを立てるので PDH_PER_UNIT_LOG 等の
        # env はそのまま継承させる。
        p = subprocess.Popen(cmd, stdout=fout, stderr=subprocess.STDOUT,
                             env=env, cwd=_ROOT)
        procs.append((p, fout))
        logs.append(logf)
        print(f"[parallel] worker {i} 起動 PID={p.pid} seed={seed} "
              f"n_trials={per[i]} log={logf}", flush=True)

    # 全 worker 完了待ち。並列中 coordinator が沈黙すると進捗が見えないため、DB を
    # 定期(既定 60s)読みして集約進捗を 1 行表示する (詳細な per-trial ログは各
    # _worker*.log に出続ける)。env PDH_PARALLEL_PROGRESS_SEC で間隔調整。
    import time as _time
    _t0 = _time.time()
    _interval = float(os.environ.get('PDH_PARALLEL_PROGRESS_SEC', '60'))
    while any(p.poll() is None for p, _ in procs):
        _time.sleep(_interval)
        _print_parallel_progress(storage_url, study_name, n_trials_total, _t0)
    for i, (p, fout) in enumerate(procs):
        rc = p.wait()
        try:
            fout.close()
        except Exception:
            pass
        print(f"[parallel] worker {i} 終了 rc={rc}", flush=True)
    return logs


def _print_parallel_progress(storage_url: str, study_name: str,
                             n_total: int, t0: float) -> None:
    """共有 study を read-only で読み、集約進捗を 1 行表示 (失敗は黙ってスキップ)。"""
    try:
        import time
        import optuna
        st = optuna.load_study(
            study_name=study_name,
            storage=optuna.storages.RDBStorage(
                url=storage_url, engine_kwargs={"connect_args": {"timeout": 30.0}}),
        )
        comp = [t for t in st.trials if t.state.name == 'COMPLETE']
        running = sum(1 for t in st.trials if t.state.name == 'RUNNING')
        feas = [t for t in comp if t.user_attrs.get('is_feasible') and t.value is not None]
        best = f"{min(t.value for t in feas):.2f}" if feas else "----"
        done = len(comp)
        pct = 100.0 * done / max(n_total, 1)
        el = (time.time() - t0) / 60.0
        eta = ((n_total - done) / done * el) if done > 0 else 0.0
        print(f"[parallel] 進捗 {done}/{n_total} ({pct:.0f}%)  running={running}  "
              f"feasible {len(feas)}  best {best}  経過 {el:.0f}分  ETA ~{eta:.0f}分",
              flush=True)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# CLI (worker mode)
# ---------------------------------------------------------------------------

def _main_cli() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--kind", choices=["main", "sub1", "sub2"], required=True)
    ap.add_argument("--study-name", required=True)
    ap.add_argument("--storage", required=True)
    ap.add_argument("--n-trials", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-startup", type=int, required=True)
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    if args.worker:
        _run_worker(args.kind, args.study_name, args.storage,
                    args.n_trials, args.seed, args.n_startup)
        # main(HYSYS) の worker は HYSYS を確実にクローズ
        if args.kind == "main":
            try:
                from units.vle.hysys.provider import shutdown_default_provider
                shutdown_default_provider()
            except Exception:
                pass
    else:
        ap.error("--worker なしの直接実行は未対応 (coordinator は main.py / sub/sub1.py / sub/sub2.py から)")


if __name__ == "__main__":
    _main_cli()
