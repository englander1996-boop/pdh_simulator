r"""monitor_main.py — sub1(旧main)/sub2(旧final) 並列 run の「全体進捗」をざっくり表示。

別ターミナルで実行して、6 worker 合算の進捗 (X/N trial・feasible・best・ETA) を見る。
(旧 main.py→sub/sub1.py, 旧 final.py→sub/sub2.py。本 main.py=旧special は HYSYS 単一プロセスで本ツール対象外。)
共有 SQLite を read-only で読むだけなので最適化には干渉しない (WAL で並行読取り可)。
個別 worker の 1 trial 詳細を見たい時は従来どおり _worker*.log を tail。

使い方 (venv python で):
  python tools\monitor_main.py                 # 最新 run を 1 回表示
  python tools\monitor_main.py --watch 15      # 15 秒ごとに更新表示 (Ctrl+C で監視終了)
  python tools\monitor_main.py --db outputs\main_XXXX\optuna.db
"""
import os
import sys
import glob
import time
import argparse
import datetime

_REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'sub'))   # 旧 main.py→sub/sub1.py を import するため
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import optuna
optuna.logging.set_verbosity(optuna.logging.ERROR)


def find_latest_db():
    # main_* と final_* の両 run を対象 (どちらも study 名は pdh_<ts>)
    dbs = glob.glob('outputs/main_*/optuna.db') + glob.glob('outputs/final_*/optuna.db')
    return max(dbs, key=os.path.getmtime) if dbs else None


def study_name_from_db(db):
    # outputs/{main,final}_<ts>/optuna.db -> pdh_<ts>
    parent = os.path.basename(os.path.dirname(db))
    for pre in ('main_', 'final_'):
        if parent.startswith(pre):
            return 'pdh_' + parent[len(pre):]
    return None


def default_total():
    # sub1(旧 main.py)/sub2(旧 final.py) の並列 run を監視する用途なので sub1 の N_TRIALS を既定に。
    try:
        import sub1 as M
        return M.N_TRIALS
    except Exception:
        return None


def fmt_dur(sec):
    sec = int(max(0, sec))
    h, r = divmod(sec, 3600)
    m, s = divmod(r, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


def load(db):
    storage = optuna.storages.RDBStorage(
        url=f"sqlite:///{db}",
        engine_kwargs={"connect_args": {"timeout": 30.0}},
    )
    name = study_name_from_db(db)
    if name:
        try:
            return optuna.load_study(study_name=name, storage=storage)
        except Exception:
            pass
    # フォールバック: storage 内の最初の study
    sums = optuna.get_all_study_summaries(storage=storage)
    return optuna.load_study(study_name=sums[0].study_name, storage=storage)


def snapshot(db, total):
    study = load(db)
    trials = study.trials
    states = {}
    for t in trials:
        states[t.state.name] = states.get(t.state.name, 0) + 1
    done = states.get('COMPLETE', 0)
    running = states.get('RUNNING', 0)
    feas_trials = [t for t in trials
                   if t.state.name == 'COMPLETE'
                   and t.user_attrs.get('is_feasible') and t.value is not None]
    feas = len(feas_trials)

    # ETA: 完了 trial の開始時刻群から「全 worker 合算スループット」を出して線形外挿。
    # ざっくり用途。QMC→TPE で pace が変わるので厳密ではないが目安にはなる。
    comp = [t for t in trials
            if t.state.name == 'COMPLETE' and t.datetime_start and t.datetime_complete]
    eta_str, rate_str = "?", "?"
    if comp and total and done > 0:
        t0 = min(t.datetime_start for t in comp)
        now = datetime.datetime.now()
        elapsed = (now - t0).total_seconds()
        if elapsed > 0:
            rate = done / elapsed                       # trial/sec (合算)
            remaining = max(0, total - done)
            eta_sec = remaining / rate if rate > 0 else 0
            finish = now + datetime.timedelta(seconds=eta_sec)
            rate_str = f"{rate*60:.1f} trial/分"
            eta_str = f"残り~{fmt_dur(eta_sec)} → 完了予定 {finish.strftime('%H:%M')}"

    # best は「is_feasible な trial の最小 effective_TAC (= t.value)」。
    # study.best_value は constraints_func の制約充足 trial の最小を返すため、
    # is_feasible=True でも constraint 違反扱いの trial が除外され値が食い違う
    # (worker ログの "best" と一致させるためこちらを採用)。
    best = f"{min(t.value for t in feas_trials):.2f}" if feas_trials else "なし"

    pct = (done / total * 100) if total else 0
    bar_n = int(pct / 100 * 28)
    bar = '█' * bar_n + '░' * (28 - bar_n)
    stamp = datetime.datetime.now().strftime('%H:%M:%S')
    print(f"[{stamp}] [{bar}] {done}/{total} ({pct:.0f}%)  "
          f"running={running} feas={feas} best={best}  {rate_str}  {eta_str}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', default=None, help='optuna.db パス (省略時は最新 main_* を自動検出)')
    ap.add_argument('--total', type=int, default=None, help='目標 trial 数 (省略時 main.N_TRIALS)')
    ap.add_argument('--watch', type=float, default=0, help='秒間隔で更新 (0=1回のみ)')
    args = ap.parse_args()

    db = args.db or find_latest_db()
    if not db or not os.path.exists(db):
        print("optuna.db が見つかりません (outputs/main_*/optuna.db)")
        return
    total = args.total or default_total()
    print(f"db = {db}  (total = {total})", flush=True)

    if args.watch > 0:
        try:
            while True:
                snapshot(db, total)
                time.sleep(args.watch)
        except KeyboardInterrupt:
            print("\n監視終了 (最適化は別プロセスなので継続中)")
    else:
        snapshot(db, total)


if __name__ == '__main__':
    main()
