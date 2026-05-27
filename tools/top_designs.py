r"""top_designs.py — main.py 並列 run の feasible 設計を「収率・Dist3段数・利益」軸で一覧。

BO の目的関数は effective_TAC だが、設計上注目したい 収率 / N_dist3 / 利益 を並べて
feasible 上位を確認する。共有 SQLite を read-only 参照 (最適化に干渉しない)。

定義:
  収率 yield    = production_kmol_h / F_C3H8_fresh_used_kmol_h * 100  (= worker ログの yield)
  利益 profit   = profit_hi_okuyen  (Stage1/HI 適用後、億円/年。BO ループは Stage2 非実行)
  Dist3 段数    = params['N_dist3']

使い方 (venv python で):
  python tools\top_designs.py                 # TAC 順で feasible 上位 10
  python tools\top_designs.py --k 15 --watch 20
  python tools\top_designs.py --sort profit   # 利益順  (tac|profit|yield|n3)
"""
import os
import sys
import glob
import time
import argparse
import datetime

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import optuna
optuna.logging.set_verbosity(optuna.logging.ERROR)


def find_latest_db():
    # main_* と final_* の両 run を対象
    dbs = glob.glob('outputs/main_*/optuna.db') + glob.glob('outputs/final_*/optuna.db')
    return max(dbs, key=os.path.getmtime) if dbs else None


def load(db):
    storage = optuna.storages.RDBStorage(
        url=f"sqlite:///{db}",
        engine_kwargs={"connect_args": {"timeout": 30.0}},
    )
    parent = os.path.basename(os.path.dirname(db))
    for pre in ('main_', 'final_'):
        if parent.startswith(pre):
            try:
                return optuna.load_study(study_name='pdh_' + parent[len(pre):], storage=storage)
            except Exception:
                pass
    sums = optuna.get_all_study_summaries(storage=storage)
    return optuna.load_study(study_name=sums[0].study_name, storage=storage)


def metrics(t):
    ua = t.user_attrs
    prod = ua.get('production_kmol_h')
    fresh = ua.get('F_C3H8_fresh_used_kmol_h')
    yld = (prod / fresh * 100) if (prod and fresh) else None
    return {
        'tac':    t.value,
        'profit': ua.get('profit_hi_okuyen'),
        'yield':  yld,
        # Dist3 段数: main は 'N_dist3'、final は 'col3_n_stages'
        'n3':     t.params.get('N_dist3', t.params.get('col3_n_stages')),
        'prod':   prod,
        'target': ua.get('target_kmol_h'),
    }


def show(db, k, sort):
    study = load(db)
    feas = [t for t in study.trials if t.state.name == 'COMPLETE'
            and t.value is not None and t.user_attrs.get('is_feasible')]
    stamp = datetime.datetime.now().strftime('%H:%M:%S')
    if not feas:
        print(f"[{stamp}] feasible まだ無し")
        return
    rows = [(t, metrics(t)) for t in feas]
    if sort == 'profit':
        rows.sort(key=lambda r: -(r[1]['profit'] if r[1]['profit'] is not None else -9e9))
    elif sort == 'yield':
        rows.sort(key=lambda r: -(r[1]['yield'] if r[1]['yield'] is not None else -9e9))
    elif sort == 'n3':
        rows.sort(key=lambda r: (r[1]['n3'] if r[1]['n3'] is not None else 9e9))
    else:  # tac
        rows.sort(key=lambda r: r[1]['tac'])

    print(f"[{stamp}] feasible {len(feas)} 件  (sort={sort}、上位 {min(k, len(rows))})")
    print(f"  {'rank':>4} {'trial':>6} {'TAC':>9} {'利益':>9} {'収率':>7} {'N_d3':>5} {'生産/目標':>13}")
    for i, (t, m) in enumerate(rows[:k], 1):
        tac = f"{m['tac']:.1f}"
        pf  = f"{m['profit']:+.1f}" if m['profit'] is not None else "?"
        yl  = f"{m['yield']:.1f}%" if m['yield'] is not None else "?"
        n3  = m['n3'] if m['n3'] is not None else '?'
        pt  = f"{m['prod']:.0f}/{m['target']:.0f}" if (m['prod'] and m['target']) else "?"
        print(f"  {i:>4} #{t.number:>5} {tac:>9} {pf:>9} {yl:>7} {str(n3):>5} {pt:>13}")

    ys = [m['yield'] for _, m in rows if m['yield'] is not None]
    ns = [m['n3'] for _, m in rows if m['n3'] is not None]
    ps = [m['profit'] for _, m in rows if m['profit'] is not None]
    if ys and ns and ps:
        print(f"  ── feasible 全体レンジ: 収率 {min(ys):.1f}〜{max(ys):.1f}%  / "
              f"N_d3 {min(ns)}〜{max(ns)}  / 利益 {min(ps):+.0f}〜{max(ps):+.0f}億")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', default=None)
    ap.add_argument('--k', type=int, default=10)
    ap.add_argument('--sort', choices=['tac', 'profit', 'yield', 'n3'], default='tac')
    ap.add_argument('--watch', type=float, default=0)
    args = ap.parse_args()
    db = args.db or find_latest_db()
    if not db or not os.path.exists(db):
        print("optuna.db が見つかりません")
        return
    print(f"db = {db}")
    if args.watch > 0:
        try:
            while True:
                show(db, args.k, args.sort)
                time.sleep(args.watch)
        except KeyboardInterrupt:
            print("\n監視終了 (最適化は継続中)")
    else:
        show(db, args.k, args.sort)


if __name__ == '__main__':
    main()
