"""parallel.py (main 並列パス) の小規模実テスト。

2 worker × 2 trial = 実フローシート 4 評価を共有 SQLite study で実行し、
worker が DB に trial を書き込み coordinator が読めることを確認する。
現走行の main/special とは別 DB・別 study なので無干渉。
"""
import os, sys, datetime
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
from optimization.study import make_sampler
from optimization.parallel import spawn_workers

ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
out_dir = 'outputs'
db = os.path.join(out_dir, f'_test_parallel_{ts}.db')
storage = f'sqlite:///{db}'
study_name = f'_test_par_{ts}'

# coordinator: 空の study を作成 (sampler は worker 側が各自持つので任意)
optuna.create_study(study_name=study_name, storage=storage, direction='minimize',
                    sampler=make_sampler('tpe', 42, 2), load_if_exists=True)
print(f"study 作成: {study_name}  db={db}", flush=True)

# 2 worker × 2 trial = 計 4 trial
spawn_workers(kind='main', study_name=study_name, storage_url=storage, db_path=db,
              n_workers=2, n_trials_total=4, n_startup=2, base_seed=999, out_dir=out_dir)

# coordinator が DB を読み直して全 trial を見られるか
study2 = optuna.load_study(study_name=study_name, storage=storage)
comp = [t for t in study2.trials if t.state.name == 'COMPLETE']
print(f"\n=== 検証 ===")
print(f"trials total = {len(study2.trials)}  complete = {len(comp)}")
# どの worker(seed) が何 trial 入れたか (params に seed は無いが trial 数で確認)
feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
print(f"feasible = {len(feas)}")
try:
    print(f"best_value = {study2.best_value:.2f}")
except Exception as e:
    print(f"best 取得不可 (feasible 無し等): {type(e).__name__}")
# 各 trial の effective_TAC を列挙
for t in study2.trials:
    fr = t.user_attrs.get('failure_unit', '?')
    print(f"  trial #{t.number} state={t.state.name} value={t.value} unit={fr}")
print("=== parallel main パス: trial が DB に書かれ coordinator が読めれば OK ===")
