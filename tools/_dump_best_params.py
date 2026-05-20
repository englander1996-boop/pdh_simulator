"""SQLite から best trial の params + user_attrs を抽出。"""
import sqlite3
import json
import sys

db = r'z:\pdh_simulator\outputs\main_20260519_202617\optuna.db'
con = sqlite3.connect(db)
cur = con.cursor()

# Best trial (min value)
cur.execute(
    'SELECT t.trial_id, t.number, v.value '
    'FROM trials t JOIN trial_values v ON t.trial_id=v.trial_id '
    'WHERE t.state="COMPLETE" ORDER BY v.value LIMIT 10'
)
top = cur.fetchall()
print("--- Top 10 trials (BO value) ---")
for tid, num, val in top:
    print(f'  trial #{num} (id={tid}): value = {val:.4f}')

best_tid, best_num, best_val = top[0]
print(f'\n=== Best trial: #{best_num} (id={best_tid}), value = {best_val:.4f} ===')

# Get params
cur.execute(
    'SELECT param_name, param_value, distribution_json '
    'FROM trial_params WHERE trial_id=?', (best_tid,)
)
print("\n--- Params ---")
for name, val_raw, dist_json in cur.fetchall():
    # param_value is internal representation; for log distributions need conversion
    # We can also read distribution_json to be safe
    import math
    try:
        dist = json.loads(dist_json)
    except Exception:
        dist = {}
    # param_value is the actual numerical value for most distributions
    print(f'  {name:30s} = {val_raw}')

# Get user_attrs
cur.execute(
    "SELECT key, value_json FROM trial_user_attributes WHERE trial_id=?", (best_tid,)
)
print("\n--- user_attrs ---")
for k, v in cur.fetchall():
    print(f'  {k:30s} = {v}')

con.close()
