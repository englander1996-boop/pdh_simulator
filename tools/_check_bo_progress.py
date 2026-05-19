"""現在進行中の BO 結果を SQLite から覗く一時診断スクリプト。"""
import sqlite3
import sys

db_path = r'z:\pdh_simulator\outputs\main_20260519_202617\optuna.db'
con = sqlite3.connect(db_path)
cur = con.cursor()

cur.execute('SELECT state, COUNT(*) FROM trials GROUP BY state')
print("states:", cur.fetchall())

cur.execute(
    'SELECT t.trial_id, v.value FROM trials t JOIN trial_values v ON t.trial_id=v.trial_id '
    'WHERE t.state="COMPLETE" ORDER BY v.value LIMIT 10'
)
print("best 10 values:", cur.fetchall())

for thr in (9999, 5000, 2000, 1000, 500, 300, 200):
    cur.execute(
        f'SELECT COUNT(*) FROM trial_values WHERE value < {thr}'
    )
    print(f"trials with value < {thr}:", cur.fetchone()[0])

# proxy_penalty_total_okuyen が user_attr に入っているか
cur.execute(
    "SELECT trial_id, value_json FROM trial_user_attributes "
    "WHERE key='proxy_penalty_total_okuyen' LIMIT 5"
)
print("proxy_penalty samples:", cur.fetchall())

cur.execute(
    "SELECT COUNT(*) FROM trial_user_attributes WHERE key='proxy_penalty_total_okuyen'"
)
print("trials with proxy_penalty recorded:", cur.fetchone()[0])

# is_feasible
cur.execute(
    "SELECT value_json, COUNT(*) FROM trial_user_attributes "
    "WHERE key='is_feasible' GROUP BY value_json"
)
print("is_feasible counts:", cur.fetchall())

# failure_reason 上位サンプル
cur.execute(
    "SELECT value_json, COUNT(*) FROM trial_user_attributes "
    "WHERE key='failure_reason' AND value_json != '\"\"' "
    "GROUP BY value_json ORDER BY COUNT(*) DESC LIMIT 10"
)
print("failure_reason top:")
for row in cur.fetchall():
    s = row[0][:200] if row[0] else ''
    print(f"  {row[1]:4d}: {s}")

con.close()
