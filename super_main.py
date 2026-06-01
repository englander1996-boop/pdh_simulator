# -*- coding: utf-8 -*-
r"""super_main.py — BO 検証バッチの単一エントリ (UI)。

これを動かすだけで検証が完結する:
  main.py を seed/sampler を散らして多数回 (直列) まわす → 集計 → 可視化 まで自動。

成果物は main.py と同様、**1 回の実行 = 1 個のタイムスタンプ付きフォルダ**に内包される:
  outputs/super_main_<ts>/
    ├ runs/            各 main 実行の成果物 (main_<ts>/ … trials.csv, best.json, top*.txt)
    ├ logs/            各 run の BO ライブログ (run NN_<sampler>_seed<seed>.log)
    ├ plots/           可視化 PNG (best-TAC 分布 / 収束 / ヒスト / params 収束)
    ├ manifest.jsonl   完了 run 記録 (resume 用)
    ├ summary.txt      集計テキスト
    └ summary.csv      1 run 1 行の一覧

検証の中身 (詳細は superbatch/ 各モジュール):
  - tpe   : 本命。best-TAC 分布と best params の収束一致 (大域性の傍証) を見る。
  - cmaes : 別手法クロスチェック (同じ最適に行くか)。
  - random: 対照群 (BO がランダムより良いか = BO 正当性)。
  - champion は exp3 再評価で要再現確認 (ノイズの楽観バイアス除去)。

⚠ HYSYS は単一インスタンス ([[project_hysys_no_parallel]])。直列実行・手動 main.py との同時不可。
   resume: 中断しても再実行すれば未完了の super_main_<ts>/ を続きから埋める。

計画の編集は superbatch/config.py のみ (PLAN / N_TRIALS_PER_RUN / seed)。

使い方:
  .\.venv\Scripts\python.exe -u super_main.py
"""
import os
import time
import datetime

from superbatch import config, manifest, runner, aggregate, plots
from superbatch.util import log, fmt_dur


_README = """# super_main 検証バッチ {ts}

main.py を seed/sampler を散らして直列に回し、BO の検証 (ロバスト性・大域性の傍証・正当性) と
最良 TAC 探索を行ったバッチの成果物一式。

## 中身
- `runs/`   … 各 main 実行の成果物 (main_<ts>/ に trials.csv, best.json, top*.txt, README.md)
- `logs/`   … 各 run の BO ライブログ
- `plots/`  … 可視化 PNG (best-TAC 分布 / 収束 / ヒスト / params 収束)
- `manifest.jsonl` … 完了 run 記録 (resume 用)
- `summary.txt` / `summary.csv` … 集計

## 見方
1. `summary.txt` … sampler 別 best-TAC 分布 + champion + BO 正当性 + params 収束。
2. `plots/best_tac_by_sampler.png` … tpe/cmaes が random より下なら BO が効いている。
3. `plots/param_convergence.png` … TPE best が 1 点に集まれば大域の傍証。
4. champion は ★ exp3 で 3-5 回再評価して TAC 再現を確認 (単一 min を結果にしない)。
"""


def _write_readme():
    path = os.path.join(config.BATCH_DIR, 'README.md')
    try:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(_README.format(ts=os.path.basename(config.BATCH_DIR)))
    except Exception:
        pass


def main():
    # --- バッチ親フォルダを確定 (resume 対応) ---
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    batch_dir, resumed = config.find_or_create_batch_dir(ts)
    config.set_batch_dir(batch_dir)
    for d in (config.BATCH_DIR, config.LOG_DIR, config.RUNS_DIR, config.PLOTS_DIR):
        os.makedirs(d, exist_ok=True)
    _write_readme()

    # --- HYSYS 競合ガード (起動時のみ) ---
    act = runner.active_main_run()
    if act:
        log(f"⚠ 別の main.py が稼働中の可能性 (直近更新: {act})。")
        log("  HYSYS は単一インスタンスのため同時実行不可。完了を待って再実行してください。")
        return

    done = manifest.load_done()
    todo = [(s, sd) for (s, sd) in config.PLAN if (s, sd) not in done]
    durations = [r['dur_sec'] for r in done.values()
                 if isinstance(r.get('dur_sec'), (int, float))]

    log(f"super_main {'再開' if resumed else '開始'}: {os.path.basename(batch_dir)}")
    log(f"  計画 {len(config.PLAN)} run "
        f"(TPE {len(config.TPE_SEEDS)} + cmaes {len(config.CMAES_SEEDS)} + "
        f"random {len(config.RANDOM_SEEDS)})、完了済 {len(config.PLAN) - len(todo)}、"
        f"今回 {len(todo)} run。")
    log(f"  N_TRIALS/run={config.N_TRIALS_PER_RUN}、1 run ~1.8h 想定 → "
        f"残り ~{fmt_dur(len(todo) * config.EST_RUN_SEC)}。")
    log(f"  成果物: {batch_dir}  (Ctrl-C 中断でも済んだ run は保持・再開可)")

    batch_start = time.time()
    remaining = len(todo)
    for idx, (sampler, seed) in enumerate(config.PLAN, 1):
        if (sampler, seed) in done:
            continue
        rec, dur = runner.run_one(idx, len(config.PLAN), sampler, seed, remaining, durations)
        durations.append(dur)
        done[(sampler, seed)] = rec
        remaining -= 1
        # 各 run 完了ごとに集計/可視化を更新 (途中でも最新が見られる)
        aggregate.print_summary()
        plots.generate_all()
        log(f"    batch 経過 {fmt_dur(time.time() - batch_start)}、残り {remaining} run")

    if not todo:
        # 完了済みでも集計・可視化は最新化しておく
        aggregate.print_summary()
        plots.generate_all()

    log(f"super_main 終了。成果物一式 → {batch_dir}")


if __name__ == '__main__':
    main()
