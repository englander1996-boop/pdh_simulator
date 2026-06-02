r"""
tools/run_feasibility.py — 既存 Optuna study に対する事後 feasibility 解析

main.py の SQLite storage (outputs/main_*.db) を後から読み込んで分類学習する。
main.py 実行時に解析を回し損ねた場合や、target_type / model を差し替えて
再解析する場合に使う。

使い方:
  .\.venv\Scripts\python.exe tools/run_feasibility.py outputs/main_<timestamp>.db
  .\.venv\Scripts\python.exe tools/run_feasibility.py outputs/main_<timestamp>.db --target spec
  .\.venv\Scripts\python.exe tools/run_feasibility.py outputs/main_<timestamp>.db --model logreg
"""

import argparse
import os
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

import optuna

from optimization.feasibility import analyze_feasibility


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('db_path', type=str,
                        help='Optuna SQLite ファイルパス (例: outputs/main_xxx.db)')
    parser.add_argument('--study-name', type=str, default=None,
                        help='study 名 (省略時は DB 内で最も新しい study を採用)')
    parser.add_argument('--target', type=str, default='convergence',
                        choices=['convergence', 'spec', 'both'],
                        help="ラベル定義 (default: 'convergence')")
    parser.add_argument('--model', type=str, default='rf',
                        choices=['rf', 'logreg'],
                        help="分類器 (default: 'rf')")
    parser.add_argument('--output-dir', type=str, default=None,
                        help='出力先ディレクトリ (default: DB ファイルと同じ subdir)')
    parser.add_argument('--prefix', type=str, default=None,
                        help='出力ファイルのプレフィックス (default: "feasibility_<target>_<model>")')
    args = parser.parse_args()

    db_path = Path(args.db_path).resolve()
    if not db_path.exists():
        print(f"ERROR: DB ファイルが見つかりません: {db_path}", file=sys.stderr)
        sys.exit(1)
    storage_url = f'sqlite:///{db_path.as_posix()}'

    # main.py は run ごと subdir 構造を採用 (outputs/main_<ts>/)。
    # output_dir 未指定なら DB と同じ subdir に出す。
    output_dir = Path(args.output_dir) if args.output_dir else db_path.parent

    # study 名の解決
    if args.study_name is None:
        names = optuna.get_all_study_names(storage=storage_url)
        if not names:
            print(f"ERROR: DB 内に study がありません: {db_path}", file=sys.stderr)
            sys.exit(1)
        study_name = names[-1]   # 最新を採用
        print(f"[info] study_name 自動選択: {study_name}")
    else:
        study_name = args.study_name

    study = optuna.load_study(study_name=study_name, storage=storage_url)
    print(f"[info] study '{study_name}' をロード ({len(study.trials)} trial)")

    prefix = args.prefix or f'feasibility_{args.target}_{args.model}'
    analyze_feasibility(
        study       = study,
        output_dir  = output_dir,
        prefix      = prefix,
        target_type = args.target,
        model       = args.model,
    )


if __name__ == '__main__':
    main()
