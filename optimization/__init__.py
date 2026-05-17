"""
optimization モジュール — 多変数最適化 + top-k 候補の re-evaluation。

構成:
  - hen_synthesis.py: HEN (Heat Exchanger Network) の Pinch Design Method 合成
                      Stage 2 (= Stage 1 ピンチ targeting に対する詳細設計)
  - search_space.py:  SEARCH_SPACE スキーマと params → FlowsheetDesignVars 変換
  - objective.py:     Optuna objective 関数ファクトリ
  - study.py:         Optuna Study 生成と最適化ループ  (optuna 必須)
  - topk.py:          BO 上位候補の rigorous + Stage 2 再評価  (optuna 必須)
  - reporting.py:     CSV / JSON / txt 比較レポート出力       (optuna 必須)

main.py は本モジュールを薄くオーケストレーションする。

依存:
  - optuna が未インストールでも `hen_synthesis` / `search_space` /
    `objective` (の make_objective ファクトリ) は使用可能。
  - 未インストール時は `study` / `topk` / `reporting` の import がスキップされ、
    本パッケージから当該シンボルは公開されない (= main.py は動かない)。
    その場合、import 時に ImportWarning が発行される。
"""

import warnings as _warnings

# ---- optuna 可用性チェック ----
try:
    import optuna as _optuna  # noqa: F401
    _HAS_OPTUNA = True
except ImportError:
    _HAS_OPTUNA = False
    # 設計判断: ImportWarning はデフォルトで非表示。UserWarning を使って必ず表示させる。
    _warnings.warn(
        "optuna が未インストールです。多変数最適化機能 (main.py / optimization.study など) は\n"
        "使用できません。インストール方法:\n"
        "    .\\.venv\\Scripts\\python.exe -m pip install -r requirements.txt\n"
        "  または\n"
        '    .\\.venv\\Scripts\\python.exe -m pip install "optuna>=3.5"',
        UserWarning,
        stacklevel=2,
    )

# ---- optuna 不要なシンボル (常に load) ----
from optimization.hen_synthesis import (
    HEMatch,
    HENResult,
    synthesize_hen,
    apply_synthesis_to_economics,
)

from optimization.search_space import (
    VarSpec,
    EXPECTED_KEYS,
    DEFAULT_BASELINE,
    P_L_FIXED_PA,
    validate_search_space,
    suggest_params,
    build_design,
)

from optimization.objective import make_objective


__all__ = [
    # hen_synthesis
    'HEMatch', 'HENResult', 'synthesize_hen', 'apply_synthesis_to_economics',
    # search_space
    'VarSpec', 'EXPECTED_KEYS', 'DEFAULT_BASELINE', 'P_L_FIXED_PA',
    'validate_search_space', 'suggest_params', 'build_design',
    # objective
    'make_objective',
]

# ---- optuna 必須シンボル (条件付き load) ----
if _HAS_OPTUNA:
    from optimization.study import (
        make_sampler,
        create_study,
        run_optimization,
    )
    from optimization.topk import (
        TopKEntry,
        select_topk,
        reevaluate_topk,
        best_entry,
    )
    from optimization.reporting import (
        save_trials_csv,
        save_best_json,
        save_topk_report,
    )
    __all__ += [
        # study
        'make_sampler', 'create_study', 'run_optimization',
        # topk
        'TopKEntry', 'select_topk', 'reevaluate_topk', 'best_entry',
        # reporting
        'save_trials_csv', 'save_best_json', 'save_topk_report',
    ]

    # ---- L1: 収束/feasibility 分類学習 (sklearn 必須、optional) ----
    try:
        import sklearn  # noqa: F401
        from optimization.feasibility import (
            FeasibilityAnalysis,
            extract_training_data,
            train_classifier,
            feature_importance_report,
            plot_feasibility_2d,
            analyze_feasibility,
        )
        __all__ += [
            'FeasibilityAnalysis',
            'extract_training_data', 'train_classifier',
            'feature_importance_report', 'plot_feasibility_2d',
            'analyze_feasibility',
        ]
        _HAS_SKLEARN = True
    except ImportError:
        _HAS_SKLEARN = False
        _warnings.warn(
            "scikit-learn が未インストールです。Feasibility 分類解析 (L1) は使用不可。\n"
            '    .\\.venv\\Scripts\\python.exe -m pip install "scikit-learn>=1.4"',
            UserWarning,
            stacklevel=2,
        )
