"""
optimization/feasibility.py — 収束/feasibility の分類学習 (L1: 事後解析)

BO 終了後の Optuna study から (params, is_feasible) を抽出し、Random Forest で
二値分類器を学習する。出力:
  - CV AUC スコア (5-fold stratified)
  - 特徴量重要度ランキング
  - 2D 散布図 (上位 2 特徴量で feasible/infeasible を色分け)

方針:
  - スクリーニング先行 (L2) や制約付き BO (L3) の input として利用想定。
  - sklearn が未インストールでも import を遅延 (ImportError は呼び出し時)。
  - 「feasible」の定義は target_type で切替:
      * 'convergence' : effective_TAC < SOLVER_FAILURE_THRESHOLD (= solver 収束)
      * 'spec'        : is_feasible == True (spec も満たした)
      * 'both'        : 両方を満たした (より厳格)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# optuna 必須 (study を受け取るため)
import optuna


# SOLVER_FAILURE_THRESHOLD は config.penalty.solver_failure_okuyen から導出する。
# 閾値 = penalty 値 - 1.0 で「ペナルティ確定」値より僅かに小さい値を取り、float 比較で
# 取りこぼしを防ぐ。operating.toml の値変更時の同期忘れを避けるため lazy load 化。
def _get_solver_failure_threshold() -> float:
    """operating.toml の penalty.solver_failure_okuyen から閾値を計算する (lazy)。"""
    from config.load import load_operating_config
    cfg = load_operating_config()
    return cfg.penalty.solver_failure_okuyen - 1.0


# Backward compat: 旧 API (module-level 定数) もサポート。
# 初回 import 時に config を 1 回だけロードする。
SOLVER_FAILURE_THRESHOLD = _get_solver_failure_threshold()


@dataclass
class FeasibilityAnalysis:
    """analyze_feasibility() の戻り値。"""
    n_total:            int
    n_feasible:         int
    feasible_rate:      float           # n_feasible / n_total
    feature_names:      List[str]
    importances:        List[Tuple[str, float]]   # 降順
    cv_auc_mean:        float
    cv_auc_std:         float
    classifier:         object          # sklearn 分類器 (再利用可)
    target_type:        str
    top2_features:      Tuple[str, str]


def extract_training_data(
    study:       optuna.Study,
    target_type: str = 'convergence',
):
    """Optuna study から (X, y, feature_names) を抽出。

    Parameters
    ----------
    target_type : str
        'convergence' : value < SOLVER_FAILURE_THRESHOLD なら feasible (y=1)
        'spec'        : user_attrs['is_feasible'] == True なら feasible
        'both'        : 両方を満たすなら feasible (より厳格)
    """
    import numpy as np

    completed = [t for t in study.trials
                 if t.state == optuna.trial.TrialState.COMPLETE
                 and t.value is not None]
    if not completed:
        raise ValueError("Study に完了 trial がありません。")

    # 特徴量名 (params の順序を保持、全 trial 共通と仮定)
    feature_names = list(completed[0].params.keys())

    X_list: List[List[float]] = []
    y_list: List[int] = []
    for t in completed:
        row = [float(t.params.get(k, np.nan)) for k in feature_names]
        if any(np.isnan(v) for v in row):
            continue   # 一部 params 欠落 trial は除外
        X_list.append(row)

        # ラベル決定
        is_converged = t.value < SOLVER_FAILURE_THRESHOLD
        is_spec_ok   = bool(t.user_attrs.get('is_feasible', False))
        if target_type == 'convergence':
            y_list.append(int(is_converged))
        elif target_type == 'spec':
            y_list.append(int(is_spec_ok))
        elif target_type == 'both':
            y_list.append(int(is_converged and is_spec_ok))
        else:
            raise ValueError(f"未知の target_type: {target_type!r}")

    X = np.asarray(X_list, dtype=float)
    y = np.asarray(y_list, dtype=int)
    return X, y, feature_names


def train_classifier(
    X, y,
    *,
    model:       str  = 'rf',
    n_estimators: int = 200,
    max_depth:    Optional[int] = None,
    random_state: int = 42,
    cv:           int = 5,
):
    """Random Forest 等の二値分類器を学習し、CV AUC を返す。

    Returns
    -------
    clf : sklearn 分類器 (全データで再学習済)
    cv_auc_mean : float
    cv_auc_std  : float
    """
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos == 0 or n_neg == 0:
        # クラス不均衡で分類不可、ダミー結果を返す
        if model == 'rf':
            clf = RandomForestClassifier(n_estimators=n_estimators,
                                         max_depth=max_depth,
                                         random_state=random_state)
        else:
            clf = Pipeline([('scaler', StandardScaler()),
                            ('lr', LogisticRegression(max_iter=2000,
                                                      random_state=random_state))])
        return clf, float('nan'), float('nan')

    if model == 'rf':
        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=random_state,
            n_jobs=-1,
            class_weight='balanced',
        )
    elif model == 'logreg':
        clf = Pipeline([
            ('scaler', StandardScaler()),
            ('lr', LogisticRegression(
                max_iter=2000, random_state=random_state, class_weight='balanced',
            )),
        ])
    else:
        raise ValueError(f"未知の model: {model!r} (許容: 'rf' | 'logreg')")

    cv_folds = min(cv, n_pos, n_neg)
    if cv_folds < 2:
        cv_mean, cv_std = float('nan'), float('nan')
    else:
        skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
        scores = cross_val_score(clf, X, y, cv=skf, scoring='roc_auc', n_jobs=-1)
        cv_mean, cv_std = float(scores.mean()), float(scores.std())

    clf.fit(X, y)
    return clf, cv_mean, cv_std


def feature_importance_report(clf, feature_names: List[str]) -> List[Tuple[str, float]]:
    """学習済み分類器の特徴量重要度を降順リストで返す。

    RF は `feature_importances_`、LogReg は `|coef_|` を使う。
    """
    import numpy as np
    if hasattr(clf, 'feature_importances_'):
        imp = clf.feature_importances_
    elif hasattr(clf, 'named_steps'):
        # Pipeline (logreg)
        lr = clf.named_steps.get('lr')
        if lr is None or not hasattr(lr, 'coef_'):
            raise ValueError("分類器から重要度を取り出せません。")
        imp = np.abs(lr.coef_[0])
        imp = imp / imp.sum() if imp.sum() > 0 else imp
    else:
        raise ValueError("未対応の分類器型")

    pairs = list(zip(feature_names, imp.tolist()))
    pairs.sort(key=lambda kv: kv[1], reverse=True)
    return pairs


def plot_feasibility_2d(
    X, y,
    feature_names: List[str],
    var_x:    str,
    var_y:    str,
    path:     Path | str,
    title:    str = "Feasibility 2D scatter",
) -> None:
    """指定 2 変数で散布図、feasible (y=1) を青、infeasible (y=0) を赤。"""
    import matplotlib
    matplotlib.use('Agg')   # GUI なし backend
    import matplotlib.pyplot as plt
    import numpy as np

    ix = feature_names.index(var_x)
    iy = feature_names.index(var_y)
    feasible    = (y == 1)
    infeasible  = (y == 0)

    fig, ax = plt.subplots(figsize=(7, 6))
    if infeasible.any():
        ax.scatter(X[infeasible, ix], X[infeasible, iy],
                   c='tab:red', alpha=0.4, s=18, label=f'infeasible (n={infeasible.sum()})')
    if feasible.any():
        ax.scatter(X[feasible, ix], X[feasible, iy],
                   c='tab:blue', alpha=0.7, s=28, label=f'feasible (n={feasible.sum()})',
                   edgecolors='white', linewidths=0.5)
    ax.set_xlabel(var_x)
    ax.set_ylabel(var_y)
    ax.set_title(title)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def analyze_feasibility(
    study:       optuna.Study,
    output_dir:  Path | str,
    prefix:      str = 'feasibility',
    target_type: str = 'convergence',
    model:       str = 'rf',
) -> Optional[FeasibilityAnalysis]:
    """L1 事後解析の最上位関数。study → CSV/PNG/txt を出力。

    Parameters
    ----------
    study : optuna.Study
        BO 終了後の study (SQLite or in-memory)。
    output_dir : Path | str
        出力先ディレクトリ。
    prefix : str
        出力ファイルのプレフィックス (例: 'main_20260514_213319_feasibility')。
    target_type : str
        'convergence' | 'spec' | 'both' (extract_training_data の定義参照)
    model : str
        'rf' | 'logreg'

    Returns
    -------
    FeasibilityAnalysis | None
        完了 trial が無い・分類不可ならログ出力して None。
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    txt_path = output_dir / f'{prefix}.txt'
    png_path = output_dir / f'{prefix}_2d.png'

    try:
        X, y, feature_names = extract_training_data(study, target_type=target_type)
    except ValueError as e:
        msg = f"[feasibility] スキップ: {e}\n"
        txt_path.write_text(msg, encoding='utf-8')
        print(msg.rstrip())
        return None

    n_total = len(y)
    n_pos   = int((y == 1).sum())
    rate    = n_pos / max(1, n_total)

    if n_pos == 0:
        msg = (f"[feasibility] 全 {n_total} trial が y=0 ({target_type})、"
               f"分類不可。探索空間の見直しが必要。\n")
        txt_path.write_text(msg, encoding='utf-8')
        print(msg.rstrip())
        return None
    if n_pos == n_total:
        msg = (f"[feasibility] 全 {n_total} trial が y=1 ({target_type})、"
               f"分類は自明 (全 feasible)。\n")
        txt_path.write_text(msg, encoding='utf-8')
        print(msg.rstrip())
        return None

    clf, cv_auc_mean, cv_auc_std = train_classifier(X, y, model=model)
    importances = feature_importance_report(clf, feature_names)
    top2 = (importances[0][0], importances[1][0])

    # txt レポート
    lines = []
    lines.append("=" * 80)
    lines.append(f"Feasibility 分類解析 ({prefix})")
    lines.append("=" * 80)
    lines.append(f"  target_type      = {target_type}")
    lines.append(f"  model            = {model}")
    lines.append(f"  完了 trial 数     = {n_total}")
    lines.append(f"  feasible 数       = {n_pos} ({rate*100:.1f}%)")
    lines.append(f"  infeasible 数     = {n_total - n_pos}")
    lines.append(f"  CV AUC (5-fold)  = {cv_auc_mean:.4f} ± {cv_auc_std:.4f}")
    lines.append("")
    lines.append("特徴量重要度 (降順):")
    lines.append(f"  {'rank':>4} {'feature':<22} {'importance':>12}")
    lines.append(f"  {'-'*4} {'-'*22} {'-'*12}")
    for rank, (name, imp) in enumerate(importances, start=1):
        lines.append(f"  {rank:>4} {name:<22} {imp:>12.4f}")
    lines.append("")
    lines.append(f"2D 散布図: {png_path.name} (x={top2[0]}, y={top2[1]})")
    lines.append("=" * 80)
    txt_path.write_text("\n".join(lines) + "\n", encoding='utf-8')

    # 2D 散布図
    try:
        plot_feasibility_2d(
            X, y, feature_names, top2[0], top2[1], png_path,
            title=f"{prefix}: {target_type} feasibility ({n_pos}/{n_total})",
        )
    except Exception as e:
        print(f"[feasibility] 2D 散布図作成失敗: {type(e).__name__}: {e}")

    print(f"[feasibility] レポート → {txt_path}")
    print(f"[feasibility] 2D 散布図 → {png_path}")
    print(f"[feasibility] CV AUC = {cv_auc_mean:.3f} ± {cv_auc_std:.3f}, "
          f"feasible rate = {rate*100:.1f}%")
    print(f"[feasibility] Top-3 重要特徴量: "
          f"{', '.join(f'{n}({i:.3f})' for n, i in importances[:3])}")

    return FeasibilityAnalysis(
        n_total       = n_total,
        n_feasible    = n_pos,
        feasible_rate = rate,
        feature_names = feature_names,
        importances   = importances,
        cv_auc_mean   = cv_auc_mean,
        cv_auc_std    = cv_auc_std,
        classifier    = clf,
        target_type   = target_type,
        top2_features = top2,
    )
