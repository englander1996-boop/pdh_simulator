"""
special.py — HYSYS バックエンドでの蒸留塔最適化 (Optuna BO)

exp3.py を評価関数として、各塔の運転条件・段数・フィード段・主スペックを
Bayesian Optimization で探索する。

設計判断 (2026-05-22):
  - main.py / optimization/study.py とは独立スクリプト。
    探索変数が HYSYS の主スペック (Comp Fraction-2 / Reflux Ratio / Draw Rate) を
    使うため、main の ColumnTunables.reflux_ratio / recovery_* 系の探索空間とは
    別物として扱う。
  - 反応器系 (Swing/PSA/膜) は exp1 の best 値 (trial #115) で固定。
    「上流ノータッチ」のユーザー方針に従う。
  - HYSYS 非収束は effective_TAC = PENALTY_VALUE で BO の探索方向から除外。
  - HI/Stage 2 は exp3 と同じく True で評価 (公平比較)。

探索変数 (12 変数):
  Dist1: col1_p_kpa, col1_n_stages, col1_feed_stage, col1_comp_frac_2
  Dist2: col2_p_kpa, col2_n_stages, col2_feed_stage, col2_reflux_ratio
  Dist3: col3_p_kpa, col3_n_stages, col3_feed_stage, col3_draw_rate_kmolh

使い方:
  N_TRIALS を変更し、.\.venv\Scripts\python.exe special.py で実行。
  storage を有効にすると中断・再開可能 (SQLite ファイル)。
"""

import os
import sys
import math
import time
from typing import Optional

# Windows コンソール対応
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import optuna

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars


# ===========================================================================
#  BO 設定
# ===========================================================================

N_TRIALS                = 2         # 試行回数 (まず動かす目的、慣れたら増やす)
PENALTY_VALUE_OKUYEN    = 1.0e6     # 非収束時の effective_TAC [億円/年]
SAMPLER_SEED            = 42
USE_SQLITE_STORAGE      = False     # True で SQLite 中断/再開対応
SQLITE_PATH             = "special_study.db"
STUDY_NAME              = "pdh_hysys_special"
SHOW_HYSYS_GUI          = False     # True にすると HYSYS を可視化 (デバッグ用)


# ===========================================================================
#  上流 (反応器・PSA・膜) は exp1 best (trial #115) 固定
# ===========================================================================
T_in_K        = 955.6260
z_cat_m       = 21.0550
t_cyc_min     = 14.2111
D_reactor_m   = 8.4363

D_psa_col_m       = 3.3187
L_psa_bed_m       = 25.6086
desorption_target = 0.2805

P_H_Pa     = 7.5403e5
P_L_Pa     = 1.0e5
A_mem_m2   = 1.2085e5

F_C3H8_fresh_kmol_h = 1647.5800

APPLY_HI     = True
APPLY_STAGE2 = True
HI_DT_MIN_K  = 10.0


# ===========================================================================
#  探索範囲 (lhs_column*.py の範囲を参考に。HSC 存在段数で制約)
# ===========================================================================
SEARCH_SPACE = {
    # Dist1
    "col1_p_kpa":         (1700.0, 2200.0),   # kPa
    "col1_n_stages":      (30, 60),           # HSC 30〜60
    "col1_feed_stage":    (5, 50),            # 5..min(N-2, 50)
    "col1_comp_frac_2":   (0.90, 0.999),

    # Dist2
    "col2_p_kpa":         (700.0, 1000.0),
    "col2_n_stages":      (15, 80),           # HSC 15〜80
    "col2_feed_stage":    (3, 70),
    "col2_reflux_ratio":  (5.0, 15.0),

    # Dist3
    "col3_p_kpa":         (1600.0, 2200.0),
    "col3_n_stages":      (69, 200),          # HSC 69〜200
    "col3_feed_stage":    (10, 180),
    "col3_draw_rate_kmolh": (800.0, 1400.0),
}


# ===========================================================================
#  Objective
# ===========================================================================

def _build_design_from_trial(trial: optuna.trial.Trial) -> FlowsheetDesignVars:
    """trial から探索変数をサンプリングして FlowsheetDesignVars を作る。"""
    s = SEARCH_SPACE

    # Dist1
    p1 = trial.suggest_float("col1_p_kpa", *s["col1_p_kpa"])
    n1 = trial.suggest_int("col1_n_stages", *s["col1_n_stages"])
    # feed_stage は段数より小さくする
    fs1_hi = min(s["col1_feed_stage"][1], n1 - 2)
    fs1_lo = min(s["col1_feed_stage"][0], fs1_hi)
    fs1 = trial.suggest_int("col1_feed_stage", fs1_lo, fs1_hi)
    cf1 = trial.suggest_float("col1_comp_frac_2", *s["col1_comp_frac_2"])

    # Dist2
    p2 = trial.suggest_float("col2_p_kpa", *s["col2_p_kpa"])
    n2 = trial.suggest_int("col2_n_stages", *s["col2_n_stages"])
    fs2_hi = min(s["col2_feed_stage"][1], n2 - 2)
    fs2_lo = min(s["col2_feed_stage"][0], fs2_hi)
    fs2 = trial.suggest_int("col2_feed_stage", fs2_lo, fs2_hi)
    rr2 = trial.suggest_float("col2_reflux_ratio", *s["col2_reflux_ratio"])

    # Dist3
    p3 = trial.suggest_float("col3_p_kpa", *s["col3_p_kpa"])
    n3 = trial.suggest_int("col3_n_stages", *s["col3_n_stages"])
    fs3_hi = min(s["col3_feed_stage"][1], n3 - 2)
    fs3_lo = min(s["col3_feed_stage"][0], fs3_hi)
    fs3 = trial.suggest_int("col3_feed_stage", fs3_lo, fs3_hi)
    dr3 = trial.suggest_float("col3_draw_rate_kmolh", *s["col3_draw_rate_kmolh"])

    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=T_in_K, z_cat=z_cat_m, t_cyc=t_cyc_min, D=D_reactor_m),
        psa=PSADesignVars(D_col=D_psa_col_m, L_bed=L_psa_bed_m,
                          desorption_target=desorption_target),
        mem=MemDesignVars(P_H=P_H_Pa, P_L=P_L_Pa, A_mem=A_mem_m2,
                          P_dist=p3 * 1000.0),
        dist1=ColumnTunables(
            P_col=p1 * 1000.0, N_stages=n1, N_feed=1, reflux_ratio=2.0,
            solver_method='hysys',
            hysys_spec_value=cf1, hysys_feed_stage=fs1,
        ),
        dist2=ColumnTunables(
            P_col=p2 * 1000.0, N_stages=n2, N_feed=1, reflux_ratio=rr2,
            solver_method='hysys',
            hysys_spec_value=rr2, hysys_feed_stage=fs2,
        ),
        dist3=ColumnTunables(
            P_col=p3 * 1000.0, N_stages=n3, N_feed=1, reflux_ratio=12.0,
            solver_method='hysys',
            # adapter で 0-1 範囲なら recovery_LK_top (動的 Draw Rate)、1.0 超なら絶対量。
            hysys_spec_value=(dr3 if 0 < dr3 <= 1.0 else dr3 / 3600.0),
            hysys_feed_stage=fs3,
        ),
    )


_CONFIG = None


def _get_config():
    global _CONFIG
    if _CONFIG is None:
        _CONFIG = load_operating_config()
    return _CONFIG


def objective(trial: optuna.trial.Trial) -> float:
    """1 trial の評価 — effective_TAC [億円/年] を返す。

    HYSYS 非収束 / フローシート失敗時は PENALTY_VALUE_OKUYEN を返し、BO 探索から
    自動的に外れる方向に学習させる。
    """
    design = _build_design_from_trial(trial)
    t0 = time.perf_counter()
    try:
        result = evaluate(
            design, _get_config(), verbose=False,
            apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
            apply_stage2=APPLY_STAGE2,
            F_C3H8_override=F_C3H8_fresh_kmol_h,
        )
    except Exception as e:
        trial.set_user_attr("failure_reason", f"evaluate 例外: {type(e).__name__}: {e}")
        trial.set_user_attr("wallclock_sec", time.perf_counter() - t0)
        return PENALTY_VALUE_OKUYEN

    trial.set_user_attr("wallclock_sec", time.perf_counter() - t0)
    trial.set_user_attr("failure_reason", result.failure_reason)
    trial.set_user_attr("is_feasible", result.is_feasible)

    tac = float(result.effective_TAC)
    if not math.isfinite(tac):
        return PENALTY_VALUE_OKUYEN
    return min(tac, PENALTY_VALUE_OKUYEN)


# ===========================================================================
#  Optuna study
# ===========================================================================

def build_study() -> optuna.Study:
    sampler = optuna.samplers.TPESampler(
        seed=SAMPLER_SEED,
        n_startup_trials=max(5, N_TRIALS // 4),
        multivariate=True,
    )
    storage: Optional[str] = None
    if USE_SQLITE_STORAGE:
        storage = f"sqlite:///{SQLITE_PATH}"
    return optuna.create_study(
        study_name=STUDY_NAME,
        sampler=sampler,
        storage=storage,
        load_if_exists=USE_SQLITE_STORAGE,
        direction="minimize",
    )


def main():
    print(f"==== special.py: PDH HYSYS BO ({N_TRIALS} trials) ====")
    study = build_study()
    study.optimize(objective, n_trials=N_TRIALS, gc_after_trial=True,
                   show_progress_bar=False)

    # 結果サマリ
    print("\n==== 結果 ====")
    feasible = [t for t in study.trials
                if t.state == optuna.trial.TrialState.COMPLETE
                and t.value is not None
                and t.value < PENALTY_VALUE_OKUYEN]
    print(f"完了 trial 数: {len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])}")
    print(f"feasible trial 数: {len(feasible)}")
    if feasible:
        best = min(feasible, key=lambda t: t.value)
        print(f"best TAC = {best.value:.2f} 億円/年 (trial #{best.number})")
        for k, v in best.params.items():
            print(f"  {k} = {v}")

    # 終了時に HYSYS default provider を確実にクローズ
    try:
        from units.vle.hysys.provider import shutdown_default_provider
        shutdown_default_provider()
    except Exception:
        pass


if __name__ == "__main__":
    main()
