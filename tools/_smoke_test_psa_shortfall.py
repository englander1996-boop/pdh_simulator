"""5 trial の mini-BO を回して、psa_*_shortfall が trial.user_attrs まで届くか確認する。

実行例:
    py -3 tools\_smoke_test_psa_shortfall.py

設計判断 (2026-05-21): 一回限りの確認スクリプト。
"""
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import optuna

from config.load import load_operating_config
from optimization.objective import make_objective
from optimization.study import _default_constraints_func


SEARCH_SPACE = {
    'T_in_K':            (900.0,  970.0,  'linear', 'float'),
    'z_cat_m':           (15.0,   40.0,   'linear', 'float'),
    't_cyc_min':         (12.0,   25.0,   'linear', 'float'),
    'D_reactor_m':       (7.0,    10.0,   'linear', 'float'),
    'D_psa_col_m':       (2.5,    5.0,    'linear', 'float'),
    'L_psa_bed_m':       (15.0,   30.0,   'linear', 'float'),
    'desorption_target': (0.15,   0.40,   'linear', 'float'),
    'P_H_Pa':            (7.5e5,  9.5e5,  'linear', 'float'),
    'A_mem_m2':          (3.0e4,  3.0e5,  'log',    'float'),
    'P_dist1_Pa':        (12.0e5, 25.0e5, 'linear', 'float'),
    'N_dist1':           (16,     30,     'linear', 'int'),
    'reflux_dist1':      (1.5,    3.0,    'linear', 'float'),
    'P_dist2_Pa':        (5.0e5,  7.0e5,  'linear', 'float'),
    'N_dist2':           (20,     40,     'linear', 'int'),
    'reflux_dist2':      (6.0,    10.0,   'linear', 'float'),
    'P_dist3_Pa':        (15.0e5, 25.0e5, 'linear', 'float'),
    'N_dist3':           (80,    250,    'linear', 'int'),
    'reflux_dist3':      (11.0,   20.0,   'linear', 'float'),
    'F_C3H8_fresh_kmol_h': (1200.0, 1700.0, 'linear', 'float'),
    'rec_LK_top_dist2':  (0.95, 0.999, 'linear', 'float'),
    'rec_HK_bot_dist2':  (0.998, 0.9999, 'linear', 'float'),
}
SOLVER_BO = {'dist1': 'fug', 'dist2': 'rigorous', 'dist3': 'fug'}


def main() -> None:
    warnings.simplefilter("once")

    config = load_operating_config()
    objective = make_objective(
        search_space=SEARCH_SPACE,
        solver_assignment=SOLVER_BO,
        config=config,
        apply_hi=True,
        apply_stage2=False,
        strict_recovery_check=False,
        recovery_tolerance=0.10,
        baseline=None,
        n_trials_total=5,
    )
    # _PhaseSwitchSampler 経由で QMC → TPE 切替を試す
    from optimization.study import make_sampler
    sampler = make_sampler(name='tpe', seed=42, n_startup=3)
    study = optuna.create_study(direction='minimize', sampler=sampler)
    study.optimize(objective, n_trials=5, catch=(Exception,), show_progress_bar=False)

    print("\n" + "=" * 70)
    print("Mini-BO 結果 (5 trial)")
    print("=" * 70)
    for t in study.trials:
        attrs = t.user_attrs
        psa_t = attrs.get('psa_t_abs_shortfall', 0.0)
        psa_u = attrs.get('psa_u_0_shortfall', 0.0)
        psa_f = attrs.get('psa_feed_shortfall', 0.0)
        d1_N  = attrs.get('dist1_N_shortfall', 0.0)
        d2_dT = attrs.get('dist2_dT_shortfall', 0.0)
        fr    = attrs.get('failure_reason', '')[:60]
        print(f"  trial {t.number}: value={t.value:.2f}, "
              f"psa_t={psa_t:.3f} psa_u={psa_u:.3f} psa_f={psa_f:.3f} "
              f"d1_N={d1_N:.3f} d2_dT={d2_dT:.3f}")
        print(f"    failure_reason: {fr}")

        # constraint vector を確認
        c = _default_constraints_func(t)
        print(f"    constraints (len={len(c)}): {[f'{x:.3f}' for x in c]}")


if __name__ == "__main__":
    main()
