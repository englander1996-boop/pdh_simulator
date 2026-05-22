"""
optimization/search_space.py — 探索空間スキーマと params ↔ FlowsheetDesignVars 変換

main.py の SEARCH_SPACE 辞書を解釈し、Optuna trial から params を suggest する。
suggest された params (+ 欠落キーの baseline 補完) から FlowsheetDesignVars を構築する。

設計判断:
  - SEARCH_SPACE のキーをコメントアウトすれば、その変数は baseline 固定となる。
    変数の追加・削除は main.py のみ編集すれば完結する (本モジュール変更不要)。
  - N_feed (蒸留塔フィード段) は探索対象外。core 側で Kirkbride 推奨を自動採用するため、
    ColumnTunables.N_feed には _N_FEED_PLACEHOLDER (=1) を渡す (実値は無視される)。
  - P_L (膜透過側圧力) は大気圧 1 atm に固定 (透過側真空ポンプ無し前提)。
  - P_dist3 と mem.P_dist は同期 (1 変数として扱う)。
"""

from typing import Dict, Tuple, Any, Optional

from flowsheet import FlowsheetDesignVars
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars


# (low, high, scale, type) のタプル形式
#   scale: 'linear' | 'log'
#   type:  'float'  | 'int'
VarSpec = Tuple[float, float, str, str]


# 25 変数の標準キー一覧 (SEARCH_SPACE で許容されるキー)
# 内訳:
#   18 (元): 反応器 4 + PSA 3 + 膜 2 + Dist1/2/3 各 3 (P/N/R)
#    7 (新, 2026-05-14): F_fresh_C3H8 + 各塔の recovery_{LK_top, HK_bot} × 3 塔
EXPECTED_KEYS = (
    # 反応器
    'T_in_K', 'z_cat_m', 't_cyc_min', 'D_reactor_m',
    # PSA
    'D_psa_col_m', 'L_psa_bed_m', 'desorption_target',
    # 膜
    'P_H_Pa', 'A_mem_m2',
    # Dist1
    'P_dist1_Pa', 'N_dist1', 'reflux_dist1',
    # Dist2
    'P_dist2_Pa', 'N_dist2', 'reflux_dist2',
    # Dist3
    'P_dist3_Pa', 'N_dist3', 'reflux_dist3',
    # F_fresh (BO で直接指定、外側ループ skip)
    'F_C3H8_fresh_kmol_h',
    # 塔別 recovery (旧 0.99 ハードコードを変数化)
    'rec_LK_top_dist1', 'rec_HK_bot_dist1',
    'rec_LK_top_dist2', 'rec_HK_bot_dist2',
    'rec_LK_top_dist3', 'rec_HK_bot_dist3',
)


# 膜透過側圧力 (固定、大気圧 1 atm)
P_L_FIXED_PA = 1.0e5

# ColumnTunables.N_feed のプレースホルダ。
# rigorous/sm では core 側で Kirkbride 推奨を自動採用、FUG では post-hoc 出力にのみ使う。
_N_FEED_PLACEHOLDER = 1


# ---------------------------------------------------------------------------
# 依存サンプリング (2026-05-22 E-plan)
# ---------------------------------------------------------------------------
# 設計判断: 物理的に必須の不等式 (例: P_H > feed.P_in = design.dist2.P_col) を
# search bounds 越しに保証する仕組み。Mem の silent penalty 経路 ph_le_pfeed を
# 物理的に発生不可にできる。
#
# 実装方式: post-clip
#   - Optuna の suggest_float は静的 bounds (low, high) を持ち、TPE の KDE/Sobol の
#     基準としてそのまま使用 (分布変更警告を避ける)。
#   - suggest 直後に依存値以下なら静かに max(suggested, dep_low) でクリップ。
#   - TPE 視点では「同じ bounds で動き続ける」ため学習が安定。クリップは
#     constraint_func の代替ではないが、ハード fail を確実に防ぐ最後の門番。
#
# 現状の SEARCH_SPACE (P_H_min=7.5e5, P_dist2_max=7.0e5, margin=0.5e5) では
# 既に P_H ≥ P_dist2 + margin が常に成立しているため、本機構は実質 no-op。
# 将来 P_dist2 上限を 8e5 等に広げたり P_H 下限を 6e5 に下げたら自動で発火する保険。
#
# margin の根拠: フィード圧縮機が有意な仕事をする最低圧縮比 ~1.07 を
# P_dist2=7bar の場合に 0.5 bar = 50 kPa として確保。
_DEPENDENT_MIN_MARGIN_PA = {
    # key: P_H_Pa の依存先 'P_dist2_Pa', 必要マージン 0.5 bar
    'P_H_Pa': ('P_dist2_Pa', 5.0e4),
}


def _apply_dependent_lows(params: Dict[str, Any]) -> Dict[str, Any]:
    """suggest 後の params に対し、依存サンプリングのクリップを適用する。

    例: P_H_Pa が P_dist2_Pa + 0.5bar を下回っていれば持ち上げる。
    依存先が params に無い (SEARCH_SPACE から除外されている = baseline 固定値) 場合は
    DEFAULT_BASELINE を参照、それも無ければ 0 として扱う (実質 no-op)。
    """
    for target_key, (dep_key, margin_pa) in _DEPENDENT_MIN_MARGIN_PA.items():
        if target_key not in params:
            continue
        dep_value = params.get(dep_key)
        if dep_value is None:
            dep_value = DEFAULT_BASELINE.get(dep_key, 0.0)
        floor = float(dep_value) + float(margin_pa)
        if params[target_key] < floor:
            params[target_key] = floor
    return params


# exp1 baseline 値 (suggest 対象外のキーを補完するため)
DEFAULT_BASELINE: Dict[str, Any] = {
    # 反応器
    'T_in_K':            950.0,
    'z_cat_m':           30.0,
    't_cyc_min':         15.0,
    'D_reactor_m':       7.0,
    # PSA
    'D_psa_col_m':       3.0,
    'L_psa_bed_m':       20.0,
    'desorption_target': 0.35,
    # 膜
    'P_H_Pa':            9.5e5,
    'A_mem_m2':          1.0e5,
    # Dist1
    'P_dist1_Pa':        17.0e5,
    'N_dist1':           20,
    'reflux_dist1':      1.5,
    # Dist2
    'P_dist2_Pa':        8.5e5,
    'N_dist2':           20,
    'reflux_dist2':      6.0,
    # Dist3
    'P_dist3_Pa':        20.0e5,
    'N_dist3':           200,
    'reflux_dist3':      12.0,
    # F_fresh: None → 外側ループに任せる (exp1 baseline 互換)。
    #          float → BO が直接指定、外側ループ skip。
    'F_C3H8_fresh_kmol_h': None,
    # 各塔 recovery: None → ラッパ既定値 0.99 (exp1 baseline 互換)、float → BO で振る
    'rec_LK_top_dist1':  None,
    'rec_HK_bot_dist1':  None,
    'rec_LK_top_dist2':  None,
    'rec_HK_bot_dist2':  None,
    'rec_LK_top_dist3':  None,
    'rec_HK_bot_dist3':  None,
}


def validate_search_space(search_space: Dict[str, VarSpec]) -> None:
    """SEARCH_SPACE 辞書のキー・値が妥当か検査。

    main.py 起動時に呼んで、誤ったキー・型・範囲を早期に弾く。
    """
    extras = set(search_space) - set(EXPECTED_KEYS)
    if extras:
        raise ValueError(
            f"SEARCH_SPACE に未知の変数: {sorted(extras)}\n"
            f"許容キーは {EXPECTED_KEYS}"
        )
    for name, spec in search_space.items():
        if not isinstance(spec, tuple) or len(spec) != 4:
            raise ValueError(
                f"SEARCH_SPACE['{name}']: spec は (low, high, scale, type) の 4-tuple 必要。"
                f"\n受け取り: {spec!r}"
            )
        low, high, scale, vtype = spec
        if not (isinstance(low, (int, float)) and isinstance(high, (int, float))):
            raise ValueError(f"SEARCH_SPACE['{name}']: low/high は数値必要 ({low!r}, {high!r})")
        if low >= high:
            raise ValueError(f"SEARCH_SPACE['{name}']: low ({low}) >= high ({high})")
        if scale not in ('linear', 'log'):
            raise ValueError(f"SEARCH_SPACE['{name}']: scale は 'linear' | 'log' (受け取り: {scale!r})")
        if vtype not in ('float', 'int'):
            raise ValueError(f"SEARCH_SPACE['{name}']: type は 'float' | 'int' (受け取り: {vtype!r})")
        if scale == 'log' and low <= 0:
            raise ValueError(f"SEARCH_SPACE['{name}']: log scale には low > 0 が必要 (low={low})")


def suggest_params(trial, search_space: Dict[str, VarSpec]) -> Dict[str, Any]:
    """Optuna trial から SEARCH_SPACE で定義された params を suggest する。

    Parameters
    ----------
    trial : optuna.Trial
        Optuna の trial オブジェクト。
    search_space : dict
        main.py の SEARCH_SPACE 辞書。

    Returns
    -------
    dict
        変数名 → suggest 値の辞書。SEARCH_SPACE に無い変数は含まれない。
        suggest 直後に `_apply_dependent_lows` で物理依存の下限クリップを適用する
        (詳細は _DEPENDENT_MIN_MARGIN_PA のコメント参照)。
    """
    params: Dict[str, Any] = {}
    for name, spec in search_space.items():
        low, high, scale, vtype = spec
        if vtype == 'int':
            params[name] = trial.suggest_int(
                name, int(low), int(high), log=(scale == 'log'),
            )
        elif vtype == 'float':
            params[name] = trial.suggest_float(
                name, float(low), float(high), log=(scale == 'log'),
            )
        else:
            raise ValueError(f"未知の type {vtype!r} (許容: 'int' | 'float')")
    # 依存サンプリングのクリップ (2026-05-22 E-plan)
    _apply_dependent_lows(params)
    return params


def build_design(
    params:            Dict[str, Any],
    solver_assignment: Dict[str, str],
    baseline:          Dict[str, Any] | None = None,
) -> FlowsheetDesignVars:
    """params 辞書から FlowsheetDesignVars を組み立てる。

    Parameters
    ----------
    params : dict
        suggest_params() で suggest された設計変数値、または手動指定値。
        SEARCH_SPACE で suggest 対象外だったキーは baseline で補完される。
    solver_assignment : dict
        {'dist1': 'fug'|'rigorous'|'sm', 'dist2': ..., 'dist3': ...}
        BO ループでは SOLVER_BO、top-k 再評価では SOLVER_TOPK を渡す想定。
    baseline : dict | None
        suggest 対象外の変数のデフォルト値。None なら DEFAULT_BASELINE を使用。
    """
    if baseline is None:
        baseline = DEFAULT_BASELINE
    p = {**baseline, **params}
    P_dist3 = float(p['P_dist3_Pa'])

    def _opt_float(key: str):
        """None なら None を返す (= ラッパ既定値採用)、それ以外は float 化。"""
        v = p.get(key)
        return None if v is None else float(v)

    return FlowsheetDesignVars(
        swing=SwingDesign(
            T_in=float(p['T_in_K']),
            z_cat=float(p['z_cat_m']),
            t_cyc=float(p['t_cyc_min']),
            D=float(p['D_reactor_m']),
        ),
        psa=PSADesignVars(
            D_col=float(p['D_psa_col_m']),
            L_bed=float(p['L_psa_bed_m']),
            desorption_target=float(p['desorption_target']),
        ),
        mem=MemDesignVars(
            P_H=float(p['P_H_Pa']),
            P_L=P_L_FIXED_PA,
            A_mem=float(p['A_mem_m2']),
            P_dist=P_dist3,                     # Dist3 と同期
        ),
        dist1=ColumnTunables(
            P_col=float(p['P_dist1_Pa']),
            N_stages=int(p['N_dist1']),
            N_feed=_N_FEED_PLACEHOLDER,        # Kirkbride 自動採用 (core 側)
            reflux_ratio=float(p['reflux_dist1']),
            solver_method=solver_assignment['dist1'],
            recovery_LK_top=_opt_float('rec_LK_top_dist1'),
            recovery_HK_bot=_opt_float('rec_HK_bot_dist1'),
        ),
        dist2=ColumnTunables(
            P_col=float(p['P_dist2_Pa']),
            N_stages=int(p['N_dist2']),
            N_feed=_N_FEED_PLACEHOLDER,
            reflux_ratio=float(p['reflux_dist2']),
            solver_method=solver_assignment['dist2'],
            recovery_LK_top=_opt_float('rec_LK_top_dist2'),
            recovery_HK_bot=_opt_float('rec_HK_bot_dist2'),
        ),
        dist3=ColumnTunables(
            P_col=P_dist3,
            N_stages=int(p['N_dist3']),
            N_feed=_N_FEED_PLACEHOLDER,
            reflux_ratio=float(p['reflux_dist3']),
            solver_method=solver_assignment['dist3'],
            recovery_LK_top=_opt_float('rec_LK_top_dist3'),
            recovery_HK_bot=_opt_float('rec_HK_bot_dist3'),
        ),
    )


def extract_F_fresh_override(
    params:    Dict[str, Any],
    baseline:  Dict[str, Any] | None = None,
) -> Optional[float]:
    """params から F_C3H8_fresh_kmol_h を取り出して evaluate に渡す形に整える。

    None なら外側ループに任せる (= 従来動作)、float なら外側ループ skip。
    """
    if baseline is None:
        baseline = DEFAULT_BASELINE
    v = params.get('F_C3H8_fresh_kmol_h', baseline.get('F_C3H8_fresh_kmol_h'))
    return None if v is None else float(v)
