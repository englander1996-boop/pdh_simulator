r"""
comparing.shared.space — 21 設計変数の探索空間と design 構築。

special.py (= 最適化の本丸) の SEARCH_SPACE / _build_design を **移植** したもの。
import で special.py に依存させない (special.py→main.py 改名・削除予定に巻き込まれないため)。
BO も素朴手法も「同じ変数集合・同じ design 構築」を使うことで、手法だけの差を切り出す。

NOTE: special.py を更新したらここも手で同期する (単一の真実をここに置く方針)。
"""

from typing import Optional

# import 順序は special.py と一致させる: flowsheet を先に完全初期化してから
# src.distillation_core を引く。逆順だと distillation_core 初期化途中に
# flowsheet.design が ColumnTunables を引いて循環 import になる。
from flowsheet import FlowsheetDesignVars
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars


# 膜透過側圧力 (大気圧固定、special.py と同じ)
P_L_Pa = 1.0e5

# 探索範囲 — 全 21 変数  形式: (low, high, scale, type)。special.py SEARCH_SPACE と一致。
SEARCH_SPACE: dict = {
    # ----- 反応器 (Swing) -----
    "T_in_K":              (880.0,  940.0,  'linear', 'float'),
    "z_cat_m":             (15.0,   40.0,   'linear', 'float'),
    "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
    "D_reactor_m":         (7.0,    10.0,   'linear', 'float'),
    # ----- PSA -----
    "D_psa_col_m":         (2.9,    5.0,    'linear', 'float'),
    "L_psa_bed_m":         (22.0,   30.0,   'linear', 'float'),
    "desorption_target":   (0.22,   0.40,   'linear', 'float'),
    # ----- 膜 (P_L 固定、P_dist は Dist3 圧と同期) -----
    "P_H_Pa":              (7.5e5,  9.5e5,  'linear', 'float'),
    "A_mem_m2":            (5.0e4,  3.0e5,  'log',    'float'),
    # ----- 原料 (外側ループ skip、override) -----
    "F_C3H8_fresh_kmol_h": (1500.0, 1750.0, 'linear', 'float'),
    # ----- Dist1 (SM) -----
    "col1_p_kpa":          (1600.0, 2000.0, 'linear', 'float'),
    "col1_n_stages":       (30,     60,     'linear', 'int'  ),
    "col1_feed_stage":     (22,     28,     'linear', 'int'  ),
    "col1_comp_frac_2":    (0.90,   0.999,  'linear', 'float'),
    # ----- Dist2 (HYSYS 脱エタン塔) -----
    "col2_p_kpa":          (500.0,  700.0,  'linear', 'float'),
    "col2_n_stages":       (44,     80,     'linear', 'int'  ),
    "col2_feed_ratio":     (0.40,   0.60,   'linear', 'float'),
    "col2_reflux_ratio":   (8.0,    10.5,   'linear', 'float'),
    # ----- Dist3 (SM, spec なし) -----
    "col3_p_kpa":          (1600.0, 1900.0, 'linear', 'float'),
    "col3_n_stages":       (115,    160,     'linear', 'int' ),
    "col3_feed_ratio":     (0.60,   0.90,   'linear', 'float'),
}

# SM/HYSYS の feed_stage 絶対 bounds (ratio→絶対段変換後に clamp)。special.py と一致。
_FEED_STAGE_ABS = {
    "col2": (2, 9999),    # HYSYS は N-2 のみ
    "col3": (70, 180),    # model3 feas 域 (feed_stage >=70)
}

# 蒸留塔バックエンド (column → solver_method)。special.py = Dist2 HYSYS。
# 'fug' に差し替えれば高速スモークテストが可能。
DEFAULT_BACKEND = {"dist1": "sm", "dist2": "hysys", "dist3": "sm"}

PARAM_NAMES = list(SEARCH_SPACE.keys())


def is_int(name: str) -> bool:
    return SEARCH_SPACE[name][3] == 'int'


def bounds(name: str) -> tuple:
    low, high, _scale, _typ = SEARCH_SPACE[name]
    return (low, high)


def clamp(name: str, value: float):
    """value を変数 name の範囲内に収め、int 変数は丸める。"""
    low, high, _scale, typ = SEARCH_SPACE[name]
    v = max(low, min(value, high))
    return int(round(v)) if typ == 'int' else float(v)


def midpoint_params() -> dict:
    """各変数を範囲中央に置いた params dict。素朴手法の「初期エンジニアリング推算」既定値。

    log scale 変数 (A_mem) は幾何中央。int は丸め。
    """
    import math
    p: dict = {}
    for name, (low, high, scale, typ) in SEARCH_SPACE.items():
        if scale == 'log':
            mid = math.sqrt(low * high)
        else:
            mid = 0.5 * (low + high)
        p[name] = int(round(mid)) if typ == 'int' else float(mid)
    return p


def grid_points(name: str, k: int) -> list:
    """変数 name を範囲内で k 点に等間隔サンプリング (int は丸めて重複除去、log は対数等間隔)。"""
    import math
    low, high, scale, typ = SEARCH_SPACE[name]
    if k < 2:
        return [clamp(name, 0.5 * (low + high))]
    if scale == 'log':
        lo, hi = math.log(low), math.log(high)
        raw = [math.exp(lo + (hi - lo) * i / (k - 1)) for i in range(k)]
    else:
        raw = [low + (high - low) * i / (k - 1) for i in range(k)]
    if typ == 'int':
        seen, out = set(), []
        for r in raw:
            iv = int(round(r))
            if iv not in seen:
                seen.add(iv); out.append(iv)
        return out
    return [float(r) for r in raw]


def _feed_stage_from_ratio(ratio: float, n: int, lo: int, hi: int) -> int:
    fs = int(round(ratio * n))
    hi_eff = min(hi, n - 2)
    lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))


def build_design(p: dict, backend: Optional[dict] = None) -> FlowsheetDesignVars:
    """params dict (21 変数) から FlowsheetDesignVars を構築。special.py._build_design と等価。

    backend で各塔の solver_method を差し替え可能 (既定 = special.py 同等 Dist2=HYSYS)。
    """
    be = backend or DEFAULT_BACKEND
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages']); fs2 = _feed_stage_from_ratio(p['col2_feed_ratio'], n2, *_FEED_STAGE_ABS['col2'])
    n3 = int(p['col3_n_stages']); fs3 = _feed_stage_from_ratio(p['col3_feed_ratio'], n3, *_FEED_STAGE_ABS['col3'])
    p3_kpa = float(p['col3_p_kpa'])
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=p['T_in_K'], z_cat=p['z_cat_m'],
                          t_cyc=p['t_cyc_min'], D=p['D_reactor_m']),
        psa=PSADesignVars(D_col=p['D_psa_col_m'], L_bed=p['L_psa_bed_m'],
                          desorption_target=p['desorption_target']),
        mem=MemDesignVars(P_H=p['P_H_Pa'], P_L=P_L_Pa, A_mem=p['A_mem_m2'],
                          P_dist=p3_kpa * 1000.0),
        dist1=ColumnTunables(
            P_col=float(p['col1_p_kpa']) * 1000.0, N_stages=n1, N_feed=1, reflux_ratio=2.0,
            solver_method=be['dist1'], hysys_spec_value=float(p['col1_comp_frac_2']),
            hysys_feed_stage=fs1),
        dist2=ColumnTunables(
            P_col=float(p['col2_p_kpa']) * 1000.0, N_stages=n2, N_feed=1,
            reflux_ratio=float(p['col2_reflux_ratio']),
            solver_method=be['dist2'], hysys_spec_value=float(p['col2_reflux_ratio']),
            hysys_feed_stage=fs2),
        dist3=ColumnTunables(
            P_col=p3_kpa * 1000.0, N_stages=n3, N_feed=1, reflux_ratio=12.0,
            solver_method=be['dist3'], hysys_spec_value=0.99, hysys_feed_stage=fs3),
    )
