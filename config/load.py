"""
operating.toml を型付き dataclass にロードする。

最適化器・ランナーから `load_operating_config()` を呼ぶことで、
「設計者が決め打ちしている運転条件」と「最適化対象の設計変数」を
明確に分離する。
"""

import os
import tomllib
from dataclasses import dataclass


@dataclass(frozen=True)
class ProductSpec:
    target_mta: float
    mw_kg_per_kmol: float


@dataclass(frozen=True)
class FeedSpec:
    lpg_c3h8_mol_fraction: float
    yield_assumed: float
    T_K: float
    P_Pa: float


@dataclass(frozen=True)
class PressureSpec:
    pump1_out_Pa: float       # 原料を液送で 17 bar (Dist1 入口) まで昇圧
    comp2_out_Pa: float
    reactor_inlet_Pa: float


@dataclass(frozen=True)
class TemperatureSpec:
    cooler_after_reactor_K: float
    mem_feed_K: float


@dataclass(frozen=True)
class InnerSolverSpec:
    """内側ソルバ設定。
    判定式: max( |Δtear[k]| / max(|tear[k]|, tol_floor_kmol_h) ) < tol_relative
    method:
      "successive_substitution": 旧来の α-relaxation (relax フィールドのみ使用)
      "wegstein": 履歴2点から最適加速係数を自動計算 (wegstein_q_* も使用)
    """
    max_iter: int
    tol_relative: float        # 相対許容誤差 (例: 0.001 = 0.1%)
    tol_floor_kmol_h: float    # 相対計算 denominator の下限 [kmol/h]
    relax: float               # 旧 SS 法の α / Wegstein 第1反復の α
    recycle_guard_ratio: float
    method: str                # "successive_substitution" or "wegstein"
    wegstein_q_min: float      # wegstein 用: 加速係数下限 (負ほど強加速)
    wegstein_q_max: float      # wegstein 用: 加速係数上限


@dataclass(frozen=True)
class OuterSolverSpec:
    """外側ソルバ設定。
    判定式: max(0, target - actual) < target × tol_relative   (片側、undershoot のみ違反)
    """
    max_iter: int
    tol_relative: float        # 片側相対許容誤差 (例: 0.01 = 1%)
    relax: float


@dataclass(frozen=True)
class InitSpec:
    T_d3_K: float
    T_mem_K: float
    tear_d3_A_per_1500_fresh: float
    tear_d3_B_per_1500_fresh: float
    tear_mem_A_per_1500_fresh: float
    tear_mem_B_per_1500_fresh: float


@dataclass(frozen=True)
class SolverSpec:
    inner: InnerSolverSpec
    outer: OuterSolverSpec
    init: InitSpec


@dataclass(frozen=True)
class SpecSpec:
    """製品仕様の最低限度 (feasibility 制約)。"""
    c3h6_min_wtfrac:         float  # C3H6 製品の質量分率下限
    h2_min_molfrac:          float  # H2 製品のモル分率下限
    production_min_relative: float  # 生産量の許容相対不足 (target × これ未満は違反)


@dataclass(frozen=True)
class PenaltySpec:
    """最適化器が見る effective_TAC への加算ロジック。
    詳細な階層構造は config/operating.toml の [penalty] セクションのコメント参照。
    """
    solver_failure_okuyen: float   # solver-level 失敗時の TAC 固定値 [億円/年]
    spec_base_okuyen:      float   # spec 違反 1 つにつき固定加算 [億円/年]
    spec_coef_okuyen:      float   # 違反 %pt 当たりの加算係数 [億円/(年・%pt)]


@dataclass(frozen=True)
class OperatingConfig:
    product:     ProductSpec
    feed:        FeedSpec
    pressure:    PressureSpec
    temperature: TemperatureSpec
    solver:      SolverSpec
    spec:        SpecSpec
    penalty:     PenaltySpec


_DEFAULT_PATH = os.path.normpath(
    os.path.join(os.path.dirname(__file__), 'operating.toml')
)


def load_operating_config(path: str | None = None) -> OperatingConfig:
    """operating.toml をロードして OperatingConfig を返す。

    Parameters
    ----------
    path : str | None
        TOML ファイルパス。None のときは config/operating.toml を読む。
    """
    p = path or _DEFAULT_PATH
    with open(p, 'rb') as f:
        d = tomllib.load(f)
    return OperatingConfig(
        product    =ProductSpec(**d['product']),
        feed       =FeedSpec(**d['feed']),
        pressure   =PressureSpec(**d['pressure']),
        temperature=TemperatureSpec(**d['temperature']),
        solver     =SolverSpec(
            inner=InnerSolverSpec(**d['solver']['inner']),
            outer=OuterSolverSpec(**d['solver']['outer']),
            init =InitSpec(**d['solver']['init']),
        ),
        spec       =SpecSpec(**d['spec']),
        penalty    =PenaltySpec(**d['penalty']),
    )
