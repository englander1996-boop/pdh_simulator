"""
製品仕様 (純度・生産量) の compliance チェック。

設計判断 (2026-05-08, 相談時の合意):
  - 「収束反復ノイズ (tear stream の前後反復差)」と「製品品質 (純度)」を
    別の判定基準として完全に分離する。
  - solver は流量の収束だけ見て、品質はここで独立に評価する。
  - 「収束しているけど純度が出ていない」状態を最適化器が認識できるようにする。

判定対象:
  (1) C3H6 純度 (Dist3 塔頂の質量分率) — contest.md §2-1 で 99.5 wt% 必須
  (2) H2 純度 (PSA 製品のモル分率)   — H2 を製品販売する想定で独自に 99.9 mol%
  (3) 生産量 (Dist3 塔頂の C3H6 流量) — target × (1 - tol) 以上で OK (片側)

オフガス (CH4 等) は contest 規定通り純度制約なし。
contest §3-3 の C2 留分 99.8 mol% は本フローでは該当ストリーム不在のため対象外。
"""

from dataclasses import dataclass

from src.cost_parameters import OPERATING_HOURS_PER_YEAR
from config.load import OperatingConfig


# 成分分子量 [kg/kmol]
# 設計判断: contest.md §3-2 の反応式に登場する成分の標準分子量。
# C3H8/C3H6/H2/C2H4/CH4/C2H6/C4H10 を A〜F, Z にマッピング (ProcessStream 仕様)。
_MW = {
    'A': 44.10,   # C3H8
    'B': 42.08,   # C3H6
    'C':  2.02,   # H2
    'D': 28.05,   # C2H4
    'E': 16.04,   # CH4
    'F': 30.07,   # C2H6
    'Z': 58.12,   # C4H10
}


@dataclass
class SpecComplianceResult:
    """各 spec の判定結果と数値。最適化器のペナルティ計算で使う。

    violation_pp_* は spec を %ポイント単位で正規化した「不足量」。
    (例: spec 99.5%, actual 99.0% → violation_pp = 0.5)
    生産量だけは元々比率なので「target に対する相対不足 × 100」で揃える。
    全制約を %pt スケールに揃えることで、ペナルティ係数を1つで管理可能。
    """
    c3h6_purity_wtfrac:  float
    h2_purity_molfrac:   float
    production_kmol_h:   float
    target_kmol_h:       float

    c3h6_pass:           bool
    h2_pass:             bool
    production_pass:     bool

    c3h6_violation_pp:       float   # %pt (>0 のとき不足)
    h2_violation_pp:         float   # %pt
    production_violation_pp: float   # %pt (target に対する不足率 × 100)

    @property
    def all_pass(self) -> bool:
        return self.c3h6_pass and self.h2_pass and self.production_pass


def _target_kmol_h(config: OperatingConfig) -> float:
    """contest.md §2-1 の年産から時間当たり kmol/h に換算。"""
    return (config.product.target_mta * 1000.0
            / config.product.mw_kg_per_kmol
            / OPERATING_HOURS_PER_YEAR)


def check_specs(one_pass: dict, config: OperatingConfig) -> SpecComplianceResult:
    """run_one_pass の戻り値から純度・生産量を計算して spec 判定。"""
    spec = config.spec

    # ---- (1) C3H6 純度 (Dist3 塔頂、wt%) ----
    top_F = one_pass['r3'].top.F_in
    mass_total = sum(top_F.get(k, 0.0) * _MW[k] for k in _MW)
    mass_c3h6  = top_F.get('B', 0.0) * _MW['B']
    c3h6_wt    = mass_c3h6 / mass_total if mass_total > 0 else 0.0

    # ---- (2) H2 純度 (PSA 製品、mol%) ----
    psa_prod      = one_pass['r_psa'].product
    mol_total_psa = sum(psa_prod.values())
    mol_h2        = psa_prod.get('C', 0.0)
    h2_mol        = mol_h2 / mol_total_psa if mol_total_psa > 0 else 0.0

    # ---- (3) 生産量 (Dist3 塔頂の C3H6 流量) ----
    production_kmol_h = top_F.get('B', 0.0)
    target_kmol_h     = _target_kmol_h(config)

    # ---- 判定 (片側: 規格値以上で pass) ----
    c3h6_pass       = c3h6_wt >= spec.c3h6_min_wtfrac
    h2_pass         = h2_mol  >= spec.h2_min_molfrac
    production_pass = production_kmol_h >= target_kmol_h * (1.0 - spec.production_min_relative)

    # ---- 違反量を %pt 単位で正規化 ----
    # 設計判断: 異なる種類の制約 (質量分率 / モル分率 / 流量比率) を %pt スケールに
    # 揃えることで、ペナルティ係数 spec_coef_okuyen を1つで管理できるようにする。
    c3h6_violation_pp = max(0.0, (spec.c3h6_min_wtfrac - c3h6_wt) * 100.0)
    h2_violation_pp   = max(0.0, (spec.h2_min_molfrac  - h2_mol)  * 100.0)
    if production_pass:
        production_violation_pp = 0.0
    else:
        # 不足率 = (target × (1-tol) - actual) / target × 100  → %pt
        threshold = target_kmol_h * (1.0 - spec.production_min_relative)
        production_violation_pp = (threshold - production_kmol_h) / target_kmol_h * 100.0

    return SpecComplianceResult(
        c3h6_purity_wtfrac     =c3h6_wt,
        h2_purity_molfrac      =h2_mol,
        production_kmol_h      =production_kmol_h,
        target_kmol_h          =target_kmol_h,
        c3h6_pass              =c3h6_pass,
        h2_pass                =h2_pass,
        production_pass        =production_pass,
        c3h6_violation_pp      =c3h6_violation_pp,
        h2_violation_pp        =h2_violation_pp,
        production_violation_pp=production_violation_pp,
    )
