"""
1 パス分の結果 (run_one_pass の戻り値) から CAPEX・OPEX・TAC を集計する。

CAPEX/OPEX 単価は src/cost_parameters.py に集約されているため、
最終単価のチューニングはそちら 1 ファイルで完結する。
"""

from dataclasses import dataclass, field

from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
    CATALYST_PTSN_JPY_PER_KG, CATALYST_PTSN_LIFE_YEARS,
    OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS,
)


@dataclass
class Economics:
    """経済計算結果。最適化器はこの TAC (または unit_jpy_per_t) を目的関数に使う。"""
    capex:          dict   # [億円]
    opex:           dict   # [億円/年]
    total_capex:    float  # [億円]
    total_opex:    float   # [億円/年]
    TAC:            float  # [億円/年]
    annual_kg_C3H6: float  # [kg/年]
    unit_jpy_per_t: float  # [円/ton]


def _ele(W_kW: float) -> float:
    """電力 [kW] → [億円/年]"""
    return W_kW * ELECTRICITY_JPY_PER_KWH * OPERATING_HOURS_PER_YEAR / 1.0e8


def _heat(Q_kW: float, jpy_per_GJ: float) -> float:
    """熱量 [kW] → [億円/年]  (1 kW × 1 h = 3.6 MJ = 3.6e-3 GJ)"""
    return Q_kW * 3.6e-3 * OPERATING_HOURS_PER_YEAR * jpy_per_GJ / 1.0e8


def collect_capex_opex(one_pass: dict) -> tuple[dict, dict]:
    """run_one_pass の戻り値から CAPEX・OPEX 内訳を抽出。"""
    R = one_pass

    capex = {
        'Comp1':       R['comp1'].equipment.CAPEX,
        'Dist1':       R['r1'].equipment.CAPEX,
        'Reactor':     R['r_rx'].equipment.Reactor_CAPEX,
        'Cooler':      R['cooled'].equipment.CAPEX,
        'Comp2':       R['comp2'].equipment.CAPEX,
        'Dist2':       R['r2'].equipment.CAPEX,
        'PSA容器':     R['r_psa'].equipment.CAPEX_vessels,
        'PSA活性炭':   R['r_psa'].equipment.CAPEX_adsorbent,
        'MemPrecool':  R['mem_precool'].equipment.CAPEX,
        'Mem気化器':   R['r_mem'].equipment.CAPEX_vap,
        'Mem F圧縮機': R['r_mem'].equipment.CAPEX_comp_feed,
        'Mem P圧縮機': R['r_mem'].equipment.CAPEX_comp_prod,
        'Mem冷却器':   R['r_mem'].equipment.CAPEX_cond,
        'Mem膜本体':   R['r_mem'].equipment.CAPEX_mem,
        'Dist3':       R['r3'].equipment.CAPEX,
    }

    opex = {
        'Comp1電力':         _ele(R['comp1'].equipment.W_kW),
        'Comp2電力':         _ele(R['comp2'].equipment.W_kW),
        'MemF圧縮機電力':    _ele(R['r_mem'].equipment.W_feed_kW),
        'MemP圧縮機電力':    _ele(R['r_mem'].equipment.W_prod_kW),
        'Dist1リボイラ蒸気': _heat(R['r1'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Dist2リボイラ蒸気': _heat(R['r2'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Dist3リボイラ蒸気': _heat(R['r3'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Mem気化器蒸気':     _heat(R['r_mem'].equipment.Q_vap_kW, LP_STEAM_JPY_PER_GJ),
    }

    # 反応器プリヒーター (GJ/h → kW 換算)
    Q_preheat_kW = R['r_rx'].effluent.Q_preheat * 1.0e9 / 3600.0 / 1000.0
    opex['Reactor予熱燃料'] = _heat(Q_preheat_kW, FUEL_JPY_PER_GJ)

    # 設計判断 (2026-05-08): cooler.py が utility_selector で選んだユーティリティ名・
    # 単価を equipment に格納するようになったため、それを直接使う。
    # ハードコードの COOLING_WATER_JPY_PER_GJ はここでは使わない。
    cooled_eq      = R['cooled'].equipment
    mem_precool_eq = R['mem_precool'].equipment
    opex[f"Cooler({cooled_eq.utility_name})"] = _heat(
        abs(cooled_eq.Q_duty_kW), cooled_eq.utility_jpy_per_GJ)
    opex[f"MemPrecool({mem_precool_eq.utility_name})"] = _heat(
        abs(mem_precool_eq.Q_duty_kW), mem_precool_eq.utility_jpy_per_GJ)

    # 設計判断: 蒸留塔のコンデンサと膜の冷却器は内部実装が COOLING_WATER 前提のまま
    # (fake_columnX, membrane_system.py)。本実装フェーズで utility_selector に
    # 統合する予定だが、現状は冷却水単独で計上 (旧版踏襲)。
    opex['Dist1コンデンサ冷水'] = _heat(R['r1'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['Dist2コンデンサ冷水'] = _heat(R['r2'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['Dist3コンデンサ冷水'] = _heat(R['r3'].equipment.Q_cond, COOLING_WATER_JPY_PER_GJ)
    opex['Mem冷却器冷水']      = _heat(R['r_mem'].equipment.Q_cond_kW,
                                        COOLING_WATER_JPY_PER_GJ)

    opex['Reactor触媒交換'] = (R['r_rx'].equipment.Catalyst_Weight_Total
                                * CATALYST_PTSN_JPY_PER_KG
                                / CATALYST_PTSN_LIFE_YEARS / 1.0e8)
    opex['PSA活性炭交換']   = R['r_psa'].equipment.OPEX_adsorbent_okuyen_per_year

    return capex, opex


def calculate_economics(one_pass: dict, mw_C3H6_kg_per_kmol: float) -> Economics:
    """CAPEX/OPEX 集計と TAC・製品単価を計算。"""
    capex, opex = collect_capex_opex(one_pass)

    # ペナルティ装置 (CAPEX >= 1e8 億円相当) は集計から除外
    total_capex = sum(v for v in capex.values() if v < 1e6)
    total_opex  = sum(opex.values())
    TAC         = total_capex / DEPRECIATION_YEARS + total_opex

    # 製品単価
    C3H6_product_kmol_h = one_pass['r3'].top.F_in.get('B', 0.0)
    annual_kg = C3H6_product_kmol_h * OPERATING_HOURS_PER_YEAR * mw_C3H6_kg_per_kmol
    unit_jpy_per_t = (TAC * 1.0e8 / (annual_kg / 1000.0)) if annual_kg > 0 else float('inf')

    return Economics(
        capex          =capex,
        opex           =opex,
        total_capex    =total_capex,
        total_opex     =total_opex,
        TAC            =TAC,
        annual_kg_C3H6 =annual_kg,
        unit_jpy_per_t =unit_jpy_per_t,
    )
