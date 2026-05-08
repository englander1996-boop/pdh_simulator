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
    LPG_FEED_JPY_PER_KG, C3H6_PRODUCT_JPY_PER_KG, H2_PRODUCT_JPY_PER_KG,
    HHV_MJ_PER_KMOL,
)
from src.component_data import MW


@dataclass
class Economics:
    """経済計算結果。

    定義 (化工テキスト標準):
      TAC      = CAPEX/年 + OPEX (= utility + 触媒 + 吸着剤 + 原料費)
      Revenue  = C3H6 売上 + H2 売上 + オフガス燃料クレジット (全て正)
      Profit   = Revenue - TAC   (正なら黒字)

    最適化器は effective_TAC = TAC - Revenue + soft_penalty を最小化する
    (利益最大化と等価)。runner.py 側で計算。
    """
    capex:          dict    # [億円]
    opex:           dict    # [億円/年]   utility + 触媒 + 吸着剤 + 原料費 (全て正)
    revenue:        dict    # [億円/年]   売上 + 燃料クレジット (全て正)
    total_capex:    float   # [億円]
    total_opex:     float   # [億円/年]
    total_revenue:  float   # [億円/年]
    TAC:            float   # [億円/年]   = total_capex/DEPR + total_opex
    profit:         float   # [億円/年]   = total_revenue - TAC (正=黒字)
    annual_kg_C3H6: float   # [kg/年]
    unit_jpy_per_t: float   # [円/ton]   TAC ベース (製造原単位)


def _ele(W_kW: float) -> float:
    """電力 [kW] → [億円/年]"""
    return W_kW * ELECTRICITY_JPY_PER_KWH * OPERATING_HOURS_PER_YEAR / 1.0e8


def _heat(Q_kW: float, jpy_per_GJ: float) -> float:
    """熱量 [kW] → [億円/年]  (1 kW × 1 h = 3.6 MJ = 3.6e-3 GJ)"""
    return Q_kW * 3.6e-3 * OPERATING_HOURS_PER_YEAR * jpy_per_GJ / 1.0e8


def _annual_okuyen(F_kmol_h: float, MW_kg_per_kmol: float, jpy_per_kg: float) -> float:
    """流量 [kmol/h] × MW × 単価 → [億円/年]。"""
    return (F_kmol_h * MW_kg_per_kmol
            * OPERATING_HOURS_PER_YEAR * jpy_per_kg / 1.0e8)


def _offgas_GJ_per_h(offgas: dict) -> float:
    """PSA オフガス組成 (dict {A,B,...: kmol/h}) を高位発熱量ベースで GJ/h 換算。"""
    return sum(
        offgas.get(c, 0.0) * HHV_MJ_PER_KMOL.get(c, 0.0) / 1000.0
        for c in offgas
    )


def collect_capex_opex(one_pass: dict) -> tuple[dict, dict, dict]:
    """run_one_pass の戻り値から CAPEX・OPEX・Revenue 内訳を抽出。

    Returns
    -------
    capex   : 装置別 CAPEX [億円]
    opex    : OPEX (utility + 触媒 + 吸着剤 + 原料費) [億円/年]   全て正
    revenue : 売上 (C3H6 + H2) + 燃料クレジット (オフガス) [億円/年]   全て正
    """
    R = one_pass

    capex = {
        'Pump1':       R['pump1'].equipment.CAPEX,
        'Dist1':       R['r1'].equipment.CAPEX,
        'Reactor':     R['r_rx'].equipment.Reactor_CAPEX,
        'Cooler':      R['cooled'].equipment.CAPEX,
        'Comp2a':      R['comp2a'].equipment.CAPEX,
        'Intercool':   R['intercool'].equipment.CAPEX,
        'Comp2b':      R['comp2b'].equipment.CAPEX,
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
        'Pump1電力':         _ele(R['pump1'].equipment.W_kW),
        'Comp2a電力':        _ele(R['comp2a'].equipment.W_kW),
        'Comp2b電力':        _ele(R['comp2b'].equipment.W_kW),
        'MemF圧縮機電力':    _ele(R['r_mem'].equipment.W_feed_kW),
        'MemP圧縮機電力':    _ele(R['r_mem'].equipment.W_prod_kW),
        'Dist1リボイラ蒸気': _heat(R['r1'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Dist2リボイラ蒸気': _heat(R['r2'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Dist3リボイラ蒸気': _heat(R['r3'].equipment.Q_reb, LP_STEAM_JPY_PER_GJ),
        'Mem気化器蒸気':     _heat(R['r_mem'].equipment.Q_vap_kW, LP_STEAM_JPY_PER_GJ),
        # 設計判断 (2026-05-08): 蒸留塔フィード予熱を独立計上 (旧版は抜けていた)。
        # distillation_core.py の DistEquipment.Q_feed_preheat_kW から読み取る。
        'Dist1フィード予熱蒸気': _heat(R['r1'].equipment.Q_feed_preheat_kW, LP_STEAM_JPY_PER_GJ),
        'Dist2フィード予熱蒸気': _heat(R['r2'].equipment.Q_feed_preheat_kW, LP_STEAM_JPY_PER_GJ),
        'Dist3フィード予熱蒸気': _heat(R['r3'].equipment.Q_feed_preheat_kW, LP_STEAM_JPY_PER_GJ),
    }

    # 反応器プリヒーター (GJ/h → kW 換算)
    Q_preheat_kW = R['r_rx'].effluent.Q_preheat * 1.0e9 / 3600.0 / 1000.0
    opex['Reactor予熱燃料'] = _heat(Q_preheat_kW, FUEL_JPY_PER_GJ)

    # 設計判断 (2026-05-08): cooler.py が utility_selector で選んだユーティリティ名・
    # 単価を equipment に格納するようになったため、それを直接使う。
    # ハードコードの COOLING_WATER_JPY_PER_GJ はここでは使わない。
    cooled_eq      = R['cooled'].equipment
    intercool_eq   = R['intercool'].equipment
    mem_precool_eq = R['mem_precool'].equipment
    opex[f"Cooler({cooled_eq.utility_name})"] = _heat(
        abs(cooled_eq.Q_duty_kW), cooled_eq.utility_jpy_per_GJ)
    opex[f"Intercool({intercool_eq.utility_name})"] = _heat(
        abs(intercool_eq.Q_duty_kW), intercool_eq.utility_jpy_per_GJ)
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

    # ---- 原料費を OPEX に追加 (TAC に含める標準的な定義) ----
    # 設計判断 (2026-05-09): TAC には化工標準で原料費を含める (Sinnott §6.5,
    # Turton §8.2)。utility と並べて opex dict に置く。
    fresh_F = R['pump1'].outlet.F_in        # Pump1 inlet/outlet 同じ組成
    F_C3H8_feed  = fresh_F.get('A', 0.0)
    F_C4H10_feed = fresh_F.get('Z', 0.0)
    opex['Fresh LPG 原料費'] = (
        _annual_okuyen(F_C3H8_feed,  MW['A'], LPG_FEED_JPY_PER_KG)
        + _annual_okuyen(F_C4H10_feed, MW['Z'], LPG_FEED_JPY_PER_KG)
    )

    # ---- Revenue (売上 + 燃料クレジット、全て正値で計上) ----
    revenue: dict = {}

    F_C3H6_prod = R['r3'].top.F_in.get('B', 0.0)
    revenue['C3H6 製品売上'] = _annual_okuyen(
        F_C3H6_prod, MW['B'], C3H6_PRODUCT_JPY_PER_KG,
    )

    F_H2_prod = R['r_psa'].product.get('C', 0.0)
    revenue['H2 製品売上'] = _annual_okuyen(
        F_H2_prod, MW['C'], H2_PRODUCT_JPY_PER_KG,
    )

    # オフガスは反応器プリヒーター燃料として利用 → equivalent な FUEL_JPY_PER_GJ で
    # 燃料費が浮く。Revenue 側に正値として計上。
    offgas_GJ_per_h = _offgas_GJ_per_h(R['r_psa'].offgas)
    revenue['PSA オフガス燃料クレジット'] = (
        offgas_GJ_per_h * OPERATING_HOURS_PER_YEAR * FUEL_JPY_PER_GJ / 1.0e8
    )

    return capex, opex, revenue


def calculate_economics(one_pass: dict, mw_C3H6_kg_per_kmol: float) -> Economics:
    """CAPEX/OPEX/Revenue 集計と TAC・Profit・製品単価を計算。"""
    capex, opex, revenue = collect_capex_opex(one_pass)

    # ペナルティ装置 (CAPEX >= 1e8 億円相当) は集計から除外
    total_capex   = sum(v for v in capex.values() if v < 1e6)
    total_opex    = sum(opex.values())
    total_revenue = sum(revenue.values())
    TAC           = total_capex / DEPRECIATION_YEARS + total_opex
    profit        = total_revenue - TAC

    # 製品単価 (TAC ベース、製造原単位)
    C3H6_product_kmol_h = one_pass['r3'].top.F_in.get('B', 0.0)
    annual_kg = C3H6_product_kmol_h * OPERATING_HOURS_PER_YEAR * mw_C3H6_kg_per_kmol
    if annual_kg > 0:
        unit_jpy_per_t = TAC * 1.0e8 / (annual_kg / 1000.0)
    else:
        unit_jpy_per_t = float('inf')

    return Economics(
        capex          =capex,
        opex           =opex,
        revenue        =revenue,
        total_capex    =total_capex,
        total_opex     =total_opex,
        total_revenue  =total_revenue,
        TAC            =TAC,
        profit         =profit,
        annual_kg_C3H6 =annual_kg,
        unit_jpy_per_t =unit_jpy_per_t,
    )
