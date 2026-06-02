"""
1 パス分の結果 (run_one_pass の戻り値) から CAPEX・OPEX・TAC を集計する。

CAPEX/OPEX 単価は src/cost_parameters.py に集約されているため、
最終単価のチューニングはそちら 1 ファイルで完結する。

OPEX 集計は Hasebe (プロセス設計R08-3.pdf) §3.4 式 (10) の経験式に準拠:
    製造コスト = 0.180·C_TM + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM) + 減価償却費
本実装では「減価償却費」を CAPEX/DEPRECIATION_YEARS として TAC 側で別途加算し、
opex dict には Hasebe 式 (10) のうち減価償却費以外の項を反映する:
    total_opex = 0.180·C_TM + 2.73·C_OL + 1.23·C_UT + 1.23·C_WT + 1.23·C_RM
                 + 触媒・吸着剤交換費 (Hasebe 式外、別建て)
    TAC        = CAPEX / DEPRECIATION_YEARS + total_opex
opex dict は装置別の生 (1.00 倍) 値に加え、Hasebe 集計項 (`[Hasebe-集計]` プレフィクス)
を末尾に持つ。これにより `sum(opex.values()) == total_opex` の不変条件が成り立つ。
"""

import math
import os
from dataclasses import dataclass, field

from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ, FURNACE_EFFICIENCY,
    CATALYST_JPY_PER_KG, CATALYST_LIFE_YEARS,
    MEM_LIFETIME_YEARS,
    OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS,
    LPG_C3H8_JPY_PER_KG, LPG_C4H10_JPY_PER_KG,
    C3H6_PRODUCT_JPY_PER_KG, H2_PRODUCT_JPY_PER_KG,
    HHV_MJ_PER_KMOL,
    HASEBE_NOL_COEFF, HASEBE_SHIFT_MULTIPLIER, OPERATOR_ANNUAL_SALARY_JPY,
    HASEBE_COEFF_C_TM, HASEBE_COEFF_C_OL, HASEBE_COEFF_C_UT_WT_RM,
    HASEBE_C_WT_OKUYEN_PER_YEAR,
    PENALTY_CAPEX_THRESHOLD_OKUYEN,
)
from src.component_data import MW

# opex dict 内で Hasebe 集計項を識別するためのキー prefix。
# HI/HEN 等の下流処理はこの prefix を見て「再計算が必要な集計項」だと判断する。
HASEBE_AGGR_PREFIX = '[Hasebe-集計] '


@dataclass
class Economics:
    """経済計算結果。

    定義 (Hasebe §3.4 式 (10) 準拠):
      TAC      = CAPEX/年 + OPEX
      OPEX     = 0.180·C_TM + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM)
                 + 触媒・吸着剤交換 (Hasebe 式外、別建て)
                 ※ 減価償却費は CAPEX/年 として上で別途加算済
      Revenue  = C3H6 売上 + H2 売上 + オフガス燃料クレジット (全て正)
      Profit   = Revenue - TAC   (正なら黒字)

    opex dict は装置別の生エントリと `[Hasebe-集計]` プレフィクス付きの集計項
    (0.180·C_TM, 2.73·C_OL, 0.23·C_UT delta, 0.23·C_RM delta) を持つ。
    `sum(opex.values()) == total_opex` の不変条件が成り立つ。

    最適化器は effective_TAC を最小化する。runner.py 側で計算。
    """
    capex:          dict    # [億円]
    opex:           dict    # [億円/年]   装置別 + Hasebe 集計項 (全て正)
    revenue:        dict    # [億円/年]   売上 + 燃料クレジット (全て正)
    total_capex:    float   # [億円]
    total_opex:     float   # [億円/年]   sum(opex.values())
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
    """ガス組成 (dict {A,B,...: kmol/h}) を高位発熱量ベースで GJ/h 換算。

    PSA オフガス・Dist1 塔底 (C4H10 主) などの「廃棄→燃料 CR」評価に共用。
    """
    return sum(
        offgas.get(c, 0.0) * HHV_MJ_PER_KMOL.get(c, 0.0) / 1000.0
        for c in offgas
    )


# ---------------------------------------------------------------------------
# Hasebe 式 (9)(10) 関連ヘルパー
# ---------------------------------------------------------------------------

def _classify_opex_term(key: str) -> str:
    """opex dict のキーから Hasebe 式 (10) の区分を返す。

    Returns
    -------
    'UT'             : 用役費 (Hasebe 式の C_UT, 1.23× 適用対象)
                       電力・蒸気・冷却水・冷媒・燃料 (Reactor予熱炉燃料含む)・PSA予熱 等
    'RM'             : 原料費 (Hasebe 式の C_RM, 1.23× 適用対象)
                       Fresh LPG 等
    'CATALYST_OUT'   : 触媒・吸着剤交換費。Hasebe 式 (10) の枠外で別建て計上
                       (1.23× は掛けない)。本プロセス固有の消耗品扱い。
    'HASEBE_AGGR'    : Hasebe 集計項 (再計算時に剥がして再算出するためマーク)
    """
    if key.startswith(HASEBE_AGGR_PREFIX):
        return 'HASEBE_AGGR'
    # 膜の定期交換費 ('Mem膜交換') も触媒・吸着剤と同じ Hasebe 式 (10) 枠外の
    #   消耗品交換費として扱う (1.23× は掛けない)。
    if '触媒' in key or '吸着剤' in key or '活性炭交換' in key or '膜交換' in key:
        return 'CATALYST_OUT'
    if '原料費' in key or 'Fresh LPG' in key:
        return 'RM'
    return 'UT'   # default: 装置別の用役費 (電力・熱・燃料・HI/Stage2: ... 等)


def _count_main_equipment(capex: dict) -> int:
    """Hasebe 式 (9) 用の主要機器数 N_eq を数える。

    Hasebe §3.3 本文「コンプレッサー、塔、反応器、加熱器、熱交換器など」に従う。

    数え方の判断 (本プロセス固定構成 = 18 単位想定):
      - 蒸留塔 (Dist1/2/3) は付帯 HE (リボイラ・コンデンサ・フィード予熱)
        を 1 塔に含めて 1 単位扱い (Hasebe/Turton 慣行)。
      - Reactor 予熱炉は CAPEX 上は反応器バンドルだが、運転監視は独立した
        加熱器として N_eq に +1 計上。
      - PSA は容器 + ローテーション弁系を 1 単位 (吸着剤エントリ 'PSA活性炭'
        は機器ではなく内容物なので除外)。
      - 膜モジュールは補機含めて 1 単位 ('Mem気化器' 'Mem F圧縮機' 'Mem P圧縮機'
        'Mem冷却器' 'Mem膜本体' は別個に計上)。
      - スイング 2 基反応器は 1 系統運用として 1 単位 (Hasebe 例題に倣う)。
    """
    # capex dict から「物理機器 1 単位 = 1 カウント」を抽出
    # 'PSA活性炭' は機器ではなく内容物なので除外。'PSA容器' を PSA 1 単位として扱う。
    physical_unit_keys = {
        k for k in capex
        if k != 'PSA活性炭' and not k.startswith('[Hasebe')
    }
    n = len(physical_unit_keys)
    # CAPEX 上は反応器バンドルだが運転監視は独立した加熱炉として +1
    n += 1   # Reactor予熱炉
    return n


def _compute_labor_cost_okuyen(N_eq: int) -> float:
    """Hasebe 式 (9) で 1 班人数を算出し、4 直で総運転員 → 年間労務費 [億円/年]。

    式: 1 班の人数 = ceil( sqrt(6.29 + HASEBE_NOL_COEFF · N_eq) )
        総運転員  = 1 班の人数 × HASEBE_SHIFT_MULTIPLIER (=4)
        C_OL      = 総運転員 × OPERATOR_ANNUAL_SALARY_JPY / 1e8
    """
    n_per_shift = math.ceil(math.sqrt(6.29 + HASEBE_NOL_COEFF * N_eq))
    total_operators = n_per_shift * HASEBE_SHIFT_MULTIPLIER
    return total_operators * OPERATOR_ANNUAL_SALARY_JPY / 1.0e8


def apply_hasebe_aggregation(
    opex_raw: dict,
    total_capex: float,
    N_eq: int,
) -> dict:
    """装置別の生 opex dict に Hasebe 式 (10) の集計項を追加した新 dict を返す。

    既存の Hasebe 集計項 (`HASEBE_AGGR_PREFIX` で始まるキー) は事前に剥がして
    再計算する。これにより HI/HEN 後にこの関数を再呼び出しすると、変更された
    C_UT を反映した最新の Hasebe 集計が得られる。

    Hasebe 式 (10) のうち本実装が反映する項:
      - 0.180·C_TM            (保全 + 諸経費 + 税保険) ─ 集計項として追加
      - 2.73·C_OL             (労務費)                ─ 集計項として追加
      - 0.23·C_UT             (用役費の 23% 上乗せ)   ─ 集計項として追加
                              (装置別の生エントリが既に 1.00·C_UT を表すため、
                               残り 0.23·C_UT を delta として加算)
      - 0.23·C_RM             (原料費の 23% 上乗せ)   ─ 集計項として追加
      - 1.23·C_WT             (廃棄物処理費)          ─ PDH では 0 で集計項なし
      - 触媒・吸着剤交換費      (式外、別建て)          ─ 生エントリのまま残す
      - 減価償却費             (CAPEX/n)               ─ TAC 側で別途加算
    """
    # 既存集計項を除去
    raw = {
        k: v for k, v in opex_raw.items()
        if _classify_opex_term(k) != 'HASEBE_AGGR'
    }

    C_UT = sum(v for k, v in raw.items() if _classify_opex_term(k) == 'UT')
    C_RM = sum(v for k, v in raw.items() if _classify_opex_term(k) == 'RM')

    C_TM = total_capex
    C_OL = _compute_labor_cost_okuyen(N_eq)
    C_WT = HASEBE_C_WT_OKUYEN_PER_YEAR

    # 用役費・原料費は装置別エントリが 1.00× で計上済みなので、追加分のみを delta 計上
    multiplier_delta = HASEBE_COEFF_C_UT_WT_RM - 1.0   # = 0.23

    out = dict(raw)
    out[HASEBE_AGGR_PREFIX + f'{HASEBE_COEFF_C_TM:.3f}·C_TM (保全+諸経費)'] = (
        HASEBE_COEFF_C_TM * C_TM
    )
    out[HASEBE_AGGR_PREFIX + f'{HASEBE_COEFF_C_OL:.2f}·C_OL (労務費, N_eq={N_eq})'] = (
        HASEBE_COEFF_C_OL * C_OL
    )
    out[HASEBE_AGGR_PREFIX + f'{multiplier_delta:.2f}·C_UT (用役費 上乗せ分)'] = (
        multiplier_delta * C_UT
    )
    out[HASEBE_AGGR_PREFIX + f'{multiplier_delta:.2f}·C_RM (原料費 上乗せ分)'] = (
        multiplier_delta * C_RM
    )
    # PDH は気相反応で水処理対象廃棄物が実質発生しない (使用済み触媒は触媒交換費に
    # 内包) ため C_WT = 0 億円/年で固定。将来 C_WT > 0 になった場合でも if 判定で
    # 項を落とすと OPEX 集計から漏れるため、本実装では常に加算する
    # (C_WT = 0 のときは加算値も 0 なので副作用なし)。
    out[HASEBE_AGGR_PREFIX + f'{HASEBE_COEFF_C_UT_WT_RM:.2f}·C_WT (廃棄物処理)'] = (
        HASEBE_COEFF_C_UT_WT_RM * C_WT
    )
    return out


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
        'Desuper':     R['desuper'].equipment.CAPEX,    # Comp2b → Dist2 desuperheater
        'Dist2':       R['r2'].equipment.CAPEX,
        'PSA容器':     R['r_psa'].equipment.CAPEX_vessels,
        'PSA活性炭':   R['r_psa'].equipment.CAPEX_adsorbent,
        'MemPrecool':  R['mem_precool'].equipment.CAPEX,
        'Mem気化器':   R['r_mem'].equipment.CAPEX_vap,
        'Mem F圧縮機': R['r_mem'].equipment.CAPEX_comp_feed,
        'Mem P圧縮機': R['r_mem'].equipment.CAPEX_comp_prod,
        'Mem冷却器':   R['r_mem'].equipment.CAPEX_cond,
        # 膜圧縮機の多段化に伴う段間冷却器。単段 (n=1) なら 0。
        'Mem段間冷却器': R['r_mem'].equipment.CAPEX_intercool,
        'Mem膜本体':   R['r_mem'].equipment.CAPEX_mem,
        'Dist3':       R['r3'].equipment.CAPEX,
    }

    # 蒸留塔のリボイラ熱媒は equipment.reb_utility_name/単価から動的に選択
    # (utility_selector 経由)。
    r1_eq, r2_eq, r3_eq = R['r1'].equipment, R['r2'].equipment, R['r3'].equipment
    opex = {
        'Pump1電力':         _ele(R['pump1'].equipment.W_kW),
        'Comp2a電力':        _ele(R['comp2a'].equipment.W_kW),
        'Comp2b電力':        _ele(R['comp2b'].equipment.W_kW),
        'MemF圧縮機電力':    _ele(R['r_mem'].equipment.W_feed_kW),
        'MemP圧縮機電力':    _ele(R['r_mem'].equipment.W_prod_kW),
        f'Dist1リボイラ ({r1_eq.reb_utility_name})':
            _heat(r1_eq.Q_reb, r1_eq.reb_utility_jpy_per_GJ),
        f'Dist2リボイラ ({r2_eq.reb_utility_name})':
            _heat(r2_eq.Q_reb, r2_eq.reb_utility_jpy_per_GJ),
        f'Dist3リボイラ ({r3_eq.reb_utility_name})':
            _heat(r3_eq.Q_reb, r3_eq.reb_utility_jpy_per_GJ),
        'Mem気化器蒸気':     _heat(R['r_mem'].equipment.Q_vap_kW, LP_STEAM_JPY_PER_GJ),
        # 蒸留塔フィード予熱を独立計上。単価は塔のリボイラと同じ utility を流用
        # (フィード予熱は塔体内ではなく専用 HE 想定だが熱源は同レンジで近似)。
        f'Dist1フィード予熱 ({r1_eq.reb_utility_name})':
            _heat(r1_eq.Q_feed_preheat_kW, r1_eq.reb_utility_jpy_per_GJ),
        f'Dist2フィード予熱 ({r2_eq.reb_utility_name})':
            _heat(r2_eq.Q_feed_preheat_kW, r2_eq.reb_utility_jpy_per_GJ),
        f'Dist3フィード予熱 ({r3_eq.reb_utility_name})':
            _heat(r3_eq.Q_feed_preheat_kW, r3_eq.reb_utility_jpy_per_GJ),
    }

    # 反応器プリヒーター (GJ/h → kW 換算)
    # 燃料消費 = Q_preheat / η_furnace (LHV ベース、化工便覧 18·4·3 表 18·11)
    Q_preheat_kW = R['r_rx'].effluent.Q_preheat * 1.0e9 / 3600.0 / 1000.0
    Q_fuel_kW    = Q_preheat_kW / FURNACE_EFFICIENCY
    opex['Reactor予熱燃料'] = _heat(Q_fuel_kW, FUEL_JPY_PER_GJ)

    # cooler.py が utility_selector で選んだユーティリティ名・単価を equipment に
    # 格納するため、それを直接使う (ハードコードの COOLING_WATER_JPY_PER_GJ は使わない)。
    cooled_eq      = R['cooled'].equipment
    intercool_eq   = R['intercool'].equipment
    desuper_eq     = R['desuper'].equipment
    mem_precool_eq = R['mem_precool'].equipment
    opex[f"Cooler({cooled_eq.utility_name})"] = _heat(
        abs(cooled_eq.Q_duty_kW), cooled_eq.utility_jpy_per_GJ)
    opex[f"Intercool({intercool_eq.utility_name})"] = _heat(
        abs(intercool_eq.Q_duty_kW), intercool_eq.utility_jpy_per_GJ)
    opex[f"Desuper({desuper_eq.utility_name})"] = _heat(
        abs(desuper_eq.Q_duty_kW), desuper_eq.utility_jpy_per_GJ)
    opex[f"MemPrecool({mem_precool_eq.utility_name})"] = _heat(
        abs(mem_precool_eq.Q_duty_kW), mem_precool_eq.utility_jpy_per_GJ)

    # 蒸留塔のコンデンサ単価を equipment.cond_utility_* から動的選択
    # (utility_selector 経由)。Dist2 のように T_top が氷点下なら自動的に
    # エチレン冷媒 (~9000 円/GJ) が選ばれる。
    opex[f'Dist1コンデンサ ({r1_eq.cond_utility_name})'] = _heat(
        r1_eq.Q_cond, r1_eq.cond_utility_jpy_per_GJ,
    )
    opex[f'Dist2コンデンサ ({r2_eq.cond_utility_name})'] = _heat(
        r2_eq.Q_cond, r2_eq.cond_utility_jpy_per_GJ,
    )
    opex[f'Dist3コンデンサ ({r3_eq.cond_utility_name})'] = _heat(
        r3_eq.Q_cond, r3_eq.cond_utility_jpy_per_GJ,
    )
    # 膜の冷却器は内部実装が COOLING_WATER 前提。
    opex['Mem冷却器冷水']      = _heat(R['r_mem'].equipment.Q_cond_kW,
                                        COOLING_WATER_JPY_PER_GJ)

    # 膜圧縮機の多段化に伴う段間冷却 (冷却水) を計上。多段化で圧縮仕事 (= 電力 OPEX)
    # が下がる代わりに段間冷却器の CAPEX・冷却水 OPEX が生じるため両方を計上する。
    # HI 統合済: 本キーは heat_integration.classify_heat_opex_key で 'cold' 登録され、
    #   extract_streams の 'H_mem_intercool' ホットストリームと対応する。apply_hi=True 時は
    #   本フラット値が削除され HI tier 別 OPEX へ差し替わる (= 段間排熱がピンチ回収対象、
    #   二重計上なし)。apply_hi=False (生 OPEX) 時は本値が冷却水費として 1 回だけ計上される。
    Q_mem_intercool = R['r_mem'].equipment.Q_intercool_kW
    if Q_mem_intercool > 0:
        opex['Mem段間冷却 冷水'] = _heat(Q_mem_intercool, COOLING_WATER_JPY_PER_GJ)

    opex['Reactor触媒交換'] = (R['r_rx'].equipment.Catalyst_Weight_Total
                                * CATALYST_JPY_PER_KG
                                / CATALYST_LIFE_YEARS / 1.0e8)
    opex['PSA活性炭交換']   = R['r_psa'].equipment.OPEX_adsorbent_okuyen_per_year

    # 膜モジュール定期交換費 (レポート §7.6 の C_mem/τ_mem)。
    # 膜本体 CAPEX を耐用年数 (!仮置き MEM_LIFETIME_YEARS) で割った年均等費。
    # 触媒・PSA 活性炭交換と同じく Hasebe 式枠外の消耗品交換費 (1.23× なし)。
    # 膜耐用年数は感度解析用に env で上書き可 (PDH_MEM_LIFETIME_YEARS)。未設定なら !仮置き 既定。
    _mem_life = float(os.environ.get('PDH_MEM_LIFETIME_YEARS', str(MEM_LIFETIME_YEARS)))
    _capex_mem = R['r_mem'].equipment.CAPEX_mem
    if _capex_mem is not None and _capex_mem < PENALTY_CAPEX_THRESHOLD_OKUYEN and _mem_life > 0:
        opex['Mem膜交換'] = _capex_mem / _mem_life

    # PSA 予熱を OPEX に計上。T_in (Dist2 塔頂、partial cond で ~-20°C) →
    # T_abs (25°C) へ温度調整。加熱方向: LP Steam、冷却方向: 冷却水 (実機ではほぼ加熱方向)。
    Q_psa_preheat = R['r_psa'].equipment.Q_preheat_kW
    if Q_psa_preheat > 0:
        opex['PSA予熱 (LP Steam)']   = _heat(Q_psa_preheat,      LP_STEAM_JPY_PER_GJ)
    elif Q_psa_preheat < 0:
        opex['PSA予熱冷却 (冷却水)'] = _heat(abs(Q_psa_preheat), COOLING_WATER_JPY_PER_GJ)

    # ---- 原料費を OPEX に追加 (TAC に含める標準的な定義) ----
    # TAC には化工標準で原料費を含める (Sinnott §6.5, Turton §8.2)。
    # utility と並べて opex dict に置く。
    fresh_F = R['pump1'].outlet.F_in        # Pump1 inlet/outlet 同じ組成
    F_C3H8_feed  = fresh_F.get('A', 0.0)
    F_C4H10_feed = fresh_F.get('Z', 0.0)
    # C3H8 と C4H10 で個別単価採用 (n-Butane は約 10% 高い実勢)
    opex['Fresh LPG 原料費'] = (
        _annual_okuyen(F_C3H8_feed,  MW['A'], LPG_C3H8_JPY_PER_KG)
        + _annual_okuyen(F_C4H10_feed, MW['Z'], LPG_C4H10_JPY_PER_KG)
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

    # Dist1 塔底 (C4H10 富化) も燃料系へ送って CR 計上。
    # LPG (C3H8+C4H10) は HHV ~50 MJ/kg と高位、製品 spec から外れた C4 を
    # burning するのは工業的に標準 (refinery fuel gas header)。
    # FUEL_JPY_PER_GJ (1500 円/GJ ≒ LNG 相当) で評価する保守的計上。
    dist1_bot_GJ_per_h = _offgas_GJ_per_h(R['r1'].bottom.F_in)
    revenue['Dist1 塔底 燃料クレジット'] = (
        dist1_bot_GJ_per_h * OPERATING_HOURS_PER_YEAR * FUEL_JPY_PER_GJ / 1.0e8
    )

    return capex, opex, revenue


def calculate_economics(one_pass: dict, mw_C3H6_kg_per_kmol: float) -> Economics:
    """CAPEX/OPEX/Revenue 集計と TAC・Profit・製品単価を計算。

    opex dict は Hasebe 式 (10) の集計項 (HASEBE_AGGR_PREFIX で始まるキー) を
    含む形で返される。装置別の生エントリと Hasebe 集計エントリの合計が
    total_opex (= 製造コスト − 減価償却費) に一致する。
    """
    capex, opex_raw, revenue = collect_capex_opex(one_pass)

    # ペナルティ装置 (CAPEX >= PENALTY_CAPEX_THRESHOLD_OKUYEN = 1e8 億円相当) は集計から除外
    total_capex   = sum(v for v in capex.values() if v < PENALTY_CAPEX_THRESHOLD_OKUYEN)
    N_eq          = _count_main_equipment(capex)
    opex          = apply_hasebe_aggregation(opex_raw, total_capex, N_eq)
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
