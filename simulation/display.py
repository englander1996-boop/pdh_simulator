"""
exp 結果表示モジュール。

各セクションを独立した関数として持ち、`display_full_results()` で
全部まとめて出すか、個別に呼ぶ。
"""

from stream.stream import ProcessStream
from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
    CATALYST_PTSN_JPY_PER_KG, CATALYST_PTSN_LIFE_YEARS,
    OPERATING_HOURS_PER_YEAR, DEPRECIATION_YEARS,
    LPG_FEED_JPY_PER_KG, C3H6_PRODUCT_JPY_PER_KG, H2_PRODUCT_JPY_PER_KG,
)


_COMP_NAMES = {
    'A': 'C3H8', 'B': 'C3H6', 'C': 'H2',  'D': 'C2H4',
    'E': 'CH4',  'F': 'C2H6', 'Z': 'C4H10',
}


# ---------------------------------------------------------------------------
# 共通ヘルパ
# ---------------------------------------------------------------------------

def hdr(title: str) -> None:
    print(f"\n{'=' * 64}")
    print(f"  {title}")
    print('=' * 64)


def show_stream(label: str, stream) -> None:
    parts = [f"{_COMP_NAMES.get(k, k)}:{v:.1f}"
             for k, v in sorted(stream.F_in.items()) if v > 0.01]
    print(f"  {label}: {', '.join(parts)}")
    print(f"  {' ' * len(label)}  T={stream.T_in - 273.15:.0f}°C  "
          f"P={stream.P_in / 1e5:.1f}bar  F={stream.total_flow():.1f}kmol/h")


# ---------------------------------------------------------------------------
# 各セクション
# ---------------------------------------------------------------------------

def show_streams_overview(R: dict, F_C3H8_feed: float,
                          F_C4H10_feed: float, config) -> None:
    """ユニット間ストリーム遷移の表示 (各ユニットが次に何を渡すか)。"""
    hdr("ストリーム遷移 (各ユニット入出力)")

    # Fresh LPG (Pump1 入口) と Pump1 → Dist1
    fresh = ProcessStream(
        F_in={'A': F_C3H8_feed, 'Z': F_C4H10_feed},
        T_in=config.feed.T_K, P_in=config.feed.P_Pa,
    )
    show_stream("[in] Fresh LPG (→ Pump1)", fresh)
    show_stream("[Pump1 → Dist1]", R['pump1'].outlet)

    # Dist1 出口 (反応器系 / 廃棄)
    show_stream("[Dist1 → 膨張弁] 塔頂 (C3 富化)", R['r1'].top)
    show_stream("[Dist1 → 廃棄] 塔底 (C4 富化)",  R['r1'].bottom)
    show_stream("[膨張弁 → Reactor 系] 0.5bar 膨張後", R['dist1_top_rx'])

    # Recycle (Mem 保留 + Dist3 塔底)
    tear_d3  = R['tear_dist3_new']
    tear_mem = R['tear_mem_new']
    recycle_total = tear_d3['A'] + tear_d3['B'] + tear_mem['A'] + tear_mem['B']
    print(f"  [Recycle 合計] Mem 保留(A={tear_mem['A']:.1f}, B={tear_mem['B']:.1f})"
          f" + Dist3 塔底(A={tear_d3['A']:.1f}, B={tear_d3['B']:.1f})"
          f"  → {recycle_total:.1f} kmol/h")

    # Reactor 入口・出口
    show_stream("[Mixer → Reactor] Fresh + Recycle 合流", R['reactor_inlet'])
    show_stream("[Reactor → Cooler]", R['rx_out'])

    # Comp2 系列
    show_stream("[Cooler → Comp2a]",     R['cooled'].outlet)
    show_stream("[Comp2a → Intercool]",  R['comp2a'].outlet)
    show_stream("[Intercool → Comp2b]",  R['intercool'].outlet)
    show_stream("[Comp2b → Desuper]",    R['comp2b'].outlet)
    show_stream("[Desuper → Dist2]",     R['desuper'].outlet)

    # Dist2 → PSA / Mem
    show_stream("[Dist2 → PSA] 塔頂", R['r2'].top)
    show_stream("[Dist2 → MemPrecool] 塔底", R['r2'].bottom)
    show_stream("[MemPrecool → Mem] (気化フィード)", R['mem_precool'].outlet)

    # Mem 透過 → Dist3 (ProcessStream に再構築)
    mem_to_dist3 = ProcessStream(
        F_in={'A': R['r_mem'].product.F_C3H8, 'B': R['r_mem'].product.F_C3H6,
              'C': 0., 'D': 0., 'E': 0., 'F': 0.},
        T_in=R['r_mem'].product.T_out, P_in=R['r_mem'].product.P_out,
    )
    show_stream("[Mem 透過 → Dist3]", mem_to_dist3)

    # 製品/排出
    show_stream("[Dist3 → C3H6 製品] 塔頂", R['r3'].top)
    show_stream("[Dist3 → Recycle] 塔底",   R['r3'].bottom)

    # PSA / Mem 排出側 (ProcessStream でないので個別表示)
    rmem = R['r_mem']
    print(f"  [Mem 保留 → Recycle]    "
          f"C3H6:{rmem.retentate.F_C3H6:.1f}, C3H8:{rmem.retentate.F_C3H8:.1f}"
          f"  T={rmem.retentate.T_out-273.15:.0f}°C  P={rmem.retentate.P_out/1e5:.1f}bar")
    rpsa = R['r_psa']
    psa_p = rpsa.product
    psa_o = rpsa.offgas
    print(f"  [PSA → 製品 H2]         "
          f"C3H8:{psa_p.get('A',0):.1f}, C3H6:{psa_p.get('B',0):.1f}, "
          f"H2:{psa_p.get('C',0):.1f}, C2H4:{psa_p.get('D',0):.1f}, "
          f"CH4:{psa_p.get('E',0):.1f}, C2H6:{psa_p.get('F',0):.1f}  "
          f"→ {sum(psa_p.values()):.1f} kmol/h")
    print(f"  [PSA → 燃料系 オフガス] "
          f"C3H8:{psa_o.get('A',0):.1f}, C3H6:{psa_o.get('B',0):.1f}, "
          f"H2:{psa_o.get('C',0):.1f}, C2H4:{psa_o.get('D',0):.1f}, "
          f"CH4:{psa_o.get('E',0):.1f}, C2H6:{psa_o.get('F',0):.1f}  "
          f"→ {sum(psa_o.values()):.1f} kmol/h")


def show_unit_details(R: dict) -> None:
    """各ユニットの装置設計値を 1 行ずつ表示 (分析用)。"""
    hdr("ユニット詳細 (装置設計)")

    # ---- Pump1 ----
    eq = R['pump1'].equipment
    print(f"  [Pump1]      W={eq.W_kW:6.2f} kW  ρ_liq={eq.rho_liq:.0f} kg/m³  "
          f"V_dot={eq.V_dot*3600:.1f} m³/h")

    # ---- Dist1 ----
    eq = R['r1'].equipment
    print(f"  [Dist1]      D={eq.D_col:5.2f}m  H={eq.H_col:5.1f}m  "
          f"T_top={eq.T_top-273.15:+6.1f}°C  T_bot={eq.T_bot-273.15:+6.1f}°C  "
          f"N_min={eq.N_min:5.1f}  R_min={eq.R_min:5.2f}  feasible={eq.feasible}")
    print(f"               Cond: Q={eq.Q_cond:6.0f}kW  A={eq.A_cond_m2:6.0f}m²  "
          f"util={eq.cond_utility_name} ({eq.cond_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_cond:.3f}億円")
    print(f"               Reb : Q={eq.Q_reb:6.0f}kW  A={eq.A_reb_m2:6.0f}m²  "
          f"util={eq.reb_utility_name} ({eq.reb_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_reb:.3f}億円  (N_feed_kirkbride={eq.N_feed_kirkbride})")

    # ---- Reactor ----
    eq   = R['r_rx'].equipment
    perf = R['r_rx'].performance
    eff  = R['r_rx'].effluent
    print(f"  [Reactor]    {eq.N_parallel}並列 × {eq.N_swing_sets} swing = "
          f"{eq.N_reactors_total} 基  V/基={eq.V_vessel_actual:.0f} m³  "
          f"W_cat={eq.Catalyst_Weight_Total/1000:.1f} t")
    print(f"               転化率 X={perf.Conversion:5.1f}%  選択率 S={perf.Selectivity:5.1f}%  "
          f"Q_preheat={eff.Q_preheat:6.2f} GJ/h  T_out_avg={eff.T_out_avg-273.15:.0f}°C")

    # ---- Cooler ----
    eq = R['cooled'].equipment
    print(f"  [Cooler]     Q={eq.Q_duty_kW:+7.0f}kW  A={eq.A_est_m2:6.0f}m²  "
          f"utility={eq.utility_name} ({eq.utility_jpy_per_GJ:.0f}円/GJ)")

    # ---- Comp2a / Intercool / Comp2b ----
    eq = R['comp2a'].equipment
    print(f"  [Comp2a]     W={eq.W_kW:6.0f}kW  T_out={eq.T_out-273.15:.0f}°C")
    eq = R['intercool'].equipment
    print(f"  [Intercool]  Q={eq.Q_duty_kW:+7.0f}kW  A={eq.A_est_m2:6.0f}m²  "
          f"utility={eq.utility_name} ({eq.utility_jpy_per_GJ:.0f}円/GJ)")
    eq = R['comp2b'].equipment
    print(f"  [Comp2b]     W={eq.W_kW:6.0f}kW  T_out={eq.T_out-273.15:.0f}°C")
    eq = R['desuper'].equipment
    print(f"  [Desuper]    Q={eq.Q_duty_kW:+7.0f}kW  A={eq.A_est_m2:6.0f}m²  "
          f"utility={eq.utility_name} ({eq.utility_jpy_per_GJ:.0f}円/GJ)")

    # ---- Dist2 ----
    eq = R['r2'].equipment
    print(f"  [Dist2]      D={eq.D_col:5.2f}m  H={eq.H_col:5.1f}m  "
          f"T_top={eq.T_top-273.15:+6.1f}°C  T_bot={eq.T_bot-273.15:+6.1f}°C  "
          f"N_min={eq.N_min:5.1f}  R_min={eq.R_min:5.2f}  feasible={eq.feasible}")
    print(f"               Cond: Q={eq.Q_cond:6.0f}kW  A={eq.A_cond_m2:6.0f}m²  "
          f"util={eq.cond_utility_name} ({eq.cond_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_cond:.3f}億円")
    print(f"               Reb : Q={eq.Q_reb:6.0f}kW  A={eq.A_reb_m2:6.0f}m²  "
          f"util={eq.reb_utility_name} ({eq.reb_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_reb:.3f}億円  (N_feed_kirkbride={eq.N_feed_kirkbride})")

    # ---- PSA ----
    eq = R['r_psa'].equipment
    rpsa = R['r_psa']
    print(f"  [PSA]        {eq.N_total_columns} 基"
          f" ({eq.N_abs_parallel}並列×{eq.N_cycle_sets}セット)  "
          f"t_abs={eq.t_abs_sec:5.0f}s  t_des={eq.t_des_sec:5.0f}s  "
          f"u_0={eq.u_0:5.2f}m/s  W_ads={eq.W_adsorbent_kg/1000:.1f}t")
    print(f"               H2 回収率={rpsa.H2_recovery*100:5.1f}%  "
          f"CH4 捕捉率={rpsa.CH4_capture*100:5.1f}%  "
          f"Q_preheat={eq.Q_preheat_kW:+.0f}kW  "
          f"H2_loss(BD/Purge)={eq.H2_loss_blowdown_kmolh:.1f}/{eq.H2_loss_purge_kmolh:.1f}")

    # ---- MemPrecool ----
    eq = R['mem_precool'].equipment
    print(f"  [MemPrecool] Q={eq.Q_duty_kW:+7.0f}kW (顕熱 {eq.Q_sensible_kW:+.0f} + "
          f"潜熱 {eq.Q_latent_kW:+.0f})  A={eq.A_est_m2:.0f}m²  "
          f"utility={eq.utility_name} ({eq.utility_jpy_per_GJ:.0f}円/GJ)")

    # ---- Membrane ----
    eq = R['r_mem'].equipment
    rmem = R['r_mem']
    print(f"  [Membrane]   A_mem={eq.A_mem:.0f}m²  n_modules={eq.n_modules}  "
          f"stage_cut={rmem.stage_cut*100:5.1f}%  perm_purity(C3H6)={rmem.perm_purity*100:5.1f}%  "
          f"ret_purity(C3H6)={rmem.ret_purity*100:5.1f}%")
    print(f"               W_feed={eq.W_feed_kW:6.0f}kW  W_prod={eq.W_prod_kW:6.0f}kW  "
          f"Q_vap={eq.Q_vap_kW:6.0f}kW  Q_cond={eq.Q_cond_kW:6.0f}kW  "
          f"A_vap={eq.A_vap:.0f}m²  A_cond={eq.A_cond:.0f}m²")

    # ---- Dist3 ----
    eq = R['r3'].equipment
    print(f"  [Dist3]      D={eq.D_col:5.2f}m  H={eq.H_col:5.1f}m  "
          f"T_top={eq.T_top-273.15:+6.1f}°C  T_bot={eq.T_bot-273.15:+6.1f}°C  "
          f"N_min={eq.N_min:5.1f}  R_min={eq.R_min:5.2f}  feasible={eq.feasible}")
    print(f"               Cond: Q={eq.Q_cond:6.0f}kW  A={eq.A_cond_m2:6.0f}m²  "
          f"util={eq.cond_utility_name} ({eq.cond_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_cond:.3f}億円")
    print(f"               Reb : Q={eq.Q_reb:6.0f}kW  A={eq.A_reb_m2:6.0f}m²  "
          f"util={eq.reb_utility_name} ({eq.reb_utility_jpy_per_GJ:.0f}円/GJ)  "
          f"CAPEX={eq.CAPEX_reb:.3f}億円  (N_feed_kirkbride={eq.N_feed_kirkbride})")


def show_production(R: dict, F_C3H8_feed: float,
                    F_C4H10_feed: float, config) -> None:
    """生産・収率の表示。"""
    target_kmol_h = (config.product.target_mta * 1000.0
                     / config.product.mw_kg_per_kmol / OPERATING_HOURS_PER_YEAR)
    C3H6_product = R['r3'].top.F_in.get('B', 0.0)
    H2_product   = R['r_psa'].product.get('C', 0.0)
    yield_pct    = C3H6_product / F_C3H8_feed * 100.0 if F_C3H8_feed > 0 else 0.0
    tear_d3  = R['tear_dist3_new']
    tear_mem = R['tear_mem_new']
    recycle_C3H8 = tear_d3['A'] + tear_mem['A']
    recycle_C3H6 = tear_d3['B'] + tear_mem['B']

    hdr("生産・収率(収束時)")
    print(f"  C3H6 目標      : {target_kmol_h:7.2f} kmol/h"
          f"  ({config.product.target_mta:.0f} t/年 @ {OPERATING_HOURS_PER_YEAR:.0f} h/年)")
    print(f"  Fresh C3H8     : {F_C3H8_feed:7.2f} kmol/h   (外側ループで決定)")
    print(f"  Fresh C4H10    : {F_C4H10_feed:7.2f} kmol/h"
          f"   (LPG 組成 C3H8={config.feed.lpg_c3h8_mol_fraction} より)")
    print(f"  C3H6 製品      : {C3H6_product:7.2f} kmol/h"
          f"   (実収率 {yield_pct:5.2f}%, 仮定 {config.feed.yield_assumed*100:.1f}%)")
    print(f"  目標との差     : {C3H6_product - target_kmol_h:+7.3f} kmol/h")
    print(f"  H2 副産物      : {H2_product:7.2f} kmol/h")
    print(f"  Recycle C3H8   : {recycle_C3H8:7.2f} kmol/h"
          f"   (Fresh比 {recycle_C3H8/F_C3H8_feed*100:.0f}%)")
    print(f"  Recycle C3H6   : {recycle_C3H6:7.2f} kmol/h")
    print(f"  Reactor 入口   : {sum(R['reactor_inlet'].F_in.values()):7.2f} kmol/h"
          f"  (Fresh の {sum(R['reactor_inlet'].F_in.values())/F_C3H8_feed:.2f} 倍)")
    print(f"  Reactor 転化率 : {R['r_rx'].performance.Conversion:5.1f}%")
    print(f"  Reactor 選択率 : {R['r_rx'].performance.Selectivity:5.1f}%")


def show_capex(econ) -> None:
    hdr("CAPEX 内訳 [億円]")
    for n, v in econ.capex.items():
        if v < 1e6:
            print(f"  {n:<14}: {v:8.4f}")
        else:
            print(f"  {n:<14}:   ペナルティ")
    print(f"  {'-' * 26}")
    print(f"  {'合計':<14}: {econ.total_capex:8.4f}")


def show_opex(econ) -> None:
    hdr("OPEX 内訳 [億円/年]  (utility + 触媒 + 吸着剤 + 原料費 — 全て TAC に含む)")
    for n, v in econ.opex.items():
        print(f"  {n:<24}: {v:9.4f}")
    print(f"  {'-' * 36}")
    print(f"  {'OPEX 合計':<24}: {econ.total_opex:9.4f}")


def show_revenue(econ) -> None:
    hdr("Revenue 内訳 [億円/年]  (売上 + オフガス燃料クレジット)")
    for n, v in econ.revenue.items():
        print(f"  {n:<24}: {v:9.4f}")
    print(f"  {'-' * 36}")
    print(f"  {'Revenue 合計':<24}: {econ.total_revenue:9.4f}")


def show_specs(specs, failure_reason: str) -> None:
    hdr("製品仕様 compliance")
    print(f"  C3H6 純度 (Dist3 塔頂) : {specs.c3h6_purity_wtfrac*100:6.3f} wt%"
          f"   {'✓' if specs.c3h6_pass else '✗'}"
          f" (spec ≥ {99.5:.1f} wt%)")
    print(f"  H2 純度 (PSA 製品)     : {specs.h2_purity_molfrac*100:6.3f} mol%"
          f"   {'✓' if specs.h2_pass else '✗'}"
          f" (spec ≥ {99.9:.1f} mol%)")
    print(f"  生産量                 : {specs.production_kmol_h:7.2f} kmol/h"
          f" {'✓' if specs.production_pass else '✗'}"
          f" (片側 spec ≥ {specs.target_kmol_h * 0.99:.2f})")
    if not specs.all_pass:
        print(f"  違反内訳 : {failure_reason}")


def show_hi_summary(result) -> None:
    """Heat Integration (pinch targeting) 結果と HI 前後の TAC/Profit 比較を表示。

    HI 未適用 (result.hi_result is None) なら何もしない。
    """
    if result.hi_result is None or result.economics_hi is None:
        return

    hr      = result.hi_result
    econ    = result.economics
    econ_hi = result.economics_hi

    hdr("Heat Integration (Pinch Targeting)")
    print(f"  ΔT_min        : {(hr.T_pinch_hot_K - hr.T_pinch_cold_K):.0f} K")
    print(f"  Q_H_min       : {hr.Q_H_min_kW:8.0f} kW = "
          f"{hr.Q_H_min_kW/1000:6.2f} MW (HI 後の必要加熱量)")
    print(f"  Q_C_min       : {hr.Q_C_min_kW:8.0f} kW = "
          f"{hr.Q_C_min_kW/1000:6.2f} MW (HI 後の必要冷却量)")
    print(f"  ピンチ温度    : hot {hr.T_pinch_hot_K-273.15:6.1f}°C "
          f"/ cold {hr.T_pinch_cold_K-273.15:6.1f}°C")
    print(f"  最少 HE 数    : {hr.N_HE_min} (Linnhoff 簡易式)")
    print(f"  総伝熱面積    : {hr.A_total_m2:7.0f} m² (Bath 式 targeting)")
    print(f"  feasible      : {hr.feasible}"
          + (f"   ({hr.message})" if hr.message else ""))

    # HI 適用後の utility tier 別 OPEX 内訳
    hi_breakdown = {k: v for k, v in econ_hi.opex.items() if k.startswith('HI:')}
    if hi_breakdown:
        print(f"\n  HI 後 utility 内訳 [億円/年]:")
        for k, v in sorted(hi_breakdown.items(), key=lambda kv: -abs(kv[1])):
            print(f"    {k:<32}: {v:9.4f}")
        print(f"    {'-'*36}")
        print(f"    {'(HI utility 合計)':<32}: {sum(hi_breakdown.values()):9.4f}")


def show_hi_comparison(result) -> None:
    """HI 前後の TAC/Profit 比較表を表示。"""
    if result.economics_hi is None:
        return
    econ    = result.economics
    econ_hi = result.economics_hi
    hdr("HI 前 vs HI 後 比較 [億円/年]")

    # 熱関連 OPEX を集計 (元 economics 側で classify)
    from flowsheet.heat_integration import classify_heat_opex_key
    heat_opex_before = sum(
        v for k, v in econ.opex.items()
        if classify_heat_opex_key(k) is not None
    )
    nonheat_opex = sum(
        v for k, v in econ.opex.items()
        if classify_heat_opex_key(k) is None
    )
    heat_opex_after = sum(
        v for k, v in econ_hi.opex.items() if k.startswith('HI:')
    )

    print(f"  {'項目':<28} | {'HI なし':>10} | {'HI 後':>10} | {'差':>10}")
    print(f"  {'-'*28} | {'-'*10} | {'-'*10} | {'-'*10}")
    rows = [
        ("CAPEX/年 (償却)",      econ.total_capex/8,  econ_hi.total_capex/8),
        ("OPEX 熱系 (utility)",  heat_opex_before,    heat_opex_after),
        ("OPEX 非熱系 (触媒+原料等)", nonheat_opex,    nonheat_opex),
        ("OPEX 合計",            econ.total_opex,     econ_hi.total_opex),
        ("TAC",                  econ.TAC,            econ_hi.TAC),
        ("Revenue",              econ.total_revenue,  econ_hi.total_revenue),
        ("Profit",               econ.profit,         econ_hi.profit),
    ]
    for name, before, after in rows:
        diff = after - before
        sign = '+' if diff >= 0 else ''
        print(f"  {name:<28} | {before:>10.3f} | {after:>10.3f} "
              f"| {sign}{diff:>9.3f}")
    print()
    print(f"  C3H6 製造原単価:")
    print(f"    HI なし: {econ.unit_jpy_per_t:8.0f} 円/ton  ({econ.unit_jpy_per_t/1000:5.1f} 円/kg)")
    print(f"    HI 後  : {econ_hi.unit_jpy_per_t:8.0f} 円/ton  ({econ_hi.unit_jpy_per_t/1000:5.1f} 円/kg)")


def show_stage2_synthesis(result) -> None:
    """Stage 2 (HEN synthesis) の結果を表示。

    apply_stage2=True で評価したときのみ動作。HE matching 内訳と追加 CAPEX、
    Stage 1 vs Stage 2 の TAC/Profit 比較。
    """
    if result.hen_result is None or result.economics_synth is None:
        return

    hr = result.hen_result

    hdr("Stage 2: HEN Synthesis (Pinch Design Method, greedy + tick-off)")
    print(f"  feasible          : {hr.feasible}"
          + (f"   ({hr.message})" if hr.message else ""))
    print(f"  process-process HE: {hr.n_process_HE} 機 (= 追加 recovery exchanger)")
    print(f"  Q_recovered       : {hr.Q_recovered_kW:8.0f} kW = "
          f"{hr.Q_recovered_kW/1000:6.2f} MW (内部熱回収量)")
    print(f"  Q_hot_utility 残  : {hr.Q_hot_utility_kW:8.0f} kW "
          f"(synthesis 後の hot utility 必要量)")
    print(f"  Q_cold_utility 残 : {hr.Q_cold_utility_kW:8.0f} kW")
    print(f"  追加 HE CAPEX     : {hr.CAPEX_added_okuyen:.4f} 億円")

    if hr.matches:
        print(f"\n  HE matching 内訳:")
        print(f"    {'name':<24} {'side':<7} {'hot ↔ cold':<60} "
              f"{'Q[kW]':>8} {'A[m²]':>8} {'CAPEX':>8}")
        print(f"    {'-'*24} {'-'*7} {'-'*60} {'-'*8} {'-'*8} {'-'*8}")
        for m in hr.matches:
            pair = f"{m.hot_label[:28]:<28} ↔ {m.cold_label[:28]:<28}"
            print(f"    {m.name:<24} {m.side:<7} {pair[:60]} "
                  f"{m.Q_kW:>8.0f} {m.A_m2:>8.0f} {m.CAPEX_okuyen:>8.3f}")


def show_stage2_comparison(result) -> None:
    """Stage 1 (HI targeting) vs Stage 2 (synthesis) の比較。"""
    if result.economics_synth is None or result.economics_hi is None:
        return

    econ_hi    = result.economics_hi
    econ_syn   = result.economics_synth
    hr         = result.hen_result

    hdr("Stage 1 (targeting) vs Stage 2 (synthesis) 比較 [億円/年]")
    print(f"  {'項目':<28} | {'Stage 1':>10} | {'Stage 2':>10} | {'差':>10}")
    print(f"  {'-'*28} | {'-'*10} | {'-'*10} | {'-'*10}")

    rows = [
        ("CAPEX/年 (償却)",      econ_hi.total_capex/8,  econ_syn.total_capex/8),
        ("OPEX 合計",            econ_hi.total_opex,      econ_syn.total_opex),
        ("TAC",                  econ_hi.TAC,             econ_syn.TAC),
        ("Revenue",              econ_hi.total_revenue,   econ_syn.total_revenue),
        ("Profit",               econ_hi.profit,          econ_syn.profit),
    ]
    for name, before, after in rows:
        diff = after - before
        sign = '+' if diff >= 0 else ''
        print(f"  {name:<28} | {before:>10.3f} | {after:>10.3f} "
              f"| {sign}{diff:>9.3f}")

    print()
    print(f"  C3H6 製造原単価:")
    print(f"    Stage 1: {econ_hi.unit_jpy_per_t:8.0f} 円/ton  ({econ_hi.unit_jpy_per_t/1000:5.1f} 円/kg)")
    print(f"    Stage 2: {econ_syn.unit_jpy_per_t:8.0f} 円/ton  ({econ_syn.unit_jpy_per_t/1000:5.1f} 円/kg)")

    print(f"\n  ★ 解釈:")
    print(f"    Stage 1 (BO): 理論限界 OPEX、既存 HE CAPEX のみ → 楽観的 lower bound")
    print(f"    Stage 2 (top-k): 実 HEN 構成、追加 HE CAPEX 加算、実 OPEX → realistic")
    if econ_syn.profit < econ_hi.profit:
        diff = econ_hi.profit - econ_syn.profit
        print(f"    → Stage 2 で Profit が {diff:.2f} 億円/年 悪化 (top-k 評価で現実化)")


def show_tac_summary(result, C3H6_product: float) -> None:
    hdr("TAC・Revenue・Profit")
    econ = result.economics
    TAC     = econ.TAC
    revenue = econ.total_revenue
    profit  = econ.profit
    print(f"  CAPEX/{DEPRECIATION_YEARS}年(償却)   : "
          f"{econ.total_capex/DEPRECIATION_YEARS:9.4f}  億円/年")
    print(f"  OPEX 合計           : {econ.total_opex:9.4f}  億円/年"
          f"  (utility + 触媒 + 原料費)")
    print(f"  ────────────────────────────────────")
    print(f"  TAC                 : {TAC:9.4f}  億円/年  (年間総費用)")
    print(f"  Revenue             : {revenue:9.4f}  億円/年  (売上 + 燃料CR)")
    print(f"  ────────────────────────────────────")
    print(f"  Profit = Rev − TAC  : {profit:+9.4f}  億円/年  (正なら黒字)")
    penalty_amount = result.effective_TAC - (TAC - revenue)
    if penalty_amount > 0:
        print(f"  + spec 違反ペナルティ : {penalty_amount:+9.4f}  億円/年")
    print(f"  ────────────────────────────────────")
    print(f"  effective_TAC       : {result.effective_TAC:+9.4f}  億円/年"
          f"  (= TAC − Revenue + ペナルティ、最適化器の最小化対象)")
    print()
    if C3H6_product > 0:
        print(f"  C3H6 年間生産量       : {econ.annual_kg_C3H6/1000.0:.0f}  ton/年")
        print(f"  C3H6 製造原単価 (TAC) : {econ.unit_jpy_per_t:9.0f}  円/ton"
              f"  ({econ.unit_jpy_per_t/1000:5.1f} 円/kg)")
    print()
    print(f"  仮単価 (★全て src/cost_parameters.py に集約):")
    print(f"    LPG {LPG_FEED_JPY_PER_KG}円/kg  /  C3H6 売価 {C3H6_PRODUCT_JPY_PER_KG}円/kg"
          f"  /  H2 売価 {H2_PRODUCT_JPY_PER_KG}円/kg")
    print(f"    電力 {ELECTRICITY_JPY_PER_KWH}円/kWh  /  LP蒸気 {LP_STEAM_JPY_PER_GJ}円/GJ"
          f"  /  冷水 {COOLING_WATER_JPY_PER_GJ}円/GJ  /  燃料 {FUEL_JPY_PER_GJ}円/GJ")
    print(f"    PtSn触媒 {CATALYST_PTSN_JPY_PER_KG:.0f}円/kg / 寿命{CATALYST_PTSN_LIFE_YEARS:.0f}年")
    print(f"    稼働 {OPERATING_HOURS_PER_YEAR:.0f}h/年  /  償却 {DEPRECIATION_YEARS}年")


# ---------------------------------------------------------------------------
# ワンショット表示
# ---------------------------------------------------------------------------

def display_full_results(result, design, config) -> None:
    """フローシート評価結果の全セクションを順次表示する。

    Parameters
    ----------
    result  : flowsheet.runner.FlowsheetResult
    design  : flowsheet.design.FlowsheetDesignVars  (現状未使用、将来の比較表示用)
    config  : config.load.OperatingConfig
    """
    # solver-level 失敗 (economics=None) は早期表示で打ち切り
    if result.economics is None:
        print(f"\n設計NG: {result.failure_reason}")
        print(f"  effective_TAC = {result.effective_TAC:.0f} 億円/年 (固定打ち切り値)")
        print("  設計変数を見直してから再実行してください。")
        return

    R = result.solver.one_pass
    F_C3H8_feed  = result.solver.fresh_C3H8
    F_C4H10_feed = result.solver.fresh_C4H10
    C3H6_product = R['r3'].top.F_in.get('B', 0.0)

    show_streams_overview(R, F_C3H8_feed, F_C4H10_feed, config)
    show_unit_details(R)
    show_production(R, F_C3H8_feed, F_C4H10_feed, config)
    show_capex(result.economics)
    show_opex(result.economics)
    show_revenue(result.economics)
    show_specs(result.specs, result.failure_reason)
    show_tac_summary(result, C3H6_product)
    # HI を適用した場合のみ表示 (apply_hi=False では何も出ない)
    show_hi_summary(result)
    show_hi_comparison(result)
    # Stage 2 (HEN synthesis、apply_stage2=True) の結果も同様に optional 表示
    show_stage2_synthesis(result)
    show_stage2_comparison(result)
