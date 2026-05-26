"""
exp 結果表示モジュール。

各セクションを独立した関数として持ち、`display_full_results()` で
全部まとめて出すか、個別に呼ぶ。
"""

from stream.stream import ProcessStream
from src.cost_parameters import (
    ELECTRICITY_JPY_PER_KWH, LP_STEAM_JPY_PER_GJ,
    COOLING_WATER_JPY_PER_GJ, FUEL_JPY_PER_GJ,
    CATALYST_JPY_PER_KG, CATALYST_LIFE_YEARS,
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


def show_input_snapshot(design, config=None, eval_kwargs: dict = None) -> None:
    """実験条件のスナップショットを出力 (PDF 冒頭・log 先頭に入る)。

    Parameters
    ----------
    design : FlowsheetDesignVars
        全 unit の設計変数。dist1/2/3 の solver_method も含む。
    config : OperatingConfig, optional
        現状未使用 (将来 spec 値などを表示する場合に拡張)。
    eval_kwargs : dict, optional
        evaluate() に渡された kwargs (apply_hi, apply_stage2, hi_dT_min_K 等)。
        None なら省略表示。
    """
    from datetime import datetime as _dt
    print("=" * 72)
    print(f"  実行スナップショット ({_dt.now().strftime('%Y-%m-%d %H:%M:%S')})")
    print("=" * 72)
    print()
    sw = design.swing
    print("[反応器 (Swing)]")
    print(f"  T_in = {sw.T_in} K, z_cat = {sw.z_cat} m, t_cyc = {sw.t_cyc} min, D = {sw.D} m")
    print()
    psa = design.psa
    print("[PSA]")
    print(f"  D_col = {psa.D_col} m, L_bed = {psa.L_bed} m, desorption_target = {psa.desorption_target}")
    print()
    mem = design.mem
    print("[Membrane]")
    print(f"  P_H = {mem.P_H/1e5:.2f} bar, P_L = {mem.P_L/1e5:.2f} bar, A_mem = {mem.A_mem:.2e} m²")
    print()
    for label, col in [('Dist1 (脱ブタン塔)', design.dist1),
                       ('Dist2 (脱エタン塔, partial cond)', design.dist2),
                       ('Dist3 (C3 スプリッタ)', design.dist3)]:
        print(f"[{label}]")
        # N_feed は探索対象外で core 側 Kirkbride 推奨を自動採用。
        # 入力値表示は誤解を招くため省略 (実値は results.equipment.N_feed_kirkbride)。
        print(f"  P = {col.P_col/1e5:.1f} bar, N = {col.N_stages}, "
              f"R = {col.reflux_ratio}  (N_feed: Kirkbride 自動採用)")
        print(f"  solver = {col.solver_method}")
        print()

    if eval_kwargs:
        print("[評価オプション (evaluate 引数)]")
        for k, v in eval_kwargs.items():
            print(f"  {k:<22} = {v}")
        print()

    print("[使用モデル・出典]")
    print(f"  EOS              : Peng-Robinson 1976 (src/eos.py)")
    print(f"  bubble_point_T   : thermo (CalebBell, MIT v0.6.0) 内部実装、API は src/eos.py")
    print(f"  蒸留塔 (FUG)      : Fenske-Underwood-Gilliland shortcut")
    print(f"  蒸留塔 (rigorous) : Wang-Henke MESH + Newton + Wegstein "
          f"(Seader/Henley/Roper Ch.10.4)")
    print(f"  HE U 値           : 第17回プロセスデザイン学生コンテスト Ver.2.0 §4-4 表")
    print(f"  蒸留塔径 G*       : contest §4-2 (SF=0.8, K=0.05 m/s)")
    print(f"  Furnace 効率 0.85: 化工便覧 18·4·3 表 18·11")
    print(f"  Compressor η_poly 0.75 (polytropic): 化工便覧 p.333")
    print(f"  Pump 効率 0.70   : 化工便覧 5·6·4 例題 5·8")
    print()
    print("=" * 72)


def show_stream(label: str, stream) -> None:
    """ストリーム表示: 成分流量 [kmol/h] + mol% を併記、T/P/総流量を 2 行目に。"""
    total = stream.total_flow()
    parts = []
    for k, v in sorted(stream.F_in.items()):
        if v > 0.01:
            pct = v / total * 100.0 if total > 0 else 0.0
            parts.append(f"{_COMP_NAMES.get(k, k)}:{v:.1f}({pct:.1f}%)")
    print(f"  {label}: {', '.join(parts)}")
    print(f"  {' ' * len(label)}  T={stream.T_in - 273.15:.0f}°C  "
          f"P={stream.P_in / 1e5:.1f}bar  F={total:.1f}kmol/h")


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
    # V_cat per vessel (= 触媒だけの体積、200 m³ 制約対象) と
    # V_vessel per vessel (= 容器の物理体積、V_cat / (1-eps)) を両方表示
    eps_assumed = 0.5  # FixedParams.eps と整合
    v_cat_per_vessel = eq.V_vessel_actual * (1.0 - eps_assumed)
    print(f"  [Reactor]    {eq.N_parallel}並列 × {eq.N_swing_sets} swing = "
          f"{eq.N_reactors_total} 基  "
          f"V_cat/基={v_cat_per_vessel:.0f} m³ (≤200制約)  "
          f"V_vessel/基={eq.V_vessel_actual:.0f} m³  "
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
    """CAPEX 内訳 [億円] — share% 付き、寄与降順 で見やすく。"""
    hdr("CAPEX 内訳 [億円]  (share% 降順、合計に対する寄与率)")
    total = econ.total_capex
    items = [(n, v) for n, v in econ.capex.items() if v < 1e6]
    items.sort(key=lambda kv: -kv[1])
    penalty_items = [(n, v) for n, v in econ.capex.items() if v >= 1e6]
    print(f"  {'装置':<14} {'CAPEX [億円]':>13} {'share':>8}")
    print(f"  {'-'*14} {'-'*13} {'-'*8}")
    for n, v in items:
        pct = v / total * 100.0 if total > 0 else 0.0
        bar = '■' * max(1, int(pct / 5))   # 5% 刻みのバー、視覚化
        print(f"  {n:<14} {v:>13.4f} {pct:>6.1f}%  {bar}")
    for n, _ in penalty_items:
        print(f"  {n:<14}     ペナルティ")
    print(f"  {'-'*14} {'-'*13} {'-'*8}")
    print(f"  {'合計':<14} {econ.total_capex:>13.4f} {100.0:>6.1f}%")


def show_opex(econ) -> None:
    """OPEX 内訳を二段表示 (装置別 raw + Hasebe 式 (10) 集計項)。

    装置別: 電力・蒸気・冷却水・燃料・原料費・触媒/吸着剤交換 (全て 1.00 倍の生値)
    Hasebe 集計: 0.180·C_TM, 2.73·C_OL, 0.23·C_UT 上乗せ, 0.23·C_RM 上乗せ 等
                 (C_WT≈0, 減価償却費は CAPEX/n として TAC 側で別途加算)
    """
    from flowsheet.economics import HASEBE_AGGR_PREFIX

    raw_items    = {k: v for k, v in econ.opex.items() if not k.startswith(HASEBE_AGGR_PREFIX)}
    hasebe_items = {k: v for k, v in econ.opex.items() if k.startswith(HASEBE_AGGR_PREFIX)}

    hdr("OPEX 内訳 [億円/年]  Hasebe §3.4 式(10) ─ share% 降順、減価償却費は CAPEX/n で別計算")

    total = econ.total_opex
    raw_sorted = sorted(raw_items.items(), key=lambda kv: -kv[1])
    print(f"  ▼ 装置別 (用役・原料・触媒/吸着剤の 1.00 倍生値)")
    print(f"    {'項目':<28} {'OPEX [億円/年]':>15} {'share':>8}")
    print(f"    {'-'*28} {'-'*15} {'-'*8}")
    for n, v in raw_sorted:
        pct = v / total * 100.0 if total > 0 else 0.0
        bar = '■' * max(1, int(pct / 5))
        print(f"    {n:<28} {v:>15.4f} {pct:>6.1f}%  {bar}")
    raw_sum = sum(raw_items.values())
    print(f"    {'-'*28} {'-'*15} {'-'*8}")
    print(f"    {'(装置別 小計)':<28} {raw_sum:>15.4f} "
          f"{raw_sum/total*100 if total>0 else 0:>6.1f}%")

    if hasebe_items:
        print()
        print(f"  ▼ Hasebe 集計項 (式 (10) の 0.180·C_TM, 2.73·C_OL, 0.23 倍率上乗せ)")
        hasebe_sorted = sorted(hasebe_items.items(), key=lambda kv: -kv[1])
        for n, v in hasebe_sorted:
            label = n[len(HASEBE_AGGR_PREFIX):]
            pct = v / total * 100.0 if total > 0 else 0.0
            bar = '■' * max(1, int(pct / 5))
            print(f"    {label:<28} {v:>15.4f} {pct:>6.1f}%  {bar}")
        hasebe_sum = sum(hasebe_items.values())
        print(f"    {'-'*28} {'-'*15} {'-'*8}")
        print(f"    {'(Hasebe 小計)':<28} {hasebe_sum:>15.4f} "
              f"{hasebe_sum/total*100 if total>0 else 0:>6.1f}%")

    print(f"  {'='*55}")
    print(f"  {'OPEX 合計':<28} {econ.total_opex:>15.4f}  100.0%")


def show_revenue(econ) -> None:
    """Revenue 内訳 [億円/年] — share% 付き、寄与降順。"""
    hdr("Revenue 内訳 [億円/年]  (売上 + オフガス燃料クレジット、share% 降順)")
    total = econ.total_revenue
    items = sorted(econ.revenue.items(), key=lambda kv: -kv[1])
    print(f"  {'項目':<24} {'Revenue [億円/年]':>17} {'share':>8}")
    print(f"  {'-'*24} {'-'*17} {'-'*8}")
    for n, v in items:
        pct = v / total * 100.0 if total > 0 else 0.0
        bar = '■' * max(1, int(pct / 5))
        print(f"  {n:<24} {v:>17.4f} {pct:>6.1f}%  {bar}")
    print(f"  {'-'*24} {'-'*17} {'-'*8}")
    print(f"  {'Revenue 合計':<24} {econ.total_revenue:>17.4f}  100.0%")


def show_specs(specs, failure_reason: str) -> None:
    hdr("製品仕様 compliance")
    print(f"  C3H6 純度 (Dist3 塔頂) : {specs.c3h6_purity_wtfrac*100:6.3f} wt%"
          f"   {'✓' if specs.c3h6_pass else '✗'}"
          f" (spec ≥ {99.5:.1f} wt%)")
    print(f"  H2 純度 (PSA 製品)     : {specs.h2_purity_molfrac*100:6.3f} mol%"
          f"   {'✓' if specs.h2_pass else '✗'}"
          f" (spec ≥ {99.9:.1f} mol%)")
    if specs.threshold_high_kmol_h == float('inf'):
        _spec_range = f"≥ {specs.threshold_low_kmol_h:.2f}"
    else:
        _spec_range = f"{specs.threshold_low_kmol_h:.2f}〜{specs.threshold_high_kmol_h:.2f}"
    print(f"  生産量                 : {specs.production_kmol_h:7.2f} kmol/h"
          f" {'✓' if specs.production_pass else '✗'}"
          f" (spec {_spec_range} kmol/h)")
    if not specs.all_pass:
        print(f"  違反内訳 : {failure_reason}")


def show_hi_summary(result) -> None:
    """Heat Integration (pinch targeting) 結果を簡潔表示。"""
    if result.hi_result is None or result.economics_hi is None:
        return

    hr = result.hi_result

    hdr("Heat Integration (Pinch Targeting)")
    print(f"  ΔT_min={hr.T_pinch_hot_K - hr.T_pinch_cold_K:.0f}K, "
          f"ピンチ hot/cold={hr.T_pinch_hot_K-273.15:.1f}/{hr.T_pinch_cold_K-273.15:.1f}°C")
    print(f"  Q_H_min={hr.Q_H_min_kW/1000:.2f} MW, Q_C_min={hr.Q_C_min_kW/1000:.2f} MW, "
          f"N_HE_min={hr.N_HE_min}, A_total={hr.A_total_m2:.0f} m²")
    if not hr.feasible:
        print(f"  infeasible: {hr.message}")


def show_stage2_synthesis(result) -> None:
    """Stage 2 (HEN synthesis) の結果: 熱交換の組み合わせ (hot↔cold) と移動熱量を表示。

    設計判断 (2026-05-26 ユーザー要望): 旧版は 1 行サマリのみで matching 詳細を省略して
    いたが、「どの hot ストリームとどの cold ストリームを、何 kW 交換するか」(= HEN の
    実構成) を結果に出す。データは HENResult.matches (各 HEMatch に hot_label/cold_label/
    Q_kW/温度/A/CAPEX) に揃っているので表で展開する。exp1/exp3/special/main top-k の
    display_full_results は全て本関数を通るので一括で反映される。
    """
    if result.hen_result is None or result.economics_synth is None:
        return
    hr = result.hen_result

    hdr("Heat Integration Stage 2 (HEN 合成: 熱交換の組み合わせ・移動熱量)")
    print(f"  追加 process-process HE: {hr.n_process_HE} 機   "
          f"内部熱回収 = {hr.Q_recovered_kW/1000:6.2f} MW   "
          f"追加 CAPEX = {hr.CAPEX_added_okuyen:.3f} 億円")
    print(f"  残ユーティリティ: 加熱 hot = {hr.Q_hot_utility_kW/1000:6.2f} MW   "
          f"冷却 cold = {hr.Q_cold_utility_kW/1000:6.2f} MW")
    if not hr.feasible:
        print(f"  ⚠ HEN synthesis infeasible: {hr.message}")

    if hr.matches:
        print()
        print("  熱交換の組み合わせ (hot ストリーム → cold ストリーム / 移動熱量):")
        print(f"  {'#':>2} {'hot stream':<18} {'cold stream':<18} {'移動熱量':>9} "
              f"{'hot T(in>out)':>15} {'cold T(in>out)':>15} {'A[m2]':>7} {'CAPEX':>7} {'pinch':>6}")
        print(f"  {'-'*2} {'-'*18} {'-'*18} {'-'*9} {'-'*15} {'-'*15} {'-'*7} {'-'*7} {'-'*6}")
        for i, m in enumerate(hr.matches, 1):
            h_lab = m.hot_label if len(m.hot_label) <= 18 else m.hot_label[:17] + '…'
            c_lab = m.cold_label if len(m.cold_label) <= 18 else m.cold_label[:17] + '…'
            print(f"  {i:>2} {h_lab:<18} {c_lab:<18} {m.Q_kW/1000:>6.2f}MW "
                  f"{m.T_h_in_K-273.15:>6.0f}>{m.T_h_out_K-273.15:<6.0f} "
                  f"{m.T_c_in_K-273.15:>6.0f}>{m.T_c_out_K-273.15:<6.0f} "
                  f"{m.A_m2:>7.0f} {m.CAPEX_okuyen:>7.3f} {m.side:>6}")
    else:
        print("  (process-process マッチなし: 全ストリームを utility で処理)")


def show_tac_summary(result, C3H6_product: float) -> None:
    """raw / HI 後 / Stage 2 後の TAC・Profit を 1 つの表に統合して表示。

    最も実態に近い stage (= Stage 2 > HI > raw) を **最終値** としてマーカー付き表示。
    Revenue は通常 stage 間で不変なので 1 回だけ表示。
    """
    hdr("TAC・Revenue・Profit (全 stage 統合)")
    econ      = result.economics
    econ_hi   = result.economics_hi
    econ_syn  = result.economics_synth

    # 最終 stage を決定 (Stage 2 > HI > raw)
    if econ_syn is not None:
        final_econ, final_label = econ_syn,  "Stage 2 後"
    elif econ_hi is not None:
        final_econ, final_label = econ_hi,   "HI 後"
    else:
        final_econ, final_label = econ,      "raw"

    print(f"  CAPEX/{DEPRECIATION_YEARS}年(償却)  : "
          f"{econ.total_capex/DEPRECIATION_YEARS:9.3f} 億円/年")
    print(f"  Revenue           : {econ.total_revenue:9.3f} 億円/年  (売上 + 燃料CR)")
    print()
    # 3 stage を 1 表で比較
    print(f"  {'stage':<18} {'TAC':>10} {'Profit':>11} {'原単価':>12}")
    print(f"  {'-'*18} {'-'*10} {'-'*11} {'-'*12}")

    def _stage_row(label: str, e, is_final: bool):
        marker = "  ← 最終値" if is_final else ""
        print(f"  {label:<18} {e.TAC:>10.3f} {e.profit:>+11.3f} "
              f"{e.unit_jpy_per_t/1000:>9.1f} 円/kg{marker}")

    _stage_row("raw (HI 無し)",       econ,     final_label == "raw")
    if econ_hi is not None:
        _stage_row("HI 後 (Stage 1)", econ_hi,  final_label == "HI 後")
    if econ_syn is not None:
        _stage_row("Stage 2 後",      econ_syn, final_label == "Stage 2 後")

    print()
    # 最終 Profit と effective_TAC を強調表示
    sign_word = "黒字" if final_econ.profit >= 0 else "赤字"
    print(f"  ▶ 最終 Profit ({final_label})    : {final_econ.profit:+9.3f} 億円/年"
          f"  ({sign_word})")
    # 設計判断 (2026-05-22): 旧表示は `effective_TAC = TAC − Revenue + ペナルティ` の
    # 旧式を前提に `penalty_amount = effective_TAC - (TAC - Revenue)` で逆算していたが、
    # flowsheet/runner.py:413 で 2026-05-21 に純 TAC 最小化 (`eff_econ.TAC + soft_penalty`)
    # に変更済 → 旧式逆算では Revenue 項が紛れ込んだ過大値 (例: 27 億の実 penalty が
    # 839 億として表示) になっていた。現コード式に合わせて単純差分にする。
    penalty_amount = result.effective_TAC - final_econ.TAC
    if penalty_amount > 0.001:
        print(f"  + spec 違反ペナルティ        : {penalty_amount:+9.3f} 億円/年")
    print(f"  ▶ effective_TAC (最適化器)   : {result.effective_TAC:+9.3f} 億円/年"
          f"  (= TAC + ペナルティ)")
    print()
    if C3H6_product > 0:
        print(f"  C3H6 年間生産量              : {econ.annual_kg_C3H6/1000.0:.0f} ton/年")
        print(f"  C3H6 製造原単価 ({final_label}) "
              f": {final_econ.unit_jpy_per_t:.0f} 円/ton")
    print()
    print(f"  単価 (src/cost_parameters.py): "
          f"LPG {LPG_FEED_JPY_PER_KG}/C3H6 {C3H6_PRODUCT_JPY_PER_KG}/H2 {H2_PRODUCT_JPY_PER_KG} 円/kg, "
          f"電力 {ELECTRICITY_JPY_PER_KWH} 円/kWh, "
          f"LP蒸気/冷水/燃料 {LP_STEAM_JPY_PER_GJ}/{COOLING_WATER_JPY_PER_GJ}/{FUEL_JPY_PER_GJ} 円/GJ")
    print(f"  触媒 Cr2O3-Al2O3 (Catofin相当) {CATALYST_JPY_PER_KG:.0f}円/kg × {CATALYST_LIFE_YEARS:.1f}年, "
          f"稼働 {OPERATING_HOURS_PER_YEAR:.0f}h/年, 償却 {DEPRECIATION_YEARS}年")


def show_final_summary_box(result, C3H6_product: float) -> None:
    """ファイル末尾に置く目立つ最終サマリ。

    全角文字 (日本語・★・記号) の幅計算問題を避けるため、右端を閉じない
    片開き枠 (heavy double line) で構成。CAPEX/OPEX/TAC/Revenue/Profit を
    最終 stage (Stage 2 > HI > raw) で集約表示。
    """
    econ      = result.economics
    econ_hi   = result.economics_hi
    econ_syn  = result.economics_synth

    if econ_syn is not None:
        final_econ, final_label = econ_syn,  "Stage 2 後 (HEN synthesis)"
    elif econ_hi is not None:
        final_econ, final_label = econ_hi,   "HI 後 (Stage 1, pinch targeting)"
    else:
        final_econ, final_label = econ,      "raw (HI 無し)"

    profit_sign = "★ 黒字" if final_econ.profit >= 0 else "× 赤字"
    bar = "═" * 64
    sep = "─" * 64

    print()
    print(bar)
    print("                       ★ 最  終  サ  マ  リ ★")
    print(bar)
    print(f"  ベース stage   : {final_label}")
    print(sep)
    print(f"  CAPEX (合計)            : {final_econ.total_capex:>12.3f} 億円")
    print(f"  CAPEX 年償却 ({DEPRECIATION_YEARS}年)     : "
          f"{final_econ.total_capex/DEPRECIATION_YEARS:>12.3f} 億円/年")
    print(f"  OPEX  (合計)            : {final_econ.total_opex:>12.3f} 億円/年")
    print(f"  TAC   (= 償却 + OPEX)   : {final_econ.TAC:>12.3f} 億円/年")
    print(f"  Revenue                 : {final_econ.total_revenue:>12.3f} 億円/年")
    print(sep)
    print(f"  Profit = Rev − TAC      : {final_econ.profit:>+12.3f} 億円/年   "
          f"[ {profit_sign} ]")
    if C3H6_product > 0:
        unit_yen_kg = final_econ.unit_jpy_per_t / 1000.0
        prod_kt = final_econ.annual_kg_C3H6 / 1.0e6
        print(f"  C3H6 年産               : {prod_kt:>12.2f} kt/年")
        print(f"  C3H6 製造原単価         : {unit_yen_kg:>12.2f} 円/kg")
    print(f"  effective_TAC           : {result.effective_TAC:>+12.3f} 億円/年   "
          f"(最適化器の最小化対象)")
    print(bar)


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
    # HI と Stage 2 の物理メタデータ (pinch / HEN 構成) を先に表示
    show_hi_summary(result)
    show_stage2_synthesis(result)
    # raw/HI/Stage 2 を 1 表で比較 (中間サマリ)
    show_tac_summary(result, C3H6_product)
    # ★ 最終の目立つサマリボックス (一番下に必ず表示)
    show_final_summary_box(result, C3H6_product)
