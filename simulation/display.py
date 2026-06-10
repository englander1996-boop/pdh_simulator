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
    # 反応器モデルの種別で表示を分岐 (Catofin=既定 / 径方向流 / 軸流 の 3 種)。
    # Catofin は N_online、径方向流は bed_thickness を持つ (どちらも無ければ軸流)。
    if hasattr(sw, 'N_online'):
        print("[反応器 (Catofin 浅床軸流スイング)]")
        print(f"  T_in = {sw.T_in} K, t_cyc = {sw.t_cyc} min, "
              f"D = {sw.D} m, L_bed = {sw.L_bed} m, N_online = {sw.N_online} 基 (並列), "
              f"d_p = {sw.d_p * 1000:.3f} mm")
    elif hasattr(sw, 'bed_thickness'):
        print("[反応器 (Radial flow 径方向流)]")
        print(f"  T_in = {sw.T_in} K, t_cyc = {sw.t_cyc} min, "
              f"D_inner = {sw.D_inner} m, bed_thickness = {sw.bed_thickness} m, H = {sw.H} m "
              f"(r_i={sw.r_i:.2f}m, r_o={sw.r_o:.2f}m)")
    else:
        print("[反応器 (Axial swing 軸流)]")
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
        # フィード段の出どころは backend で異なる:
        #   sm / hysys     : 最適化されたフィード段 (col.hysys_feed_stage) を実際に使用。
        #   fug / rigorous : core 側で Kirkbride 推奨を自動採用 (入力 N_feed は無視)。
        if col.solver_method in ('sm', 'hysys') and col.hysys_feed_stage is not None:
            _feed_note = f"feed段 = {col.hysys_feed_stage} (最適化値)"
        else:
            _feed_note = "feed段 = Kirkbride 自動採用 (理論式)"
        print(f"  P = {col.P_col/1e5:.1f} bar, N = {col.N_stages}, "
              f"R = {col.reflux_ratio}  ({_feed_note})")
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


def show_unit_details(R: dict, design=None) -> None:
    """各ユニットの装置設計値を 1 行ずつ表示 (分析用)。

    design を渡すと各塔の solver_method に応じてフィード段ラベルを正しく出す。
    DistEquipment.N_feed_kirkbride フィールドの中身は backend で意味が異なる:
      - sm / hysys     : 入力した最適化フィード段 (provider/SM がそのまま格納)。
      - fug / rigorous : core が _kirkbride_feed_stage で計算した推奨段。
    フィールド名が "kirkbride" 固定で誤解を招くため、表示ラベルを backend で分岐する。
    """
    hdr("ユニット詳細 (装置設計)")

    def _feed_label(solver, n_feed):
        if solver in ('sm', 'hysys'):
            return f"feed段={n_feed} (最適化値)"
        if solver in ('fug', 'rigorous'):
            return f"feed段={n_feed} (Kirkbride推奨)"
        return f"feed段={n_feed}"
    _sm1 = design.dist1.solver_method if design is not None else None
    _sm2 = design.dist2.solver_method if design is not None else None
    _sm3 = design.dist3.solver_method if design is not None else None

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
          f"CAPEX={eq.CAPEX_reb:.3f}億円  ({_feed_label(_sm1, eq.N_feed_kirkbride)})")

    # ---- Reactor ----
    eq   = R['r_rx'].equipment
    perf = R['r_rx'].performance
    eff  = R['r_rx'].effluent
    # V_cat per vessel (= 触媒だけの体積、200 m³ 制約対象) と
    # V_vessel per vessel (= 容器の物理体積) を両方表示。
    # 軸流は V_cat=V_vessel·(1-eps) だが、径方向流は V_vessel=πr_o²H が中心捕集管
    #   void を含むため V_cat/V_vessel≠(1-eps)。どちらの幾何でも正しい W_cat から
    #   逆算する (W_cat=V_cat_total×N_swing×ρ_b)。
    rho_b = 900.0  # FixedParams.rho_b と整合
    # V_cat/基 (= 1 缶あたり触媒体積、200 m³ 制約対象): Catalyst_Weight_Total は
    #   全段(全 n_beds)合計なので、総基数 N_reactors_total で割る。
    #   N_reactors_total = N_parallel × N_swing_sets × n_beds (軸流/単段は n_beds=1)。
    #   Oleflex 多段は各段同一ジオメトリ ([[radial_flow]] 434) なので per-vessel は厳密に等価。
    if eq.N_reactors_total > 0:
        v_cat_per_vessel = eq.Catalyst_Weight_Total / (rho_b * eq.N_reactors_total)
    else:
        v_cat_per_vessel = 0.0
    _denom = eq.N_parallel * eq.N_swing_sets
    n_beds = max(eq.N_reactors_total // _denom, 1) if _denom > 0 else 1
    config_str = (f"{eq.N_parallel}並列 × {n_beds}段 × {eq.N_swing_sets}swing"
                  if n_beds > 1 else f"{eq.N_parallel}並列 × {eq.N_swing_sets} swing")
    print(f"  [Reactor]    {config_str} = "
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
          f"CAPEX={eq.CAPEX_reb:.3f}億円  ({_feed_label(_sm2, eq.N_feed_kirkbride)})")

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
    print(f"  [MemPrecool] Q={eq.Q_duty_kW:+7.0f}kW (顕熱 {eq.Q_sensible_kW:+.0f}kW + "
          f"潜熱 {eq.Q_latent_kW:+.0f}kW)  A={eq.A_est_m2:.0f}m²  "
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
          f"CAPEX={eq.CAPEX_reb:.3f}億円  ({_feed_label(_sm3, eq.N_feed_kirkbride)})")


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


def show_specs(specs, failure_reason: str, config=None) -> None:
    hdr("製品仕様 compliance")
    # spec 閾値は判定に使った config から取る (main/exp3 は SM Dist3 の mol↔wt 差吸収で
    # 99.45 wt% に緩和)。ハードコードだと c3h6_pass=✓ なのに「spec ≥ 99.5」と矛盾表示になる。
    # config 未指定時は contest 名目値にフォールバック。
    c3h6_thr = config.spec.c3h6_min_wtfrac * 100.0 if config is not None else 99.5
    h2_thr   = config.spec.h2_min_molfrac  * 100.0 if config is not None else 99.9
    print(f"  C3H6 純度 (Dist3 塔頂) : {specs.c3h6_purity_wtfrac*100:6.3f} wt%"
          f"   {'✓' if specs.c3h6_pass else '✗'}"
          f" (spec ≥ {c3h6_thr:.2f} wt%)")
    print(f"  H2 純度 (PSA 製品)     : {specs.h2_purity_molfrac*100:6.3f} mol%"
          f"   {'✓' if specs.h2_pass else '✗'}"
          f" (spec ≥ {h2_thr:.1f} mol%)")
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

    「どの hot ストリームとどの cold ストリームを、何 kW 交換するか」(= HEN の
    実構成) を表で展開する。データは HENResult.matches (各 HEMatch に hot_label/
    cold_label/Q_kW/温度/A/CAPEX) に揃っている。
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


def show_heat_balance(result, design) -> None:
    """プロセス各ストリームの詳細熱収支 + HI ユーティリティ tier 別配分を表示。

    HI (pinch targeting) の入力となる全 hot/cold ストリームを、顕熱・潜熱・相・
    温度域つきで一覧化し、内部熱回収量とユーティリティ tier 別の熱量・費用を集計する。
    show_hi_summary (ピンチ点・GCC 指標のみ) を補完する「細かい熱収支」の詳細版。

    extract_streams は run_one_pass の戻り値 (各装置の equipment) から直接 Q を読み取る
    ため、apply_hi=False でも per-stream の熱収支は表示できる (tier 配分は hi_result が
    あるときのみ)。
    """
    if result.solver is None or result.solver.one_pass is None:
        return
    from flowsheet.heat_integration import (
        extract_streams, get_default_utility_tiers, calc_hi_opex_okuyen,
    )

    streams = extract_streams(result.solver.one_pass, design.swing.T_in)
    if not streams:
        return

    hdr("熱収支 詳細 (プロセスストリーム別 + ユーティリティ tier 配分) [MW]")

    hot  = [s for s in streams if s.is_hot]
    cold = [s for s in streams if not s.is_hot]

    def _stream_row(s) -> None:
        # 潜熱区間 (T_in≈T_out の相変化) は温度域を 1 点表記にする
        if abs(s.T_in_K - s.T_out_K) < 0.05:
            t_range = f"{s.T_in_K-273.15:>6.1f}(相変化)"
        else:
            t_range = f"{s.T_in_K-273.15:>6.1f}→{s.T_out_K-273.15:<6.1f}"
        fcp = f"{s.F_Cp_kW_per_K:>7.2f}" if s.F_Cp_kW_per_K > 1e-6 else f"{'-':>7}"
        print(f"    {s.name:<24} {t_range:<14} "
              f"{s.Q_sensible_kW/1000:>8.3f} {s.Q_latent_kW/1000:>8.3f} "
              f"{s.Q_total_kW/1000:>8.3f}  {fcp}  {s.phase}")

    hdr_cols = (f"    {'stream':<24} {'T_in→T_out[°C]':<14} "
                f"{'Q_顕熱[MW]':>8} {'Q_潜熱[MW]':>8} {'Q_合計[MW]':>8}  {'F·Cp':>7}  相")
    sep_cols = f"    {'-'*24} {'-'*14} {'-'*8} {'-'*8} {'-'*8}  {'-'*7}"

    # ---- 与熱流体 (Hot) ----
    print(f"  ▼ 与熱流体 Hot (冷却が必要、{len(hot)} 本)  ※F·Cp は kW/K")
    print(hdr_cols)
    print(sep_cols)
    hot_sens = hot_lat = 0.0
    for s in sorted(hot, key=lambda s: -s.Q_total_kW):
        _stream_row(s)
        hot_sens += s.Q_sensible_kW
        hot_lat  += s.Q_latent_kW
    hot_total = hot_sens + hot_lat
    print(sep_cols)
    print(f"    {'(Hot 小計 = 総冷却要求)':<24} {'':<14} "
          f"{hot_sens/1000:>8.3f} {hot_lat/1000:>8.3f} {hot_total/1000:>8.3f}")

    # ---- 受熱流体 (Cold) ----
    print()
    print(f"  ▼ 受熱流体 Cold (加熱が必要、{len(cold)} 本)")
    print(hdr_cols)
    print(sep_cols)
    cold_sens = cold_lat = 0.0
    for s in sorted(cold, key=lambda s: -s.Q_total_kW):
        _stream_row(s)
        cold_sens += s.Q_sensible_kW
        cold_lat  += s.Q_latent_kW
    cold_total = cold_sens + cold_lat
    print(sep_cols)
    print(f"    {'(Cold 小計 = 総加熱要求)':<24} {'':<14} "
          f"{cold_sens/1000:>8.3f} {cold_lat/1000:>8.3f} {cold_total/1000:>8.3f}")

    # ---- HI 内部熱回収 + ユーティリティ tier 別配分 ----
    hr = result.hi_result
    if hr is None:
        print()
        print(f"  (HI 未適用: 総加熱要求 {cold_total/1000:.3f} MW / 総冷却要求 "
              f"{hot_total/1000:.3f} MW、内部熱回収なしの生 OPEX で計上)")
        return

    # 内部熱回収量 (pinch targeting): 回収 = 総要求 − 残ユーティリティ
    Q_rec_hot  = cold_total - hr.Q_H_min_kW   # cold 側で内部回収された分
    Q_rec_cold = hot_total  - hr.Q_C_min_kW   # hot  側で内部回収された分 (理論上 Q_rec_hot と一致)
    heating_tiers, cooling_tiers = get_default_utility_tiers()
    opex_calc = calc_hi_opex_okuyen(
        hr, heating_tiers, cooling_tiers, OPERATING_HOURS_PER_YEAR,
    )

    print()
    print(f"  ▼ ヒートインテグレーション (pinch targeting, ΔT_min="
          f"{hr.T_pinch_hot_K - hr.T_pinch_cold_K:.0f}K)")
    print(f"    総加熱要求       : {cold_total/1000:>8.3f} MW")
    print(f"    総冷却要求       : {hot_total/1000:>8.3f} MW")
    print(f"    内部熱回収 (pinch): {Q_rec_hot/1000:>8.3f} MW "
          f"(冷却側回収 {Q_rec_cold/1000:.3f} MW)")
    print(f"    残 加熱 Q_H_min  : {hr.Q_H_min_kW/1000:>8.3f} MW (外部熱媒で供給)")
    print(f"    残 冷却 Q_C_min  : {hr.Q_C_min_kW/1000:>8.3f} MW (外部冷媒で除去)")
    if not hr.feasible:
        print(f"    ⚠ HI infeasible: {hr.message}")

    # tier 別配分 (Q [MW] と費用 [億円/年])
    bd_Q    = hr.utility_breakdown
    bd_cost = opex_calc.get('breakdown', {})
    if bd_Q:
        print()
        print(f"    {'utility tier':<24} {'Q [MW]':>9} {'費用 [億円/年]':>15}")
        print(f"    {'-'*24} {'-'*9} {'-'*15}")
        # 加熱 tier (is_heating) を先、冷却 tier を後で、各々 Q 降順
        heat_names = {t.name for t in heating_tiers}
        heat_rows = [(n, q) for n, q in bd_Q.items() if n in heat_names]
        cool_rows = [(n, q) for n, q in bd_Q.items() if n not in heat_names]
        for group_label, rows in (('加熱', heat_rows), ('冷却', cool_rows)):
            for name, q_kW in sorted(rows, key=lambda kv: -kv[1]):
                cost = bd_cost.get(name, 0.0)
                print(f"    {('['+group_label+'] '+name):<24} "
                      f"{q_kW/1000:>9.3f} {cost:>15.4f}")
        print(f"    {'-'*24} {'-'*9} {'-'*15}")
        print(f"    {'(熱系 OPEX 合計, HI後)':<24} {'':>9} {opex_calc.get('total', 0.0):>15.4f}")
        unmatched = opex_calc.get('unmatched', {})
        if unmatched:
            print(f"    ⚠ tier 未マッチ (費用未計上): "
                  f"{', '.join(f'{k}={v/1000:.3f}MW' for k, v in unmatched.items())}")


# 完全酸化物基準エンタルピー計算の共有ヘルパ (show_oxide_enthalpy_flows /
# show_compression_detail / show_recycle_mixing で再利用)。datum = CO2+液H2O @298.15K。
_OXIDE_TREF_K = 298.15
_ox_thermo = None


def _oxide_enthalpy_parts(F_in: dict, T_K: float):
    """組成 F_in [kmol/h]・温度 T_K の完全酸化物基準エンタルピー流量。

    Returns (化学[GJ/h], 顕熱[GJ/h], 合計[GJ/h])。
      化学 = Σ n_i × HHV_i              (datum = CO2 + 液H2O @ 298.15K)
      顕熱 = Σ n_i × ∫_{298.15K}^{T} Cp_i dT
    """
    global _ox_thermo
    from src.cost_parameters import HHV_MJ_PER_KMOL
    if _ox_thermo is None:
        from src.thermo import PDHThermo
        _ox_thermo = PDHThermo()
    chem = sens = 0.0
    for k, n_kmolh in F_in.items():
        n = float(n_kmolh or 0.0)
        if n <= 0:
            continue
        chem += n * HHV_MJ_PER_KMOL.get(k, 0.0) / 1000.0                              # GJ/h
        try:
            sens += n * _ox_thermo.calc_enthalpy_change(k, _OXIDE_TREF_K, T_K) / 1.0e6  # GJ/h
        except KeyError:
            pass
    return chem, sens, chem + sens


def show_oxide_enthalpy_flows(result, design=None, config=None) -> None:
    """完全酸化物基準のエンタルピー流量を詳細表示 [GJ/h]。

    基準状態 (datum): 各元素の完全酸化物 (CO2 ガス + H2O 液) を 298.15 K で 0 とする
    「完全酸化物基準」。各成分のモルエンタルピーは

        h_i(T) = HHV_i + ∫_{298.15K}^{T} Cp_i dT   [J/mol]

    HHV_i (高位発熱量) = 完全燃焼で放出される熱 = 「酸化物を datum にしたときの化学
    エンタルピー」。これにストリーム温度までの顕熱を足したものが酸化物基準エンタルピー。
    ストリーム流量 H = Σ_i n_i × h_i を「化学 (HHV)」「顕熱」「合計」に分けて表示する。

    この基準では反応熱が HHV 差として自動的に含まれる (例: PDH 脱水素
    C3H8→C3H6+H2 は HHV 2220→2058+286=2344 kJ/mol で +124 kJ/mol = ΔH_rxn と一致、
    吸熱)。よって反応器・加熱炉を跨いだ系全体のエネルギー収支を Q = ΔH として一貫して
    追える。プラント境界 (Fresh 入 / 製品・オフガス・燃料 出) の収支も併記する。

    出典: HHV は cost_parameters.HHV_MJ_PER_KMOL (!仮置き 出典確認中、KNOWN_PLACEHOLDERS
    §A.3 参照)、Cp は THERMO_DATA の多項式 (化工便覧 改訂六版)。
    """
    if result.solver is None or result.solver.one_pass is None:
        return
    R = result.solver.one_pass
    _parts = _oxide_enthalpy_parts

    rmem = R['r_mem']
    rpsa = R['r_psa']
    # (グループ, ラベル, F_in[kmol/h], T_K)。
    # 非 ProcessStream の出口は組成 dict と代表温度を組み立てる。
    #   PSA 製品/オフガス: PSA 操作温度 T_abs=25°C (=298.15K) で出るため顕熱≈0 で近似。
    rows = [
        ('入口',    'Fresh LPG (Pump1出口)',        R['pump1'].outlet.F_in,   R['pump1'].outlet.T_in),
        ('内部',    'Reactor 入口 (Fresh+Recycle)', R['reactor_inlet'].F_in,  R['reactor_inlet'].T_in),
        ('内部',    'Reactor 出口',                 R['rx_out'].F_in,         R['rx_out'].T_in),
        ('内部',    'Dist2 塔頂 (→PSA)',            R['r2'].top.F_in,         R['r2'].top.T_in),
        ('内部',    'Dist2 塔底 (→Mem)',            R['r2'].bottom.F_in,      R['r2'].bottom.T_in),
        ('内部',    'Mem 透過 (→Dist3)',
            {'A': rmem.product.F_C3H8, 'B': rmem.product.F_C3H6},            rmem.product.T_out),
        ('出口',    'C3H6 製品 (Dist3塔頂)',        R['r3'].top.F_in,         R['r3'].top.T_in),
        ('出口',    'H2 製品 (PSA, ~25°C近似)',     rpsa.product,             298.15),
        ('出口',    'PSA オフガス (→燃料, ~25°C近似)', rpsa.offgas,           298.15),
        ('出口',    'Dist1 塔底 (→燃料)',           R['r1'].bottom.F_in,      R['r1'].bottom.T_in),
        ('Recycle', 'Dist3 塔底 → 反応器',          R['tear_dist3_new'],      R.get('T_d3_new', 298.15)),
        ('Recycle', 'Mem 保留 → 反応器',            R['tear_mem_new'],        R.get('T_mem_new', 298.15)),
    ]

    hdr("完全酸化物基準 エンタルピー流量 [GJ/h]  (datum: CO2 + 液H2O @ 298.15K = 0)")
    print("  h_i(T) = HHV_i + ∫Cp dT(298.15K→T)、H_flow = Σ n_i × h_i  (反応熱は HHV 差に内包)")
    print()
    print(f"    {'ストリーム':<26} {'T[°C]':>7} {'F[kmol/h]':>10} "
          f"{'H_化学[GJ/h]':>10} {'H_顕熱[GJ/h]':>10} {'H_合計[GJ/h]':>11}")
    print(f"    {'-'*26} {'-'*7} {'-'*10} {'-'*10} {'-'*10} {'-'*11}")

    last_group = None
    for group, label, F_in, T_K in rows:
        if group != last_group:
            print(f"  ▼ {group}")
            last_group = group
        chem, sens, tot = _parts(F_in, T_K)
        F_tot = sum(float(v or 0.0) for v in F_in.values())
        print(f"    {label:<26} {T_K-273.15:>7.1f} {F_tot:>10.1f} "
              f"{chem:>10.2f} {sens:>+10.2f} {tot:>11.2f}")

    # ---- プラント境界収支 (Recycle は内部ループなので除外) ----
    H_in = _parts(R['pump1'].outlet.F_in, R['pump1'].outlet.T_in)
    out_streams = [
        ('C3H6 製品', R['r3'].top.F_in,  R['r3'].top.T_in),
        ('H2 製品',   rpsa.product,      298.15),
        ('オフガス',  rpsa.offgas,       298.15),
        ('Dist1塔底', R['r1'].bottom.F_in, R['r1'].bottom.T_in),
    ]
    out_chem = out_sens = out_tot = 0.0
    for _, F_in, T_K in out_streams:
        c, s, t = _parts(F_in, T_K)
        out_chem += c; out_sens += s; out_tot += t

    print()
    print(f"  ─ プラント境界 エンタルピー収支 (完全酸化物基準、Recycle 除く) ─")
    print(f"    {'':<22} {'H_化学':>10} {'H_顕熱':>9} {'H_合計':>11}  [GJ/h]")
    print(f"    {'入  (Fresh LPG)':<22} {H_in[0]:>10.2f} {H_in[1]:>+9.2f} {H_in[2]:>11.2f}")
    print(f"    {'出  (製品+H2+OG+塔底)':<22} {out_chem:>10.2f} {out_sens:>+9.2f} {out_tot:>11.2f}")
    print(f"    {'ΔH (出 − 入)':<22} {out_chem-H_in[0]:>+10.2f} "
          f"{out_sens-H_in[1]:>+9.2f} {out_tot-H_in[2]:>+11.2f}")
    print(f"      ・ΔH_化学 > 0 = 反応で化学エンタルピー増 (PDH 脱水素は吸熱、加熱炉が供給)")
    print(f"      ・ΔH_合計 = 系へ正味で投入される熱+仕事に相当 (Q_利用・圧縮仕事と対応)")
    print(f"      ・HHV は !仮置き 値のため絶対量は暫定 (相対比較・収支検算に使用)")


def show_process_stream_table(result, design, config=None) -> None:
    """全ノード材料ストリーム表: 状態が変わる全ての点を T・P・組成・相・エンタルピー流量で記録。

    フローシートを工程順に辿り、状態変化が起きる**すべてのノード**を 1 行ずつ出す
    (HYSYS の material stream table 相当)。特に Dist2 塔底→膜 の区間は JT 減圧・気化・
    フィード圧縮機・膜・製品圧縮機・製品冷却器と状態変化が連続するため、その中間状態も
    すべて展開する。膜内部の圧縮機段の中間状態は equipment が保持する温度から復元する。

    列: ノード / 相 / T[°C] / P[bar] / F合計 / 各成分流量[kmol/h] / H_合計[GJ/h] (完全酸化物基準)。
    相は工程位置から既知のものをラベル付け (* = 膨張弁の気相モデル前提。液/二相の可能性あり)。
    """
    if result.solver is None or result.solver.one_pass is None:
        return
    R = result.solver.one_pass
    ORDER = ['A', 'B', 'C', 'D', 'E', 'F', 'Z']
    SHORT = {'A': 'C3H8', 'B': 'C3H6', 'C': 'H2', 'D': 'C2H4',
             'E': 'CH4', 'F': 'C2H6', 'Z': 'C4H10'}

    rmem = R['r_mem']
    meq  = rmem.equipment
    mem  = design.mem
    P_psa = R['r2'].top.P_in   # PSA は Dist2 塔頂圧で運転 (製品/オフガスの代表圧)

    # 膜内部ノードの復元 (membrane_system は中間ストリームを返さないので equipment から再構成)
    #   フィード圧縮機 出口  : 組成=retentate+product(C3), T=retentate.T_out, P=P_H
    #   膜 非透過(retentate) : 組成=retentate,            T=retentate.T_out, P=P_H
    #   膜 透過(permeate)    : 組成=product,              T=retentate.T_out(等温), P=P_L
    #   製品圧縮機 出口      : 組成=product,              T=T_cond_in_K,     P=P_dist
    #   製品冷却器 出口(=製品): 組成=product,             T=product.T_out,   P=P_dist
    mem_feedcomp_F = {'A': rmem.retentate.F_C3H8 + rmem.product.F_C3H8,
                      'B': rmem.retentate.F_C3H6 + rmem.product.F_C3H6}
    mem_ret_F      = {'A': rmem.retentate.F_C3H8, 'B': rmem.retentate.F_C3H6}
    mem_perm_F     = {'A': rmem.product.F_C3H8,   'B': rmem.product.F_C3H6}

    # Fresh LPG (Pump1 入口) は config から復元。config 無ければ Pump1 出口で代用。
    if config is not None:
        fresh_node = ('原料', 'Fresh LPG (Pump1入口)', '液',
                      R['pump1'].outlet.F_in, config.feed.T_K, config.feed.P_Pa)
    else:
        fresh_node = ('原料', 'Fresh LPG (≈Pump1入口)', '液',
                      R['pump1'].outlet.F_in, R['pump1'].outlet.T_in, R['pump1'].outlet.P_in)

    # (セクション, ラベル, 相, F_in, T_K, P_Pa)
    nodes = [
        fresh_node,
        ('原料',    'Pump1 出口 (→Dist1)',      '液',   R['pump1'].outlet.F_in,   R['pump1'].outlet.T_in,  R['pump1'].outlet.P_in),
        ('Dist1',   'Dist1 塔頂 (C3)',          '液',   R['r1'].top.F_in,         R['r1'].top.T_in,        R['r1'].top.P_in),
        ('Dist1',   'Dist1 塔頂 膨張後 (→Rx)',  '気*',  R['dist1_top_rx'].F_in,   R['dist1_top_rx'].T_in,  R['dist1_top_rx'].P_in),
        ('Dist1',   'Dist1 塔底 (→燃料)',       '液',   R['r1'].bottom.F_in,      R['r1'].bottom.T_in,     R['r1'].bottom.P_in),
        ('Recycle', 'Dist3 塔底 recycle 膨張後','気*',  R['recycle_dist3'].F_in,  R['recycle_dist3'].T_in, R['recycle_dist3'].P_in),
        ('Recycle', 'Mem 保留 recycle 膨張後',  '気',   R['recycle_mem'].F_in,    R['recycle_mem'].T_in,   R['recycle_mem'].P_in),
        ('Reactor', 'Reactor 入口 (合流後)',    '気',   R['reactor_inlet'].F_in,  R['reactor_inlet'].T_in, R['reactor_inlet'].P_in),
        ('Reactor', 'Reactor 出口',             '気',   R['rx_out'].F_in,         R['rx_out'].T_in,        R['rx_out'].P_in),
        ('Comp2',   'Cooler 出口',              '気',   R['cooled'].outlet.F_in,    R['cooled'].outlet.T_in,    R['cooled'].outlet.P_in),
        ('Comp2',   'Comp2a 出口',              '気',   R['comp2a'].outlet.F_in,    R['comp2a'].outlet.T_in,    R['comp2a'].outlet.P_in),
        ('Comp2',   'Intercool 出口',           '気',   R['intercool'].outlet.F_in, R['intercool'].outlet.T_in, R['intercool'].outlet.P_in),
        ('Comp2',   'Comp2b 出口',              '気',   R['comp2b'].outlet.F_in,    R['comp2b'].outlet.T_in,    R['comp2b'].outlet.P_in),
        ('Comp2',   'Desuper 出口 (→Dist2)',    '気',   R['desuper'].outlet.F_in,   R['desuper'].outlet.T_in,   R['desuper'].outlet.P_in),
        ('Dist2',   'Dist2 塔頂 (→PSA)',        '気',   R['r2'].top.F_in,         R['r2'].top.T_in,        R['r2'].top.P_in),
        ('PSA',     'PSA 製品 H2 (~25°C近似)',  '気',   R['r_psa'].product,       298.15,                  P_psa),
        ('PSA',     'PSA オフガス (~25°C近似)', '気',   R['r_psa'].offgas,        298.15,                  P_psa),
        ('Dist2→Mem', 'Dist2 塔底',             '液',   R['r2'].bottom.F_in,      R['r2'].bottom.T_in,     R['r2'].bottom.P_in),
        ('Dist2→Mem', '塔底 JT減圧 (→P_H)',     '液',   R['mem_feed_letdown'].F_in, R['mem_feed_letdown'].T_in, R['mem_feed_letdown'].P_in),
        ('Dist2→Mem', 'MemPrecool 出口 (気化)', '気',   R['mem_precool'].outlet.F_in, R['mem_precool'].outlet.T_in, R['mem_precool'].outlet.P_in),
        ('Dist2→Mem', 'Mem フィード圧縮機 出口','気',   mem_feedcomp_F,           rmem.retentate.T_out,    mem.P_H),
        ('Dist2→Mem', 'Mem 非透過 (→Recycle)',  '気',   mem_ret_F,                rmem.retentate.T_out,    mem.P_H),
        ('Dist2→Mem', 'Mem 透過 (P_L)',         '気',   mem_perm_F,               rmem.retentate.T_out,    mem.P_L),
        ('Dist2→Mem', 'Mem 製品圧縮機 出口',    '気',   mem_perm_F,               meq.T_cond_in_K,         mem.P_dist),
        ('Dist2→Mem', 'Mem 製品冷却器 出口(→D3)','液',  mem_perm_F,               rmem.product.T_out,      mem.P_dist),
        ('Dist3',   'Dist3 塔頂 (C3H6 製品)',   '液',   R['r3'].top.F_in,         R['r3'].top.T_in,        R['r3'].top.P_in),
        ('Dist3',   'Dist3 塔底 (→Recycle)',    '液',   R['r3'].bottom.F_in,      R['r3'].bottom.T_in,     R['r3'].bottom.P_in),
    ]

    hdr("全ノード 材料ストリーム表 (状態変化を全て記録)  成分流量 [kmol/h]、H=完全酸化物基準 [GJ/h]")
    comp_hdr = ' '.join(f"{SHORT[k]:>7}" for k in ORDER)
    print(f"  {'ノード':<26} {'相':>3} {'T[°C]':>7} {'P[bar]':>7} {'F計':>8}  {comp_hdr}  {'H[GJ/h]':>9}")
    print(f"  {'-'*26} {'-'*3} {'-'*7} {'-'*7} {'-'*8}  {'-'*(8*len(ORDER))}  {'-'*9}")

    last_sec = None
    for sec, label, phase, F_in, T_K, P_Pa in nodes:
        if sec != last_sec:
            print(f"  ── {sec} " + "─" * max(0, 60 - len(sec)))
            last_sec = sec
        F_tot = sum(float(F_in.get(k, 0.0) or 0.0) for k in ORDER)
        comps = ' '.join(f"{float(F_in.get(k, 0.0) or 0.0):>7.1f}" for k in ORDER)
        _, _, h = _oxide_enthalpy_parts(F_in, T_K)
        print(f"  {label:<26} {phase:>3} {T_K-273.15:>7.1f} {P_Pa/1e5:>7.3f} {F_tot:>8.1f}  "
              f"{comps}  {h:>9.1f}")

    print(f"  {'-'*26}")
    print("  注) 相: * = 膨張弁を気相モデルで通過 (実際は液/二相で出る可能性、別途協議中)。")
    print("      PSA 製品/オフガスは出口温度を持たないため T_abs=25°C で近似。")
    print("      Mem 内部 (フィード/製品圧縮機・膜) は equipment 保持温度から復元 (中間段の per-stage 状態は非保持)。")


def show_compression_detail(result, design) -> None:
    """多段圧縮系の詳細: 各段前後の温度・圧力・エンタルピー流量と、段間/前後の熱・仕事。

    対象:
      ▼ Comp2 系列 (反応器出口 0.5bar → Dist2 供給): Cooler→Comp2a→Intercool→Comp2b→
        Desuper の各ノードの T/P/H_合計 と、各ユニットの Q(冷却) / W(圧縮)。
        段間冷却 (Intercool) を挟む 2 段圧縮なので「圧縮機の間及び前後」が一望できる。
      ▼ 膜 圧縮機: フィード圧縮機 (feed.P→P_H) と 製品圧縮機 (P_L→P_dist)。各々 in/out の
        T/P/H・W・段数、段間冷却 (多段時) の Q/A/温度域。膜内部の多段は集約値で表示
        (per-stage の中間状態はモデルが保持しないため n_stages と合計で記述)。

    エンタルピー流量は完全酸化物基準 [GJ/h]。組成一定区間では H_合計差 = 顕熱差 で、
    圧縮機は ΔH ≈ +W、冷却器は ΔH ≈ −Q の検算ができる (1 kW = 3.6e-3 GJ/h)。
    """
    if result.solver is None or result.solver.one_pass is None:
        return
    R = result.solver.one_pass

    hdr("多段圧縮系 詳細 (T・P・エンタルピー流量・熱/仕事)")

    # ---- Comp2 系列 (反応器出口 → Dist2) ----
    print("  ▼ Comp2 系列 (反応器出口 → Dist2 供給、段間冷却付き 2 段圧縮)")
    print(f"    {'ノード(出口)':<22} {'T[°C]':>7} {'P[bar]':>7} {'H_合計[GJ/h]':>13}   "
          f"{'直前ユニット':<13} {'Q/W[kW]':>10}")
    print(f"    {'-'*22} {'-'*7} {'-'*7} {'-'*13}   {'-'*13} {'-'*10}")

    def _comp_node(label, stream, unit, duty_kW):
        _, _, h = _oxide_enthalpy_parts(stream.F_in, stream.T_in)
        duty_s = '' if unit is None else f"{duty_kW:>+10.0f}"
        unit_s = '' if unit is None else unit
        print(f"    {label:<22} {stream.T_in-273.15:>7.1f} {stream.P_in/1e5:>7.3f} "
              f"{h:>13.2f}   {unit_s:<13} {duty_s}")

    _comp_node('Reactor 出口',         R['rx_out'],           None,          0.0)
    _comp_node('Cooler 出口',          R['cooled'].outlet,    'Cooler (冷)',   R['cooled'].equipment.Q_duty_kW)
    _comp_node('Comp2a 出口',          R['comp2a'].outlet,    'Comp2a (圧)',   R['comp2a'].equipment.W_kW)
    _comp_node('Intercool 出口',       R['intercool'].outlet, 'Intercool (冷)', R['intercool'].equipment.Q_duty_kW)
    _comp_node('Comp2b 出口',          R['comp2b'].outlet,    'Comp2b (圧)',   R['comp2b'].equipment.W_kW)
    _comp_node('Desuper 出口(→Dist2)', R['desuper'].outlet,   'Desuper (冷)',  R['desuper'].equipment.Q_duty_kW)
    W_comp2 = R['comp2a'].equipment.W_kW + R['comp2b'].equipment.W_kW
    P_in  = R['rx_out'].P_in
    P_mid = R['comp2a'].outlet.P_in
    P_out = R['desuper'].outlet.P_in
    print(f"    → Comp2 合計圧縮仕事 W = {W_comp2:,.0f} kW、総圧縮比 {P_out/P_in:.1f}:1 "
          f"(段間 {P_in/1e5:.2f}→{P_mid/1e5:.2f}→{P_out/1e5:.2f} bar)")

    # ---- 膜 圧縮機 ----
    rmem = R['r_mem']
    meq  = rmem.equipment
    mem  = design.mem
    mp_out = R['mem_precool'].outlet
    print()
    print("  ▼ 膜 圧縮機 (フィード圧縮機 + 製品圧縮機、多段時は段間冷却)")

    # フィード圧縮機: mem_precool 出口 (≈P_H へ JT 減圧済の C3 ガス) → P_H
    feed_F = {'A': rmem.retentate.F_C3H8 + rmem.product.F_C3H8,
              'B': rmem.retentate.F_C3H6 + rmem.product.F_C3H6}
    _, _, h_fin  = _oxide_enthalpy_parts(feed_F, mp_out.T_in)
    _, _, h_fout = _oxide_enthalpy_parts(feed_F, rmem.retentate.T_out)
    print(f"    [フィード圧縮機] {meq.n_stages_feed} 段  W={meq.W_feed_kW:,.0f} kW")
    print(f"      in : T={mp_out.T_in-273.15:>7.1f}°C  P={mp_out.P_in/1e5:>7.3f} bar  H={h_fin:>10.2f} GJ/h")
    print(f"      out: T={rmem.retentate.T_out-273.15:>7.1f}°C  P={mem.P_H/1e5:>7.3f} bar  H={h_fout:>10.2f} GJ/h")

    # 製品圧縮機: 透過ガス P_L (1atm) → P_dist (Dist3 圧)。膜は等温 → in T = フィード圧縮機 out T。
    prod_F = {'A': rmem.product.F_C3H8, 'B': rmem.product.F_C3H6}
    _, _, h_pin  = _oxide_enthalpy_parts(prod_F, rmem.retentate.T_out)
    _, _, h_pout = _oxide_enthalpy_parts(prod_F, meq.T_cond_in_K)
    print(f"    [製品圧縮機]   {meq.n_stages_prod} 段  W={meq.W_prod_kW:,.0f} kW")
    print(f"      in : T={rmem.retentate.T_out-273.15:>7.1f}°C  P={mem.P_L/1e5:>7.3f} bar  H={h_pin:>10.2f} GJ/h")
    print(f"      out: T={meq.T_cond_in_K-273.15:>7.1f}°C  P={mem.P_dist/1e5:>7.3f} bar  H={h_pout:>10.2f} GJ/h")

    if meq.Q_intercool_kW > 1e-6:
        print(f"    [段間冷却] Q={meq.Q_intercool_kW:,.0f} kW  A={meq.A_intercool_m2:.0f} m²  "
              f"T {meq.T_intercool_in_K-273.15:.1f}→{meq.T_intercool_out_K-273.15:.1f}°C (feed+prod 合計)")
    else:
        print(f"    [段間冷却] なし (単段)")


def show_recycle_mixing(result, design=None) -> None:
    """リサイクル合流前後の圧力・温度・エンタルピー流量。

    反応器入口 mixer の「前」(Dist1 塔頂 C3 + Dist3 塔底 recycle + Mem 保留 recycle、
    いずれも膨張弁で反応器圧 0.5bar に減圧済) と「後」(reactor_inlet) を
    T/P/F/H_合計 で対比する。エンタルピー流量は完全酸化物基準 [GJ/h]。mixer は断熱合流
    (ΔH≈0、温度は成分別 Cp 加重平均) なので ΔH≈0 が収支検算になる。
    """
    if result.solver is None or result.solver.one_pass is None:
        return
    R = result.solver.one_pass
    # recycle_dist3/recycle_mem は run_one_pass 成功時のみ result に格納される。
    rd3  = R.get('recycle_dist3')
    rmem = R.get('recycle_mem')
    if rd3 is None or rmem is None:
        return

    hdr("リサイクル合流 詳細 (反応器入口 mixer 前後: T・P・エンタルピー流量)")
    print(f"    {'ストリーム':<26} {'T[°C]':>7} {'P[bar]':>7} {'F[kmol/h]':>10} {'H_合計[GJ/h]':>13}")
    print(f"    {'-'*26} {'-'*7} {'-'*7} {'-'*10} {'-'*13}")

    def _row(label, s):
        _, _, h = _oxide_enthalpy_parts(s.F_in, s.T_in)
        F = sum(float(v or 0.0) for v in s.F_in.values())
        print(f"    {label:<26} {s.T_in-273.15:>7.1f} {s.P_in/1e5:>7.3f} {F:>10.1f} {h:>13.2f}")

    print("  ▼ 合流前 (各々 膨張弁で反応器圧 0.5bar に減圧済)")
    before = [
        ('Dist1 塔頂 C3 (Fresh 分)', R['dist1_top_rx']),
        ('Dist3 塔底 recycle',       rd3),
        ('Mem 保留 recycle',         rmem),
    ]
    for label, s in before:
        _row(label, s)
    h_before = sum(_oxide_enthalpy_parts(s.F_in, s.T_in)[2] for _, s in before)
    F_before = sum(sum(float(v or 0.0) for v in s.F_in.values()) for _, s in before)
    print(f"    {'(合流前 合計)':<26} {'':>7} {'':>7} {F_before:>10.1f} {h_before:>13.2f}")

    print("  ▼ 合流後")
    _row('Reactor 入口 (mixed)', R['reactor_inlet'])
    _, _, h_after = _oxide_enthalpy_parts(R['reactor_inlet'].F_in, R['reactor_inlet'].T_in)
    print(f"    ΔH (後 − 前) = {h_after - h_before:+.3f} GJ/h  "
          f"(mixer 断熱合流のため理論上 ≈ 0、組成も保存)")


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
    # effective_TAC は純 TAC 最小化 (eff_econ.TAC + soft_penalty) なので、penalty は
    # effective_TAC と TAC の単純差分で求まる。
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
    # 全ノード材料ストリーム表 (状態変化を全て記録: T/P/組成/相/エンタルピー流量)
    show_process_stream_table(result, design, config)
    show_unit_details(R, design)
    show_production(R, F_C3H8_feed, F_C4H10_feed, config)
    show_capex(result.economics)
    show_opex(result.economics)
    show_revenue(result.economics)
    show_specs(result.specs, result.failure_reason, config)
    # HI と Stage 2 の物理メタデータ (pinch / HEN 構成) を先に表示
    show_hi_summary(result)
    # 細かい熱収支 (プロセスストリーム別の顕熱/潜熱 + ユーティリティ tier 別配分)
    show_heat_balance(result, design)
    # 完全酸化物基準のエンタルピー流量 (化学(HHV)+顕熱、反応熱込みの系全体収支)
    show_oxide_enthalpy_flows(result, design, config)
    # 多段圧縮系 (Comp2 系列・膜圧縮機) の段前後・段間の T/P/H/熱/仕事
    show_compression_detail(result, design)
    # リサイクル合流前後の T/P/エンタルピー流量
    show_recycle_mixing(result, design)
    show_stage2_synthesis(result)
    # raw/HI/Stage 2 を 1 表で比較 (中間サマリ)
    show_tac_summary(result, C3H6_product)
    # ★ 最終の目立つサマリボックス (一番下に必ず表示)
    show_final_summary_box(result, C3H6_product)
