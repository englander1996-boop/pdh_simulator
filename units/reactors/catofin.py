r"""units/reactors/catofin.py — Catofin型 浅床・大断面・多基並列スイング固定床 PDH 反応器。

設計思想 (商用 Catofin 準拠):
  - 軸流「深床」は 0.5 bar 低圧で Ergun 圧損・空塔速度が成立しない。捨てるのは軸流ではなく
    "深床"。→ **浅床 (L_bed 0.3-1.0 m)** にして圧損を桁で下げる。
  - 大流量を低空塔速度で通すため **N_online 基を並列** (= 触媒体積分割でなく「断面積確保数」、
    明示的設計変数)。各基は F_total / N_online を処理する。
  - 1基を軸流 ODE で解き、出口を N_online 倍する (F_1 = F_total/N_online → 解 → ×N_online)。
    velocity も ΔP も per-vessel になり、N_online を増やすほど低速・低圧損・高転化になる
    (Catofin の本質的トレードオフを BO が探索可能)。
  - 熱管理は **HGM 等価補償 (案A)**: 床温が T_in-ΔT_max 以下に下がろうとしたら HGM/再生蓄熱が
    反応吸熱を相殺して床温を維持する等価モデル。HGM 量は陽に解かず ΔT_max で表現 (一次近似)。
    供給熱 Q_HGM はエネルギー収支から後計算し燃料 OPEX に計上 ("無料の熱"にしない)。HGM が床
    体積を占める分は φ_cat<1 で有効触媒を減らす。
  - **設計思想 = 低速・大断面・多基**: 0.5bar 大流量では、空塔速度を下げる(u_design_min)ほど
    接触時間 τ=L/u が伸び (床ΔP∝uL 一定なら τ∝1/u²)、転化率↑・床ΔP↓。低速化で余った床ΔP
    予算は L_bed に戻せる (L↔N_online 同時最適化)。代償は基数増と流量分配の難化。

実機リスク (本モデル未計上、レポート明記事項):
  - 低速ほど床が整流抵抗として効かず流量分配がシビア → 分配板/マニホールド設計が主役 (偏差±5-10%)
  - 多基(数十基)スイング: 高温切替弁・再生空気/パージ弁・マニホールド・計装が支配的 (K_swing 感度)
  - 再生時ホットスポット: T_regen,max < T_limit を詳細設計制約とする (HGM 補償が大きいほど再生熱大)
  - 粒内拡散(Weisz-Prater)・外部物質移動(Mears): 3mm 粒子・低速で要確認 (低速で外部境膜厚化)
  - 圧損は **2段制約**: 床 Ergun ΔP/P_in < 5%、総 ΔP (床×K_internals)/P_in < 10%
    (内部品=分配板/支持グリッド/スクリーン/高温弁/ノズル/マニホールドを K_internals で一括計上)。
  - 総基数 **N_total = N_online × N_swing** (反応・再生・パージをカバーするセット数)。

物理 (反応速度・失活・断熱・Ergun ODE) は swing.py を per-vessel 流量で流用する。
swing._ode_axial は velocity を N_parallel で割らない (u = Q_vol/A_cross) ため、per-vessel
流量を渡すだけで SV・ΔP が自動的に per-vessel になる。
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

import os

from src.cost_calculator import calc_reactor_capex_okuyen
from units.reactors.swing import (
    DesignVars as _SwingDesign, FeedStream, FixedParams,
    SimulationResult, EffluentStream, EquipmentCost, PerformanceMetrics,
    _simulate_one_time, _penalty_result, _PENALTY_CAPEX, _COMPS, _thermo,
    THERMO_DATA, _T_REF,
)

# 床単独 ΔP/P_in の上限 (浅床なら容易に満たすが、0.5 bar abs では小さな圧損も効く)
_DP_BED_MAX = 0.05
# 総 ΔP = 床ΔP × K_internals。内部品 (分配板/支持グリッド/スクリーン/高温弁/ノズル/
# マニホールド/再生・パージライン) ぶんを一括計上。swing 既定 1.4 より厳しめの 2.0。
_K_INTERNALS = 2.0

# ---- HGM/再生蓄熱 熱補償モデル (案A) のパラメータ ----
# ΔT_max: 床内の許容温度降下 [K]。HGM/再生蓄熱が反応吸熱を相殺し、床温を T_in-ΔT_max 以上に
#   保つと仮定する等価モデル。HGM 量を陽に解かず「温度降下上限」で表現 (防御的な一次近似)。
#   感度: 30(楽観)/50(標準)/80(保守)。env PDH_CATOFIN_DTMAX で切替。
_DT_MAX_K = float(os.environ.get('PDH_CATOFIN_DTMAX', '50'))
# φ_cat: 有効触媒分率 [-]。HGM/蓄熱材が床体積の一部を占めるため、有効触媒は減る。
#   反応速度・転化率に反映 (φ_cat<1 で転化↓)。感度 0.80/0.85/0.90。env PDH_CATOFIN_PHI_CAT。
_PHI_CAT = float(os.environ.get('PDH_CATOFIN_PHI_CAT', '0.85'))

# D_eff: 触媒粒子内の有効拡散係数 [m²/s]。粒内拡散 η (Thiele) の計算に使う。多孔質・高温・
#   低圧で不確かさ大 (Knudsen/屈曲度/コーク)。標準 1e-5。感度 3e-6(保守)/1e-5/3e-5。env で上書き。
_D_EFF = float(os.environ.get('PDH_CATOFIN_DEFF', '1e-5'))
# 分配板/マニホールドの最小圧損 [bar]。低速・大粒径で床ΔPが小さくなりすぎると流量分配が
#   悪化するため、内部品ΔPに下限を設ける (実機は流量均一化のため意図的に圧損をつける)。
#   標準 0.01 bar (0.5bar運転で ~2%)。感度 0.005/0.015。env PDH_CATOFIN_DPDIST_MIN_BAR。
_DP_DIST_MIN_PA = float(os.environ.get('PDH_CATOFIN_DPDIST_MIN_BAR', '0.01')) * 1e5

# u_design_min: 空塔速度(superficial velocity)の設計下限 [m/s]。※ "SV"(空間速度 GHSV/WHSV)
#   ではなく空塔速度。Catofin の設計思想 = 低圧大流量を「大断面・低速・多基」で通すこと。
#   低速ほど接触時間(W/F)↑で転化↑・床ΔP↓だが、流量分配がシビアになる(分配板設計課題)。
#   標準 0.15。感度 0.10(技術上限)/0.20/0.30。env PDH_CATOFIN_UMIN。
#   注: フルイダイゼーション等の明確な物理下限ではなく、分配不良/偏流/外部物質移動低下を
#   避けるための設計判断値。0.10 は感度ケース扱い。
_U_DESIGN_MIN = float(os.environ.get('PDH_CATOFIN_UMIN', '0.15'))


def _enthalpy_J_per_mol(comp: str, T: float) -> float:
    """成分 comp の温度 T [K] における生成エンタルピー込みモルエンタルピー [J/mol]。
    h_i(T) = ΔHf_298 + ∫Cp dT (T_ref→T)。エネルギー収支で Q_HGM を求めるのに使う。"""
    return THERMO_DATA[comp].dHf_298 + _thermo.calc_enthalpy_change(comp, _T_REF, T)


@dataclass
class CatofinDesignVars:
    """Catofin型 浅床軸流スイング反応器の設計変数 (6個)。"""
    T_in:     float   # 反応器入口温度 [K]
    t_cyc:    float   # 反応フェーズ運転時間 [min]
    D:        float   # 1基あたり反応器内径 [m]
    L_bed:    float   # 浅床厚み [m] (0.3-1.0 を想定)
    N_online: int     # 反応中の並列基数 (大流量を低速で通すための断面積確保数)
    d_p:      float = 0.003   # 触媒粒径 [m] (2-6mm)。d_p↑ で ΔP↓ だが粒内拡散 η↓ のトレードオフ


def _catofin_fixed(d_p: float = 0.003) -> FixedParams:
    """Catofin 用 FixedParams。総ΔPマージンを K_internals に、空塔速度下限を低速運転用に下げ、
    触媒粒径を設計変数 d_p にする。

    FixedParams.SV_min_m_per_s は名前に反して「空塔速度 [m/s] の下限」(GHSV ではない)。
    Catofin は大断面・低速・多基運転を採るため、軸流深床用の 0.5 から u_design_min に下げる。
    """
    return FixedParams(dP_margin_factor=_K_INTERNALS, SV_min_m_per_s=_U_DESIGN_MIN,
                       d_p_m=float(d_p))


def simulate_catofin_reactor_system(
    design: CatofinDesignVars,
    feed: FeedStream,
    fixed: FixedParams | None = None,
    n_time_samples: int = 20,
) -> SimulationResult:
    """Catofin型 浅床・多基並列スイング反応器をシミュレーションする。

    1基あたり F_total/N_online を軸流 ODE (swing._simulate_one_time) で解き、出口を
    N_online 倍する。SV・ΔP は per-vessel で評価し、CAPEX/触媒量は N_total=N_online×N_swing。
    """
    if fixed is None:
        fixed = _catofin_fixed(getattr(design, 'd_p', 0.003))

    N = int(design.N_online)

    # ---- 入力バリデーション ----
    if design.t_cyc <= 0 or design.L_bed <= 0 or design.D <= 0 or N < 1:
        return _penalty_result(reason='input_invalid')
    if design.T_in <= 0 or feed.T_feed <= 0 or feed.P_in <= 0:
        return _penalty_result(reason='input_invalid')
    if any(v < 0 for v in feed.F_in.values()):
        return _penalty_result(reason='input_invalid')
    F_total_in = sum(feed.F_in.values())
    if F_total_in <= 0:
        return _penalty_result(reason='input_invalid')

    # ---- 1基あたり触媒量上限 (コンテスト §3-3: 最大 200 m³/基) ----
    # 深床許可 (L_bed≤5m) で V_cat/基 が 200m³ を超えうるため、ODE 前に早期厳守する。
    #   旧軸流 swing は V_cat_max_per_vessel で並列分割していたが、catofin は N_online を
    #   自由変数化した際にこの上限を落としていた。再導入して非現実な巨大単基を排除する。
    #   超過は _compute_reactor_shortfall の else 経由で reactor_other_shortfall=1.0 となり、
    #   BO に「この反応器幾何は不可」と伝わる (D↓/L_bed↓ 方向)。
    _A_cross_chk = math.pi / 4.0 * design.D ** 2
    _V_cat_chk = _A_cross_chk * design.L_bed * (1.0 - fixed.eps)
    if _V_cat_chk > fixed.V_cat_max_per_vessel:
        warnings.warn(
            f"catofin: V_cat/基={_V_cat_chk:.1f} m³ > 上限 "
            f"{fixed.V_cat_max_per_vessel:.0f} m³ "
            f"(D={design.D:.2f}m, L_bed={design.L_bed:.2f}m) — infeasible",
            UserWarning, stacklevel=2)
        return _penalty_result(reason='vcat_above_max')

    # ---- per-vessel フィード (流量を N_online で分割) ----
    feed_1 = FeedStream(
        F_in={c: feed.F_in.get(c, 0.0) / N for c in feed.F_in},
        T_feed=feed.T_feed,
        P_in=feed.P_in,
    )
    # 浅床軸流 = z_cat を L_bed に置いた swing 設計。N_parallel=1 (per-vessel) で解く。
    swing_design = _SwingDesign(
        T_in=design.T_in, z_cat=design.L_bed, t_cyc=design.t_cyc, D=design.D,
    )

    # ---- 時間方向サンプリング + 空間積分 (per-vessel) ----
    t_samples = np.linspace(0.0, design.t_cyc, n_time_samples)
    F_list, T_list, P_list = [], [], []
    # HGM 等価熱補償: 床温を T_in - ΔT_max 以上に保つ (t_floor_K)。HGM 占有で有効触媒 φ_cat。
    _t_floor = design.T_in - _DT_MAX_K
    for t in t_samples:
        F_out, T_out, P_out = _simulate_one_time(
            swing_design, feed_1, fixed, float(t), 1,
            t_floor_K=_t_floor, phi_cat=_PHI_CAT, d_eff=_D_EFF)
        if F_out is None:
            return _penalty_result(reason='sim_failure')
        F_list.append(F_out)
        T_list.append(T_out)
        P_list.append(P_out)

    F_arr = np.array(F_list)   # (n, 6) [mol/s] per-vessel
    T_arr = np.array(T_list)
    P_arr = np.array(P_list)
    _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
    F_avg_mol_s = _trapz(F_arr, t_samples, axis=0) / design.t_cyc   # per-vessel 時間平均
    T_out_avg = float(_trapz(T_arr, t_samples) / design.t_cyc)
    P_out_bed = float(_trapz(P_arr, t_samples) / design.t_cyc)

    # ---- 2段 ΔP 制約 (床 < 5% / 総 < 10%) ----
    # 総ΔP = 床ΔP + 内部品ΔP。内部品ΔP = max((K_internals-1)·床ΔP, 分配最小ΔP)。
    #   低速・大粒径で床ΔPが小さくなりすぎても、分配板/マニホールドには流量均一化のための
    #   最小圧損 (_DP_DIST_MIN_PA) が要る。これを入れないと BO が "超低ΔP" を選びすぎる。
    dP_bed_pa = max(feed.P_in - P_out_bed, 0.0)
    dP_internals_pa = max((fixed.dP_margin_factor - 1.0) * dP_bed_pa, _DP_DIST_MIN_PA)
    dP_total_pa = dP_bed_pa + dP_internals_pa
    dP_bed = dP_bed_pa / feed.P_in if feed.P_in > 0 else 0.0
    dP_total = float(np.clip(dP_total_pa / feed.P_in, 0.0, 1.0)) if feed.P_in > 0 else 0.0
    P_out_total = feed.P_in * (1.0 - dP_total)   # 内部品込み出口圧 (→下流圧縮機へ伝播)
    if dP_bed > _DP_BED_MAX:
        warnings.warn(
            f"catofin: 床ΔP/P_in={dP_bed*100:.1f}% > {_DP_BED_MAX*100:.0f}% "
            f"(L_bed={design.L_bed:.2f}m, D={design.D:.2f}m, N_online={N}) — infeasible",
            UserWarning, stacklevel=2)
        return _penalty_result(reason='dP_excess', dP_over_P=dP_total)
    if dP_total > fixed.dP_over_P_max:
        warnings.warn(
            f"catofin: 総ΔP/P_in={dP_total*100:.1f}% (床×{fixed.dP_margin_factor}) > "
            f"{fixed.dP_over_P_max*100:.0f}% — infeasible",
            UserWarning, stacklevel=2)
        return _penalty_result(reason='dP_excess', dP_over_P=dP_total)

    # ---- 空塔速度 (per-vessel) ----
    A_cross = math.pi / 4.0 * design.D ** 2
    n_inlet_mol_s = (F_total_in / N) * 1000.0 / 3600.0   # per-vessel [mol/s]
    Q_vol = n_inlet_mol_s * 8.314 * design.T_in / feed.P_in if feed.P_in > 0 else 0.0
    SV = Q_vol / A_cross if A_cross > 0 else 0.0
    if not (fixed.SV_min_m_per_s <= SV <= fixed.SV_max_m_per_s):
        warnings.warn(
            f"catofin: SV={SV:.2f} m/s が範囲 [{fixed.SV_min_m_per_s}, {fixed.SV_max_m_per_s}] 外 "
            f"(D={design.D:.2f}m, N_online={N}, T_in={design.T_in:.1f}K) — infeasible",
            UserWarning, stacklevel=2)
        return _penalty_result(reason='sv_out_of_range', SV_actual=SV)

    # ---- 出口を N_online 倍してシステム全体出口に ----
    F_out_total_kmolh = {
        c: float(F_avg_mol_s[i]) * 3600.0 / 1000.0 * N for i, c in enumerate(_COMPS)
    }

    # ---- 入口予熱 Q_preheat [GJ/h] (全流量, T_feed→T_in) ----
    q_w = 0.0
    for comp in _COMPS:
        F_mol_s = feed.F_in.get(comp, 0.0) * 1000.0 / 3600.0
        q_w += F_mol_s * _thermo.calc_enthalpy_change(comp, feed.T_feed, design.T_in)
    Q_preheat_GJh = q_w * 3600.0 / 1e9

    # ---- HGM 補償熱 Q_HGM [GJ/h] をエネルギー収支から後計算し、予熱炉燃料 OPEX に加算 ----
    # 床のエネルギー収支: Σ F_in·h(T_in) + Q_HGM = Σ F_out·h(T_out)  (h は生成エンタルピー込み)。
    #   → Q_HGM = Σ F_out·h(T_out) − Σ F_in·h(T_in)。完全断熱なら ~0、HGM が床温を維持した分だけ
    #   正になる。これを Q_preheat に足して "無料の熱" にしない (= 再生/燃料コストとして TAC 反映)。
    #   商用 Catofin では HGM 酸化還元/コーク燃焼/再生で供給される熱に対応 (本モデルは燃料等価)。
    Q_hgm_W = 0.0
    for comp in _COMPS:
        f_out_mol_s = F_out_total_kmolh.get(comp, 0.0) * 1000.0 / 3600.0
        f_in_mol_s = feed.F_in.get(comp, 0.0) * 1000.0 / 3600.0
        Q_hgm_W += (f_out_mol_s * _enthalpy_J_per_mol(comp, T_out_avg)
                    - f_in_mol_s * _enthalpy_J_per_mol(comp, design.T_in))
    Q_hgm_GJh = max(Q_hgm_W, 0.0) * 3600.0 / 1e9
    Q_preheat_GJh += Q_hgm_GJh

    # ---- 装置・基数 (共有プール フリート・サイジング) ----
    # 各反応器は固定床のまま (触媒は動かない=移動床ではない)、t_cyc 反応 → t_off オフライン
    # (再生+パージ+切替+再加熱) を位相をずらして循環する。定常で N_online 基を常時オンラインに
    # 保つ総数は時間占有率から:
    #   N_total = ceil( N_online × (t_cyc + t_off) / t_cyc ) = ceil( N_online × (1 + t_off/t_cyc) )
    # 旧「ブロック複製式」 N_online × (ceil(t_regen/t_cyc)+1) は N_online 基を1ブロックとして
    # 一斉切替する保守近似で、per-set の切り上げを ×N_online 増幅し過大計上していた
    # (例: N_online=24, t_cyc=14, t_off=30 → 旧96 基)。実機は1基/小グループ単位で位相を
    # ずらせるので、切り上げは総数に対し最後に1回だけ。
    #
    # ★t_off = t_cyc を採用 (balanced cycle): 実機 Catofin は反応と再生をほぼ同尺
    #   (各~10-15min) でバランスさせるため、t_off (再生+パージ+切替+再加熱) ≈ t_cyc と置く。
    #   → fleet 倍率 = 1 + t_off/t_cyc = 2.0 (t_cyc の値によらず)、N_total = 2 × N_online。
    #   旧 !仮置き t_regen=30min は反応の2倍超で倍率3.14 となり台数を過大にしていた。
    #   (fixed.t_regen は swing/radial 用に残置、catofin では未使用。)
    _t_off = design.t_cyc                          # balanced: 再生≈反応 (実機 Catofin)
    _fleet_mult = 1.0 + _t_off / design.t_cyc      # = 2.0
    N_total = math.ceil(N * _fleet_mult)
    N_swing_sets = math.ceil(_fleet_mult)          # = 2 (表示用: (反応+休止)/反応 倍率)
    V_cat_pervessel = A_cross * design.L_bed * (1.0 - fixed.eps)   # 1基床体積 [m³]
    V_vessel_pervessel = A_cross * design.L_bed                    # 1基容器体積 [m³]
    catalyst_weight_total = V_cat_pervessel * N_total * fixed.rho_b  # [kg]
    if V_vessel_pervessel <= 0:
        return _penalty_result(reason='volume_zero')
    try:
        reactor_capex = calc_reactor_capex_okuyen(
            V_vessel_m3=V_vessel_pervessel, P_abs_pa=feed.P_in,
            D_m=design.D, N_reactors_total=N_total,
        )
    except Exception:
        reactor_capex = _PENALTY_CAPEX

    # ---- 性能 (intensive: per-vessel = system 全体で同一) ----
    F_A_in = feed.F_in.get('A', 0.0); F_A_out = F_out_total_kmolh['A']
    F_B_in = feed.F_in.get('B', 0.0); F_B_out = F_out_total_kmolh['B']
    conversion = (F_A_in - F_A_out) / F_A_in * 100.0 if F_A_in > 0 else 0.0
    delta_A = F_A_in - F_A_out
    selectivity = (F_B_out - F_B_in) / delta_A * 100.0 if delta_A > 0 else 0.0
    conversion = float(np.clip(conversion, 0.0, 100.0))
    selectivity = float(np.clip(selectivity, 0.0, 100.0))

    return SimulationResult(
        effluent=EffluentStream(
            F_out_avg=F_out_total_kmolh,
            T_out_avg=T_out_avg,
            Q_preheat=Q_preheat_GJh,
            P_out=P_out_total,
        ),
        equipment=EquipmentCost(
            V_vessel_actual=V_vessel_pervessel,
            N_parallel=N,                 # = N_online (断面積確保数)
            N_swing_sets=N_swing_sets,
            N_reactors_total=N_total,
            Catalyst_Weight_Total=catalyst_weight_total,
            Reactor_CAPEX=reactor_capex,
            dP_over_P_actual=dP_total,
        ),
        performance=PerformanceMetrics(
            Conversion=conversion,
            Selectivity=selectivity,
        ),
    )
