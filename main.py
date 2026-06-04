r"""
main.py — HYSYS + SM バックエンドでの PDH プロセス全変数最適化 (制約付き Optuna BO)

exp3.py を評価関数として、**全設計変数**(反応器・PSA・膜・原料 + 蒸留塔 3 つ。反応器は
REACTOR_KIND で Catofin23変数/径方向流22変数/軸流21変数を切替、既定=Catofin23)を
Bayesian Optimization で探索する。**BO を実行する本丸ファイル**。

BO インフラ: QMC→TPE 2相サンプラ・constraints_func(連続制約)・penalty_scale・
  _store_diagnostics・compact callback(flush 付きライブログ, callbacks.py)・
  詳細レポート保存(display_full_results)。
構成: Dist1/Dist3=SM(学習済み GPR, ~瞬時)、Dist2=HYSYS。
  純度は SM Dist3 が 99.5 mol%=99.497 wt% 固定 → spec を 99.45 wt% に緩和。
  塔本体 CAPEX は provider 側で FUG と同式で計算(N/還流が CAPEX に効く)。

変数 (反応器は REACTOR_KIND で Catofin/径方向流/軸流を切替。既定=catofin):
  反応器(Catofin 6): T_in_K, t_cyc_min, D_reactor_m, L_bed_m, N_online, d_p_mm   (← 既定)
  反応器(径方向流 5): T_in_K, t_cyc_min, D_inner_m, bed_thickness_m, H_m   (棄却、比較用に残置)
  反応器(軸流 4)    : T_in_K, z_cat_m, t_cyc_min, D_reactor_m              (旧、比較用に残置)
  PSA(3)   : D_psa_col_m, L_psa_bed_m, desorption_target
  膜(2)    : P_H_Pa, A_mem_m2   (P_L=1atm 固定、mem.P_dist=Dist3圧 同期)
  原料(1)  : F_C3H8_fresh_kmol_h
  Dist1(4) : col1_p_kpa, col1_n_stages, col1_feed_stage, col1_comp_frac_2  (SM)
  Dist2(4) : col2_p_kpa, col2_n_stages, col2_feed_ratio, col2_reflux_ratio (HYSYS)
  Dist3(3) : col3_p_kpa, col3_n_stages, col3_feed_ratio                    (SM, spec なし)

出力 (run ごとの subdir へ集約):
  outputs/main_<ts>/README.md            : 結果の見方ガイド (最初に開く)
  outputs/main_<ts>/trials.csv           : 全 trial の params + 診断
  outputs/main_<ts>/best.json            : best trial 要約
  outputs/main_<ts>/top{1..N}_trial*.txt : 上位候補の詳細レポート (CAPEX/OPEX/spec/HI 内訳)
  outputs/main_<ts>/optuna.db            : Optuna SQLite (USE_SQLITE_STORAGE or 並列時のみ)
  stdout(リダイレクト推奨): compact callback による trial 毎ライブログ

使い方:  .\.venv\Scripts\python.exe main.py > outputs\main_run.log 2>&1
  (Python は 3.13。flush 付きログなのでリダイレクトでも live に書き出される。)
"""

import os
import sys
import csv
import json
import time
import datetime
import contextlib
import io
from typing import Optional

os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '300')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import optuna

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate, FlowsheetResult
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.reactors.radial_flow import RadialDesignVars
from units.reactors.catofin import CatofinDesignVars
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

from optimization.study import make_sampler, run_optimization, _default_constraints_func
from optimization.objective import _store_diagnostics
from optimization.penalty_scale import set_scale, default_schedule
from simulation import display_full_results, show_input_snapshot


# ===========================================================================
# § 1. BO 設定
# ===========================================================================
N_TRIALS    = int(os.environ.get('PDH_N_TRIALS', '300'))   # env で短縮可 (検証用)。既定 300
N_STARTUP   = int(os.environ.get('PDH_N_STARTUP', '50'))   # QMC Sobol 広域カバレッジ (以降 TPE)
# SEED/SAMPLER は env で上書き可能(既定 tpe/42)。単体起動は env 未設定なので既定で動く。
#   seed散らし/対照群バッチ (verification/run_seed_robustness.py) が run ごとに
#   PDH_SEED・PDH_SAMPLER('tpe' or 'random') を設定して切り替える。
#   狙い: 初期引きによらず同じ解に収束するか(再現性)と、BO が random 探索を
#   安定して上回るか(BO の正当性)の検証。
SEED        = int(os.environ.get('PDH_SEED', '42'))
SAMPLER     = os.environ.get('PDH_SAMPLER', 'tpe')   # 'tpe'(既定) | 'random'(学習なし対照群)
N_JOBS      = 1              # HYSYS COM + penalty_scale global のため 1 (スレッド並列は不可)
# マルチプロセス並列 worker 数。>1 で N プロセスが共有 SQLite study を分担し、各 worker が
# 自前の HYSYS インスタンスを起動する (重い: 各 ~数百MB + 起動 ~1分)。各 worker 単スレッドで
# penalty_scale/GIL 問題なし、constant_liar=True で冗長サンプリング抑制。1 で単一プロセス。
# 本 main (HYSYS+SM) は N_WORKERS=1 のみ正しい。HYSYS.Application はユーザーセッションに単一
# インスタンスの COM サーバーで、複数 worker が Dispatch しても同じ HYSYS を共有 →
# ケース切替(Close+Open)が跨プロセスで割り込み結果を汚染する。並列化できるのは HYSYS 不使用の
# sub1/sub2 のみ。本 main は SM 化で元々高速。
N_WORKERS   = 1
N_TOPK      = 3              # 詳細レポートを出す上位候補数

USE_SQLITE_STORAGE = False
STUDY_NAME         = "pdh_hysys_sm_main"


# ===========================================================================
# § 2. 固定値 / 評価オプション
# ===========================================================================
P_L_Pa       = 1.0e5         # 膜透過側圧力 (大気圧固定)
APPLY_HI     = True
# Stage2(HEN greedy 合成) を外し HI(targeting) のみで評価する。
#   (1) HI の経済的本体=OPEX 削減は Stage1 (Q_H_min/Q_C_min = MER 理論最少) が既に正確に捉える。
#       現 greedy Stage2 は under-recover (95%許容) で OPEX をむしろ過大評価する。
#   (2) Stage2 が追加するのは「回収網 CAPEX」= ΔTmin の energy-capital トレードオフの capital 側。
#       本構成は ΔTmin=10K 固定で最適化しないため、その CAPEX は ~一定オフセットで BO ランキングを
#       歪めず、省いても実害が小さい (回収 CAPEX を課す意味が出るのは ΔTmin を最適化する時のみ)。
#   (3) 目的関数が滑らかになり TPE/GP に優しい。greedy HEN (optimization/hen_synthesis.py) は残置・実害なし。
#   注: HI-only は回収網 CAPEX を計上しないターゲティング水準の評価。レポートにその旨を明記すること。
APPLY_STAGE2 = False
HI_DT_MIN_K  = 10.0


# ===========================================================================
# § 3. 探索範囲 — 反応器 REACTOR_KIND 依存 (軸流21/径方向流22)  形式: (low, high, scale, type)
#   bounds は forensic 調整値 ∩ SM 分類器 feasible 領域。
#   SM 崖 (予測不能域) は除外、解の縁 (収束ぎりぎり) は含める。
# ===========================================================================
# ===========================================================================
# 反応器モデル選択
#   'radial' = 径方向流 (薄い環状床を半径方向に通す。0.5bar 低圧でも圧損が桁で小さく、
#              現実的な 3mm 触媒で feasible。実機 Oleflex/Catofin の設計思想)。本既定。
#   'axial'  = 旧 軸流深床 (units/reactors/swing.py)。Ergun 圧損を入れると 0.5bar では
#              現探索域 z_cat 15-40m が全 infeasible になる (= 比較用に残置)。
#   詳細: monitor/reactor_pressure_drop_and_geometry.ipynb, units/reactors/SPEC_swing.md
# ===========================================================================
#   'catofin' = 浅床・大断面・多基並列スイング (units/reactors/catofin.py)。商用 Catofin 準拠。
#              深床軸流の圧損問題を「浅床(L_bed≤1m)+多基並列(N_online)」で回避。本既定。
#              捨てたのは"軸流"でなく"深床"。径方向流は内部品/分配/再生温度の不確かさで撤去。
REACTOR_KIND = 'catofin'   # 'catofin'(既定) | 'radial' | 'axial'

_REACTOR_SPACE = {
    # Catofin 型 浅床軸流スイング (既定)。浅床+多基並列で 0.5bar 低圧の圧損を成立させる。
    'catofin': {
        # T_in 下限を 930 に引上げ (旧 880)。catofin の per-pass 転化は低く、低 T_in では
        # クラッキング由来の C2 (CH4/C2H6/C2H4) 副生が乏しい (<3.5%) → 脱エタン塔 (Dist2) の
        # 塔頂が H2 偏重で −97〜−105℃ に冷え、エチレン冷媒 −100℃ で凝縮不能 (cold-top)。
        # T_in≥928 で C2≳4% になり塔頂 ≥−93℃ で凝縮可。930 で −88℃級の余裕を確保。
        # (検証: tools/_diag_feasible_search.py / _probe_dist2_levers_c2lean.py。低 R/低 N/高 col2_p
        #  では C2-lean フィードを暖められず、C2 生成=高 T_in が唯一の有効レバーと確認。)
        "T_in_K":              (930.0,  940.0,  'linear', 'float'),
        "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),  # 反応フェーズ時間
        "D_reactor_m":         (9.0,    12.0,   'linear', 'float'),  # 1基あたり内径 (下限9: N_online≤8では小径=流速過大で床ΔP死、探索条件付け)
        # L_bed と N_online は同時最適化が肝。基数上限を 8 に絞り (N_total=N_online×N_swing≤32)、
        # 旧解の「17m³級×96基」という非現実設計 (= 1基あたりコンテスト上限 200m³ の 8.6% しか
        # 使わず弁/マニホールド/計装が支配的になる) を排除する。少数基に必要触媒を収めるため
        # 浅床縛り (旧 L_bed≤1.0) を外して深床を許可 (上限 5.0m)。深床化による床ΔP増は
        # d_p↑ (粒径拡大) と既存 ΔP 制約 (床<5%/総<10%) のトレードオフに管理させる。
        # 1基あたり触媒量はコンテスト §3-3 の 200m³ で上限を厳守 (catofin.py 側)。
        "L_bed_m":             (0.3,    3.0,    'linear', 'float'),  # 床厚み (深床許可だが上限3.0: L>3は0.5barでΔP死域。D≥10で十分な触媒量)
        "N_online":            (6,      8,      'linear', 'int'  ),  # 並列基数 (≤8 → balanced で N_total=2×N_online≤16, 実機妥当性)
        "d_p_mm":              (4.0,    6.0,    'linear', 'float'),  # 触媒粒径 (下限4: 小粒径は床ΔP過大で死。4-6では転化率ほぼ不変=損失なし)
    },
    # 軸流深床 (旧、圧損で 0.5bar では成立しない。比較・回帰用)
    'axial': {
        "T_in_K":              (880.0,  940.0,  'linear', 'float'),
        "z_cat_m":             (15.0,   40.0,   'linear', 'float'),
        "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
        "D_reactor_m":         (7.0,    10.0,   'linear', 'float'),
    },
    # 径方向流 (圧損レビュー後の既定)。bounds は _smoke_test_ergun / 圧損 nb の feasible 帯から。
    'radial': {
        "T_in_K":              (880.0,  940.0,  'linear', 'float'),
        "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
        "D_inner_m":           (6.0,    12.0,   'linear', 'float'),   # 中心捕集管径 (r_i=3-6m)
        #   太い径=低空塔速度=低ΔP で可行域が広がる。容器CAPEX増は原料費比で軽微。
        #   下限 6: 多段(3段)化で累積ΔPが ~3倍になり、Di<6(小流路=高流速)は H/dr によらず
        #   ほぼ ΔP>10% で棄却される。Di≥6 で可行域に入る。
        "bed_thickness_m":     (0.3,    0.8,    'linear', 'float'),   # 環状床厚 Δr (薄い=低圧損)。
        #   ΔP ∝ 床厚。上限 0.8 で ΔP/P を 10% 制約内に抑える(床厚 1.4m で ΔP 16%)。
        #   触媒量は H↑/D_inner↑ で補う。
        "H_m":                 (22.0,   30.0,   'linear', 'float'),   # 床高 (触媒量を稼ぐ)
        #   下限 22: 多段(3段)化で累積ΔP≈3倍。H≥22・Di≥6 なら床厚0.8でも累積ΔP<0.10 に収まる。
    },
}[REACTOR_KIND]

SEARCH_SPACE = {
    # ----- 反応器 (REACTOR_KIND で軸流/径方向流を切替。上の _REACTOR_SPACE 参照) -----
    **_REACTOR_SPACE,

    # ----- PSA -----
    "D_psa_col_m":         (2.9,    5.0,    'linear', 'float'),
    "L_psa_bed_m":         (22.0,   30.0,   'linear', 'float'),
    "desorption_target":   (0.22,   0.40,   'linear', 'float'),

    # ----- 膜 (P_L 固定、P_dist は Dist3 圧と同期) -----
    "P_H_Pa":              (7.5e5,  9.5e5,  'linear', 'float'),
    "A_mem_m2":            (5.0e4,  3.0e5,  'log',    'float'),

    # ----- 原料 (外側ループ skip、override) -----
    # HYSYS+SM フローの実効収率は ~72-78%。中央収率で target 1188 kmol/h を満たすには
    #   F ≈ 1188/0.72 ≈ 1645 が必要。F_fresh は TAC の第1ドライバ。生産量下限 target×0.95
    #   =1128.6 kmol/h は best 収率78.4%なら F≈1440 でも充足可能なので、下限を 1450 まで開けて
    #   より安い解を探索する(下げ過ぎ=production_under 増を避けるため 1450 で留める)。
    "F_C3H8_fresh_kmol_h": (1450.0, 1750.0, 'linear', 'float'),

    # ----- Dist1 (SM model1: N30-60/P1600-2200/feed_stage10-39/CF0.9-0.999) -----
    # feed_stage は SM feas ≥22 (プローブ: <21 で 0%)。範囲は (22,28) 固定: N 下限 30 でも
    # N-2=28 以下で常に有効になり、動的 search space を避けて TPE(multivariate) をフル活用。
    # (SM feas 域 ≤39 のうち 22-28 を採用。Dist1 は cost driver でないため探索損失は許容。)
    "col1_p_kpa":          (1600.0, 2000.0, 'linear', 'float'),
    "col1_n_stages":       (30,     60,     'linear', 'int'  ),
    "col1_feed_stage":     (22,     28,     'linear', 'int'  ),
    "col1_comp_frac_2":    (0.90,   0.999,  'linear', 'float'),

    # ----- Dist2 (HYSYS 脱エタン塔)。収束 envelope 狭、縁を含む (N≈45 安/N=75-80 頑健) -----
    # 塔頂を浅冷化(高P→塔頂温度↑)して深冷コンデンサ費(エチレン-100C 17731円/GJ)を下げるレバー。
    # 低 P 側は Dist2 塔頂が cold-top (H2 希釈で -104℃ 級、エチレン-100℃ で凝縮不能) になりやすく
    #   feasible 0 の最多要因。col2_p↑ が塔頂を暖める方向。HSC は 950kPa まで実走確認済。
    #   col2_p > 膜 P_H のケースは膜前 JT let-down (run_one_pass) で P_H へ減圧して吸収する。
    # Catofin 低C2フィードでは Dist2 が分岐収束し、warm 物理解(塔頂-85〜-87℃)が安定するのは
    #   col2_p ≈ 740-800kPa の帯 (probe 確認: 全 N=70-80 で feasible)。P≥830 で時々 cold 分岐
    #   (-128/-150℃)に飛ぶため上限を 820 に絞る。※radial の "高P=暖かい" とは逆 (低C2のため)。
    #   反応器は 0.5bar 固定で col2_p と競合しない (Dist2 feed はどのみち圧縮)。
    # 下限を 770 に引上げ (旧 740)。T_in=930 (下限) では col2_p≤750 だと依然 cold-top
    #   (検証: _diag_feasible_search.py で 930/750=不可・930/820=可・935/780=可・940/760=可)。
    #   低 T_in × 低 col2_p の死に隅を削る。塔頂を暖める方向 (高 col2_p) に寄せる。
    "col2_p_kpa":          (770.0,  820.0,  'linear', 'float'),
    # Dist2 塔頂温度は N で強く決まり (N44→-125℃, N60→-101℃, N80→-98℃)、低 N は塔頂が冷えすぎて
    #   -100℃ エチレン冷媒でも凝縮不能 (cold-top)。feasible 帯は高 N 側 (N≥~70 で -98℃ 級)。
    #   HSC は 80 段まで (hysys_cases/column2)。
    "col2_n_stages":       (60,     80,     'linear', 'int'  ),
    "col2_feed_ratio":     (0.40,   0.60,   'linear', 'float'),
    # Dist2(HYSYS) 深冷コンデンサ(−83°C, エチレン-100C)が本 main TAC の最大コスト。
    #   還流比↓ = 凝縮 duty↓ = 直接削減だが、Dist2 が feasible になる還流比は R≥11
    #   (N=80 で R9.5→塔頂-98.1℃/dT1.9K=不可, R11→-96.1℃/dT3.9K=凝縮可)。下限を下回ると
    #   永遠に cold-top/不成立。低 R(≤9) は HYSYS 非収束(空出力)も多発。
    #   reboiler OPEX 増は BO が最小 R を探って均衡させる。
    "col2_reflux_ratio":   (10.0,   15.0,   'linear', 'float'),

    # ----- Dist3 (SM model3: N69-200/P1600-2200, spec なし)。feas: N≥115, P≤1900 -----
    "col3_p_kpa":          (1600.0, 1900.0, 'linear', 'float'),
    "col3_n_stages":       (115,    160,    'linear', 'int'  ),
    "col3_feed_ratio":     (0.60,   0.90,   'linear', 'float'),
}

# SM/HYSYS の feed_stage 絶対 bounds (ratio→絶対段変換後に clamp)。
_FEED_STAGE_ABS = {
    "col2": (2, 9999),    # HYSYS は N-2 のみ
    "col3": (70, 180),    # model3 feas 域 (feed_stage ≥70)
}


# ===========================================================================
# § 4. 出力
# ===========================================================================
OUTPUT_DIR = 'outputs'


# ===========================================================================
# ↓↓↓ 以下はパイプライン (通常触らない) ↓↓↓
# ===========================================================================

import dataclasses as _dc
_CONFIG = load_operating_config()
# SM Dist3 は 99.5 mol%=99.497 wt% 固定。mol↔wt 差 0.003pp を吸収するため
# purity 閾値を 99.45 wt% に緩和。
_CONFIG = _dc.replace(_CONFIG, spec=_dc.replace(_CONFIG.spec, c3h6_min_wtfrac=0.9945))


def _feed_stage_from_ratio(ratio: float, n: int, lo: int, hi: int) -> int:
    fs = int(round(ratio * n))
    hi_eff = min(hi, n - 2)
    lo_eff = min(lo, hi_eff)
    return max(lo_eff, min(fs, hi_eff))


def _suggest_params(trial: optuna.trial.Trial) -> dict:
    """SEARCH_SPACE 全変数を suggest して params dict を返す。

    全変数が固定範囲なので multivariate TPE をフル活用できる (動的 search space なし)。
    col2/col3 の feed は ratio で suggest し、_build_design で絶対段に変換 (N 依存 clamp)。
    """
    p: dict = {}
    for name, (low, high, scale, typ) in SEARCH_SPACE.items():
        if typ == 'int':
            p[name] = trial.suggest_int(name, int(low), int(high), log=(scale == 'log'))
        else:
            p[name] = trial.suggest_float(name, float(low), float(high), log=(scale == 'log'))
    return p


def _build_design(p: dict) -> FlowsheetDesignVars:
    """params dict (Catofin 23 / 径方向流 22 / 軸流 21 変数) から FlowsheetDesignVars を構築。trial 非依存 (best 再評価でも使用)。"""
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages']); fs2 = _feed_stage_from_ratio(p['col2_feed_ratio'], n2, *_FEED_STAGE_ABS['col2'])
    n3 = int(p['col3_n_stages']); fs3 = _feed_stage_from_ratio(p['col3_feed_ratio'], n3, *_FEED_STAGE_ABS['col3'])
    p3_kpa = float(p['col3_p_kpa'])
    # 反応器: REACTOR_KIND に応じて軸流 (SwingDesign) / 径方向流 (RadialDesignVars) を構築。
    # run_one_pass が型でディスパッチする。
    if REACTOR_KIND == 'catofin':
        reactor = CatofinDesignVars(T_in=p['T_in_K'], t_cyc=p['t_cyc_min'],
                                    D=p['D_reactor_m'], L_bed=p['L_bed_m'],
                                    N_online=int(p['N_online']),
                                    d_p=float(p['d_p_mm']) / 1000.0)
    elif REACTOR_KIND == 'radial':
        reactor = RadialDesignVars(T_in=p['T_in_K'], t_cyc=p['t_cyc_min'],
                                   D_inner=p['D_inner_m'], bed_thickness=p['bed_thickness_m'],
                                   H=p['H_m'])
    else:
        reactor = SwingDesign(T_in=p['T_in_K'], z_cat=p['z_cat_m'],
                              t_cyc=p['t_cyc_min'], D=p['D_reactor_m'])
    return FlowsheetDesignVars(
        swing=reactor,
        psa=PSADesignVars(D_col=p['D_psa_col_m'], L_bed=p['L_psa_bed_m'],
                          desorption_target=p['desorption_target']),
        mem=MemDesignVars(P_H=p['P_H_Pa'], P_L=P_L_Pa, A_mem=p['A_mem_m2'],
                          P_dist=p3_kpa * 1000.0),
        dist1=ColumnTunables(
            P_col=float(p['col1_p_kpa']) * 1000.0, N_stages=n1, N_feed=1, reflux_ratio=2.0,
            solver_method='sm', hysys_spec_value=float(p['col1_comp_frac_2']), hysys_feed_stage=fs1),
        dist2=ColumnTunables(
            P_col=float(p['col2_p_kpa']) * 1000.0, N_stages=n2, N_feed=1,
            reflux_ratio=float(p['col2_reflux_ratio']),
            solver_method='hysys', hysys_spec_value=float(p['col2_reflux_ratio']), hysys_feed_stage=fs2),
        dist3=ColumnTunables(
            P_col=p3_kpa * 1000.0, N_stages=n3, N_feed=1, reflux_ratio=12.0,
            solver_method='sm', hysys_spec_value=0.99, hysys_feed_stage=fs3),  # spec は SM 未使用
    )


def objective(trial: optuna.trial.Trial) -> float:
    scale = default_schedule(trial.number, N_TRIALS)
    set_scale(scale)
    trial.set_user_attr('penalty_scale', scale)

    params = _suggest_params(trial)
    design = _build_design(params)
    F_fresh = float(params['F_C3H8_fresh_kmol_h'])

    t0 = time.perf_counter()
    result: FlowsheetResult = evaluate(
        design, _CONFIG, verbose=False,
        apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K, apply_stage2=APPLY_STAGE2,
        F_C3H8_override=F_fresh,
    )
    trial.set_user_attr('wallclock_sec', time.perf_counter() - t0)
    _store_diagnostics(trial, result)
    trial.set_user_attr('F_C3H8_fresh_used_kmol_h', F_fresh)
    return result.effective_TAC


def constraints_func(trial: optuna.trial.FrozenTrial):
    """連続制約 (feas / 生産量方向 / 反応器SV / PSA / 膜 / 塔 shortfall)。

    全変数最適化なので反応器/PSA/膜の shortfall 信号も活性化し、TPE が上流の方向も学習する。
    純度は SM で不変のため制約化しない (定数=無意味)。
    """
    return _default_constraints_func(trial)


# ===========================================================================
# ライブログ用 compact callback (flush 付き)
# ===========================================================================
from collections import Counter as _Counter, deque as _deque
import time as _time
from optimization.callbacks import _fmt_dur, _fmt_reason_from_trial, _fmt_tally

_BAR_W = 30

# main は in-memory study で完走まで trials.csv を書かないため、run 中に失敗内訳を確認できるよう
# trial ごとに 1 行 JSONL を増分追記する (tools/_scan_run_log.py で mid-run 解析可)。
# 失敗は failure_unit より細かい「カテゴリ」に正規化して live tally + JSONL に出す (1 行/trial で軽量)。
_LIVE_JSONL = {'path': None}   # main() で run subdir のパスをセット


def _failure_category(failure_unit: str, failure_reason: str) -> str:
    """失敗を普遍的なカテゴリに正規化 (想定外エラーも必ずどれかに落ちる)。

    failure_unit (r1/r2/.../spec_*/timeout/exception:*) と failure_reason 文字列から、
    cold-top・HYSYS 空出力(COM)・反応器 ΔP/SV・spec 違反・例外・タイムアウト等に分類。
    """
    r = failure_reason or ''
    if not failure_unit or failure_unit == 'success':
        return 'success'
    if 'HysysEmptyOutput' in r or 'COM エラー' in r or 'empty' in r:
        return 'hysys_empty(COM)'
    if 'condenser' in r and ('ΔT' in r or '不成立' in r):
        return 'dist2_coldtop'
    if 'reboiler' in r and '不成立' in r:
        return 'reboiler_dT'
    if 'dP_excess' in r or 'ΔP/P_in' in r:
        return 'reactor_dP'
    if 'sv_out_of_range' in r or 'SV=' in r:
        return 'reactor_SV'
    if 'タイムアウト' in r or failure_unit == 'timeout':
        return 'timeout'
    if str(failure_unit).startswith('exception'):
        return 'exception'
    if str(failure_unit).startswith('spec_'):
        return failure_unit            # spec_production_under / spec_c3h6_purity 等
    if 'Wang-Henke' in r or 'rigorous' in r:
        return 'rigorous_fail'
    return failure_unit                # r1 / r2 / r3 / r_rx / r_psa / r_mem の素失敗


def _append_live_jsonl(trial, category: str, is_feas: bool, dur: float) -> None:
    """trial 1 件を JSONL に増分追記 (run 中の mid-run 解析用、~1 行/trial で軽量)。"""
    path = _LIVE_JSONL.get('path')
    if not path:
        return
    p = trial.params
    rec = {
        'n': trial.number,
        'v': (round(trial.value, 2) if trial.value is not None else None),
        'feas': is_feas, 'cat': category,
        'reason': (trial.user_attrs.get('failure_reason', '') or '')[:100],
        'sec': round(dur, 1),
        # 失敗診断に効く主要パラメータのみ (全変数は出さない=軽量)
        'col2_p': p.get('col2_p_kpa'), 'col2_n': p.get('col2_n_stages'),
        'col2_R': p.get('col2_reflux_ratio'), 'T_in': p.get('T_in_K'),
        'bed_dr': p.get('bed_thickness_m'), 'F': p.get('F_C3H8_fresh_kmol_h'),
        'A_mem': p.get('A_mem_m2'), 'P_H': p.get('P_H_Pa'),
    }
    try:
        import json as _json
        with open(path, 'a', encoding='utf-8') as f:
            f.write(_json.dumps(rec, ensure_ascii=False, default=str) + '\n')
    except Exception:
        pass   # ログ失敗は本計算を止めない


def _make_main_callback(n_total: int):
    state = {'start': None, 'prev_best': float('inf'), 'n_feas': 0,
             'n_done': 0, 'recent': _deque(maxlen=20), 'tally': _Counter()}

    def cb(study, trial):
        if state['start'] is None:
            state['start'] = _time.monotonic()
        dur = 0.0
        if trial.datetime_start is not None and trial.datetime_complete is not None:
            dur = (trial.datetime_complete - trial.datetime_start).total_seconds()
        state['recent'].append(dur); state['n_done'] += 1
        a = trial.user_attrs
        is_feas = bool(a.get('is_feasible', False))
        if is_feas:
            state['n_feas'] += 1
        # failure_unit より細かい「カテゴリ」で集計 (HYSYS 空出力 COM エラー・cold-top・
        # 反応器 ΔP/SV 等を live tally で区別)。さらに JSONL に増分追記して mid-run 解析可に。
        fu = a.get('failure_unit', '') or ('legacy' if not is_feas else '')
        category = _failure_category(fu, a.get('failure_reason', ''))
        if not is_feas:
            state['tally'][category] += 1
        _append_live_jsonl(trial, category, is_feas, dur)
        val = trial.value if trial.value is not None else float('inf')
        new_best = is_feas and val < state['prev_best']
        delta = state['prev_best'] - val if new_best else None
        if new_best:
            state['prev_best'] = val
        badge = '★ BEST  ' if new_best else ('✓ feas  ' if is_feas else '✗ infeas')
        val_s = '   ----' if val >= 9999.0 else f"{val:9.2f}"
        delta_s = f" (Δ-{delta:.1f})" if (delta is not None and delta < 1e9) else ""
        reason_s = ("  reason=" + _fmt_reason_from_trial(trial)) if not is_feas else ""
        line0 = f"[#{trial.number:03d}] {badge}  TAC={val_s}{delta_s}{reason_s}  {dur:5.1f}s"

        p = trial.params
        if 'z_cat_m' in p:        # 軸流 (旧)
            rx = (f"Rx(axial): T={p.get('T_in_K',0):.0f}K z={p.get('z_cat_m',0):.1f} "
                  f"t={p.get('t_cyc_min',0):.1f} D={p.get('D_reactor_m',0):.2f}")
        elif 'N_online' in p:     # Catofin 浅床・多基並列 (既定)
            rx = (f"Rx(catofin): T={p.get('T_in_K',0):.0f}K t={p.get('t_cyc_min',0):.1f} "
                  f"D={p.get('D_reactor_m',0):.2f} L_bed={p.get('L_bed_m',0):.2f} "
                  f"N={p.get('N_online',0)} dp={p.get('d_p_mm',0):.1f}mm")
        else:                     # 径方向流 (棄却、比較用)
            rx = (f"Rx(radial): T={p.get('T_in_K',0):.0f}K t={p.get('t_cyc_min',0):.1f} "
                  f"Di={p.get('D_inner_m',0):.1f} dr={p.get('bed_thickness_m',0):.2f} H={p.get('H_m',0):.1f}")
        v0 = (f"{rx} | PSA: D={p.get('D_psa_col_m',0):.2f} L={p.get('L_psa_bed_m',0):.1f} "
              f"des={p.get('desorption_target',0):.3f} | Mem: P_H={p.get('P_H_Pa',0)/1e5:.2f}bar "
              f"A={p.get('A_mem_m2',0):.2e} | F={p.get('F_C3H8_fresh_kmol_h',0):.0f}")
        v1 = (f"Dist1(SM): P={p.get('col1_p_kpa',0):.0f}kPa N={p.get('col1_n_stages',0)} "
              f"feed={p.get('col1_feed_stage',0)} cf={p.get('col1_comp_frac_2',0):.3f}")
        v2 = (f"Dist2(HY): P={p.get('col2_p_kpa',0):.0f}kPa N={p.get('col2_n_stages',0)} "
              f"fr={p.get('col2_feed_ratio',0):.2f} R={p.get('col2_reflux_ratio',0):.1f}")
        v3 = (f"Dist3(SM): P={p.get('col3_p_kpa',0):.0f}kPa N={p.get('col3_n_stages',0)} "
              f"fr={p.get('col3_feed_ratio',0):.2f}")
        pur = a.get('c3h6_purity_wtfrac'); prod = a.get('production_kmol_h')
        fused = a.get('F_C3H8_fresh_used_kmol_h')
        outs = ""
        if pur and prod:
            outs = f"       -> purity={float(pur)*100:.2f}wt% prod={float(prod):.0f}kmol/h"
            if fused:
                outs += f" yield={float(prod)/float(fused)*100:.1f}%"

        elapsed = _time.monotonic() - state['start']; n = state['n_done']
        med = (sorted(state['recent'])[len(state['recent']) // 2] if state['recent'] else 0.0)
        eta = max(n_total - n, 0) * med
        pct = 100.0 * n / max(n_total, 1)
        filled = int(_BAR_W * n / max(n_total, 1))
        bar = '█' * filled + '░' * (_BAR_W - filled)
        feas_pct = 100.0 * state['n_feas'] / max(n, 1)
        best_s = f"{state['prev_best']:.2f}" if state['prev_best'] < 1e9 else '----'
        prog = (f"       [{bar}] {n}/{n_total} ({pct:.0f}%)  feas {state['n_feas']}/{n} "
                f"({feas_pct:.0f}%)  elapsed {_fmt_dur(elapsed)} ETA {_fmt_dur(eta)} "
                f"pace {med:.1f}s best {best_s}")
        tally_s = _fmt_tally(state['tally'], 5)

        print(line0, flush=True)
        print("       " + v0, flush=True)
        print("       " + v1, flush=True)
        print("       " + v2, flush=True)
        print("       " + v3, flush=True)
        if outs:
            print(outs, flush=True)
        print(prog, flush=True)
        if tally_s:
            print(f"       top fails: {tally_s}", flush=True)
        print('', flush=True)

    return cb


# ===========================================================================
# レポート / 保存
# ===========================================================================

def _save_trials_csv(study: optuna.Study, path: str) -> None:
    param_keys: list = []
    attr_keys: list = []
    for t in study.trials:
        for k in t.params:
            if k not in param_keys:
                param_keys.append(k)
        for k in t.user_attrs:
            if k not in attr_keys:
                attr_keys.append(k)
    header = ['number', 'value', 'state'] + param_keys + [f'attr.{k}' for k in attr_keys]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        for t in study.trials:
            row = [t.number, (t.value if t.value is not None else ''), t.state.name]
            row += [t.params.get(k, '') for k in param_keys]
            row += [t.user_attrs.get(k, '') for k in attr_keys]
            w.writerow(row)


def _save_best_reports(study: optuna.Study, out_dir: str, top_n: int) -> list:
    """上位 top_n 候補を再評価して exp3 形式の詳細レポート (CAPEX/OPEX/spec/HI 内訳) を保存。

    feasible 優先、無ければ TAC 最小。レポートは out_dir 直下に top{rank}_trial{N}.txt として書き出す。
    """
    comp = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
    # display_full_results は economics 前提なので feasible のみ対象。
    # feasible が無ければスキップ (penalty result でのクラッシュ回避)。
    if not feas:
        print("  feasible trial が無いため詳細レポートはスキップ (CSV/JSON は保存済み)", flush=True)
        return []
    cand = sorted(feas, key=lambda t: t.value)[:top_n]
    eval_kwargs = dict(apply_hi=APPLY_HI, apply_stage2=APPLY_STAGE2, hi_dT_min_K=HI_DT_MIN_K)
    saved = []
    for rank, t in enumerate(cand, 1):
        try:
            design = _build_design(t.params)
            F_fresh = float(t.params.get('F_C3H8_fresh_kmol_h'))
            res = evaluate(design, _CONFIG, verbose=False,
                           apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
                           apply_stage2=APPLY_STAGE2, F_C3H8_override=F_fresh)
            # レポート本文を一度 StringIO に組み立て、ファイル保存とコンソール出力で共用する。
            # top1 はコンソールにも全文を出す。
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                print(f"# main.py top{rank}  trial #{t.number}  "
                      f"effective_TAC(BO)={t.value:.2f} 億円/年  "
                      f"feasible={t.user_attrs.get('is_feasible')}")
                print("#" + "=" * 70)
                show_input_snapshot(design, _CONFIG, eval_kwargs)
                display_full_results(res, design, _CONFIG)
            report_text = buf.getvalue()
            path = os.path.join(out_dir, f"top{rank}_trial{t.number}.txt")
            with open(path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            saved.append(path)
            if rank == 1:
                # top1 のみコンソールにも全文出力。
                print("\n" + "=" * 72, flush=True)
                print(f"  ★ ベスト候補 詳細レポート (top1 / trial #{t.number}) ─ コンソール表示",
                      flush=True)
                print("=" * 72, flush=True)
                print(report_text, flush=True)
            print(f"  top{rank} 詳細レポート(CAPEX/OPEX/spec内訳): {path}", flush=True)
        except Exception as e:
            print(f"  top{rank} レポート生成失敗 (trial #{t.number}): {type(e).__name__}: {e}", flush=True)
    return saved


def _write_readme(out_dir: str, ts: str, study: optuna.Study, best, saved_reports: list) -> None:
    """run subdir に README.md を出力 (結果の見方ガイド)。"""
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    top1_name = os.path.basename(saved_reports[0]) if saved_reports else None

    L = []
    L.append(f"# main.py (HYSYS+SM 全{len(SEARCH_SPACE)}変数 BO, 反応器={REACTOR_KIND}) run — {ts}")
    L.append("")
    L.append("Dist1/Dist3 = SM (学習済み GPR), Dist2 = HYSYS / 反応器・PSA・膜・F_fresh も変数化。")
    L.append("")
    L.append("## まず見るべきファイル (推奨順)")
    L.append("")
    if top1_name:
        L.append(f"1. **`{top1_name}`** ─ ★最終結果。ベスト候補の再評価詳細 "
                 "(CAPEX/OPEX/spec/HI 内訳 + 入力スナップショット)。")
        L.append(f"   - `top2_*` / `top3_*` は次点候補。同じ形式で比較できる。")
    else:
        L.append("1. **`top*_trial*.txt`** ─ ★最終結果 (今回は feasible 無しで未生成)。")
    L.append("2. **`best.json`** ─ BO 単体ベスト trial の params + 診断値。再現・簡易確認用。")
    L.append("3. **`trials.csv`** ─ 全 trial 履歴。Excel/pandas で散布図・統計解析。")
    if (USE_SQLITE_STORAGE or N_WORKERS > 1):
        L.append("4. **`optuna.db`** ─ Optuna SQLite。可視化: "
                 "`optuna-dashboard sqlite:///optuna.db`")
    L.append("")
    L.append("## この run の設定")
    L.append("")
    L.append(f"- N_TRIALS = {N_TRIALS}, N_STARTUP(QMC) = {N_STARTUP}, N_TOPK = {N_TOPK}")
    L.append(f"- SAMPLER = {SAMPLER}, SEED = {SEED}, N_WORKERS = {N_WORKERS}")
    L.append(f"- 探索変数数 = {len(SEARCH_SPACE)} "
             f"(反応器[{REACTOR_KIND}]{len(_REACTOR_SPACE)} + PSA3 + 膜2 + 原料1 + Dist1/2/3 各4/4/3)")
    L.append(f"- Dist1 = SM, Dist2 = HYSYS, Dist3 = SM (spec なし)")
    L.append(f"- purity 閾値 = 99.45 wt% (SM Dist3 の 99.5 mol%=99.497 wt% を尊重した緩和)")
    L.append("")
    L.append("## ベスト要約")
    L.append("")
    L.append(f"- 完了 trial = {len(complete)} / feasible = {len(feasible)}")
    if best is not None:
        tag = "feasible ✓" if best.user_attrs.get('is_feasible', False) else "infeasible ✗ (feasible 無し、TAC 最小)"
        L.append(f"- ベスト: **trial #{best.number}** ({tag})")
        L.append(f"- effective_TAC = **{best.value:.3f}** 億円/年")
        try:
            _pur  = float(best.user_attrs.get('c3h6_purity_wtfrac'))
            _prod = float(best.user_attrs.get('production_kmol_h'))
            _ff   = float(best.params.get('F_C3H8_fresh_kmol_h'))
            L.append(f"- purity = {_pur*100:.2f} wt%, 生産量 = {_prod:.1f} kmol/h, "
                     f"F_fresh = {_ff:.1f} kmol/h, 収率 = {_prod/_ff*100:.1f}%")
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    else:
        L.append("- 完了 trial なし")
    L.append("")
    with open(os.path.join(out_dir, 'README.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


def main():
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    print(f"==== main.py: PDH HYSYS+SM 全{len(SEARCH_SPACE)}変数 制約付き BO (反応器={REACTOR_KIND}) ====", flush=True)
    print(f"  N_TRIALS={N_TRIALS}, N_STARTUP(QMC)={N_STARTUP}, sampler={SAMPLER}, seed={SEED}, top-k report={N_TOPK}", flush=True)
    print(f"  Dist1/Dist3 = SM, Dist2 = HYSYS / 上流(反応器・PSA・膜)・F_fresh も変数化", flush=True)
    # run ごとの subdir へ成果物を集約。
    # outputs/main_<ts>/{trials.csv, best.json, top*.txt, optuna.db, README.md}
    out_dir = os.path.join(OUTPUT_DIR, f'main_{ts}')
    os.makedirs(out_dir, exist_ok=True)
    # 増分 JSONL ログのパスをセット。trial ごとに 1 行追記され、
    # run 中に `tools/_scan_run_log.py outputs/main_<ts>/trials_live.jsonl` で失敗内訳を即解析できる。
    _LIVE_JSONL['path'] = os.path.join(out_dir, 'trials_live.jsonl')

    sampler = make_sampler(SAMPLER, SEED, N_STARTUP, constraints_func=constraints_func)
    # 並列(N_WORKERS>1)時は worker 間で study 共有のため SQLite 必須。
    _use_sqlite = USE_SQLITE_STORAGE or N_WORKERS > 1
    _db_path = os.path.join(out_dir, 'optuna.db')
    storage = f"sqlite:///{_db_path}" if _use_sqlite else None
    study = optuna.create_study(
        study_name=f"{STUDY_NAME}_{ts}" if _use_sqlite else STUDY_NAME,
        sampler=sampler, direction='minimize',
        storage=storage, load_if_exists=bool(storage),
    )

    if N_WORKERS > 1 and storage is not None:
        # N worker プロセスで共有 study を分担。各 worker は自前の HYSYS を起動 (重い)。
        # 完了後は study.trials(DB) を読んでサマリ/レポートを生成する。
        print(f"  並列実行: {N_WORKERS} worker (各々 HYSYS 起動)。worker ログ: outputs/_worker*.log", flush=True)
        from optimization.parallel import spawn_workers
        spawn_workers(
            kind='main', study_name=study.study_name, storage_url=storage,
            db_path=_db_path, n_workers=N_WORKERS, n_trials_total=N_TRIALS,
            n_startup=N_STARTUP, base_seed=SEED, out_dir=out_dir,
        )
    else:
        run_optimization(
            study, objective, n_trials=N_TRIALS,
            show_progress_bar=False, n_jobs=N_JOBS,
            callbacks=[_make_main_callback(N_TRIALS)],
        )

    # ---- 結果サマリ ----
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    print("\n==== 結果 ====", flush=True)
    print(f"  完了 trial: {len(complete)} / feasible: {len(feasible)}", flush=True)
    best = None
    if feasible:
        best = min(feasible, key=lambda t: t.value); tag = "feasible best"
    elif complete:
        best = min(complete, key=lambda t: t.value); tag = "best (feasible 無し)"
    if best is not None:
        print(f"  {tag}: trial #{best.number}  effective_TAC={best.value:.2f} 億円/年", flush=True)
        # 見やすい入力スナップショットは top1 詳細レポート (show_input_snapshot) に出力され、
        # 再現用の生 params は best JSON に保存される。ここでは要点 1 行のみ。
        try:
            _pur  = float(best.user_attrs.get('c3h6_purity_wtfrac'))
            _prod = float(best.user_attrs.get('production_kmol_h'))
            _ff   = float(best.params.get('F_C3H8_fresh_kmol_h'))
            print(f"    purity={_pur*100:.2f}wt%  prod={_prod:.1f}kmol/h  "
                  f"F_fresh={_ff:.1f}kmol/h", flush=True)
        except Exception:
            pass

    # ---- 保存 (run subdir に集約) ----
    trials_csv = os.path.join(out_dir, 'trials.csv')
    _save_trials_csv(study, trials_csv)
    print(f"\n  trial 履歴 CSV: {trials_csv}", flush=True)
    saved_reports: list = []
    if best is not None:
        with open(os.path.join(out_dir, 'best.json'), 'w', encoding='utf-8') as f:
            json.dump({'number': best.number, 'effective_TAC': best.value,
                       'params': best.params,
                       'user_attrs': {k: v for k, v in best.user_attrs.items()}},
                      f, ensure_ascii=False, indent=2, default=str)
        print(f"  best JSON: {os.path.join(out_dir, 'best.json')}", flush=True)
        # 上位候補の詳細レポート (CAPEX/OPEX/spec/HI 内訳)
        print(f"\n  上位 {N_TOPK} 候補の詳細レポートを生成中...", flush=True)
        saved_reports = _save_best_reports(study, out_dir, N_TOPK)

    # ---- README (結果の見方ガイド) ----
    _write_readme(out_dir, ts, study, best, saved_reports)

    # ---- 成果物サマリ ----
    print()
    print("=" * 72, flush=True)
    print(f"成果物: {os.path.abspath(out_dir)}/", flush=True)
    print("=" * 72, flush=True)
    print(f"  📌 README.md       … 結果の見方ガイド (最初に開いて)", flush=True)
    if saved_reports:
        print(f"  ★ top1_*.txt      … ベスト候補の詳細 (★最終結果はここ)", flush=True)
    print(f"  ・ best.json       … BO ベスト trial (JSON、簡易)", flush=True)
    print(f"  ・ trials.csv      … 全 {N_TRIALS} trial 履歴 (Excel/pandas で解析)", flush=True)
    if (USE_SQLITE_STORAGE or N_WORKERS > 1):
        print(f"  ・ optuna.db       … SQLite (中断・再開・dashboard 用)", flush=True)
    print("=" * 72, flush=True)

    try:
        from units.vle.hysys.provider import shutdown_default_provider
        shutdown_default_provider()
    except Exception:
        pass


if __name__ == "__main__":
    main()
