r"""
exp3.py — リサイクルあり PDH プロセス全体フロー シミュレーション (HYSYS + SM バックエンド)

exp1.py の HYSYS 版。蒸留塔は塔ごとにバックエンドを選択:
  - Dist1 / Dist3 : SM (学習済み GPR サロゲート、2026-05-25 導入)。HYSYS 解とほぼ一致を検証済み。
  - Dist2         : HYSYS COM (SM 化できなかったため HYSYS 継続)。
それ以外 (反応器・PSA・膜・リサイクル合流) は exp1 と同じく既存実装で動かす。

高速化の経緯 (2026-05-25):
  290s (全塔 HYSYS) → 159s (Dist1 メモ化 + swap_case sleep 短縮) → 26s (Dist1/Dist3 SM 化)。
  SM 化で HYSYS 塔が Dist2 のみになり HSC swap が消滅。Dist2 は開きっぱなしで warm 再解に
  なるが、force-cold (PDH_HYSYS_FORCE_COLD=1) 版と結果完全一致で経路依存なしを実証済み。

設計判断 (2026-05-22):
  HYSYS の探索変数は exp1 (FUG/rigorous) と異なる:
    - exp1   : ColumnTunables.reflux_ratio / recovery_LK_top / recovery_HK_bot
    - exp3   : ColumnTunables.hysys_spec_value (column1=Comp Fraction-2,
                                                column2=Reflux Ratio,
                                                column3=Draw Rate kgmol/s)
              + hysys_feed_stage
  そのため main.py / optimization/study.py には組み込まず、独立スクリプトとして運用。
  HI / Stage 2 は exp1 と同じく True で実行 (HYSYS 出力からも Q_cond/Q_reb が取れる)。

使い方:
  下の "実験で振る設計変数" を編集して .\.venv\Scripts\python.exe exp/exp3.py
  N_stages は HSC 存在範囲内のみ:
    column1: 30〜60
    column2: 15〜80
    column3: 69〜200
"""

import os
import sys

# 設計判断 (2026-05-22): HYSYS は 1 塔あたり ~15秒、3塔 × リサイクル iter で
# 既定の 60s 予算では 2 iter で打ち切られる。HYSYS バックエンド向けに 30分予算に拡大。
os.environ.setdefault('PDH_TRIAL_TIME_BUDGET_SEC', '1800')
os.environ.setdefault('PDH_PER_UNIT_LOG', '1')   # 各 iter で penalty 発火ユニットを stderr へ

# Windows コンソール (cp932) でも記号を表示
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from config.load import load_operating_config
from flowsheet import FlowsheetDesignVars, evaluate
from src.distillation_core import ColumnTunables
from units.reactors.swing import DesignVars as SwingDesign
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars
from simulation import display_full_results, hdr, show_input_snapshot, run_exp


# ===========================================================================
#  実験で振る設計変数 (ここを書き換えて再実行)
# ===========================================================================

# 設計判断 (2026-05-26): special.py の BO best = trial #294 の最適化済み 21 変数を
# そのまま投入。出典 outputs/special_20260526_172608_best.json
#   (BO 記録 effective_TAC=1305.59 億円/年、feasible=True)。
# ユーザーが exp3 でフル詳細出力を見て手動検証するため。旧テンプレ値 (exp1 best #115
# 系) は git 履歴で復元可能。
# 注意: special では col2/col3 の feed を「比率」で suggest し _feed_stage_from_ratio で
#   絶対段に変換している。exp3 は絶対段を直接指定するため、変換後の値を記載:
#     Dist2: feed_ratio 0.553435 × N67 → FEED_STAGE_dist2 = 37  (レポート kirkbride 37 と一致)
#     Dist3: feed_ratio 0.715715 × N156 → FEED_STAGE_dist3 = 112 (同 112 と一致)

# === 反応器 (Swing) ============================================================
T_in_K        = 936.9962277
z_cat_m       = 29.1345707
t_cyc_min     = 19.0206554
D_reactor_m   = 9.4751160

# === PSA =====================================================================
D_psa_col_m       = 4.4154197
L_psa_bed_m       = 27.8226005
desorption_target = 0.2974322

# === 膜分離 ===================================================================
P_H_Pa     = 8.0711387e5
P_L_Pa     = 1.0e5
A_mem_m2   = 1.1854758e5

# === Dist1 (脱ブタン塔) SM ===================================================
P_dist1_kPa            = 1996.1098
N_dist1                = 31
FEED_STAGE_dist1       = 28
COMP_FRAC_2_dist1      = 0.9083715

# === Dist2 (脱エタン塔) HYSYS ===============================================
# P_dist2=645.4 kPa < P_H=807.1 kPa (Mem の ph_le_pfeed 回避を満たす)。
P_dist2_kPa            = 645.4137
N_dist2                = 67
FEED_STAGE_dist2       = 37       # special: feed_ratio 0.553435 × N67 → 37
REFLUX_RATIO_dist2     = 8.9916122

# === Dist3 (C3 スプリッタ) SM ===============================================
P_dist3_kPa            = 1677.1573
N_dist3                = 156
FEED_STAGE_dist3       = 112      # special: feed_ratio 0.715715 × N156 → 112
DRAW_RATE_dist3_kmolh  = 0.99     # SM では未使用 (Dist3 は spec なし)

# === Fresh LPG ==============================================================
F_C3H8_fresh_kmol_h = 1523.3995


# ===========================================================================
#  HI オプション (exp1 と同じ)
# ===========================================================================
APPLY_HI     = True
APPLY_STAGE2 = True
HI_DT_MIN_K  = 10.0


# ===========================================================================
#  出力モード
# ===========================================================================
SAVE_OUTPUT = True


# ===========================================================================
#  フローシート設計変数の組み立て
# ===========================================================================
#  ColumnTunables の P_col は Pa 単位、kPa から変換する。
#  hysys 経路では reflux_ratio (PDH 側) と recovery_LK_top/recovery_HK_bot は
#  参照されないので dummy 値で OK。HYSYS の主スペックは hysys_spec_value で渡す。
#  draw_rate は HYSYS では kgmol/s 単位で書込むため kmolh → kgmol/s に変換。
design = FlowsheetDesignVars(
    swing=SwingDesign(
        T_in=T_in_K, z_cat=z_cat_m, t_cyc=t_cyc_min, D=D_reactor_m,
    ),
    psa=PSADesignVars(
        D_col=D_psa_col_m, L_bed=L_psa_bed_m,
        desorption_target=desorption_target,
    ),
    mem=MemDesignVars(
        P_H=P_H_Pa, P_L=P_L_Pa, A_mem=A_mem_m2, P_dist=P_dist3_kPa * 1000.0,
    ),
    dist1=ColumnTunables(
        P_col=P_dist1_kPa * 1000.0,
        N_stages=N_dist1, N_feed=1,
        reflux_ratio=2.0,                       # dummy (SM/HYSYS 経路では未使用)
        # 設計判断 (2026-05-25): Dist1 を SM (学習済み GPR) に置換。
        # In_CompFraction2 = hysys_spec_value を流用。HYSYS 解とほぼ完全一致を検証済み。
        solver_method='sm',
        hysys_spec_value=COMP_FRAC_2_dist1,
        hysys_feed_stage=FEED_STAGE_dist1,
    ),
    dist2=ColumnTunables(
        P_col=P_dist2_kPa * 1000.0,
        N_stages=N_dist2, N_feed=1,
        reflux_ratio=REFLUX_RATIO_dist2,        # exp1 系の dummy 互換、未使用
        solver_method='hysys',
        hysys_spec_value=REFLUX_RATIO_dist2,
        hysys_feed_stage=FEED_STAGE_dist2,
    ),
    dist3=ColumnTunables(
        P_col=P_dist3_kPa * 1000.0,
        N_stages=N_dist3, N_feed=1,
        reflux_ratio=12.0,                      # dummy
        # 設計判断 (2026-05-25): Dist3 を SM (学習済み GPR) に置換 (速度の本命)。
        # model3 は spec 入力なし → 分配(回収率)は SM 予測をそのまま設計値として採用
        # (ユーザー決定)。製品純度 99.5% は満たす。hysys_spec_value は SM では未使用。
        solver_method='sm',
        # 設計判断 (2026-05-22): adapter で 0-1 範囲なら recovery_LK_top (動的 Draw Rate
        # 計算)、1.0 超なら絶対量 (kgmol/s) として扱う。
        # DRAW_RATE_dist3_kmolh が 0.99 (純度モード) なら /3600 せずそのまま渡す。
        hysys_spec_value=(DRAW_RATE_dist3_kmolh
                          if 0 < DRAW_RATE_dist3_kmolh <= 1.0
                          else DRAW_RATE_dist3_kmolh / 3600.0),
        hysys_feed_stage=FEED_STAGE_dist3,
    ),
)


# ===========================================================================
#  実行 + 結果表示 (exp1 と共通の run_exp ラッパ)
# ===========================================================================
#  設計判断 (2026-05-25): Wegstein 減衰 (q_min=-5→-2) も検討したが、振動抑制で
#  反復は 10→9 に減るものの収束経路が変わり 1% 許容球内の別点に着地して結果が
#  わずかに動く (生産量 1248→1258) 割に効果が限定的だった。結果中立な高速化
#  (Dist1 メモ化 + swap_case 待ち時間削減) を優先し、減衰は不採用とした。
config = load_operating_config()


def _run_simulation():
    eval_kwargs = {
        'apply_hi':     APPLY_HI,
        'apply_stage2': APPLY_STAGE2,
        'hi_dT_min_K':  HI_DT_MIN_K,
    }
    show_input_snapshot(design, config, eval_kwargs)
    hdr("外側ループ: 製品流量厳密化 (Fresh を調整) — HYSYS バックエンド")
    res = evaluate(
        design, config, verbose=True,
        apply_hi=APPLY_HI, hi_dT_min_K=HI_DT_MIN_K,
        apply_stage2=APPLY_STAGE2,
        F_C3H8_override=F_C3H8_fresh_kmol_h,
    )
    display_full_results(res, design, config)
    return res


_out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'outputs'))
result = run_exp(
    label                 = 'exp3',
    eval_callable         = _run_simulation,
    output_dir            = _out_dir,
    save_output           = SAVE_OUTPUT,
    expected_outer_iters  = 6,
)

# シミュレーション終了時に HYSYS default provider を確実にクローズ
# (プロセス残留対策、Python 終了時の __del__ より早く)
try:
    from units.vle.hysys.provider import shutdown_default_provider
    shutdown_default_provider()
except Exception:
    pass

if result.economics is None:
    sys.exit(1)
