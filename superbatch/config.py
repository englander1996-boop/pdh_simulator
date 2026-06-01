# -*- coding: utf-8 -*-
"""superbatch.config — 検証計画とパス・定数。**編集はこのファイルだけ**。"""
import os
import glob

# --- 固定パス (repo 直下基準) ---
# 設計方針 (2026-06-01 ユーザー指示): main.py と同様に、super_main の 1 回の実行 =
# **タイムスタンプ付き親フォルダ 1 個 (outputs/super_main_<ts>/) に全部内包**する。
# その中に runs/(各 main 実行の成果物) logs/ plots/ manifest.jsonl summary.* を置く。
# outputs 直下には super_main_<ts>/ が 1 個できるだけで散乱しない。
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(REPO, '.venv', 'Scripts', 'python.exe')   # main.py を回す python
OUTPUTS = os.path.join(REPO, 'outputs')
SUPER_PREFIX = 'super_main_'    # 親フォルダ名の接頭辞 (outputs/super_main_<ts>/)

# --- バッチ dir 相対パス (実行時に set_batch_dir() で確定。未設定なら None) ---
BATCH_DIR = None
LOG_DIR = None
PLOTS_DIR = None
RUNS_DIR = None
MANIFEST = None
SUMMARY_TXT = None
SUMMARY_CSV = None


def set_batch_dir(path):
    """このバッチの親フォルダを確定し、配下パスを全モジュールへ反映する。"""
    global BATCH_DIR, LOG_DIR, PLOTS_DIR, RUNS_DIR, MANIFEST, SUMMARY_TXT, SUMMARY_CSV
    BATCH_DIR = path
    LOG_DIR = os.path.join(path, 'logs')             # 各 run の BO ライブログ
    PLOTS_DIR = os.path.join(path, 'plots')          # 可視化 PNG
    RUNS_DIR = os.path.join(path, 'runs')            # 各 main 実行の成果物 (移動先)
    MANIFEST = os.path.join(path, 'manifest.jsonl')  # 完了 run 記録 (resume)
    SUMMARY_TXT = os.path.join(path, 'summary.txt')  # 集計テキスト
    SUMMARY_CSV = os.path.join(path, 'summary.csv')  # 1 run 1 行の一覧


def _manifest_done_count(batch_dir):
    mpath = os.path.join(batch_dir, 'manifest.jsonl')
    if not os.path.exists(mpath):
        return 0
    return sum(1 for ln in open(mpath, encoding='utf-8') if ln.strip())


def find_or_create_batch_dir(ts):
    """resume 対応のバッチ dir 選択。

    最新の super_main_<ts> が未完了 (manifest の done < len(PLAN)) ならそれを resume、
    無ければ outputs/super_main_<ts>/ を新規作成する。戻り値 (path, resumed)。
    """
    existing = sorted(glob.glob(os.path.join(OUTPUTS, SUPER_PREFIX + '*')))
    for d in reversed(existing):
        if os.path.isdir(d) and _manifest_done_count(d) < len(PLAN):
            return d, True
    return os.path.join(OUTPUTS, SUPER_PREFIX + ts), False


# ===========================================================================
# 検証計画 — ここを編集する
# ===========================================================================
# 1 run あたりの BO 試行数 (PDH_N_TRIALS)。
# 設計判断 (2026-06-01): 300→400。run2(新バウンド)分析で best が末尾 #299 に出る等
# 300 では広いバウンドに対し未収束/ノイズ blip が起きた。400 で収束余裕を持たせる。
N_TRIALS_PER_RUN = 400

# 各 run = (sampler, seed)。seed は QMC/TPE/CMA-ES/Random すべてに渡る (make_sampler)
# ので run ごとに探索経路が変わる。
#   tpe    : 本命 (best-TAC 分布 / ロバスト性)
#   cmaes  : 別手法クロスチェック (同じ最適に行くか = 大域性の傍証)
#   random : 対照群 (BO がランダム探索より良いか = BO 正当性)
#
# 設計判断 (2026-06-01): 30 run (TPE 20 + cmaes 5 + random 5)。1 run ~2.4h(400 trial)
# → 約 72h (~3 日)。分散推定+対照群はこの規模で十分。
TPE_SEEDS    = list(range(1001, 1021))   # 20
CMAES_SEEDS  = list(range(3001, 3006))   # 5
RANDOM_SEEDS = list(range(2001, 2006))   # 5
# ===========================================================================

PLAN = ([('tpe', s) for s in TPE_SEEDS]
        + [('cmaes', s) for s in CMAES_SEEDS]
        + [('random', s) for s in RANDOM_SEEDS])

SAMPLERS = ('tpe', 'cmaes', 'random')
COLORS = {'tpe': 'tab:blue', 'cmaes': 'tab:green', 'random': 'tab:gray'}

# 1 run の想定所要 (ETA 計算の初期値、実測が貯まれば中央値で上書きされる)
# 400 trial: 前回 300 trial 実測 1.74h を線形外挿 → ~2.3h。
EST_RUN_SEC = 2.4 * 3600

# params 収束/strip プロットで見る主要変数 (反応器・原料・Dist2 系)
KEY_VARS = ['F_C3H8_fresh_kmol_h', 'D_inner_m', 'bed_thickness_m', 'H_m',
            'T_in_K', 'col2_p_kpa', 'col2_reflux_ratio']
