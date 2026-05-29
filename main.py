r"""
main.py — HYSYS + SM バックエンドでの PDH プロセス全変数最適化 (制約付き Optuna BO)

exp3.py を評価関数として、**全 21 設計変数**(反応器・PSA・膜・原料 + 蒸留塔 3 つ)を
Bayesian Optimization で探索する。**BO を実行する本丸ファイル**。
  (旧 main.py = FUG/rigorous 全フローシート版は sub/sub1.py、
   旧 final.py = SM/rigorous/SM 版は sub/sub2.py に退避済み。
   本 main は旧 special.py を改名したもの = HYSYS+SM 版。)

設計判断 (2026-05-25): sub1 (旧 main.py) の BO 成功インフラを移植 + HYSYS/SM 特性に適応。
  借用 (sub1 = 旧 main.py): QMC→TPE 2相サンプラ・constraints_func(連続制約)・penalty_scale・
    _store_diagnostics・compact callback(flush 付きライブログ, callbacks.py)・
    詳細レポート保存(display_full_results)。
  適応 (HYSYS/SM): Dist1/Dist3=SM(学習済み GPR, ~瞬時)、Dist2=HYSYS。
    探索 bounds は sub1 (旧 main.py) の forensic 値 ∩ SM 分類器 feasible 領域。
    純度は SM Dist3 が 99.5 mol%=99.497 wt% 固定 → spec を 99.45 wt% に緩和(決定A)。
    塔本体 CAPEX は provider 側で FUG と同式で計算済み(N/還流が CAPEX に効く)。

21 変数:
  反応器(4): T_in_K, z_cat_m, t_cyc_min, D_reactor_m
  PSA(3)   : D_psa_col_m, L_psa_bed_m, desorption_target
  膜(2)    : P_H_Pa, A_mem_m2   (P_L=1atm 固定、mem.P_dist=Dist3圧 同期)
  原料(1)  : F_C3H8_fresh_kmol_h
  Dist1(4) : col1_p_kpa, col1_n_stages, col1_feed_stage, col1_comp_frac_2  (SM)
  Dist2(4) : col2_p_kpa, col2_n_stages, col2_feed_ratio, col2_reflux_ratio (HYSYS)
  Dist3(3) : col3_p_kpa, col3_n_stages, col3_feed_ratio                    (SM, spec なし)

出力 (sub1 = 旧 main.py 流に run ごとの subdir へ集約):
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
from units.separators.psa.psa_system import PSADesignVars
from units.separators.membrane.membrane_system import MemDesignVars

from optimization.study import make_sampler, run_optimization, _default_constraints_func
from optimization.objective import _store_diagnostics
from optimization.penalty_scale import set_scale, default_schedule
from simulation import display_full_results, show_input_snapshot


# ===========================================================================
# § 1. BO 設定
# ===========================================================================
N_TRIALS    = 300            # sub1(旧main.py) 準拠。全21変数なので 300 推奨 (~1-1.5h 目安)
N_STARTUP   = 50             # QMC Sobol 広域カバレッジ (以降 TPE)
SEED        = 42
N_JOBS      = 1              # HYSYS COM + penalty_scale global のため 1 (スレッド並列は不可)
# 設計判断 (2026-05-26): マルチプロセス並列 worker 数。>1 で N プロセスが共有 SQLite study を
# 分担し、各 worker が自前の HYSYS インスタンスを起動する (重い: 各 ~数百MB + 起動 ~1分。
# RAM 次第で 3 worker 程度が現実的)。各 worker 単スレッドで penalty_scale/GIL 問題なし、
# constant_liar=True で冗長サンプリング抑制。1 で従来の単一プロセス。
# 設計判断 (2026-05-26): 本 main (HYSYS+SM) は N_WORKERS=1 のみ正しい。HYSYS.Application は
# ユーザーセッションに単一インスタンスの COM サーバーで、複数 worker が Dispatch しても
# 同じ HYSYS を共有 → ケース切替(Close+Open)が跨プロセスで割り込み結果汚染する (実証済)。
# 並列化できるのは HYSYS 不使用の sub1(旧main)/sub2(旧final) のみ。本 main は SM 化で元々高速(~30分)。
N_WORKERS   = 1
N_TOPK      = 3              # 詳細レポートを出す上位候補数

USE_SQLITE_STORAGE = False
STUDY_NAME         = "pdh_hysys_sm_main"


# ===========================================================================
# § 2. 固定値 / 評価オプション
# ===========================================================================
P_L_Pa       = 1.0e5         # 膜透過側圧力 (大気圧固定、sub1 旧main と同じ)
APPLY_HI     = True
APPLY_STAGE2 = True
HI_DT_MIN_K  = 10.0


# ===========================================================================
# § 3. 探索範囲 — 全 21 変数  形式: (low, high, scale, type)
#   bounds は sub1(旧main) の forensic 調整値 ∩ SM 分類器 feasible 領域 (2026-05-25 プローブ)。
#   SM 崖 (予測不能域) は除外、解の縁 (収束ぎりぎり) は含める。
# ===========================================================================
SEARCH_SPACE = {
    # ----- 反応器 (Swing) — sub1(旧main.py) 準拠 -----
    # 設計判断 (2026-05-26): 940→925 に下げたら feas 0/61 に激減 → 940 へ戻す。
    # 理由: 本 main は Dist3 が固定スペック(99.5mol%/~1200製品)で、高い C3H6 throughput を
    # 維持するには高転化=高 T_in が必須。前回 feasible は全て T_in 929-940。925 で全除外された。
    # sub1(旧main) は自由Dist3で低Tでも1129生産できる(選択率↑で収率79%)が、本 main では低T=production不足。
    # → 本 main の収率は固定1200 Dist3 demand により高T(~76%)で構造的に頭打ち。T_in は下げない。
    "T_in_K":              (880.0,  940.0,  'linear', 'float'),
    "z_cat_m":             (15.0,   40.0,   'linear', 'float'),
    "t_cyc_min":           (12.0,   25.0,   'linear', 'float'),
    "D_reactor_m":         (7.0,    10.0,   'linear', 'float'),

    # ----- PSA — sub1(旧main.py) 準拠 -----
    "D_psa_col_m":         (2.9,    5.0,    'linear', 'float'),
    "L_psa_bed_m":         (22.0,   30.0,   'linear', 'float'),
    "desorption_target":   (0.22,   0.40,   'linear', 'float'),

    # ----- 膜 (P_L 固定、P_dist は Dist3 圧と同期) — sub1(旧main.py) 準拠 -----
    "P_H_Pa":              (7.5e5,  9.5e5,  'linear', 'float'),
    "A_mem_m2":            (5.0e4,  3.0e5,  'log',    'float'),

    # ----- 原料 (外側ループ skip、override) -----
    # 設計判断 (2026-05-26): (1380,1500) → (1500,1750) に上げる。
    # 根拠: HYSYS+SM フローの実効収率は ~72% (観測 66-77%) と FUG/sub1(旧main) より低い。
    #   旧上限 1500 では F×yield = 1500×0.72 ≈ 1083 で生産量下限 1128.6 (target 1188 の -5%)
    #   すら割り、special_run_20260526_125938 で Dist2 を通過した 28 trial の 27 件が
    #   prod_under、feasible は 1/154 に留まった。中央収率で target 1188 を満たすには
    #   F ≈ 1188/0.72 ≈ 1645 が必要。sub1(旧main) の縮小前 (1200,1700) 寄りに戻す方向で、
    #   feasible 生産量帯 [1128.6, 1247.4] を F×yield の範囲で到達可能にする。
    "F_C3H8_fresh_kmol_h": (1500.0, 1750.0, 'linear', 'float'),

    # ----- Dist1 (SM model1: N30-60/P1600-2200/feed_stage10-39/CF0.9-0.999) -----
    # feed_stage は SM feas ≥22 (プローブ: <21 で 0%)。範囲は (22,28) 固定: N 下限 30 でも
    # N-2=28 以下で常に有効になり、動的 search space を避けて TPE(multivariate) をフル活用。
    # (SM feas 域 ≤39 のうち 22-28 を採用。Dist1 は cost driver でないため探索損失は許容。)
    "col1_p_kpa":          (1600.0, 2000.0, 'linear', 'float'),
    "col1_n_stages":       (30,     60,     'linear', 'int'  ),
    "col1_feed_stage":     (22,     28,     'linear', 'int'  ),
    "col1_comp_frac_2":    (0.90,   0.999,  'linear', 'float'),

    # ----- Dist2 (HYSYS 脱エタン塔)。収束 envelope 狭、縁を含む (N≈45 安/N=75-80 頑健) -----
    # 設計判断 (2026-05-26): 上限 620→700。Dist2(HYSYS) 塔頂を浅冷化(高P→塔頂温度↑)して
    # 深冷コンデンサ費(エチレン-100C 17731円/GJ)を下げるレバー。P_H 下限 750kPa 未満を維持(膜 ph_le_pfeed 回避)。
    "col2_p_kpa":          (500.0,  700.0,  'linear', 'float'),  # 膜 P_H 未満 (上限 620→700, 浅冷化)
    "col2_n_stages":       (44,     80,     'linear', 'int'  ),
    "col2_feed_ratio":     (0.40,   0.60,   'linear', 'float'),
    # 設計判断 (2026-05-26): 上限 13→10.5。Dist2(HYSYS) 深冷コンデンサ(−83°C, エチレン-100C)
    # が 本 main TAC の ~25%(427億)の最大コスト。還流比↓ = 凝縮 duty↓ = 直接削減。
    # 高還流(13)を切って低還流帯へ誘導。下限 8 は維持(HYSYS 収束 envelope + C3 封じ込め)。
    "col2_reflux_ratio":   (8.0,    10.5,   'linear', 'float'),

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
# 決定A (2026-05-25): SM Dist3 は 99.5 mol%=99.497 wt% 固定。mol↔wt 差 0.003pp を吸収するため
# purity 閾値を 99.45 wt% に緩和 (SM の実力を尊重)。詳細は project_sm_integration メモ。
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
    """params dict (21 変数) から FlowsheetDesignVars を構築。trial 非依存 (best 再評価でも使用)。"""
    n1 = int(p['col1_n_stages']); fs1 = int(p['col1_feed_stage'])
    n2 = int(p['col2_n_stages']); fs2 = _feed_stage_from_ratio(p['col2_feed_ratio'], n2, *_FEED_STAGE_ABS['col2'])
    n3 = int(p['col3_n_stages']); fs3 = _feed_stage_from_ratio(p['col3_feed_ratio'], n3, *_FEED_STAGE_ABS['col3'])
    p3_kpa = float(p['col3_p_kpa'])
    return FlowsheetDesignVars(
        swing=SwingDesign(T_in=p['T_in_K'], z_cat=p['z_cat_m'],
                          t_cyc=p['t_cyc_min'], D=p['D_reactor_m']),
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
    """sub1(旧main) 由来の連続制約 (feas / 生産量方向 / 反応器SV / PSA / 膜 / 塔 shortfall)。

    全変数最適化なので反応器/PSA/膜の shortfall 信号も活性化し、TPE が上流の方向も学習する。
    純度は SM で不変のため制約化しない (定数=無意味)。
    """
    return _default_constraints_func(trial)


# ===========================================================================
# ライブログ用 compact callback (sub1=旧main の make_compact_callback 相当、flush 付き)
# ===========================================================================
from collections import Counter as _Counter, deque as _deque
import time as _time
from optimization.callbacks import _fmt_dur, _fmt_reason_from_trial, _fmt_tally

_BAR_W = 30


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
        fu = a.get('failure_unit', '') or ('legacy' if not is_feas else '')
        if fu:
            state['tally'][fu] += 1
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
        v0 = (f"Rx: T={p.get('T_in_K',0):.0f}K z={p.get('z_cat_m',0):.1f} t={p.get('t_cyc_min',0):.1f} "
              f"D={p.get('D_reactor_m',0):.2f} | PSA: D={p.get('D_psa_col_m',0):.2f} L={p.get('L_psa_bed_m',0):.1f} "
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

    sub1(旧main) の display_best_full / top-k レポート相当。feasible 優先、無ければ TAC 最小。
    レポートは out_dir 直下に top{rank}_trial{N}.txt として書き出す。
    """
    comp = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    feas = [t for t in comp if t.user_attrs.get('is_feasible', False)]
    # 設計判断: display_full_results は economics 前提なので feasible のみ対象。
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
            # 設計判断 (2026-05-26): レポート本文を一度 StringIO に組み立て、ファイル保存と
            # コンソール出力で共用する。top1 は sub1(旧main) (pipeline._display_best_full) と同様に
            # コンソールにも全文を出す (ファイルだけだと分析しづらいという指摘に対応)。
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
                # top1 のみコンソールにも全文出力 (sub1=旧main の最終詳細レポート相当)。
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
    """run subdir に README.md を出力 (結果の見方ガイド、sub1=旧main の _write_readme 相当)。"""
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    feasible = [t for t in complete if t.user_attrs.get('is_feasible', False)]
    top1_name = os.path.basename(saved_reports[0]) if saved_reports else None

    L = []
    L.append(f"# main.py (HYSYS+SM 全21変数 BO) run — {ts}")
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
    L.append(f"- SAMPLER = tpe, SEED = {SEED}, N_WORKERS = {N_WORKERS}")
    L.append(f"- 探索変数数 = {len(SEARCH_SPACE)} (反応器4 + PSA3 + 膜2 + 原料1 + Dist1/2/3 各4/4/3)")
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
    print(f"==== main.py: PDH HYSYS+SM 全21変数 制約付き BO ====", flush=True)
    print(f"  N_TRIALS={N_TRIALS}, N_STARTUP(QMC)={N_STARTUP}, seed={SEED}, top-k report={N_TOPK}", flush=True)
    print(f"  Dist1/Dist3 = SM, Dist2 = HYSYS / 上流(反応器・PSA・膜)・F_fresh も変数化", flush=True)
    # 設計判断 (2026-05-28): sub1(旧main) 流に run ごとの subdir へ成果物を集約。
    # outputs/main_<ts>/{trials.csv, best.json, top*.txt, optuna.db, README.md}
    out_dir = os.path.join(OUTPUT_DIR, f'main_{ts}')
    os.makedirs(out_dir, exist_ok=True)

    sampler = make_sampler('tpe', SEED, N_STARTUP, constraints_func=constraints_func)
    # 設計判断 (2026-05-26): 並列(N_WORKERS>1)時は worker 間で study 共有のため SQLite 必須。
    _use_sqlite = USE_SQLITE_STORAGE or N_WORKERS > 1
    _db_path = os.path.join(out_dir, 'optuna.db')
    storage = f"sqlite:///{_db_path}" if _use_sqlite else None
    study = optuna.create_study(
        study_name=f"{STUDY_NAME}_{ts}" if _use_sqlite else STUDY_NAME,
        sampler=sampler, direction='minimize',
        storage=storage, load_if_exists=bool(storage),
    )

    if N_WORKERS > 1 and storage is not None:
        # 設計判断 (2026-05-26): N worker プロセスで共有 study を分担。各 worker は自前の
        # HYSYS を起動 (重い)。完了後は study.trials(DB) を読んでサマリ/レポートを生成する。
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
        # 設計判断 (2026-05-26): 生 params の羅列 (float 21 行) は撤去。見やすい入力
        # スナップショットは下の top1 詳細レポート (show_input_snapshot) に出力され、
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

    # ---- 成果物サマリ (sub1=旧main 相当) ----
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
