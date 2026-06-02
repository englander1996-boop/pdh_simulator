"""SM (surrogate model) バックエンド — Dist1 / Dist3 を学習済み GPR で置換する。

HYSYS バックエンドの遅さは「HSC を毎パス開き直す COM オーバーヘッド」が支配項
(塔の計算自体は ~0.3s と高速)。Dist1/Dist3 を学習済みサロゲートに置換すると
当該塔の HYSYS 呼び出しが消え、HYSYS 塔が Dist2 のみになって swap が原理的に消滅する。
(SM は入力→出力の決定的写像なので経路依存はない。)

モデル (models/column{1,3}_sm.pkl):
  dict 形式。`regressors`=ターゲット別 GaussianProcessRegressor、
  `classifier`=RandomForest Pipeline (feasibility 判定)、
  `x_scaler`/`y_scalers`=StandardScaler。`target_scales` は y_scaler.scale_ と同値で冗長。

予測手順 (HYSYS 実測との突合せで確定):
  Xs = x_scaler.transform(X)
  各 target: y = y_scalers[t].inverse_transform(regressors[t].predict(Xs))
  feasibility: classifier.predict(X)  ← raw X (Pipeline 内蔵 scaler)

入力 (input_columns):
  Dist1: In_Total_Stages, In_Feed_Stage, In_Column_P[kPa], In_Flow[kgmol/h], In_CompFraction2(spec)
  Dist3: In_Total_Stages, In_Feed_Stage, In_Column_P[kPa], In_Flow[kgmol/s], In_Propane(feed C3H8 mol frac)

分配 (回収率) の扱い:
  model3 は spec 入力を持たず、サンプリング時の固定 spec 戦略を埋め込んでいる。
  Dist3 SM の分配はそのまま設計値として「正」とする。
  製品純度 (Out_TopPropene) は spec を満たすので製品仕様上 OK、回収率は SM 準拠。

出力 → DistResult は HYSYS 経路と同じ provider._column_result_to_dist_result を再利用し、
utility 選択・HE 面積・CAPEX 計算を完全に揃える。組成は binary feed 前提
(Dist1=A/Z, Dist3=A/B) で mass balance (bottom = feed − top) により復元し保存則を守る。
"""
from __future__ import annotations

import os
import pickle
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from stream.stream import ProcessStream
from src.distillation_core import ColumnTunables, DistFixedParams, DistResult

# DistResult 組立・penalty は HYSYS 経路と共通ロジックを再利用 (utility/HE/CAPEX 一致)。
# provider の import は win32com を引くが、exp3 は Dist2 を HYSYS で解くため問題なし。
from units.vle.hysys.provider import (
    _column_result_to_dist_result, _failure_dist_result,
)
from units.vle.hysys.adapters.types import ColumnResult


_MODELS_DIR = Path(__file__).resolve().parents[1] / "models"
_MODEL_FILES = {
    "column1": _MODELS_DIR / "column1_sm.pkl",
    "column3": _MODELS_DIR / "column3_sm.pkl",
}

# 学習済みモデルのプロセス内キャッシュ (pkl ロードは数百 ms かかるので 1 回だけ)
_MODEL_CACHE: Dict[str, dict] = {}


def _load_model(column_id: str) -> dict:
    if column_id not in _MODEL_CACHE:
        path = _MODEL_FILES[column_id]
        if not path.exists():
            raise FileNotFoundError(f"SM モデルが見つからない: {path}")
        with open(path, "rb") as f:
            _MODEL_CACHE[column_id] = pickle.load(f)
    return _MODEL_CACHE[column_id]


def _predict(model: dict, x: Dict[str, float]) -> Tuple[bool, Dict[str, float], bool]:
    """SM 予測。

    Returns
    -------
    (feasible, outputs, clamped)
        feasible : classifier の判定 (True=実行可能)
        outputs  : target_columns → 予測値 [物理単位]
        clamped  : 入力が bounds 外で clip されたら True
    """
    cols = model["input_columns"]
    bounds = model.get("bounds", {})
    clamped = False
    vals = []
    for c in cols:
        v = float(x[c])
        lo, hi = bounds.get(c, (None, None))
        if lo is not None and v < lo:
            v = lo; clamped = True
        elif hi is not None and v > hi:
            v = hi; clamped = True
        vals.append(v)
    X = np.asarray([vals], dtype=float)

    with warnings.catch_warnings():
        # scaler/classifier は DataFrame (feature 名付き) で学習されているため、
        # ndarray を渡すと sklearn が feature-name 警告を出す。予測値には影響しないので抑制。
        warnings.simplefilter("ignore")
        Xs = model["x_scaler"].transform(X)
        try:
            feasible = bool(model["classifier"].predict(X)[0])
        except Exception:
            feasible = True  # 分類器が無い/失敗時は feasible 扱い (回帰のみで判断)
        out: Dict[str, float] = {}
        for t in model["target_columns"]:
            ps = np.asarray(model["regressors"][t].predict(Xs)).reshape(-1, 1)
            ys = model["y_scalers"][t]
            try:
                val = float(ys.inverse_transform(ps)[0, 0])
            except Exception:
                # y_scaler が非スケール (scale_=None) のターゲットは恒等
                val = float(ps[0, 0])
            out[t] = val
    return feasible, out, clamped


def _binary_split_streams(
    feed: ProcessStream,
    top_total: float,
    top_key: str,
    top_key_frac: float,
    other_key: str,
) -> Tuple[Dict[str, float], float, Dict[str, float], float]:
    """binary feed (top_key / other_key) を top/bottom に分配 (mass balance 厳守)。

    top は SM 予測の (top_total, top_key_frac) を採用。bottom = feed − top で保存則を守る。
    成分が feed を超えないよう clip し、top_total を再計算する。

    Returns (top_comp_frac, top_flow, bot_comp_frac, bot_flow)
    """
    feed_top_key   = max(feed.F_in.get(top_key, 0.0),   0.0)
    feed_other_key = max(feed.F_in.get(other_key, 0.0), 0.0)

    top_kk = top_total * top_key_frac
    top_ok = top_total * (1.0 - top_key_frac)
    # feed を超えないよう clip (SM の微小なマスバランス誤差・外挿対策)
    top_kk = min(max(top_kk, 0.0), feed_top_key)
    top_ok = min(max(top_ok, 0.0), feed_other_key)
    top_flow = top_kk + top_ok

    bot_kk = feed_top_key - top_kk
    bot_ok = feed_other_key - top_ok
    bot_flow = bot_kk + bot_ok

    top_frac = ({top_key: top_kk / top_flow, other_key: top_ok / top_flow}
                if top_flow > 0 else {top_key: 0.0, other_key: 0.0})
    bot_frac = ({top_key: bot_kk / bot_flow, other_key: bot_ok / bot_flow}
                if bot_flow > 0 else {top_key: 0.0, other_key: 0.0})
    return top_frac, top_flow, bot_frac, bot_flow


def _build_column_result(
    top_frac: Dict[str, float], top_flow: float, top_temp_c: float,
    bot_frac: Dict[str, float], bot_flow: float, bot_temp_c: float,
    qc_mw: float, qr_mw: float, p_col_kpa: float,
    reflux_ratio: Optional[float], message: str,
) -> ColumnResult:
    cr = ColumnResult()
    cr.success = True
    cr.message = message
    cr.qc_mw = float(qc_mw)
    cr.qr_mw = float(qr_mw)
    cr.top_flow_kmolh = float(top_flow)
    cr.top_temp_c = float(top_temp_c)
    cr.top_pressure_kpa = float(p_col_kpa)
    cr.top_composition = dict(top_frac)
    cr.bot_flow_kmolh = float(bot_flow)
    cr.bot_temp_c = float(bot_temp_c)
    cr.bot_pressure_kpa = float(p_col_kpa)
    cr.bot_composition = dict(bot_frac)
    cr.reflux_ratio = reflux_ratio
    return cr


# ---------------------------------------------------------------------------
# 塔別エントリ (column1.py / column3.py の wrapper から呼ばれる)
# ---------------------------------------------------------------------------

def solve_column1_via_sm(
    feed: ProcessStream,
    tunables: ColumnTunables,
    fixed: Optional[DistFixedParams] = None,
) -> DistResult:
    """Dist1 (脱ブタン塔) を SM で解く。

    In_CompFraction2 = tunables.hysys_spec_value (Comp Fraction-2 spec)。
    top = C3H8(A) 富化、bottom = C4H10(Z) 富化 (binary feed)。
    """
    if tunables.hysys_feed_stage is None or tunables.hysys_spec_value is None:
        return _failure_dist_result(tunables, "SM: hysys_feed_stage / hysys_spec_value が None")
    model = _load_model("column1")
    feed_total_kmolh = sum(max(v, 0.0) for v in feed.F_in.values())
    x = {
        "In_Total_Stages":  int(tunables.N_stages),
        "In_Feed_Stage":    int(tunables.hysys_feed_stage),
        "In_Column_P":      float(tunables.P_col) / 1000.0,   # Pa → kPa
        "In_Flow":          feed_total_kmolh,                  # kgmol/h
        "In_CompFraction2": float(tunables.hysys_spec_value),
    }
    feasible, o, clamped = _predict(model, x)
    if not feasible:
        return _failure_dist_result(tunables, "SM(Dist1): classifier が infeasible 判定")
    # Out_TopPropane = 塔頂 C3H8(A) mol 分率。top_key=A, other=Z(C4H10)。
    top_frac, top_flow, bot_frac, bot_flow = _binary_split_streams(
        feed, top_total=o["Out_TopFlow"], top_key="A",
        top_key_frac=o["Out_TopPropane"], other_key="Z",
    )
    cr = _build_column_result(
        top_frac, top_flow, o["Out_TopTemp_C"],
        bot_frac, bot_flow, o["Out_BtmTemp_C"],
        o["Out_Qc_MW"], o["Out_Qr_MW"], x["In_Column_P"],
        o.get("Out_RefluxRatio"),
        message="sm" + (" (input clamped to bounds)" if clamped else ""),
    )
    return _column_result_to_dist_result(cr, tunables, feed, "column1")


def solve_column3_via_sm(
    feed: ProcessStream,
    tunables: ColumnTunables,
    fixed: Optional[DistFixedParams] = None,
) -> DistResult:
    """Dist3 (C3 スプリッタ) を SM で解く。

    In_Propane = feed の C3H8(A) mol 分率 (spec 入力ではない)。
    top = C3H6(B) 製品、bottom = C3H8(A) 富化リサイクル (binary feed)。
    分配 (回収率) は SM 予測 (Out_TopFlow, Out_TopPropene) をそのまま設計値として採用。
    """
    if tunables.hysys_feed_stage is None:
        return _failure_dist_result(tunables, "SM: hysys_feed_stage が None")
    model = _load_model("column3")
    feed_total_kmolh = sum(max(v, 0.0) for v in feed.F_in.values())
    if feed_total_kmolh <= 0:
        return _failure_dist_result(tunables, "SM(Dist3): feed flow ≤ 0")
    feed_A = max(feed.F_in.get("A", 0.0), 0.0)   # C3H8
    propane_frac = feed_A / feed_total_kmolh
    x = {
        "In_Total_Stages": int(tunables.N_stages),
        "In_Feed_Stage":   int(tunables.hysys_feed_stage),
        "In_Column_P":     float(tunables.P_col) / 1000.0,    # Pa → kPa
        "In_Flow":         feed_total_kmolh / 3600.0,          # kgmol/s
        "In_Propane":      propane_frac,
    }
    feasible, o, clamped = _predict(model, x)
    if not feasible:
        return _failure_dist_result(tunables, "SM(Dist3): classifier が infeasible 判定")
    # Out_TopPropene = 塔頂 C3H6(B) mol 分率。top_key=B, other=A(C3H8)。
    top_frac, top_flow, bot_frac, bot_flow = _binary_split_streams(
        feed, top_total=o["Out_TopFlow"], top_key="B",
        top_key_frac=o["Out_TopPropene"], other_key="A",
    )
    cr = _build_column_result(
        top_frac, top_flow, o["Out_TopTemp_C"],
        bot_frac, bot_flow, o["Out_BtmTemp_C"],
        o["Out_Qc_MW"], o["Out_Qr_MW"], x["In_Column_P"],
        o.get("Out_RefluxRatio"),
        message="sm" + (" (input clamped to bounds)" if clamped else ""),
    )
    return _column_result_to_dist_result(cr, tunables, feed, "column3")


# ---------------------------------------------------------------------------
# generic dispatch 互換 (simulate_distillation_column の 'sm' 経路用)
# ---------------------------------------------------------------------------

def simulate_sm(design, feed, fixed=None):
    """src.distillation_core の 'sm' dispatch 用。LK/HK で塔を識別して委譲。

    ただし DistDesignVars は hysys_spec_value/hysys_feed_stage を持たないため、
    通常は column1.py / column3.py の wrapper 側で solve_column*_via_sm を直接呼ぶ。
    本関数は後方互換のスタブ (LK/HK から判別できない構成では NotImplementedError)。
    """
    raise NotImplementedError(
        "simulate_sm: generic dispatch は未対応。column{1,3}.py の "
        "solve_column{1,3}_via_sm を使うこと (hysys_spec_value/feed_stage が必要)。"
    )
