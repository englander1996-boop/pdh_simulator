r"""
comparing.shared.simulator — sim 評価の薄いラッパ + 設定保持 (問題点レポート §4.2 IntegratedSimulator 相当)。

special.py (最適化の本丸) と同じ CONFIG (SM Dist3 の 99.5mol%=99.497wt% を尊重した
purity 99.45wt% 緩和) と評価オプション (apply_hi / hi_dT_min_K / apply_stage2) を 1 箇所に集約。

harness はここの raw_evaluate で FlowsheetResult を取得し、_store_diagnostics と
display_full_results にそのまま渡す (= special.py と同じ忠実な診断・レポート経路)。

バックエンドは space.DEFAULT_BACKEND (Dist2=HYSYS) が既定。HYSYS は単一 COM インスタンス・
並列不可なので、評価は必ず単一プロセス・直列で行うこと (special.py が N_WORKERS=1 の理由と同じ)。
"""

import os
import sys
import dataclasses as _dc

# --- repo root を import path に載せる (どこから import されても動くように) ---
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from config.load import load_operating_config
from flowsheet import evaluate as _flowsheet_evaluate, FlowsheetResult, FlowsheetDesignVars


# 設定は 1 度だけロード。special.py と同じ purity 緩和 (決定A 2026-05-25)。
CONFIG = load_operating_config()
CONFIG = _dc.replace(CONFIG, spec=_dc.replace(CONFIG.spec, c3h6_min_wtfrac=0.9945))

# special.py の既定評価オプション。
APPLY_HI_DEFAULT     = True
HI_DT_MIN_K_DEFAULT  = 10.0
APPLY_STAGE2_DEFAULT = True

# top-k 詳細レポート / best 再評価に渡す既定 kwargs (special.py の eval_kwargs と同形)。
EVAL_KWARGS_DEFAULT = dict(
    apply_hi=APPLY_HI_DEFAULT,
    hi_dT_min_K=HI_DT_MIN_K_DEFAULT,
    apply_stage2=APPLY_STAGE2_DEFAULT,
)


def raw_evaluate(
    design: FlowsheetDesignVars,
    *,
    apply_hi: bool = APPLY_HI_DEFAULT,
    hi_dT_min_K: float = HI_DT_MIN_K_DEFAULT,
    apply_stage2: bool = APPLY_STAGE2_DEFAULT,
    F_fresh: float = None,
) -> FlowsheetResult:
    """design を本シミュレータで評価し FlowsheetResult を返す (special.objective と同じ呼び方)。

    返り値は raw な FlowsheetResult。effective_TAC / is_feasible / economics / specs /
    economics_hi / economics_synth を持ち、_store_diagnostics・display_full_results に
    そのまま渡せる。
    """
    return _flowsheet_evaluate(
        design, CONFIG, verbose=False,
        apply_hi=apply_hi, hi_dT_min_K=hi_dT_min_K, apply_stage2=apply_stage2,
        F_C3H8_override=F_fresh,
    )


def shutdown() -> None:
    """HYSYS provider を後片付け (special.py 終了処理と同じ)。HYSYS 未使用なら無害。"""
    try:
        from units.vle.hysys.provider import shutdown_default_provider
        shutdown_default_provider()
    except Exception:
        pass
