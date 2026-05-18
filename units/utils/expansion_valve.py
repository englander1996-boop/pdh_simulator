"""
膨張弁 (Joule-Thomson 等エンタルピー膨張)

設計判断 (2026-05-08):
  ヒートインテグレーション設計のために、配管中の減圧操作で起こる JT 効果
  (温度低下) を陽に扱う。これまで run_one_pass.py では P を書き換えるだけで
  T を維持していたが、それでは反応器入口プレヒート Q_preheat が過小評価
  される問題があった。

  本モジュールは PR EOS (src/eos.py) を使い、等エンタルピー条件
       H(T_in, P_in) = H(T_out, P_out)
  を満たす T_out を数値的に解く。

仮定 (2026-05-08, 文献根拠なし):
  - 「vapor 相のまま膨張する」(部分気化を扱わない)
    根拠: PDH プロセスでの対象 3 流体 (dist1_top_rx, recycle_dist3,
          recycle_mem) は名目上ガス相。
    限界: 飽和近い流体 (例: 17 bar の dist1_top_rx) は実際には
          17 → 0.5 bar 等エンタルピー膨張で部分気化する可能性があり、
          そのときは潜熱吸収で温度低下が緩和される。本実装は温度低下を
          過大評価する側に偏る。
    将来課題: VLE フラッシュ計算 (Rachford-Rice) を src/eos.py に追加して
          二相膨張に対応する。

依存:
  - src/eos.py: z_factor, residual_enthalpy, _dh_ig
  - src/config.py: THERMO_DATA (PR パラメータ Tc, Pc, omega)
  - 新しい数値定数の導入は **なし** (PR EOS の Ω_a, Ω_b, κ 式は既存)
"""

import os
import sys
import warnings
from typing import Tuple, List

from scipy.optimize import brentq

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from stream.stream import ProcessStream
from src.eos import z_factor, residual_enthalpy, _dh_ig


# 標準成分キー順序 (ProcessStream.F_in のキーに対応)
_DEFAULT_KEYS = ['A', 'B', 'C', 'D', 'E', 'F', 'Z']


def _composition_from_stream(stream: ProcessStream) -> Tuple[List[float], List[str]]:
    """ProcessStream から PR EOS 用の (x, keys) を作る。

    流量 0 の成分は除外する (PR EOS は x_i > 0 のみを期待)。
    """
    F_total = sum(max(F, 0.0) for F in stream.F_in.values())
    if F_total <= 0.0:
        raise ValueError("expansion_valve: 全成分流量がゼロです。")

    keys: List[str] = []
    x:    List[float] = []
    for k in _DEFAULT_KEYS:
        F = stream.F_in.get(k, 0.0)
        if F > 0.0:
            keys.append(k)
            x.append(F / F_total)

    # 微小数値誤差で sum(x) != 1 になる可能性 → 再正規化
    s = sum(x)
    if s > 0:
        x = [xi / s for xi in x]

    return x, keys


def simulate_jt_expansion(
    stream: ProcessStream,
    P_out:  float,
    T_search_lower_K: float = 100.0,
) -> ProcessStream:
    """ProcessStream を等エンタルピー (Joule-Thomson) 膨張させる。

    Parameters
    ----------
    stream : ProcessStream
        入口ストリーム (T_in [K], P_in [Pa], F_in [kmol/h])
    P_out : float
        出口圧力 [Pa]。stream.P_in より低い必要がある (膨張)。
    T_search_lower_K : float
        T_out の探索下限 [K]。デフォルト 100 K (ほぼ全成分の沸点を下回る安全側)。

    Returns
    -------
    ProcessStream
        F_in は同じ (物質収支)、P_in は P_out に、T_in は等エンタルピー条件を
        満たす T_out に置き換わる。

    Raises
    ------
    ValueError
        全成分流量がゼロ、または P_out >= stream.P_in の場合。

    Notes
    -----
    P_out > stream.P_in (= 圧縮) のときは膨張弁モデルとして無効。エラーを返す。
    P_out == stream.P_in のときは入力をそのまま返す (no-op)。
    """
    # ---- 入力チェック ----
    if P_out <= 0:
        raise ValueError(f"expansion_valve: P_out={P_out} は正値でなければなりません。")
    if P_out > stream.P_in:
        raise ValueError(
            f"expansion_valve: P_out ({P_out/1e5:.2f} bar) > P_in "
            f"({stream.P_in/1e5:.2f} bar) は圧縮方向です。膨張弁モデルでは扱えません。"
        )
    if abs(P_out - stream.P_in) < 1e-3:
        # no-op: 圧力差がほぼないので入力をそのまま返す
        return ProcessStream(F_in=dict(stream.F_in), T_in=stream.T_in, P_in=P_out)

    # ---- モル分率と成分キー ----
    x, keys = _composition_from_stream(stream)

    T1 = stream.T_in
    P1 = stream.P_in

    # ---- 入口エンタルピー (vapor 相) ----
    # 設計判断 (2026-05-08): vapor 相を仮定 (上記 docstring の仮定参照)。
    try:
        Z1  = z_factor(T1, P1, x, keys, 'vapor')
        Hr1 = residual_enthalpy(T1, P1, x, keys, Z1)
    except Exception as e:
        warnings.warn(
            f"expansion_valve: 入口エンタルピー計算失敗 ({e})。"
            f" 理想気体近似で T 維持を返す。",
            UserWarning, stacklevel=2,
        )
        return ProcessStream(F_in=dict(stream.F_in), T_in=T1, P_in=P_out)

    # ---- 等エンタルピー条件: H(T2, P_out) = H(T1, P1) ----
    # _dh_ig(T1, T2) は ∫_{T1}^{T2} Cp_ig dT
    # H(T2, P_out) - H(T1, P_in) = _dh_ig(T1, T2) + Hr(T2, P_out) - Hr(T1, P_in)
    def enthalpy_balance(T2: float) -> float:
        try:
            Z2  = z_factor(T2, P_out, x, keys, 'vapor')
            Hr2 = residual_enthalpy(T2, P_out, x, keys, Z2)
        except Exception:
            return float('nan')
        return _dh_ig(T1, T2, x, keys) + (Hr2 - Hr1)

    # JT 効果は通常 T2 < T1 (温度低下)。極稀に逆転気体 (H2 など) で T2 > T1 だが
    # 探索範囲を広めに取る。
    T_lo = max(T_search_lower_K, T1 - 300.0)   # 最大 300 K の温度低下まで想定
    T_hi = T1 + 50.0                            # 逆転 JT 効果のための余裕

    try:
        T2 = brentq(enthalpy_balance, T_lo, T_hi, xtol=0.05, maxiter=200)
    except (ValueError, RuntimeError) as e:
        # 設計判断 (2026-05-18): JT 膨張の brentq 失敗時は理想気体仮定 (T2=T1) で
        # 返す。実気体の JT 冷却効果が無視されるため反応器入口プレヒート Q_preheat
        # が過小推算になる可能性あり (定量はユーザー検証要)。warning は
        # flowsheet/run_one_pass.py の _capture_warnings 経由で BO log に残るため、
        # 頻発時は探索範囲または PR EOS 適用範囲を見直すこと。
        warnings.warn(
            f"expansion_valve: brentq 収束失敗 ({e}, T1={T1:.1f}K, "
            f"P1={P1/1e5:.2f}→P2={P_out/1e5:.2f} bar)。"
            f"理想気体仮定 (T2=T1) で fallback、JT 冷却効果は無視される。",
            UserWarning, stacklevel=2,
        )
        T2 = T1

    return ProcessStream(F_in=dict(stream.F_in), T_in=T2, P_in=P_out)
