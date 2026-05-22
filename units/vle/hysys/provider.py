"""HYSYS VLE プロバイダ。

役割:
  - ColumnTunables (PDH 側設計変数) + ProcessStream (PDH 側 feed) を
    Column{1,2,3}Input に組み立て、HysysSession + Column{1,2,3}Adapter で実行
  - HYSYS の出口値 (ColumnResult) を DistResult に変換 (CAPEX は cost_calculator
    で後計算)
  - 非収束/エラーは feasible=False の penalty DistResult で吸収 → run_one_pass
    側の penalty_hit 経路に乗せる (BO 探索から自動的に除外)

特徴:
  - column*.py の simulate_column* から呼ばれる。tunables.solver_method == 'hysys'
    のときに dispatch される
  - 段数解決は HsysRegistry (デフォルトは hysys_cases/ 配下)
  - 同じ段数 HSC が連続して呼ばれる場合に備えて LRU 1スロットのセッションキャッシュ
    を備える (BO 中に同段数連続のとき HYSYS 再起動を回避)
"""
from __future__ import annotations

import math
import warnings
from typing import Dict, Optional

from stream.stream import ProcessStream
from src.cost_calculator import calc_he_capex_okuyen
from src.cost_parameters import PENALTY_CAPEX_THRESHOLD_OKUYEN
from src.utility_selector import select_cooling_utility, select_heating_utility
from flowsheet.heat_integration import StreamPhase, lookup_U, utility_phase
from src.distillation_core import (
    ColumnTunables, DistDesignVars, DistFixedParams, DistResult, DistEquipment,
    _supply_T_for_utility,
)

from units.vle.hysys.session import HysysSession
from units.vle.hysys.registry import HsysRegistry, StageNotAvailableError
from units.vle.hysys.adapters import (
    Column1Adapter, Column2Adapter, Column3Adapter,
    Column1Input, Column2Input, Column3Input, ColumnResult,
)

_ALL_PDH_KEYS = ('A', 'B', 'C', 'D', 'E', 'F', 'Z')
_ZERO_FLOW = {k: 0.0 for k in _ALL_PDH_KEYS}
_T_COOLING_FLOOR_K = 173.15 + 10.0 + 1.0  # = -89°C: -100°C 冷媒 + 10K approach + 1K margin
_SENTINEL_CAPEX = max(PENALTY_CAPEX_THRESHOLD_OKUYEN, 1e8)


# ---------------------------------------------------------------------------
# セッションキャッシュ (LRU 1 スロット)
# ---------------------------------------------------------------------------

class _SessionCache:
    """3 セッション同時保持型キャッシュ (app 1 つで複数 case 保持)。

    設計判断 (2026-05-22): フローシート全体評価で Dist1/Dist2/Dist3 を順次呼ぶとき、
    case swap (Close + Open) するだけでも sleep 含めて 1塔 ~10秒。3塔で 30秒/iter。
    → 3 つの case を全て開きっぱなしにして、Activate で切替えるだけにすれば、
    各 iter で HSC Open/Close 不要になり大幅高速化。

    仕様:
      - app は HysysSession 1 つで保持
      - case は (column_id, n_stages) ごとに別オブジェクトで dict 保持
      - acquire() で既知の key なら Activate のみ、新規なら Open + Activate
      - release() で全 case を Close + app Quit
      - 異なる n_stages が要求されたら別 case として追加で Open (蓄積)
        → BO で N_stages が散らばっても問題ないが、メモリは増える。
          max_cases を超えたら LRU で古い case を Close。
    """
    def __init__(self, visible: bool = False, keep_open: bool = False,
                 max_cases: int = 6):
        self.visible = visible
        self.keep_open = keep_open
        self.max_cases = max_cases
        self._app_session: Optional[HysysSession] = None
        # OrderedDict 相当: 古い順に並ぶ
        self._cases: "dict[tuple, object]" = {}

    def acquire(self, column_id: str, n_stages: int, hsc_path) -> HysysSession:
        import time as _t
        key = (column_id, int(n_stages))

        # 初回: app を起動して最初の case を Open
        if self._app_session is None:
            self._app_session = HysysSession(
                hsc_path, visible=self.visible, keep_open=self.keep_open,
            ).__enter__()
            self._cases[key] = self._app_session._case
            return self._app_session

        # 設計判断 (2026-05-22 改訂): 3 セッション同時保持で Dist3 が iter 2 以降で
        # Solver 走らない症状が出たので、swap_case 経由に戻す (app 使い回しのみ維持)。
        # 同 key なら現セッションそのまま、異なる key なら HSC を swap。
        if key == next(iter(self._cases), None):
            return self._app_session

        # 異なる key: swap で切替 (debug で動作実績ある経路)
        try:
            self._app_session.swap_case(hsc_path)
            # _cases を最新の key だけに更新 (古い key は捨てる、case は Close 済み)
            self._cases.clear()
            self._cases[key] = self._app_session._case
            return self._app_session
        except Exception:
            # swap 失敗時は app ごと再起動
            self.release()
            return self.acquire(column_id, n_stages, hsc_path)

    def release(self) -> None:
        if self._app_session is None:
            return
        # 全 case を Close
        for case in self._cases.values():
            try:
                case.Close(False)
            except Exception:
                pass
        self._cases = {}
        try:
            self._app_session.__exit__(None, None, None)
        except Exception:
            pass
        self._app_session = None


# ---------------------------------------------------------------------------
# HysysVLEProvider
# ---------------------------------------------------------------------------

class HysysVLEProvider:
    """3塔を HYSYS で解く高水準 API。

    通常は process-global な singleton として利用する (_default_provider)。
    """
    def __init__(
        self,
        registry: Optional[HsysRegistry] = None,
        visible: bool = False,
        timeout_sec: float = 180.0,
        use_cache: bool = True,
        keep_open: bool = False,
    ):
        """
        keep_open : True のとき close() / __exit__ でも HYSYS を残す (デバッグ用)。
        """
        self.registry = registry if registry is not None else HsysRegistry()
        self.visible = visible
        self.timeout_sec = timeout_sec
        self.keep_open = keep_open
        self._cache = (
            _SessionCache(visible=visible, keep_open=keep_open) if use_cache else None
        )
        self._run_counter = 0

    # ---- ライフサイクル ----
    def close(self) -> None:
        if self._cache is not None:
            self._cache.release()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False

    # ---- 内部: セッション取得 (キャッシュ経由 or 単発) ----
    def _open_session(self, column_id: str, n_stages: int):
        """(session, owns_session) を返す。owns_session=True なら呼出側で閉じる。"""
        try:
            path = self.registry.get_path(column_id, n_stages)
        except StageNotAvailableError as e:
            raise
        if self._cache is not None:
            return self._cache.acquire(column_id, n_stages, path), False
        sess = HysysSession(
            path, visible=self.visible, keep_open=self.keep_open,
        ).__enter__()
        return sess, True

    def _next_run_id(self) -> int:
        self._run_counter += 1
        return self._run_counter

    # ---- 各塔の solve API ----
    def solve_column1(self, feed: ProcessStream, tunables: ColumnTunables) -> DistResult:
        return _solve_via_hysys(self, feed, tunables, column_id="column1")

    def solve_column2(self, feed: ProcessStream, tunables: ColumnTunables) -> DistResult:
        return _solve_via_hysys(self, feed, tunables, column_id="column2")

    def solve_column3(self, feed: ProcessStream, tunables: ColumnTunables) -> DistResult:
        return _solve_via_hysys(self, feed, tunables, column_id="column3")


# ---------------------------------------------------------------------------
# プロセス global な default provider
# ---------------------------------------------------------------------------

_default_provider: Optional[HysysVLEProvider] = None


def get_default_provider() -> HysysVLEProvider:
    """プロセス内で再利用される default の HysysVLEProvider。"""
    global _default_provider
    if _default_provider is None:
        _default_provider = HysysVLEProvider()
    return _default_provider


def shutdown_default_provider() -> None:
    """default provider をクローズ (テスト・終了処理用)。"""
    global _default_provider
    if _default_provider is not None:
        _default_provider.close()
        _default_provider = None


# ---------------------------------------------------------------------------
# column*.py から呼ばれる薄いエントリ
# ---------------------------------------------------------------------------

def solve_column1_via_hysys(
    feed: ProcessStream,
    tunables: ColumnTunables,
    fixed: Optional[DistFixedParams] = None,
) -> DistResult:
    """column1.py の simulate_column1 から呼ばれる HYSYS 経路。"""
    return get_default_provider().solve_column1(feed, tunables)


def solve_column2_via_hysys(
    feed: ProcessStream,
    tunables: ColumnTunables,
    fixed: Optional[DistFixedParams] = None,
) -> DistResult:
    """column2.py の simulate_column2 から呼ばれる HYSYS 経路。"""
    return get_default_provider().solve_column2(feed, tunables)


def solve_column3_via_hysys(
    feed: ProcessStream,
    tunables: ColumnTunables,
    fixed: Optional[DistFixedParams] = None,
) -> DistResult:
    """column3.py の simulate_column3 から呼ばれる HYSYS 経路。"""
    return get_default_provider().solve_column3(feed, tunables)


# ---------------------------------------------------------------------------
# 共通実行ロジック (3塔分岐)
# ---------------------------------------------------------------------------

def _solve_via_hysys(
    provider: HysysVLEProvider,
    feed: ProcessStream,
    tunables: ColumnTunables,
    column_id: str,
) -> DistResult:
    # 1. HYSYS 用入力を組み立て
    try:
        hysys_input = _build_input(column_id, feed, tunables)
    except Exception as e:
        return _failure_dist_result(tunables, f"HYSYS input 組立失敗: {e}")

    # 2. HSC オープン
    try:
        sess, owns = provider._open_session(column_id, int(tunables.N_stages))
    except StageNotAvailableError as e:
        return _failure_dist_result(tunables, f"段数解決失敗: {e}")
    except Exception as e:
        return _failure_dist_result(tunables, f"HSC オープン失敗: {e}")

    # 3. アダプタ実行
    try:
        adapter_cls = _ADAPTER_BY_COLUMN_ID[column_id]
        run_id = provider._next_run_id()
        cr = adapter_cls.run(
            sess, hysys_input,
            timeout_sec=provider.timeout_sec, run_id=run_id,
        )
    except Exception as e:
        if owns:
            try:
                sess.__exit__(None, None, None)
            except Exception:
                pass
        return _failure_dist_result(tunables, f"HYSYS 実行失敗: {e}")
    finally:
        if owns:
            try:
                sess.__exit__(None, None, None)
            except Exception:
                pass

    # 4. ColumnResult → DistResult 変換
    return _column_result_to_dist_result(cr, tunables, feed, column_id)


_ADAPTER_BY_COLUMN_ID = {
    "column1": Column1Adapter,
    "column2": Column2Adapter,
    "column3": Column3Adapter,
}


# ---------------------------------------------------------------------------
# 入力組立 (ColumnTunables + ProcessStream → Column{1,2,3}Input)
# ---------------------------------------------------------------------------

def _build_input(column_id: str, feed: ProcessStream, t: ColumnTunables):
    if t.hysys_feed_stage is None:
        raise ValueError("ColumnTunables.hysys_feed_stage が None。HYSYS 経路では必須。")
    if t.hysys_spec_value is None:
        raise ValueError("ColumnTunables.hysys_spec_value が None。HYSYS 経路では必須。")

    # PDH 側の P [Pa] → HYSYS 側の kPa
    col_p_kpa = float(t.P_col) / 1000.0
    # PDH 側の feed は kmol/h で持っている
    feed_total_kmolh = float(sum(max(v, 0.0) for v in feed.F_in.values()))

    # 設計判断 (2026-05-22 改訂): 温度は HSC のデフォルト (Calculated 設定) を維持。
    # 圧力は lhs_column*.py の運用通り「塔圧 + 50 kPa」で書込む (Independent 設定で書ける)。
    # アダプタ側で feed_pressure_kpa=None なら自動的に col_p_kpa + FEED_PRESSURE_MARGIN_KPA
    # を採用するので、ここでは None のまま渡してアダプタに委ねる。
    #   - 温度: 入力しない (HSC 内 Calculated、ユーザー方針 2026-05-22)
    #   - 圧力: 書込む (lhs_column*.py 実績、Independent)
    feed_temp_c = None
    feed_p_kpa  = None    # アダプタが col_p_kpa + 50 で計算する

    if column_id == "column1":
        return Column1Input(
            feed_composition=dict(feed.F_in),
            feed_temperature_c=feed_temp_c,
            feed_pressure_kpa=feed_p_kpa,
            col_p_kpa=col_p_kpa,
            feed_flow_kmolh=feed_total_kmolh,
            comp_frac_2=float(t.hysys_spec_value),
            feed_stage=int(t.hysys_feed_stage),
        )
    if column_id == "column2":
        return Column2Input(
            feed_composition=dict(feed.F_in),
            feed_temperature_c=feed_temp_c,
            feed_pressure_kpa=feed_p_kpa,
            col_p_kpa=col_p_kpa,
            feed_flow_kmolh=feed_total_kmolh,
            reflux_ratio=float(t.hysys_spec_value),
            feed_stage=int(t.hysys_feed_stage),
        )
    if column_id == "column3":
        return Column3Input(
            feed_composition=dict(feed.F_in),
            feed_temperature_c=feed_temp_c,
            feed_pressure_kpa=feed_p_kpa,
            col_p_kpa=col_p_kpa,
            feed_flow_kgmols=feed_total_kmolh / 3600.0,
            propane_frac=None,                    # feed_composition で全成分上書き
            draw_rate_kgmols=float(t.hysys_spec_value),
            feed_stage=int(t.hysys_feed_stage),
        )
    raise ValueError(f"未知の column_id: {column_id}")


# ---------------------------------------------------------------------------
# ColumnResult → DistResult 変換
# ---------------------------------------------------------------------------

def _failure_dist_result(t: ColumnTunables, msg: str) -> DistResult:
    """非収束/エラー時の penalty DistResult。

    run_one_pass の solver_failure 検出経路に乗せるため CAPEX に sentinel を入れる
    (column1.py 内の _penalty_result と整合)。
    """
    return DistResult(
        top=ProcessStream(F_in=dict(_ZERO_FLOW), T_in=298.15, P_in=t.P_col),
        bottom=ProcessStream(F_in=dict(_ZERO_FLOW), T_in=298.15, P_in=t.P_col),
        equipment=DistEquipment(
            D_col=0.0, H_col=0.0, V_col=0.0,
            CAPEX_vessel=_SENTINEL_CAPEX, CAPEX_trays=0.0, CAPEX=_SENTINEL_CAPEX,
            Q_cond=0.0, Q_reb=0.0,
            feasible=False, message=f"hysys: {msg}",
        ),
    )


def _column_result_to_dist_result(
    cr: ColumnResult,
    t: ColumnTunables,
    feed: ProcessStream,
    column_id: str,
) -> DistResult:
    if not cr.success:
        return _failure_dist_result(t, cr.message or "HYSYS 失敗 (理由不明)")

    # mol fraction × 流量 → 各成分流量 [kmol/h]
    top_flows: Dict[str, float] = {
        k: float(cr.top_composition.get(k, 0.0)) * float(cr.top_flow_kmolh)
        for k in _ALL_PDH_KEYS
    }
    bot_flows: Dict[str, float] = {
        k: float(cr.bot_composition.get(k, 0.0)) * float(cr.bot_flow_kmolh)
        for k in _ALL_PDH_KEYS
    }

    T_top_K = float(cr.top_temp_c) + 273.15
    T_bot_K = float(cr.bot_temp_c) + 273.15

    Q_cond_kW = float(cr.qc_mw) * 1000.0
    Q_reb_kW  = float(cr.qr_mw)  * 1000.0

    # utility 選択 + A 計算 + HE CAPEX
    # 設計判断: distillation_core._simulate_rigorous の Step 6-8 と同じロジックで揃える。
    # 失敗時は penalty (sentinel CAPEX) で BO 探索から除外させる。
    try:
        cond_utility = select_cooling_utility(max(T_top_K, _T_COOLING_FLOOR_K), mode='continuous')
        reb_utility  = select_heating_utility(T_bot_K, mode='continuous')
    except ValueError as e:
        return _failure_dist_result(t, f"utility 選択失敗 (T_top={T_top_K-273.15:.1f}°C, T_bot={T_bot_K-273.15:.1f}°C): {e}")

    try:
        A_cond, dT_lm_cond = _compute_he_area(
            T_proc=T_top_K, Q_kW=Q_cond_kW, utility_name=cond_utility.name,
            phase=StreamPhase.CONDENSING, utility_first=True,
        )
        A_reb, dT_lm_reb = _compute_he_area(
            T_proc=T_bot_K, Q_kW=Q_reb_kW, utility_name=reb_utility.name,
            phase=StreamPhase.EVAPORATING, utility_first=False,
        )
    except _HEAreaError as e:
        return _failure_dist_result(t, str(e))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        capex_cond = calc_he_capex_okuyen(A_cond)
        capex_reb  = calc_he_capex_okuyen(A_reb)

    # 塔本体 CAPEX: 暫定で 0。塔径計算 (ρ_v 必要) は後続フェーズで詰める。
    # CAPEX_vessel + CAPEX_trays = 0 でも HE CAPEX + OPEX で BO は方向感を持てる。
    capex_vessel = 0.0
    capex_trays  = 0.0
    capex_total  = capex_vessel + capex_trays + capex_cond + capex_reb

    return DistResult(
        top    = ProcessStream(
            F_in=top_flows, T_in=T_top_K,
            P_in=float(cr.top_pressure_kpa) * 1000.0,
        ),
        bottom = ProcessStream(
            F_in=bot_flows, T_in=T_bot_K,
            P_in=float(cr.bot_pressure_kpa) * 1000.0,
        ),
        equipment = DistEquipment(
            D_col=0.0, H_col=0.0, V_col=0.0,
            CAPEX_vessel=capex_vessel, CAPEX_trays=capex_trays, CAPEX=capex_total,
            Q_cond=Q_cond_kW, Q_reb=Q_reb_kW,
            Q_feed_preheat_kW=0.0,
            N_min=0.0, R_min=0.0,
            N_feed_kirkbride=int(t.hysys_feed_stage) if t.hysys_feed_stage else 0,
            feasible=True,
            message=f"hysys: {cr.message}".strip(),
            A_cond_m2=A_cond, A_reb_m2=A_reb,
            CAPEX_cond=capex_cond, CAPEX_reb=capex_reb,
            cond_utility_name=cond_utility.name,
            cond_utility_jpy_per_GJ=cond_utility.jpy_per_GJ,
            reb_utility_name=reb_utility.name,
            reb_utility_jpy_per_GJ=reb_utility.jpy_per_GJ,
            T_top=T_top_K, T_bot=T_bot_K,
        ),
    )


class _HEAreaError(RuntimeError):
    pass


def _compute_he_area(T_proc: float, Q_kW: float, utility_name: str,
                     phase: StreamPhase, utility_first: bool) -> tuple:
    """HE 伝熱面積 [m²] と LMTD [K] を計算。distillation_core 流用のロジック。

    utility_first=True (condenser): T_proc > T_util, dT1=T_proc-T_in, dT2=T_proc-T_out
    utility_first=False (reboiler): T_util > T_proc, dT_lm = T_util_supply - T_proc
    """
    T_util_in = _supply_T_for_utility(utility_name)
    if utility_first:
        T_util_out = T_util_in + 10.0
        dT1 = T_proc - T_util_in
        dT2 = T_proc - T_util_out
        if dT1 <= 0 or dT2 <= 0:
            raise _HEAreaError(
                f"condenser ΔT 不成立 (T_top={T_proc-273.15:.1f}°C, T_util={T_util_in-273.15:.1f}°C)"
            )
        if abs(dT1 - dT2) < 1e-3:
            dT_lm = 0.5 * (dT1 + dT2)
        else:
            dT_lm = (dT1 - dT2) / math.log(dT1 / dT2)
        U = lookup_U(phase, utility_phase(utility_name))
    else:
        dT_lm = T_util_in - T_proc
        if dT_lm <= 0:
            raise _HEAreaError(
                f"reboiler ΔT 不成立 (T_bot={T_proc-273.15:.1f}°C, T_util={T_util_in-273.15:.1f}°C)"
            )
        U = lookup_U(utility_phase(utility_name), phase)
    A = max(Q_kW * 1000.0 / (U * dT_lm), 10.0)
    return A, dT_lm
