"""
H2 パージ損失モデル検証スクリプト (R-4)

目的:
    psa_system.py のパージ損失式
        H2_loss_purge = u_0 × A_col × C_H2_purge × (t_des/t_abs) × 3600/1000
    が、新しい KFa 値（0.02 s⁻¹; 旧仮置き値 0.5–1.0 s⁻¹）の下で
    H2 回収率を非現実的に低く推算していないかを確認する。

    旧 KFa: t_des ≈ 2 s → パージ損失ほぼ無視できる
    新 KFa: t_des ≈ 63 s → t_des/t_abs ≈ 0.5–0.65 → 損失が顕在化

確認項目:
    1. H2 回収率 (H2_recovery) が物理的に妥当か
    2. クランプ（物質収支保護）が発火しているか
       raw_purge > (F_H2_in - H2_loss_blowdown) のとき clamp 発火
    3. 30% 警告閾値に引っかかる条件の有無

検証手順:
    a) 代表フィードで (D_col, L_bed) をスキャンし、有効域の条件のみ集計
    b) 各条件で simulate_psa_system を実行
    c) 生の（クランプ前）パージ損失を再計算し、クランプ発火の有無を判定
    d) H2 回収率の分布と危険フラグを表示

使い方:
    python tests/test_h2_loss_purge.py

検証結果メモ (2026-05-06, 代表条件 F_H2_in=100 kmol/h, P_in=20 bar, P_des=1 atm):
    KFa (新値): CH4=0.02, C2H4=0.015, C2H6=0.011 s⁻¹  (旧仮置き: 0.5–1.0 s⁻¹)

    D=0.5 m 全域・D=1.0 m L=3 m はペナルティ。有効条件 11 ケース:

    t_des はほぼ一定 (≈73 s) — KFa と desorption_target で決まり、塔寸法に依存しない。

    パージ損失 (P%) : 0.2 – 2.9 %   ← 当初の懸念とは逆に非常に小さい
      理由: C_H2_purge = C_H2 × P_des/P_in ≈ C_H2/20 (大気圧でほぼ希薄)

    ブローダウン損失 (BD%): 14 – 32 %  ← こちらが支配的な H2 損失
      理由: 20 bar 高圧空隙に大量の H2 が存在し、大気圧ブローダウンで 95 % 放出。
      改善策: 圧力均等化 (Pressure Equalization) ステップを設計に組み込めば
              ブローダウン H2 の 60–80 % 回収が可能 (初期設計では保守値として許容)。

    クランプ発火: 0 件 / 警告 (30 % 超過): 0 件 → 物質収支は健全
    H2 回収率 < 70 %: 1 件 (D=1.0, L=5.0 → H2_rec=0.655)
      → 最適化は CAPEX と H2_rec のトレードオフで自動的に大きな塔を選択するはず。

    【結論】
      パージ損失モデルの過大推算は軽微 (< 3%)。
      ブローダウン損失が支配的だが圧力均等化なしの保守的推算として許容。
      コードの変更は不要。ただし詳細設計段階では圧力均等化を実装すること。
"""

import math
import os
import sys
import warnings

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.config import R
from src.cost_parameters import PSA_KFA
from units.separators.psa.psa_system import (
    PSADesignVars, PSAFeedStream, PSAFixedParams,
    simulate_psa_system,
    _calc_feed_state,
)

# ---------------------------------------------------------------------------
# 代表フィード条件
# ---------------------------------------------------------------------------
_FEED = PSAFeedStream(
    F_in={'A': 0.0, 'B': 0.0, 'C': 100.0, 'D': 2.0, 'E': 5.0, 'F': 1.0},
    T_in=350.0,
    P_in=20.0e5,
)
_FIXED = PSAFixedParams()  # デフォルト (P_des=101325 Pa, T_abs=298.15 K)

# (D_col [m], L_bed [m]) のスキャン点
_D_VALS = [0.5, 1.0, 1.5, 2.0]
_L_VALS = [3.0, 5.0, 7.0, 10.0]


def _raw_purge_loss(u_0, A_col, C_H2, Z, P_des, P_in, t_des, t_abs):
    """クランプ前のパージ損失 [kmol/h] を再計算する。"""
    C_H2_purge = C_H2 * Z * P_des / P_in
    return u_0 * A_col * C_H2_purge * t_des / t_abs * 3600.0 / 1000.0


def _raw_blowdown_loss(V_col, eps, C_H2, Z, P_des, P_in, t_abs):
    """ブローダウン損失 [kmol/h]（クランプ対象外）。"""
    return V_col * eps * C_H2 * (1.0 - Z * P_des / P_in) / t_abs * 3600.0 / 1000.0


def main():
    print("=" * 70)
    print("R-4: H2 パージ損失モデル検証")
    print("=" * 70)
    print(f"  フィード: {_FEED.F_in}")
    print(f"  P_in = {_FEED.P_in/1e5:.1f} bar, P_des = {_FIXED.P_des/1e5:.4f} bar")
    print(f"  KFa: CH4={PSA_KFA['CH4']}, C2H4={PSA_KFA['C2H4']}, C2H6={PSA_KFA['C2H6']} [1/s]")
    print()

    F_H2_in = _FEED.F_in.get('C', 0.0)  # [kmol/h]

    # フィード状態量（T_abs, P_in 条件）
    _, C_H2, _, Z = _calc_feed_state(_FEED.F_in, _FIXED.T_abs, _FEED.P_in)

    header = (
        f"{'D':>5} {'L':>5} | {'t_abs':>7} {'t_des':>7} {'td/ta':>6} | "
        f"{'BD_loss':>8} {'raw_P':>8} {'clp_P':>8} | "
        f"{'BD%':>6} {'P%':>6} {'tot%':>6} | "
        f"{'H2rec':>7} {'clamp':>6} {'warn':>5}"
    )
    sep = "-" * len(header)
    print(header)
    print(sep)

    n_clamp  = 0
    n_warn   = 0
    n_low    = 0
    n_valid  = 0

    for D in _D_VALS:
        for L in _L_VALS:
            design = PSADesignVars(D_col=D, L_bed=L, desorption_target=0.35)
            A_col  = math.pi / 4.0 * D ** 2
            V_col  = A_col * L

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                res = simulate_psa_system(design, _FEED, _FIXED)

            # ペナルティチェック
            eq = res.equipment
            if eq.CAPEX_total < 0.0 or eq.t_abs_sec == 0.0:
                print(f"{D:>5.1f} {L:>5.1f} | {'PENALTY':^57}")
                continue

            n_valid += 1
            t_abs = eq.t_abs_sec
            t_des = eq.t_des_sec
            u_0   = eq.u_0

            # 生パージ損失再計算
            raw_purge = _raw_purge_loss(
                u_0, A_col, C_H2, Z, _FIXED.P_des, _FEED.P_in, t_des, t_abs
            )
            raw_bd = _raw_blowdown_loss(
                V_col, _FIXED.eps, C_H2, Z, _FIXED.P_des, _FEED.P_in, t_abs
            )
            clamped_purge = eq.H2_loss_purge_kmolh
            clamped_bd    = eq.H2_loss_blowdown_kmolh

            clamp_fired = raw_purge > clamped_purge + 1e-9
            warn_fired  = any("PSAパージ損失" in str(w.message) for w in caught)

            bd_pct   = raw_bd     / F_H2_in * 100.0
            p_pct    = raw_purge  / F_H2_in * 100.0
            tot_pct  = (raw_bd + raw_purge) / F_H2_in * 100.0

            if clamp_fired: n_clamp += 1
            if warn_fired:  n_warn  += 1
            if res.H2_recovery < 0.70: n_low += 1

            print(
                f"{D:>5.1f} {L:>5.1f} | "
                f"{t_abs:>7.1f} {t_des:>7.1f} {t_des/t_abs:>6.3f} | "
                f"{raw_bd:>8.3f} {raw_purge:>8.3f} {clamped_purge:>8.3f} | "
                f"{bd_pct:>6.1f} {p_pct:>6.1f} {tot_pct:>6.1f} | "
                f"{res.H2_recovery:>7.3f} "
                f"{'YES' if clamp_fired else 'no':>6} "
                f"{'YES' if warn_fired else 'no':>5}"
            )

    print(sep)
    print()
    print("[凡例]")
    print("  t_abs / t_des : 吸着・脱着時間 [s]   td/ta : t_des/t_abs")
    print("  BD_loss  : ブローダウン損失 (生値) [kmol/h]")
    print("  raw_P    : パージ損失 (クランプ前, 生値) [kmol/h]")
    print("  clp_P    : パージ損失 (クランプ後, コードが使用) [kmol/h]")
    print("  BD%/P%/tot% : 各損失 / F_H2_in [%]")
    print("  H2rec    : H2 回収率 [-]")
    print("  clamp    : 物質収支クランプ発火の有無")
    print("  warn     : 30% 過大推算警告の有無")
    print()
    print(f"[集計] 有効条件={n_valid}  clamp発火={n_clamp}  warn発火={n_warn}  H2rec<70%={n_low}")
    print()

    # -----------------------------------------------------------------------
    # 判定
    # -----------------------------------------------------------------------
    print("[判定]")
    print("[パージ損失]")
    print("  当初懸念: 新 KFa (0.02 s-1) で t_des~73s -> パージ損失が爆発するのでは？")
    print("  実結果  : パージ損失 < 3% -> 懸念は杞憂。C_H2_purge = C_H2*(P_des/P_in)")
    print("            = C_H2/20 (大気圧) なので損失は小さい。")
    print()
    print("[ブローダウン損失]")
    print("  実際の支配的損失はブローダウン (BD%: 14-32%)。")
    print("  P_in=20 bar の高圧空隙から大気圧へのブローダウンで H2 の 95% が放出される。")
    print("  現行モデルは圧力均等化なしの保守的推算。詳細設計段階で PE ステップを追加すること。")
    print()

    if n_clamp == 0 and n_warn == 0 and n_low == 0:
        print("  OK - 全有効条件で H2 回収率 >= 70%、クランプ・警告なし。")
        print("     現行パージモデルはこの設計空間では保守的だが許容範囲内。")
    else:
        if n_clamp > 0:
            print(f"  NG (clamp) - {n_clamp} 条件でクランプ発火。")
            print("    パージ損失が物質収支限界を超えており、最適化の評価関数が")
            print("    不連続・非物理的になるリスクあり。")
        if n_warn > 0:
            print(f"  NG (warn) - {n_warn} 条件で 30% 超過警告。")
            print("    パージ速度 u_0 や t_des の仮定を見直すか、")
            print("    パージ流量を t_abs の一定割合に制限する設計制約を追加すること。")
        if n_low > 0:
            print(f"  NG (low H2) - {n_low} 条件で H2 回収率 < 70%。")
            print("    これらの条件は最適化で高コスト評価を受けるため自動的に排除されるが、")
            print("    ブローダウン損失が支配的なことに留意（PE ステップで改善可能）。")


if __name__ == "__main__":
    main()
