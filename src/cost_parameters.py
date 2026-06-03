"""
Bare Module Cost 法の定数定義

【一次出典】
  資料名 : プロセス設計R08-3.pdf
  タイトル: 「プロセス設計(No. 3) 建設費と運転費の推算」
  著者   : 長谷部 伸治・外輪 健一郎

  各数値の引用元ページ:
    p.9  : 付録A Table A.1  — K1, K2, K3（プロセス容器・熱交換器・圧縮機）
    p.10 : 付録B Table A.4  — B1, B2（熱交換器）
    p.11 : 付録C Table A.2  — 圧縮機の圧力補正係数 C1, C2, C3
    p.13–14: Table A.6 および Bare Module Factor グラフ — 圧縮機 FBM

【OPEX式の係数の出典】
  授業資料 §3.4 式 (10): 製造コスト = 0.180·C_TM + 2.73·C_OL
                                    + 1.23·(C_UT + C_WT + C_RM) + 減価償却費
  本実装では「減価償却費」を CAPEX/DEPRECIATION_YEARS として別途計算するため、
  Hasebe 式 (10) のうち減価償却費以外の項を OPEX 集計に反映する。

収録装置:
  - 縦型プロセス容器 (Process vessels, Vertical)
  - 熱交換器・固定管板式 (Shell & tube, Fixed tube sheet)
  - 遠心式圧縮機 (Centrifugal compressor, CS)
"""

# ---------------------------------------------------------------------------
# 基本コスト算出用係数（K1, K2, K3）
# 対象装置: 縦型プロセス容器 (Process vessels, Vertical)
# 出典: 授業資料 プロセス設計R08-3.pdf 付録A Table A.1
# ---------------------------------------------------------------------------

K1: float = 3.4974
K2: float = 0.4485
K3: float = 0.1074

# 推算式の適用可能範囲 [m³]
A_RANGE_MIN: float =   0.3
A_RANGE_MAX: float = 520.0

# ---------------------------------------------------------------------------
# モジュール係数算出用定数（B1, B2）
# 対象: 縦型プロセス容器 (Process vessels, Vertical including towers)
# 出典: プロセス設計R08-3.pdf 付録B Table A.4
# 注意: 横型容器 (Horizontal) は B1=1.49, B2=1.52 であり別値。混同しないこと。
# ---------------------------------------------------------------------------

B1: float = 2.25
B2: float = 1.82

# ---------------------------------------------------------------------------
# 材質係数（FM）: 炭素鋼（Carbon Steel）を前提
# ---------------------------------------------------------------------------

FM: float = 1.0

# ---------------------------------------------------------------------------
# コストインデックス（CEPCI: Chemical Engineering Plant Cost Index）
# ---------------------------------------------------------------------------

CEPCI_BASE:    float = 397.0    # 2001年基準値 (R08-3.pdf 掲載値)
CEPCI_CURRENT: float = 800.0    # !仮置き 推定値 2026年
                                # 推定根拠: 既存値 397→544 (2001→2016) の trend 外挿。
                                #   線形外挿 642、複利 670、近年インフレ加速考慮で 800 採用。
                                # ⚠️  Chemical Engineering 誌 PCI archive 等で実値を確認後、
                                #    確定値で置換要 (CEPCI 2026 公表値)。
                                # 影響: COOLING_WATER, AIR_COOLING, PROPYLENE/ETHYLENE_REFRIG_*,
                                #   LP/MP/HP_STEAM 等が escalation factor (CEPCI/544) で連動。

# スイング操作ペナルティ係数: K_SWING (標準 1.5)
# Bare Module Cost法の F_BM は連続定常操作を前提としており、スイング操作固有のコストを含まない。
# 本プロセスは600℃高温下で数十分ごとに流路を切り替えるため、高温用自動切替バルブ群・
# 配管マニホールド・安全インターロック機構が追加で必要となる。
# Catofin型は多基(数十基)スイングで弁数・マニホールド・計装が容器本体より支配的になりうるため、
# 旧 1.2 は楽観的。標準を 1.5 とし、感度解析で 1.2(楽観)/1.5(標準)/2.0(保守) を振る。
# env PDH_K_SWING で上書き可。
import os as _os
K_SWING: float = float(_os.environ.get('PDH_K_SWING', '1.5'))

# ---------------------------------------------------------------------------
# 経済パラメータ
# ---------------------------------------------------------------------------

USD_TO_JPY: float = 158.8595        # 為替レート [JPY/USD] 2026-05-18 06:35:00 UTC
                                    # 出典: Google Finance USD-JPY
                                    # https://www.google.com/finance/beta/quote/USD-JPY

# 減価償却年数: 8年
# 根拠: 国税庁「減価償却資産の耐用年数等に関する省令」別表第二
#   大分類「8 化学工業用設備」の「その他の設備」（法定耐用年数 8年）を適用。
#   臭素・塩化りん・活性炭等の特定設備には非該当。
#   平成20年改正以降の総合償却方式に基づき、反応器・圧縮機・蒸留塔・膜分離器等の
#   個別機器ではなくプラント一式を「化学工業用設備」として一律 8年で償却する。
DEPRECIATION_YEARS: int = 8

PLANT_INDIRECT_FACTOR: float = 1.18  # 据付・間接費拡張係数


# ---------------------------------------------------------------------------
# 熱交換器（固定管板式 / Shell & tube, Fixed tube sheet）
#
# 【機器タイプ選択理由】
#   プロピレン/プロパン系の操作温度範囲（−50〜+200°C 程度）では
#   熱膨張が管理可能であり固定管板式で問題ない。
#   遊動頭式より構造が単純でコストも低い。
#
# 【K1, K2, K3】
#   出典: プロセス設計R08-3.pdf p.9「付録A Table A.1」
#   Equipment Type    : Heat exchangers / Fixed tube sheet
#   サイジングパラメータ: 伝熱面積 A [m²]
#   K1 = 4.3247 / K2 = -0.3030 / K3 = 0.1634
#
# 【B1, B2】
#   出典: プロセス設計R08-3.pdf p.10「付録B Table A.4」
#   Equipment Type       : Heat exchangers
#   Equipment Description: Fixed tube sheet, floating head, U-tube,
#                          bayonet, kettle reboiler, and Teflon tube
#   B1 = 1.63 / B2 = 1.66
# ---------------------------------------------------------------------------

K1_HE: float =  4.3247
K2_HE: float = -0.3030
K3_HE: float =  0.1634
A_HE_MIN: float =  10.0  # 適用範囲下限 [m²]
A_HE_MAX: float = 1000.0  # 適用範囲上限 [m²]

B1_HE: float = 1.63
B2_HE: float = 1.66

# FM = 1.0: 炭素鋼/炭素鋼（C3炭化水素サービスの初期設計値）
FM_HE: float = 1.0

# Fp = 1.0: プロセス設計R08-3.pdf に HE 向け Fp 式の記載なし。
# 膜分離システムの操作圧力（< 25 bar abs）は中低圧領域であり、
# 初期設計として Fp = 1.0 を仮定する（やや保守的な見積もり）。
FP_HE_DEFAULT: float = 1.0


# ---------------------------------------------------------------------------
# 圧縮機（遠心式 / Centrifugal compressor, 炭素鋼）
#
# 【K1, K2, K3】
#   出典: プロセス設計R08-3.pdf p.9「付録A Table A.1」
#   Equipment Type    : Centrifugal, axial, and reciprocating compressors
#   サイジングパラメータ: 流体動力（Fluid power） W [kW]
#   K1 = 2.2897 / K2 = 1.3604 / K3 = -0.1027
#
# 【FBM = 2.15】
#   出典: プロセス設計R08-3.pdf p.13–14「Table A.6」および
#         「Bare Module Factor, FBM グラフ」
#   Equipment Type: Centrifugal compressor or blower
#   Material      : CS → ID番号 1 → グラフ横軸 ID=1 より FBM = 2.15 を読み取り
#
# 【Fp = 1.0（圧力補正なし）】
#   出典: プロセス設計R08-3.pdf p.11「付録C Table A.2」
#   Equipment Type: Compressors → C1 = 0, C2 = 0, C3 = 0 → Fp = 1.0
# ---------------------------------------------------------------------------

K1_COMP: float =  2.2897
K2_COMP: float =  1.3604
K3_COMP: float = -0.1027
W_COMP_MIN: float =  450.0  # 適用範囲下限 [kW]
W_COMP_MAX: float = 3000.0  # 適用範囲上限 [kW]
# 注意: 膜分離システムのテストケース（100 kmol/h フィード, P_H=10 bar, P_dist=20 bar）では
# W_feed≈182 kW, W_prod≈108 kW となり適用範囲を下回る。範囲外の外挿となるため
# コスト推算精度が低下する。最終設計流量が確定した時点で再確認すること。

FBM_COMP: float = 2.15


# ---------------------------------------------------------------------------
# ポンプ（遠心式 / Centrifugal pump, 炭素鋼）
#
# 【K1, K2, K3】
#   出典: Turton et al., "Analysis, Synthesis, and Design of Chemical Processes",
#         4th ed. (2013), Appendix A, Table A.1.
#   Equipment Type    : Centrifugal pump
#   サイジングパラメータ: 流体動力 W [kW]
#   K1 = 3.3892 / K2 = 0.0536 / K3 = 0.1538 (適用範囲 W = 1〜300 kW)
#
# 【B1, B2】
#   出典: Turton 4th ed., Table A.4. Centrifugal pump.
#   B1 = 1.89 / B2 = 1.35
#
# 【Fp 計算式】
#   出典: Turton 4th ed., Table A.2. Centrifugal pump.
#   log10(Fp) = C1 + C2*log10(P_g) + C3*(log10(P_g))^2  (P_g [barg], 10〜100 barg)
#   C1 = -0.3935, C2 = 0.3957, C3 = -0.00226
#   P_g < 10 barg では Fp = 1.0 (補正なし)。
#
# !仮置き — 要出典: 授業資料 R08-3.pdf に Pump entry の記載がない場合は
#                 Turton 4th ed. の値を流用。確認後にコメントを更新する。
# ---------------------------------------------------------------------------

K1_PUMP: float =  3.3892
K2_PUMP: float =  0.0536
K3_PUMP: float =  0.1538
W_PUMP_MIN: float =   1.0   # 適用範囲下限 [kW]
W_PUMP_MAX: float = 300.0   # 適用範囲上限 [kW]

B1_PUMP: float = 1.89
B2_PUMP: float = 1.35
FM_PUMP: float = 1.0   # 炭素鋼 (LPG サービスは構造材的に CS で十分)

C1_PUMP_FP: float = -0.3935
C2_PUMP_FP: float =  0.3957
C3_PUMP_FP: float = -0.00226


# ---------------------------------------------------------------------------
# 蒸留塔トレイ (Sieve trays / 多孔板トレイ、炭素鋼)
#
# 【機器タイプ選択理由】
#   PDH の蒸留塔は中低圧 (8.5〜20 bar) かつ C3 系炭化水素サービスで
#   流量変動が小さい。Sieve tray は最も一般的・低コストで、本プロジェクトの
#   設計レンジで標準。Valve tray (turn-down あり) は K1=3.3322 とコスト約 1.5 倍。
#
# 【K1, K2, K3 — Sieve trays】
#   出典: プロセス設計授業資料 (R08-3.pdf) 付録A Table A.1
#         (ユーザー提供の式と整合する Turton 系の標準値)
#   Equipment Type     : Sieve trays
#   サイジングパラメータ: トレイ面積 A [m²] (= 塔断面積)
#   K1 = 2.9949 / K2 = 0.4465 / K3 = 0.3961
#   適用面積範囲      : 0.07 ～ 12.30 m²
#
# 【トレイ段数係数 Fq の式】
#   出典: 同上 (R08-3.pdf 付録)
#   N < 20: log10(Fq) = 0.4771 + 0.08516*log10(N) - 0.3473*(log10(N))^2
#   N ≥ 20: Fq = 1
#   小段数では「規模の不経済」で割高になる補正項。
#
# 【FBM_TRAY】
#   出典: Turton 系標準値、CS/Sieve trays
#   FBM_TRAY = 1.0 (炭素鋼 + 圧力補正なし: トレイは塔本体ほど圧力影響を受けない)
#   材質変更時 (例 SS 等) は要見直し。
#
# 【ベアモジュールコスト】
#   C_BM,trays = Cp0_tray × N × FBM_TRAY × Fq
#   塔本体の Fp は適用しない (上記出典明記)。
# ---------------------------------------------------------------------------

K1_TRAY_SIEVE: float = 2.9949
K2_TRAY_SIEVE: float = 0.4465
K3_TRAY_SIEVE: float = 0.3961
A_TRAY_MIN:    float = 0.07     # 適用面積下限 [m²]
A_TRAY_MAX:    float = 12.30    # 適用面積上限 [m²]
FBM_TRAY:      float = 1.0      # CS / Sieve tray

# トレイ段数係数 Fq の係数 (N < 20 のとき)
#   log10(Fq) = FQ_C0 + FQ_C1*log10(N) + FQ_C2*(log10(N))^2
FQ_C0: float =  0.4771
FQ_C1: float =  0.08516
FQ_C2: float = -0.3473


# ---------------------------------------------------------------------------
# 膜モジュール単価 [USD/m²]
#
# !仮置き — 要出典: 膜モジュール単価
#   - Hua et al. (2024) または ZIF-8 膜 TEA 論文で実際の単価を確認
#   - 参考: 高分子膜の文献値 $50〜100/m² (Baker & Lokhandwala 2008)
#           ZIF-8 系 MMM は上記より高コストになる可能性あり
#   - CEPCI 基準年も要確認 (現状 CEPCI_BASE=2001年 と仮定)
# ---------------------------------------------------------------------------
MEM_UNIT_PRICE_USD_PER_M2: float = 50.0  # [USD/m²]

# 膜モジュール耐用年数 [年]
# !仮置き — 要メーカーカタログ確認。ZIF-8 系混合マトリックス膜の実機交換サイクルは
#   未確定。高分子膜の一般的なレンジ 1〜5 年の中央付近として 3.0 年を仮置きする。
# 用途: economics.py で膜の定期交換費 OPEX = CAPEX_mem / MEM_LIFETIME_YEARS を計上
#   (レポート §7.6 の C_mem/τ_mem に対応。PSA 活性炭・触媒の交換費と同じ別建て扱い)。
# 注: 本値は感度解析対象。exp/exp_membrane_sensitivity.py で 1〜5 年を振って TAC 影響を見る。
MEM_LIFETIME_YEARS: float = 3.0  # [年] !仮置き


# ==============================================================================
# PSA システム固有パラメータ
# ==============================================================================

# ---------------------------------------------------------------------------
# 活性炭に対する Langmuir 吸着定数 (293-298 K 付近)
#
# 出典: Choi, B.-U. et al., "Adsorption equilibria of methane, ethane,
#       ethylene, nitrogen, and hydrogen onto activated carbon",
#       J. Chem. Eng. Data, 48, 603-607 (2003).
#
# データ構造:
#   q_s : 飽和吸着量 [mol/kg-adsorbent]
#   a   : Langmuir 親和定数 [m³/mol]
#         文献の b [kPa⁻¹] から a = b × R × T (T=298K) で換算。
#         変換式の導出: Langmuir 等温線の圧力形 q = q_s·b·P/(1+b·P) と
#         濃度形 q = q_s·a·C/(1+a·C) を理想気体 C=P/(RT) で接続すると a=b·RT。
# ---------------------------------------------------------------------------
PSA_LANGMUIR_PARAMS: dict = {
    "CH4":  {"q_s": 5.23,  "a": 0.00332},  # [mol/kg], [m³/mol]
    "C2H4": {"q_s": 6.20,  "a": 0.0327},   # [mol/kg], [m³/mol]
    "C2H6": {"q_s": 5.66,  "a": 0.0533},   # [mol/kg], [m³/mol]
}

# ---------------------------------------------------------------------------
# 総括物質移動容量係数 K_Fa [1/s]
#
# 【設計パラメータ設定根拠】
#
# ■ CH4 = 0.02 s⁻¹ (文献実測値ベース)
#   根拠: T. E. Rufford et al., "Adsorption Equilibria and Kinetics of
#         Methane + Nitrogen Mixtures on the Activated Carbon Norit RB3",
#         Ind. Eng. Chem. Res., 2013.
#   データ: 302.3 K における Norit RB3 活性炭への CH4 の LDF 係数 (MTC) の
#           実測範囲 0.004〜0.052 s⁻¹ の代表中央値として 0.02 を採用。
#
# ■ C2H6 = 0.011 s⁻¹、C2H4 = 0.015 s⁻¹ (エンジニアリング推算値)
#   参照: X. Hu and D. D. Do, "Multicomponent adsorption kinetics of
#         hydrocarbons onto activated carbon", Chem. Eng. Sci., 1993.
#   計算: C2H6 の表面拡散係数 D_μ0 = 4.86×10⁻⁹ m²/s (Table 1, Ajax炭)
#         Glueckauf 式 K_Fa = 15·D_μ0/Rp² (Rp = 1 mm) → 0.073 s⁻¹
#   調整理由: 上記計算値 0.073 s⁻¹ は別炭種（Ajax炭）の値のため、CH4 の
#         実測値 0.02 s⁻¹ (Norit RB3) をアンカーとして整合させた。
#         一般に動的直径が大きい分子ほど細孔内拡散抵抗が大きいため、
#         CH4(3.8 Å) > C2H4(4.2 Å) > C2H6(4.4 Å) の順で減少させた。
#         C2H4・C2H6 は文献直接値ではなくエンジニアリング判断による推算値。
#
# !仮置き — 要出典: 実機データで K_Fa (特に C2H4, C2H6) を直接測定して更新
# ---------------------------------------------------------------------------
PSA_KFA: dict = {
    "CH4":  0.02,   # [1/s] Rufford et al. (2013) 実測値ベース
    "C2H4": 0.015,  # [1/s] CH4 をアンカーとしたエンジニアリング推算値
    "C2H6": 0.011,  # [1/s] 同上
}

# !仮置き — 要出典: 工業用活性炭単価
#   試薬メーカーカタログ or プラント設計文献から工業グレード単価を調査
# ---------------------------------------------------------------------------
ACTIVATED_CARBON_PRICE_USD_PER_KG: float = 5.0  # [USD/kg]

# !仮置き — 要出典: PSA 用活性炭の交換サイクル
#   典型 2〜5 年 (Ruthven, "Principles of Adsorption Processes", 1984)
#   詳細設計段階でメーカーデータ・運転条件 (圧縮強度・コーキング・不純物吸着)
#   に基づいて更新すること。中央値 4 年を仮置き。
# ---------------------------------------------------------------------------
ADSORBENT_LIFETIME_YEARS: float = 4.0  # [年]


# ==============================================================================
# ユーティリティ単価 / OPEX パラメータ
# !仮置き — 要出典: 全てコンテスト課題 Ver.2.0 のサイト仕様に置換予定
#
# 使用箇所: exp1.py のリサイクル付き TAC 計算 (flowsheet/economics.py 経由)
#
# 単位の取り決め:
#   電力     : 円/kWh
#   蒸気・燃料・冷水/冷媒 : 円/GJ  (kW 連続消費を年換算: kW × 3.6e-3 × OP_HOURS = GJ/年)
#   触媒交換  : 単価 [円/kg] ÷ 寿命 [年] × 全量 [kg] / 1e8 → 億円/年
# ==============================================================================

# 電力単価 [円/kWh]  大規模 PDH プラント = 特別高圧契約想定
# 出典: 新電力ネット「法人・家庭の電気料金の平均単価の推移」2026年1月
#   特別高圧 16.72 円/kWh (燃調費・再エネ賦課金込み)
#   URL: https://pps-net.org/unit (accessed 2026-05-18)
# 参照: 経産省 電力調査統計 https://www.enecho.meti.go.jp/statistics/electric_power/ep002/
# 参考: 再エネ賦課金 FY2026 = 4.18 円/kWh (https://enemanex.jp/saienefukakin/)
ELECTRICITY_JPY_PER_KWH: float = 17.0

# ---------------------------------------------------------------------------
# スチーム階層
# 温度の出典: PDH 学生コンテスト要綱 Ver.2.0 (ア)(イ)(ウ) より:
#   - HP Steam: 230°C 飽和蒸気で供給、戻り 230°C 飽和液 or 3bar 液
#   - MP Steam: 186°C 飽和蒸気で供給、戻り 186°C 飽和液 or 3bar 液
#   - LP Steam: 160°C 飽和蒸気で供給、戻り 160°C 飽和液 or 3bar 液
#   コンデンセートの顕熱利用も可 (本実装では潜熱のみ計上で保守)。
# 加熱箇所の温度に応じて階層から自動選択する (utility_selector.py)。
# ---------------------------------------------------------------------------
# スチーム単価 [円/GJ]
# ⚠️ !仮置き (Turton 書籍未入手 + 日本補正定数の根拠も簡易) — 入手後/精緻化後に更新要
# 推定根拠: Turton et al. "Analysis, Synthesis, and Design of Chemical Processes"
#   5th ed. Appendix Table 8.3 (CEPCI 2016 基準) の推定転載値:
#   LP (5 barg、160°C): 4.5 USD/GJ
#   MP (10 barg、184°C): 4.8 USD/GJ
#   HP (41 barg、254°C): 5.7 USD/GJ
#   Google Books ID: kWXyhVXztZ8C (Pearson Education, ISBN 9780132459181)
#
# 補正フロー:
#   円/GJ = Turton USD/GJ × CEPCI escalation (1.471) × USD/JPY (158.8595) × 日本補正 (2.0)
#                          ─────────────────         ───────────────────   ──────────────
#                          一般インフレ (CEPCI)       為替                  US 安価 NG ⇒ 日本高 LNG 差
#
# 日本補正定数 (JAPAN_STEAM_FUEL_CORRECTION)
#   日本/US 燃料差の確たる実勢値が無い段階では Turton US 基準値そのまま (補正なし、1.0)
#   を採用する。Turton 基準値自体が書籍未入手の推定なので !仮置き マーカーは継続。
#
# 計算 (補正 1.0):
#   LP: 4.5 × 1.471 × 158.8595 × 1.0 ≈ 1051 → 1050
#   MP: 4.8 × 1.471 × 158.8595 × 1.0 ≈ 1122 → 1120
#   HP: 5.7 × 1.471 × 158.8595 × 1.0 ≈ 1332 → 1330
# 温度の出典: PDH 学生コンテスト要綱 Ver.2.0 (ア)(イ)(ウ)
JAPAN_STEAM_FUEL_CORRECTION: float = 1.0   # !仮置き 日本補正なし (Turton US 基準値そのまま)
LP_STEAM_JPY_PER_GJ: float = 1050.0   # 160°C 飽和 (!仮置き Turton×CEPCI、日本補正なし)
MP_STEAM_JPY_PER_GJ: float = 1120.0   # 186°C 飽和 (!仮置き Turton×CEPCI、日本補正なし)
HP_STEAM_JPY_PER_GJ: float = 1330.0   # 230°C 飽和 (!仮置き Turton×CEPCI、日本補正なし)

# ---------------------------------------------------------------------------
# 冷媒階層 (contest.md §2-3 のユーティリティ仕様)
#   - 冷却水・空冷を超えて低温が必要な箇所では冷凍冷媒を使う必要がある。
#   - 温度が低いほど製造動力が大きく、単価が指数的に増加。
#   - エチレン冷媒の製造にプロピレン冷媒を使うため、同一温度なら
#     エチレン冷媒の方がプロピレン冷媒より高コスト (本来は重複しない温度域)。
# ---------------------------------------------------------------------------
# 冷却水・空冷単価 [円/GJ]
# ⚠️  !仮置き (Turton 書籍未入手のため数値推定) — 入手後に厳密値で更新要
# 推定根拠: Turton et al. (2018) "Analysis, Synthesis, and Design of Chemical Processes"
#   5th ed. Appendix Table 8.3 (CEPCI 2016 基準) の推定転載値
#   Google Books ID: kWXyhVXztZ8C (Pearson Education, ISBN 9780132459181)
# Escalation factor: CEPCI_CURRENT / CEPCI_2016 = 800 / 544 ≈ 1.471 (CEPCI_CURRENT も推定値)
# 換算: USD/GJ × 1.471 × USD/JPY (158.8595) = 円/GJ
#   冷却水 (towers, ΔT=10°C) ≈ 0.354 USD/GJ × 1.471 × 158.8595 ≈ 83 → 85
#   空冷 (fan power) ≈ 0.40 USD/GJ × 1.471 × 158.8595 ≈ 94 → 95
AIR_COOLING_JPY_PER_GJ:           float =    95.0   # 空冷 (!仮置き Turton 推定値、書籍未入手)
COOLING_WATER_JPY_PER_GJ:         float =    85.0   # 30→40°C (!仮置き Turton 推定値、書籍未入手)

# 冷凍冷媒単価 [円/GJ] — contest §2-3 準拠「製造所要動力」ベース
#   contest §2-3 原文要旨: 冷凍冷媒は隣接エチレンプラントから供給を受けられる (= 購入費ゼロ、
#     飽和液で受けプロセス冷却後に飽和蒸気で返す)。コストは「製造するための所要動力」で、温度が
#     低いほど大きく、エチレン冷媒はプロピレン冷媒を用いたカスケードのためプロピレンより大きい。
#     冷媒の評価方法は設計者が考察して決める旨も明記。→ 動力(電力)ベースが要綱に忠実。
#
# モデル: 冷媒製造 = 冷凍サイクルの圧縮機電力。COP = η × Carnot = η × T_c / (T_h - T_c)。
#   単価[円/GJ_cold] = 電力[円/GJ] / COP。 電力 = ELECTRICITY_JPY_PER_KWH / 3.6e-3 (=4722円/GJ)。
#   プロピレン冷媒 (+15〜-40℃): 単段、冷却水へ排熱。
#   エチレン冷媒 (-60〜-100℃): 2段カスケード — エチレン段が T_c→プロピレン-40℃ で吸熱、
#     プロピレン段が (Q + W_eth) を -40℃→冷却水 で排熱。→ 同温でもエチレン>プロピレン動力 (要綱通り)。
#   η, T_h は !仮置き の設計仮定。方針「現実にあり得る範囲でぎりぎり安く」→ η=0.60 (良好な大型
#   冷凍機の対Carnot効率の上端。これ以上は非現実的)。要調整。
REFRIG_CARNOT_EFF: float = 0.60     # !仮置き 冷凍機の対Carnot効率 (現実的 best 端、要調整)
_T_H_REJECT_K:     float = 308.15   # !仮置き 排熱先 (冷却水30℃供給 + 接近余裕 ≒ 35℃)
_T_CASCADE_MID_K:  float = 233.15   # エチレン冷媒カスケードの中間段 = プロピレン-40℃ (contest 記述)
_ELEC_JPY_PER_GJ:  float = ELECTRICITY_JPY_PER_KWH / 3.6e-3   # 17円/kWh → 4722円/GJ


def _refrig_single_jpy_per_gj(T_c_K: float) -> float:
    """単段冷凍 (T_c で吸熱 → 冷却水へ排熱) の所要動力ベース単価 [円/GJ_cold]。"""
    cop = REFRIG_CARNOT_EFF * T_c_K / (_T_H_REJECT_K - T_c_K)
    return _ELEC_JPY_PER_GJ / cop


def _refrig_cascade_jpy_per_gj(T_c_K: float) -> float:
    """2段カスケード (エチレン段 T_c→中間, プロピレン段 中間→冷却水) の単価 [円/GJ_cold]。

    エチレン冷媒製造がプロピレンより大きくなる contest §2-3 記述を物理で反映。
    W は Q_cold 基準の比仕事 (per unit Q_cold)。
    """
    w_eth  = (_T_CASCADE_MID_K - T_c_K) / (REFRIG_CARNOT_EFF * T_c_K)
    w_prop = (1.0 + w_eth) * (_T_H_REJECT_K - _T_CASCADE_MID_K) / (REFRIG_CARNOT_EFF * _T_CASCADE_MID_K)
    return _ELEC_JPY_PER_GJ * (w_eth + w_prop)


PROPYLENE_REFRIG_15C_JPY_PER_GJ:  float = _refrig_single_jpy_per_gj(288.15)   # +15°C
PROPYLENE_REFRIG_5C_JPY_PER_GJ:   float = _refrig_single_jpy_per_gj(278.15)   # +5°C
PROPYLENE_REFRIG_M25C_JPY_PER_GJ: float = _refrig_single_jpy_per_gj(248.15)   # -25°C
PROPYLENE_REFRIG_M40C_JPY_PER_GJ: float = _refrig_single_jpy_per_gj(233.15)   # -40°C
ETHYLENE_REFRIG_M60C_JPY_PER_GJ:  float = _refrig_cascade_jpy_per_gj(213.15)  # -60°C
ETHYLENE_REFRIG_M75C_JPY_PER_GJ:  float = _refrig_cascade_jpy_per_gj(198.15)  # -75°C
ETHYLENE_REFRIG_M100C_JPY_PER_GJ: float = _refrig_cascade_jpy_per_gj(173.15)  # -100°C

# 燃料 (反応器プリヒーター LNG) 単価 [円/GJ]
# 出典: JOGMEC 月例 LNG 市況調査 + 財務省 貿易統計 LNG 輸入 CIF 価格
#   URL: https://oilgas-info.jogmec.go.jp/info_reports/1009848/index.html
#   URL: https://www.customs.go.jp/toukei/suii/html/time_latest.htm
# 値: 2026-05-18 baseline、約 100,000 円/トン (1 トン = 約 54.6 GJ) → 1,830 円/GJ
# 換算: $12.0/MMBtu × 158.8595 / 1.055 GJ/MMBtu ≈ 1,807 円/GJ で近似一致
FUEL_JPY_PER_GJ: float = 1830.0

# 加熱炉熱効率 [-] (LHV ベース、燃焼計算 熱勘定の損失率より)
# Ref: 化工便覧 改訂六版 18·4·3 項 表 18·11 (れき青炭燃焼計算で全熱損失 15% = 効率 85%)
# 反応器プリヒーター燃料消費は Q_preheat / η_furnace で計算する。
FURNACE_EFFICIENCY: float = 0.85

# 触媒モデル: Cr2O3-Al2O3 (Catofin プロセス相当) をモデルとする。
# コンテスト要項 §3-3 で提供される a (失活係数) データは
#   架空触媒のものだが、その挙動 (700°C で 30 min で a≈0.04 と分単位で急速失活、
#   400°C で a≈0.81) は工業 Cr2O3-Al2O3 触媒 (Catofin プロセス) の物理と整合する。
#   Pt-Sn/Al2O3 (UOP Oleflex) は本来 CCR (連続再生) 方式で日〜週オーダーの緩慢
#   失活であり、本実装のスイング (t_regen=30 min) + a データの分単位失活とは
#   物理整合しない。よって経済パラメータも Cr2O3 系の文献値を採用する。
# 環境懸念: Cr2O3 触媒は六価 Cr の生成/排出規制が厳しい (大気汚染防止法 特定物質、
#   水質汚濁防止法 有害物質)。実プラント設計では Pt 系 (Pt-Sn/Oleflex, Pt-Ga/STARplus)
#   への置換を検討すべきだが、本シミュレータでは a データとの物理整合性を優先。
#
# 触媒単価
# 出典: Mangalindan, J. R. et al. (2025) "Tandem Cu/ZnO/ZrO2-SAPO-34 System
#   for Dimethyl Ether Synthesis from CO2 and H2: Catalyst Optimization,
#   Techno-Economic, and Carbon-Footprint Analyses". ACS Engineering Au.
#   Table 2 (TEA Parameters) にて、白金を含まない遷移金属酸化物ベース触媒の
#   単価として $23/kg を採用。CrOx/Al2O3 直接の単価表は機密情報のため、
#   検証可能なベースメタル系酸化物触媒リファレンスとして引用。
# 参照 PDF: report_for_processdesign/references/Mangalindan_2025_DME_TEA_ACS_Eng_Au.pdf
# 換算: $23/kg × USD_TO_JPY (158.8595 円/USD) ≒ 3,654 円/kg
CATALYST_USD_PER_KG: float = 23.0
CATALYST_JPY_PER_KG: float = CATALYST_USD_PER_KG * USD_TO_JPY

# 触媒寿命 [年]
# 出典: Ni, L. (2022). "Propane dehydrogenation on highly active and selective
#   Ga/BEA and ethanol conversion to butadiene on zincosilicate BEA"
#   (Doctoral dissertation, Technical University of Munich, Fakultät für Chemie;
#   accepted 2022-02-24; mediaTUM repository, doc id 1638095).
#   PDF p.20 (Chapter 1, 工業 PDH プロセス解説節) で Catofin について
#   "stable dehydrogenation performance (2-3 years lifetime)" と明記。
# 参照 PDF: report_for_processdesign/references/Ni_2022_TUM_PDH_thesis.pdf
# 採用値: レンジ 2-3 年の中央値 2.5 年。
CATALYST_LIFE_YEARS: float = 2.5

# 年間稼働時間 [h/年]  333 日連続稼働相当
# 出典: PDH 学生コンテスト要綱 (Ver.2.0) で指定された値
OPERATING_HOURS_PER_YEAR: float = 8000.0


# ==============================================================================
# ペナルティ sentinel — solver 失敗装置の CAPEX 値判定閾値
#
# 装置 CAPEX が物理計算で求まらないとき、各 unit は
# _PENALTY_CAPEX = 1e9 億円 (= 1 京円) を sentinel として返す (units/reactors/
# swing.py, units/separators/psa/psa_system.py, units/separators/membrane/
# membrane_system.py)。下流側 (flowsheet/economics.py の集計、flowsheet/
# solver.py のペナルティ判定) はこの sentinel を識別するため閾値で篩い分ける。
#
# 正常 CAPEX の最大は数十億円程度のため、1e8 億円を境界に置けば誤判定なし。
# 本定数で閾値を一元化する。
# ==============================================================================
PENALTY_CAPEX_THRESHOLD_OKUYEN: float = 1.0e8


# ==============================================================================
# Hasebe 式 (9)(10) — 労務費 / 製造コスト経験式 用パラメータ
#
# 出典: プロセス設計R08-3.pdf §3.3, §3.4
#
# Hasebe 式 (10) 製造コスト = 0.180·C_TM + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM)
#                                                                + 減価償却費
# Hasebe 式 (9)  1 班の人数 = sqrt(6.29 + 0.21·N_eq)
#                (N_eq = 主要機器数: コンプレッサー、塔、反応器、加熱器、熱交換器など)
# 4 直 3 交替: 総運転員 = 1 班の人数 × 4
# ==============================================================================

# Hasebe 式 (9) の N_eq 係数。
# 注: Turton 4th ed. では同じ式の係数を 0.23 と記載しているが、本実装は
#     一次出典である長谷部資料 (§3.3 式 (9)) に従い 0.21 を採用する。
HASEBE_NOL_COEFF: float = 0.21

# 4 直 3 交替体制の交替倍率。Hasebe §3.3 本文の例題に従い ×4 (欠勤考慮なし)。
# 注: Turton 標準は ×4.5 (欠勤・休暇を考慮した割増)。長谷部資料準拠で ×4 を採用。
HASEBE_SHIFT_MULTIPLIER: int = 4

# !仮置き — 要出典: 化学プラント運転員の年俸 [円/人/年]
#   厚労省賃金構造基本統計調査の化学工業 大企業中堅水準 (40 代男性、賞与込み)
#   600〜800 万円レンジの下〜中央値として 600 万円を仮置き。
#   コンテストの site fixed cost 仕様や直接的な業界統計で更新すること。
OPERATOR_ANNUAL_SALARY_JPY: float = 6_000_000.0

# Hasebe 式 (10) の固定係数 (上記引用ブロック参照)。可読性のためここで定数化。
HASEBE_COEFF_C_TM:    float = 0.180   # 保全 + 諸経費 + 税保険
HASEBE_COEFF_C_OL:    float = 2.73    # 労務費 (運転員 + 監督者 + 福利厚生など)
HASEBE_COEFF_C_UT_WT_RM: float = 1.23 # 用役費・廃棄物処理費・原料費の間接費込み倍率

# 廃棄物処理費 C_WT [億円/年]。
# PDH は気相反応で水処理対象廃棄物が実質発生しない (PSA オフガスは
# 燃料クレジットで回収済み、使用済み触媒は CATALYST_JPY_PER_KG の
# 回収費込みで別建て計上) ため 0 を採用する。レポート §3.3 にも注記する。
HASEBE_C_WT_OKUYEN_PER_YEAR: float = 0.0


# ==============================================================================
# 原料・製品単価 / 副産物クレジット
# !仮置き — 要出典: コンテスト課題 Ver.2.0 のサイト仕様 / 業界統計で確定単価に置換
#
# 使用箇所: flowsheet/economics.py
#   TAC = CAPEX/年 + ユーティリティ OPEX + 原料費 - 製品売上 - オフガス燃料クレジット
#
# 単位の取り決め: 全て 円/kg。
#   年量 [kg/年]      = kmol/h × MW [kg/kmol] × OPERATING_HOURS_PER_YEAR
#   年額 [億円/年]    = 年量 × 単価 [円/kg] / 1e8
#
# 単価は全て 円/kg で統一。kmol/h ベースの流量と MW で
#   年量に変換する一貫した式 (OPEX_utilities と同様の 1e8 スケール) で扱う。
# ==============================================================================

# LPG プロパン (C3H8) 単価 [円/kg]  CIF 基準 (輸入基地隣接想定)
# 大規模 PDH プラント (年産 40 万 t C3H6) は港湾コンビナート
#   立地が一般的なため、国内配送費を含まない CIF 価格を採用。
# 出典: 財務省 貿易統計 液化プロパン 輸入 CIF 価格
#   URL: https://www.customs.go.jp/toukei/suii/html/time_latest.htm (accessed 2026-05-18)
# 換算: CP プロパン $600/t × 158.8595 円/USD / 1000 = 95.3 円/kg
# 参照値 (採用しない、内陸プラント想定の場合):
#   大阪ガス CIF+フレート 119 円/kg (https://ene.osakagas.co.jp/gas/lpg/pdf/lpgkakaku_CP.pdf)
#   ChemAnalyst CFR (Q1 2026) 〜120 円/kg
LPG_C3H8_JPY_PER_KG: float = 95.0

# LPG n-ブタン (C4H10) 単価 [円/kg]  CIF 基準 (C3H8 と同シナリオ)
# C3H8 と同じ CIF ベースで揃える。実勢でも C3H8 ± 5% 以内。
# 出典: 財務省 貿易統計 液化ブタン 輸入 CIF 価格
#   URL: https://www.customs.go.jp/toukei/suii/html/time_latest.htm (accessed 2026-05-18)
# 換算: CP n-ブタン $580/t × 158.8595 円/USD / 1000 ≒ 92.5 円/kg、丸めて 95 円/kg
# 注: C3H8 と同価で扱う (ペナルティ的差別化が必要なら個別更新)
LPG_C4H10_JPY_PER_KG: float = 95.0

# Backward compat: LPG_FEED_JPY_PER_KG を参照しているコード用 (主にプロパン基準)。
# C3H8 と C4H10 を個別単価化したため、共通参照名としてプロパン値を割当。
LPG_FEED_JPY_PER_KG: float = LPG_C3H8_JPY_PER_KG

# ポリマーグレード C3H6 (≥99.5 wt%) 出荷価格 [円/kg]
# 出典: ChemAnalyst Q1 2026 Northeast Asia propylene 0.93 USD/kg ≈ 147 円/kg
#   URL: https://www.chemanalyst.com/Pricing-data/polypropylene-10
# 参照: IMARC Q3 2025 Japan CFR Nagoya 962.67 USD/t ≈ 153 円/kg
#   URL: https://www.imarcgroup.com/propylene-pricing-report
# 採用: 両者中央値 150 円/kg
C3H6_PRODUCT_JPY_PER_KG: float = 150.0

# 化学グレード H2 (≥99.9 mol%) 出荷価格 [円/kg]
# 出典: ChemAnalyst Q4 2025 Japan industrial hydrogen (SMR-based merchant gas)
#   URL: https://www.chemanalyst.com/Pricing-data/hydrogen-1165
#   値: 2,130 USD/t = 2.13 USD/kg × 158.8595 = 338 円/kg
# 注: PDH 副生 H2 を merchant SMR グレードとして売却想定。
#     圧縮トレーラー渡し (NEDO 100-120 円/Nm³ = 1100-1300 円/kg) は買値であり該当しない。
#     売値として 400 円/kg を採用 (販売マージン込み)
H2_PRODUCT_JPY_PER_KG: float = 400.0


# ==============================================================================
# HHV (高位発熱量) — オフガス燃料クレジット計算用
# !仮置き — 出典確認要: プロジェクト既存物性 (THERMO_DATA の dHf_298 等) は
#                       化工便覧 改訂六版 もしくは thermo パッケージ初期値を一次出典と
#                       しているため、HHV も同系列で確認推奨。値は標準値で問題ないが
#                       一次出典の確定をユーザーが行うこと。
#
# PSA オフガス (主に CH4 + 残留 H2/C2 系 + C3 トレース) を反応器プリヒーター燃料として
# 戻すと、equivalent FUEL_JPY_PER_GJ 分だけ燃料費が浮く。economics.py で:
#   credit [億円/年] = Σ (offgas[i] [kmol/h] × HHV[i] / 1000)  [GJ/h]
#                      × OPERATING_HOURS_PER_YEAR × FUEL_JPY_PER_GJ / 1e8
#
# 単位: [MJ/kmol]   (= kJ/mol × 1)
# ==============================================================================
HHV_MJ_PER_KMOL: dict = {
    'A': 2220.0,  # C3H8
    'B': 2058.0,  # C3H6
    'C':  286.0,  # H2  (液水基準 HHV)
    'D': 1411.0,  # C2H4
    'E':  891.0,  # CH4
    'F': 1561.0,  # C2H6
    'Z': 2877.0,  # C4H10
}
