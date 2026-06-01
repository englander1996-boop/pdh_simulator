# -*- coding: utf-8 -*-
"""superbatch — super_main.py 用の検証バッチ部品 (モジュール化)。

構成:
  config    : 検証計画 (PLAN / 試行数 / seed) とパス・定数。**編集はここだけ**。
  util      : ログ・時間整形の小物。
  manifest  : 完了 run の記録 (resume) と best.json/params の読み出し。
  runner    : main.py を 1 本ずつ直列起動する (HYSYS 競合ガード込み)。
  aggregate : 全 run の集計 (best-TAC 分布 / champion / BO 正当性 / params 収束)。
  plots     : matplotlib による可視化 (PNG 出力)。

エントリは repo 直下の super_main.py (薄い UI)。
"""
