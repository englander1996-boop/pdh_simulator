"""PDH 側成分キーと HYSYS 内コンポーネント名のマッピング。

PDH 側キー (stream.stream.ProcessStream.F_in のキー):
  'A'=C3H8, 'B'=C3H6, 'C'=H2, 'D'=C2H4, 'E'=CH4, 'F'=C2H6, 'Z'=C4H10

HYSYS 側コンポーネント名は HSC ファイルの fluid package に依存するが、
慣例的に下記キーワードを名前に含む。マッチ判定は部分一致 (substring) で行う
(lhs_column2.py の実装と同じ)。例: "Propene" は "Propylene" も拾うため別途扱う。
"""
from __future__ import annotations

from typing import Dict, List, Optional

# PDH キー → HYSYS 成分名に含まれそうなキーワード一覧。
# 各キーワード集合のいずれかが HYSYS 成分名に含まれていればマッチとみなす。
# 順序重要: "Propane" は "Propanol" 等とは別、また "Propene" と "Propylene" は別物
# (HYSYS では Propylene 表記が一般的) として両方フォールバックを並べる。
PDH_TO_HYSYS_KEYWORDS: Dict[str, List[str]] = {
    "A": ["Propane"],
    "B": ["Propene", "Propylene"],
    "C": ["Hydrogen"],
    "D": ["Ethylene", "Ethene"],
    "E": ["Methane"],
    "F": ["Ethane"],
    "Z": ["n-Butane", "Butane"],
}


def build_pdh_to_hysys_index(comps) -> Dict[str, int]:
    """HYSYS Components コレクションを走査し、PDH キー → HYSYS index の辞書を返す。

    見つからなかったキーは戻り値の辞書に含めない (キーが落ちる)。
    呼び出し側で必要に応じて警告・エラーを出すこと。
    """
    n = comps.Count
    names = [comps.Item(j).Name for j in range(n)]
    result: Dict[str, int] = {}
    for key, kws in PDH_TO_HYSYS_KEYWORDS.items():
        for j, name in enumerate(names):
            if any(kw in name for kw in kws):
                result[key] = j
                break
    return result


def hysys_index_of(comps, pdh_key: str) -> Optional[int]:
    """単独 PDH キーの HYSYS index を返す (見つからなければ None)。"""
    kws = PDH_TO_HYSYS_KEYWORDS.get(pdh_key)
    if not kws:
        return None
    for j in range(comps.Count):
        name = comps.Item(j).Name
        if any(kw in name for kw in kws):
            return j
    return None
