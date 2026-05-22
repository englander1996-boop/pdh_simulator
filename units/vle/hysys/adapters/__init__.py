"""各塔の HYSYS 入出力アダプタ。

アダプタの責務:
  - 塔ごとに異なる HYSYS オブジェクト名 (Topv 等)、スペック名 (Comp Fraction-2 /
    Reflux Ratio / Draw Rate)、流量単位 (kgmole/h vs kgmol/s) を吸収
  - 入力契約 (Column{1,2,3}Input) を受けて HYSYS に書込み → ソルバ → 出力回収
  - 結果は塔非依存の ColumnResult で返す
"""

from units.vle.hysys.adapters.types import (
    Column1Input, Column2Input, Column3Input, ColumnResult,
)
from units.vle.hysys.adapters.column1 import Column1Adapter
from units.vle.hysys.adapters.column2 import Column2Adapter
from units.vle.hysys.adapters.column3 import Column3Adapter
from units.vle.hysys.adapters.components import (
    PDH_TO_HYSYS_KEYWORDS, build_pdh_to_hysys_index,
)

__all__ = [
    "Column1Input", "Column2Input", "Column3Input", "ColumnResult",
    "Column1Adapter", "Column2Adapter", "Column3Adapter",
    "PDH_TO_HYSYS_KEYWORDS", "build_pdh_to_hysys_index",
]
