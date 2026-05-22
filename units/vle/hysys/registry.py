"""段数 → HSC ファイルパス解決レジストリ。

PDH シミュレーター内 (hysys_cases/column{1,2,3}/{N}.hsc) を走査し、
要求段数に対する HSC ファイルパスを返す。
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# プロジェクト直下を起点に HSC ルートを決定 (本ファイル: units/vle/hysys/registry.py)。
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_HSC_ROOT = _PROJECT_ROOT / "hysys_cases"

# 塔識別子 → サブディレクトリ名
COLUMN_SUBDIR = {
    "column1": "column1",
    "column2": "column2",
    "column3": "column3",
}


class StageNotAvailableError(FileNotFoundError):
    """要求された段数の HSC が見つからない。"""


@dataclass
class HsysRegistry:
    """段数 → HSC パスの解決。

    Parameters
    ----------
    root : Path
        HSC ルートディレクトリ。配下に column1/, column2/, column3/ を期待。
    """
    root: Path = field(default_factory=lambda: DEFAULT_HSC_ROOT)
    # 塔識別子 → {段数: パス} のキャッシュ
    _index: Dict[str, Dict[int, Path]] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        self.root = Path(self.root)
        if not self.root.exists():
            raise FileNotFoundError(
                f"HSC ルートが見つからない: {self.root}\n"
                f"先に hysys_cases/ を作成して HSC をコピーすること。"
            )
        self._reindex()

    # -------- 構築 --------
    def _reindex(self) -> None:
        """配下を走査して段数 → パスの辞書を作る。"""
        self._index.clear()
        for col_id, subdir in COLUMN_SUBDIR.items():
            d = self.root / subdir
            stage_map: Dict[int, Path] = {}
            if d.is_dir():
                for f in d.iterdir():
                    if not f.is_file() or f.suffix.lower() != ".hsc":
                        continue
                    base = f.stem
                    try:
                        n = int(base)
                    except ValueError:
                        # 数字以外のファイル名はスキップ (例: backup ファイル等)
                        continue
                    stage_map[n] = f
            self._index[col_id] = stage_map

    # -------- 公開 API --------
    def list_stages(self, column_id: str) -> List[int]:
        """その塔で利用可能な段数を昇順で返す。"""
        if column_id not in self._index:
            raise KeyError(f"未知の塔識別子: {column_id} (有効: {list(COLUMN_SUBDIR)})")
        return sorted(self._index[column_id].keys())

    def get_path(self, column_id: str, n_stages: int) -> Path:
        """指定段数の HSC パスを返す。無ければ StageNotAvailableError。"""
        if column_id not in self._index:
            raise KeyError(f"未知の塔識別子: {column_id} (有効: {list(COLUMN_SUBDIR)})")
        stage_map = self._index[column_id]
        if n_stages not in stage_map:
            avail = sorted(stage_map.keys())
            raise StageNotAvailableError(
                f"{column_id}: 段数 {n_stages} の HSC が無い。"
                f" 利用可能: {avail[0]}〜{avail[-1]} ({len(avail)}件)"
                if avail else
                f"{column_id}: 段数 {n_stages} の HSC が無い (利用可能段数なし)"
            )
        return stage_map[n_stages]

    def resolve(
        self,
        column_id: str,
        n_stages: int,
        mode: str = "exact",
    ) -> Tuple[int, Path]:
        """段数解決。

        mode='exact'   : 完全一致のみ、無ければ StageNotAvailableError
        mode='nearest' : 最寄り段に丸めて (実段数, パス) を返す
        """
        if mode == "exact":
            return n_stages, self.get_path(column_id, n_stages)
        if mode == "nearest":
            avail = self.list_stages(column_id)
            if not avail:
                raise StageNotAvailableError(
                    f"{column_id}: 利用可能段数なし"
                )
            n_best = min(avail, key=lambda x: abs(x - n_stages))
            return n_best, self._index[column_id][n_best]
        raise ValueError(f"未対応 mode: {mode} (使えるのは 'exact'|'nearest')")

    def summary(self) -> str:
        """デバッグ用の人間可読サマリ。"""
        lines = [f"HsysRegistry root={self.root}"]
        for col_id in COLUMN_SUBDIR:
            stages = self.list_stages(col_id)
            if stages:
                lines.append(
                    f"  {col_id}: {len(stages)}件, {stages[0]}〜{stages[-1]}"
                )
            else:
                lines.append(f"  {col_id}: なし")
        return "\n".join(lines)
