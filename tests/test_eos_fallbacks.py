"""
src/eos.py のフォールバック動作テスト

テスト対象:
    1. Z=1 フォールバック
       z_factor 内: _cubic_z が有効実根を返せない場合に Z=1 + UserWarning を返す

    2. 理想気体フォールバック
       compress_isentropic 内: T2_actual の brentq が ValueError の場合に
       Cp 一定の理想気体近似で T2 を推算する

各フォールバックについて以下の 2 点を確認する:
    A) 通常の物理条件では発動しない (正常稼働確認)
    B) 意図して発動させたとき正しく動作する (フォールバック動作確認)
"""

import os
import sys
import warnings
import unittest
from unittest.mock import patch

import scipy.optimize as sp_opt

# プロジェクトルートを sys.path に追加してから import
_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import src.eos as eos_module
from src.eos import z_factor, compress_isentropic

# ---------- テスト共通パラメータ ----------
# C3H6 (キー 'B') / C3H8 (キー 'A')  ← membrane_system.py と同じ順序
_KEYS   = ["B", "A"]
_X_HALF = [0.5, 0.5]   # 等モル混合
_X_C3H6 = [0.9, 0.1]   # C3H6 リッチ（膜透過側の典型組成）


# ===========================================================================
# 1-A. Z=1 フォールバック — 通常条件では発動しないこと
# ===========================================================================

class TestZFactorNormal(unittest.TestCase):
    """現実的な C3H6/C3H8 膜分離条件で Z=1 フォールバックが発動しないことを確認する"""

    def _check(self, T: float, P: float, x: list) -> float:
        """UserWarning が出ないこと & Z が物理的な範囲内であることを確認"""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            Z = z_factor(T, P, x, _KEYS, "vapor")
        self.assertGreater(Z, 0.3, "Z が物理的下限を下回った")
        self.assertLess(Z, 1.5,   "Z が物理的上限を超えた")
        return Z

    def test_atmospheric(self):
        """低圧（1 bar, 320 K）でフォールバックなし"""
        self._check(320.0, 1.0e5, _X_HALF)

    def test_feed_10bar(self):
        """膜フィード典型圧（10 bar, 320 K）でフォールバックなし、かつ Z < 1"""
        Z = self._check(320.0, 10.0e5, _X_HALF)
        self.assertLess(Z, 1.0, "高圧 C3 混合は Z < 1 になるはず（引力項優位）")

    def test_feed_20bar(self):
        """膜フィード最大圧（20 bar, 330 K）でフォールバックなし"""
        self._check(330.0, 20.0e5, _X_HALF)

    def test_pure_c3h6(self):
        """純 C3H6 でフォールバックなし"""
        self._check(320.0, 10.0e5, [1.0, 0.0])

    def test_pure_c3h8(self):
        """純 C3H8 でフォールバックなし"""
        self._check(320.0, 10.0e5, [0.0, 1.0])

    def test_c3h6_rich(self):
        """C3H6 90% 混合でフォールバックなし"""
        self._check(310.0, 10.0e5, _X_C3H6)


# ===========================================================================
# 1-B. Z=1 フォールバック — 発動時の動作確認
# ===========================================================================

class TestZFactorFallbackZ1(unittest.TestCase):
    """_cubic_z が [] を返す状況を人工的に作り、z_factor の挙動を確認する"""

    def test_returns_one_with_warning(self):
        """`_cubic_z` が [] のとき Z=1.0 を返し UserWarning を発行する"""
        with patch("src.eos._cubic_z", return_value=[]):
            with self.assertWarns(UserWarning) as cm:
                Z = z_factor(320.0, 10.0e5, _X_HALF, _KEYS, "vapor")
        self.assertEqual(Z, 1.0)
        self.assertIn("Z=1", str(cm.warning))

    def test_warning_message_contains_T_and_P(self):
        """警告メッセージに温度・圧力の数値が含まれること"""
        with patch("src.eos._cubic_z", return_value=[]):
            with self.assertWarns(UserWarning) as cm:
                z_factor(350.0, 15.0e5, _X_HALF, _KEYS, "vapor")
        msg = str(cm.warning)
        self.assertIn("350.0", msg,  "温度が警告メッセージに含まれていない")
        self.assertIn("15.00", msg,  "圧力（bar）が警告メッセージに含まれていない")

    def test_cubic_z_returns_roots_in_normal_conditions(self):
        """通常条件では _cubic_z 自体が非空リストを返すことを確認"""
        from src.eos import _cubic_z, _mix
        A, B, *_ = _mix(320.0, 10.0e5, _X_HALF, _KEYS)
        roots = _cubic_z(A, B)
        self.assertGreater(len(roots), 0, "_cubic_z が通常条件で空リストを返した")
        # 最大根（気相 Z）は B より大きいはず
        self.assertGreater(max(roots), B)


# ===========================================================================
# 2-A. 理想気体フォールバック — 通常条件では発動しないこと
# ===========================================================================

class TestCompressIsentropicNormal(unittest.TestCase):
    """現実的な圧縮条件で brentq フォールバックが発動しないことを確認する"""

    def _run(self, T1: float, P1: float, P2: float,
             x: list, eta: float = 0.75):
        """UserWarning なしで T2 > T1、W > 0 を確認"""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            T2, W = compress_isentropic(T1, P1, P2, x, _KEYS, eta=eta)
        self.assertGreater(T2, T1,   "圧縮後 T2 が T1 以下になっている")
        self.assertGreater(W,  0.0,  "圧縮仕事が非正")
        self.assertLess(T2,  1000.0, "T2 が物理的上限（1000 K）を超えた")
        return T2, W

    def test_feed_compressor(self):
        """フィード圧縮（1 bar → 10 bar, 310 K）でフォールバックなし"""
        self._run(310.0, 1.0e5, 10.0e5, _X_HALF)

    def test_product_compressor(self):
        """製品圧縮（0.5 bar → 15 bar, C3H6 リッチ）でフォールバックなし"""
        self._run(310.0, 0.5e5, 15.0e5, _X_C3H6)

    def test_mild_compression(self):
        """軽度圧縮（5 bar → 10 bar, 320 K）でフォールバックなし"""
        self._run(320.0, 5.0e5, 10.0e5, _X_HALF)

    def test_high_eta(self):
        """高効率 η=0.95 でフォールバックなし"""
        self._run(310.0, 1.0e5, 10.0e5, _X_HALF, eta=0.95)


# ===========================================================================
# 2-B. 理想気体フォールバック — 発動時の動作確認
# ===========================================================================

class TestCompressIsentropicFallbackIdealGas(unittest.TestCase):
    """
    T2_actual の brentq を人工的に失敗させ、理想気体近似フォールバックを検証する。

    compress_isentropic 内の brentq 呼び出し順序:
        呼び出し 1: entropy_balance から T2s を探索 (通常条件では必ず成功)
        呼び出し 2: enthalpy_balance から T2_actual を探索
                    ← ここを ValueError にして理想気体フォールバックを発動させる
    """

    _T1, _P1, _P2 = 310.0, 1.0e5, 10.0e5
    _x = _X_HALF
    _eta = 0.75

    def _make_mock_brentq(self):
        """1 回目は本物の brentq に通し、2 回目は ValueError を発生させるモック"""
        real_brentq = sp_opt.brentq
        call_count = [0]

        def _mock(f, a, b, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                # T2s 探索: 本物の brentq に通す
                return real_brentq(f, a, b, **kwargs)
            # T2_actual 探索: 強制 ValueError → 理想気体フォールバックへ
            raise ValueError("テスト用強制失敗: T2_actual brentq")

        return _mock, call_count

    def test_fallback_returns_valid_T2(self):
        """フォールバック後も T2 > T1、W > 0 かつ上限以内であること"""
        mock_brentq, call_count = self._make_mock_brentq()

        with patch("src.eos.brentq", side_effect=mock_brentq):
            T2, W = compress_isentropic(
                self._T1, self._P1, self._P2, self._x, _KEYS, eta=self._eta
            )

        self.assertGreater(T2, self._T1, "フォールバック T2 が T1 以下")
        self.assertGreater(W,  0.0,      "フォールバック後の仕事量が非正")
        self.assertLess(T2,   1000.0,    "フォールバック T2 が 1000 K 超")
        self.assertEqual(call_count[0], 2, "brentq の呼び出し回数が想定 (2) と異なる")

    def test_W_actual_is_unchanged_by_fallback(self):
        """
        W_actual は T2s から計算済みのため、フォールバックの有無に関係なく
        正常値と一致するはず。
        """
        # 正常時の W を取得
        _, W_normal = compress_isentropic(
            self._T1, self._P1, self._P2, self._x, _KEYS, eta=self._eta
        )

        mock_brentq, _ = self._make_mock_brentq()
        with patch("src.eos.brentq", side_effect=mock_brentq):
            _, W_fallback = compress_isentropic(
                self._T1, self._P1, self._P2, self._x, _KEYS, eta=self._eta
            )

        # フォールバックは T2_actual のみに影響し、W_actual は変わらない
        self.assertAlmostEqual(W_fallback, W_normal, places=3,
                               msg="フォールバック時の W_actual が正常値と異なる")

    def test_fallback_T2_in_reasonable_range(self):
        """
        理想気体近似は粗い近似だが、正常値の ±30% 以内に収まることを確認する。
        """
        T2_normal, _ = compress_isentropic(
            self._T1, self._P1, self._P2, self._x, _KEYS, eta=self._eta
        )

        mock_brentq, _ = self._make_mock_brentq()
        with patch("src.eos.brentq", side_effect=mock_brentq):
            T2_fallback, _ = compress_isentropic(
                self._T1, self._P1, self._P2, self._x, _KEYS, eta=self._eta
            )

        delta = T2_normal * 0.30
        self.assertAlmostEqual(
            T2_fallback, T2_normal, delta=delta,
            msg=f"フォールバック T2={T2_fallback:.1f} K が正常値 {T2_normal:.1f} K から ±30% 超"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
