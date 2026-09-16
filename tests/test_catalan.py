"""Unit tests for the Catalan's Constant calculator.

The module lives in a file with a space and apostrophe in its name, so it is
loaded via importlib rather than a normal import.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parent.parent / "Catalan's Constant.py"
_spec = importlib.util.spec_from_file_location("catalan_impl", MODULE_PATH)
catalan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(catalan)

# Published decimal expansion of Catalan's constant (OEIS A006752).
KNOWN_DIGITS = "0915965594177219015054603514932384110774149374281"


class TestSeriesMath:
    """The Zucker-Bradley series machinery."""

    def test_series_matches_direct_sum(self):
        """bs_series() must equal the term-by-term sum of 1/((2n+1)^2 * C(2n,n))."""
        import math
        from fractions import Fraction

        num_terms = 40
        T, Q = catalan.bs_series(num_terms)

        direct = Fraction(0)
        for n in range(num_terms):
            direct += Fraction(1, (2 * n + 1) ** 2 * math.comb(2 * n, n))

        assert abs(Fraction(int(T), int(Q)) - direct) < Fraction(1, 10 ** 30)

    def test_terms_needed_scales(self):
        """More digits require proportionally more terms (~1.7 per digit)."""
        for digits in (10, 50, 200):
            n = catalan.series_terms_needed(digits)
            assert 1.5 * digits <= n <= 2 * digits + 12


class TestComputeCatalan:
    """End-to-end digit generation."""

    def test_known_digits(self):
        """First 49 digits match the published expansion of Catalan's constant."""
        assert catalan.compute_catalan(len(KNOWN_DIGITS)) == KNOWN_DIGITS

    def test_prefix_property(self):
        """Digits computed at lower precision are a prefix of higher precision."""
        short = catalan.compute_catalan(15)
        long = catalan.compute_catalan(150)
        assert long.startswith(short)

    def test_requested_length_returned(self):
        for n in (1, 2, 10, 77):
            assert len(catalan.compute_catalan(n)) == n

    def test_rejects_nonpositive_digits(self):
        with pytest.raises(ValueError):
            catalan.compute_catalan(0)
        with pytest.raises(ValueError):
            catalan.compute_catalan(-5)


class TestOeisFiles:
    """OEIS output file generation."""

    def test_files_written_correctly(self, tmp_path, monkeypatch):
        """Both output files appear with the expected names, content, and b-file format."""
        monkeypatch.chdir(tmp_path)
        digits = catalan.compute_catalan(20)
        catalan.save_oeis_files("Catalan", digits, 20)

        raw = tmp_path / "Catalan_20_digits.txt"
        bfile = tmp_path / "b_file_Catalan_20.txt"

        assert raw.read_text() == digits
        lines = bfile.read_text().splitlines()
        assert len(lines) == 20
        assert lines[0] == f"1 {digits[0]}"
        assert lines[19] == f"20 {digits[19]}"

    def test_extra_digits_truncated(self, tmp_path, monkeypatch):
        """A digit string longer than target_digits is truncated."""
        monkeypatch.chdir(tmp_path)
        catalan.save_oeis_files("Test", "0123456789", 5)
        assert (tmp_path / "Test_5_digits.txt").read_text() == "01234"


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
