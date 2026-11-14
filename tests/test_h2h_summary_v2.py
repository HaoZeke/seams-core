#!/usr/bin/env python3
"""Power-law fit and the adaptive-CNA cause labels."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from h2h_summary_v2 import disagreement_cause, extrapolate_ms, fit_power


def test_square_law():
    sizes = [64, 216, 512, 1000, 4096]
    times = [size * size * 1.0e-4 for size in sizes]
    fit = fit_power(sizes, times)
    assert abs(fit["exponent"] - 2.0) < 1e-8
    assert abs(fit["r2"] - 1.0) < 1e-8
    assert abs(extrapolate_ms(fit, 8000) - 8000 * 8000 * 1.0e-4) < 1e-6


def test_causes():
    assert disagreement_cause("hexagonal_ice", 0, 4) == "four_neighbor"
    assert disagreement_cause("hydrate", 0, 5) == "structured_other_shell"
    assert disagreement_cause("other", 0, 2) == "agree"
    assert disagreement_cause("other", 1, 12) == "cna_only"
    assert disagreement_cause("cubic_ice", 1, 4) == "both_named"


if __name__ == "__main__":
    test_square_law()
    test_causes()
    print("ok")
