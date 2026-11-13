#!/usr/bin/env python3
"""Assert carried rings stay at or under a bare ringNetwork call."""

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.argv = [
    str(ROOT / "scripts" / "h2h_carry_fraction.py"),
    sys.argv[1],
    "--samples",
    "8",
]
runpy.run_path(str(ROOT / "scripts" / "h2h_carry_fraction.py"), run_name="__main__")
