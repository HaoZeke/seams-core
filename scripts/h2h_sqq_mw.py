#!/usr/bin/env python3
"""215 K mW liquids for the SQQ head-to-head.

input/mw215/in.equilibrate on one OpenMP thread with velocity seed 215.
16 cells on an edge is 4096 molecules and 48 is 110592. The public
d-SEAMS/mw-liquid-215K download answers 401, so the deck is rerun with
`lmp` from PATH. Frames land under MW215_DIR, outside the tree. The log
records the LAMMPS banner and the sha256 of every dump.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "input" / "mw215"))

import regenerate  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def lammps_banner() -> str:
    proc = subprocess.run(["lmp", "-h"], capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        if "Large-scale Atomic/Molecular" in line:
            return line.strip()
    return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", nargs="*", type=int, default=[16, 48])
    parser.add_argument("--log", required=True)
    args = parser.parse_args()
    lines = [lammps_banner(), "deck input/mw215/in.equilibrate seed 215 OMP_NUM_THREADS 1"]
    for cells in args.cells:
        path = regenerate.frame_path(cells)
        if not path.is_file():
            regenerate.generate(cells)
        lines.append(f"mw-liquid-{cells ** 3} cells {cells} sha256 {sha256(path)}")
    Path(args.log).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
