"""Ideal ice frames and the shipped mW liquid, as LAMMPS dumps.

Ice Ic and Ih use the lattices of paper/tex/figures/worked_example.py
(O-O distance 2.75 angstrom). The liquid is the 512-molecule mW dump
converted from dseams2_repro/data/mW_liquid_T215_512.data. The comparator
campaign's liquid dump is not in the seams-core checkout.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT

OUT = ROOT / "benchmarks" / "cpc" / "frames"
LIQUID_DATA = ROOT.parent / "dseams2_repro" / "data" / "mW_liquid_T215_512.data"


def _worked():
    spec = importlib.util.spec_from_file_location(
        "worked_example", ROOT / "paper" / "tex" / "figures" / "worked_example.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_liquid(path: Path, source: Path) -> None:
    lines = source.read_text().splitlines()
    lo = [0.0, 0.0, 0.0]
    hi = [0.0, 0.0, 0.0]
    for line in lines:
        for axis, key in enumerate(("xlo xhi", "ylo yhi", "zlo zhi")):
            if line.strip().endswith(key):
                lo[axis], hi[axis] = (float(part) for part in line.split()[:2])
    start = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 2
    atoms = []
    for line in lines[start:]:
        if not line.strip():
            if atoms:
                break
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        atoms.append((int(parts[0]), 1, float(parts[2]), float(parts[3]), float(parts[4])))
    atoms.sort()
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n" % len(atoms))
        handle.write("ITEM: BOX BOUNDS pp pp pp\n")
        for origin, upper in zip(lo, hi):
            handle.write("%.12f %.12f\n" % (origin, upper))
        handle.write("ITEM: ATOMS id type x y z\n")
        for index, (_, typ, x, y, z) in enumerate(atoms, start=1):
            handle.write("%d %d %.10f %.10f %.10f\n" % (index, typ, x, y, z))


def main() -> None:
    worked = _worked()
    OUT.mkdir(parents=True, exist_ok=True)
    for cells in (3, 4, 5, 6, 8):
        points, box = worked.cubic_diamond(cells)
        worked.write_lammps(OUT / f"ic-{len(points)}.lammpstrj", points, box)
        print(f"ic-{len(points)} {len(points)}")
    points, box = worked.lonsdaleite(5)
    worked.write_lammps(OUT / f"ih-{len(points)}.lammpstrj", points, box)
    print(f"ih-{len(points)} {len(points)}")
    if LIQUID_DATA.is_file():
        write_liquid(OUT / "mw-liquid-512.lammpstrj", LIQUID_DATA)
        print("mw-liquid-512")
    else:
        print("mw-liquid-512 absent")


if __name__ == "__main__":
    main()
