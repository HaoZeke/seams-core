#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Regenerate the 4096-molecule bulk TIP4P/2005 liquid.

Molecules start on a 16 x 16 x 16 simple-cubic grid at 0.997 g/cm^3
with random orientations, the bulk builder of input/brine. in.bulk
minimises with harmonic bonds, then holds the geometry with Shake,
melts the grid in NVT at 300 K for 20 ps and equilibrates in NPT at
300 K and 1 atm for 100 ps on a 2 fs step. Only the oxygens of the
last step are kept, as tip4p2005-4096.lammpstrj under TIP4P2005_DIR
or `tip4p2005-frames`. The dump stays outside git; the log records
its sha256, the LAMMPS banner, the thread count and the last thermo
rows. tip4p2005-4096.sha256 is the frame the benchmarks used; check a
regenerated one with `sha256sum -c` inside TIP4P2005_DIR. OpenMP
reductions make the trajectory, and so the hash, machine-dependent.

Usage:
  python3 input/tip4p2005/regenerate.py --log benchmarks/h2h-sqq-tip4p2005.txt
"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DECK = Path(__file__).resolve().parent
ROOT = DECK.parents[1]
sys.path.insert(0, str(ROOT / "input" / "hop"))
import regenerate_tip4p as hop  # noqa: E402

brine = hop.brine
N_WATER = 4096
NAME = f"tip4p2005-{N_WATER}.lammpstrj"


def out_dir() -> Path:
    return Path(os.environ.get("TIP4P2005_DIR", "tip4p2005-frames"))


def frame_path() -> Path:
    return out_dir() / NAME


def build(seed: int):
    rng = random.Random(seed)
    origins, box = brine.liquid_origins(N_WATER, rng)
    waters = [hop.whole_molecule(*brine.tip4p_at(o, rng), box) for o in origins]
    hop.molecules_intact(waters)
    if hop.cell_pairs([w[0] for w in waters], box, 2.2):
        raise SystemExit("two oxygens are closer than 2.2 A")
    return waters, box


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def oxygen_count(path: Path) -> int:
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                return int(next(handle))
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--seed", type=int, default=2005)
    parser.add_argument("--temp", type=float, default=300.0)
    parser.add_argument("--nmelt", type=int, default=10000)
    parser.add_argument("--nequil", type=int, default=50000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()
    if shutil.which(args.lmp) is None and not Path(args.lmp).is_file():
        raise SystemExit(f"lmp not found: {args.lmp}")

    waters, box = build(args.seed)
    out_dir().mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    with tempfile.TemporaryDirectory(prefix="tip4p2005-") as tmp:
        work = Path(tmp)
        hop.write_data(work / "data.tip4p", waters, box)
        shutil.copy(DECK / "in.bulk", work / "in.bulk")
        cmd = [args.lmp, "-sf", "omp", "-pk", "omp", str(args.threads), "-nocite",
               "-in", "in.bulk", "-log", "log.bulk", "-screen", "none",
               "-var", "seed", str(args.seed), "-var", "temp", f"{args.temp:.2f}",
               "-var", "nmelt", str(args.nmelt), "-var", "nequil", str(args.nequil)]
        subprocess.run(cmd, cwd=work, env=env, check=True)
        banner = hop.lammps_banner(work / "log.bulk")
        tail = hop.thermo_tail(work / "log.bulk")
        shutil.copy(work / "bulk.lammpstrj", frame_path())

    dest = frame_path()
    lines = [
        banner,
        f"deck input/tip4p2005/in.bulk seed {args.seed} temp {args.temp:.2f} K "
        f"nmelt {args.nmelt} nequil {args.nequil} step 2 fs threads {args.threads}",
        f"start grid {N_WATER} molecules box {box[0]:.4f} A",
    ]
    lines += [f"thermo step {s} temp {t:.2f} density {d:.4f}" for s, t, d in tail]
    lines.append(f"{dest.stem} oxygens {oxygen_count(dest)} sha256 {sha256(dest)}")
    args.log.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
