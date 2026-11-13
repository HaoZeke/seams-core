#!/usr/bin/env python3
"""Orthorhombic supercell of the tilted sH dump, for SQQ.

The genice sH cell has xy = -lx/4 and no xz or yz tilt. Four images
along b tile the rectangle (lx) by (4*ly) by lz. SQQ accepts that
dump. Cage counts are compared with the seams command on the same
file. Eight samples, one thread and eight threads.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from h2h_cells import _parse_sH  # noqa: E402

N = 8
OUT = ROOT / "bbdir-release" / "h2h-sh-ortho.lammpstrj"


def ortho_images():
    pos, origin, (lx, ly, lz), (xy, xz, yz) = _parse_sH()
    if abs(xz) > 1e-8 or abs(yz) > 1e-8 or abs(xy * 4 + lx) > 1e-6:
        raise SystemExit(f"sH tilt is not xy=-lx/4: {lx} {xy} {xz} {yz}")
    pos = pos - np.asarray(origin)
    blocks = []
    for image in range(4):
        shifted = pos + np.array([image * xy, image * ly, 0.0])
        shifted[:, 0] -= lx * np.floor(shifted[:, 0] / lx)
        shifted[:, 1] -= (4.0 * ly) * np.floor(shifted[:, 1] / (4.0 * ly))
        shifted[:, 2] -= lz * np.floor(shifted[:, 2] / lz)
        blocks.append(shifted)
    stacked = np.vstack(blocks)
    # Drop a duplicate if periodic wrap landed on the same point.
    key = np.round(stacked, 5)
    _, index = np.unique(key, axis=0, return_index=True)
    stacked = stacked[np.sort(index)]
    return stacked, (lx, 4.0 * ly, lz)


def write_dump(path: Path, pos, lengths) -> None:
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(pos)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {lengths[0]:.8f}",
        f"0 {lengths[1]:.8f}",
        f"0 {lengths[2]:.8f}",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def sqq_counts(pos, lengths, width: int) -> dict:
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame

    pos_nm = np.asarray(pos, dtype=np.float64) / 10.0
    box = np.asarray(lengths, dtype=np.float64) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(len(pos_nm))
    ]
    frame = Frame(name="sH-ortho", atoms=atoms, box=box)
    cfg = load_config({
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
        "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": True},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": width, "math_threads": 1},
    }, engine="py")
    result = analyze_frame(frame, cfg, engine="py")
    counts = {}
    for cage in result.cages:
        counts[cage.cage_type] = counts.get(cage.cage_type, 0) + 1
    return counts


def seams_counts(path: Path, width: int) -> str:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
    env["OMP_NUM_THREADS"] = str(width)
    cmd = [
        "taskset", "-c", "0" if width == 1 else "0-7",
        str(ROOT / "bbdir-release" / "seams"),
        "--jobs", str(width), "--cutoff", "3.5", "--graph", "cutoff",
        "--signature", "512+51268+sh", "cages", str(path),
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        raise RuntimeError((done.stderr or done.stdout)[-500:])
    return done.stdout


def main() -> None:
    pos, lengths = ortho_images()
    write_dump(OUT, pos, lengths)
    print(f"cores {os.cpu_count()}", flush=True)
    print(f"samples {N}", flush=True)
    print(
        f"sH orthorhombic supercell molecules {len(pos)} "
        f"edges {lengths[0]:.6f} {lengths[1]:.6f} {lengths[2]:.6f}",
        flush=True,
    )
    from sqq.io.input.lammps import inspect_lammps_dump, normalize_lammps_config
    cfg = normalize_lammps_config(
        {"units": "real", "atom_style": "molecular", "type_map": {"1": ["SOL", "OW"]}}
    )
    inspect_lammps_dump(OUT, cfg)
    print("sqq accepts the orthorhombic supercell", flush=True)
    text = seams_counts(OUT, 1)
    print("seams", " ".join(line.strip() for line in text.splitlines() if line.strip()), flush=True)
    ours = sqq_counts(pos, lengths, 1)
    print(f"sqq cages {ours}", flush=True)
    for width in (1, 8):
        os.sched_setaffinity(0, set(range(width)))

        def sqq_once(width=width):
            return sqq_counts(pos, lengths, width)

        def seams_once(width=width):
            seams_counts(OUT, width)

        sqq_once()
        seams_once()
        sqq_samples = []
        seams_samples = []
        for _ in range(N):
            t0 = time.perf_counter_ns()
            sqq_once()
            sqq_samples.append((time.perf_counter_ns() - t0) / 1e6)
            t0 = time.perf_counter_ns()
            seams_once()
            seams_samples.append((time.perf_counter_ns() - t0) / 1e6)
        for label, samples in (("sqq", sqq_samples), ("seams", seams_samples)):
            ordered = sorted(samples)
            med = 0.5 * (ordered[N // 2 - 1] + ordered[N // 2])
            print(
                f"{label} sH-ortho threads {width} molecules {len(pos)} "
                f"median {med:.4f} min {ordered[0]:.4f} samples {N}",
                flush=True,
            )


if __name__ == "__main__":
    main()
