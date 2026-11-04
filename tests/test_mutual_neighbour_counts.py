#!/usr/bin/env python3
"""The skewed 512-molecule frame has one mutual four-neighbour count.

d-SEAMS, freud and OVITO each keep a neighbour only when it is among the
four nearest of both molecules. The directed edge counts are equal.
"""

import os
import subprocess
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bench_h2h_neighbours import (  # noqa: E402
    make_skew_narrow_frame,
    mutual_directed_count,
    mutual_from_slots,
    write_lammps_dump,
)

EXPECTED = 2016
SEED = 7


def freud_count(frame) -> int:
    import freud

    freud.parallel.set_num_threads(1)
    box = freud.box.Box.from_matrix(frame.H)
    points = box.wrap(frame.positions - frame.origin)
    neighbours = (
        freud.locality.AABBQuery(box, points)
        .query(points, {"mode": "nearest", "num_neighbors": 4, "exclude_ii": True})
        .toNeighborList()
    )
    return mutual_directed_count(neighbours.query_point_indices, neighbours.point_indices)


def ovito_count(frame) -> int:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    os.environ["OVITO_THREAD_COUNT"] = "1"
    from ovito.data import DataCollection, NearestNeighborFinder

    data = DataCollection()
    data.create_cell(frame.H, pbc=(True, True, True))
    parts = data.create_particles(count=frame.positions.shape[0])
    parts.create_property("Position", data=frame.positions - frame.origin)
    slots = NearestNeighborFinder(4, data).find_all()[0]
    return mutual_from_slots(slots)


def main() -> int:
    frame = make_skew_narrow_frame(SEED)
    freud_edges = freud_count(frame)
    ovito_edges = ovito_count(frame)
    if freud_edges != EXPECTED or ovito_edges != EXPECTED:
        print(f"freud {freud_edges} ovito {ovito_edges} expected {EXPECTED}", file=sys.stderr)
        return 1
    binary = Path("bbdir/tests/bench_neigh_only")
    if binary.is_file():
        dump = Path("input/traj/skew_triclinic_width_lt_cutoff.lammpstrj")
        write_lammps_dump(str(dump), frame, atom_type=1)
        proc = subprocess.run(
            [
                "taskset", "-c", "0",
                str(binary),
                "--traj", str(dump),
                "--frame", "1",
                "--type", "1",
                "--reps", "1",
                "--warmup", "0",
                "--cutoff", "3.5",
                "--k", "4",
                "--threads", "1",
            ],
            check=False,
            capture_output=True,
            text=True,
            env={**os.environ, "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
        )
        if proc.returncode != 0:
            print(proc.stderr, file=sys.stderr)
            return proc.returncode
        edges = None
        for line in proc.stdout.splitlines():
            if line.startswith("knn_directed_edges "):
                edges = int(line.split()[1])
        if edges != EXPECTED:
            print(f"d-SEAMS {edges} expected {EXPECTED}", file=sys.stderr)
            return 1
    print(f"mutual directed edges {EXPECTED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
