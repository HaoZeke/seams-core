"""Time OVITO's ring finder against the d-SEAMS ring driver on shared frames.

Rebuilds the frames of ``h2h_campaign.py prepare`` (ideal ice Ic and Ih,
liquid water, confined water at four scales), then times on each frame

* ``dseams-rings``: the ``rings`` mode of ``tests/bench_paper`` (3.5 A
  cutoff graph and primitive rings of at most six members), and
* ``ovito-rings``: OVITO's ``CreateBondsModifier`` at 3.5 A followed by
  ``FindRingsModifier`` for ring sizes 3 to 6, on a fresh data collection
  built outside the timer.

Run on one pinned core with ``OVITO_THREAD_COUNT=1``. Both methods read the
same LAMMPS file, each in its own child process, so the OVITO runtime is never
resident while the d-SEAMS driver runs on the shared core. Writes one row per
method and frame with every sample and the ring census.

Usage: python3 h2h_ovito_rings.py BENCH LIQUID_DUMP CONFINED_DUMP OUT.tsv
"""

import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("h2h", HERE / "h2h_campaign.py")
h2h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h2h)

# Warm-up and timed calls of both methods at each scale.
PLAN = {"1e3": (5, 30), "1e4": (5, 30), "1e5": (5, 30), "1e6": (1, 10)}


def frame(system, scale, liquid, confined):
    if system == "ic":
        return h2h.cubic_diamond(h2h.IDEAL_REPEATS[scale])
    if system == "ih":
        return h2h.lonsdaleite(h2h.IDEAL_REPEATS[scale])
    if system == "liquid":
        rep = h2h.LIQUID_REPEATS[scale]
        return h2h.replicate(*liquid, (rep, rep, rep))
    return h2h.replicate(*confined, h2h.CONFINED_REPEATS[scale])


def dseams(bench, trajectory, warmup, repeats, work):
    sidecar = work / "census.txt"
    process = subprocess.run(
        ["taskset", "-c", str(h2h.CORE), str(bench), str(trajectory),
         str(h2h.DRIVER_FRAME), str(warmup), str(repeats), "rings", str(sidecar)],
        check=True, capture_output=True, text=True,
    )
    samples = h2h.parse_driver_samples(process.stdout)
    census = dict(line.split() for line in sidecar.read_text().splitlines())
    return samples, [int(census.get(str(size), 0)) for size in h2h.RING_SIZES]


def ovito(positions, box, warmup, repeats):
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier

    samples = []
    data = None
    for repetition in range(warmup + repeats):
        data = h2h.ovito_data(positions, box)
        start = time.perf_counter_ns()
        data.apply(CreateBondsModifier(cutoff=h2h.CUTOFF_A))
        data.apply(FindRingsModifier(minimum_ring_size=h2h.MIN_RING,
                                     maximum_ring_size=h2h.MAX_RING))
        elapsed = (time.perf_counter_ns() - start) / h2h.NS_PER_MS
        if repetition >= warmup:
            samples.append(elapsed)
    census = [0] * len(h2h.RING_SIZES)
    for size in h2h.RING_SIZES:
        table = data.tables.get(f"{size}-rings")
        census[size - h2h.MIN_RING] = 0 if table is None else table.count
    return samples, census


def ovito_child(trajectory, warmup, repeats):
    process = subprocess.run(
        ["taskset", "-c", str(h2h.CORE), sys.executable, str(pathlib.Path(__file__).resolve()),
         "--ovito-case", str(trajectory), str(warmup), str(repeats)],
        check=True, capture_output=True, text=True,
    )
    result = json.loads(process.stdout.strip().splitlines()[-1])
    return result["samples"], result["census"]


def ovito_case(trajectory, warmup, repeats):
    positions, box = h2h.read_lammps(pathlib.Path(trajectory), frame=h2h.DRIVER_FRAME)
    samples, census = ovito(positions, box, int(warmup), int(repeats))
    sys.stdout.write(json.dumps({"samples": samples, "census": census}) + "\n")


def main():
    if sys.argv[1] == "--ovito-case":
        ovito_case(*sys.argv[2:5])
        return
    bench, liquid_dump, confined_dump, out = map(pathlib.Path, sys.argv[1:5])
    liquid = h2h.read_lammps(liquid_dump, frame=h2h.SOURCE_FRAME)
    confined = h2h.read_lammps(confined_dump, frame=h2h.SOURCE_FRAME,
                               atom_type=h2h.CONFINED_OXYGEN_TYPE)
    rows = ["system\tscale\tmolecules\tmethod\trepeats\trings_3_to_6\tsamples_ms"]
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        for system in h2h.SYSTEMS:
            for scale in h2h.SCALES:
                positions, box = frame(system, scale, liquid, confined)
                trajectory = work / f"{system}-{scale}.lammpstrj"
                h2h.write_lammps(trajectory, positions, box)
                warmup, repeats = PLAN[scale]
                for method in ("dseams-rings", "ovito-rings"):
                    if method == "dseams-rings":
                        samples, census = dseams(bench, trajectory, warmup, repeats, work)
                    else:
                        samples, census = ovito_child(trajectory, warmup, repeats)
                    rows.append("\t".join((
                        system, scale, str(len(positions)), method, str(len(samples)),
                        ",".join(map(str, census)),
                        ",".join(f"{x:.6f}" for x in samples),
                    )))
                    sys.stdout.write(rows[-1][:120] + "\n")
                    sys.stdout.flush()
                trajectory.unlink()
    out.write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
