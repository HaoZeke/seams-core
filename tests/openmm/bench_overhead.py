# Copyright (c) d-SEAMS developers
"""Serial wall time of plain OpenMM against the reporter.

Two TIP4P/Ew boxes, report intervals 10, 100 and 1000, n paired trials.
One warmup of each leg is discarded per box. Trials run one after another
on one thread. The overhead of a trial is the in-situ time minus the plain
time. Dumps stay in a temporary directory.
"""

import os
import statistics
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "adapters" / "python"))
sys.path.insert(0, str(ROOT / "tests" / "openmm"))

from dseams_insitu import DseamsReporter, Engine  # noqa: E402
from dseams_insitu.engine import load  # noqa: E402

from water import build_simulation, restore, save_state  # noqa: E402

LIB = Path(
    os.environ.get("SEAMS_LIB", ROOT / "bbdir-eessi" / "bindings" / "libseams_lang.so")
)
BOXES_NM = (2.0, 4.93)
INTERVALS = (10, 100, 1000)
STEPS = 1000
N = 6


def hinges(values):
    """Median, and the Tukey hinge (upper-half median minus lower-half median)."""
    ordered = sorted(values)
    half = len(ordered) // 2
    lower = ordered[:half]
    upper = ordered[half:] if len(ordered) % 2 == 0 else ordered[half + 1 :]
    return statistics.median(ordered), statistics.median(upper) - statistics.median(lower)


def time_plain(simulation, state):
    simulation.reporters.clear()
    restore(simulation, state)
    start = time.perf_counter()
    simulation.step(STEPS)
    return time.perf_counter() - start


def time_insitu(simulation, state, lib, dump, record, indices, interval):
    simulation.reporters.clear()
    restore(simulation, state)
    with Engine(lib, record, cutoff=3.5) as engine:
        simulation.reporters.append(
            DseamsReporter(engine, dump, interval=interval, indices=indices)
        )
        start = time.perf_counter()
        simulation.step(STEPS)
        elapsed = time.perf_counter() - start
    simulation.reporters.clear()
    return elapsed


def reports_in(path):
    text = Path(path).read_text(encoding="ascii")
    return sum(1 for line in text.splitlines() if line.startswith("step "))


def run_box(box_nm, lib, tmp):
    simulation, oxygens = build_simulation(box_nm)
    waters = simulation.topology.getNumResidues()
    state = save_state(simulation)
    warm = Path(tmp) / f"warm-{box_nm}"
    time_plain(simulation, state)
    time_insitu(
        simulation, state, lib, str(warm.with_suffix(".lammpstrj")),
        str(warm.with_suffix(".txt")), oxygens, INTERVALS[-1],
    )
    rows = []
    for interval in INTERVALS:
        plain = []
        live = []
        for trial in range(N):
            plain.append(time_plain(simulation, state))
            dump = Path(tmp) / f"{box_nm}-{interval}-{trial}.lammpstrj"
            record = Path(tmp) / f"{box_nm}-{interval}-{trial}.txt"
            elapsed = time_insitu(
                simulation, state, lib, str(dump), str(record), oxygens, interval
            )
            if trial == 0:
                got = reports_in(record)
                expect = STEPS // interval
                if got != expect:
                    raise SystemExit(
                        f"{waters} waters interval {interval}: "
                        f"{got} reports, expected {expect}"
                    )
            dump.unlink(missing_ok=True)
            record.unlink(missing_ok=True)
            live.append(elapsed)
            print(
                f"waters {waters} interval {interval} trial {trial + 1} "
                f"plain {plain[-1]:.6f} insitu {elapsed:.6f}",
                flush=True,
            )
        overhead = [b - a for a, b in zip(plain, live)]
        rows.append((interval, plain, live, overhead))
    return waters, rows


def emit(waters, rows):
    print(f"box waters {waters} steps {STEPS} n {N}")
    for interval, plain, live, overhead in rows:
        for name, values in ("plain", plain), ("insitu", live), ("overhead", overhead):
            median, iqr = hinges(values)
            print(
                f"waters {waters} interval {interval} {name} "
                f"median {median:.6f} iqr {iqr:.6f}"
            )
            print(" ".join(f"{v:.6f}" for v in values))


def main():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["NO_COLOR"] = "1"
    lib = load(str(LIB))
    with tempfile.TemporaryDirectory(prefix="dseams-openmm-bench-") as tmp:
        for box_nm in BOXES_NM:
            waters, rows = run_box(box_nm, lib, tmp)
            emit(waters, rows)


if __name__ == "__main__":
    main()
