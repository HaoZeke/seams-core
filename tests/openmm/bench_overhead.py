# Copyright (c) d-SEAMS developers
"""Serial wall time of plain OpenMM against the same run with the reporter.

One warmup of each leg is discarded. The recorded sample is n paired
trials, run one after another on one thread. The overhead of a trial is
the in-situ time minus the plain time.
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
STEPS = 100
INTERVAL = 10
N = 6


def hinges(values):
    """Median, and the Tukey hinge (upper-half median minus lower-half median)."""
    ordered = sorted(values)
    half = len(ordered) // 2
    lower = ordered[:half]
    upper = ordered[half:] if len(ordered) % 2 == 0 else ordered[half + 1 :]
    q1 = statistics.median(lower)
    q3 = statistics.median(upper)
    return statistics.median(ordered), q3 - q1


def time_plain(simulation, state):
    simulation.reporters.clear()
    restore(simulation, state)
    start = time.perf_counter()
    simulation.step(STEPS)
    return time.perf_counter() - start


def time_insitu(simulation, state, lib, dump, record, indices):
    simulation.reporters.clear()
    restore(simulation, state)
    with Engine(lib, record, cutoff=3.5) as engine:
        simulation.reporters.append(
            DseamsReporter(engine, dump, interval=INTERVAL, indices=indices)
        )
        start = time.perf_counter()
        simulation.step(STEPS)
        elapsed = time.perf_counter() - start
    simulation.reporters.clear()
    return elapsed


def main():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["NO_COLOR"] = "1"
    simulation, oxygens = build_simulation()
    state = save_state(simulation)
    lib = load(str(LIB))
    with tempfile.TemporaryDirectory(prefix="dseams-openmm-bench-") as tmp:
        warm_dump = str(Path(tmp) / "warm.lammpstrj")
        warm_record = str(Path(tmp) / "warm.txt")
        time_plain(simulation, state)
        time_insitu(simulation, state, lib, warm_dump, warm_record, oxygens)
        plain = []
        live = []
        for trial in range(N):
            plain.append(time_plain(simulation, state))
            dump = str(Path(tmp) / f"run-{trial}.lammpstrj")
            record = str(Path(tmp) / f"run-{trial}.txt")
            live.append(time_insitu(simulation, state, lib, dump, record, oxygens))
            print(f"trial {trial + 1} plain {plain[-1]:.6f} insitu {live[-1]:.6f}", flush=True)
    overhead = [b - a for a, b in zip(plain, live)]
    for name, values in ("plain", plain), ("insitu", live), ("overhead", overhead):
        median, iqr = hinges(values)
        print(f"{name} n={len(values)} median {median:.6f} iqr {iqr:.6f}")
        print(" ".join(f"{v:.6f}" for v in values))
    print(f"waters {simulation.topology.getNumResidues()} steps {STEPS} interval {INTERVAL}")


if __name__ == "__main__":
    main()
