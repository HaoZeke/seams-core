#!/usr/bin/env python3
"""q6 at 110592 molecules and 8 threads stays at least 1.5x mdapy.

The values stay within 1e-12 of mdapy's cutoff Steinhardt parameter.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.setdefault("LD_LIBRARY_PATH", str(ROOT / "bbdir-release" / "src"))
os.environ["OMP_NUM_THREADS"] = "8"
os.environ["OMP_PROC_BIND"] = "close"
os.environ["OMP_PLACES"] = "cores"
os.environ["OMP_WAIT_POLICY"] = "PASSIVE"

import mdapy  # noqa: E402
import numpy as np  # noqa: E402

import h2h_mdapy_pyscal3_vesin as h  # noqa: E402


def main() -> None:
    os.sched_setaffinity(0, set(range(8)))
    lib = h.load_seams(h.library_path())
    pos, lengths = h.ice_positions(24)
    n = len(pos)
    if n != 110592:
        raise SystemExit(f"expected 110592 molecules, got {n}")
    held = h.frame_system(lib, pos, lengths)
    h.require(lib, lib.seams_set_num_threads(8), "threads")
    import ctypes

    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    h.dseams_q6(lib, held, q, bar)
    system = mdapy.System(pos, box=lengths)
    mdapy.set_num_threads(8)
    built = mdapy.NeighborList.build(system, cutoff=h.CUTOFF)
    local = mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=False).compute(
        system, neighbors=built
    )
    averaged = mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=True).compute(
        system, neighbors=built
    )
    mq = np.asarray(local.q, dtype=np.float64).reshape(n, -1)[:, 0]
    mbar = np.asarray(averaged.q, dtype=np.float64).reshape(n, -1)[:, 0]
    diff = max(h.max_abs(np.ctypeslib.as_array(q), mq), h.max_abs(np.ctypeslib.as_array(bar), mbar))
    if diff > 1e-12:
        raise SystemExit(f"q6 max_abs {diff:.3e} exceeds 1e-12")

    def dseams():
        h.dseams_q6(lib, held, q, bar)

    def other():
        neigh = mdapy.NeighborList.build(system, cutoff=h.CUTOFF)
        mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=False).compute(
            system, neighbors=neigh
        )
        mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=True).compute(
            system, neighbors=neigh
        )

    dseams()
    other()
    ds, md = [], []
    for i in range(8):
        if i % 2 == 0:
            t0 = time.perf_counter()
            dseams()
            ds.append((time.perf_counter() - t0) * 1e3)
            t0 = time.perf_counter()
            other()
            md.append((time.perf_counter() - t0) * 1e3)
        else:
            t0 = time.perf_counter()
            other()
            md.append((time.perf_counter() - t0) * 1e3)
            t0 = time.perf_counter()
            dseams()
            ds.append((time.perf_counter() - t0) * 1e3)
    ds.sort()
    md.sort()
    dmed = 0.5 * (ds[3] + ds[4])
    mmed = 0.5 * (md[3] + md[4])
    ratio = mmed / dmed
    print(
        f"q6 110592 threads 8 dseams {dmed:.4f} mdapy {mmed:.4f} ratio {ratio:.3f} max_abs {diff:.3e}",
        flush=True,
    )
    if ratio < 1.5 or dmed >= mmed or min(ds) >= min(md):
        raise SystemExit(f"q6 ratio {ratio:.3f} is below 1.5")


if __name__ == "__main__":
    main()
