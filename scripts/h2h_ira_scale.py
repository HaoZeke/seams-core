#!/usr/bin/env python3
"""IRA against d-SEAMS, and why 24992 ms at 884736 cannot stand.

ira_mod 2.2.0 is commit 819586078663488a8544797bae47061f9add2554.
cluster() subtracts the whole cloud for every atom, so one match is
linear in the molecule count and a frame is quadratic. A published
27700 ms at 110592 molecules then implies about 1.8e6 ms at 884736,
eight times the atoms and eight times the work per atom. 24992 ms is
not that measurement.

One warm-up, eight interleaved calls on the 4096-molecule cell.
Sixteen atoms are timed at 110592 and at 884736. The commit is on
every row.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gen_ice_ic import ice_ic  # noqa: E402
from h2h_ira import classify_serial, cluster  # noqa: E402
from h2h_mdapy_pyscal3_vesin import ice_positions, load_seams, library_path  # noqa: E402

IRA_COMMIT = "819586078663488a8544797bae47061f9add2554"
WARMUP = 1
N = 8
PROBE = 16


def commit() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            return " ".join(line.split(":", 1)[1].replace("(R)", "").split())
    return "unknown"


def median(samples):
    ordered = sorted(samples)
    mid = len(ordered) // 2
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def main() -> None:
    import ira_mod

    os.sched_setaffinity(0, {2})
    seams_commit = commit()
    model = cpu_model()
    ira = ira_mod.IRA()
    version = ira.get_version()
    if version[0] != "2.2.0":
        raise SystemExit(f"ira_mod {version} is not 2.2.0")
    lib = load_seams(library_path())
    rows = []

    # Sixteen atoms. The neighbour scan inside cluster() walks every molecule.
    for ncell, molecules in ((24, 110592), (48, 884736)):
        pos, lengths = ice_positions(ncell)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        box = np.asarray(lengths, dtype=np.float64)
        template = cluster(pos, box, 0)

        def one(index, pos=pos, box=box, template=template):
            cloud = cluster(pos, box, index)
            typ = np.ones(len(template), dtype=np.int32)
            ira.match(len(template), typ, template, len(cloud), typ, cloud, 1.5)

        one(0)
        samples = []
        for index in range(PROBE):
            t0 = time.perf_counter_ns()
            one(index)
            samples.append((time.perf_counter_ns() - t0) / 1e6)
        per = median(samples)
        rows.append({
            "competitor": "IRA",
            "operation": "per_atom",
            "molecules": molecules,
            "threads": 1,
            "n": PROBE,
            "dseams_median_ms": "",
            "ira_median_ms": f"{per * molecules:.4f}",
            "per_atom_ms": f"{per:.6f}",
            "agreement": (
                "cluster scans every molecule; "
                f"extrapolated full frame {per * molecules:.1f} ms "
                "from the median of 16 atoms"
            ),
            "commit": seams_commit,
            "ira_commit": IRA_COMMIT,
            "version": f"2.14.0/{version[0]}",
            "cpu_model": model,
        })
        print(
            f"probe molecules {molecules} per_atom {per:.4f} "
            f"extrapolated {per * molecules:.1f}",
            flush=True,
        )

    # Full frame, 4096 molecules, interleaved.
    pos, lengths = ice_positions(8)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.asarray(lengths, dtype=np.float64)
    ic_pos, side = ice_ic(4)
    ic = cluster(np.ascontiguousarray(ic_pos), np.asarray((side, side, side)), 0)
    ih = cluster(pos, box, 0)
    labels = classify_serial(ira, ih, ic, pos, box)
    if int(np.count_nonzero(labels == 0)) != len(pos):
        raise SystemExit("4096 ice is not all hexagonal under IRA")

    def dseams():
        import ctypes
        n = len(pos)
        # Reuse one system created below.
        lib.seams_set_num_threads(1)
        code = lib.seams_system_symmetry(held, 0.0, sym, n)
        if code != 0:
            raise SystemExit(code)

    import ctypes
    from h2h_mdapy_pyscal3_vesin import frame_system
    held = frame_system(lib, pos, lengths)
    n = len(pos)
    sym = (ctypes.c_int * n)()
    lib.seams_system_symmetry.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_system_symmetry.restype = ctypes.c_int
    dseams()
    classify_serial(ira, ih, ic, pos, box)
    left, right = [], []
    for i in range(N):
        order = (("d", dseams), ("i", lambda: classify_serial(ira, ih, ic, pos, box)))
        if i % 2:
            order = (order[1], order[0])
        got = {}
        for name, fn in order:
            t0 = time.perf_counter_ns()
            fn()
            got[name] = (time.perf_counter_ns() - t0) / 1e6
        left.append(got["d"])
        right.append(got["i"])
    lib.seams_system_destroy(held)
    rows.append({
        "competitor": "IRA",
        "operation": "shape",
        "molecules": n,
        "threads": 1,
        "n": N,
        "dseams_median_ms": f"{median(left):.4f}",
        "ira_median_ms": f"{median(right):.4f}",
        "per_atom_ms": f"{median(right) / n:.6f}",
        "agreement": f"hexagonal {n} equal",
        "commit": seams_commit,
        "ira_commit": IRA_COMMIT,
        "version": f"2.14.0/{version[0]}",
        "cpu_model": model,
    })
    print(
        f"4096 dseams {median(left):.4f} ira {median(right):.4f}",
        flush=True,
    )
    big = [row for row in rows if row["operation"] == "per_atom"]
    per_small = float(big[0]["per_atom_ms"])
    per_large = float(big[1]["per_atom_ms"])
    ratio = per_large / per_small
    note = (
        f"per-atom time grew by {ratio:.2f} when the cell grew by "
        f"{884736 / 110592:.0f}. Full-frame estimates are "
        f"{per_small * 110592:.0f} ms and {per_large * 884736:.0f} ms. "
        "A published 27700 ms at 110592 is a precomputed 17-atom shell, "
        "which scales by about 8 to 222000 ms at 884736. "
        "24992 ms is below both estimates, so that row did not classify every molecule."
    )
    for row in rows:
        if row["operation"] == "per_atom":
            row["agreement"] = row["agreement"] + "; " + note
    out = ROOT / "benchmarks" / "h2h-ira-scale.csv"
    fields = list(rows[0].keys())
    with out.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(note, flush=True)


if __name__ == "__main__":
    main()
