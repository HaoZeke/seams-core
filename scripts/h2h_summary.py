#!/usr/bin/env python3
"""One sitting: vesin, freud, pyscal3, and mdapy against d-SEAMS.

Neighbours are the 3.5 A cutoff list. Pair counts must match.
q6 is the cutoff Steinhardt parameter and its neighbour average.
CHILL+ is mdapy's modifier against d-SEAMS ovito mode, the same
neighbour rule. pyscal3 has no CHILL+.

Cells are orthorhombic ice Ih at 4096, 110592, and 884736 molecules.
Eight alternating samples. Serial is core 0. Eight threads are cores
0-7. A win is a strictly smaller median and a strictly smaller minimum.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = 8
CELLS = (8, 24, 48)  # 4096, 110592, 884736 molecules
FIELDS = [
    "competitor", "operation", "molecules", "threads", "cores", "cpus", "n",
    "dseams_median_ms", "dseams_min_ms", "other_median_ms", "other_min_ms",
    "ratio", "verdict", "agreement", "budget_s", "largest_finished",
    "version", "cpu_model",
]


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def pin(width: int) -> None:
    os.sched_setaffinity(0, {0} if width == 1 else set(range(8)))


def median(samples: list[float]) -> float:
    ordered = sorted(samples)
    mid = len(ordered) // 2
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def timed(fn) -> float:
    t0 = time.perf_counter_ns()
    fn()
    return (time.perf_counter_ns() - t0) / 1e6


def abab(ours, other) -> tuple[list[float], list[float]]:
    ours()
    other()
    left, right = [], []
    for i in range(SAMPLES):
        order = (("d", ours), ("o", other)) if i % 2 == 0 else (("o", other), ("d", ours))
        got = {}
        for name, fn in order:
            got[name] = timed(fn)
        left.append(got["d"])
        right.append(got["o"])
    return left, right


def verdict(left: list[float], right: list[float]) -> str:
    return "WIN" if median(left) < median(right) and min(left) < min(right) else "LOSE"


def row_of(competitor, operation, molecules, threads, cores, left, right, agreement, version, model):
    other_med = median(right)
    our_med = median(left)
    ratio = other_med / our_med if our_med > 0 else ""
    return {
        "competitor": competitor,
        "operation": operation,
        "molecules": molecules,
        "threads": threads,
        "cores": cores,
        "cpus": "0" if threads == 1 else "0-7",
        "n": SAMPLES,
        "dseams_median_ms": f"{our_med:.4f}",
        "dseams_min_ms": f"{min(left):.4f}",
        "other_median_ms": f"{other_med:.4f}",
        "other_min_ms": f"{min(right):.4f}",
        "ratio": f"{ratio:.4f}" if ratio != "" else "",
        "verdict": verdict(left, right),
        "agreement": agreement,
        "budget_s": "",
        "largest_finished": "",
        "version": version,
        "cpu_model": model,
    }


def run_width(threads: int) -> list[dict]:
    pin(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    sys.path.insert(0, str(ROOT / "scripts"))
    import ctypes

    import freud
    import mdapy
    import numpy as np
    import pyscal3
    import vesin
    from ase import Atoms

    import h2h_mdapy_pyscal3_vesin as h

    pin(threads)
    lib = h.load_seams(h.library_path())
    lib.seams_version_string.restype = ctypes.c_char_p
    seams = lib.seams_version_string().decode()
    cores = os.cpu_count() or 1
    model = cpu_model()
    versions = {
        "vesin": getattr(vesin, "__version__", "vesin"),
        "freud": freud.__version__,
        "mdapy": mdapy.__version__,
        "pyscal3": getattr(pyscal3, "__version__", "pyscal3"),
    }
    rows = []
    # pyscal3 exports no CHILL or CHILL+ function. The reason cites the API.
    from h2h_summary_v2 import pyscal_chill_reason

    rows.append({
        "competitor": "pyscal3",
        "operation": "chill_plus",
        "molecules": "",
        "threads": threads,
        "cores": cores,
        "cpus": "0" if threads == 1 else "0-7",
        "n": 0,
        "dseams_median_ms": "",
        "dseams_min_ms": "",
        "other_median_ms": "",
        "other_min_ms": "",
        "ratio": "",
        "verdict": "absent",
        "agreement": pyscal_chill_reason(),
        "version": f"{seams}/{versions['pyscal3']}",
        "cpu_model": model,
    })

    for ncell in CELLS:
        pin(threads)
        pos, lengths = h.ice_positions(ncell)
        n = int(pos.shape[0])
        lx, ly, lz = lengths
        box = np.diag(np.asarray(lengths, dtype=np.float64))
        fbox = freud.box.Box(lx, ly, lz)
        system = mdapy.System(pos, box=lengths)
        atoms = Atoms(positions=pos, cell=lengths, pbc=True)
        held = h.frame_system(lib, pos, lengths)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        d_off, d_col = h.dseams_cutoff(lib, held)
        d_pairs = h.pairs_from_csr(d_off, d_col)
        v_pairs = h.vesin_pairs(pos, box, threads)
        if not h.same_edges(d_pairs, v_pairs):
            raise SystemExit(
                f"vesin pairs {v_pairs.shape[0]} differ from d-SEAMS {d_pairs.shape[0]} at {n}"
            )
        found = freud.locality.AABBQuery(fbox, pos).query(
            pos, {"r_max": h.CUTOFF, "exclude_ii": True}
        ).toNeighborList()
        f_pairs = h.undirected(np.asarray(found[:, 0]), np.asarray(found[:, 1]))
        if not h.same_edges(d_pairs, f_pairs):
            raise SystemExit(
                f"freud pairs {f_pairs.shape[0]} differ from d-SEAMS {d_pairs.shape[0]} at {n}"
            )
        pair_note = f"pairs {d_pairs.shape[0]} equal"

        q = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()
        h.dseams_q6(lib, held, q, bar)
        q_np = np.ctypeslib.as_array(q).copy()
        bar_np = np.ctypeslib.as_array(bar).copy()
        mdapy.set_num_threads(threads)
        shell = mdapy.NeighborList.build(system, cutoff=h.CUTOFF)
        local = mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=False).compute(
            system, neighbors=shell
        )
        averaged = mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=True).compute(
            system, neighbors=shell
        )
        mq = np.asarray(local.q, dtype=np.float64).reshape(n, -1)[:, 0]
        mbar = np.asarray(averaged.q, dtype=np.float64).reshape(n, -1)[:, 0]
        q_diff = max(h.max_abs(q_np, mq), h.max_abs(bar_np, mbar))
        if q_diff > 1e-8:
            raise SystemExit(f"mdapy q6 max_abs {q_diff:.3e} at {n}")
        pyscal3.set_num_threads(threads)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=h.CUTOFF, store_rows=False)
        pyscal3.steinhardt_parameter(atoms, 6, averaged=True)
        pq = np.asarray(atoms.arrays["pyscal_q6"], dtype=np.float64)
        pbar = np.asarray(atoms.arrays["pyscal_avg_q6"], dtype=np.float64)
        p_diff = max(h.max_abs(q_np, pq), h.max_abs(bar_np, pbar))
        if p_diff > 1e-8:
            raise SystemExit(f"pyscal3 q6 max_abs {p_diff:.3e} at {n}")

        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()
        h.dseams_chill_plus(lib, held, labels, counts, h.SEAMS_CHILL_OVITO)
        ours = np.array([h.DSEAMS_AS_OVITO[int(v)] for v in np.ctypeslib.as_array(labels)])
        ice = mdapy.identify.ChillPlus(cutoff=h.CUTOFF).compute(system)
        md_names = np.array([tuple(ice.names)[int(v)] for v in np.asarray(ice.structure)])
        if not np.array_equal(ours, md_names):
            differ = int(np.count_nonzero(ours != md_names))
            raise SystemExit(f"mdapy CHILL+ {differ} atoms differ at {n}")
        chill_note = f"labels {n} equal q6_max_abs {q_diff:.3e}"
        q_note = f"q6_max_abs {q_diff:.3e} q6bar_max_abs {h.max_abs(bar_np, mbar):.3e}"
        pq_note = f"q6_max_abs {p_diff:.3e} q6bar_max_abs {h.max_abs(bar_np, pbar):.3e}"

        nnz = ctypes.c_size_t()
        n_off = ctypes.c_size_t()

        def dseams_neighbours():
            pin(threads)
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            nb = ctypes.c_void_p()
            h.require(
                lib,
                lib.seams_neighbours_build(held, h.CUTOFF, 1, h.SEAMS_GRAPH_CUTOFF, ctypes.byref(nb)),
                "cutoff",
            )
            h.require(
                lib,
                lib.seams_neighbours_copy(
                    nb, None, 0, ctypes.byref(n_off), None, 0, ctypes.byref(nnz)
                ),
                "cutoff length",
            )
            lib.seams_neighbours_destroy(nb)
            if nnz.value == 0:
                raise SystemExit("empty cutoff")

        def vesin_neighbours():
            pin(threads)
            calc = vesin.NeighborList(cutoff=h.CUTOFF, full_list=True, n_threads=threads)
            left, _right = calc.compute(
                points=pos, box=box, periodic=True, quantities="ij", copy=False
            )
            if int(left.shape[0]) < n:
                raise SystemExit("vesin short")

        def freud_neighbours():
            pin(threads)
            freud.parallel.set_num_threads(threads)
            query = freud.locality.AABBQuery(fbox, pos).query(
                pos, {"r_max": h.CUTOFF, "exclude_ii": True}
            ).toNeighborList()
            if int(query.neighbor_counts.shape[0]) != n:
                raise SystemExit("freud short")

        def dseams_q():
            pin(threads)
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            h.dseams_q6(lib, held, q, bar)

        def mdapy_q():
            pin(threads)
            mdapy.set_num_threads(threads)
            built = mdapy.NeighborList.build(system, cutoff=h.CUTOFF)
            mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=False).compute(
                system, neighbors=built
            )
            mdapy.analysis.SteinhardtParameters((6,), cutoff=h.CUTOFF, average=True).compute(
                system, neighbors=built
            )

        def pyscal_q():
            pin(threads)
            pyscal3.set_num_threads(threads)
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff=h.CUTOFF, store_rows=False)
            pyscal3.steinhardt_parameter(atoms, 6, averaged=True)

        def dseams_ice():
            pin(threads)
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            h.dseams_chill_plus(lib, held, labels, counts, h.SEAMS_CHILL_OVITO)

        def mdapy_ice():
            pin(threads)
            mdapy.set_num_threads(threads)
            result = mdapy.identify.ChillPlus(cutoff=h.CUTOFF).compute(system)
            if int(result.structure.shape[0]) != n:
                raise SystemExit("mdapy chill short")

        jobs = (
            ("vesin", "neighbours", dseams_neighbours, vesin_neighbours, pair_note),
            ("freud", "neighbours", dseams_neighbours, freud_neighbours, pair_note),
            ("mdapy", "q6", dseams_q, mdapy_q, q_note),
            ("pyscal3", "q6", dseams_q, pyscal_q, pq_note),
            ("mdapy", "chill_plus", dseams_ice, mdapy_ice, chill_note),
        )
        for competitor, operation, ours, other, note in jobs:
            left, right = abab(ours, other)
            made = row_of(
                competitor, operation, n, threads, cores, left, right, note,
                f"{seams}/{versions[competitor]}", model,
            )
            rows.append(made)
            print(
                f"threads {threads} {competitor} {operation} molecules {n} "
                f"dseams {made['dseams_median_ms']} other {made['other_median_ms']} "
                f"{made['verdict']}",
                file=sys.stderr,
                flush=True,
            )
        from h2h_summary_extra import bind, extra_rows
        bind(lib)
        rows.extend(extra_rows(
            lib, held, pos, lengths, threads, cores, model, seams, pin, abab, row_of, h,
        ))
        lib.seams_system_destroy(held)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--worker", type=int, choices=(1, 8))
    args = ap.parse_args()
    if not args.worker:
        from build_ptm_ice import build
        build()
    if args.worker:
        rows = run_width(args.worker)
        writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        if any(row["verdict"] == "LOSE" for row in rows):
            raise SystemExit(1)
        return

    rows = []
    for threads, cores in ((1, "0"), (8, "0-7")):
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = (
            "/usr/lib/x86_64-linux-gnu" + os.pathsep + str(ROOT / "bbdir-release" / "src")
        )
        env["OVITO_THREAD_COUNT"] = str(threads)
        env["QT_QPA_PLATFORM"] = "offscreen"
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_DYNAMIC"] = "FALSE"
        env["OMP_PROC_BIND"] = "close"
        env["OMP_PLACES"] = "cores"
        env["OMP_WAIT_POLICY"] = "PASSIVE"
        env["QT_QPA_PLATFORM"] = "offscreen"
        env["PYTHONUNBUFFERED"] = "1"
        done = subprocess.run(
            ["taskset", "-c", cores, sys.executable, "-u", str(Path(__file__).resolve()),
             "--csv", args.csv, "--worker", str(threads)],
            cwd=ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
        )
        if done.returncode not in (0, 1):
            sys.stderr.write(done.stdout[-2000:])
            raise SystemExit(done.returncode or 1)
        lines = done.stdout.splitlines()
        start = next((i for i, line in enumerate(lines) if line.startswith("competitor,")), None)
        if start is None:
            sys.stderr.write(done.stdout[-2000:])
            raise SystemExit(f"no csv header from threads {threads}")
        parsed = list(csv.DictReader(lines[start:]))
        if not parsed:
            sys.stderr.write(done.stdout[-2000:])
            raise SystemExit(f"no rows from threads {threads}")
        rows.extend(parsed)
    if any(row["verdict"] == "LOSE" for row in rows):
        _write(Path(args.csv), rows)
        raise SystemExit(1)
    cage_path = ROOT / "benchmarks" / "h2h-clavis-cages.csv"
    if cage_path.is_file():
        with cage_path.open(encoding="ascii", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    _write(Path(args.csv), rows)
    print(f"wrote {Path(args.csv).name} rows {len(rows)}", file=sys.stderr, flush=True)


def _write(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
