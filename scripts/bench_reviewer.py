#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""One sitting for the paper timings the reviewer asked to see together.

Serial and eight threads, one warm-up then eight calls, from one release
build. The host name is not read. The core count and the CPU model are.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paper_numbers as pn  # noqa: E402

FIELDS = pn.FIELDS
CUTOFF = pn.CUTOFF


def pin(threads: int, args: list[str]) -> str:
    """One core for the serial arm, the first eight for the wide arm.

    Core 0 is left free while another sitting holds it: the serial arm
    uses core 1, and that affinity is what the row records.
    """
    import paper_stats as stats
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    cores = "1" if threads == 1 else stats.taskset_cores(threads)
    env["PAPER_TASKSET"] = cores
    cmd = ["taskset", "-c", cores, sys.executable, *args]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise SystemExit("\n".join(tail[-40:]) if tail else f"exit {done.returncode}")
    return done.stdout


def freud_rows(info) -> list[dict]:
    import freud
    import numpy as np
    import paper_ptm_templates as ptm
    threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
    freud.parallel.set_num_threads(threads)
    rows = []
    for molecules, cells in ((4096, 8), (110592, 24)):
        pos, lattice = ptm.ice_ic(cells)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        side = float(np.asarray(lattice).reshape(-1)[0])
        box = freud.box.Box.cube(side)
        pts = np.ascontiguousarray(pos - side / 2.0)

        def once(box=box, pts=pts, n=len(pos)):
            found = freud.locality.AABBQuery(box, pts).query(
                pts, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
            if int(found.neighbor_counts.shape[0]) != n:
                raise RuntimeError("freud neighbours")

        once()
        samples = [pn.timed(once) for _ in range(pn.REPEATS)]
        rows.append(pn.row_of(
            info, "reviewer", f"ice-Ic-{molecules}", "neighbours", "freud",
            molecules, threads, samples,
            f"freud.AABBQuery n={pn.REPEATS} after one warm-up; "
            f"get_num_threads={freud.parallel.get_num_threads()}",
            version=freud.__version__,
        ))
    return rows


def sqq_rows(info) -> list[dict]:
    import numpy as np
    import networkx as nx
    import freud
    import paper_ptm_templates as ptm
    from sqq.core.sqq_py.ring import find_rings
    import sqq
    threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
    pos, lattice = ptm.ice_ic(8)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    n = len(pos)
    side = float(np.asarray(lattice).reshape(-1)[0])
    query = freud.locality.AABBQuery(freud.box.Box.cube(side), pos - side / 2.0).query(
        pos - side / 2.0, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    graph.add_edges_from(
        (int(i), int(j)) for i, j in zip(query[:, 0], query[:, 1]) if int(i) < int(j))
    adj = {i: set(graph[i]) for i in range(n)}

    def once():
        find_rings(adj, [4, 5, 6])

    once()
    samples = [pn.timed(once) for _ in range(pn.REPEATS)]
    rows = [pn.row_of(
        info, "reviewer", "ice-Ic-4096", "rings", "SQQ", n, threads, samples,
        "sqq find_rings sizes 4-6; bond graph built outside the timer",
        version=getattr(sqq, "__version__", "0.5.8"),
    )]
    return rows


def scale_rows(info) -> list[dict]:
    """d-SEAMS 2.0 size sweep, the current-code side of the version figure."""
    import paper_stats as stats
    bench = ROOT / "bbdir-release/tests/bench_strong"
    threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
    rows = []
    for n_atoms in (1000, 2000, 4000, 8000):
        def once(n_atoms=n_atoms):
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(threads)
            env["OMP_WAIT_POLICY"] = "PASSIVE"
            done = subprocess.run(
                ["taskset", "-c", os.environ.get(
                    "PAPER_TASKSET", "1" if threads == 1 else stats.taskset_cores(threads)),
                 str(bench), str(n_atoms), "1"],
                cwd=ROOT, env=env, text=True, capture_output=True)
            if done.returncode != 0:
                raise RuntimeError((done.stderr or done.stdout or "bench_strong")[-400:])
            lines = [line for line in done.stdout.splitlines() if line and not line.startswith("#")]
            fields = dict(zip(lines[0].split(), lines[1].split()))
            return float(fields["total/ms"])

        once()
        samples = [once() for _ in range(pn.REPEATS)]
        rows.append(pn.row_of(
            info, "fig8", f"n-{n_atoms}", "pipeline", "d-SEAMS-2", n_atoms, threads, samples,
            "bench_strong total; d-SEAMS 2.0 on this CPU; v1 binary is a separate row",
        ))
    return rows


def v1_note(info) -> dict:
    binary = ROOT / "build/baseline-v1.0.1/driver-build/baseline_dseams1"
    return {
        "source_commit": info["source_commit"],
        "version": info["version"],
        "cpu_model": info["cpu_model"],
        "cpu_cores": info["cpu_cores"],
        "affinity": "1",
        "processes": 1,
        "family": "fig8",
        "phase": "v1",
        "operation": "cages",
        "code": "d-SEAMS-1",
        "molecules": 0,
        "threads": 1,
        "warmup": pn.WARMUP,
        "n": 0,
        "median_ms": "",
        "iqr_ms": "",
        "q1_ms": "",
        "q3_ms": "",
        "samples_ms": "",
        "detail": (
            "present" if binary.is_file() else
            "v1 driver is not built in this checkout; the 2.0 rows carry this CPU model"
        ),
    }


def macros(rows: list[dict], build: str) -> list[dict]:
    """Paper-macro keys. Threads stay a column, so the key names the arm."""
    named = []

    def add(key, row):
        if not row or not row.get("median_ms"):
            return
        named.append({
            "key": key,
            "value": row["median_ms"],
            "n": row["n"],
            "threads": row["threads"],
            "build": build,
            "cpu_model": row["cpu_model"],
        })

    def pick(family, phase, operation, code, threads):
        hits = [
            row for row in rows
            if row.get("family") == family and row.get("phase") == phase
            and row.get("operation") == operation and row.get("code") == code
            and str(row.get("threads")) == str(threads) and row.get("median_ms")
        ]
        return hits[0] if hits else None

    for threads, suffix in ((1, "t1"), (8, "t8")):
        add(f"HringDsInside.{suffix}", pick("rings", "ice-Ic-110592", "rings-inside", "d-SEAMS", threads))
        add(f"HringDsPre.{suffix}", pick("rings", "ice-Ic-110592", "rings-prebuilt", "d-SEAMS", threads))
        add(f"HringOvInside.{suffix}", pick("rings", "ice-Ic-110592", "rings-inside", "OVITO", threads))
        add(f"HringOvPre.{suffix}", pick("rings", "ice-Ic-110592", "rings-prebuilt", "OVITO", threads))
        add(f"HclassDs.{suffix}", pick("sitting-large", "ice", "classifier", "d-SEAMS", threads))
        add(f"HptmDefault.{suffix}", pick("sitting-large", "ice", "ptm-default", "OVITO", threads))
        add(f"HptmIce.{suffix}", pick("sitting-large", "ice", "ptm-ice", "OVITO", threads))
        add(f"Hira.{suffix}", pick("sitting-large", "ice", "ira", "IRA", threads))
        add(f"HfreudNeigh4096.{suffix}", pick("reviewer", "ice-Ic-4096", "neighbours", "freud", threads))
        add(f"HfreudNeigh110592.{suffix}", pick("reviewer", "ice-Ic-110592", "neighbours", "freud", threads))
        add(f"Hsqq.{suffix}", pick("reviewer", "ice-Ic-4096", "rings", "SQQ", threads))
        full = pick("carry-sweep", "nucleation", "rings-full", "d-SEAMS", threads)
        inc = pick("carry-sweep", "nucleation", "rings-carry", "d-SEAMS", threads)
        add(f"HcarryFull.{suffix}", full)
        add(f"HcarryInc.{suffix}", inc)
        if full and inc and full.get("median_ms") and inc.get("median_ms"):
            frac = float(inc["median_ms"]) / float(full["median_ms"])
            named.append({
                "key": f"HcarryFrac.{suffix}",
                "value": f"{frac:.6f}",
                "n": full["n"],
                "threads": threads,
                "build": build,
                "cpu_model": full["cpu_model"],
            })
        for n_atoms in (1000, 2000, 4000, 8000):
            add(f"Hfig8Ds2.n{n_atoms}.{suffix}",
                pick("fig8", f"n-{n_atoms}", "pipeline", "d-SEAMS-2", threads))
    return named


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/reviewer")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--core", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--worker", default="")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    pn.stamp_env()
    pn.prepare_libs(args.lib)
    info = pn.identity()
    if args.worker:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
        dispatch = {
            "rings": lambda: pn.run_rings(info, args.lib),
            "carry": lambda: pn.run_carry(info),
            "freud": lambda: freud_rows(info),
            "sqq": lambda: sqq_rows(info),
            "large": lambda: pn.run_large(info, args.core),
            "scale": lambda: scale_rows(info),
        }
        aff = os.environ.get("PAPER_TASKSET", "")
        for row in dispatch[args.worker]():
            if aff:
                row["affinity"] = aff
            print(json.dumps(row), flush=True)
        return 0
    out = args.out
    script = str(Path(__file__).resolve())
    rows: list[dict] = []
    for threads in (1, 8):
        if threads == 8:
            # Eight threads share the cores. Wait out a cycless point that
            # still owns one of them so the wide arm is not sharing a core.
            while True:
                listed = subprocess.run(["pgrep", "-af", "cycles_iter"], text=True, capture_output=True)
                if "cycles_iter" not in (listed.stdout or ""):
                    break
                time.sleep(20)
        for worker in ("rings", "freud", "sqq", "carry", "scale", "large"):
            text = pin(threads, [
                script, "--worker", worker, "--threads", str(threads),
                "--lib", args.lib, "--core", args.core,
            ])
            rows.extend(json.loads(line) for line in text.splitlines() if line.startswith("{"))
            print(f"reviewer {worker} threads {threads}", flush=True)
        ira = pin(threads, [
            str(ROOT / "scripts/paper_numbers.py"), "--worker", "ira",
            "--threads", str(threads), "--lib", args.lib, "--core", args.core,
        ])
        got = [json.loads(line) for line in ira.splitlines() if line.startswith("{")]
        for row in got:
            row["affinity"] = "1" if threads == 1 else "0-7"
        rows.extend(got)
        print(f"reviewer ira threads {threads}", flush=True)
    rows.append(v1_note(info))
    write_csv(out / "samples.csv", rows)
    build = info["source_commit"]
    table = macros(rows, build)
    dest = out / "MACROS.tsv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("key", "value", "n", "threads", "build", "cpu_model"),
            delimiter="\t")
        writer.writeheader()
        writer.writerows(table)
    manifest = out / "manifest.txt"
    manifest.write_text(
        "\n".join([
            f"source_commit {build}",
            f"version {info['version']}",
            f"cpu_model {info['cpu_model']}",
            f"cpu_cores {info['cpu_cores']}",
            "warmup 1",
            "n 8",
            "one release build",
            "host name is not recorded",
            "serial affinity is core 1",
            "eight-thread affinity is cores 0-7",
        ]) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {dest} rows {len(table)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
