#!/usr/bin/env python3
"""Ice Ic rings, the same work on both sides.

Two variants, each interleaved, eight timed calls after one warm-up:

  bonds_and_rings  neighbour list or CreateBonds, then the ring search
  rings            ring search only, graph already built

Ideal ice Ic, bond 2.75 A, cutoff 3.5 A, rings of size at most 6.
Serial is core 0. Eight threads are cores 0-7. cycless is cycles_iter,
which takes no thread count, so it is recorded only at one thread.

    python3 scripts/ice_ic_session.py --csv benchmarks/h2h-ice-ic-rings.csv
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gen_ice_ic import ice_ic  # noqa: E402

CUTOFF = 3.5
MAX_RING = 6
WARMUP = 1
REPEATS = 8
BENCH = ROOT / "bbdir-release" / "tests" / "bench_rings_h2h"
LIB = ROOT / "bbdir-release" / "src" / "libseams.so"


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def median(samples):
    ordered = sorted(samples)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def version_string() -> str:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
    out = subprocess.check_output(
        [
            sys.executable, "-c",
            "import ctypes,sys; lib=ctypes.CDLL(sys.argv[1]); "
            "lib.seams_version_string.restype=ctypes.c_char_p; "
            "print(lib.seams_version_string().decode())",
            str(LIB),
        ],
        env=env,
        text=True,
    )
    return out.strip()


def write_dump(path: Path, points, side: float) -> None:
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(points)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0.0 {side:.8f}",
        f"0.0 {side:.8f}",
        f"0.0 {side:.8f}",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(points, start=1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


class Pipe:
    def __init__(self, argv, env):
        self.proc = subprocess.Popen(
            argv,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )

    def ready(self) -> str:
        assert self.proc.stdout is not None
        line = self.proc.stdout.readline()
        if not line.startswith("ready "):
            raise RuntimeError(f"worker failed: {line!r}")
        return line[len("ready ") :].strip()

    def ask(self, word: str) -> float:
        assert self.proc.stdin is not None and self.proc.stdout is not None
        self.proc.stdin.write(word + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline().strip()
        if not line.startswith("ms "):
            raise RuntimeError(f"bad sample {line!r}")
        return float(line.split()[1])

    def close(self) -> None:
        if self.proc.stdin is not None:
            self.proc.stdin.write("quit\n")
            self.proc.stdin.flush()
        self.proc.wait(timeout=120)


def dseams_env(threads: int) -> dict:
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    return env


def ovito_env(threads: int) -> dict:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu"
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    return env


def start_dseams(traj: Path, threads: int) -> tuple[Pipe, dict]:
    cores = "0" if threads == 1 else "0-7"
    pipe = Pipe(
        [
            "taskset", "-c", cores, str(BENCH), "--repl", str(traj),
            "1", "1", f"{CUTOFF:.8f}", str(MAX_RING),
        ],
        dseams_env(threads),
    )
    text = pipe.ready()
    parts = text.split()
    hist = {}
    if "hist" in parts:
        blob = " ".join(parts[parts.index("hist") + 1 :])
        if blob != "none":
            for item in blob.split():
                size, count = item.split(":")
                hist[size] = int(count)
    return pipe, hist


def start_ovito(python: str, traj: Path, threads: int) -> tuple[Pipe, dict]:
    cores = "0" if threads == 1 else "0-7"
    pipe = Pipe(
        [
            "taskset", "-c", cores, python, str(Path(__file__).resolve()),
            "--ovito-worker", str(traj),
        ],
        ovito_env(threads),
    )
    payload = json.loads(pipe.ready())
    return pipe, payload


def ovito_worker(traj: str) -> None:
    import time

    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource
    import ovito

    text = Path(traj).read_text(encoding="ascii").splitlines()
    n = int(text[3])
    side = float(text[5].split()[1])
    pos = np.empty((n, 3), dtype=np.float64)
    for i, line in enumerate(text[9 : 9 + n]):
        parts = line.split()
        pos[i, 0] = float(parts[-3])
        pos[i, 1] = float(parts[-2])
        pos[i, 2] = float(parts[-1])

    def fresh():
        data = DataCollection()
        data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=pos)
        return data

    bond = fresh()
    bond.apply(CreateBondsModifier(cutoff=CUTOFF))
    probe = Pipeline(source=StaticSource(data=bond))
    probe.modifiers.append(FindRingsModifier(maximum_ring_size=MAX_RING))
    table = probe.compute().tables["ring-size-histogram"]
    hist = {
        str(int(size)): int(count)
        for size, count in zip(table["Ring Size"][...], table["Count"][...])
        if int(count) != 0
    }
    sys.stdout.write(
        "ready "
        + json.dumps(
            {
                "hist": hist,
                "version": ".".join(str(part) for part in ovito.version),
                "molecules": n,
            }
        )
        + "\n"
    )
    sys.stdout.flush()
    for line in sys.stdin:
        word = line.strip()
        if word.startswith("q"):
            break
        if word not in ("go", "all"):
            continue
        t0 = time.perf_counter()
        if word == "all":
            data = fresh()
            data.apply(CreateBondsModifier(cutoff=CUTOFF))
            data.apply(FindRingsModifier(minimum_ring_size=3, maximum_ring_size=MAX_RING))
        else:
            pipe = Pipeline(source=StaticSource(data=bond))
            pipe.modifiers.append(FindRingsModifier(maximum_ring_size=MAX_RING))
            pipe.compute()
        t1 = time.perf_counter()
        sys.stdout.write(f"ms {(t1 - t0) * 1000.0:.6f}\n")
        sys.stdout.flush()


def cycless_block(points, side: float) -> list[dict]:
    """Serial walk. cycles_iter has no thread argument."""
    import time

    import networkx as nx
    import numpy as np
    from cycless.cycles import cycles_iter

    held = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {0})
    try:
        pos = np.asarray(points, dtype=np.float64)
        n = len(pos)

        def graph_of():
            ncell = max(1, int(side / CUTOFF))
            cell = side / ncell
            buckets = [[] for _ in range(ncell ** 3)]

            def cell_of(x):
                c = int(x / cell)
                return 0 if c < 0 else ncell - 1 if c >= ncell else c

            for i, p in enumerate(pos):
                buckets[(cell_of(p[0]) * ncell + cell_of(p[1])) * ncell + cell_of(p[2])].append(i)
            edges = []
            cut2 = CUTOFF * CUTOFF
            for i, pi in enumerate(pos):
                cx, cy, cz = cell_of(pi[0]), cell_of(pi[1]), cell_of(pi[2])
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            key = (
                                ((cx + dx) % ncell) * ncell + ((cy + dy) % ncell)
                            ) * ncell + ((cz + dz) % ncell)
                            for j in buckets[key]:
                                if j <= i:
                                    continue
                                d = pi - pos[j]
                                d -= side * np.round(d / side)
                                if float(d @ d) <= cut2:
                                    edges.append((i, j))
            graph = nx.Graph()
            graph.add_nodes_from(range(n))
            graph.add_edges_from(edges)
            return graph

        built = graph_of()

        def count(graph):
            return sum(1 for _ in cycles_iter(graph, MAX_RING))

        rows = []
        for variant, fn in (
            ("bonds_and_rings", lambda: count(graph_of())),
            ("rings", lambda: count(built)),
        ):
            samples = []
            census = None
            for i in range(WARMUP + REPEATS):
                t0 = time.perf_counter()
                census = fn()
                dt = (time.perf_counter() - t0) * 1000.0
                if i >= WARMUP:
                    samples.append(dt)
            rows.append(
                {
                    "variant": variant,
                    "code": "cycless",
                    "version": "0.7",
                    "threads": 1,
                    "cpus": "0",
                    "median_ms": median(samples),
                    "min_ms": min(samples),
                    "samples_ms": samples,
                    "census": census,
                    "note": "cycles_iter has no thread argument; the walk is serial",
                }
            )
        return rows
    finally:
        os.sched_setaffinity(0, held)


def paired(dseams: Pipe, ovito: Pipe, word: str) -> tuple[list, list]:
    left, right = [], []
    for i in range(WARMUP + REPEATS):
        order = ("d", "o") if i % 2 == 0 else ("o", "d")
        got = {}
        for arm in order:
            got[arm] = dseams.ask(word) if arm == "d" else ovito.ask(word)
        if i >= WARMUP:
            left.append(got["d"])
            right.append(got["o"])
    return left, right


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=False)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--ovito-worker", metavar="TRAJ")
    args = ap.parse_args()
    if args.ovito_worker:
        ovito_worker(args.ovito_worker)
        return
    if not args.csv:
        raise SystemExit("--csv is required")
    if not BENCH.is_file():
        raise SystemExit(f"missing {BENCH.name}")

    seams_version = version_string()
    cores = os.cpu_count() or 1
    model = cpu_model()
    rows = []
    work = ROOT / "bbdir-release" / "h2h-ice-ic"
    work.mkdir(parents=True, exist_ok=True)
    for reps in (5, 11, 24):
        points, side = ice_ic(reps)
        traj = work / f"ice-ic-{len(points)}.lammpstrj"
        write_dump(traj, points, side)
        for threads in (1, 8):
            dseams, dhist = start_dseams(traj, threads)
            ovito, payload = start_ovito(args.python, traj, threads)
            ohist = {str(k): int(v) for k, v in payload["hist"].items()}
            if {k: v for k, v in dhist.items() if int(k) >= 4} != ohist:
                raise SystemExit(f"census d-SEAMS {dhist} OVITO {ohist}")
            try:
                for variant, word in (("rings", "go"), ("bonds_and_rings", "all")):
                    left, right = paired(dseams, ovito, word)
                    for code, version, samples in (
                        ("d-SEAMS", seams_version, left),
                        ("OVITO", payload["version"], right),
                    ):
                        rows.append(
                            {
                                "phase": "Ic",
                                "molecules": len(points),
                                "variant": variant,
                                "code": code,
                                "version": version,
                                "threads": threads,
                                "cores": cores,
                                "cpus": "0" if threads == 1 else "0-7",
                                "n": len(samples),
                                "median_ms": f"{median(samples):.4f}",
                                "min_ms": f"{min(samples):.4f}",
                                "samples_ms": " ".join(f"{v:.4f}" for v in samples),
                                "census": json.dumps(payload["hist"], sort_keys=True),
                                "cpu_model": model,
                                "note": "",
                            }
                        )
                    print(
                        f"Ic {len(points)} {variant} threads {threads} "
                        f"dseams {median(left):.4f} ovito {median(right):.4f}",
                        flush=True,
                    )
            finally:
                dseams.close()
                ovito.close()
        if len(points) == 1000:
            for block in cycless_block(points, side):
                rows.append(
                    {
                        "phase": "Ic",
                        "molecules": 1000,
                        "variant": block["variant"],
                        "code": "cycless",
                        "version": block["version"],
                        "threads": block["threads"],
                        "cores": cores,
                        "cpus": block["cpus"],
                        "n": len(block["samples_ms"]),
                        "median_ms": f"{block['median_ms']:.4f}",
                        "min_ms": f"{block['min_ms']:.4f}",
                        "samples_ms": " ".join(f"{v:.4f}" for v in block["samples_ms"]),
                        "census": str(block["census"]),
                        "cpu_model": model,
                        "note": block["note"],
                    }
                )
                print(
                    f"Ic 1000 {block['variant']} cycless threads 1 "
                    f"median {block['median_ms']:.4f}",
                    flush=True,
                )

    import csv

    fields = [
        "phase", "molecules", "variant", "code", "version", "threads", "cores",
        "cpus", "n", "median_ms", "min_ms", "samples_ms", "census", "cpu_model",
        "note",
    ]
    out = Path(args.csv)
    with out.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out.name} version {seams_version} cores {cores}", flush=True)


if __name__ == "__main__":
    main()
