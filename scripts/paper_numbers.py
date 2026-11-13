#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Every paper timing from one release build.

One warm-up, then eight calls. When two codes are compared, the calls
alternate. The median of an even count is the mean of the two central
samples. Quartiles are type 7. Every row names the commit, the version
and the core count. The host name is not read.

This is the body of the Snakemake target ``paper-numbers``. A later
release commit reruns the same target; this script does not special-case
that commit.
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
sys.path.insert(0, str(ROOT / "scripts"))

WARMUP = 1
REPEATS = 8
CUTOFF = 3.5

FIELDS = [
    "source_commit", "version", "cpu_model", "cpu_cores", "affinity", "processes",
    "family", "phase", "operation", "code", "molecules", "threads",
    "warmup", "n", "median_ms", "iqr_ms", "q1_ms", "q3_ms", "samples_ms", "detail",
]


def stamp_env() -> None:
    os.environ["LJOS_SEAT"] = "dseams-paper"
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("MPLBACKEND", "Agg")
    extra = os.environ.get("OVITO_EXTRA", "").strip() or str(ROOT / ".cache" / "ovito-extra")
    os.environ["OVITO_EXTRA"] = extra


def identity():
    import paper_stats as stats
    info = stats.machine(ROOT)
    version = subprocess.check_output(
        [str(ROOT / "bbdir-release/seams"), "--version"], text=True).strip()
    info["version"] = version
    return info


def row_of(info, family, phase, operation, code, molecules, threads, samples, detail, processes=1, version=None):
    import paper_stats as stats
    summary = stats.summarize(samples)
    return {
        "source_commit": info["source_commit"],
        "version": version or info["version"],
        "cpu_model": info["cpu_model"],
        "cpu_cores": info["cpu_cores"],
        "affinity": stats.taskset_cores(int(threads)),
        "processes": processes,
        "family": family,
        "phase": phase,
        "operation": operation,
        "code": code,
        "molecules": molecules,
        "threads": threads,
        "warmup": WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.6f}",
        "iqr_ms": f"{summary['iqr_ms']:.6f}",
        "q1_ms": f"{summary['q1_ms']:.6f}",
        "q3_ms": f"{summary['q3_ms']:.6f}",
        "samples_ms": stats.fmt_samples(samples),
        "detail": detail,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    commits = {row["source_commit"] for row in rows}
    if len(commits) != 1:
        raise SystemExit(f"{path.name} mixed commits {commits}")
    print(f"wrote {path} rows {len(rows)} commit {commits.pop()[:12]}", flush=True)


def timed(fn):
    t0 = time.perf_counter()
    fn()
    return (time.perf_counter() - t0) * 1e3


def alternate(calls: dict, repeats: int = REPEATS):
    names = list(calls)
    buckets = {name: [] for name in names}
    for index in range(WARMUP + repeats):
        order = names[index % len(names):] + names[:index % len(names)]
        for name in order:
            sample = timed(calls[name])
            if index >= WARMUP:
                buckets[name].append(sample)
    return buckets


def pin_and_run(threads: int, script_args: list[str]) -> str:
    import paper_stats as stats
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OVITO_THREAD_COUNT"] = str(threads)
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    cmd = ["taskset", "-c", stats.taskset_cores(threads), sys.executable, *script_args]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise SystemExit("\n".join(tail[-30:]) if tail else f"exit {done.returncode}")
    return done.stdout


def prepare_libs(lib: str) -> None:
    import paper_sitting as sit
    build = Path(lib).resolve().parents[1]
    extra = os.environ["OVITO_EXTRA"]
    parts = [sit.library_path(lib), str(build / "subprojects/chemfiles-prefix/lib"), extra]
    prior = [
        item for item in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep)
        if item and not item.startswith("/usr/lib")
    ]
    os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(parts + prior)


def run_h2h(info) -> list[dict]:
    """d-SEAMS, freud and OVITO on the ice cells, one warm-up and eight calls."""
    env = os.environ.copy()
    env["OVITO_THREAD_COUNT"] = "8"
    done = subprocess.run(
        [sys.executable, str(ROOT / "scripts/h2h_session.py")],
        cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise SystemExit("h2h session failed\n" + "\n".join(tail[-30:]))
    rows = []
    for line in done.stdout.splitlines():
        if "samples_ms" not in line:
            continue
        parts = line.split()
        # arm CODE op OP threads T molecules N samples_ms ...
        code, op, threads, molecules = parts[1], parts[3], int(parts[5]), int(parts[7])
        samples = [float(v) for v in parts[9:]]
        if len(samples) < 8:
            raise SystemExit(f"h2h {code} {op} has {len(samples)} samples")
        rows.append(row_of(
            info, "h2h-summary", "ice-Ih", op, code, molecules, threads, samples,
            "h2h_session one warm-up then eight calls",
            version=info["version"] if code == "dseams" else code,
        ))
    if not rows:
        raise SystemExit("h2h session wrote no sample lines")
    return rows


def run_carry(info) -> list[dict]:
    """Full rebuild against the carried update, on the nucleation dump."""
    dump = ROOT / "data/hop/nucleation.lammpstrj"
    bench = ROOT / "bbdir-release/tests/bench_trajectory_incremental"
    if not dump.is_file() or not bench.is_file():
        raise SystemExit("carry sweep needs the nucleation dump and bench_trajectory_incremental")
    rows = []
    import paper_stats as stats
    for threads in (1, 8):
        full, inc = [], []
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_PROC_BIND"] = "close"
        env["OMP_PLACES"] = "cores"

        def one():
            done = subprocess.run(
                ["taskset", "-c", stats.taskset_cores(threads), str(bench),
                 str(dump), "148", "1", "6", "2.0", "cutoff"],
                cwd=ROOT, env=env, text=True, capture_output=True)
            if done.returncode != 0:
                raise RuntimeError((done.stderr or done.stdout or "carry failed")[-400:])
            steady = [line for line in done.stdout.splitlines() if line.startswith("# steady ")]
            if not steady:
                raise RuntimeError(done.stdout[-400:])
            parts = steady[-1].split()
            return float(parts[3]), float(parts[5])

        one()
        for _ in range(REPEATS):
            ring_full, ring_inc = one()
            full.append(ring_full)
            inc.append(ring_inc)
        rows.append(row_of(
            info, "carry-sweep", "nucleation", "rings-full", "d-SEAMS", 0, threads, full,
            "bench_trajectory_incremental ringFull 148 frames",
        ))
        rows.append(row_of(
            info, "carry-sweep", "nucleation", "rings-carry", "d-SEAMS", 0, threads, inc,
            "bench_trajectory_incremental ringInc 148 frames",
        ))
        print(f"carry threads {threads}", flush=True)
    return rows


def run_rings(info, lib: str) -> list[dict]:
    """Ice Ic 4096 and 110592, and structure H at 136, bonds prebuilt and inside."""
    import ctypes
    import numpy as np
    import paper_ptm_templates as ptm
    import paper_sitting as sit
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource
    import ovito

    lang = sit.open_lang(lib)
    ovito_version = ".".join(str(part) for part in tuple(ovito.version)[:3])
    frames = []
    for molecules, cells in ((4096, 8), (110592, 24)):
        pos, lattice = ptm.ice_ic(cells)
        side = float(lattice[0][0])
        box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
        frames.append((f"ice-Ic-{molecules}", np.ascontiguousarray(pos), box, np.asarray(lattice)))
    rows = []
    for phase, pos, box, lattice in frames:
        n = len(pos)

        def fresh():
            return sit.frame_from(lang, pos, box)[0]

        def six_on(handle):
            counts = (ctypes.c_int * 7)()
            sit.fail_lang(lang, lang.seams_frame_rings(
                handle, CUTOFF, 1, b"cutoff", 6, counts, 7), "rings")
            return int(counts[6])

        def dseams(inside: bool):
            handle = fresh()
            try:
                if not inside:
                    lang.seams_frame_neighbours(handle, CUTOFF, 1, b"cutoff")
                t0 = time.perf_counter()
                six_on(handle)
                return (time.perf_counter() - t0) * 1e3
            finally:
                lang.seams_frame_free(handle)

        def ovito(inside: bool):
            data = DataCollection()
            data.create_cell(np.asarray(lattice, dtype=float), pbc=(True, True, True))
            parts = data.create_particles(count=n)
            parts.create_property("Position", data=np.ascontiguousarray(pos))
            if not inside:
                data.apply(CreateBondsModifier(cutoff=CUTOFF))
            pipe = Pipeline(source=StaticSource(data=data))
            if inside:
                pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
            pipe.modifiers.append(FindRingsModifier(minimum_ring_size=3, maximum_ring_size=6))
            t0 = time.perf_counter()
            pipe.compute()
            return (time.perf_counter() - t0) * 1e3

        calls = {
            "d-SEAMS prebuilt": lambda: dseams(False),
            "d-SEAMS inside": lambda: dseams(True),
            "OVITO prebuilt": lambda: ovito(False),
            "OVITO inside": lambda: ovito(True),
        }
        names = list(calls)
        buckets = {name: [] for name in names}
        for index in range(WARMUP + REPEATS):
            order = names[index % len(names):] + names[:index % len(names)]
            for name in order:
                sample = calls[name]()
                if index >= WARMUP:
                    buckets[name].append(sample)
        for name, samples in buckets.items():
            code, bonds = name.split()
            rows.append(row_of(
                info, "rings", phase, f"rings-{bonds}", code, n,
                int(os.environ.get("OMP_NUM_THREADS", "1")), samples,
                "alternating calls; bonds prebuilt or built inside the sample",
                version=info["version"] if code == "d-SEAMS" else ovito_version,
            ))
        print(f"rings {phase}", flush=True)
    return rows


def run_symmetry(info) -> list[dict]:
    """Classifier at 110592, the shell match bench_v2 already times."""
    import paper_stats as stats
    import paper_symmetry_sitting as sym
    rows = []
    binary = ROOT / "bbdir-release/tests/bench_v2"
    for threads in (1, 8):
        samples = []
        env = sym.bench_env(threads)
        cmd = [
            "taskset", "-c", stats.taskset_cores(threads),
            str(binary), "symmetry", "110592", str(threads), "1",
        ]
        subprocess.run(cmd, cwd=ROOT, env=env, check=True, text=True, capture_output=True)
        for _ in range(REPEATS):
            done = subprocess.run(cmd, cwd=ROOT, env=env, check=True, text=True, capture_output=True)
            line = [row for row in done.stdout.splitlines() if row.startswith("mode=")][-1]
            samples.append(float(line.split("samples_ms=", 1)[1].split(",")[0]))
        rows.append(row_of(
            info, "symmetry", "ice-Ic-110592", "classifier", "d-SEAMS",
            110592, threads, samples,
            "bench_v2 symmetry; neighbour list and diamondFromCsr inside the call; one process",
        ))
        print(f"symmetry threads {threads}", flush=True)
    return rows


def run_large(info, core: str) -> list[dict]:
    """884736 classifier, PTM and IRA on one build, into sitting-large rows."""
    import paper_ptm_templates as ptm
    pos, lattice = ptm.ice_ic(48)
    n = len(pos)
    if n != 884736:
        raise SystemExit(f"ice Ic 48^3 is {n}, need 884736")
    lib = ptm.open_core(core)
    os.environ["OVITO_THREAD_COUNT"] = os.environ.get("OMP_NUM_THREADS", "1")

    def classify():
        ptm.classify(lib, pos, lattice)

    def ptm_kind(kind):
        def run():
            ptm.ovito_labels(pos, lattice, kind)
        return run

    calls = {
        "classifier": classify,
        "ptm-default": ptm_kind("default"),
        "ptm-ice": ptm_kind("ice"),
    }
    buckets = alternate(calls)
    threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
    rows = []
    for operation, samples in buckets.items():
        code = "d-SEAMS" if operation == "classifier" else "OVITO"
        version = info["version"] if code == "d-SEAMS" else "3.17.2"
        rows.append(row_of(
            info, "sitting-large", "ice", operation, code, n, threads, samples,
            "884736 ideal ice Ic; one process; alternating with the other labels",
            version=version,
        ))
    ira = ira_samples(pos, lattice, threads)
    if ira is not None:
        samples, version, processes, detail = ira
        rows.append(row_of(
            info, "sitting-large", "ice", "ira", "IRA", n, threads, samples, detail,
            processes=processes, version=version,
        ))
    print(f"large threads {threads} rows {len(rows)}", flush=True)
    return rows


def ira_samples(pos, lattice, threads: int):
    """One process over every oxygen. No row is written if IRA is absent."""
    import numpy as np
    try:
        import ira_mod
    except ImportError:
        return None
    import vesin
    version = str(getattr(ira_mod, "version", "ira"))
    points = np.ascontiguousarray(pos, dtype=np.float64)
    n = len(points)
    side = float(np.asarray(lattice).reshape(-1)[0])
    box = np.eye(3) * side
    finder = vesin.NeighborList(cutoff=5.0, full_list=True, n_threads=1)
    left, dist, vec = finder.compute(
        points=points, box=box, periodic=True, quantities="idD", copy=True)
    order = np.argsort(left, kind="mergesort")
    left, dist, vec = left[order], dist[order], vec[order]
    cuts = np.searchsorted(left, np.arange(n + 1))
    shells = np.full((n, 17, 3), np.nan)
    for i in range(n):
        start, stop = int(cuts[i]), int(cuts[i + 1])
        if stop - start < 16:
            continue
        block = dist[start:stop]
        pick = np.argpartition(block, 15)[:16]
        pick = pick[np.argsort(block[pick])]
        shells[i, 0] = points[i]
        shells[i, 1:] = points[i] + vec[start:stop][pick]
    types = np.ones(17, dtype=np.int32)
    templates = [shells[i] for i in range(min(3, n)) if np.isfinite(shells[i, 0, 0])]
    if not templates:
        return None
    ira = ira_mod.IRA()

    def one():
        for shell in shells:
            if not np.isfinite(shell[0, 0]):
                continue
            best = 1e9
            for templ in templates:
                _r, _t, _p, dist_h = ira.match(17, types, templ, 17, types, shell, 1.5)
                best = min(best, float(dist_h))
            if best > 1e8:
                raise RuntimeError("IRA match failed")

    one()
    samples = [timed(one) for _ in range(REPEATS)]
    return (
        samples, version, 1,
        f"one process; shells built once; OMP_NUM_THREADS={threads}; no process pool",
    )


def run_figs(info) -> tuple[list[dict], list[dict]]:
    """Figs 8 and 9: size sweep and the 65536-atom pipeline, 1 and 8 threads."""
    import paper_stats as stats
    bench = ROOT / "bbdir-release/tests/bench_strong"
    scale_sizes = (1000, 2000, 4000, 8000, 16000, 32000)
    strong_n = 65536

    def once(n_atoms: int, threads: int) -> dict[str, float]:
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)
        env.pop("OMP_PROC_BIND", None)
        env.pop("OMP_PLACES", None)
        done = subprocess.run(
            ["taskset", "-c", stats.taskset_cores(threads), str(bench), str(n_atoms), "1"],
            cwd=ROOT, env=env, text=True, capture_output=True)
        if done.returncode != 0:
            raise RuntimeError((done.stderr or done.stdout or "bench_strong failed")[-400:])
        lines = [line for line in done.stdout.splitlines() if line and not line.startswith("#")]
        header = lines[0].split()
        values = lines[1].split()
        fields = dict(zip(header, values))
        return {
            "cutoff-neighbours": float(fields["neigh/ms"]),
            "steinhardt": float(fields["ql/ms"]),
            "index-neighbours": float(fields["index/ms"]),
            "rings": float(fields["rings/ms"]),
            "pipeline": float(fields["total/ms"]),
        }

    scale_rows, strong_rows = [], []
    for threads in (1, 8):
        for n_atoms in scale_sizes:
            once(n_atoms, threads)
            buckets = {name: [] for name in (
                "cutoff-neighbours", "index-neighbours", "rings")}
            for _ in range(REPEATS):
                sample = once(n_atoms, threads)
                for name in buckets:
                    buckets[name].append(sample[name])
            for operation, samples in buckets.items():
                scale_rows.append(row_of(
                    info, "fig-scaling", f"scale-{n_atoms}", operation, "d-SEAMS",
                    n_atoms, threads, samples,
                    "bench_strong reps=1; one warm-up then eight calls; this host",
                ))
            print(f"fig8 n={n_atoms} threads={threads}", flush=True)
        once(strong_n, threads)
        buckets = {name: [] for name in (
            "pipeline", "rings", "steinhardt", "cutoff-neighbours", "index-neighbours")}
        for _ in range(REPEATS):
            sample = once(strong_n, threads)
            for name in buckets:
                buckets[name].append(sample[name])
        for operation, samples in buckets.items():
            strong_rows.append(row_of(
                info, "fig-strong", "strong-65536", operation, "d-SEAMS",
                strong_n, threads, samples,
                "bench_strong reps=1; one warm-up then eight calls; this host",
            ))
        print(f"fig9 threads={threads}", flush=True)
    return scale_rows, strong_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-numbers")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--core", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--worker", default="")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    stamp_env()
    prepare_libs(args.lib)
    if args.worker == "rings":
        info = identity()
        for item in run_rings(info, args.lib):
            import json
            print(json.dumps(item), flush=True)
        return 0
    if args.worker == "large":
        info = identity()
        for item in run_large(info, args.core):
            import json
            print(json.dumps(item), flush=True)
        return 0
    info = identity()
    out = args.out
    h2h = run_h2h(info)
    write_csv(out / "h2h-summary.csv", h2h)
    carry = run_carry(info)
    write_csv(out / "carry-sweep.csv", carry)
    ring_rows = []
    for threads in (1, 8):
        text = pin_and_run(threads, [
            str(Path(__file__).resolve()), "--worker", "rings",
            "--threads", str(threads), "--lib", args.lib,
        ])
        import json
        ring_rows.extend(json.loads(line) for line in text.splitlines() if line.startswith("{"))
    write_csv(out / "h2h-ice-ic-rings.csv", ring_rows)
    log = out / "ice-ic-session.txt"
    log.write_text(
        "\n".join([
            f"version {info['version']}",
            f"source_commit {info['source_commit']}",
            f"cpu_model {info['cpu_model']}",
            f"cpu_cores {info['cpu_cores']}",
            "warmup 1",
            "n 8",
            "protocol alternating calls, one process",
            "this log replaces the 2.13.0 ice-ic session",
            *[
                f"{row['phase']} {row['code']} {row['operation']} threads {row['threads']} "
                f"median_ms={row['median_ms']} samples_ms={row['samples_ms']}"
                for row in ring_rows
            ],
        ]) + "\n",
        encoding="utf-8",
    )
    write_csv(out / "symmetry-rows.csv", run_symmetry(info))
    large = []
    for threads in (1, 8):
        text = pin_and_run(threads, [
            str(Path(__file__).resolve()), "--worker", "large",
            "--threads", str(threads), "--lib", args.lib, "--core", args.core,
        ])
        import json
        large.extend(json.loads(line) for line in text.splitlines() if line.startswith("{"))
    write_csv(out / "sitting-large.csv", large)
    scale, strong = run_figs(info)
    write_csv(out / "fig-scaling.csv", scale)
    write_csv(out / "fig-strong.csv", strong)
    manifest = out / "manifest.txt"
    manifest.write_text(
        "\n".join([
            f"source_commit {info['source_commit']}",
            f"version {info['version']}",
            f"cpu_model {info['cpu_model']}",
            f"cpu_cores {info['cpu_cores']}",
            "warmup 1",
            "n 8",
            "one release build",
            "host name is not recorded",
        ]) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {manifest}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
