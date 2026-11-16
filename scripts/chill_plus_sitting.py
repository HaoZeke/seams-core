#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""OVITO-rule CHILL+ of two d-SEAMS builds beside mdapy ChillPlus, one sitting.

Ice Ic at the literature sizes, one thread and eight. Every row is a
long-lived worker of scripts/lit_one_sitting.py: the d-SEAMS timer is
seams_chill_plus with the OVITO rule on the cutoff graph, the mdapy timer
is ChillPlus.compute(). Rows run one block at a time and every other
worker is stopped while a block runs.

source_tree is the tree of the commit that measured; build_tree is the
tree each library was built from. No commit hash is written.

    python3 scripts/chill_plus_sitting.py \\
        --before bbdir-before/src --before-tree "$(git rev-parse ea8b9c55^{tree})" \\
        --after bbdir-release/src --out benchmarks/chill-plus-1t
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import lit_one_sitting as lit  # noqa: E402
import paper_stats as stats  # noqa: E402

OPERATION = "chill_plus"
MIN_REPEATS = 6
BUILDS = ("before", "after", "mdapy")


def code_of(build: str) -> str:
    return "mdapy" if build == "mdapy" else "d-SEAMS"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--before", type=Path, required=True,
                        help="build directory holding libyodaLib.so and libseams.so")
    parser.add_argument("--before-tree", required=True, help="source tree of --before")
    parser.add_argument("--after", type=Path, default=ROOT / "bbdir-release" / "src")
    parser.add_argument("--after-tree", default="",
                        help="source tree of --after (default: the measuring tree)")
    parser.add_argument("--sizes", default=",".join(map(str, lit.SIZES)))
    parser.add_argument("--widths", default=",".join(map(str, lit.WIDTHS)))
    parser.add_argument("--warmup", type=int, default=lit.WARMUP)
    parser.add_argument("--repeats", type=int, default=lit.REPEATS)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < MIN_REPEATS:
        parser.error(f"--repeats must be at least {MIN_REPEATS}")
    sizes = tuple(int(v) for v in args.sizes.split(","))
    widths = tuple(int(v) for v in args.widths.split(","))
    tree = stats.source_tree(ROOT)
    trees = {"before": args.before_tree, "after": args.after_tree or tree, "mdapy": ""}
    libs = {"before": args.before.resolve(), "after": args.after.resolve()}
    for build, lib_dir in libs.items():
        for name in ("libyodaLib.so", "libseams.so"):
            if not (lib_dir / name).is_file():
                raise SystemExit(f"{build}: no {name} in {lib_dir}")
    model = stats.machine(ROOT)["cpu_model"]
    cores, affinity = stats.cores_now()
    state = stats.worktree_state(ROOT, ("src", "bindings", "scripts"))

    frames = Path(tempfile.mkdtemp(prefix="chill-sitting-"))
    lit.write_frames(frames, sizes)
    script = ROOT / "scripts" / "lit_one_sitting.py"
    workers: dict[tuple[str, int], lit.Worker] = {}
    samples: dict[tuple[str, int, int], list[float]] = {}
    checks: dict[tuple[str, int, int], str] = {}
    ntasks: dict[tuple[str, int, int], str] = {}
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        for build in BUILDS:
            for threads in widths:
                code = code_of(build)
                lib_dir = libs.get(build, libs["after"])
                env = lit.child_env(code, threads, lib_dir, lib_dir.parent / "bindings")
                env["LIT_FRAMES"] = str(frames)
                env["SEAMS_CORE_LIB"] = str(lib_dir / "libyodaLib.so")
                env["SEAMS_API_LIB"] = str(lib_dir / "libseams.so")
                print(f"start {build} t{threads}", flush=True)
                workers[(build, threads)] = lit.Worker(
                    code, threads, lit.taskset_cores(threads), env, script)
        with (args.out / "raw-chill.txt").open("w", encoding="utf-8") as handle:
            handle.write("\n".join([
                f"cpu_model={model}",
                f"cores={cores}",
                f"affinity={affinity}",
                f"source_tree={tree}",
                f"worktree={state}",
                *(f"build_tree {build}={trees[build]}" for build in ("before", "after")),
                f"warmup={args.warmup}",
                f"n={args.repeats}",
                f"cutoff_A={lit.CUTOFF}",
                "timer=d-SEAMS seams_chill_plus OVITO rule, cutoff graph; mdapy ChillPlus.compute()",
                "layout=one process per build and width; every other process is SIGSTOP'd",
            ]) + "\n")
            for (build, threads), worker in workers.items():
                handle.write(f"worker build={build} threads={threads} version={worker.version} "
                             f"omp_max={worker.omp_max} affinity={worker.affinity} "
                             f"env={worker.env}\n")
            for molecules in sizes:
                for threads in widths:
                    for build in BUILDS:
                        key = (build, threads, molecules)
                        samples[key] = []
                        print(f"{build} {molecules} t{threads}", flush=True)
                        worker = workers[(build, threads)]
                        worker.begin()
                        try:
                            for k in range(args.warmup + args.repeats):
                                payload = worker.call(OPERATION, molecules)
                                if k >= args.warmup:
                                    samples[key].append(float(payload["ms"]))
                                checks[key] = payload["check"]
                                ntasks[key] = str(payload.get("ntask", ""))
                                handle.write(
                                    f"round={k} molecules={molecules} build={build} "
                                    f"threads={threads} ms={payload['ms']:.6f} "
                                    f"check={payload['check']} ntask={payload.get('ntask', '')} "
                                    f"warmup={int(k < args.warmup)}\n")
                        finally:
                            worker.end()
                        handle.flush()
    finally:
        for worker in workers.values():
            worker.close()

    fields = [
        "build", "code", "version", "molecules", "threads", "median_ms", "min_ms", "max_ms",
        "q1_ms", "q3_ms", "iqr_ms", "samples_ms", "check", "agreement", "ntask",
        "source_tree", "build_tree", "worktree", "warmup", "n", "cpus", "env", "cpu_model",
    ]
    failures = []
    with (args.out / "chill-sitting.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for molecules in sizes:
            for threads in widths:
                ref = checks[("before", threads, molecules)]
                for build in BUILDS:
                    key = (build, threads, molecules)
                    if build == "mdapy":
                        ok, note = lit.agree(OPERATION, ref, checks[key])
                    else:
                        ok, note = checks[key] == ref, ""
                    if not ok:
                        failures.append(f"{build} {molecules} t{threads}: {checks[key]} "
                                        f"against {ref} {note}")
                    summary = stats.summarize(samples[key])
                    worker = workers[(build, threads)]
                    writer.writerow({
                        "build": build, "code": code_of(build), "version": worker.version,
                        "molecules": molecules, "threads": threads,
                        "median_ms": f"{summary['median_ms']:.6f}",
                        "min_ms": f"{min(samples[key]):.6f}",
                        "max_ms": f"{max(samples[key]):.6f}",
                        "q1_ms": f"{summary['q1_ms']:.6f}",
                        "q3_ms": f"{summary['q3_ms']:.6f}",
                        "iqr_ms": f"{summary['iqr_ms']:.6f}",
                        "samples_ms": stats.fmt_samples(samples[key]),
                        "check": checks[key], "agreement": "agree" if ok else "differ",
                        "ntask": ntasks[key], "source_tree": tree,
                        "build_tree": trees[build], "worktree": state,
                        "warmup": args.warmup, "n": summary["n"],
                        "cpus": lit.taskset_cores(threads), "env": worker.env,
                        "cpu_model": model,
                    })
                    print(f"{build:7} {molecules:7} t{threads} {summary['median_ms']:10.3f} ms  "
                          f"{checks[key]}", flush=True)
    for line in failures:
        print(line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
