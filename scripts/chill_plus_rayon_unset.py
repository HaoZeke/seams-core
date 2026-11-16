#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""One-thread OVITO-rule CHILL+ under taskset -c 0 with RAYON_NUM_THREADS unset.

The case the sittings no longer measure: OMP_NUM_THREADS=1 and one cpu,
but no rayon width, so a library that widens its mask past the launch
cpu lets the linkcell pool start a worker per allowed cpu. Each round
records the time, the worker's mask after the call and its thread count
(/proc/<pid>/task). This is a regression record, not a benchmark row.

    python3 scripts/chill_plus_rayon_unset.py \\
        --before bbdir-before/src --before-tree "$(git rev-parse ca7fce9d^{tree})" \\
        --after bbdir-release/src --out benchmarks/chill-plus-1t/affinity
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import lit_one_sitting as lit  # noqa: E402
import paper_stats as stats  # noqa: E402

OPERATION = "chill_plus"
MOLECULES = 110592
CPUS = "0"
MIN_REPEATS = 6


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--before-tree", required=True)
    parser.add_argument("--after", type=Path, default=ROOT / "bbdir-release" / "src")
    parser.add_argument("--warmup", type=int, default=lit.WARMUP)
    parser.add_argument("--repeats", type=int, default=lit.REPEATS)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < MIN_REPEATS:
        parser.error(f"--repeats must be at least {MIN_REPEATS}")
    tree = stats.source_tree(ROOT)
    trees = {"before": args.before_tree, "after": tree}
    libs = {"before": args.before.resolve(), "after": args.after.resolve()}
    cores, affinity = stats.cores_now()
    frames = Path(tempfile.mkdtemp(prefix="chill-rayon-unset-"))
    lit.write_frames(frames, (MOLECULES,))
    script = ROOT / "scripts" / "lit_one_sitting.py"
    args.out.mkdir(parents=True, exist_ok=True)
    medians = {}
    with (args.out / "rayon-unset.txt").open("w", encoding="utf-8") as handle:
        handle.write("\n".join([
            f"cpu_model={stats.machine(ROOT)['cpu_model']}",
            f"cores={cores}",
            f"affinity={affinity}",
            f"source_tree={tree}",
            f"worktree={stats.worktree_state(ROOT, ('src', 'bindings', 'scripts'))}",
            *(f"build_tree {build}={trees[build]}" for build in libs),
            f"warmup={args.warmup}",
            f"n={args.repeats}",
            f"taskset={CPUS}",
            "timer=d-SEAMS seams_chill_plus OVITO rule, cutoff graph",
            "threads=entries of /proc/<worker>/task after the call",
        ]) + "\n")
        for build, lib_dir in libs.items():
            env = lit.child_env("d-SEAMS", 1, lib_dir, lib_dir.parent / "bindings")
            env.pop("RAYON_NUM_THREADS", None)
            env["LIT_FRAMES"] = str(frames)
            env["SEAMS_CORE_LIB"] = str(lib_dir / "libyodaLib.so")
            env["SEAMS_API_LIB"] = str(lib_dir / "libseams.so")
            worker = lit.Worker("d-SEAMS", 1, CPUS, env, script, confine=False)
            try:
                handle.write(f"worker build={build} version={worker.version} "
                             f"omp_max={worker.omp_max} env={worker.env}\n")
                worker.begin()
                samples = []
                for k in range(args.warmup + args.repeats):
                    payload = worker.call(OPERATION, MOLECULES)
                    threads = len(os.listdir(f"/proc/{worker.proc.pid}/task"))
                    if k >= args.warmup:
                        samples.append(float(payload["ms"]))
                    handle.write(
                        f"round={k} molecules={MOLECULES} build={build} "
                        f"ms={payload['ms']:.6f} check={payload['check']} "
                        f"cpus_effective={payload['cpus_effective']} cpus={payload['cpus']} "
                        f"escaped={int(payload['escaped'])} threads={threads} "
                        f"warmup={int(k < args.warmup)}\n")
                worker.end()
            finally:
                worker.close()
            summary = stats.summarize(samples)
            medians[build] = summary
            handle.write(f"summary build={build} median_ms={summary['median_ms']:.6f} "
                         f"q1_ms={summary['q1_ms']:.6f} q3_ms={summary['q3_ms']:.6f} "
                         f"n={summary['n']}\n")
            handle.flush()
    for build, summary in medians.items():
        print(f"{build:7} {MOLECULES} t1 rayon unset {summary['median_ms']:10.3f} ms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
