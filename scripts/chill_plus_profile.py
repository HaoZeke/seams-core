#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Where one-thread OVITO-rule CHILL+ spends its time, from perf samples.

Times seams_chill_plus (OVITO rule, cutoff graph) on 110592-molecule ice
Ic under `perf record`, pinned to core 0 with OpenMP and rayon at one
thread, then sums the self time of each symbol into the neighbour build,
the q3 harmonics, the shell sort and classification, allocation, and the
API copy. Shares are scaled to the one-thread 110592 median of the same
build in chill-sitting.csv (or to --ms).

    python3 scripts/chill_plus_profile.py --api bbdir-release/src/libseams.so \\
        --sitting benchmarks/chill-plus-1t/chill-sitting.csv --build after \\
        --out benchmarks/chill-plus-1t/profile/after.txt
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

MOLECULES = 110592
CATEGORIES = (
    ("neighbour build",
     r"linkcell|minimage|vesin|cutoffRows|rowsByID|pairedRows|orthoFullCsr|sortRow|neighListO|^round$"),
    ("q3 harmonics", r"realYlm3|realL3Dot|addRealYlm3"),
    ("shell sort, deltas, c_ij and classification",
     r"classifyBondsLeft|__adjust_heap|__sort|__introsort|__insertion|fracDelta|"
     r"chillPlusOvitoFromCsr|assignIceTypePlus|cutoffShellCount|classifyBond|correlationOf"),
    ("allocation and free",
     r"malloc|_int_free|^free$|cfree|operator new|operator delete|unlink_chunk|"
     r"memset|memcpy|memmove|release_rows"),
    ("API copy and labels", r"kernel_cloud|write_ice_labels|normalize_types|seams_chill_plus"),
)


def drive(api_path: str, reps: int) -> None:
    import ctypes
    import numpy as np
    import lit_one_sitting as lit
    api = ctypes.CDLL(api_path)
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    vp = ctypes.c_void_p
    api.seams_system_create.restype = vp
    api.seams_system_set_ortho.argtypes = [vp, dp, ctypes.c_int] + [ctypes.c_double] * 6 + [ip, ip]
    api.seams_chill_plus.argtypes = [vp, ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_int, ip, ip]
    pos, _box9, _cell, side = lit.ice_ic(MOLECULES)
    pos = np.ascontiguousarray(pos)
    n = len(pos)
    handle = api.seams_system_create()
    if api.seams_system_set_ortho(handle, pos.ctypes.data_as(dp), n, 0.0, 0.0, 0.0,
                                  side, side, side, None, None) != 0:
        raise SystemExit("seams_system_set_ortho failed")
    labels = np.zeros(n, dtype=np.int32)
    counts = np.zeros(9, dtype=np.int32)
    samples = []
    for k in range(3 + reps):
        t0 = time.perf_counter()
        api.seams_chill_plus(handle, 3.5, 0, 0, 1, labels.ctypes.data_as(ip),
                             counts.ctypes.data_as(ip))
        if k >= 3:
            samples.append((time.perf_counter() - t0) * 1e3)
    samples.sort()
    print(f"under perf: median {samples[len(samples) // 2]:.3f} ms, counts {counts.tolist()}, "
          f"ntask {len(os.listdir('/proc/self/task'))}")


def group(report: str, ms: float) -> list[str]:
    totals = {name: 0.0 for name, _ in CATEGORIES}
    other = driver = 0.0
    for line in report.splitlines():
        m = re.match(r"\s*([\d.]+)%\s+(\S+)\s+\[.\]\s+(.*)", line)
        if not m:
            continue
        pct, dso = float(m.group(1)), m.group(2)
        sym = re.sub(r"^\(anonymous namespace\)::", "", m.group(3).strip())
        for name, pattern in CATEGORIES:
            if re.search(pattern, sym):
                totals[name] += pct
                break
        else:
            if dso.startswith("python"):
                driver += pct
            else:
                other += pct
    inside = sum(totals.values()) + other
    lines = [f"{name:46s} {100 * totals[name] / inside:5.1f}%  {ms * totals[name] / inside:6.1f} ms"
             for name, _ in CATEGORIES]
    lines.append(f"{'other':46s} {100 * other / inside:5.1f}%  {ms * other / inside:6.1f} ms")
    lines.append(f"python driver: {driver:.1f}% of samples, left out of the shares")
    return lines


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "--drive":
        drive(sys.argv[2], int(sys.argv[3]))
        return 0
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--api", type=Path, required=True, help="libseams.so to profile")
    parser.add_argument("--ms", type=float, help="sitting median of the call")
    parser.add_argument("--sitting", type=Path, help="chill-sitting.csv to take the median from")
    parser.add_argument("--build", choices=("before", "after"), help="row of --sitting")
    parser.add_argument("--build-tree", default="", help="source tree of --api")
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.ms is None:
        if args.sitting is None or args.build is None:
            parser.error("give --ms, or --sitting and --build")
        import csv
        with args.sitting.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if (row["build"], row["molecules"], row["threads"]) == (
                        args.build, str(MOLECULES), "1"):
                    args.ms = float(row["median_ms"])
                    args.build_tree = args.build_tree or row["build_tree"]
        if args.ms is None:
            raise SystemExit(f"no {args.build} {MOLECULES} t1 row in {args.sitting}")
    import paper_stats as stats
    env = dict(os.environ, OMP_NUM_THREADS="1", RAYON_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    with tempfile.TemporaryDirectory() as tmp:
        data = Path(tmp) / "perf.data"
        run = subprocess.run(
            stats.pinned("0", 1, None, [
                "perf", "record", "-q", "-F", "2999", "-o", str(data),
                sys.executable, str(Path(__file__).resolve()), "--drive",
                str(args.api.resolve()), str(args.reps)]),
            env=env, check=True, capture_output=True, text=True)
        report = subprocess.run(
            ["perf", "report", "-i", str(data), "--no-children", "--sort", "dso,symbol",
             "--stdio", "-q", "--percent-limit", "0"],
            check=True, capture_output=True, text=True).stdout
    top = [line.rstrip()[:160] for line in report.splitlines() if line.strip()][:25]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join([
        f"source_tree={stats.source_tree(ROOT)}",
        f"build_tree={args.build_tree}",
        f"library={os.path.relpath(args.api.resolve(), ROOT)}",
        f"molecules={MOLECULES} rule=ovito graph=cutoff threads=1 cpus=0 reps={args.reps}",
        f"sitting_median_ms={args.ms}",
        run.stdout.strip(),
        stats.stderr_record(run.stderr),
        "",
        *group(report, args.ms),
        "",
        "top symbols by self time:",
        *top,
    ]) + "\n", encoding="utf-8")
    print(args.out.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
