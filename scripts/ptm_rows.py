#!/usr/bin/env python3
"""Time the PTM pair and CHILL+ on one frame, five warmups and 30 repeats.

The d-SEAMS side is tests/bench_paper. seeded-adj and chill both build
their neighbour lists inside the timer. PTM is OVITO in this process.
Samples are the paper_rows protocol: CORE, WARMUP and REPS from
timing_protocol.py, and PTM_RMSD from paper_rows.py.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

from paper_rows import PTM_RMSD, ptm_pair, run_bench
from render_timings import order_statistic_interval
from timing_protocol import CORE, WARMUP, REPS

ROOT = Path(__file__).resolve().parents[1]


def host_lines():
    cpu = "unknown"
    family = model = ""
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name") and cpu == "unknown":
            cpu = line.split(":", 1)[1].strip()
        if line.startswith("cpu family") and not family:
            family = line.split(":", 1)[1].strip()
        if line.startswith("model\t") and not model:
            model = line.split(":", 1)[1].strip()
    rev = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    return [
        f"cpu {cpu}",
        f"cpuid family {family} model {model}",
        f"commit {rev}",
        f"python {sys.version.split()[0]}",
        f"protocol warmup {WARMUP} reps {REPS} core {CORE}",
        f"rmsd_cutoff {PTM_RMSD}",
        "neighbour lists are inside the bench_paper timer for chill and seeded-adj",
        "chill is d-SEAMS getCorrelPlus and getIceTypePlusNoPrint on neighListO at 3.5 A",
    ]


def stats_line(name, samples):
    values = [row[0] for row in samples]
    stats = order_statistic_interval(values)
    return (
        f"{name} n={stats['n']} median={stats['median']:.6f} "
        f"low={stats['low']:.6f} high={stats['high']:.6f} "
        f"min={stats['minimum']:.6f} max={stats['maximum']:.6f}"
    )


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: ptm_rows.py traj log")
    traj = Path(sys.argv[1])
    log = Path(sys.argv[2])
    os.environ["OMP_NUM_THREADS"] = "1"
    os.sched_setaffinity(0, {CORE})
    load_start = Path("/proc/loadavg").read_text().strip()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    side = log.with_suffix(".chill-labels")
    chill_samples, chill_cmd = run_bench("chill", traj, side)
    pair = ptm_pair(traj)
    finished = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    load_end = Path("/proc/loadavg").read_text().strip()
    lines = host_lines()
    lines.append(f"started {started}")
    lines.append(f"finished {finished}")
    lines.append(f"load_start {load_start}")
    lines.append(f"load_end {load_end}")
    lines.append("chill_cmd " + " ".join(chill_cmd))
    lines.append("dseams_cmd " + " ".join(str(x) for x in pair[0]["command"]))
    lines.append("ovito_cmd " + " ".join(str(x) for x in pair[1]["command"]))
    lines.append(stats_line("chill", chill_samples))
    lines.append("chill_ms")
    lines.extend(f"{row[0]:.6f}" for row in chill_samples)
    for session in pair:
        lines.append(stats_line(session["id"], session["samples"]))
        lines.append(f"agreement {session['agreement']}")
        lines.append(session["id"] + "_ms")
        lines.extend(f"{row[0]:.6f}" for row in session["samples"])
    chill_med = order_statistic_interval([row[0] for row in chill_samples])["median"]
    seeded_med = order_statistic_interval([row[0] for row in pair[0]["samples"]])["median"]
    lines.append(f"ratio seeded-adj/chill {seeded_med / chill_med:.6f}")
    log.write_text("\n".join(lines) + "\n")
    print(log)


if __name__ == "__main__":
    main()
