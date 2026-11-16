#!/usr/bin/env python3
"""Carried-ring crossover. Eight samples, one thread and eight threads.

The partial search is used while the changed vertices stay under 0.55 of
the frame on one thread and under 0.33 with a team. Above that, and when
the partial search would not fill the team, the update calls ringNetwork.
Each run goes through scripts/pin_run.py, which sets OMP_NUM_THREADS and
RAYON_NUM_THREADS to the width; its cpus line follows the run's output.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import paper_stats as stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "bbdir-release" / "tests" / "bench_ring_carry"


def main() -> None:
    print(f"cores {os.cpu_count()}", flush=True)
    print("samples 8", flush=True)
    print("crossover one_thread 0.55 eight_threads 0.33", flush=True)
    if not BENCH.is_file():
        raise SystemExit(f"missing {BENCH.name}")
    specs = [("", "churn"), ("local", "local"), ("spread 2", "spread 2"),
             ("spread 6", "spread 6"), ("spread 15", "spread 15")]
    for width in (1, 8):
        for arg, name in specs:
            cmd = [str(BENCH), "4096", "8"]
            if arg:
                cmd.extend(arg.split())
            cmd = stats.pinned("0" if width == 1 else "0-7", width, None, cmd)
            env = os.environ.copy()
            env["OMP_PROC_BIND"] = "close"
            env["OMP_PLACES"] = "cores"
            env["OMP_WAIT_POLICY"] = "PASSIVE"
            env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
            print(f"width {width} spec {name}", flush=True)
            done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
            sys.stdout.write(done.stdout)
            print(stats.stderr_record(done.stderr), flush=True)
            if done.returncode != 0:
                sys.stderr.write(done.stderr)
                raise SystemExit(done.returncode)


if __name__ == "__main__":
    main()
