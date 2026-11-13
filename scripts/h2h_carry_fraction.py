#!/usr/bin/env python3
"""Carried rings at or under a bare ringNetwork call.

Moved fractions 0, 0.05, ..., 1 on a 4096-molecule cell. Eight samples.
One thread on core 0, then eight threads on cores 0-7. A row passes when
the carried median and the carried minimum are at most the bare search.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def version_of(lib: Path) -> str:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(lib.parent)
    out = subprocess.check_output(
        [
            sys.executable, "-c",
            "import ctypes,sys; lib=ctypes.CDLL(sys.argv[1]); "
            "lib.seams_version_string.restype=ctypes.c_char_p; "
            "print(lib.seams_version_string().decode())",
            str(lib),
        ],
        env=env,
        text=True,
    )
    return out.strip()


def run_width(bench: Path, threads: int, samples: int) -> list[dict]:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(bench.resolve().parents[1] / "src")
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    cores = "0" if threads == 1 else "0-7"
    done = subprocess.run(
        ["taskset", "-c", cores, str(bench), "4096", str(samples), "sweep"],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
    )
    if done.returncode not in (0, 1):
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode)
    rows = list(csv.DictReader(done.stdout.splitlines()))
    if not rows:
        sys.stderr.write(done.stderr)
        raise SystemExit("sweep produced no rows")
    for row in rows:
        row["threads"] = str(threads)
        row["cpus"] = cores
    if done.returncode != 0 or any(row["within"] != "yes" for row in rows):
        for row in rows:
            if row["within"] != "yes":
                print(
                    f"over threads {threads} fraction {row['moved_fraction']} "
                    f"carry {row['carry_median_ms']} bare {row['bare_median_ms']}",
                    flush=True,
                )
        raise SystemExit(1)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("bench")
    ap.add_argument("--csv")
    ap.add_argument("--samples", type=int, default=8)
    args = ap.parse_args()
    bench = Path(args.bench)
    if not bench.is_file():
        raise SystemExit(f"missing {bench}")
    lib = bench.resolve().parents[1] / "src" / "libseams.so"
    seams = version_of(lib) if lib.is_file() else ""
    cores = os.cpu_count() or 1
    model = cpu_model()
    rows = []
    for threads in (1, 8):
        for row in run_width(bench, threads, args.samples):
            row["cores"] = str(cores)
            row["version"] = seams
            row["cpu_model"] = model
            rows.append(row)
            print(
                f"threads {threads} fraction {row['moved_fraction']} "
                f"carry {row['carry_median_ms']} bare {row['bare_median_ms']} "
                f"within {row['within']}",
                flush=True,
            )
    if args.csv:
        fields = [
            "threads", "cores", "cpus", "moved_fraction", "recomputed",
            "molecules", "samples", "carry_median_ms", "carry_min_ms",
            "bare_median_ms", "bare_min_ms", "within", "version", "cpu_model",
        ]
        out = Path(args.csv)
        with out.open("w", encoding="ascii", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {out.name} rows {len(rows)} cores {cores}", flush=True)


if __name__ == "__main__":
    main()
