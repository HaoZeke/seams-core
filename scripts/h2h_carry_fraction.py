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


def four_of_five(bench: Path, threads: int, samples: int) -> dict:
    """Every molecule whose index is not a multiple of 5."""
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(bench.resolve().parents[1] / "src")
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    cores = "0" if threads == 1 else "0-7"
    done = subprocess.run(
        ["taskset", "-c", cores, str(bench), "4096", str(samples)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        sys.stderr.write(done.stdout)
        raise SystemExit(done.returncode or 1)
    got = {}
    for line in done.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            got[parts[0]] = parts
    carry = float(got["carry_median"][1])
    carry_min = float(got["carry_median"][3])
    bare = float(got["full_median"][1])
    bare_min = float(got["full_median"][3])
    within = carry <= bare and carry_min <= bare_min and got["identical"][1] == "yes"
    if not within:
        raise SystemExit(
            f"four-of-five threads {threads} carry {carry:.4f} bare {bare:.4f}"
        )
    return {
        "threads": str(threads),
        "cpus": cores,
        "moved_fraction": "0.8000",
        "pattern": "four_of_five",
        "recomputed": got["recomputed"][1],
        "molecules": got["molecules"][1],
        "samples": got["molecules"][3],
        "carry_median_ms": f"{carry:.4f}",
        "carry_min_ms": f"{carry_min:.4f}",
        "bare_median_ms": f"{bare:.4f}",
        "bare_min_ms": f"{bare_min:.4f}",
        "within": "yes",
    }


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
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
        for row in rows:
            row["commit"] = commit
            row["pattern"] = "i_mod_20"
        # The four-of-five frame is every molecule whose index is not a
        # multiple of 5. The modulus-20 row at 0.8 is a different set.
        archival = {
            "threads": "8", "cores": str(cores), "cpus": "0-7",
            "moved_fraction": "0.8000", "pattern": "four_of_five",
            "recomputed": "4096", "molecules": "4096", "samples": "8",
            "carry_median_ms": "4.8771", "carry_min_ms": "4.7689",
            "bare_median_ms": "4.6538", "bare_min_ms": "4.4966",
            "within": "no", "version": seams, "cpu_model": model,
            "commit": "c38ed3c592caaf48493fc78bb66644d21e2b96bb",
        }
        live = four_of_five(bench, 8, args.samples)
        live["cores"] = str(cores)
        live["version"] = seams
        live["cpu_model"] = model
        live["commit"] = commit
        rows.append(live)
        rows.append(archival)
        fields = [
            "threads", "cores", "cpus", "moved_fraction", "pattern", "recomputed",
            "molecules", "samples", "carry_median_ms", "carry_min_ms",
            "bare_median_ms", "bare_min_ms", "within", "version", "commit",
            "cpu_model",
        ]
        out = Path(args.csv)
        with out.open("w", encoding="ascii", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {out.name} rows {len(rows)} cores {cores}", flush=True)


if __name__ == "__main__":
    main()
