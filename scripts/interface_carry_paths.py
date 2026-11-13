#!/usr/bin/env python3
"""Paired carry and rebuild timings for the nucleation path switch.

The driver times each transition as a carried update followed by a
rebuild, repeated for every replicate in one session. A paired
bootstrap of 10000 draws, seed 20261103, gives a 95% interval on the
speedup of each transition and of each path class. A class whose
interval lies above 1 clears 1; an interval that contains 1 is parity.
On a full pass the overhead is the paired difference, carry minus
rebuild, in milliseconds: the rebuild plus the row count.

The per-transition table is data/interface-carry-paths.csv. The class
intervals are data/interface-carry-classes.csv.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261103
DRAWS = 10000
COVERAGE = 0.95

PATH_FIELDS = [
    "graph",
    "frame",
    "path",
    "recomputed",
    "carry_ms",
    "rebuild_ms",
    "speedup",
    "ci_lo",
    "ci_hi",
    "commit",
]
CLASS_FIELDS = [
    "graph",
    "path",
    "transitions",
    "speedup",
    "ci_lo",
    "ci_hi",
    "verdict",
    "overhead_ms",
    "overhead_lo",
    "overhead_hi",
    "cores",
    "omp_threads",
    "reps",
    "commit",
]


def type7(values: list[float], p: float) -> float:
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise ValueError("empty sample")
    if n == 1:
        return ordered[0]
    h = (n - 1) * p
    lo = int(h)
    hi = int(h) if h == int(h) else int(h) + 1
    if hi >= n:
        hi = n - 1
    weight = h - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def median(values: list[float]) -> float:
    return type7(values, 0.5)


def percentile_ends(draws: list[float]) -> tuple[float, float]:
    tail = (1.0 - COVERAGE) / 2.0
    return type7(draws, tail), type7(draws, 1.0 - tail)


def mix_seed(*parts: object) -> int:
    acc = SEED
    for part in parts:
        for ch in str(part):
            acc = (acc * 131 + ord(ch)) & 0xFFFFFFFF
    return acc or 1


def verdict(lo: float, hi: float) -> str:
    if lo > 1.0:
        return "clears"
    if hi < 1.0:
        return "slower"
    return "parity"


def frame_speedup(pairs: list[tuple[float, float]], seed: int) -> tuple[float, float, float]:
    carry = median([c for c, _r in pairs])
    rebuild = median([r for _c, r in pairs])
    point = rebuild / carry if carry > 0.0 else float("nan")
    rng = random.Random(seed)
    n = len(pairs)
    draws = []
    for _ in range(DRAWS):
        cs = []
        rs = []
        for _k in range(n):
            c, r = pairs[rng.randrange(n)]
            cs.append(c)
            rs.append(r)
        mc = median(cs)
        draws.append(median(rs) / mc if mc > 0.0 else float("nan"))
    lo, hi = percentile_ends(draws)
    return point, lo, hi


def class_speedup(frames: list[list[tuple[float, float]]], seed: int) -> tuple[float, float, float]:
    n = len(frames[0])

    def ratio(indices: list[int]) -> float:
        carry = 0.0
        rebuild = 0.0
        for pairs in frames:
            carry += median([pairs[i][0] for i in indices])
            rebuild += median([pairs[i][1] for i in indices])
        return rebuild / carry if carry > 0.0 else float("nan")

    point = ratio(list(range(n)))
    rng = random.Random(seed)
    draws = []
    for _ in range(DRAWS):
        draws.append(ratio([rng.randrange(n) for _k in range(n)]))
    lo, hi = percentile_ends(draws)
    return point, lo, hi


def class_overhead(frames: list[list[tuple[float, float]]], seed: int) -> tuple[float, float, float]:
    """Median extra milliseconds of the carried update over the rebuild."""
    n = len(frames[0])

    def extra(indices: list[int]) -> float:
        per_frame = []
        for pairs in frames:
            per_frame.append(median([pairs[i][0] - pairs[i][1] for i in indices]))
        return median(per_frame)

    point = extra(list(range(n)))
    rng = random.Random(seed)
    draws = []
    for _ in range(DRAWS):
        draws.append(extra([rng.randrange(n) for _k in range(n)]))
    lo, hi = percentile_ends(draws)
    return point, lo, hi


def child_env(bench: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", os.environ.get("SCIENCE_THREADS", "2"))
    env["OMP_DYNAMIC"] = "FALSE"
    origin = bench.resolve().parent
    libs = [
        Path("/cvmfs/software.eessi.io/versions/2026.06/compat/linux/x86_64/lib64"),
        origin.parent / "src",
        origin.parent / "chemfiles-prefix" / "lib",
        origin.parent / "subprojects" / "readcon-core",
        origin.parent / "subprojects" / "readcon-db",
        Path("/lib/x86_64-linux-gnu"),
    ]
    parts = [str(path) for path in libs if path.is_dir()]
    old = env.get("LD_LIBRARY_PATH", "")
    if old:
        parts.append(old)
    if parts:
        env["LD_LIBRARY_PATH"] = ":".join(parts)
    return env


def run_bench(bench: Path, dump: Path, reps: int, warmup: int, last: int) -> str:
    cmd = [
        str(bench),
        "time",
        str(dump),
        "--graph",
        "both",
        "--reps",
        str(reps),
        "--warmup",
        str(warmup),
        "--last",
        str(last),
        "--type",
        "1",
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=child_env(bench), text=True, capture_output=True)
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(f"bench exited {done.returncode}")
    if "# order abab" not in done.stdout.splitlines():
        raise SystemExit("bench did not interleave carry and rebuild")
    return done.stdout


def parse_log(text: str) -> tuple[dict[str, str], dict[tuple[str, int], dict]]:
    header = {"cores": "", "omp_threads": "", "reps": ""}
    frames: dict[tuple[str, int], dict] = {}
    for line in text.splitlines():
        if line.startswith("# cores "):
            header["cores"] = line.split()[2]
        elif line.startswith("# omp_threads "):
            header["omp_threads"] = line.split()[2]
        elif line.startswith("# reps "):
            header["reps"] = line.split()[2]
        elif line.startswith("frame "):
            parts = line.split()
            # frame stride S graph G index I ... recomputed R ... path P ...
            graph = parts[4]
            index = int(parts[6])
            row = frames.setdefault((graph, index), {
                "graph": graph,
                "frame": index,
                "path": "",
                "recomputed": "",
                "pairs": [],
            })
            row["recomputed"] = int(parts[10])
            row["path"] = parts[14]
        elif line.startswith("sample "):
            parts = line.split()
            graph = parts[4]
            frame = int(parts[6])
            carry = float(parts[10])
            rebuild = float(parts[12])
            # The opening frame is primed, not a transition, so it has
            # samples and no summary line. It is dropped below.
            row = frames.setdefault((graph, frame), {
                "graph": graph,
                "frame": frame,
                "path": "",
                "recomputed": "",
                "pairs": [],
            })
            row["pairs"].append((carry, rebuild))
    if not frames:
        raise SystemExit("bench wrote no transitions")
    return header, frames


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    return done.stdout.strip()


def build_tables(header: dict[str, str], frames: dict, commit: str):
    paths = []
    by_class: dict[tuple[str, str], list] = {}
    for key in sorted(frames, key=lambda item: (item[0], item[1])):
        row = frames[key]
        if row["path"] in ("", "cold") or not row["pairs"]:
            continue
        point, lo, hi = frame_speedup(row["pairs"], mix_seed(row["graph"], row["frame"]))
        carry = median([c for c, _r in row["pairs"]])
        rebuild = median([r for _c, r in row["pairs"]])
        paths.append({
            "graph": row["graph"],
            "frame": str(row["frame"]),
            "path": row["path"],
            "recomputed": str(row["recomputed"]),
            "carry_ms": f"{carry:.6f}",
            "rebuild_ms": f"{rebuild:.6f}",
            "speedup": f"{point:.6f}",
            "ci_lo": f"{lo:.6f}",
            "ci_hi": f"{hi:.6f}",
            "commit": commit,
        })
        by_class.setdefault((row["graph"], row["path"]), []).append(row["pairs"])
    classes = []
    for (graph, path) in sorted(by_class):
        groups = by_class[(graph, path)]
        point, lo, hi = class_speedup(groups, mix_seed(graph, path, "speedup"))
        extra, elo, ehi = class_overhead(groups, mix_seed(graph, path, "overhead"))
        classes.append({
            "graph": graph,
            "path": path,
            "transitions": str(len(groups)),
            "speedup": f"{point:.6f}",
            "ci_lo": f"{lo:.6f}",
            "ci_hi": f"{hi:.6f}",
            "verdict": verdict(lo, hi),
            "overhead_ms": f"{extra:.6f}",
            "overhead_lo": f"{elo:.6f}",
            "overhead_hi": f"{ehi:.6f}",
            "cores": header.get("cores", ""),
            "omp_threads": header.get("omp_threads", ""),
            "reps": header.get("reps", ""),
            "commit": commit,
        })
    return paths, classes


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def report(classes: list[dict]) -> None:
    for row in classes:
        print(
            f"class graph {row['graph']} path {row['path']} "
            f"transitions {row['transitions']} speedup {row['speedup']} "
            f"ci_lo {row['ci_lo']} ci_hi {row['ci_hi']} verdict {row['verdict']} "
            f"overhead_ms {row['overhead_ms']} "
            f"overhead_lo {row['overhead_lo']} overhead_hi {row['overhead_hi']}"
        )


def self_test() -> None:
    equal = [(10.0, 10.0)] * 8
    point, lo, hi = frame_speedup(equal, 1)
    if not (point == 1.0 and lo <= 1.0 <= hi and verdict(lo, hi) == "parity"):
        raise SystemExit("equal pairs are not parity")
    faster = [(10.0, 12.0)] * 8
    point, lo, hi = frame_speedup(faster, 2)
    if not (abs(point - 1.2) < 1e-9 and verdict(lo, hi) == "clears"):
        raise SystemExit("a faster rebuild did not clear 1")
    slower = [(12.0, 10.0)] * 8
    _point, lo, hi = frame_speedup(slower, 3)
    if verdict(lo, hi) != "slower":
        raise SystemExit("a slower carry was not slower")
    groups = [[(10.2 + i * 0.01, 10.0 + i * 0.01) for i in range(8)]]
    extra, elo, ehi = class_overhead(groups, 4)
    if not (elo <= extra <= ehi and 0.15 < extra < 0.25):
        raise SystemExit(f"overhead {extra} not near 0.2")
    text = "\n".join([
        "# cores 8",
        "# omp_threads 2",
        "# reps 8",
        "# order abab",
        "sample stride 0 graph cutoff frame 1 rep 1 carry_ms 6.000000 rebuild_ms 6.000000",
        "sample stride 0 graph cutoff frame 2 rep 1 carry_ms 10.000000 rebuild_ms 10.000000",
        "sample stride 0 graph cutoff frame 2 rep 2 carry_ms 10.000000 rebuild_ms 11.000000",
        "frame stride 0 graph cutoff index 2 changed_bonds 1 recomputed 4 "
        "nop 4 path full carry_ms 1 rebuild_ms 1 speedup 1 identical yes",
    ])
    _header, frames = parse_log(text)
    if frames[("cutoff", 2)]["path"] != "full":
        raise SystemExit("parser missed the path")
    if len(frames[("cutoff", 2)]["pairs"]) != 2:
        raise SystemExit("parser missed a pair")
    print("self-test ok")


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path, required=True)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=8)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--last", type=int, default=8)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--classes", type=Path, required=True)
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    if args.reps < 8:
        raise SystemExit("reps must be at least 8")
    if not args.bench.is_file():
        raise SystemExit(f"missing bench {args.bench}")
    if not args.dump.is_file():
        raise SystemExit(f"missing dump {args.dump}")
    text = run_bench(args.bench, args.dump, args.reps, args.warmup, args.last)
    header, frames = parse_log(text)
    commit = args.commit or head_commit()
    paths, classes = build_tables(header, frames, commit)
    write_csv(args.out, PATH_FIELDS, paths)
    write_csv(args.classes, CLASS_FIELDS, classes)
    report(classes)
    print(f"wrote {args.out}")
    print(f"wrote {args.classes}")


if __name__ == "__main__":
    main()
