#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Paired-bootstrap ratio intervals, and the table the paper's \\num macros read.

sitting-large.csv gains ratio_vs_dseams, ratio_lo and ratio_hi on each
PTM and IRA row. paper/NUMBERS.tsv has one row per published number:
key, value, CI, source CSV, commit. The check fails when a value or a
CI endpoint is not a field of the source row.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAWS = 10000
SEED = 20261103
COVERAGE = 0.95
TAIL = (1.0 - COVERAGE) / 2.0
TSV_FIELDS = ("key", "value", "ci", "source_csv", "commit")
RATIO_FIELDS = ("ratio_vs_dseams", "ratio_lo", "ratio_hi")


def median(samples: list[float]) -> float:
    ordered = sorted(samples)
    count = len(ordered)
    if count == 0:
        raise ValueError("empty sample")
    if count % 2:
        return ordered[count // 2]
    return 0.5 * (ordered[count // 2 - 1] + ordered[count // 2])


def paired_ratio(num: list[float], den: list[float], rng: random.Random) -> tuple[float, float, float]:
    """Ratio of medians, and a paired percentile interval.

    Each draw resamples one call index and applies it to both series.
    The index convention matches paper/tex/figures/stats.py.
    """
    if len(num) != len(den) or len(num) < 2:
        raise ValueError("a paired ratio needs two samples of the same length")
    count = len(num)
    ratios = []
    for _ in range(DRAWS):
        picked = [rng.randrange(count) for _ in range(count)]
        ratios.append(median([num[i] for i in picked]) / median([den[i] for i in picked]))
    ratios.sort()
    lo = ratios[int(math.floor(TAIL * DRAWS))]
    hi = ratios[int(math.ceil((1.0 - TAIL) * DRAWS)) - 1]
    return median(num) / median(den), lo, hi


def _samples(row: dict) -> list[float]:
    return [float(part) for part in row["samples_ms"].split() if part]


def _fmt(value: float) -> str:
    return f"{value:.6f}"


def annotate_sitting_large(path: Path) -> None:
    """Fill d-SEAMS-vs-PTM and d-SEAMS-vs-IRA ratio columns from the samples."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not rows:
        raise SystemExit(f"{path} is empty")
    for name in RATIO_FIELDS:
        if name not in fields:
            fields.append(name)
    groups: dict[tuple, dict[str, dict]] = {}
    for row in rows:
        key = (row["phase"], row["molecules"], row["threads"])
        groups.setdefault(key, {})[row["operation"]] = row
    rng = random.Random(SEED)
    for group in groups.values():
        base = group.get("classifier")
        if base is None or base["code"] != "d-SEAMS":
            continue
        den = _samples(base)
        for operation in ("ptm-default", "ptm-ice", "ira"):
            other = group.get(operation)
            if other is None:
                continue
            point, lo, hi = paired_ratio(_samples(other), den, rng)
            other["ratio_vs_dseams"] = _fmt(point)
            other["ratio_lo"] = _fmt(lo)
            other["ratio_hi"] = _fmt(hi)
            note = (
                f"paired-bootstrap ratio of medians vs d-SEAMS classifier "
                f"draws={DRAWS} seed={SEED}"
            )
            if note not in other.get("detail", ""):
                other["detail"] = f"{other['detail']}; {note}"
    by_op: dict[tuple, dict[str, dict]] = {}
    for row in rows:
        if not row.get("samples_ms"):
            continue
        by_op.setdefault(
            (row.get("phase", ""), row.get("operation", ""), row.get("molecules", ""), row.get("threads", "")),
            {},
        )[row["code"]] = row
    for group in by_op.values():
        base = group.get("d-SEAMS")
        if base is None:
            continue
        den = _samples(base)
        for code in sorted(group):
            if code == "d-SEAMS":
                continue
            other = group[code]
            point, lo, hi = paired_ratio(_samples(other), den, rng)
            other["ratio_vs_dseams"] = _fmt(point)
            other["ratio_lo"] = _fmt(lo)
            other["ratio_hi"] = _fmt(hi)
            note = (
                f"paired-bootstrap ratio of medians vs d-SEAMS "
                f"draws={DRAWS} seed={SEED}"
            )
            if note not in other.get("detail", ""):
                other["detail"] = f"{other['detail']}; {note}"
    for row in rows:
        for name in RATIO_FIELDS:
            row.setdefault(name, "")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _key(row: dict, column: str) -> str:
    return ".".join((
        row.get("family", "row"),
        row.get("phase", ""),
        row.get("operation", ""),
        row.get("code", ""),
        str(row.get("molecules", "")),
        f"t{row.get('threads', '')}",
        column,
    ))


def numbers_from(csv_dir: Path) -> list[dict]:
    found = []
    for path in sorted(csv_dir.glob("*.csv")):
        rel = path.resolve().relative_to(ROOT).as_posix()
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if "median_ms" not in row or not row["median_ms"]:
                    continue
                commit = row["source_commit"]
                found.append({
                    "key": _key(row, "median_ms"),
                    "value": row["median_ms"],
                    "ci": f"{row['q1_ms']},{row['q3_ms']}" if row.get("q1_ms") else "",
                    "source_csv": rel,
                    "commit": commit,
                })
                if row.get("ratio_vs_dseams"):
                    found.append({
                        "key": _key(row, "ratio_vs_dseams"),
                        "value": row["ratio_vs_dseams"],
                        "ci": f"{row['ratio_lo']},{row['ratio_hi']}",
                        "source_csv": rel,
                        "commit": commit,
                    })
    return found


def write_numbers(csv_dir: Path, dest: Path) -> None:
    rows = numbers_from(csv_dir)
    if not rows:
        raise SystemExit(f"no numbers in {csv_dir}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {dest} rows {len(rows)}", flush=True)


def _same(left: str, right: str) -> bool:
    return abs(float(left) - float(right)) <= 1e-9 * max(1.0, abs(float(right)))


def check(tsv_path: Path) -> None:
    """Fail if a published number is not a field of its source row."""
    with tsv_path.open(newline="") as handle:
        published = list(csv.DictReader(handle, delimiter="\t"))
    if not published:
        raise SystemExit(f"{tsv_path} has no numbers")
    missing = []
    for item in published:
        source = ROOT / item["source_csv"]
        if not source.is_file():
            missing.append(f"{item['key']}: missing {item['source_csv']}")
            continue
        with source.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        parts = item["key"].split(".")
        if len(parts) < 7:
            missing.append(f"{item['key']}: key has no source identity")
            continue
        family, phase, operation, code, molecules, thread, column = parts[:7]
        threads = thread[1:] if thread.startswith("t") else thread
        hits = [
            row for row in rows
            if row.get("source_commit") == item["commit"]
            and row.get("family", family) == family
            and row.get("phase") == phase
            and row.get("operation") == operation
            and row.get("code") == code
            and str(row.get("molecules")) == molecules
            and str(row.get("threads")) == threads
            and column in row
            and row[column] != ""
            and _same(item["value"], row[column])
        ]
        if not hits:
            missing.append(f"{item['key']}: value {item['value']} has no source row")
            continue
        if item.get("ci"):
            ends = item["ci"].split(",")
            if column == "ratio_vs_dseams":
                wanted = (hits[0].get("ratio_lo", ""), hits[0].get("ratio_hi", ""))
            else:
                wanted = (hits[0].get("q1_ms", ""), hits[0].get("q3_ms", ""))
            if len(ends) != 2 or not (_same(ends[0], wanted[0]) and _same(ends[1], wanted[1])):
                missing.append(f"{item['key']}: CI {item['ci']} has no source row")
    if missing:
        raise SystemExit("numbers with no source row:\n" + "\n".join(missing))
    print(f"checked {len(published)} numbers, each has a source row", flush=True)


def publish(csv_dir: Path, tsv_path: Path) -> None:
    for path in sorted(csv_dir.glob("*.csv")):
        with path.open(newline="") as handle:
            fields = csv.DictReader(handle).fieldnames or []
        if "samples_ms" in fields and "code" in fields:
            annotate_sitting_large(path)
    write_numbers(csv_dir, tsv_path)
    check(tsv_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv-dir", type=Path, default=ROOT / "benchmarks/paper-numbers")
    parser.add_argument("--tsv", type=Path, default=ROOT / "paper/NUMBERS.tsv")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only:
        check(args.tsv)
        return 0
    publish(args.csv_dir, args.tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
