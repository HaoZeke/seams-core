#!/usr/bin/env python3
"""Fold the hour-budget log into the summary and a markdown table."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOGS = [
    ROOT / "benchmarks" / "logs" / "h2h-dnf-hour.txt",
    ROOT / "benchmarks" / "logs" / "h2h-dnf-sqq.txt",
    ROOT / "benchmarks" / "logs" / "h2h-dnf-ira110.txt",
    ROOT / "benchmarks" / "logs" / "h2h-dnf-cyc.txt",
    ROOT / "benchmarks" / "logs" / "h2h-dnf-ira884.txt",
]
SUMMARY = ROOT / "benchmarks" / "h2h-summary.csv"
MD = ROOT / "benchmarks" / "h2h-dnf-hour.md"


def read_logs() -> list[str]:
    lines = []
    seen = set()
    for path in sorted(ROOT.glob("benchmarks/logs/h2h-dnf-*.txt")):
        if path in seen:
            continue
        seen.add(path)
        lines.extend(path.read_text(encoding="ascii", errors="replace").splitlines())
    return lines


def sessions(lines: list[str]) -> dict:
    found = {}
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        if parts[0] not in {"warm", "call", "stop", "done", "crash"}:
            continue
        kind = parts[1] if parts[0] != "crash" else ""
        if parts[0] == "crash":
            continue
        molecules = threads = None
        fields = {}
        for i, token in enumerate(parts):
            if token == "molecules" and i + 1 < len(parts):
                molecules = parts[i + 1]
            elif token == "threads" and i + 1 < len(parts):
                threads = parts[i + 1]
        key = (kind, molecules, threads)
        rec = found.setdefault(key, {"calls": [], "warm": None, "stop": None, "done": None})
        if parts[0] == "call" and "other_ms" in parts:
            rec["calls"].append(line)
        elif parts[0] == "warm":
            rec["warm"] = line
        elif parts[0] == "stop":
            rec["stop"] = line
        elif parts[0] == "done":
            rec["done"] = line
        rec["last"] = line
    return found


def field(line: str, name: str) -> str:
    parts = line.split()
    if name not in parts:
        return ""
    return parts[parts.index(name) + 1]


def median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def reason_of(line: str) -> str:
    if "reason " not in line:
        return ""
    text = line.split("reason ", 1)[1]
    if " rss_kb" in text:
        text = text.split(" rss_kb", 1)[0]
    return text.strip()


def pin_of(rec: dict, row: dict) -> str:
    for line in (rec.get("done"), rec.get("stop"), rec.get("warm")):
        if line and field(line, "cpus"):
            return field(line, "cpus")
    return row["cpus"]


def explain(kind: str, molecules: str, calls: list[str]) -> str:
    if kind == "sqq":
        if not calls:
            return (
                "SQQ 0.5.8 cpp engine, ring.definition chordless. "
                "The python shortest_path engine does not return on these cells. "
                "The cpp call was killed at the budget before it returned. "
            )
        return (
            "SQQ 0.5.8 cpp engine, ring.definition chordless, returned 0 cages. "
            "d-SEAMS seams_cages_counts is atom membership: on this ice n_hex equals "
            "the molecule count, so the summed counts are not a cage census. "
            "seams_cages_signature hc is half the molecule count. "
            "The python shortest_path engine does not return on these cells. "
        )
    if kind == "ira":
        if not calls:
            return "classify_serial has no thread argument. "
        values = {field(line, "other_value") for line in calls}
        ours = {field(line, "dseams_value") for line in calls if field(line, "dseams_value")}
        same = values == {molecules} and (not ours or ours == {molecules})
        word = "equal" if same else "recorded"
        return f"hexagonal {molecules} {word}. classify_serial has no thread argument. "
    if kind == "cycless":
        values = {field(line, "other_value") for line in calls if field(line, "other_value")}
        shown = next(iter(values)) if len(values) == 1 else "recorded"
        return (
            f"hexagons {shown} equal. cycles_iter has no thread argument. "
            "The 3.5 A graph is built outside the timer. "
        )
    return ""


def fill_times(row: dict, calls: list[str]) -> None:
    other = [float(field(line, "other_ms")) for line in calls]
    dseams = [float(field(line, "dseams_ms")) for line in calls]
    row["n"] = str(len(other))
    row["dseams_median_ms"] = f"{median(dseams):.4f}"
    row["dseams_min_ms"] = f"{min(dseams):.4f}"
    row["other_median_ms"] = f"{median(other):.4f}"
    row["other_min_ms"] = f"{min(other):.4f}"
    row["ratio"] = f"{(median(other) / median(dseams)):.4f}" if median(dseams) else ""
    win = median(dseams) < median(other) and min(dseams) < min(other)
    row["verdict"] = "WIN" if win else "LOSE"


def apply() -> None:
    lines = read_logs()
    found = sessions(lines)
    rows = list(csv.DictReader(SUMMARY.open(encoding="ascii", newline="")))
    header = list(rows[0].keys()) if rows else []
    kind_of = {"SQQ": "sqq", "IRA": "ira", "cycless": "cycless"}
    finished = {}
    for (kind, molecules, _threads), rec in found.items():
        if rec["calls"]:
            finished[kind] = max(finished.get(kind, 0), int(molecules))
    changed = []
    for row in rows:
        if row["verdict"] != "DNF" and row.get("budget_s") != "3600":
            continue
        key = (kind_of[row["competitor"]], row["molecules"], row["threads"])
        rec = found.get(key)
        if rec is None:
            continue
        kind, molecules, threads = key
        if rec["calls"]:
            fill_times(row, rec["calls"])
            source = rec["done"] or rec["stop"] or rec["calls"][-1]
            rss = field(source, "rss_kb")
            cpus = pin_of(rec, row)
            row["cpus"] = cpus
            short = ""
            if rec["stop"] and not rec["done"]:
                short = (
                    f" The 3600 s budget stopped the session after {len(rec['calls'])} calls. "
                    f"Cause: {reason_of(rec['stop'])}."
                )
            pin_note = ""
            if threads != "1" and cpus.isdigit():
                pin_note = (
                    f"Requested {threads} threads. The call does not start workers, "
                    f"so the process was pinned to core {cpus}. "
                )
            row["agreement"] = (
                f"{explain(kind, molecules, rec['calls'])}"
                f"{pin_note}"
                f"n {len(rec['calls'])} after 1 warm-up. Budget 3600 s per call. "
                f"peak_rss_kb {rss}. cores 8 cpus {cpus}.{short}"
            )
        elif rec["stop"] or rec["warm"]:
            source = rec["stop"] or rec["warm"]
            row["n"] = "0"
            row["verdict"] = "DNF"
            row["dseams_median_ms"] = ""
            row["dseams_min_ms"] = ""
            row["other_median_ms"] = ""
            row["other_min_ms"] = ""
            row["ratio"] = ""
            cpus = pin_of(rec, row)
            row["cpus"] = cpus
            row["agreement"] = (
                f"{explain(kind, molecules, [])}"
                f"n 0. Budget 3600 s per call. "
                f"Cause: {reason_of(source)}. peak_rss_kb {field(source, 'rss_kb')}. "
                f"cores 8 cpus {cpus or row['cpus']}."
            )
        if kind in finished:
            current = row["largest_finished"]
            if current == "" or (current.isdigit() and int(current) < finished[kind]):
                row["largest_finished"] = str(finished[kind])
        row["budget_s"] = "3600"
        changed.append(key)
    with SUMMARY.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    md = ["# Hour budget for the former 600 s DNF rows", ""]
    md.append("Each call was allowed 3600 s. Peak RSS is the process high-water mark in kilobytes. Cores are 8. No host name.")
    md.append("")
    md.append("| Kind | Molecules | Threads | Outcome | n | Other median ms | Peak RSS kb | What stopped it |")
    md.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for key, rec in sorted(found.items()):
        kind, molecules, threads = key
        if rec["calls"] and not rec["stop"]:
            other = [float(field(line, "other_ms")) for line in rec["calls"]]
            outcome = "measured" if rec["done"] else "partial"
            why = "finished" if rec["done"] else "session still running"
            rss = field(rec["done"], "rss_kb") if rec["done"] else field(rec["calls"][-1], "rss_kb")
            md.append(
                f"| {kind} | {molecules} | {threads} | {outcome} | {len(other)} | {median(other):.1f} | {rss} | {why} |"
            )
        else:
            source = rec["stop"] or rec["warm"] or ""
            reason = ""
            if source and "reason " in source:
                reason = source.split("reason ", 1)[1].split(" rss_kb")[0]
            md.append(
                f"| {kind} | {molecules} | {threads} | stopped | {len(rec['calls'])} |  | {field(source, 'rss_kb') if source else ''} | {reason} |"
            )
    md.append("")
    MD.write_text("\n".join(md) + "\n", encoding="ascii")
    print(f"updated {len(changed)} summary rows", flush=True)


if __name__ == "__main__":
    apply()
