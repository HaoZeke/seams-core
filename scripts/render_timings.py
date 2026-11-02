#!/usr/bin/env python3
"""Regenerate Table 2 and the quoted stage timings from one timings JSON.

The only numbers this script prints are computed from the recorded samples.
It does not read a precomputed median from the JSON.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

from timing_protocol import REPS


def binomial_coverage(n: int, k: int) -> Fraction:
    """P(X_(k) <= median <= X_(n-k+1)) for an i.i.d. sample, k one-based.

    alpha = 0.05 gives the conventional 95 % interval.
    """
    if k < 1 or n - k + 1 < k:
        return Fraction(0)
    ways = 0
    # C(n, i) via multiplicative formula, walking upward.
    term = 1
    for i in range(0, n + 1):
        if k <= i <= n - k:
            ways += term
        if i == n:
            break
        term = term * (n - i) // (i + 1)
    return Fraction(ways, 2**n)


def order_statistic_interval(samples: list[float], alpha: float = 0.05):
    """Narrowest central order-statistic interval whose coverage is at least 1-alpha.

    alpha = 0.05 gives the conventional 95 % interval.
    """
    ordered = sorted(samples)
    n = len(ordered)
    if n == 0:
        raise ValueError("empty sample")
    chosen = None
    for k in range(1, n // 2 + 1):
        coverage = binomial_coverage(n, k)
        if coverage >= 1 - Fraction(str(alpha)):
            chosen = (k, coverage)
        else:
            break
    if chosen is None:
        k, coverage = 1, binomial_coverage(n, 1)
    else:
        k, coverage = chosen
    low = ordered[k - 1]
    high = ordered[n - k]
    if n % 2:
        median = ordered[n // 2]
    else:
        median = 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])
    return {
        "n": n,
        "k": k,
        "coverage": coverage,
        "alpha": alpha,
        "low": low,
        "high": high,
        "median": median,
        "minimum": ordered[0],
        "maximum": ordered[-1],
    }


def column(session: dict, name: str) -> list[float]:
    idx = session["columns"].index(name)
    return [float(row[idx]) for row in session["samples"]]


def fmt(value: float) -> str:
    return f"{value:.3f}"


def cell(stats: dict) -> str:
    return (
        f"{fmt(stats['median'])} "
        f"[{fmt(stats['low'])}, {fmt(stats['high'])}]"
    )


def marked_cell(stats: dict, session: dict) -> str:
    """n=30 keeps the order-statistic bracket. Any other n is labelled and
    prints the minimum and the maximum of the recorded samples."""
    text = cell(stats)
    if stats["n"] != REPS or session.get("protocol_note"):
        text += (
            f" (n={stats['n']}; min {fmt(stats['minimum'])}; "
            f"max {fmt(stats['maximum'])})"
        )
    return text


def cage_pair(session: dict):
    cages = session.get("counts", {}).get("cages")
    if not isinstance(cages, dict):
        return None
    return (
        int(cages.get("5:12", 0)),
        int(cages.get("5:12,6:2", 0)),
        int(cages.get("5:12,6:4", 0)),
    )


def compare_outputs(ours: dict, other: dict) -> str:
    left = cage_pair(ours)
    right = cage_pair(other)
    if left is not None and right is not None:
        if left == right:
            return f"same {left[0]} and {left[1]}"
        return (
            f"d-SEAMS {left[0]}/{left[1]}/{left[2]}; "
            f"other {right[0]}/{right[1]}/{right[2]}"
        )
    left_text = ours.get("counts", {}).get("stdout", "") or ours.get("counts", {}).get("ring_census", "")
    right_text = other.get("counts", {}).get("ring_census", "")
    if left_text and right_text:
        left_map = census_map(left_text)
        right_map = census_map(right_text)
        shared = sorted(set(left_map) & set(right_map))
        if shared and all(left_map[size] == right_map[size] for size in shared):
            shown = " ".join(f"{size}:{right_map[size]}" for size in shared)
            return f"same {shown}"
        return f"d-SEAMS {left_text}; other {right_text}"
    return other.get("agreement", "")


def census_map(text: str) -> dict[int, int]:
    out = {}
    for token in text.replace(";", " ").split():
        if ":" not in token:
            continue
        size, count = token.split(":", 1)
        if size.isdigit() and count.lstrip("-").isdigit():
            out[int(size)] = int(count)
    return out


def tex_escape(text: str) -> str:
    return (
        text.replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("_", "\\_")
        .replace("#", "\\#")
        .replace("^", "\\^{}")
    )


def session_stats(session: dict, name: str) -> dict:
    return order_statistic_interval(column(session, name))


def write_table(sessions: list[dict], path: Path) -> str:
    groups: dict[tuple, dict] = {}
    for session in sessions:
        if session.get("role") != "head-to-head":
            continue
        key = (session["comparison"], session["scope"], session["frame"])
        groups.setdefault(key, {})[session["side"]] = session
    lines = [
        "% Generated by scripts/render_timings.py from the timings JSON.",
        "% Each time is the sample median in milliseconds, with the",
        "% nonparametric order-statistic interval for the median beside it.",
        "\\begin{table}[htbp]",
        "\\caption{Single-core head-to-head timings regenerated from one",
        "JSON of raw repetitions. Each entry is the median of the repeats",
        "recorded in that session. The bracket is the order-statistic",
        "interval for the median. An entry marked n=3 is three recorded",
        "repetitions and no further warmup; its bracket is the minimum and",
        "the maximum, and the frame states the protocol. Whole-process rows",
        "time the public command, including its own input and output. Kernel",
        "rows time the d-SEAMS library call with the frame already in memory.}",
        "\\label{tab:compare-protocol}",
        "\\centering",
        "\\scriptsize",
        "\\begin{tabular}{p{3.2cm}p{2.4cm}p{1.6cm}p{2.6cm}p{2.6cm}p{3.4cm}}",
        "computation & frame & other code & d-SEAMS & other & agreement\\\\",
        "\\hline",
    ]
    notes = []
    for key in sorted(groups):
        sides = groups[key]
        if "dseams" not in sides:
            continue
        comparison, scope, frame = key
        ours = sides["dseams"]
        our_time = session_stats(ours, ours["columns"][0])
        for side, other in sorted(sides.items()):
            if side == "dseams":
                continue
            other_time = session_stats(other, other["columns"][0])
            agreement = compare_outputs(ours, other)
            lines.append(
                " & ".join(
                    [
                        tex_escape(f"{comparison} ({scope})"),
                        tex_escape(frame),
                        tex_escape(other.get("code_name", side)),
                        marked_cell(our_time, ours),
                        marked_cell(other_time, other),
                        tex_escape(agreement),
                    ]
                )
                + " \\\\"
            )
            notes.append(
                f"{comparison} {scope} {frame} {side}: "
                f"n={our_time['n']} k={our_time['k']} "
                f"coverage={our_time['coverage']:.6f} "
                f"dseams {cell(our_time)} other {cell(other_time)} "
                f"agreement={agreement}"
            )
    lines += [
        "\\end{tabular}",
        "\\end{table}",
        "",
        "\\begin{table}[htbp]",
        "\\caption{d-SEAMS library kernel on the same frames, median and",
        "order-statistic interval. This timer includes the cutoff neighbour",
        "list. It is not paired with a competitor kernel because TRACE, HTR,",
        "SQQ, GRADE and R.I.N.G.S.\\ are timed as whole processes.}",
        "\\centering",
        "\\scriptsize",
        "\\begin{tabular}{p{4cm}p{3cm}p{4cm}p{4cm}}",
        "computation & frame & median [interval] & output\\\\",
        "\\hline",
    ]
    paired = {}
    for session in sessions:
        if session.get("role") != "head-to-head":
            continue
        key = (session.get("comparison"), session.get("scope"), session.get("frame"))
        paired.setdefault(key, set()).add(session.get("side"))
    for session in sessions:
        if session.get("scope") != "kernel" or session.get("role") == "quoted-path":
            continue
        key = (session.get("comparison"), session.get("scope"), session.get("frame"))
        if len(paired.get(key, ())) > 1:
            continue
        stats = session_stats(session, session["columns"][0])
        lines.append(
            " & ".join(
                [
                    tex_escape(session.get("comparison", "")),
                    tex_escape(session.get("frame", "")),
                    cell(stats),
                    tex_escape(session.get("agreement", "")),
                ]
            )
            + " \\\\"
        )
        notes.append(
            f"kernel {session['id']} n={stats['n']} {cell(stats)} "
            f"{session.get('agreement', '')}"
        )
    lines += ["\\end{tabular}", "\\end{table}", ""]
    path.write_text("\n".join(lines))
    return "\n".join(notes)


def write_quoted(sessions: list[dict], path: Path) -> str:
    by_id = {session["id"]: session for session in sessions}
    lines = [
        "% Quoted stage timings from the same JSON.",
        "% A paired comparison below uses one session for every figure in the pair.",
        "\\begin{table}[htbp]",
        "\\caption{Ring and cage stage medians on the 4096-molecule mW cubic",
        "frame, one session for the current engine and one session for the",
        "baseline. The path column is the end-to-end timer around ring",
        "detection, the ring index, \\texttt{findHC} and \\texttt{findDDC}.",
        "The combined column is one timer around \\texttt{findHC} and",
        "\\texttt{findDDC}. Brackets are order-statistic intervals.}",
        "\\label{tab:path-protocol}",
        "\\centering",
        "\\scriptsize",
        "\\begin{tabular}{lrrrrr}",
        "engine & rings & findHC & findDDC & combined & path\\\\",
        "\\hline",
    ]
    summary = []
    for sid, label in (
        ("path-current", "2.0"),
        ("path-baseline", "baseline"),
    ):
        session = by_id.get(sid)
        if session is None:
            lines.append(f"% missing session {sid}")
            continue
        stats = {name: session_stats(session, name) for name in session["columns"]}
        lines.append(
            " & ".join(
                [
                    label,
                    cell(stats["ring_ms"]),
                    cell(stats["hc_ms"]),
                    cell(stats["ddc_ms"]),
                    cell(stats["combined_ms"]),
                    cell(stats["path_ms"]),
                ]
            )
            + " \\\\"
        )
        # Per-repeat sums, still this session: do not add the medians.
        rings = column(session, "ring_ms")
        hc = column(session, "hc_ms")
        ddc = column(session, "ddc_ms")
        summed = [a + b + c for a, b, c in zip(rings, hc, ddc)]
        alone = [b + c for b, c in zip(hc, ddc)]
        sum_stats = order_statistic_interval(summed)
        alone_stats = order_statistic_interval(alone)
        summary.append(
            f"{sid} n={stats['path_ms']['n']} k={stats['path_ms']['k']} "
            f"coverage={stats['path_ms']['coverage']:.6f}\n"
            f"  rings {cell(stats['ring_ms'])}\n"
            f"  findHC {cell(stats['hc_ms'])}\n"
            f"  findDDC {cell(stats['ddc_ms'])}\n"
            f"  combined {cell(stats['combined_ms'])}\n"
            f"  path {cell(stats['path_ms'])}\n"
            f"  median of per-repeat (ring+HC+DDC) {cell(sum_stats)}\n"
            f"  median of per-repeat (HC+DDC alone timers) {cell(alone_stats)}\n"
            f"  counts rings={session.get('counts', {}).get('rings')} "
            f"hc={session.get('counts', {}).get('hc_rings')} "
            f"ddc={session.get('counts', {}).get('ddc_rings')} "
            f"cages={session.get('counts', {}).get('cages')}"
        )
    lines += ["\\end{tabular}", "\\end{table}", ""]
    lines += rings_quote(by_id, summary)
    path.write_text("\n".join(lines))
    return "\n".join(summary)


def rings_quote(by_id: dict, summary: list[str]) -> list[str]:
    """n=3 mW rows. Every printed figure is an order statistic of recorded samples."""
    pairs = (
        ("rings-kernel-mW-n3", "rings-stage-mW-n3", "ring stage"),
        ("rings-fp-mW-n3", "rings-wall-mW-n3", "whole process"),
    )
    if any(by_id.get(left) is None or by_id.get(right) is None for left, right, _ in pairs):
        return []
    lines = [
        "",
        "% mW RINGS rows from the recorded n=3 samples.",
        "\\begin{table}[htbp]",
        "\\caption{Primitive rings on the 4096-molecule mW frame.",
        "n=3 recorded repetitions and no further warmup.",
        "The three repetitions ran at the same time, one pinned to each of",
        "cores 3, 4 and 5. The d-SEAMS repetitions used those same cores",
        "and finished before RINGS started.",
        "Each time is the median. The bracket is the minimum and the maximum.",
        "The ratio is the RINGS median divided by the d-SEAMS median.}",
        "\\label{tab:mw-rings-n3}",
        "\\centering",
        "\\scriptsize",
        "\\begin{tabular}{llrr}",
        "row & side & median [min, max] & core 3, 4, 5\\\\",
        "\\hline",
    ]
    for left_id, right_id, label in pairs:
        for sid, side in ((left_id, "d-SEAMS"), (right_id, "RINGS")):
            session = by_id[sid]
            stats = session_stats(session, session["columns"][0])
            values = column(session, session["columns"][0])
            cores = session.get("sample_cores") or []
            listed = ", ".join(
                f"{core}: {fmt(value)}" for core, value in zip(cores, values)
            )
            lines.append(
                " & ".join(
                    [
                        tex_escape(label),
                        tex_escape(side),
                        marked_cell(stats, session),
                        tex_escape(listed),
                    ]
                )
                + " \\\\"
            )
            summary.append(
                f"{sid} n={stats['n']} {marked_cell(stats, session)} cores {listed}"
            )
        left = by_id[left_id]
        right = by_id[right_id]
        left_stats = session_stats(left, left["columns"][0])
        right_stats = session_stats(right, right["columns"][0])
        shown_left = float(fmt(left_stats["median"]))
        shown_right = float(fmt(right_stats["median"]))
        shown_ratio = shown_right / shown_left
        exact_ratio = right_stats["median"] / left_stats["median"]
        lines.append(
            f"% {label} ratio of three-decimal medians "
            f"{fmt(right_stats['median'])} / {fmt(left_stats['median'])} "
            f"= {shown_ratio:.6g}; unrounded medians {exact_ratio:.6g}"
        )
        summary.append(
            f"ratio {label} {fmt(right_stats['median'])} / {fmt(left_stats['median'])} "
            f"= {shown_ratio:.6g} (unrounded {exact_ratio:.6g})"
        )
    lines += ["\\end{tabular}", "\\end{table}", ""]
    return lines


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: render_timings.py timings.json table2.tex quoted.tex", file=sys.stderr)
        return 2
    payload = json.loads(Path(sys.argv[1]).read_text())
    sessions = [s for s in payload["sessions"] if not s.get("superseded")]
    table_notes = write_table(sessions, Path(sys.argv[2]))
    quoted = write_quoted(sessions, Path(sys.argv[3]))
    print(quoted)
    print("--- table ---")
    print(table_notes)
    # Coverage check for the protocol size, from the first session that has samples.
    for session in sessions:
        if session["samples"]:
            stats = order_statistic_interval(column(session, session["columns"][0]))
            print(
                f"interval check n={stats['n']} k={stats['k']} "
                f"coverage={stats['coverage']:.10f}"
            )
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
