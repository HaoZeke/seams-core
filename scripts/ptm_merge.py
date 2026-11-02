#!/usr/bin/env python3
"""Build the PTM and CHILL+ sessions from an unedited ptm_rows.py log.

Uses paper_rows.session and paper_rows.merge. The samples, commands,
load and timestamps are the lines the log printed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import paper_rows
from timing_protocol import WARMUP


def blocks(text: str):
    lines = text.splitlines()
    out = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.endswith("_ms") and " " not in line:
            name = line[: -len("_ms")]
            i += 1
            vals = []
            while i < len(lines) and lines[i][:1].isdigit():
                vals.append([float(lines[i])])
                i += 1
            out[name] = vals
            continue
        i += 1
    return out


def field(text: str, key: str) -> str:
    for line in text.splitlines():
        if line.startswith(key + " "):
            return line[len(key) + 1 :]
    raise SystemExit(f"log has no {key}")


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: ptm_merge.py log timings.json")
    text = Path(sys.argv[1]).read_text()
    samples = blocks(text)
    agree = field(text, "agreement")
    started = field(text, "started")
    finished = field(text, "finished")
    load_start = field(text, "load_start")
    load_end = field(text, "load_end")

    def one(sid, side, code, command, key, comparison, frame):
        made = paper_rows.session(
            sid, comparison, frame, side, code, command.split(),
            samples[key], {}, agree,
        )
        made["started"] = started
        made["finished"] = finished
        made["load_start"] = load_start
        made["load_end"] = load_end
        made["warmup"] = WARMUP
        return made

    sessions = [
        one(
            "ptm-ic-dseams-0.1", "dseams", "d-SEAMS 2.0",
            field(text, "dseams_cmd"), "ptm-ic-dseams",
            "completed seeded cages vs PTM", "jittered Ic, 512",
        ),
        one(
            "ptm-ic-ovito-0.1", "ovito", "OVITO 3.15.5 PTM",
            field(text, "ovito_cmd"), "ptm-ic-ovito",
            "completed seeded cages vs PTM", "jittered Ic, 512",
        ),
        one(
            "chill-ic512-0.1", "dseams", "d-SEAMS 2.0 CHILL+",
            field(text, "chill_cmd"), "chill",
            "CHILL+ cutoff 3.5 A, list inside timer", "jittered Ic, 512",
        ),
    ]
    note = "rmsd_cutoff 0.15; replaced by the 0.1 session in logs/ptm-0.1.txt"
    paper_rows.JSON_PATH = Path(sys.argv[2])
    payload_path = paper_rows.JSON_PATH
    import json
    payload = json.loads(payload_path.read_text())
    for session in payload["sessions"]:
        if session.get("id") in ("ptm-ic-dseams", "ptm-ic-ovito"):
            session["superseded"] = note
    payload_path.write_text(json.dumps(payload, indent=2) + "\n")
    paper_rows.merge(sessions)


if __name__ == "__main__":
    main()
