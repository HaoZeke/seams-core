# Copyright (c) d-SEAMS developers
"""In-situ OpenMM output matches seams on the frames the reporter wrote."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "adapters" / "python"))

from dseams_insitu import DseamsReporter, Engine  # noqa: E402

from water import build_simulation  # noqa: E402

SEAMS = Path(os.environ.get("SEAMS", ROOT / "bbdir-eessi" / "seams"))
LIB = Path(
    os.environ.get("SEAMS_LIB", ROOT / "bbdir-eessi" / "bindings" / "libseams_lang.so")
)
CUTOFF = "3.5"
TYPE = "1"


def seams(dump, command, frame):
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.run(
        [
            str(SEAMS),
            command,
            str(dump),
            "--cutoff",
            CUTOFF,
            "--type",
            TYPE,
            "--graph",
            "cutoff",
            "--frame",
            str(frame),
            "--jobs",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return proc.stdout.strip()


def blocks(text):
    groups = []
    current = []
    for line in text.splitlines():
        if line.startswith("step ") and current:
            groups.append(current)
            current = [line]
        else:
            current.append(line)
    if current and any(line.strip() for line in current):
        groups.append(current)
    return groups


def section(lines, head, tail=None):
    start = lines.index(head)
    if tail is None:
        return lines[start]
    stop = lines.index(tail, start)
    return "\n".join(lines[start : stop + 1])


def match_box(box_nm, steps, interval, min_waters):
    simulation, oxygens = build_simulation(box_nm)
    n_mol = simulation.topology.getNumResidues()
    if n_mol < min_waters:
        raise SystemExit(f"expected at least {min_waters} waters, got {n_mol}")
    with tempfile.TemporaryDirectory(prefix="dseams-openmm-") as tmp:
        dump = str(Path(tmp) / "frames.lammpstrj")
        record = str(Path(tmp) / "census.txt")
        with Engine(str(LIB), record, cutoff=float(CUTOFF)) as engine:
            simulation.reporters.append(
                DseamsReporter(engine, dump, interval=interval, indices=oxygens)
            )
            simulation.step(steps)
        parsed = blocks(Path(record).read_text(encoding="ascii"))
        expect = steps // interval
        if len(parsed) != expect:
            raise SystemExit(f"{n_mol} waters: {len(parsed)} frames, expected {expect}")
        for frame, lines in enumerate(parsed, start=1):
            live = {
                "neighbours": section(lines, "neighbours", "end_neighbours"),
                "rings": section(lines, next(line for line in lines if line.startswith("rings "))),
                "chill": next(line for line in lines if line.startswith("chill ")).split(" ", 1)[1],
                "chill-plus": next(
                    line for line in lines if line.startswith("chill-plus ")
                ).split(" ", 1)[1],
            }
            offline = {
                "neighbours": seams(dump, "neighbours", frame),
                "rings": seams(dump, "rings", frame),
                "chill": seams(dump, "chill", frame),
                "chill-plus": seams(dump, "chill-plus", frame),
            }
            for name in offline:
                if live[name] != offline[name]:
                    raise SystemExit(
                        f"{n_mol} waters frame {frame} {name} differs\n"
                        f"in-situ:\n{live[name]}\n"
                        f"seams:\n{offline[name]}"
                    )
        print(
            f"matched rings and chill-plus on {len(parsed)} frames, "
            f"{n_mol} TIP4P/Ew waters"
        )


def test_insitu_matches_seams():
    if not SEAMS.is_file():
        raise SystemExit("seams binary is missing: " + str(SEAMS))
    if not LIB.is_file():
        raise SystemExit("libseams_lang is missing: " + str(LIB))
    os.environ["OMP_NUM_THREADS"] = "1"
    match_box(2.0, steps=100, interval=50, min_waters=200)
    match_box(4.93, steps=200, interval=100, min_waters=3000)


if __name__ == "__main__":
    test_insitu_matches_seams()
