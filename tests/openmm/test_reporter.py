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


def test_insitu_matches_seams():
    if not SEAMS.is_file():
        raise SystemExit("seams binary is missing: " + str(SEAMS))
    if not LIB.is_file():
        raise SystemExit("libseams_lang is missing: " + str(LIB))
    os.environ["OMP_NUM_THREADS"] = "1"
    simulation, oxygens = build_simulation()
    n_mol = simulation.topology.getNumResidues()
    if n_mol < 200:
        raise SystemExit(f"expected a few hundred waters, got {n_mol}")
    with tempfile.TemporaryDirectory(prefix="dseams-openmm-") as tmp:
        dump = str(Path(tmp) / "frames.lammpstrj")
        record = str(Path(tmp) / "census.txt")
        with Engine(str(LIB), record, cutoff=float(CUTOFF)) as engine:
            simulation.reporters.append(
                DseamsReporter(engine, dump, interval=50, indices=oxygens)
            )
            simulation.step(100)
        text = Path(record).read_text(encoding="ascii")
        parsed = blocks(text)
        if len(parsed) != 2:
            raise SystemExit(f"expected 2 analysed steps, got {len(parsed)}")
        for frame, lines in enumerate(parsed, start=1):
            got_neighbours = section(lines, "neighbours", "end_neighbours")
            got_rings = section(lines, next(l for l in lines if l.startswith("rings ")))
            got_chill = next(l for l in lines if l.startswith("chill ")).split(" ", 1)[1]
            got_plus = next(l for l in lines if l.startswith("chill-plus ")).split(
                " ", 1
            )[1]
            offline = {
                "neighbours": seams(dump, "neighbours", frame),
                "rings": seams(dump, "rings", frame),
                "chill": seams(dump, "chill", frame),
                "chill-plus": seams(dump, "chill-plus", frame),
            }
            live = {
                "neighbours": got_neighbours,
                "rings": got_rings,
                "chill": got_chill,
                "chill-plus": got_plus,
            }
            for name in offline:
                if live[name] != offline[name]:
                    raise SystemExit(
                        f"frame {frame} {name} differs\n"
                        f"in-situ:\n{live[name]}\n"
                        f"seams:\n{offline[name]}"
                    )
        print(f"matched {len(parsed)} frames, {n_mol} TIP4P/Ew waters")


if __name__ == "__main__":
    test_insitu_matches_seams()
