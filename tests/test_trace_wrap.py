#!/usr/bin/env python3
"""TRACE misses the sI cages whose bonds wrap the cell.

On the 46-molecule GenIce sI cell, d-SEAMS finds 2 of 5^12 and 6 of
5^12 6^2. TRACE returns one 5^12 and no 5^12 6^2. The seven cages it
drops are the rows in benchmarks/h2h-clavis-cage-ids.csv, and each of
them has a bond that crosses the periodic boundary. The 5^12 both codes
return does not.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from h2h_clavis_cages import (  # noqa: E402
    dseams_vertices,
    read_ortho,
    trace_parse,
    trace_payload,
    wrapping_bonds,
)

IDS = ROOT / "benchmarks" / "h2h-clavis-cage-ids.csv"


def listed_cages():
    found = {}
    with IDS.open(encoding="ascii", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["phase"] != "sI" or row["side"] != "dseams":
                continue
            vertices = tuple(int(item) for item in row["vertices"].split())
            found.setdefault(row["signature"], set()).add(vertices)
            assert row["cage_id"]
    return found


def test_wrap_miss():
    import clavis_engine as engine

    pos, lengths = read_ortho(ROOT / "input" / "traj" / "genice_sI.lammpstrj")
    assert len(pos) == 46
    ours = dseams_vertices(pos, lengths)
    assert len(ours["512"]) == 2
    assert len(ours["51262"]) == 6
    assert len(ours["51264"]) == 0
    text, box_nm, cutoff_nm, _note = trace_payload(pos, lengths)
    raw = engine.analyze_waters(text, "[]", box_nm, "trace", cutoff_nm, 35.0)
    theirs = trace_parse(raw)
    assert len(theirs["512"]) == 1
    assert theirs["51262"] == set()
    assert theirs["51264"] == set()
    missed = {
        "512": ours["512"] - theirs["512"],
        "51262": ours["51262"] - theirs["51262"],
    }
    assert missed == listed_cages()
    assert len(missed["512"]) == 1
    assert len(missed["51262"]) == 6
    for cages in missed.values():
        for cage in cages:
            wrap, bonds = wrapping_bonds(cage, pos, lengths)
            assert bonds > 0
            assert wrap > 0
    shared = ours["512"] & theirs["512"]
    assert len(shared) == 1
    wrap, bonds = wrapping_bonds(next(iter(shared)), pos, lengths)
    assert bonds == 30
    assert wrap == 0


if __name__ == "__main__":
    test_wrap_miss()
    print("ok")
