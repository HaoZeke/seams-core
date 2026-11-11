# Copyright (c) d-SEAMS developers
"""Classify every CHILL+ label that disagrees with OVITO on the mW liquid.

Reads the three bond dumps written by explain_chill and the OVITO
structure-type column. Does not change a threshold.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

HERE = Path(__file__).resolve().parent / "records" / "explain"
OVITO = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial-hydrate",
}
MERGE = {
    "cubic": "cubic",
    "reCubic": "cubic",
    "hexagonal": "hexagonal",
    "reHex": "hexagonal",
    "water": "other",
    "unclassified": "other",
    "interfacial": "interfacial",
    "clathrate": "hydrate",
    "interClathrate": "interfacial-hydrate",
}


def load(path: Path):
    rows = []
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        rows.append({
            "index": int(parts[0]),
            "shell": int(parts[1]),
            "nbonds": int(parts[2]),
            "stag": int(parts[3]),
            "ecl": int(parts[4]),
            "other": int(parts[5]),
            "label": parts[6],
            "cij": tuple(round(float(value), 8) for value in parts[7:]),
        })
    return rows


def merged(row):
    return MERGE[row["label"]]


def classify(left, right):
    """Shell gate on main is a bug. The rest is the published q3 rule."""
    if left["shell"] != 4:
        return (
            "shell_gate",
            "BUG",
            "src/bop.cpp:1118",
        )
    return (
        "q3_definition",
        "DEFINITION",
        "e8610637:src/include/internal/chill_plus_scheme.hpp:28-36",
    )


def main():
    main_rows = load(HERE / "main.tsv")
    published = load(HERE / "published.tsv")
    ovito_rule = load(HERE / "ovito-rule.tsv")
    ovito = [
        OVITO[int(line)]
        for line in (HERE.parent / "a-chill" / "mW-liquid-T215-ovito.txt").read_text().splitlines()
        if line.strip()
    ]
    lines = [
        "index shell comparison left right step verdict source"
    ]
    for kind, rows in (("main", main_rows), ("published", published)):
        for left, right, ov in zip(rows, ovito_rule, ovito):
            if merged(left) == ov:
                continue
            which, verdict, source = classify(left, right)
            # Published code already applies the shell gate, so its
            # remaining misses are the q3 definition, not the main bug.
            if kind == "published":
                which, verdict, source = (
                    "q3_definition",
                    "DEFINITION",
                    "e8610637:src/include/internal/chill_plus_scheme.hpp:28-36",
                )
            lines.append(
                f"{left['index']} {left['shell']} {kind}-vs-ovito "
                f"{merged(left)} {ov} {which} {verdict} {source}"
            )
    out = HERE / "liquid-steps.tsv"
    out.write_text("\n".join(lines) + "\n")
    print(out.name, len(lines) - 1)


if __name__ == "__main__":
    main()
