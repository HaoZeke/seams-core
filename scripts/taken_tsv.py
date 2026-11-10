#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write TAKEN.tsv from one merged gcovr JSON summary.

The branch bar is branches taken. gcovr is run with
--exclude-throw-branches. Branches executed is not this table.

usage: taken_tsv.py SUMMARY.json TAKEN.tsv
"""

import json
import sys


def main() -> int:
    if len(sys.argv) != 3:
        sys.stderr.write("usage: taken_tsv.py SUMMARY.json TAKEN.tsv\n")
        return 2
    data = json.load(open(sys.argv[1]))
    rows = ["file\ttaken\ttotal\ttaken_percent\tlines_percent\tbar"]
    for item in data.get("files", []):
        name = item["filename"]
        if not name.startswith("src/") or not name.endswith(".cpp"):
            continue
        taken = item.get("branch_covered") or 0
        total = item.get("branch_total") or 0
        percent = item.get("branch_percent")
        if percent is None:
            percent = -1.0
        lines = float(item.get("line_percent") or 0)
        if not total:
            bar = "no-branches"
        elif percent < 98:
            bar = "below"
        else:
            bar = "met"
        rows.append(
            "%s\t%s\t%s\t%.2f\t%.2f\t%s" % (name, taken, total, percent, lines, bar)
        )
    open(sys.argv[2], "w").write("\n".join(rows) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
