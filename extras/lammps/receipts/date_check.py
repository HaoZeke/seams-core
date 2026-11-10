# Copyright (c) d-SEAMS developers
"""Compare each commit date with the first log that cites that commit.

A log is a file whose text contains a start: timestamp. The check
passes when every cited commit is dated before that earliest start.
"""

import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

START = re.compile(r"^start:\s*(\S+)", re.M)
HASH = re.compile(r"\b[0-9a-f]{40}\b")


def commits(repo, base):
    text = subprocess.check_output(
        ["git", "-C", repo, "log", "--reverse", "--format=%H %aI", f"{base}..HEAD"],
        text=True,
    )
    rows = []
    for line in text.splitlines():
        digest, stamp = line.split()
        rows.append((digest, datetime.fromisoformat(stamp)))
    return rows


def logs(roots):
    found = []
    for root in roots:
        root = Path(root)
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix not in {".log", ".txt", ".md"} and path.name not in {
                "STATUS.txt", "TIP.txt",
            }:
                continue
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            match = START.search(text)
            if not match:
                continue
            found.append((path, datetime.fromisoformat(match.group(1)), text))
    return found


def main():
    seams = Path(os.environ["SRC_SEAMS"])
    lammps = Path(os.environ["SRC_LAMMPS"])
    roots = sys.argv[1:]
    series = [
        ("seams", commits(seams, "7f310225b7ea6c9e9e79d48d2b3547ef0f66610c")),
        ("lammps", commits(lammps, "098ffdefa525cabc82272925de02e9884f58c657")),
    ]
    cited = logs(roots)
    failed = 0
    for name, rows in series:
        print("series %s commits %d" % (name, len(rows)), flush=True)
        for digest, when in rows:
            hits = []
            for path, start, text in cited:
                if digest in text:
                    hits.append((start, path))
            if not hits:
                print("uncited %s %s %s" % (name, digest, when.isoformat()), flush=True)
                continue
            hits.sort(key=lambda item: item[0])
            first, path = hits[0]
            state = "before" if when < first else "after"
            if state != "before":
                failed += 1
            print(
                "%s %s %s commit %s log %s %s"
                % (state, name, digest, when.isoformat(), first.isoformat(), path),
                flush=True,
            )
    print("failures %d" % failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
