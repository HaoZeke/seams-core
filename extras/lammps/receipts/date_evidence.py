# Copyright (c) d-SEAMS developers
"""Pair each rewritten commit with the clock reading that dated it.

The readings come from the rewrite driver. This script does not invent
a timestamp and does not add an offset to a start time. A commit whose
author date differs from its recorded reading is an error.
"""

import os
import subprocess
import sys


def git(repo, *args):
    return subprocess.check_output(["git", "-C", repo, *args], text=True)


def commits(repo, base, tip):
    text = git(repo, "rev-list", "--reverse", "%s..%s" % (base, tip))
    return [line for line in text.split() if line]


def main():
    clock_path, seams, seams_base, lammps, lammps_base = sys.argv[1:6]
    rows = []
    current = None
    for line in open(clock_path, encoding="utf-8"):
        parts = line.rstrip("\n").split(" ", 1)
        kind = parts[0]
        rest = parts[1] if len(parts) > 1 else ""
        if kind == "repo":
            current = {"repo": rest, "reads": []}
        elif kind == "read" and current is not None:
            current["reads"].append(rest)
        elif kind == "use" and current is not None:
            current["use"] = rest
        elif kind == "old" and current is not None:
            current["old"] = rest
        elif kind == "new" and current is not None:
            current["new"] = rest
        elif kind == "subject" and current is not None:
            current["subject"] = rest
            rows.append(current)
            current = {"repo": rows[-1]["repo"], "reads": []}
    print("commit clock author committer subject", flush=True)
    seen = set()
    for row in rows:
        repo = seams if row["repo"] == seams else lammps
        author = git(repo, "log", "-1", "--format=%aI", row["new"]).strip()
        committer = git(repo, "log", "-1", "--format=%cI", row["new"]).strip()
        signed = git(repo, "log", "-1", "--format=%G?", row["new"]).strip()
        if signed != "N":
            raise SystemExit("signed %s" % row["new"])
        if author != committer:
            raise SystemExit("dates differ %s" % row["new"])
        pretty = git(repo, "log", "-1", "--format=%ad", "--date=format:%Y-%m-%d %H:%M:%S %z", row["new"]).strip()
        if pretty != row["use"]:
            raise SystemExit("clock %s stored %s" % (row["use"], pretty))
        if row["use"] in seen:
            raise SystemExit("shared second %s" % row["use"])
        seen.add(row["use"])
        print("commit %s" % row["new"], flush=True)
        print("clock %s" % row["use"], flush=True)
        print("reads %d" % len(row["reads"]), flush=True)
        print("author %s" % author, flush=True)
        print("committer %s" % committer, flush=True)
        print("subject %s" % row["subject"], flush=True)
        print("---", flush=True)
    seam_n = len(commits(seams, seams_base, git(seams, "rev-parse", "HEAD").strip()))
    lammps_n = len(commits(lammps, lammps_base, git(lammps, "rev-parse", "HEAD").strip()))
    if len(rows) != seam_n + lammps_n:
        raise SystemExit("rows %d commits %d" % (len(rows), seam_n + lammps_n))
    print("commits %d" % len(rows), flush=True)
    print("exit 0", flush=True)


if __name__ == "__main__":
    main()
