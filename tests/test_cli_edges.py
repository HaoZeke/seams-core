# Copyright (c) d-SEAMS developers
"""Exercise seams error paths on real arguments.

usage: test_cli_edges.py SEAMS TRAJECTORY
"""

import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]
traj = sys.argv[2]


def run(args):
    return subprocess.run([seams, *args], capture_output=True, text=True)


def check(name, cond):
    if not cond:
        raise SystemExit(name)


def expect(args, code, fragment):
    got = run(args)
    check(" ".join(args) + " exit", got.returncode == code)
    if fragment:
        check(" ".join(args) + " message", fragment in got.stderr + got.stdout)


def main():
    expect([], 2, "required")
    expect(["no-such-command", traj], 2, "unknown command")
    expect(["--format", "yaml", "read", traj], 2, "bad --format")
    expect(["--axis", "q", "density-z", traj], 2, "bad --axis")
    expect(["ions", traj], 2, "ions needs --ion-types")
    expect(["pairs", traj], 2, "pairs needs --site")
    expect(["domains", traj], 2, "domains needs --site")
    expect(["cn", "--ions", traj], 2, "cn --ions needs --site")
    expect(["--subset", "sideways", "--site", "1=polar", "domains", traj], 2, "bad --subset")
    expect(["--cutoff", "not-a-number", "rings", traj], 2, "")
    expect(["--frame", "1", "--last", "3", "read", traj + ".missing-suffix"], 2, "")
    with tempfile.TemporaryDirectory() as tmp:
        empty = Path(tmp) / "empty.xyz"
        empty.write_text("")
        expect(["read", str(empty)], 2, "XYZ read failed: the file has no frames")
        short = Path(tmp) / "short.lammpstrj"
        short.write_text(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        )
        expect(
            ["read", str(short)],
            2,
            "LAMMPS dump frame is truncated: expected 3 atoms, read 1.",
        )
    expect(["--features"], 0, "OpenMP")
    expect(["read", traj], 0, "nop ")


if __name__ == "__main__":
    main()
