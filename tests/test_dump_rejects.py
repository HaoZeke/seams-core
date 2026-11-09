"""seams exits 2 and names a short xyz, a junk header, a bad coordinate, and a short dump.

usage: test_dump_rejects.py SEAMS
"""

import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]


def run(path):
    return subprocess.run([seams, "read", str(path)], capture_output=True, text=True)


def check(name, cond):
    if not cond:
        raise SystemExit(name)


def expect(path, text, message):
    got = run(path)
    check(text + " exit", got.returncode == 2)
    check(text + " message", got.stderr.strip() == message)
    check(text + " frame", "nop " not in got.stdout)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        short = root / "short.xyz"
        short.write_text("3\ncomment\nO 0 0 0\n")
        expect(short, "short", "XYZ frame is truncated: expected 3 atoms, read 1.")
        junk = root / "junk.xyz"
        junk.write_text("this is not an xyz file\n")
        expect(junk, "junk", "XYZ header is not an atom count.")
        bad = root / "bad.xyz"
        bad.write_text("1\ncomment\nO a b c\n")
        expect(bad, "bad", "XYZ coordinate is not a number.")
        dump = root / "short.lammpstrj"
        dump.write_text(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        )
        expect(dump, "dump", "LAMMPS dump frame is truncated: expected 3 atoms, read 1.")


if __name__ == "__main__":
    main()
