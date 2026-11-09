"""seams exits 2 on an empty xyz and reads a multi-frame xyz.

usage: test_xyz_cli.py SEAMS
"""

import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]


def run(args):
    return subprocess.run([seams, *args], capture_output=True, text=True)


def check(name, cond):
    if not cond:
        raise SystemExit(name)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        empty = root / "empty.xyz"
        empty.write_text("")
        got = run(["read", str(empty)])
        check("empty exit", got.returncode == 2)
        check("empty message", got.stderr.strip() == "XYZ read failed: the file has no frames")
        check("empty frame", "nop " not in got.stdout)

        multi = root / "multi.xyz"
        multi.write_text("1\nframe\nO 1 2 3\n1\nframe\nO 4 5 6\n")
        feat = run(["--features"])
        got = run(["--frame", "1", "--last", "2", "--jobs", "1", "read", str(multi)])
        check("multi exit", got.returncode == 0)
        want = 2 if "chemfiles: enabled" in feat.stdout else 1
        check("multi frames", got.stdout.count("nop 1 frame ") == want)


if __name__ == "__main__":
    main()
