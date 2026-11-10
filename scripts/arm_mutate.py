#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Mutate one source string, rebuild, run a test, and restore the file.

The receipt is the raw compiler and test output. A closed arm fails
this test after the mutation and passes before it.

usage: arm_mutate.py FILE OLD NEW NINJA_DIR TARGETS... -- TEST_CMD...
"""

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if "--" not in sys.argv:
        sys.stderr.write(
            "usage: arm_mutate.py FILE OLD NEW NINJA_DIR TARGET -- TEST...\n"
        )
        return 2
    cut = sys.argv.index("--")
    file_s, old, new, ninja_dir = sys.argv[1:5]
    targets = sys.argv[5:cut]
    test = sys.argv[cut + 1 :]
    path = Path(file_s)
    original = path.read_text()
    if original.count(old) != 1:
        sys.stderr.write("old string count %d\n" % original.count(old))
        return 2
    path.write_text(original.replace(old, new, 1))
    try:
        build = subprocess.run(
            ["ninja", "-C", ninja_dir, "-j1", *targets],
            text=True,
            capture_output=True,
        )
        sys.stdout.write(build.stdout)
        sys.stdout.write(build.stderr)
        sys.stdout.write("NINJA:%d\n" % build.returncode)
        if build.returncode != 0:
            return build.returncode
        ran = subprocess.run(test, text=True, capture_output=True)
        sys.stdout.write(ran.stdout)
        sys.stdout.write(ran.stderr)
        sys.stdout.write("TEST:%d\n" % ran.returncode)
        failed = ran.returncode != 0
        return 0 if failed else 1
    finally:
        path.write_text(original)
        subprocess.run(["ninja", "-C", ninja_dir, "-j1", *targets], check=False)


if __name__ == "__main__":
    raise SystemExit(main())
