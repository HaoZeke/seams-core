# Copyright (c) d-SEAMS developers
"""Import the built nanobind module and run its unittest. Exit 77 if it was not built."""

import os
import shutil
import sys
from pathlib import Path

try:
    import numpy  # noqa: F401
except ImportError:
    print("numpy is not installed")
    sys.exit(77)

build = Path(os.environ.get("DSEAMS_BUILD", "build"))
src = Path(os.environ.get("DSEAMS_SRC", Path(__file__).resolve().parents[2]))
# A conda python and the extension must share a stable ABI tag. A module
# built for another interpreter is a skip, not a failed import.
tag = sys.implementation.cache_tag  # cpython-312
# The extension is built at the top of the build tree. The staged installs
# below it hold copies that other tests delete and recreate.
cores = [path for path in build.glob("_core*.so") if tag in path.name]
if not cores:
    built = [path.name for path in build.glob("_core*.so")]
    print(f"dseams._core was not built for {tag}; found {built or 'nothing'}")
    sys.exit(77)

stage = build / "python-stage-unittest"
pkg = stage / "dseams"
pkg.mkdir(parents=True, exist_ok=True)
root = src / "bindings" / "python"


def place(source, dest):
    # A concurrent run may have dest mapped. Copying over it in place
    # crashes that run; a rename leaves its inode alone.
    part = dest.with_name(f".{dest.name}.{os.getpid()}")
    shutil.copy(source, part)
    os.replace(part, dest)


for name in ("__init__.py", "__main__.py"):
    place(root / "dseams" / name, pkg / name)
place(root / "frames.py", pkg / "frames.py")
for so in cores:
    dest = pkg / so.name
    if so.resolve() != dest.resolve():
        place(so, dest)

if os.environ.get("DSEAMS_PY_READY") != "1":
    env = os.environ.copy()
    lib = build / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    ctypes_src = root / "src"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(stage), str(ctypes_src), env.get("PYTHONPATH", "")]
    )
    env["DSEAMS_PY_READY"] = "1"
    os.execve(sys.executable, [sys.executable, __file__], env)

# sys.path[0] is bindings/python, whose source dseams package has no
# extension and would shadow the staged one.
sys.path[:] = [entry for entry in sys.path if Path(entry or os.curdir).resolve() != root.resolve()]

import unittest

suite = unittest.defaultTestLoader.discover(str(root / "tests"), pattern="test_*.py")
result = unittest.TextTestRunner(verbosity=1).run(suite)
sys.exit(0 if result.wasSuccessful() else 1)
