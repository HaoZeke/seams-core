#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Install the R package from this build and run its testthat file.

usage: run_r_testthat.py ROOT BUILD

The package goes into BUILD/r-test-lib on every run, so the test never
sees an install from an older checkout. R_LIBS_USER from the caller stays
on the library path after it, for cpp11 and testthat.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_r_testthat.py ROOT BUILD")
    root = Path(sys.argv[1]).resolve()
    build = Path(sys.argv[2]).resolve()
    script = shutil.which("Rscript")
    if script is None:
        print("Rscript is not on PATH")
        return 77
    dest = build / "r-test-lib"
    dest.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(root)
    env["SEAMS_INCLUDE"] = str(root / "bindings" / "cxx")
    env["SEAMS_LIBDIR"] = str(build / "bindings")
    libs = [str(build / "bindings"), str(build / "src")]
    if env.get("LD_LIBRARY_PATH"):
        libs.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(libs)
    user = env.get("R_LIBS_USER", "")
    env["R_LIBS_USER"] = os.pathsep.join(p for p in (str(dest), user) if p)
    env.setdefault("OMP_NUM_THREADS", "2")
    r_cmd = str(Path(script).with_name("R"))
    subprocess.check_call(
        [r_cmd, "CMD", "INSTALL", "-l", str(dest), str(root / "bindings" / "r" / "dseams")],
        env=env,
    )
    test = root / "bindings" / "r" / "dseams" / "tests" / "testthat" / "test-ice.R"
    return subprocess.call(
        [script, "-e", 'testthat::test_file("{}", stop_on_failure = TRUE)'.format(test)],
        env=env,
    )


if __name__ == "__main__":
    raise SystemExit(main())
