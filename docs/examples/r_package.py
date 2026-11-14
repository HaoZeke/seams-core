# Copyright (c) d-SEAMS developers
"""Install bindings/r/dseams into a library directory owned by one test.

meson runs tests in parallel. A library shared between tests lets one test's
R CMD INSTALL remove dseams while another test loads it, so each test gets
build/r-lib/<test> and reinstalls there on every run.
"""

import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def install(rscript, build, test, langdir, env):
    """Install dseams for `test` and point env's R_LIBS_USER at it first.

    rscript is the Rscript executable; R sits next to it. langdir holds
    libseams_lang. Returns the library directory.
    """
    dest = Path(build) / "r-lib" / test
    shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir(parents=True)
    prior = env.get("R_LIBS_USER", "")
    env["R_LIBS_USER"] = os.pathsep.join([str(dest)] + ([prior] if prior else []))
    install_env = dict(env)
    install_env["DSEAMS_ROOT"] = str(ROOT)
    install_env["SEAMS_INCLUDE"] = str(ROOT / "bindings" / "cxx")
    install_env["SEAMS_LIBDIR"] = str(langdir)
    install_env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(langdir), install_env.get("LD_LIBRARY_PATH", "")]
    )
    r_cmd = str(Path(rscript).with_name("R"))
    done = subprocess.run(
        [r_cmd, "CMD", "INSTALL", "-l", str(dest), str(ROOT / "bindings" / "r" / "dseams")],
        env=install_env,
        capture_output=True,
        text=True,
        check=False,
    )
    if done.returncode != 0:
        raise RuntimeError("R CMD INSTALL for {}:\n{}{}".format(test, done.stdout[-2000:], done.stderr[-2000:]))
    return dest
