# Copyright (c) d-SEAMS developers
"""Install bindings/r/dseams into a library directory owned by one test.

meson runs tests in parallel, and `meson test --repeat` runs one test beside
itself. A shared library lets one R CMD INSTALL remove dseams while another
process loads it, so each process gets build/r-lib/<test>-<pid>, installs
there, and removes it on exit.
"""

import atexit
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _alive(pid):
    """A killed run leaves its library behind; its pid no longer exists."""
    try:
        os.kill(int(pid), 0)
    except (ValueError, ProcessLookupError):
        return False
    except PermissionError:
        return True
    return True


def install(rscript, build, test, langdir, env):
    """Install dseams for `test` and point env's R_LIBS_USER at it first.

    rscript is the Rscript executable; R sits next to it. langdir holds
    libseams_lang. Returns the library directory.
    """
    lib = Path(build) / "r-lib"
    for stale in lib.glob(test + "-*"):
        if not _alive(stale.name.rsplit("-", 1)[1]):
            shutil.rmtree(stale, ignore_errors=True)
    dest = lib / "{}-{}".format(test, os.getpid())
    shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir(parents=True)
    atexit.register(shutil.rmtree, dest, True)
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
