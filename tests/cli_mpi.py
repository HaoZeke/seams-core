# Copyright (c) d-SEAMS developers
"""Exercise MPI branches of the seams command.

usage: cli_mpi.py MPIEXEC SEAMS TRAJ
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

mpiexec, seams, traj = sys.argv[1:4]


def run(args, env=None, n=2):
    cmd = [mpiexec, "-n", str(n), seams, *args]
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(cmd, capture_output=True, text=True, env=merged)


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(name + " " + detail)


def main():
    one = run(["rings", traj, "--frame", "1", "--last", "1"],
              env={"SEAMS_MPI": "0"})
    check("SEAMS_MPI=0", one.returncode == 0 and "rings " in one.stdout,
          one.stderr[-200:])

    split = run(["rings", traj, "--frame", "1", "--last", "2"])
    check("two ranks", split.returncode == 0 and "rings " in split.stdout,
          split.stderr[-240:])

    mismatch = subprocess.run(
        [mpiexec, "-n", "1", seams, "rings", traj, "--frame", "1", "--last", "1",
         ":", "-n", "1", seams, "read", traj, "--frame", "1", "--last", "1"],
        capture_output=True, text=True)
    err = mismatch.stderr
    check("mismatch", mismatch.returncode == 2 and "same command" in err, err[-300:])

    with tempfile.TemporaryDirectory() as tmp:
        per = str(Path(tmp) / "per.txt")
        both = run(["cages", traj, "--frame", "1", "--last", "1",
                    "--per-atom", per, "--graph", "cutoff", "--cutoff", "3.5"])
        check("per-atom", both.returncode == 2 and "one rank" in both.stderr,
              both.stderr[-300:])

    xyz = str(Path(traj).with_name("configuration.xyz"))
    if not Path(xyz).is_file():
        xyz = str(Path(traj).parent.parent / "xyz" / "configuration.xyz")
    other = run(["rings", xyz, "--frame", "1", "--last", "2"])
    check("xyz ranks", other.returncode == 0, other.stderr[-240:])

    axis = run(["density-z", traj, "--frame", "1", "--last", "1",
                "--axis", "2", "--bins", "2"], n=1)
    check("axis 2", axis.returncode == 0 and "# z rho" in axis.stdout,
          axis.stderr[-200:])

    held = run(["rings", traj, "--frame", "1", "--last", "2"],
               env={"OMP_NUM_THREADS": "2"})
    check("omp set", held.returncode == 0 and "rings " in held.stdout,
          held.stderr[-200:])

    # A rank variable with no world size still marks the process as launched.
    env = os.environ.copy()
    for key in list(env):
        if key.startswith(("OMPI_", "PMI_", "PMIX_", "MV2_")):
            del env[key]
    env.pop("SEAMS_MPI", None)
    env["PMI_RANK"] = "0"
    env["PMIX_RANK"] = "0"
    alone = subprocess.run(
        [seams, "rings", traj, "--frame", "1", "--last", "1"],
        capture_output=True, text=True, env=env)
    check("pmi rank", alone.returncode == 0 and "rings " in alone.stdout,
          alone.stderr[-300:])


if __name__ == "__main__":
    main()
