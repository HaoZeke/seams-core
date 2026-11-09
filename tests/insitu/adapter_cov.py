# Copyright (c) d-SEAMS developers
"""Build the LAMMPS adapter with coverage and exercise it.

Tool paths come from the environment. Sources are taken from this repo.
"""

import os
import subprocess
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def need(name):
    value = os.environ.get(name, "")
    if not value:
        raise SystemExit("set %s" % name)
    return value


def run(cmd, cwd, env=None):
    proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stdout[-1500:] + proc.stderr[-1500:])
    return proc


def lammps_cases(work, plugin, lmp, mpirun, libdir):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = libdir + ":" + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    header = textwrap.dedent(
        """\
        units real
        atom_style atomic
        boundary p p p
        lattice diamond 6.4
        region box block 0 2 0 2 0 2
        create_box 2 box
        create_atoms 1 box
        create_atoms 2 single 6.0 6.0 6.0
        mass * 18.0
        pair_style lj/cut 8.0
        pair_coeff * * 0.05 2.5
        neighbor 2.0 bin
        timestep 2.0
        velocity all set 0 0 0
        plugin load %s
        """
        % plugin
    )
    cases = {
        "bad-every": header + "fix sea all dseams record rec.txt every 0\nrun 1\n",
        "bad-cutoff": header + "fix sea all dseams record rec.txt cutoff 0\nrun 1\n",
        "bad-depth": header + "fix sea all dseams record rec.txt depth 0\nrun 1\n",
        "bad-double": header + "fix sea all dseams record rec.txt cutoff nope\nrun 1\n",
        "bad-word": header + "fix sea all dseams record rec.txt nope 1\nrun 1\n",
        "no-record": header + "fix sea all dseams every 1 cutoff 3.5\nrun 1\n",
        "bad-record": header + "fix sea all dseams record /no/such/dir/rec.txt every 1\nrun 1\n",
        "live": header
        + "fix hold all setforce 0 0 0\n"
        + "fix sea all dseams dump dump.lammpstrj record rec.txt every 1 cutoff 3.2 depth 4\n"
        + "run 2\n",
        "move": header
        + "velocity all create 30.0 7\nfix int all nve\n"
        + "fix sea oxy dseams record rec.txt every 1 cutoff 3.5\n"
        + "group oxy type 1\n"
        + "run 2\n",
    }
    # group must exist before the fix that names it
    cases["move"] = header + textwrap.dedent(
        """\
        group oxy type 1
        velocity all set 0 0 0
        fix hold all setforce 0.0 0.0 0.0
        fix sea oxy dseams record rec.txt every 1 cutoff 3.5
        run 1
        displace_atoms all move 1.0 0.0 0.0
        run 1
        """
    )
    for name, text in cases.items():
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in").write_text(text)
        which = env
        # The still run uses the host list. A second still run skips it.
        run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log"], folder, which)
        print("lammps", name, "ok", flush=True)
    both = work / "ranks"
    both.mkdir(parents=True, exist_ok=True)
    (both / "in").write_text(
        header
        + "fix sea all dseams record rec.txt every 2 cutoff 3.5\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "2", lmp, "-in", "in", "-log", "log"], both, env)
    print("lammps ranks ok", flush=True)
    still = work / "nohost"
    still.mkdir(parents=True, exist_ok=True)
    (still / "in").write_text(
        header + "fix sea all dseams record rec.txt every 1 cutoff 3.5\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log"], still, env_host)
    print("lammps nohost ok", flush=True)



def sys_python():
    return os.environ.get("DSEAMS_PYTHON", "python3")


def main():
    lmp_src = need("DSEAMS_LAMMPS_SRC")
    seams_lib = need("SEAMS_LIB")
    cxx = need("DSEAMS_MPICXX")
    lmp = need("DSEAMS_LMP")
    mpirun = need("DSEAMS_MPIRUN")
    work = Path(os.environ.get("DSEAMS_ADAPTER_WORK", "var/run/dseams-adapter-cov"))
    work.mkdir(parents=True, exist_ok=True)
    build = work / "plugin"
    if not (build / "CMakeCache.txt").exists():
        run(
            [
                "cmake", "-S", str(ROOT / "src" / "adapters" / "lammps"), "-B", str(build),
                "-DLAMMPS_SRC=" + lmp_src,
                "-DSEAMS_INCLUDE=" + str(ROOT / "src" / "include" / "internal"),
                "-DSEAMS_LIBRARY=" + seams_lib,
                "-DCMAKE_BUILD_TYPE=Debug",
                "-DCMAKE_CXX_COMPILER=" + cxx,
                "-DCMAKE_CXX_FLAGS=--coverage -O0 -fprofile-update=single",
            ],
            work,
        )
    run(["cmake", "--build", str(build), "-j", "4"], work)
    plugin = str(build / "libdseams_lammps.so")
    libdir = str(Path(seams_lib).parent)
    lammps_cases(work / "lmp", plugin, lmp, mpirun, libdir + ":" + str(Path(seams_lib).parent.parent / "subprojects" / "readcon-core"))
    print(work)


if __name__ == "__main__":
    main()
