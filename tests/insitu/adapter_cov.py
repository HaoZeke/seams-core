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


def command_text(proc, cwd):
    log_path = Path(cwd) / "log.lammps"
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    return (proc.stdout or "") + (proc.stderr or "") + log


def run(cmd, cwd, env=None):
    proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        raise SystemExit("exit %d want 0\n%s" % (proc.returncode, command_text(proc, cwd)[-1500:]))
    return proc


def expect(cmd, cwd, env, code, message):
    proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True)
    text = command_text(proc, cwd)
    if proc.returncode != code:
        raise SystemExit("exit %d want %d\n%s" % (proc.returncode, code, text[-1500:]))
    if message not in text:
        raise SystemExit("missing %s\n%s" % (message, text[-1500:]))
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
    rejects = {
        "bad-every": (header + "fix sea all dseams record rec.txt every 0\nrun 1\n", 1, "fix dseams: every must be a positive integer"),
        "bad-cutoff": (header + "fix sea all dseams record rec.txt cutoff 0\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "bad-depth": (header + "fix sea all dseams record rec.txt depth 0\nrun 1\n", 1, "fix dseams: depth must be a positive integer"),
        "bad-double": (header + "fix sea all dseams record rec.txt cutoff nope\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "bad-word": (header + "fix sea all dseams record rec.txt nope 1\nrun 1\n", 1, "Illegal fix dseams command"),
        "bad-host": (header + "fix sea all dseams record rec.txt host maybe\nrun 1\n", 1, "fix dseams: host is yes or no"),
        "no-record": (header + "fix sea all dseams every 1 cutoff 3.5\nrun 1\n", 1, "fix dseams: record is required"),
        "bad-record": (header + "fix sea all dseams record /no/such/dir/rec.txt every 1\nrun 1\n", 1, "fix dseams: record did not open"),
    }
    cases = {
        "live": header
        + "fix hold all setforce 0 0 0\n"
        + "fix sea all dseams dump dump.lammpstrj record rec.txt every 1 cutoff 3.2 depth 4 host yes\n"
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
    for name, (text, code, message) in rejects.items():
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in").write_text(text)
        expect(
            [mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"],
            folder, env, code, message,
        )
        print("lammps", name, "exit", code, flush=True)
    for name, text in cases.items():
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in").write_text(text)
        run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], folder, env)
        print("lammps", name, "exit 0", flush=True)
    both = work / "ranks"
    both.mkdir(parents=True, exist_ok=True)
    (both / "in").write_text(
        header
        + "fix sea all dseams record rec.txt every 2 cutoff 3.5\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "2", lmp, "-in", "in", "-log", "log.lammps"], both, env)
    print("lammps ranks ok", flush=True)
    still = work / "nohost"
    still.mkdir(parents=True, exist_ok=True)
    (still / "in").write_text(
        header + "fix sea all dseams record rec.txt every 1 cutoff 3.5 host no\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], still, env)
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
