"""Build the LAMMPS and PLUMED adapters with coverage and exercise them.

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
    env_host = env.copy()
    env_host["DSEAMS_NO_HOST"] = "1"
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


def plumed_cases(work, so, plumed_lib, libdir, plumed_inc, gxx):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(Path(so).parent) + ":" + libdir + ":" + str(Path(plumed_lib).parent) + ":" + env.get("LD_LIBRARY_PATH", "")
    src = work / "drive_plumed.cpp"
    src.write_text(r'''
#include "wrapper/Plumed.h"
#include <exception>
#include <iostream>
#include <string>
#include <vector>
int once(const std::string &so, const std::vector<std::string> &lines, const std::vector<int> &steps, const double *box) {
  try {
    PLMD::Plumed p;
    int nat = 4;
    double dt = 0.001;
    p.cmd("setNatoms", &nat);
    p.cmd("setTimestep", &dt);
    p.cmd("init");
    for (const auto &line : lines) {
      p.cmd("readInputLine", line.c_str());
    }
    double xyz[12] = {0,0,0, 0.12,0,0, 0,0.12,0, 0.04,0.04,0.12};
    double masses[4] = {1,1,1,1};
    double forces[12] = {};
    double virial[9] = {};
    for (int step : steps) {
      p.cmd("setStep", &step);
      p.cmd("setPositions", xyz);
      p.cmd("setMasses", masses);
      p.cmd("setForces", forces);
      p.cmd("setBox", const_cast<double *>(box));
      p.cmd("setVirial", virial);
      p.cmd("calc");
    }
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
}
int main(int argc, char **argv) {
  const std::string load = std::string("LOAD FILE=") + argv[1];
  const double good[9] = {2,0,0, 0,2,0, 0,0,2};
  const double bad[9] = {0,0,0, 0,0,0, 0,0,0};
  once(argv[1], {load, "d: DSEAMS ATOMS=1-4 RECORD=rec.txt CUTOFF=3.5 EVERY=2"}, {0, 1}, good);
  once(argv[1], {load, "d: DSEAMS ATOMS=1-4 RECORD=rec.txt DUMP=dump.lammpstrj CUTOFF=3.5 EVERY=2"}, {0, 1}, bad);
  once(argv[1], {load, "d: DSEAMS ATOMS=1-4 RECORD=/no/such/dir/rec.txt CUTOFF=3.5"}, {0}, good);
  once(argv[1], {load, "d: DSEAMS ATOMS= RECORD=rec.txt CUTOFF=3.5"}, {0}, good);
  once(argv[1], {load, "d: DSEAMS ATOMS=1-4 RECORD=rec.txt CUTOFF=-1"}, {0}, good);
  once(argv[1], {load, "d: DSEAMS ATOMS=1-4 RECORD=rec.txt EVERY=0"}, {0}, good);
  std::cout << "plumed driver ok\n";
  return 0;
}
''')
    exe = work / "drive_plumed"
    run(
        [gxx, "-std=c++17", "-O0", "-g", "-I" + plumed_inc, str(src), "-o", str(exe),
         "-L" + str(Path(plumed_lib).parent), "-lplumed",
         "-Wl,-rpath," + str(Path(plumed_lib).parent)],
        work, env,
    )
    run([str(exe), so], work, env)


def sys_python():
    return os.environ.get("DSEAMS_PYTHON", "python3")


def main():
    lmp_src = need("DSEAMS_LAMMPS_SRC")
    seams_lib = need("SEAMS_LIB")
    cxx = need("DSEAMS_MPICXX")
    lmp = need("DSEAMS_LMP")
    mpirun = need("DSEAMS_MPIRUN")
    plumed_inc = need("DSEAMS_PLUMED_INC")
    plumed_lib = need("DSEAMS_PLUMED_LIB")
    gxx = need("DSEAMS_GXX")
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
    so = work / "libDseams.so"
    run(
        [
            gxx, "-shared", "-fPIC", "-std=c++17", "-O0", "-g", "--coverage",
            "-fprofile-update=single",
            "-I" + plumed_inc, "-I" + str(ROOT / "src" / "include" / "internal"),
            str(ROOT / "src" / "adapters" / "plumed" / "Dseams.cpp"),
            "-o", str(so),
            "-L" + libdir, "-lyodaLib",
            "-L" + str(Path(plumed_lib).parent), "-lplumed",
            "-Wl,-rpath," + libdir,
            "-Wl,-rpath," + str(Path(plumed_lib).parent),
        ],
        work,
    )
    plumed_cases(work, str(so), plumed_lib, libdir, plumed_inc, gxx)
    print(work)


if __name__ == "__main__":
    main()
