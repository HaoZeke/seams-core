"""Run each engine and compare the live census with postprocessing.

The binaries and the library come from the environment.
A missing binary is an error.
"""

import os
import shutil
from pathlib import Path

import live


def need(name):
    value = os.environ.get(name, "")
    if not value:
        raise SystemExit("set %s" % name)
    return value


def lib_env(extra):
    merged = os.environ.copy()
    parts = [extra, str(Path(live.LIB).parent), merged.get("LD_LIBRARY_PATH", "")]
    merged["LD_LIBRARY_PATH"] = os.pathsep.join(part for part in parts if part)
    merged["OMP_NUM_THREADS"] = "1"
    return merged


def lammps_case(binary, plugin, folder):
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    record = folder / "live.txt"
    dump = folder / "native.lammpstrj"
    script = folder / "in"
    live.lammps_input(script, 4, dump, record, plugin or "UNUSED", every=1, ncell=2, write_dump=True)
    if not plugin:
        kept = [line for line in script.read_text().splitlines(True) if "plugin load" not in line]
        script.write_text("".join(kept))
    live.run([binary, "-in", "in", "-log", "log"], folder, env=lib_env(str(Path(binary).parent)))
    matched = live.match_lammps_cli(record, dump)
    (folder / "equal.log").write_text("frames %d\n" % matched)
    if matched < 1:
        raise RuntimeError("no matched frames in %s" % folder)
    return matched


def gromacs_case():
    binary = need("DSEAMS_GMX")
    mpirun = need("DSEAMS_GMX_MPIRUN")
    folder = live.WORK / "e2e-gromacs"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    rows = live.diamond_positions(2)
    live.write_gro(folder / "conf.gro", rows, 1.28)
    live.write_oxygen_top(folder / "topol.top", len(rows))
    (folder / "md.mdp").write_text(
        "integrator = md\n"
        "dt = 0.001\n"
        "nsteps = 4\n"
        "nstxout = 1\n"
        "nstvout = 0\n"
        "nstfout = 0\n"
        "nstlog = 1\n"
        "nstlist = 10\n"
        "cutoff-scheme = Verlet\n"
        "coulombtype = Cut-off\n"
        "rcoulomb = 0.6\n"
        "rvdw = 0.6\n"
        "pbc = xyz\n"
        "constraints = none\n"
        "gen-vel = yes\n"
        "gen-temp = 30\n"
        "gen-seed = 12345\n"
        "dseams-active = yes\n"
        "dseams-input-group = System\n"
        "dseams-nstcalc = 1\n"
        "dseams-cutoff = 3.5\n"
        "dseams-output = live.txt\n"
    )
    env = lib_env(str(Path(binary).resolve().parent.parent / "lib"))
    live.run(
        [mpirun, "--oversubscribe", "-np", "1", binary, "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "topol.tpr", "-maxwarn", "2"],
        folder,
        env=env,
    )
    live.run(
        [mpirun, "--oversubscribe", "-np", "1", binary, "mdrun", "-s", "topol.tpr", "-deffnm", "run", "-noconfout"],
        folder,
        env=env,
    )
    frames = live.trr_frames(folder / "run.trr")
    dump = folder / "from-trr.lammpstrj"
    live.trr_to_dump(frames, dump)
    matched = live.match_lammps_cli(folder / "live.txt", dump)
    (folder / "equal.log").write_text("frames %d\n" % matched)
    if matched < 1:
        raise RuntimeError("gromacs matched no frames")
    return matched


def main():
    need("SEAMS_LIB")
    need("SEAMS_CLI")
    work = Path(os.environ.get("DSEAMS_E2E_WORK", "var/run/dseams-e2e"))
    live.WORK = work
    print("openmm", live.openmm_identity(), flush=True)
    print("lammps-plugin", lammps_case(need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work / "lammps-plugin"), flush=True)
    print("lammps-package", lammps_case(need("DSEAMS_LMP_INTREE"), "", work / "lammps-package"), flush=True)
    print("gromacs", gromacs_case(), flush=True)


if __name__ == "__main__":
    main()
