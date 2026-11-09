# Copyright (c) d-SEAMS developers
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


def lammps_case(binary, plugin, folder, ranks):
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    record = folder / "live.txt"
    dump = folder / "native.lammpstrj"
    script = folder / "in"
    # 216 waters, 40 steps, census and dump every 5 steps.
    live.lammps_input(script, 40, dump, record, plugin or "UNUSED", every=5, ncell=3, write_dump=True)
    text = script.read_text().replace(
        "dump_modify native sort id format float %.17g",
        "dump_modify native sort id delay 5 format float %.17g",
    )
    if not plugin:
        text = "".join(line for line in text.splitlines(True) if "plugin load" not in line)
    script.write_text(text)
    env = lib_env(str(Path(binary).parent))
    mpirun = os.environ.get("DSEAMS_MPIRUN", "")
    if ranks > 1:
        if not mpirun:
            raise SystemExit("set DSEAMS_MPIRUN")
        cmd = [mpirun, "--oversubscribe", "-np", str(ranks), binary, "-in", "in", "-log", "log"]
    else:
        cmd = [binary, "-in", "in", "-log", "log"]
    live.run(cmd, folder, env=env)
    matched = live.match_lammps_cli(record, dump)
    (folder / "equal.log").write_text("frames %d ranks %d\n" % (matched, ranks))
    if matched != 8:
        raise RuntimeError("expected 8 frames in %s, got %d" % (folder, matched))
    return matched


def gromacs_case():
    binary = need("DSEAMS_GMX")
    mpirun = need("DSEAMS_GMX_MPIRUN")
    folder = live.WORK / "e2e-gromacs"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    rows = live.diamond_positions(3)
    live.write_gro(folder / "conf.gro", rows, 1.92)
    live.write_oxygen_top(folder / "topol.top", len(rows))
    (folder / "md.mdp").write_text(
        "integrator = md\n"
        "dt = 0.001\n"
        "nsteps = 40\n"
        "nstxout = 5\n"
        "nstvout = 0\n"
        "nstfout = 0\n"
        "nstlog = 5\n"
        "nstlist = 10\n"
        "cutoff-scheme = Verlet\n"
        "coulombtype = Cut-off\n"
        "rcoulomb = 0.9\n"
        "rvdw = 0.9\n"
        "pbc = xyz\n"
        "constraints = none\n"
        "gen-vel = yes\n"
        "gen-temp = 30\n"
        "gen-seed = 12345\n"
        "dseams-active = yes\n"
        "dseams-input-group = System\n"
        "dseams-nstcalc = 5\n"
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
    trr = folder / "run.trr"
    frames = live.trr_frames(trr)
    # seams reads the TRR through chemfiles. The atoms are oxygen, so
    # the type filter is atomic number 8.
    dump = folder / "from-trr.lammpstrj"
    live.trr_to_dump(frames, dump)
    matched = live.match_lammps_cli(folder / "live.txt", dump, type_filter=1)
    trr_rings = [
        ln
        for ln in live.cli_text("rings", trr, type_filter=1).splitlines()
        if ln.startswith("rings ")
    ]
    dump_rings = [
        ln
        for ln in live.cli_text("rings", dump, type_filter=1).splitlines()
        if ln.startswith("rings ")
    ]
    live.require_equal(trr_rings, dump_rings, "seams on the TRR does not match seams on the same frames")
    native = live.match_native_cli(folder / "live.txt", trr, type_filter=1)
    live.require_native_frames(native, matched, len(frames))
    (folder / "equal.log").write_text("frames %d\n" % matched)
    receipt = folder / "receipt.txt"
    receipt.write_text(
        "seams %s\nlib %s\ngmx %s\nmpirun %s\n"
        % (live.SEAMS, live.LIB, binary, mpirun)
    )
    return matched


def main():
    need("SEAMS_LIB")
    need("SEAMS_CLI")
    work = Path(os.environ.get("DSEAMS_E2E_WORK", "var/run/dseams-e2e"))
    live.WORK = work
    print("openmm", live.openmm_identity(), flush=True)
    print("lammps-plugin-1", lammps_case(need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work / "lammps-plugin-1", 1), flush=True)
    print("lammps-plugin-2", lammps_case(need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work / "lammps-plugin-2", 2), flush=True)
    print("lammps-package", lammps_case(need("DSEAMS_LMP_INTREE"), "", work / "lammps-package", 1), flush=True)
    print("gromacs", gromacs_case(), flush=True)


if __name__ == "__main__": main()
