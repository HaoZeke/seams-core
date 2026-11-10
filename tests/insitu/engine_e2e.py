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
    chem = Path(live.LIB).resolve().parent.parent / "chemfiles-prefix" / "lib"
    parts = [extra, str(Path(live.LIB).parent), str(chem), merged.get("LD_LIBRARY_PATH", "")]
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
    text = script.read_text()
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
    if matched != 9:
        raise RuntimeError("expected 9 frames in %s, got %d" % (folder, matched))
    return matched


def main():
    need("SEAMS_LIB")
    need("SEAMS_CLI")
    work = Path(os.environ.get("DSEAMS_E2E_WORK", "var/run/dseams-e2e"))
    live.WORK = work
    print("lammps-plugin-1", lammps_case(need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work / "lammps-plugin-1", 1), flush=True)
    print("lammps-plugin-2", lammps_case(need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work / "lammps-plugin-2", 2), flush=True)
    print("lammps-package", lammps_case(need("DSEAMS_LMP_INTREE"), "", work / "lammps-package", 1), flush=True)


if __name__ == "__main__":
    main()
