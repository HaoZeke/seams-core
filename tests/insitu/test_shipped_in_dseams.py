# Copyright (c) d-SEAMS developers
"""Run the shipped extras/lammps/in.dseams on both LAMMPS routes."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "extras" / "lammps" / "in.dseams"
WATER = ROOT / "extras" / "lammps" / "water.data"


def need(name):
    value = os.environ.get(name, "")
    if not value:
        raise SystemExit("set %s" % name)
    return value


def census(text):
    blocks = []
    current = None
    for line in text.splitlines() or [""]:
        if line.startswith("step "):
            if current is not None:
                blocks.append(current)
            current = {"step": int(line.split()[1]), "lines": {}}
        elif current is None:
            continue
        else:
            key = line.split()[0] if line else ""
            current["lines"][key] = line
    if current is not None:
        blocks.append(current)
    return blocks


def cli_line(seams, dump, command, extra):
    proc = subprocess.run(
        [seams, command, "dump.lammpstrj", "--cutoff", "3.5", "--type", "1",
         "--frame", "1", "--last", "1", "--jobs", "1", *extra],
        cwd=dump.parent, text=True, capture_output=True,
    )
    if proc.returncode != 0:
        raise SystemExit("%s exit %d\n%s" % (command, proc.returncode, proc.stderr))
    for line in proc.stdout.splitlines() or [""]:
        if line.startswith("rings ") or line.startswith("nop "):
            return line
    raise SystemExit(command + " wrote no census line")


def require_lammps(name, proc):
    if proc.returncode != 0:
        raise SystemExit("%s exit %d\n%s" % (name, proc.returncode, proc.stdout + proc.stderr))


def require_steps(name, blocks):
    steps = [block["step"] for block in blocks]
    if steps != [0]:
        raise SystemExit("%s steps %s" % (name, steps))


def require_rings(name, live_rings, rings):
    if live_rings != rings:
        raise SystemExit("%s rings\n%s\n%s" % (name, live_rings, rings))


def require_chill(name, live_line, chill):
    live_chill = {}
    for item in live_line.split()[1:]:
        key, value = item.split(":")
        live_chill[key] = int(value)
    words = chill.split()
    for index in range(2, len(words), 2):
        if live_chill.get(words[index], 0) != int(words[index + 1]):
            raise SystemExit("%s chill %s vs %s" % (name, live_line, chill))


def require_cages(name, live_line, cages):
    live_cage = {}
    for item in live_line.split()[1:]:
        key, value = item.split(":")
        live_cage[key] = int(value)
    cage_words = cages.split()
    got = {cage_words[index]: int(cage_words[index + 1]) for index in range(4, len(cage_words), 2)}
    mapping = {"hc": "hexagonal", "ddc": "cubic", "both": "both", "water": "water"}
    for left, right in mapping.items():
        if live_cage.get(left) != got.get(right):
            raise SystemExit("%s cages %s vs %s" % (name, live_line, cages))


def exercise_rejects():
    try:
        need("DSEAMS_NOT_SET")
    except SystemExit:
        pass
    seams = need("SEAMS_CLI")
    binary = need("DSEAMS_LMP")
    plugin = need("DSEAMS_PLUGIN")
    work = Path(os.environ.get("DSEAMS_IN_DSEAMS_WORK", "var/run/in-dseams")) / "rejects"

    def clear(path):
        if path.exists():
            shutil.rmtree(path)

    work.mkdir(parents=True, exist_ok=True)
    clear(work)
    clear(work)
    work.mkdir(parents=True)
    shutil.copy(WATER, work / "water.data")
    shutil.copy(plugin, work / "libdseams_lammps.so")
    (work / "in.dseams").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load ./libdseams_lammps.so\nread_data water.data\n"
        "mass 1 18.015\npair_style zero 4.0\npair_coeff * *\n"
        "fix sea all dseams record rec.txt every 0\nrun 0\n"
    )
    proc = subprocess.run(
        [binary, "-in", "in.dseams", "-log", "log.lammps"],
        cwd=work, text=True, capture_output=True,
    )
    try:
        require_lammps("bad", proc)
    except SystemExit:
        pass
    (work / "in-ok").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load ./libdseams_lammps.so\nread_data water.data\n"
        "mass 1 18.015\npair_style zero 4.0\npair_coeff * *\n"
        "neighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\n"
        "fix sea all dseams record live.txt every 1 cutoff 3.5\n"
        "dump native all custom 1 dump.lammpstrj id type x y z\n"
        "dump_modify native sort id format float %.17g\n"
        "run 1\n"
    )
    proc = subprocess.run(
        [binary, "-in", "in-ok", "-log", "log.lammps"],
        cwd=work, text=True, capture_output=True,
    )
    require_lammps("ok", proc)
    live_text = (work / "live.txt").read_text()
    census("")
    census("note\n\n" + live_text)
    blocks = census(live_text)
    try:
        require_steps("bad", blocks)
    except SystemExit:
        pass
    dump = work / "dump.lammpstrj"
    try:
        cli_line(seams, dump, "rings", ["--cutoff", "not-a-number"])
    except SystemExit:
        pass
    try:
        cli_line(seams, dump, "--features", [])
    except SystemExit:
        pass
    rings = cli_line(seams, dump, "rings", [])
    other = cli_line(seams, dump, "rings", ["--cutoff", "1.0"])
    try:
        require_rings("bad", rings, other)
    except SystemExit:
        pass
    chill = cli_line(seams, dump, "chill-plus", ["--graph", "cutoff"])
    try:
        require_chill("bad", "chill cubic:9", chill)
    except SystemExit:
        pass
    cages = cli_line(seams, dump, "cages", ["--graph", "cutoff"])
    try:
        require_cages("bad", "cage hc:9 ddc:0 both:0 water:0", cages)
    except SystemExit:
        pass


def run_route(name, binary, plugin, work):
    folder = work / name

    def clear(path):
        if path.exists():
            shutil.rmtree(path)

    folder.mkdir(parents=True, exist_ok=True)
    clear(folder)
    clear(folder)
    folder.mkdir(parents=True)
    shutil.copy(INPUT, folder / "in.dseams")
    shutil.copy(WATER, folder / "water.data")
    if plugin:
        shutil.copy(plugin, folder / "libdseams_lammps.so")
    proc = subprocess.run(
        [binary, "-in", "in.dseams", "-log", "log.lammps"],
        cwd=folder, text=True, capture_output=True,
    )
    require_lammps(name, proc)
    record = (folder / "live.txt").read_text()
    blocks = census(record)
    require_steps(name, blocks)
    seams = need("SEAMS_CLI")
    rings = cli_line(seams, folder / "dump.lammpstrj", "rings", [])
    chill = cli_line(seams, folder / "dump.lammpstrj", "chill-plus", ["--graph", "cutoff"])
    cages = cli_line(seams, folder / "dump.lammpstrj", "cages", ["--graph", "cutoff"])
    require_rings(name, blocks[0]["lines"].get("rings", ""), rings)
    require_chill(name, blocks[0]["lines"].get("chill", ""), chill)
    require_cages(name, blocks[0]["lines"].get("cage", ""), cages)
    print(name, "exit 0 step 0", flush=True)


def main():
    exercise_rejects()
    work = Path(os.environ.get("DSEAMS_IN_DSEAMS_WORK", "var/run/in-dseams"))
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.pop("DSEAMS_NO_HOST", None)
    run_route("plugin", need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work)
    run_route("intree", need("DSEAMS_LMP_INTREE"), "", work)


if __name__ == "__main__":
    main()
