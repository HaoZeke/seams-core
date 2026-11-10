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
    for line in text.splitlines():
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
        raise SystemExit(command + " failed\n" + proc.stderr)
    for line in proc.stdout.splitlines():
        if line.startswith("rings ") or line.startswith("nop "):
            return line
    raise SystemExit(command + " wrote no census line")


def run_route(name, binary, plugin, work):
    folder = work / name
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    shutil.copy(INPUT, folder / "in.dseams")
    shutil.copy(WATER, folder / "water.data")
    if plugin:
        shutil.copy(plugin, folder / "libdseams_lammps.so")
    proc = subprocess.run(
        [binary, "-in", "in.dseams", "-log", "log.lammps"],
        cwd=folder, text=True, capture_output=True,
    )
    if proc.returncode != 0:
        raise SystemExit("%s exit %d\n%s" % (name, proc.returncode, proc.stdout + proc.stderr))
    record = (folder / "live.txt").read_text()
    blocks = census(record)
    if [block["step"] for block in blocks] != [0]:
        raise SystemExit("%s steps %s" % (name, [block["step"] for block in blocks]))
    seams = need("SEAMS_CLI")
    rings = cli_line(seams, folder / "dump.lammpstrj", "rings", [])
    chill = cli_line(seams, folder / "dump.lammpstrj", "chill-plus", ["--graph", "cutoff"])
    cages = cli_line(seams, folder / "dump.lammpstrj", "cages", ["--graph", "cutoff"])
    live_rings = blocks[0]["lines"].get("rings", "")
    if live_rings != rings:
        raise SystemExit("%s rings\n%s\n%s" % (name, live_rings, rings))
    live_chill = {}
    for item in blocks[0]["lines"].get("chill", "").split()[1:]:
        key, value = item.split(":")
        live_chill[key] = int(value)
    words = chill.split()
    for index in range(2, len(words), 2):
        if live_chill.get(words[index], 0) != int(words[index + 1]):
            raise SystemExit("%s chill %s vs %s" % (name, blocks[0]["lines"].get("chill"), chill))
    live_cage = {}
    for item in blocks[0]["lines"].get("cage", "").split()[1:]:
        key, value = item.split(":")
        live_cage[key] = int(value)
    cage_words = cages.split()
    got = {cage_words[index]: int(cage_words[index + 1]) for index in range(4, len(cage_words), 2)}
    mapping = {"hc": "hexagonal", "ddc": "cubic", "both": "both", "water": "water"}
    for left, right in mapping.items():
        if live_cage.get(left) != got.get(right):
            raise SystemExit("%s cages %s vs %s" % (name, blocks[0]["lines"].get("cage"), cages))
    print(name, "exit 0 step 0", flush=True)


def main():
    work = Path(os.environ.get("DSEAMS_IN_DSEAMS_WORK", "var/run/in-dseams"))
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.pop("DSEAMS_NO_HOST", None)
    run_route("plugin", need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work)
    run_route("intree", need("DSEAMS_LMP_INTREE"), "", work)


if __name__ == "__main__":
    main()
