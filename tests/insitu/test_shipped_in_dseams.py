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
        raise SystemExit("%s exit %d\n%s" % (command, proc.returncode, proc.stderr))
    for line in proc.stdout.splitlines():
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


class Finished(object):
    def __init__(self, code, out="", err=""):
        self.returncode = code
        self.stdout = out
        self.stderr = err


def assert_reject(fn, message):
    try:
        fn()
    except SystemExit as exc:
        if str(exc) != message:
            raise SystemExit("reject was %r" % str(exc))
        return
    raise SystemExit("reject did not fire")


def exercise_rejects():
    assert_reject(lambda: need("DSEAMS_NOT_SET"), "set DSEAMS_NOT_SET")
    assert_reject(lambda: cli_line("/bin/false", Path("."), "rings", []), "rings exit 1\n")
    assert_reject(lambda: cli_line("/bin/true", Path("."), "rings", []), "rings wrote no census line")
    assert_reject(lambda: require_lammps("plugin", Finished(1, "boom")), "plugin exit 1\nboom")
    assert_reject(lambda: require_steps("plugin", [{"step": 1}]), "plugin steps [1]")
    assert_reject(lambda: require_rings("plugin", "rings 3:1", "rings 3:0"), "plugin rings\nrings 3:1\nrings 3:0")
    assert_reject(
        lambda: require_chill("plugin", "chill cubic:1", "nop 1 cubic 0"),
        "plugin chill chill cubic:1 vs nop 1 cubic 0",
    )
    assert_reject(
        lambda: require_cages("plugin", "cage hc:1 ddc:0 both:0 water:0", "nop 1 x 0 hexagonal 0 cubic 0 both 0 water 0"),
        "plugin cages cage hc:1 ddc:0 both:0 water:0 vs nop 1 x 0 hexagonal 0 cubic 0 both 0 water 0",
    )
    blocks = census("step 0\nrings a\nstep 1\nrings b\n")
    if [block["step"] for block in blocks] != [0, 1]:
        raise SystemExit("two steps")


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
    (work / "plugin").mkdir(parents=True, exist_ok=True)
    run_route("plugin", need("DSEAMS_LMP"), need("DSEAMS_PLUGIN"), work)
    run_route("intree", need("DSEAMS_LMP_INTREE"), "", work)


if __name__ == "__main__":
    main()
