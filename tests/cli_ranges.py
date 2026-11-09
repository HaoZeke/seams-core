# Copyright (c) d-SEAMS developers
"""Frame ranges, carried state, and extension fallbacks.

usage: cli_ranges.py SEAMS ROOT
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]
root = Path(sys.argv[2])
tiny = root / "input" / "traj" / "agi_water_tiny.lammpstrj"
water = root / "input" / "traj" / "exampleTraj.lammpstrj"


def run(args):
    return subprocess.run([seams, *args], capture_output=True, text=True)


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(name + " " + detail)


def expect(args, code, out_has=None, err_line=None):
    got = run(args)
    check("exit", got.returncode == code,
          " ".join(args[-2:]) + " rc " + str(got.returncode) + " " + got.stderr[:200])
    if out_has is not None:
        check("stdout " + out_has, out_has in got.stdout,
              "args " + " ".join(args[:4]) + " out " + got.stdout[:180])
    if err_line is not None:
        first = got.stderr.splitlines()[0] if got.stderr else ""
        check("stderr", first == err_line, first)
    return got


def main():
    body = tiny.read_text()
    _cleanup = tempfile.TemporaryDirectory()
    work = Path(_cleanup.name)
    two = work / "two.lammpstrj"
    two.write_text(body + body.replace("TIMESTEP\n0\n", "TIMESTEP\n1\n", 1))
    expect(["rings", str(two), "--frame", "1", "--last", "2", "--jobs", "1"],
           0, out_has="rings ")
    expect(["cages", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--graph", "cutoff", "--cutoff", "3.5"], 0, out_has="nop ")
    expect(["cages", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--graph", "knn", "-k", "3", "--rebuild"], 0, out_has="nop ")
    expect(["chill-plus", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--layers", "--axis", "2"], 0, out_has="nop ")
    expect(["fingerprint", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--hops", "2", "--format", "json"], 0, out_has='"command":"fingerprint"')
    expect(["hbonds", str(two), "--frame", "1", "--last", "2", "--type", "1",
            "--htype", "1", "--donors", "--hdist", "2.5", "--hangle", "30"],
           0, out_has="hbonds ")
    expect(["rdf", str(two), "--frame", "1", "--last", "2", "--bins", "4",
            "--types", "1, 1"], 0, out_has="# r g")
    expect(["density-z", str(two), "--frame", "1", "--last", "2", "--bins", "0",
            "--axis", "0"], 0, out_has="# z rho")
    expect(["rings", str(two), "--frame", "2", "--last", "2", "--strict-input"],
           0, out_has="rings ")

    for suffix in ("xyz", "pdb", "gro", "dcd", "lammpstrj", "con"):
        path = work / ("dir." + suffix)
        path.mkdir()
    expect(["rings", str(work / "dir.xyz"), "--frame", "1", "--last", "2"],
           2, err_line="XYZ header is not an atom count.")
    expect(["rings", str(work / "dir.pdb"), "--frame", "1", "--last", "2"],
           0, err_line="chemfiles cannot read " + str(work / "dir.pdb"))
    expect(["rings", str(work / "dir.gro"), "--frame", "1", "--last", "3"],
           0, err_line="chemfiles cannot read " + str(work / "dir.gro"))
    expect(["read", str(work / "dir.dcd"), "--frame", "1", "--last", "2"],
           0, err_line="chemfiles cannot read " + str(work / "dir.dcd"))
    expect(["rings", str(work / "dir.lammpstrj"), "--frame", "1", "--last", "2"], 0)
    expect(["read", str(work / "dir.con"), "--frame", "1", "--last", "1"],
           2, err_line="readcon-core is not in this build; CON input needs it")
    nodot = work / "nodot"
    nodot.write_text(body)
    expect(["rings", str(nodot), "--frame", "1", "--last", "1"], 0, out_has="rings ")

    expect(["pairs", str(water), "--site", "1=cationHead,2=anion",
            "--frame", "1", "--last", "1"], 0, out_has="contact-pair ")
    expect(["domains", str(water), "--site", "1=polar,2=apolar", "--cutoff", "6",
            "--subset", "polar", "--frame", "1", "--last", "1"], 0, out_has="subset ")
    expect(["domains", str(water), "--site", "1=tail,2=apolar", "--cutoff", "6",
            "--subset", "tail"], 0, out_has="subset ")
    expect(["cn", str(water), "--ions", "--site", "1=cationHead,2=anion",
            "--cutoff", "6"], 0, out_has="nCation ")
    if os.geteuid() != 0:
        locked = work / "locked.xyz"
        locked.write_text("1\ncomment\nO 0 0 0\n")
        locked.chmod(0)
        expect(["rings", str(locked), "--frame", "1", "--last", "2"], 0, out_has="rings ")
        expect(["rings", str(locked), "--frame", "2", "--last", "2"], 0, out_has="rings ")
        locked_pdb = work / "locked.pdb"
        locked_pdb.write_text("END\n")
        locked_pdb.chmod(0)
        expect(["read", str(locked_pdb), "--frame", "1", "--last", "3"], 0)
    empty_dump = work / "empty.lammpstrj"
    empty_dump.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
    )
    expect(["cages", str(empty_dump), "--frame", "1", "--last", "1",
            "--signature", "512", "--graph", "cutoff"], 0, out_has="cages 0")
    expect(["rings", str(empty_dump), "--frame", "1", "--last", "1", "--type", "0"],
           0, out_has="rings ")
    expect(["ions", str(tiny), "--ion-types", "3", "--ion-cutoff", "4",
            "--frame", "1", "--last", "1", "--complete"], 0, out_has="ions ")
    expect(["cages", str(ice_or_tiny()), "--frame", "1", "--last", "1",
            "--graph", "cutoff", "--signature", "512", "--shape",
            "--cage-edge", "0,1"], 0, out_has="nop ")


def ice_or_tiny():
    ice = root / "input" / "traj" / "genice_sI.lammpstrj"
    return ice if ice.is_file() else tiny


if __name__ == "__main__":
    main()
