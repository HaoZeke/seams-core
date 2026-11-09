# Copyright (c) d-SEAMS developers
"""Cover the remaining seams command paths on real inputs.

usage: cli_cover.py SEAMS ROOT
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
mw = root / "input" / "traj" / "mW_cubic.lammpstrj"
nucleation = root / "input" / "traj" / "nucleation_f1.lammpstrj"


def run(args, env=None):
    return subprocess.run([seams, *args], capture_output=True, text=True, env=env)


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(name + " " + detail)


def expect(args, code, out_has=None, err_line=None, env=None):
    got = run(args, env)
    check("exit " + " ".join(args[:3]), got.returncode == code,
          "rc " + str(got.returncode) + " err " + got.stderr[:240]
          + " out " + got.stdout[:180])
    if out_has is not None:
        check("stdout " + out_has, out_has in got.stdout, got.stdout[:300])
    if err_line is not None:
        first = got.stderr.splitlines()[0] if got.stderr else ""
        check("stderr " + err_line, first == err_line, first)
    return got


def dump(path, n, rows, step=0):
    path.write_text(
        "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id type x y z\n%s" % (step, n, rows)
    )


def main():
    _cleanup = tempfile.TemporaryDirectory()
    work = Path(_cleanup.name)
    one = ["--frame", "1", "--last", "1", "--jobs", "1"]

    empty = work / "empty.lammpstrj"
    dump(empty, 0, "")
    expect(["rdf", str(empty), *one, "--type", "0", "--bins", "4"], 0, out_has="# types 1 ")
    expect(["cn", str(empty), *one, "--type", "0", "--bins", "4"], 0, out_has="cn ")
    expect(["cages", str(empty), *one, "--cutoff", "3.0"], 0, out_has="both 0")

    lone = work / "lone.lammpstrj"
    dump(lone, 1, "1 1 1 1 1\n")
    expect(["cn", str(lone), "--ions", "--site", "1=cationHead",
            "--cutoff", "6", *one], 0, out_has="nCation 1")
    anion = work / "anion.lammpstrj"
    dump(anion, 1, "1 2 1 1 1\n")
    expect(["cn", str(anion), "--ions", "--site", "2=anion",
            "--cutoff", "6", *one], 0, out_has="nAnion 1")
    expect(["cn", str(water), "--ions", "--site", "1=", *one], 2,
           err_line="site spec token '1=' is not key=value")

    gap = work / "gap.lammpstrj"
    dump(gap, 1, "1 1 1 1 1\n", 0)
    gap.write_text(gap.read_text() + (
        "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id type x y z\n"
    ))
    expect(["rings", str(gap), "--frame", "1", "--last", "2", "--strict-input"], 2,
           out_has="input frame 2 contains no atoms")

    xyz_dir = work / "fallback.xyz"
    xyz_dir.mkdir()
    expect(["rings", str(xyz_dir), "--frame", "2", "--last", "3"], 2,
           err_line="frame 2 is past the last frame in this file")
    pdb_dir = work / "fallback.pdb"
    pdb_dir.mkdir()
    expect(["read", str(pdb_dir), "--frame", "1", "--last", "2"], 0,
           out_has="nop 0", err_line="chemfiles cannot read " + str(pdb_dir))

    expect(["hbonds", str(water), *one, "--type", "2", "--htype", "1",
            "--cutoff", "3.5"], 0, out_has="hbonds ")
    expect(["f4", str(mw), *one, "--type", "1", "--htype", "2", "--cutoff", "3.5"],
           0, out_has="f4 ")

    expect(["cages", str(mw), *one, "--cutoff", "3.0", "--graph", "cutoff",
            "--water-types", "1", "--per-atom", str(work / "mw-cage.txt")],
           0, out_has="cubic ")
    expect(["cages", str(nucleation), *one, "--cutoff", "3.0", "--graph", "cutoff",
            "--water-types", "1", "--per-atom", str(work / "nuc-cage.txt")],
           0, out_has="both ")
    expect(["cages", str(tiny), *one, "--cutoff", "3.5", "--signature", "512",
            "--cage-edge", "none", "--graph", "knn", "-k", "4"], 0, out_has="nop ")
    expect(["cages", str(mw), *one, "--cutoff", "3.0", "--graph", "cutoff",
            "--signature", "512", "--shape", "--water-types", "1",
            "--cage-edge", "2.5,3.2"], 0, out_has="nop ")

    expect(["ions", str(tiny), *one, "--type", "1", "--ion-types", "1,3",
            "--ion-cutoff", "5", "--cutoff", "3.5"], 0, out_has="ions ")
    expect(["ions", str(nucleation), *one, "--cutoff", "3.0", "--ion-types", "1",
            "--ion-cutoff", "3.5", "--water-types", "1",
            "--per-atom", str(work / "nuc-ions.txt")], 0, out_has="ions ")

    expect(["rings", str(water), *one, "--water-types", "1,,2"], 0, out_has="rings ")
    expect(["cages", str(tiny), *one, "--guest-types", ",3,", "--signature", "512"],
           0, out_has="nop ")
    expect(["ions", str(tiny), *one, "--ion-types", "3,,4"], 0, out_has="ions ")

    lib_a = expect(["fingerprint", str(tiny), *one, "--emit-library", "tiny",
                    "--cutoff", "3.5", "--hops", "2"], 0)
    lib_b = expect(["fingerprint", str(water), *one, "--emit-library", "water",
                    "--cutoff", "3.5", "--type", "2", "--hops", "1"], 0)
    path_a = work / "a.lib"
    path_b = work / "b.lib"
    path_a.write_text(lib_a.stdout)
    path_b.write_text(lib_b.stdout)
    per = work / "labels.txt"
    expect(["fingerprint", str(tiny), *one, "--cutoff", "3.5", "--hops", "1",
            "--library", str(path_b), "--per-atom", str(per)], 0, out_has="unmatched")
    expect(["fingerprint", str(tiny), *one, "--cutoff", "3.5", "--hops", "2",
            "--library", str(path_a) + "," + str(path_b)], 0, out_has="depth ")
    expect(["fingerprint", str(tiny), *one, "--cutoff", "3.5", "--by-atom"],
           0, out_has="atom ")

    env = dict(os.environ)
    env["SEAMS_CARRY"] = "0"
    two = work / "two.lammpstrj"
    body = tiny.read_text()
    two.write_text(body + "\n" + body.replace("TIMESTEP\n0\n", "TIMESTEP\n1\n", 1))
    expect(["cages", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--cutoff", "3.5"], 0, out_has="nop ", env=env)
    expect(["rings", str(two), "--frame", "1", "--last", "2", "--jobs", "2"],
           0, out_has="rings ")

    expect(["chill", str(mw), *one, "--cutoff", "3.0", "--graph", "cutoff",
            "--water-types", "1"], 0, out_has="cubic ")
    expect(["chill", str(mw), *one, "--cutoff", "3.0", "--graph", "mutual",
            "--water-types", "1"], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), "--frame", "1", "--last", "1", "--cutoff", "3.5",
            "--graph", "mutual"], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), *one, "--cutoff", "3.5", "--graph", ""],
           0, out_has="graph cutoff")
    expect(["chill", str(mw), *one, "--cutoff", "3.0", "--graph", "",
            "--water-types", "1"], 0, out_has="nop ")

    xyz = root / "input" / "xyz" / "configuration.xyz"
    expect(["rings", str(xyz), "--frame", "1", "--last", "2", "--jobs", "2"],
           0, out_has="rings ")

    three = work / "three.lammpstrj"
    three.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 1 1 1 1\n"
        "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
        "ITEM: TIMESTEP\n2\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
    )
    expect(["rings", str(three), "--frame", "1", "--last", "3", "--strict-input"],
           2, out_has="input frame 2 contains no atoms")

    expect(["ions", str(tiny), *one, "--ion-types", ","], 0, out_has="ions 0")

    dimer = work / "dimer.lammpstrj"
    # O-H 0.96 A, H-O-H about 104 deg, one donated H between the oxygens.
    dimer.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n6\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id mol type x y z\n"
        "1 1 2 0.00 0.00 0.00\n"
        "2 1 1 0.96 0.00 0.00\n"
        "3 1 1 -0.24 0.93 0.00\n"
        "4 2 2 2.80 0.00 0.00\n"
        "5 2 1 1.84 0.00 0.00\n"
        "6 2 1 3.50 0.75 0.00\n"
    )
    expect(["hbonds", str(dimer), *one, "--type", "2", "--htype", "1",
            "--cutoff", "3.5", "--hdist", "2.5", "--hangle", "30"],
           0, out_has="nop 2 hbonds 1")
    expect(["f4", str(dimer), *one, "--type", "2", "--htype", "1", "--cutoff", "3.5"],
           0, out_has="nop 6 f4 1 n 2")
    expect(["hbonds", str(water), *one, "--type", "2", "--htype", "1",
            "--hdist", "2.5", "--hangle", "30"], 0, out_has="nop 250 hbonds 471")
    expect(["f4", str(water), *one, "--type", "2", "--htype", "1"],
           0, out_has="nop 750 f4 0.221753 n 250")

    expect(["fingerprint", str(mw), *one, "--cutoff", "3.0", "--hops", "2",
            "--type", "1", "--library", str(path_a) + "," + str(path_b)],
           0, out_has="depth")

    shuf = work / "shuf.lammpstrj"
    shuf.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n6\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id mol type x y z\n"
        "6 2 1 3.50 0.75 0.00\n"
        "2 1 1 0.96 0.00 0.00\n"
        "5 2 2 2.80 0.00 0.00\n"
        "1 1 2 0.00 0.00 0.00\n"
        "3 1 1 -0.24 0.93 0.00\n"
        "4 2 1 1.84 0.00 0.00\n"
    )
    per = work / "shuf-per.txt"
    expect(["cages", str(shuf), *one, "--graph", "cutoff", "--per-atom", str(per)],
           0, out_has="nop 2")
    text = per.read_text()
    check("shuf id 1", "1 2 0 0 0" in text, text)
    check("shuf id 5", "5 2 2.8" in text or "5 2 2.80" in text, text)
    expect(["hbonds", str(shuf), *one, "--type", "2", "--htype", "1",
            "--hdist", "2.5", "--hangle", "30"], 0, out_has="nop 2 hbonds 1")
    gap = work / "gap.dump"
    gap.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id mol type x y z\n"
        "4 7 2 1 1 1\n"
        "1 7 1 0 0 0\n"
    )
    gper = work / "gap-per.txt"
    expect(["cages", str(gap), *one, "--graph", "cutoff", "--per-atom", str(gper)],
           0, out_has="nop 1")
    check("gap id 4", "4 2 1 1 1" in gper.read_text(), gper.read_text())
    stray = work / "stray.lammpstrj"
    dump(stray, 2, "1 1 1 1 1\n2 1 4 1 1\n")
    expect(["fingerprint", str(stray), *one, "--cutoff", "3.5", "--hops", "2",
            "--type", "1", "--library", str(path_a) + "," + str(path_b)],
           0, out_has="depth")

    zero = work / "zero-box.lammpstrj"
    zero.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 0\n0 0\n0 0\n"
        "ITEM: ATOMS id mol type x y z\n"
        "1 1 1 0 0 0\n2 2 2 0 0 0\n"
    )
    expect(["cn", str(zero), "--ions", "--site", "1=cationHead,2=anion",
            "--cutoff", "6", *one], 0, out_has="degree 0",
           err_line="The point cloud has no simulation box; cannot build a neighbour list.")

    short = work / "short-box.lammpstrj"
    short.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 1 1 1 1\n"
    )
    expect(["density-z", str(short), *one, "--bins", "0"], 0, out_has="# z rho")

    nodot = work / "nodot"
    nodot.write_text(tiny.read_text())
    expect(["rings", str(nodot), "--frame", "1", "--last", "2"], 0, out_has="rings ")
    expect(["rings"], 2, err_line="A command and a trajectory file are required")

    ice = expect(["chill", str(nucleation), *one, "--cutoff", "3.5"],
                 0, out_has="reCubic 5")
    check("reHex 33", "reHex 33" in ice.stdout, ice.stdout[:240])

    one_xyz = work / "one.xyz"
    one_xyz.write_text("1\ncomment\nO 0 0 0\n")
    expect(["rings", str(one_xyz), "--frame", "2", "--last", "2"], 2,
           err_line="frame 2 is past the last frame in this file")
    expect(["rings", str(tiny), "--frame", "5", "--last", "6"], 2,
           err_line="frame 5 is past the last frame in this file")

    empty = work / "empty.lammpstrj"
    empty.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
    )
    expect(["f4", str(empty), *one, "--type", "1", "--htype", "2"],
           0, out_has="nop 0 f4 nan")
    expect(["ions", str(empty), *one, "--ion-types", "3"], 0, out_has="nop 0")
    expect(["read", str(empty), *one, "--strict-input"], 2,
           out_has="contains no atoms")
    expect(["rings", str(tiny), *one, "--family", "moltenSalt"], 0, out_has="rings ")

    ions = work / "ions.lammpstrj"
    ions.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
        "1 1 0 0 0\n"
        "2 2 8 0 0\n"
    )
    expect(["cn", str(ions), "--ions", "--site", "1=cationHead,2=anion",
            "--cutoff", "6", *one], 0, out_has="nAnion 1")

    junk = work / "junk.pdb"
    junk.write_text("END\n")
    expect(["read", str(junk), "--frame", "1", "--last", "2"], 0, out_has="nop ")

    noid = work / "noid.lammpstrj"
    noid.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS type x y z\n"
        "1 0 0 0\n"
        "trailing\n"
    )
    expect(["read", str(noid), *one], 0, out_has="nop ")
    trunc = work / "trunc.lammpstrj"
    trunc.write_text("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
    expect(["read", str(trunc), *one], 0, out_has="nop ")
    nomol = work / "nomol.lammpstrj"
    nomol.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS type x y z id\n"
        "1 0 0 0 7\n"
        "1 0 0\n"
    )
    expect(["read", str(nomol), *one], 0, out_has="nop ")
    shortmol = work / "shortmol.lammpstrj"
    shortmol.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type mol\n"
        "1 1\n"
    )
    expect(["read", str(shortmol), *one], 2,
           err_line="LAMMPS dump frame is truncated: expected 1 atoms, read 0.")
    pdb = work / "end.pdb"
    pdb.write_text("END\n")
    expect(["read", str(pdb), "--frame", "1", "--last", "2"], 0, out_has="nop ")
    con = work / "end.con"
    con.write_text("END\n")
    expect(["read", str(con), "--frame", "1", "--last", "2"], 2)


if __name__ == "__main__":
    main()
