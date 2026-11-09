# Copyright (c) d-SEAMS developers
"""Run seams on real trajectories and check a printed count or token.

usage: cli_commands.py SEAMS ROOT
"""

import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]
root = Path(sys.argv[2])
water = root / "input" / "traj" / "exampleTraj.lammpstrj"
ice = root / "input" / "traj" / "genice_sI.lammpstrj"
tiny = root / "input" / "traj" / "agi_water_tiny.lammpstrj"
mw = root / "input" / "traj" / "mW_cubic.lammpstrj"
nucleation = root / "input" / "traj" / "nucleation_f1.lammpstrj"
xyz = root / "input" / "xyz" / "configuration.xyz"
one = ["--frame", "1", "--last", "1", "--cutoff", "3.5", "--jobs", "1"]


def run(args):
    return subprocess.run([seams, *args], capture_output=True, text=True)


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(name + " " + detail)


def expect(args, code, out_has=None, err_line=None):
    got = run(args)
    check("exit " + " ".join(args[-3:]), got.returncode == code,
          "rc " + str(got.returncode) + " " + got.stderr[:240])
    if out_has is not None:
        check("stdout " + out_has, out_has in got.stdout, got.stdout[:240])
    if err_line is not None:
        first = got.stderr.splitlines()[0] if got.stderr else ""
        check("stderr " + err_line, first == err_line, first)
    return got


def main():
    expect(["--version"], 0, out_has="seams")
    expect(["--features"], 0, out_has="OpenMP")
    expect(["--print-config"], 0, out_has="SEAMS_CARRY=")
    expect(["read", str(tiny), *one], 0, out_has="nop ")
    expect(["rings", str(tiny), *one], 0, out_has="rings ")
    expect(["chill", str(water), *one, "--type", "1"], 0, out_has="nop ")
    expect(["chill-plus", str(water), *one, "--type", "1"], 0, out_has="nop ")
    expect(["chill_plus", str(water), *one, "--layers", "--axis", "z"], 0, out_has="nop ")
    expect(["chill-plus", str(mw), *one, "--water-types", "1"], 0, out_has="cubic ")
    expect(["chill-plus", str(nucleation), *one], 0, out_has="hexagonal ")
    expect(["chill-plus", str(ice), *one, "--graph", "mutual"], 0, out_has="nop ")
    expect(["cages", str(ice), *one, "--graph", "cutoff"], 0, out_has="nop ")
    expect(["cages", str(ice), *one, "--graph", "mutual", "--complete"], 0, out_has="nop ")
    expect(["cages", str(tiny), *one, "--graph", "knn", "-k", "3"], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), *one, "--hops", "2", "--type", "1"], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), *one, "--by-atom", "--colour-types"], 0, out_has="nop ")
    expect(["f4", str(tiny), *one, "--type", "1", "--htype", "1"], 0, out_has="f4 ")
    expect(["steinhardt", str(tiny), *one, "--type", "1"], 0, out_has="nop ")
    expect(["rdf", str(tiny), *one, "--bins", "4", "--type", "1"], 0, out_has="# r g")
    expect(["rdf", str(water), *one, "--bins", "8", "--types", "1,2"], 0, out_has="# r g")
    expect(["cn", str(tiny), *one, "--bins", "4", "--type", "1"], 0, out_has="cn ")
    expect(["density-z", str(tiny), *one, "--bins", "0", "--axis", "Z"], 0, out_has="# z rho")
    expect(["density-z", str(tiny), *one, "--axis", "X", "--format", "json"], 0, out_has='"command":"density-z"')
    expect(["density-z", str(tiny), *one, "--axis", "Y"], 0, out_has="# z rho")
    expect(["density-z", str(tiny), *one, "--axis", "1"], 0, out_has="# z rho")
    expect(["hbonds", str(tiny), *one, "--type", "1", "--htype", "1", "--donors"], 0, out_has="hbonds ")
    expect(["rings", str(xyz), "--frame", "1", "--last", "1", "--jobs", "1"], 0, out_has="rings ")

    _cleanup = tempfile.TemporaryDirectory()
    two = Path(_cleanup.name) / "two.lammpstrj"
    body = tiny.read_text()
    two.write_text(body + "\n" + body.replace("TIMESTEP\n0\n", "TIMESTEP\n1\n", 1))
    expect(["cages", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--cutoff", "3.5", "--graph", "cutoff"], 0, out_has="nop ")
    expect(["chill-plus", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--cutoff", "3.5", "--rebuild"], 0, out_has="nop ")
    expect(["fingerprint", str(two), "--frame", "1", "--last", "2", "--jobs", "1",
            "--cutoff", "3.5", "--hops", "1"], 0, out_has="nop ")
    expect(["rings", str(two), "--frame", "1", "--last", "2", "--jobs", "1"], 0, out_has="rings ")

    per = two.parent / "per.txt"
    expect(["cages", str(tiny), *one, "--graph", "cutoff", "--per-atom", str(per)], 0)
    check("per-atom file", per.is_file() and per.stat().st_size > 0)

    emitted = expect(["fingerprint", str(tiny), *one, "--emit-library", "tiny"], 0)
    lib = two.parent / "lib.txt"
    lib.write_text(emitted.stdout)
    check("library text", "tiny" in lib.read_text() or lib.stat().st_size > 0)
    expect(["fingerprint", str(tiny), *one, "--library", str(lib)], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), *one, "--graph", "knn", "-k", "4"], 0, out_has="nop ")
    expect(["fingerprint", str(tiny), *one, "--graph", "seeded"], 0, out_has="nop ")
    expect(["cages", str(tiny), *one, "--graph", "knn", "--water-types", "1"], 0, out_has="nop ")
    expect(["cages", str(tiny), *one, "--graph", "cutoff", "--signature", "512",
            "--guest-types", "3", "--inside", "--incomplete"], 0, out_has="nop ")
    expect(["cages", str(mw), *one, "--graph", "cutoff", "--shape", "--cage-edge", "0,1",
            "--water-types", "1"], 0, out_has="nop ")
    expect(["ions", str(tiny), *one, "--type", "1", "--ion-types", "3,4",
            "--ion-cutoff", "5", "--complete"], 0, out_has="ions ")
    expect(["f4", str(tiny), *one, "--type", "0", "--htype", "1"], 0, out_has="f4 ")
    expect(["steinhardt", str(mw), *one, "--type", "1", "--format", "json"], 0,
           out_has='"command":"steinhardt"')
    expect(["chill", str(mw), *one, "--graph", "knn", "--water-types", "1"], 0, out_has="cubic ")
    expect(["chill-plus", str(nucleation), *one, "--graph", "cutoff", "--layers",
            "--axis", "y", "--tum-layers"], 0, out_has="nop ")
    expect(["rings", str(tiny), *one, "--tpp", "2", "--block", "64",
            "--resident", "0.4", "--cell", "2.5"], 0, out_has="rings ")
    expect(["read", str(tiny), "--frame", "1", "--last", "1", "--type", "0"], 0, out_has="nop ")
    only = two.parent / "type1.lammpstrj"
    only.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 1 1 1 1\n2 1 1.2 1 1\n"
    )
    expect(["f4", str(only), "--frame", "1", "--last", "1", "--type", "0",
            "--htype", "1", "--cutoff", "3.5"], 0, out_has="f4 ")
    ions = two.parent / "ions.lammpstrj"
    ions.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n4\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id type x y z\n"
        "1 1 1 1 1\n2 2 1.2 1 1\n3 2 8 8 8\n4 3 9 9 9\n"
    )
    expect(["cn", str(ions), "--ions", "--site", "1=cationHead,2=anion,3=polar",
            "--cutoff", "6", "--frame", "1", "--last", "1"], 0, out_has="nCation ")
    expect(["cages", str(nucleation), "--frame", "1", "--last", "1", "--jobs", "1",
            "--cutoff", "3.5", "--graph", "cutoff", "--water-types", "1"], 0, out_has="both ")
    expect(["hbonds", str(water), *one, "--type", "1", "--htype", "2"], 0, out_has="hbonds ")
    expect(["rings", str(two), "--frame", "1", "--last", "2", "--jobs", "2"], 0, out_has="rings ")
    expect(["fingerprint", str(tiny), "--frame", "1", "--last", "1", "--cutoff", "3.5"],
           0, out_has="nop ")
    expect(["cages", str(tiny), "--frame", "1", "--last", "1", "--cutoff", "3.5"],
           0, out_has="nop ")
    expect(["rdf", str(tiny), "--frame", "1", "--last", "1", "--type", "0", "--bins", "4"],
           0, out_has="# r g")
    expect(["cn", str(tiny), "--frame", "1", "--last", "1", "--type", "0", "--bins", "4"],
           0, out_has="cn ")
    expect(["domains", str(water), "--site", "1=polar", "--cutoff", "6"], 0, out_has="subset ")
    expect(["rings", str(two), "--frame", "1", "--last", "2", "--strict-input"], 0, out_has="rings ")
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["SEAMS_CARRY"] = "0"
    got = subprocess.run([seams, "cages", str(two), "--frame", "1", "--last", "2",
                          "--jobs", "1", "--cutoff", "3.5"],
                         capture_output=True, text=True, env=env)
    check("carry off", got.returncode == 0 and "nop " in got.stdout, got.stderr[:160])
    expect(["fingerprint", str(tiny), *one, "--format", "json",
            "--emit-library", "a\x01b"], 0, out_has="\\u0001")
    expect(["rings", str(water), "--water-types", "1,,2", "--frame", "1", "--last", "1"],
           0, out_has="rings ")
    expect(["cages", str(tiny), *one, "--guest-types", "1,,3", "--signature", "512"],
           0, out_has="nop ")
    expect(["ions", str(ions), "--ion-types", "3,,4", "--frame", "1", "--last", "1"],
           0, out_has="ions ")


if __name__ == "__main__":
    main()
