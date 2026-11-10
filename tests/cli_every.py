# Copyright (c) d-SEAMS developers
"""Run every seams command on a one-frame trajectory.

The process must exit 0. That is the coverage of the command bodies.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

seams = sys.argv[1]
root = Path(sys.argv[2])
water = root / "input" / "traj" / "exampleTraj.lammpstrj"
ice = root / "input" / "traj" / "genice_sI.lammpstrj"
tiny = root / "input" / "traj" / "agi_water_tiny.lammpstrj"
def resolve_dump(path):
    if not path.is_file():
        return tiny
    return path


mw = resolve_dump(root / "input" / "traj" / "tiny-selectionTest.lammpstrj")

one = ["--frame", "1", "--last", "1", "--cutoff", "3.5"]


def run(args, expect=0, message=None):
    proc = subprocess.run([seams, *args], capture_output=True, text=True)
    text = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != expect:
        label = " ".join(args[:4]) if args else "no-command"
        raise SystemExit("%s exit %d want %d" % (label, proc.returncode, expect))
    if expect != 0 and (not message or message not in text):
        raise SystemExit("missing reject message")
    return True


ok = run(["rings", str(water), "--frame", "1", "--last", "2"])
with tempfile.TemporaryDirectory() as tmp:
    ions = Path(tmp) / "ions.lammpstrj"
    ions.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n4\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
        "1 1 1.0 1.0 1.0\n2 1 1.2 1.0 1.0\n3 2 1.1 1.1 1.0\n4 2 8.0 8.0 8.0\n"
    )
    blank = Path(tmp) / "blank.con"
    blank.write_text("0\n")
    zero = Path(tmp) / "zero.lammpstrj"
    zero.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 0\n0 0\n0 0\n"
        "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
    )
    empty = Path(tmp) / "empty.lammpstrj"
    empty.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 1 1 1 1\n"
        "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n0\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n"
    )
    per = str(Path(tmp) / "per.lammpstrj")
    lib = str(Path(tmp) / "lib.txt")
    calls = [
        ["--help"],
        ["--version"],
        ["--features"],
        ["--print-config"],
        ["read", str(water), *one],
        ["rings", str(water), *one],
        ["rings", str(ice), *one, "--jobs", "1"],
        ["chill", str(water), *one, "--type", "1"],
        ["chill-plus", str(water), *one, "--type", "1"],
        ["chill_plus", str(water), *one, "--layers", "--axis", "z"],
        ["chill-plus", str(water), *one, "--tum-layers", "--graph", "cutoff"],
        ["chill-plus", str(ice), *one, "--graph", "mutual", "--family", "waterIce"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--type", "1"],
        ["cages", str(ice), *one, "--graph", "mutual", "--complete"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--shape", "--cage-edge", "0,1"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--signature", "512",
         "--guest-types", "1", "--guest-radius", "4"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--per-atom", per, "--inside"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--incomplete", "--signature", "512"],
        ["cages", str(water), *one, "--water-types", "1", "--format", "json"],
        ["fingerprint", str(water), *one, "--type", "1", "--hops", "2"],
        ["fingerprint", str(water), *one, "--by-atom", "--colour-types"],
        ["fingerprint", str(tiny), *one, "--emit-library", "tiny"],
        ["ions", str(water), *one, "--type", "1", "--ion-types", "2", "--ion-cutoff", "4"],
        ["f4", str(water), *one, "--type", "1", "--htype", "2"],
        ["steinhardt", str(water), *one, "--type", "1"],
        ["rdf", str(water), *one, "--bins", "20", "--type", "1"],
        ["rdf", str(water), *one, "--bins", "10", "--types", "1,2"],
        ["cn", str(water), *one, "--bins", "20", "--type", "1"],
        ["cn", str(water), "--ions", "--site", "1=cationHead,2=anion", "--cutoff", "6"],
        ["pairs", str(water), "--site", "1=cationHead,2=anion", *one],
        ["pairs", str(water), "--site",
         "1=tail,2=donorH,1=acceptor,2=waterO,1=waterH,2=solvent,1=unspecified,2=polar",
         *one],
        ["density-z", str(water), *one, "--type", "1", "--bins", "10", "--axis", "z"],
        ["density-z", str(water), "--bins", "8", "--axis", "x"],
        ["domains", str(water), "--site", "1=cationHead,2=anion", "--cutoff", "6",
         "--subset", "polar"],
        ["domains", str(water), "--site", "1=cationHead,2=anion", "--cutoff", "6",
         "--subset", "apolar"],
        ["hbonds", str(water), *one, "--type", "1", "--htype", "2",
         "--hdist", "2.5", "--hangle", "30"],
        ["hbonds", str(water), *one, "--type", "1", "--htype", "2", "--donors"],
        ["rings", str(tiny), *one, "--rebuild", "--strict-input"],
        ["cages", str(tiny), *one, "--graph", "knn", "--complete"],
        ["density-z", str(water), *one, "--axis", "y", "--format", "json"],
        ["density-z", str(water), *one, "--axis", "x"],
        ["density-z", str(water), *one, "--axis", "0"],
        ["chill-plus", str(water), *one, "--graph", "cutoff", "--layers", "--axis", "y"],
        ["f4", str(tiny), *one, "--type", "1", "--htype", "1"],
        ["steinhardt", str(tiny), *one, "--type", "1", "--format", "json"],
        ["rdf", str(tiny), *one, "--bins", "4", "--types", "1,1"],
        ["cn", str(tiny), *one, "--bins", "4", "--type", "1", "--format", "json"],
        ["hbonds", str(tiny), *one, "--donors", "--format", "json"],
        ["ions", str(water), *one, "--ion-types", "1", "--ion-cutoff", "3", "--complete"],
        ["domains", str(water), "--site", "1=tail,2=apolar", "--cutoff", "6", "--subset", "tail"],
        ["pairs", str(water), "--site", "1=apolar,2=polar", *one],
        ["fingerprint", str(tiny), *one, "--hops", "1", "--format", "json"],
        ["fingerprint", str(tiny), *one, "--emit-library", 'a\\b"c', "--format", "json"],
        ["fingerprint", str(tiny), *one, "--emit-library", "a\rb\x01", "--format", "json"],
        ["fingerprint", str(tiny), "--type", "9", "--frame", "1", "--last", "1"],
        ["cages", str(tiny), "--type", "9", "--signature", "512", "--frame", "1", "--last", "1"],
        ["f4", str(tiny), "--type", "9", "--htype", "2", "--frame", "1", "--last", "1"],
        ["steinhardt", str(tiny), "--type", "9", "--frame", "1", "--last", "1"],
        ["hbonds", str(tiny), "--type", "9", "--htype", "9", "--frame", "1", "--last", "1"],
        ["cages", str(ice), *one, "--graph", "cutoff", "--signature", "512+hc",
         "--shape", "--guest-types", "1", "--inside", "--incomplete"],
        ["cages", str(ice), *one, "--graph", "seeded", "--signature", "512", "--water-types", "1"],
        ["cages", str(ice), *one, "--graph", "knn", "--water-types", "1"],
        ["chill", str(ice), *one, "--graph", "knn", "--water-types", "1"],
        ["pairs", str(water), "--site", "family=ionicLiquid,1=cationHead,2=anion", *one],
        ["pairs", str(water), "--site", "1=cationHead,,2=anion", *one],
        ["rings", str(root / "input" / "xyz" / "configuration.xyz"), "--frame", "1", "--last", "1"],
        ["cn", str(ions), "--ions", "--site", "1=cationHead,2=anion", "--cutoff", "6"],
        ["pairs", str(ions), "--site", "1=cationHead,2=anion", "--cutoff", "6"],
        ["ions", str(ions), "--type", "1", "--ion-types", "2", "--ion-cutoff", "4", "--per-atom", per, "--complete"],
        ["ions", str(tiny), "--type", "9", "--ion-types", "2", "--frame", "1", "--last", "1"],
        ["rdf", str(water), *one, "--bins", "0", "--type", "1"],
        ["cn", str(water), *one, "--bins", "0", "--types", "1,2"],
        ["f4", str(tiny), *one, "--htype", "1"],
        ["chill-plus", str(root / "input" / "traj" / "iceXXI_lee2026.lammpstrj"), *one],
        ["chill", str(root / "input" / "traj" / "iceXXI_lee2026.lammpstrj"), *one],
        ["chill-plus", str(root / "input" / "traj" / "nucleation_f1.lammpstrj"), *one],
        ["chill-plus", str(root / "input" / "traj" / "mW_cubic.lammpstrj"), *one, "--water-types", "1"],
        ["cages", str(tiny), *one, "--graph", "knn", "-k", "3"],
        ["chill-plus", str(tiny), "--type", "9", "--frame", "1", "--last", "1"],
        ["chill", str(tiny), "--type", "9", "--frame", "1", "--last", "1"],
        ["f4", str(tiny), "--type", "0", "--htype", "1", "--frame", "1", "--last", "1"],
        ["chill-plus", str(tiny), *one, "--water-types", "2"],
        ["ions", str(ions), "--type", "2", "--ion-types", "2", "--water-types", "1", "--ion-cutoff", "4"],
        ["cn", str(water), "--ions", "--site", "9=cationHead", "--cutoff", "6"],
        ["cn", str(ions), "--ions", "--site", "2=cationHead", "--cutoff", "6"],
        ["cn", str(zero), "--bins", "0", "--type", "1", "--frame", "1", "--last", "1"],
        ["fingerprint", str(tiny), *one, "--per-atom", per],
        ["fingerprint", str(tiny), *one, "--emit-library", "a\tb", "--format", "json"],
        ["rings", str(tiny), *one, "--tpp", "2", "--block", "32", "--resident", "0.5",
         "--cell", "3.2"],
        ["f4", str(tiny), "--frame", "1", "--last", "1"],
        ["ions", str(tiny), "--ion-types", "2", "--frame", "1", "--last", "1"],
        ["cages", str(tiny), "--frame", "1", "--last", "1"],
        ["fingerprint", str(tiny), "--frame", "1", "--last", "1"],
        ["rings", str(empty), "--frame", "2", "--last", "2"],
        ["hbonds", str(water), *one, "--type", "1", "--htype", "2", "--donors", "--hdist", "3.5"],
        ["rings", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "2"],
        ["read", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "1"],
        ["rings", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "9"],
    ]
    junk = Path(tmp) / "junk.lammpstrj"
    junk.write_text("this is not a dump\n")
    calls.append(["rings", str(junk), "--frame", "1", "--last", "1"])
    for args in calls:
        ok = run(args) and ok
    emitted = subprocess.run(
        [seams, "fingerprint", str(tiny), *one, "--emit-library", "tiny"],
        capture_output=True, text=True)
    libp = Path(tmp) / "keys.txt"
    libp.write_text(emitted.stdout)
    run_check = emitted.returncode == 0 and libp.stat().st_size > 0
    ok = run_check and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(libp), "--per-atom", per]) and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(libp) + "," + str(libp)], expect=2, message="two key libraries share a hop count") and ok
    hop1 = subprocess.run(
        [seams, "fingerprint", str(tiny), *one, "--hops", "1", "--emit-library", "one"],
        capture_output=True, text=True)
    hop2 = subprocess.run(
        [seams, "fingerprint", str(tiny), *one, "--hops", "2", "--emit-library", "two"],
        capture_output=True, text=True)
    p1 = Path(tmp) / "hop1.txt"
    p2 = Path(tmp) / "hop2.txt"
    p1.write_text(hop1.stdout)
    p2.write_text(hop2.stdout)
    ok = hop1.returncode == 0 and hop2.returncode == 0 and ok
    ok = run(["fingerprint", str(tiny), *one, "--hops", "2", "--library", str(p1) + "," + str(p2)]) and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", "," + str(libp), "--per-atom", str(Path(tmp))], expect=2, message="cannot write per-atom dump") and ok
    ok = run(["pairs", str(water), "--site", "=cationHead", *one], expect=2, message="site spec token '=cationHead' is not key=value") and ok
    ok = run(["chill-plus", str(water), "--graph", "knn-union"], expect=2, message="CHILL and CHILL+ take --graph cutoff or mutual") and ok
    ok = run(["rings", str(empty), "--frame", "1", "--last", "2", "--strict-input"], expect=2, message="input frame 2 contains no atoms") and ok
    def text_rejected(code, text, needle):
        if code != 0 or needle not in text:
            raise SystemExit(text)

    blank_proc = subprocess.run(
        [seams, "rings", str(blank), "--frame", "1", "--last", "1"],
        capture_output=True, text=True,
    )
    blank_text = blank_proc.stdout + blank_proc.stderr
    text_rejected(blank_proc.returncode, blank_text, "Frame 1 not found")
    seen = ""
    try:
        text_rejected(0, "rings 3:0", "Frame 1 not found")
    except SystemExit as err:
        seen = str(err)
    assert seen == "rings 3:0"
    ok = run(["chill-plus", str(water), "--family", "not-a-family"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "1,2,3"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "a,b"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "1,"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", ",1"], expect=2, message="Usage:") and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(Path(tmp) / "missing-lib.txt")], expect=2, message="cannot read key library") and ok
    ver = subprocess.run([seams, "--version"], capture_output=True, text=True)
    ok = ver.returncode == 0 and "seams" in ver.stdout and ok
    ok = resolve_dump(root / "input" / "traj" / "not-deposited.lammpstrj") == tiny and ok
    ok = run(["rings", str(water), "--cutoff", "not-a-number"], expect=2, message="Usage:") and ok
    def with_offload():
        off = os.environ.get("SEAMS_OFFLOAD")
        os.environ["SEAMS_OFFLOAD"] = "1"
        try:
            return run(["cages", str(tiny), *one, "--graph", "cutoff"])
        finally:
            if off is None:
                os.environ.pop("SEAMS_OFFLOAD", None)
            else:
                os.environ["SEAMS_OFFLOAD"] = off

    os.environ.pop("SEAMS_OFFLOAD", None)
    ok = with_offload() and ok
    os.environ["SEAMS_OFFLOAD"] = "0"
    ok = with_offload() and ok
    missing = str(Path(tmp) / "no-such.lammpstrj")
    ok = run(["rings", missing]) and ok
    rejected = [
        ([], "Usage:"),
        (["no-such-command", str(water)], "unknown command: no-such-command"),
        (["rings", str(water), "--format", "yaml"], "bad --format (want text|json)"),
        (["density-z", str(water), "--axis", "q"], "Usage:"),
        (["rdf", str(water), "--bins", "-1"], "bad --bins (want N > 0)"),
        (["rdf", str(water), "--types", "1"], "Usage:"),
        (["ions", str(water)], "ions needs --ion-types"),
        (["pairs", str(water)], "pairs needs --site"),
        (["cn", str(water), "--ions"], "cn --ions needs --site"),
        (["domains", str(water), "--site", "1=cationHead", "--subset", "nope"], "bad --subset (want polar|apolar)"),
        (["cages", str(water), "--cage-edge", "1.0"], "Usage:"),
        (["pairs", str(water), "--site", "1=not-a-kind"], "unknown site kind 'not-a-kind'"),
        (["pairs", str(water), "--site", "no-equals"], "site spec token 'no-equals' is not key=value"),
        (["pairs", str(water), "--site", "abc=cationHead"], "site spec type 'abc' is not an integer"),
        (["chill", str(water), "--graph", "foo"], "unknown bond graph 'foo' (use cutoff, knn, knn-union)"),
        (["rings", str(water), "--frame", "1-2"], "Usage:"),
        (["cages", str(tiny), "--guest-types", "x"], "invalid arguments: value \"x\" is not a number"),
        (["chill", str(tiny), "--family", "ionicLiquid"], "CHILL/TUM refused for family ionicLiquid"),
        (["chill", str(tiny), "--family", "moltenSalt"], "CHILL/TUM refused for family moltenSalt"),
        (["chill", str(tiny), "--family", "des"], "CHILL/TUM refused for family des"),
        (["chill", str(tiny), "--family", "electrolyte"], "CHILL/TUM refused for family electrolyte (name a waterIce subset)"),
        (["chill", str(tiny), "--family", "confinedIL"], "CHILL/TUM refused for family confinedIL"),
        (["chill", str(tiny), "--family", "confinedWater"], "CHILL/TUM refused for family confinedWater"),
        (["chill", str(tiny), "--family", "networkFormer"], "CHILL/TUM refused for family networkFormer"),
        (["domains", str(tiny)], "domains needs --site"),
        (["rings", str(tiny), "--format", "yaml"], "bad --format (want text|json)"),
        (["rings", str(tiny), "--bins", "-3"], "bad --bins (want N > 0)"),
        (["cages", str(tiny), "--cage-edge", "9,1"], "bad cage edge window (want 0 <= lo < hi, or none)"),
        ([], "A command and a trajectory file are required"),
    ]
    for args, message in rejected:
        ok = run(args, expect=2, message=message) and ok
    xyz = Path(tmp) / "two.xyz"
    xyz.write_text("2\nbox\nO 0 0 0\nO 2 0 0\n")
    nucleation = root / "input" / "traj" / "nucleation_f1.lammpstrj"
    observed = [
        (["cages", str(nucleation), "--type", "1", "--graph", "cutoff", "--shape",
          "--per-atom", per, "--frame", "1", "--last", "1"], "both"),
        (["cages", str(nucleation), "--type", "1", "--graph", "seeded",
          "--frame", "1", "--last", "1"], "seeded"),
        (["chill-plus", str(nucleation), "--type", "1", "--frame", "1", "--last", "1"], "hexagonal"),
        (["cages", str(tiny), "--graph", "knn-union", "--frame", "1", "--last", "1"], "knn-union"),
        (["fingerprint", str(tiny), "--graph", "knn", "--by-atom", "--colour-types",
          "--frame", "1", "--last", "1"], "atom"),
        (["fingerprint", str(tiny), "--graph", "seeded", "--frame", "1", "--last", "1"], "graph"),
        (["rings", str(xyz), "--frame", "1", "--last", "1"], "rings"),
        (["read", str(zero), "--frame", "1", "--last", "1"], "nop 1 frame 1 box"),
        (["rings", str(zero), "--strict-input", "--frame", "1", "--last", "1"],
         "The point cloud has no simulation box; cannot build a neighbour list."),
        (["cages", str(tiny), "--cage-edge", "0,1", "--graph", "cutoff",
          "--frame", "1", "--last", "1"], "graph"),
        (["rings", str(water), "--jobs", "2", "--frame", "1", "--last", "2"], "rings"),
        (["rings", str(water), "--rebuild", "--jobs", "2", "--frame", "1", "--last", "2"], "rings"),
    ]
    for args, needle in observed:
        proc = subprocess.run([seams, *args], capture_output=True, text=True)
        text_rejected(proc.returncode, proc.stdout + proc.stderr, needle)

    missed = ""
    try:
        run(["rings", str(tiny), "--frame", "1", "--last", "1"], expect=2)
    except SystemExit as err:
        missed = str(err)
    assert "want 2" in missed
    missed = ""
    try:
        run(["domains", str(tiny)], expect=2, message="not-the-domains-text")
    except SystemExit as err:
        missed = str(err)
    assert missed == "missing reject message"

    anion = Path(tmp) / "anion.lammpstrj"
    anion.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
        "ITEM: ATOMS id type x y z\n1 2 1 1 1\n2 2 2 1 1\n"
    )
    far = Path(tmp) / "far.lammpstrj"
    far.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 30\n0 30\n0 30\n"
        "ITEM: ATOMS id type x y z\n1 1 1 1 1\n2 1 20 1 1\n"
    )
    hydro = Path(tmp) / "hydro.lammpstrj"
    hydro.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 2 1 1 1\n2 2 2 1 1\n"
    )
    pdb = Path(tmp) / "one.pdb"
    pdb.write_text(
        "ATOM      1  O   MOL A   1       1.000   1.000   1.000  1.00  0.00           O\n"
        "END\n"
    )
    nodot = Path(tmp) / "nodotdump"
    nodot.write_bytes(tiny.read_bytes())
    checked = [
        (["cages", str(tiny), "--signature", "nope", "--frame", "1", "--last", "1"],
         2, "bad cage signature 'nope'"),
        (["rings", "--frame", "1"], 2, "A command and a trajectory file are required"),
        (["density-z", str(tiny), "--axis", "Y", "--bins", "2", "--frame", "1", "--last", "1"],
         0, "5 0.00025"),
        (["density-z", str(tiny), "--axis", "1", "--bins", "2", "--frame", "1", "--last", "1"],
         0, "15 0.0015"),
        (["density-z", str(tiny), "--axis", "2", "--bins", "2", "--frame", "1", "--last", "1"],
         0, "15 0.00175"),
        (["density-z", str(tiny), "--axis", "Z", "--bins", "2", "--frame", "1", "--last", "1"],
         0, "5 0"),
        (["cn", str(anion), "--ions", "--site", "2=anion", "--cutoff", "6", "--frame", "1", "--last", "1"],
         0, "degree 0 nCation 0 nAnion 2"),
        (["cn", str(far), "--ions", "--site", "1=cationHead", "--cutoff", "3", "--frame", "1", "--last", "1"],
         0, "degree 0 nCation 2 nAnion 0"),
        (["f4", str(hydro), "--type", "0", "--htype", "2", "--frame", "1", "--last", "1"],
         0, "f4 nan n 0"),
        (["read", str(pdb), "--frame", "1", "--last", "1"], 0, "nop 1 frame 1 box"),
        (["rings", str(xyz), "--frame", "2", "--last", "2"],
         0, "Frame 2 does not exist in"),
        (["rings", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "0"],
         0, "rings 3:0 4:0 5:0 6:0 7:0 8:0"),
        (["rings", str(nodot), "--frame", "1", "--last", "1"], 0, "chemfiles cannot read"),
        (["domains", str(far), "--site", "1=cationHead", "--cutoff", "6", "--frame", "1", "--last", "1"],
         0, "subset polar n 2 largest 1 P 0.5"),
        (["rings", str(far), "--water-types", ",", "--frame", "1", "--last", "1"],
         0, "rings 3:0 4:0 5:0 6:0 7:0 8:0"),
        (["rings", str(xyz), "--jobs", "2", "--frame", "1", "--last", "1"],
         0, "The point cloud has no simulation box; cannot build a neighbour list."),
        (["cages", str(tiny), "--cage-edge", "none", "--signature", "512", "--frame", "1", "--last", "1"],
         0, "signature 5:12 cages 0 atoms 0"),
        (["cages", str(far), "-k", "0", "--frame", "1", "--last", "1"],
         0, "graph seeded hexagonal 0 cubic 0 both 0 water 0"),
        (["fingerprint", str(tiny), "--graph", "", "--frame", "1", "--last", "1"],
         0, "graph cutoff"),
        (["fingerprint", str(tiny), "--graph", "mutual", "--frame", "1", "--last", "1"],
         0, "graph mutual"),
        (["cages", str(root / "input" / "traj" / "genice_sI.lammpstrj"),
          "--graph", "cutoff", "--signature", "512", "--shape", "--type", "1",
          "--frame", "1", "--last", "1"],
         0, "shape 20:off 20:off"),
    ]
    def demand(proc, code, needle):
        text = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode != code or needle not in text:
            raise SystemExit("%s exit %d\n%s" % (" ".join(proc.args[:3]), proc.returncode, text[-500:]))

    caught = ""
    try:
        demand(subprocess.CompletedProcess(["seams", "read", "x"], 2, "not the phrase", ""), 0, "phrase")
    except SystemExit as err:
        caught = str(err)
    if "not the phrase" not in caught:
        raise SystemExit("demand dropped the output")
    for args, code, needle in checked:
        demand(subprocess.run([seams, *args], capture_output=True, text=True), code, needle)
    library = Path(tmp) / "keys.txt"
    emitted = subprocess.run(
        [seams, "fingerprint", str(tiny), "--emit-library", "tiny", "--frame", "1", "--last", "1"],
        capture_output=True, text=True,
    )
    demand(emitted, 0, "tiny")
    library.write_text(emitted.stdout)
    matched = subprocess.run(
        [seams, "fingerprint", str(far), "--library", str(library), "--hops", "2",
         "--frame", "1", "--last", "1"],
        capture_output=True, text=True,
    )
    demand(matched, 0, "unmatched=2")

sys.exit(0 if ok else 1)
