"""Run every seams command on a one-frame trajectory.

The process must exit 0. That is the coverage of the command bodies.
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
mw = root / "input" / "traj" / "tiny-selectionTest.lammpstrj"
if not mw.exists():
    mw = tiny

one = ["--frame", "1", "--last", "1", "--cutoff", "3.5"]


def run(args, expect=0):
    proc = subprocess.run([seams, *args], capture_output=True, text=True)
    if proc.returncode != expect:
        sys.stderr.write("fail " + " ".join(args) + " rc " + str(proc.returncode) + "\n")
        sys.stderr.write(proc.stdout[-400:] + proc.stderr[-400:] + "\n")
        return False
    return True


ok = True
with tempfile.TemporaryDirectory() as tmp:
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
    ]
    for args in calls:
        ok = run(args) and ok
    emitted = subprocess.run(
        [seams, "fingerprint", str(tiny), *one, "--emit-library", "tiny"],
        capture_output=True, text=True)
    libp = Path(tmp) / "keys.txt"
    libp.write_text(emitted.stdout)
    ok = emitted.returncode == 0 and libp.stat().st_size > 0 and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(libp), "--per-atom", per]) and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(libp) + "," + str(libp)], expect=2) and ok
    missing = str(Path(tmp) / "no-such.lammpstrj")
    ok = run(["rings", missing]) and ok
    rejected = [
        [],
        ["no-such-command", str(water)],
        ["rings", str(water), "--format", "yaml"],
        ["density-z", str(water), "--axis", "q"],
        ["rdf", str(water), "--bins", "-1"],
        ["rdf", str(water), "--types", "1"],
        ["ions", str(water)],
        ["pairs", str(water)],
        ["cn", str(water), "--ions"],
        ["domains", str(water), "--site", "1=cationHead", "--subset", "nope"],
        ["cages", str(water), "--cage-edge", "1.0"],
        ["pairs", str(water), "--site", "1=not-a-kind"],
        ["pairs", str(water), "--site", "no-equals"],
        ["pairs", str(water), "--site", "abc=cationHead"],
        ["chill", str(water), "--graph", "foo"],
        ["rings", str(water), "--frame", "1-2"],
        ["cages", str(tiny), "--guest-types", "x"],
    ]
    for args in rejected:
        ok = run(args, expect=2) and ok

sys.exit(0 if ok else 1)
