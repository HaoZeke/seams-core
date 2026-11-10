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


ok = run(["rings", str(water), "--select", "type 2", "--frame", "1", "--last", "2"])
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
        ["f4", str(tiny), "--select", "type 99", "--frame", "1", "--last", "1"],
        ["ions", str(tiny), "--select", "type 99", "--ion-types", "2", "--frame", "1", "--last", "1"],
        ["cages", str(tiny), "--select", "type 99", "--frame", "1", "--last", "1"],
        ["fingerprint", str(tiny), "--select", "name ZZ", "--frame", "1", "--last", "1"],
        ["rings", str(empty), "--frame", "2", "--last", "2"],
        ["hbonds", str(water), *one, "--type", "1", "--htype", "2", "--donors", "--hdist", "3.5"],
        ["rings", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "2"],
        ["read", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "1"],
        ["rings", str(root / "input" / "con" / "tiny_multi_cuh2.con"), "--frame", "1", "--last", "9"],
    ]
    junk = Path(tmp) / "junk.lammpstrj"
    junk.write_text("this is not a dump\n")
    calls.append(["rings", str(junk), "--frame", "1", "--last", "1"])
    calls.append(["rings", str(junk), "--select", "type 1", "--frame", "1", "--last", "1"])
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
    blank_rc = subprocess.run([seams, "rings", str(blank), "--frame", "1", "--last", "1"], capture_output=True).returncode
    ok = blank_rc in (0, 2) and ok
    ok = run(["chill-plus", str(water), "--family", "not-a-family"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "1,2,3"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "a,b"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", "1,"], expect=2, message="Usage:") and ok
    ok = run(["rdf", str(water), "--types", ",1"], expect=2, message="Usage:") and ok
    ok = run(["fingerprint", str(tiny), *one, "--library", str(Path(tmp) / "missing-lib.txt")], expect=2, message="cannot read key library") and ok
    ver = subprocess.run([seams, "--version"], capture_output=True, text=True)
    ok = ver.returncode == 0 and "seams" in ver.stdout and ok
    ok = resolve_dump(root / "input" / "traj" / "not-deposited.lammpstrj") == tiny and ok
    def demand_message(fn, want):
        try:
            fn()
        except SystemExit as exc:
            if str(exc) != want:
                raise
            return
        raise SystemExit("cutoff check did not fire")

    def captured(fn, want):
        try:
            demand_message(fn, want)
        except SystemExit as exc:
            return str(exc)
        return ""

    def same(left, right):
        return left == right

    want_code = "rings %s --cutoff not-a-number exit 2 want 0" % water
    demand_message(lambda: run(["rings", str(water), "--cutoff", "not-a-number"], expect=0), want_code)
    demand_message(
        lambda: run(["rings", str(water), "--cutoff", "not-a-number"], expect=2, message="not-the-message"),
        "missing reject message",
    )
    got_match = captured(lambda: run(["rings", str(water), "--cutoff", "not-a-number"], expect=0), want_code)
    got_bad = captured(lambda: run(["rings", str(water), "--cutoff", "not-a-number"], expect=0), "not the real message")
    got_quiet = captured(lambda: run(["--version"], expect=0), "x")
    seen = [same(got_match, ""), same(got_bad, ""), same(got_quiet, "cutoff check did not fire"), "exit 2 want 0" in got_bad]
    if seen != [True, False, True, True]:
        raise SystemExit("demand arms %s" % seen)
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
    ]
    for args, message in rejected:
        ok = run(args, expect=2, message=message) and ok

sys.exit(0 if ok else 1)
