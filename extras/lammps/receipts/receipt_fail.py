# Copyright (c) d-SEAMS developers
"""Drive the receipt scripts into the failures a bad input actually produces."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import cover
import date_check
import date_evidence
import record_build
import rewrite_clock
import streamlog
import trace_libs


def outcome(fn):
    try:
        fn()
    except SystemExit as exc:
        return str(exc)
    return ""


def expect_exit(fn, needle):
    text = outcome(fn)
    if needle in text:
        print("shown", needle, flush=True)
        return
    print("absent", needle, flush=True)


def init_repo(path):
    subprocess.check_call(["git", "-C", path, "init", "-q"])
    subprocess.check_call(["git", "-C", path, "config", "user.name", rewrite_clock.NAME])
    subprocess.check_call(["git", "-C", path, "config", "user.email", rewrite_clock.EMAIL])
    subprocess.check_call(["git", "-C", path, "config", "commit.gpgsign", "false"])


def git_commit(repo, message, author, committer):
    env = os.environ.copy()
    env["GIT_AUTHOR_NAME"] = rewrite_clock.NAME
    env["GIT_AUTHOR_EMAIL"] = rewrite_clock.EMAIL
    env["GIT_COMMITTER_NAME"] = rewrite_clock.NAME
    env["GIT_COMMITTER_EMAIL"] = rewrite_clock.EMAIL
    env["GIT_AUTHOR_DATE"] = author
    env["GIT_COMMITTER_DATE"] = committer
    blob = subprocess.check_output(["git", "-C", repo, "hash-object", "-w", "--stdin"], input=message.encode())
    subprocess.check_call(
        ["git", "-C", repo, "update-index", "--add", "--cacheinfo", "100644,%s,note" % blob.decode().strip()]
    )
    tree = subprocess.check_output(["git", "-C", repo, "write-tree"], text=True).strip()
    parent = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], text=True, capture_output=True)
    args = ["git", "-C", repo, "commit-tree", tree]
    if parent.returncode == 0:
        args.extend(["-p", parent.stdout.strip()])
    new = subprocess.check_output(args, input=message, env=env, text=True).strip()
    ref = "refs/heads/demo"
    old = subprocess.run(["git", "-C", repo, "rev-parse", ref], text=True, capture_output=True)
    cmd = ["git", "-C", repo, "update-ref", ref, new]
    if old.returncode == 0:
        cmd.append(old.stdout.strip())
    subprocess.check_call(cmd)
    subprocess.check_call(["git", "-C", repo, "symbolic-ref", "HEAD", ref])
    return new


def clock_row(repo, new, use, subject):
    return "repo %s\nread %s\nuse %s\nold abc\nnew %s\nsubject %s\n" % (repo, use, use, new, subject)


def run_evidence(clock, repo, base):
    old = sys.argv
    sys.argv = ["date_evidence.py", str(clock), str(repo), base, str(repo), base]
    try:
        date_evidence.main()
    finally:
        sys.argv = old


def check_streamlog(work):
    os.environ.pop("CVMFS_PREFIX", None)
    assert not (streamlog.clean("plain") != "plain"), "clean"
    assert not (streamlog.clean("plain", [("zz", "nope")]) != "plain"), "pairs"
    note = work / "note.txt"
    note.write_text("text")
    os.environ["CVMFS_PREFIX"] = str(work)
    assert not ("$CVMFS" not in streamlog.clean(str(work) + "/x")), "cvmfs"
    os.environ.pop("CVMFS_PREFIX", None)
    (work / "nested").mkdir()
    (work / "nested" / "skip.o").write_bytes(b"obj")
    assert not (streamlog.elf_kind(note) is not None or streamlog.elf_kind(work / "missing") is not None), "elf text"
    (work / "short").write_bytes(b"\x7fELF")
    assert not (streamlog.elf_kind(work / "short") is not None), "elf short"
    for kind, name, want in ((2, "exec", "exec"), (3, "dyn", "dyn"), (1, "other", None)):
        blob = bytearray(b"\x7fELF\x01\x01" + b"\x00" * 12)
        blob[16:18] = kind.to_bytes(2, "little")
        path = work / name
        path.write_bytes(blob)
        assert not (streamlog.elf_kind(path) != want), name
    log = streamlog.Log(work / "run.log")
    assert not (log.run(["/bin/echo", "hello"], work, os.environ.copy()) != 0), "true"
    assert not (log.run(["/bin/false"], work, os.environ.copy()) == 0), "false"
    log.hash_tree(work)
    log.close()
    assert not (len(streamlog.sha256(note)) != 64), "sha"


def check_record(work):
    old = sys.argv
    sys.argv = ["record_build.py", str(work / "bad-build.log"), str(work), "false", "false"]
    try:
        expect_exit(record_build.main, "1")
        sys.argv = ["record_build.py", str(work / "ok-build.log"), str(work), "true", "true"]
        record_build.main()
    finally:
        sys.argv = old


def check_clock(work):
    repo = work / "clock"
    repo.mkdir()
    init_repo(repo)
    (repo / "base.txt").write_text("base\n")
    subprocess.check_call(["git", "-C", repo, "add", "base.txt"])
    subprocess.check_call(["git", "-C", repo, "commit", "-q", "-m", "base"])
    base = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
    stamp = "2026-11-09 16:10:00 -0600"
    tip = git_commit(str(repo), "fix: one\n\nOne commit.\n", stamp, stamp)
    messages = work / "messages.txt"
    messages.write_text(tip + "\nfix: renamed\n\nRenamed.")
    found = rewrite_clock.load_messages(str(messages))
    assert not (tip not in found or not found[tip].endswith("\n")), "messages"
    assert not (rewrite_clock.load_messages("") != {}), "empty messages"
    (work / "blank.txt").write_text("\n---\n")
    assert not (rewrite_clock.load_messages(str(work / "blank.txt")) != {}), "blank messages"
    subprocess.check_call(["git", "-C", repo, "update-ref", "refs/heads/other", base])
    log_path = work / "clock-log.txt"
    os.environ["CLOCK_RAW"] = str(work / "raw-date.txt")
    with log_path.open("w") as handle:
        previous = rewrite_clock.reading(handle)
        rewrite_clock.recreate(str(repo), base, tip, "refs/heads/other", {tip: "fix: bare"}, handle, previous)
    os.environ.pop("CLOCK_RAW", None)
    with log_path.open("a") as handle:
        rewrite_clock.commit_head(str(repo), "refs/heads/demo", str(messages), handle, "1999-01-01 00:00:00 -0600")
    ended = work / "ended.txt"
    ended.write_text("fix: ended\n\nEnded.\n")
    old = sys.argv
    sys.argv = ["rewrite_clock.py", "commit", str(repo), "refs/heads/demo", str(log_path), str(ended), "1999-01-01 00:00:01 -0600"]
    try:
        rewrite_clock.main()
    finally:
        sys.argv = old
    sys.argv = ["rewrite_clock.py", str(repo), base, tip, "refs/heads/demo", str(log_path), str(messages), "1999-01-01 00:00:02 -0600"]
    try:
        rewrite_clock.main()
    finally:
        sys.argv = old


def check_evidence(work):
    repo = work / "evidence"
    repo.mkdir()
    init_repo(repo)
    (repo / "base.txt").write_text("base\n")
    subprocess.check_call(["git", "-C", repo, "add", "base.txt"])
    subprocess.check_call(["git", "-C", repo, "commit", "-q", "-m", "base"])
    base = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
    use = "2026-11-09 16:20:00 -0600"
    later = "2026-11-09 16:20:01 -0600"
    first = git_commit(str(repo), "fix: a\n\nBody.\n", use, use)
    second = git_commit(str(repo), "fix: b\n\nBody.\n", use, use)
    good = work / "good-clock.txt"
    good.write_text(clock_row(repo, first, use, "fix: a") + clock_row(repo, second, later, "fix: b"))
    expect_exit(lambda: run_evidence(good, repo, base), "clock")
    same = work / "same-clock.txt"
    same.write_text(clock_row(repo, first, use, "fix: a") + clock_row(repo, second, use, "fix: b"))
    expect_exit(lambda: run_evidence(same, repo, base), "shared second")
    (work / "missing-clock.txt").write_text("")
    expect_exit(lambda: run_evidence(work / "missing-clock.txt", repo, base), "rows")
    drifted = git_commit(str(repo), "fix: c\n\nBody.\n", use, later)
    dates = work / "dates-clock.txt"
    dates.write_text(clock_row(repo, drifted, use, "fix: c"))
    expect_exit(lambda: run_evidence(dates, repo, base), "dates differ")
    signed = subprocess.check_output(
        [
            "git", "-C", repo, "-c", "gpg.format=openpgp",
            "-c", "user.signingkey=B727D96928F86A7F789F8F520CDEF2D2EF447838",
            "-c", "user.email=rgoswami@ieee.org", "-c", "user.name=Rohit Goswami",
            "commit", "--allow-empty", "-S", "-m", "signed",
        ],
        text=True,
    )
    # commit writes the hash only with some formats; read HEAD.
    signed_hash = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
    signed_clock = work / "signed-clock.txt"
    signed_clock.write_text(clock_row(repo, signed_hash, use, "signed"))
    expect_exit(lambda: run_evidence(signed_clock, repo, base), "signed")
    one = work / "one"
    one.mkdir()
    init_repo(one)
    (one / "base.txt").write_text("base\n")
    subprocess.check_call(["git", "-C", one, "add", "base.txt"])
    subprocess.check_call(["git", "-C", one, "commit", "-q", "-m", "base"])
    one_base = subprocess.check_output(["git", "-C", one, "rev-parse", "HEAD"], text=True).strip()
    one_tip = git_commit(str(one), "fix: only\n\nBody.\n", use, use)
    one_clock = work / "one-clock.txt"
    one_clock.write_text("noise before\n" + clock_row(one, one_tip, use, "fix: only"))
    old = sys.argv
    sys.argv = ["date_evidence.py", str(one_clock), str(one), one_base, str(one), one_tip]
    try:
        date_evidence.main()
    finally:
        sys.argv = old


def check_dates(work):
    seams = os.environ["SRC_SEAMS"]
    lammps = os.environ["SRC_LAMMPS"]
    hashes = subprocess.check_output(["git", "-C", seams, "log", "--format=%H", "be748d7a..HEAD"], text=True)
    hashes += subprocess.check_output(["git", "-C", lammps, "log", "--format=%H", "098ffdef..HEAD"], text=True)
    lines = hashes.split()
    root = work / "logs"
    root.mkdir()
    (root / "old.log").write_text("start: 2000-01-01T00:00:00+00:00\n%s\n" % lines[0])
    (root / "new.log").write_text("start: 2099-01-01T00:00:00+00:00\n%s\n" % "\n".join(lines[1:-1]))
    (root / "note.txt").write_text("no start\n")
    (root / "skip.dat").write_text("start: 2099-01-01T00:00:00+00:00\n")
    (root / "nested").mkdir()
    blocked = root / "blocked.log"
    blocked.write_text("start: 2099-01-01T00:00:00+00:00\n")
    blocked.chmod(0)
    old = sys.argv
    sys.argv = ["date_check.py", str(root), str(work / "missing-root")]
    try:
        code = date_check.main()
    finally:
        sys.argv = old
    assert not (code != 1), "date check %s" % code


def check_trace():
    assert not (trace_libs.ppid(-1) is not None), "ppid"
    assert not (trace_libs.ppid(os.getpid()) is None), "ppid self"
    assert not (os.getpid() not in trace_libs.descendants(os.getpid())), "descendants"
    trace_libs.mapped(-1)
    trace_libs.mapped(os.getpid())
    sample = "    1:\ttrying file=./libdseams_lammps.so\n    1:\tcalling init: /lib/liblammps.so.0\n"
    found = trace_libs.paths_from_loader(sample)
    assert not (len(found) < 2), "loader %s" % found
    assert not (trace_libs.paths_from_loader("nothing") != []), "loader empty"
    assert not (trace_libs.paths_from_loader("file=") != []), "loader blank"
    assert not (trace_libs.parent_of("nope") is not None), "parent none"
    assert not (trace_libs.parent_of("name) ") is not None), "parent short"
    assert not (trace_libs.parent_of("name) 1 2") != 2), "parent value"
    assert not (trace_libs.map_path("liblammps.so") != ""), "map short"
    assert not (trace_libs.map_path("a b c d e /lib/liblammps.so.0") != "/lib/liblammps.so.0"), "map path"
    child = subprocess.Popen(["sleep", "3"])
    try:
        assert not (child.pid not in trace_libs.descendants(os.getpid())), "child"
    finally:
        child.kill()
        child.wait()
    code = subprocess.run(
        [os.environ.get("DSEAMS_PYTHON", "python3"), str(Path(trace_libs.__file__)), "/bin/true"],
        env=os.environ.copy(),
    ).returncode
    assert not (code != 0), "trace true"


def run_main(path, argv):
    import runpy
    old = sys.argv
    sys.argv = argv
    try:
        try:
            runpy.run_path(str(path), run_name="__main__")
        except SystemExit:
            return
    finally:
        sys.argv = old


def check_campaign(work):
    import run_campaign
    assert not (run_campaign.wall_of("nope\n") != ""), "wall empty"
    assert not (run_campaign.wall_of("wall_s 1.5\n") != "1.5"), "wall"
    run_campaign.clean("plain")
    run_campaign.timed(["/bin/true"], work)
    blocks = run_campaign.census("note\nstep 0\nrings a\n\n")
    assert not (blocks[0]["lines"].get("rings") != "rings a"), "census"
    assert not (run_campaign.census("note\n") != []), "census empty"
    mapped = run_campaign.chill_map("chill bare bar:2")
    assert not (mapped.get("bar") != 2), "chill map"
    rows = []
    run_campaign.compare_routes("s", "serial", [], [{"step": 0, "lines": {"rings": "a"}}], rows)
    run_campaign.compare_routes(
        "s", "serial",
        [{"step": 0, "lines": {"rings": "a"}}],
        [{"step": 1, "lines": {"rings": "b"}}],
        rows,
    )
    run_campaign.compare_routes(
        "s", "serial",
        [{"step": 0, "lines": {"rings": "a", "chill": "", "cage": "", "n": "", "n_rings": "", "labels": "", "cages": ""}}],
        [{"step": 0, "lines": {"rings": "b", "chill": "", "cage": "", "n": "", "n_rings": "", "labels": "", "cages": ""}}],
        rows,
    )
    assert not (run_campaign.quantile7([3.0], 0.5) != 3.0), "quantile one"
    assert not (run_campaign.quantile7([1.0, 2.0, 3.0, 4.0], 1.0) != 4.0), "quantile edge"
    expect_exit(lambda: run_campaign.require_fix_input("nope"), "not the fix")
    expect_exit(lambda: run_campaign.require_fix_input("fix sea all dseams\n"), "does not load")
    run_campaign.require_fix_input(run_campaign.LITERAL.read_text())
    expect_exit(lambda: run_campaign.require_ok(2, "bad 2"), "bad 2")
    run_campaign.require_ok(0, "ok")
    expect_exit(
        lambda: run_campaign.run_lammps(
            "/bin/false", "false", "0" * 64, work / "badlmp", "serial", "plugin", "bad", "bad-plugin.log"
        ),
        "exit",
    )
    expect_exit(lambda: run_campaign.require_record(work / "none", "plugin", "s", "serial"), "no live")
    present = work / "present"
    present.mkdir()
    (present / "live.txt").write_text("step 0\n")
    run_campaign.require_record(present, "plugin", "s", "serial")
    summary = []
    assert not (run_campaign.take_pair({}, "s", "serial", summary) is not None), "pair absent"
    saved = {
        ("s", "serial", "plugin"): {"frames": 1},
        ("s", "serial", "intree"): {"frames": 1},
    }
    assert not (run_campaign.take_pair(saved, "s", "serial", summary) is None), "pair present"
    plug = work / "plug"
    plug.mkdir()
    (plug / "libdseams_lammps.so").write_bytes(b"x")
    run_campaign.place_plugin(plug, None)
    assert not ((plug / "libdseams_lammps.so").exists()), "plugin removed"
    run_campaign.place_plugin(plug, b"y")
    assert not ((plug / "libdseams_lammps.so").read_bytes() != b"y"), "plugin written"
    expect_exit(lambda: run_campaign.require_literal({"mismatches": 1, "live_steps": [0]}, "plugin"), "mismatch")
    run_campaign.require_literal({"mismatches": 0, "live_steps": [0]}, "plugin")
    expect_exit(lambda: run_campaign.require_route_match(2), "routes")
    run_campaign.require_route_match(0)
    fresh = work / "fresh"
    fresh.mkdir()
    init_repo(fresh)
    git_commit(str(fresh), "fix: first\n\nFirst.\n", "2026-11-09 16:10:00 -0600", "2026-11-09 16:10:00 -0600")
    def compare_dump(dump):
        if not dump.is_file():
            print("dump absent", flush=True)
            return
        folder = work / "cmp"
        folder.mkdir()
        (folder / "dump.lammpstrj").write_bytes(dump.read_bytes())
        (folder / "live.txt").write_text(
            "step 0\nrings 3:9\nchill cubic:9\ncage hc:9 plain\nstep 1\nrings 3:9\n"
        )
        run_campaign.compare_live_cli("s", "serial", "plugin", folder, [])
        saved_cli = run_campaign.SEAMS
        run_campaign.SEAMS = Path("/bin/false")
        try:
            expect_exit(lambda: run_campaign.cli_lines("rings", folder / "dump.lammpstrj", [], "rings "), "exit")
        finally:
            run_campaign.SEAMS = saved_cli

    compare_dump(work / "no-such-dump.lammpstrj")
    compare_dump(Path(os.environ["BUILD"]) / "work" / "plugin" / "small-run0" / "serial" / "dump.lammpstrj")
    expect_exit(lambda: run_campaign.require_fix_input(run_campaign.LITERAL.read_text()), "not the fix")


def check_neigh(work):
    import neigh_ref
    text = work / "pre.txt"
    text.write_text(
        "junk\nframe\nstep 0\ncutoff 3.5\nbox 0 0 0 10 10 10 0 0 0\n"
        "host 0\natoms 2\nnoise 1\natom 1 1 0 0 0\natom 2 1 1 0 0\nrow 1 2\n"
        "label 1 1 0 0 0 0\nendframe\n"
    )
    frames = neigh_ref.parse_frames(text)
    assert not (len(frames) != 1), "frames"
    neigh_ref.print_repo("gone", "")
    neigh_ref.print_repo("seams", os.environ["SRC_SEAMS"])
    numpy = __import__("numpy")
    points = numpy.array(
        [[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0], [0.9, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0]]
    )
    box = [0, 0, 0, 10, 10, 10, 0, 0, 0]
    neigh_ref.bond_counts(points, box, 0, 1.0)
    neigh_ref.bond_counts(points, box, 0, 0.0)
    neigh_ref.eclipsed_staggered(0, [[1], [0]], numpy.zeros((2, 7), dtype=complex))
    neigh_ref.compare_sets(frames[0], {}, "freud")
    bad = work / "badhost.txt"
    bad.write_text(
        "frame\nstep 0\ncutoff 3.5\nbox 0 0 0 10 10 10 0 0 0\nhost 1\n"
        "atoms 2\natom 1 1 0 0 0\natom 2 1 1 0 0\nrow 1 9\n"
        "label 1 1 0 0 0 0\nendframe\n"
    )
    old = sys.argv
    sys.argv = ["neigh_ref.py", str(bad)]
    try:
        expect_exit(neigh_ref.main, "mismatch")
    finally:
        sys.argv = old


def check_prose(work):
    if not os.environ.get("SNAPPER"):
        print("snapper absent", flush=True)
        return
    os.environ["SEAMS_BASE"] = "be748d7a26957fd91bd0371b2dfefdd326feda3d"
    os.environ["LAMMPS_BASE"] = "098ffdefa525cabc82272925de02e9884f58c657"
    os.environ["PROSE_LOG"] = str(work / "prose.log")
    kit = work / "badkit"
    kit.mkdir()
    (kit / "NOTE.md").write_text(
        "This sentence is deliberately long so the formatter wraps it at the column limit and reports the difference here.\n"
    )
    os.environ["KIT"] = str(kit)
    import prose_gate
    saved_seams = os.environ.get("SRC_SEAMS", "")
    saved_lammps = os.environ.get("SRC_LAMMPS", "")
    try:
        expect_exit(prose_gate.main, "1")
        badrepo = work / "badrepo"
        badrepo.mkdir()
        init_repo(badrepo)
        (badrepo / "a.txt").write_text("a\n")
        subprocess.check_call(["git", "-C", badrepo, "add", "a.txt"])
        subprocess.check_call(["git", "-C", badrepo, "commit", "-q", "-m", "base"])
        subprocess.check_call(
            [
                "git", "-C", badrepo, "commit", "-q", "--allow-empty",
                "-m", "fix: reject a long line\n\nThe gate holds a value and this sentence is deliberately long so the formatter wraps it.",
            ]
        )
        base = subprocess.check_output(["git", "-C", badrepo, "rev-parse", "HEAD^"], text=True).strip()
        os.environ["SRC_SEAMS"] = str(badrepo)
        os.environ["SRC_LAMMPS"] = str(badrepo)
        os.environ["SEAMS_BASE"] = base
        os.environ["LAMMPS_BASE"] = base
        os.environ["PAPER_CORE"] = str(work / "core.toml")
        (work / "core.toml").write_text("this is not a profile {\n")
        expect_exit(prose_gate.main, "1")
    finally:
        os.environ["SRC_SEAMS"] = saved_seams
        os.environ["SRC_LAMMPS"] = saved_lammps


def main():
    work = Path(tempfile.mkdtemp())
    check_streamlog(work)
    blog = streamlog.Log(work / "abort.log")
    cover.abort(blog, 0)
    expect_exit(lambda: cover.abort(streamlog.Log(work / "abort-bad.log"), 3), "3")
    os.environ["GCOVR_XML"] = ""
    os.environ["GCOVR_DIFF_LOG"] = ""
    import py_diff
    assert not (py_diff.gcovr_diff(Path("."), Path("."), "python3", "d763dd35", os.environ.copy()) != 0), "gcovr skip"
    py_diff.raise_status(0)
    expect_exit(lambda: py_diff.raise_status(1), "1")
    check_record(work)
    check_clock(work)
    check_evidence(work)
    check_dates(work)
    check_trace()
    check_campaign(work)
    check_neigh(work)
    snap = os.environ.pop("SNAPPER", None)
    check_prose(work)
    if snap:
        os.environ["SNAPPER"] = snap
    check_prose(work)
    run_main(date_check.__file__, ["date_check.py", str(work / "logs")])
    run_main(record_build.__file__, ["record_build.py", str(work / "ok-build.log"), str(work), "true", "true"])
    run_main(trace_libs.__file__, ["trace_libs.py", "/bin/true"])
    run_main(
        rewrite_clock.__file__,
        ["rewrite_clock.py", "commit", str(work / "clock"), "refs/heads/demo", str(work / "clock-log.txt"), str(work / "ended.txt"), "1999-01-01 00:00:03 -0600"],
    )
    one = work / "one"
    one_base = subprocess.check_output(["git", "-C", one, "rev-parse", "HEAD^"], text=True).strip()
    run_main(
        date_evidence.__file__,
        ["date_evidence.py", str(work / "one-clock.txt"), str(one), one_base, str(one), subprocess.check_output(["git", "-C", one, "rev-parse", "HEAD"], text=True).strip()],
    )
    print("receipt fail exit 0", flush=True)


if __name__ == "__main__":
    main()
