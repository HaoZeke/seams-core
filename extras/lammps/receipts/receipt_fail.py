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


def expect_exit(fn, needle):
    try:
        fn()
    except SystemExit as exc:
        if needle not in str(exc):
            raise SystemExit("wanted %s got %s" % (needle, exc))
        return
    raise SystemExit("no failure")


def exercise_expect():
    expect_exit(lambda: (_ for _ in ()).throw(SystemExit("needle")), "needle")
    try:
        expect_exit(lambda: None, "needle")
    except SystemExit as exc:
        if "no failure" not in str(exc):
            raise
    try:
        expect_exit(lambda: (_ for _ in ()).throw(SystemExit("other")), "needle")
    except SystemExit as exc:
        if "wanted" not in str(exc):
            raise


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
    if streamlog.clean("plain") != "plain":
        raise SystemExit("clean")
    if streamlog.clean("plain", [("zz", "nope")]) != "plain":
        raise SystemExit("pairs")
    note = work / "note.txt"
    note.write_text("text")
    os.environ["CVMFS_PREFIX"] = str(work)
    if "$CVMFS" not in streamlog.clean(str(work) + "/x"):
        raise SystemExit("cvmfs")
    os.environ.pop("CVMFS_PREFIX", None)
    (work / "nested").mkdir()
    (work / "nested" / "skip.o").write_bytes(b"obj")
    if streamlog.elf_kind(note) is not None or streamlog.elf_kind(work / "missing") is not None:
        raise SystemExit("elf text")
    (work / "short").write_bytes(b"\x7fELF")
    if streamlog.elf_kind(work / "short") is not None:
        raise SystemExit("elf short")
    for kind, name, want in ((2, "exec", "exec"), (3, "dyn", "dyn"), (1, "other", None)):
        blob = bytearray(b"\x7fELF\x01\x01" + b"\x00" * 12)
        blob[16:18] = kind.to_bytes(2, "little")
        path = work / name
        path.write_bytes(blob)
        if streamlog.elf_kind(path) != want:
            raise SystemExit(name)
    log = streamlog.Log(work / "run.log")
    if log.run(["/bin/echo", "hello"], work, os.environ.copy()) != 0:
        raise SystemExit("true")
    if log.run(["/bin/false"], work, os.environ.copy()) == 0:
        raise SystemExit("false")
    log.hash_tree(work)
    log.close()
    if len(streamlog.sha256(note)) != 64:
        raise SystemExit("sha")


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
    if tip not in found or not found[tip].endswith("\n"):
        raise SystemExit("messages")
    if rewrite_clock.load_messages("") != {}:
        raise SystemExit("empty messages")
    (work / "blank.txt").write_text("\n---\n")
    if rewrite_clock.load_messages(str(work / "blank.txt")) != {}:
        raise SystemExit("blank messages")
    subprocess.check_call(["git", "-C", repo, "update-ref", "refs/heads/other", base])
    log_path = work / "clock-log.txt"
    with log_path.open("w") as handle:
        previous = rewrite_clock.reading(handle)
        rewrite_clock.recreate(str(repo), base, tip, "refs/heads/other", {tip: "fix: bare"}, handle, previous)
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
    if code != 1:
        raise SystemExit("date check %s" % code)


def check_trace():
    if trace_libs.ppid(-1) is not None:
        raise SystemExit("ppid")
    if trace_libs.ppid(os.getpid()) is None:
        raise SystemExit("ppid self")
    if os.getpid() not in trace_libs.descendants(os.getpid()):
        raise SystemExit("descendants")
    trace_libs.mapped(-1)
    trace_libs.mapped(os.getpid())
    sample = "    1:\ttrying file=./libdseams_lammps.so\n    1:\tcalling init: /lib/liblammps.so.0\n"
    found = trace_libs.paths_from_loader(sample)
    if len(found) < 2:
        raise SystemExit("loader %s" % found)
    if trace_libs.paths_from_loader("nothing") != []:
        raise SystemExit("loader empty")
    if trace_libs.paths_from_loader("file=") != []:
        raise SystemExit("loader blank")
    code = subprocess.run(
        [os.environ.get("DSEAMS_PYTHON", "python3"), str(Path(trace_libs.__file__)), "/bin/true"],
        env=os.environ.copy(),
    ).returncode
    if code != 0:
        raise SystemExit("trace true")


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


def check_prose(work):
    if not os.environ.get("SNAPPER"):
        return
    os.environ["SEAMS_BASE"] = "be748d7a26957fd91bd0371b2dfefdd326feda3d"
    os.environ["LAMMPS_BASE"] = "098ffdefa525cabc82272925de02e9884f58c657"
    os.environ["PROSE_LOG"] = str(work / "prose.log")
    kit = work / "badkit"
    kit.mkdir()
    (kit / "NOTE.md").write_text("x" * 200 + "\n")
    os.environ["KIT"] = str(kit)
    import prose_gate
    expect_exit(prose_gate.main, "1")


def main():
    work = Path(tempfile.mkdtemp())
    exercise_expect()
    check_streamlog(work)
    blog = streamlog.Log(work / "abort.log")
    cover.abort(blog, 0)
    expect_exit(lambda: cover.abort(streamlog.Log(work / "abort-bad.log"), 3), "3")
    os.environ["GCOVR_XML"] = ""
    os.environ["GCOVR_DIFF_LOG"] = ""
    import py_diff
    if py_diff.gcovr_diff(Path("."), Path("."), "python3", "d763dd35", os.environ.copy()) != 0:
        raise SystemExit("gcovr skip")
    py_diff.raise_status(0)
    expect_exit(lambda: py_diff.raise_status(1), "1")
    check_record(work)
    check_clock(work)
    check_evidence(work)
    check_dates(work)
    check_trace()
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
