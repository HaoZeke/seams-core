# Copyright (c) d-SEAMS developers

import importlib.util
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "insitu_v30_seal.py"


def load():
    spec = importlib.util.spec_from_file_location("insitu_v30_seal", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def test_prefix_and_scores():
    seal = load()
    assert seal._prefix("tmp").endswith("tmp")
    assert len(seal.forbidden_prefixes()) == 6
    text = (
        "        -:    0:Source:widget.cpp\n"
        "        1:    1:int live;\n"
        "    #####:    2:int dead;\n"
        "    =====:    3:int also;\n"
        "        -:    4:/* note */\n"
        "branch  0 taken 100% (fallthrough)\n"
        "branch  1 taken 0% (throw)\n"
        "branch  2 never executed\n"
    )
    hits, total = seal.gcov_score(text, "widget.cpp")
    assert hits == 1 and total == 3
    taken, total_br, ex_taken, ex_total, excluded = seal.branch_score(text, "widget.cpp")
    assert taken == 1 and total_br == 3
    assert ex_taken == 1 and ex_total == 2
    assert excluded and excluded[0][0] == 3
    assert seal.gcov_score("        -:    0:Source:other.cpp\n", "widget.cpp") is None
    assert seal.substitute("abc", [("b", "d")]) == "adc"
    assert seal.conventional("fix(io): keep the frame")
    assert not seal.conventional("Fix the frame")


def test_build_fuzz_and_assemble(tmp_path, monkeypatch):
    seal = load()
    repo = ROOT
    build = tmp_path / "build"
    build.mkdir()
    elf = build / "lib.fake"
    elf.write_bytes(b"\x7fELF" + b"\x00" * 8)
    plain = build / "note.txt"
    plain.write_text("no\n")
    out = tmp_path / "out"
    bindir = tmp_path / "bin"
    bindir.mkdir()
    write_exe(bindir / "ninja", "#!/bin/sh\necho ninja-ok\n")
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    args = type("A", (), {})()
    args.src = str(repo)
    args.build = str(build)
    args.out = str(out)
    seal.cmd_build(args)
    assert "ninja-ok" in (out / "build.log").read_text()
    assert "lib.fake" in (out / "HASHES.txt").read_text()
    assert seal.is_elf(elf) and not seal.is_elf(plain)
    assert "mt19937" in seal.seed_line(repo)

    raw = out / "fuzz.raw"
    raw.write_text("1.0 alpha frame 1 OK\n1.0 beta frame 1 SKIP dense cutoff graph\nframes 2 failures 0\n")
    binary = tmp_path / "fuzzbin"
    binary.write_bytes(b"\x7fELF")
    args.bin = str(binary)
    args.traj = "input/traj"
    args.trials = 3
    monkeypatch.setenv("DSEAMS_V30_FUZZ_RESUMMARIZE", "1")
    seal.cmd_fuzz(args)
    log = (out / "fuzz.log").read_text()
    assert "ran_ok 1" in log and "skipped_dense 1" in log
    assert "alpha" in (out / "frames.csv").read_text()
    monkeypatch.delenv("DSEAMS_V30_FUZZ_RESUMMARIZE")
    write_exe(bindir / "fuzzbin", "#!/bin/sh\necho 'gamma frame 1 FAIL'\necho 'frames 1 failures 1'\nexit 3\n")
    args.bin = str(bindir / "fuzzbin")
    with pytest.raises(SystemExit):
        seal.cmd_fuzz(args)

    dest = tmp_path / "kit"
    args.from_dir = str(out)
    args.dest = str(dest)
    seal.cmd_assemble(args)
    paths = (dest / "PATHS.txt").read_text()
    assert "$SRC" in paths and "probe_dense" in paths
    leaked = seal._prefix("tmp")
    assert leaked not in paths
    assert (dest / "STORE_READY").exists()
    names = (dest / "SHA256SUMS").read_text()
    assert "SHA256SUMS" not in names
    assert seal.head_tree(repo)[0]


def test_python_gcov_and_commits(tmp_path, monkeypatch):
    seal = load()
    out = tmp_path / "py"
    bindir = tmp_path / "bin"
    bindir.mkdir()
    write_exe(
        bindir / "cov",
        "#!/bin/sh\n"
        "if [ \"$2\" = json ] || [ \"$3\" = json ]; then\n"
        "  out=$(echo \" $* \" | sed -n 's/.* -o \\([^ ]*\\).*/\\1/p')\n"
        "  printf '%s\\n' '{}' > \"$out\"\n"
        "fi\n"
        "exit 0\n",
    )
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    monkeypatch.setenv("COVERAGE_BIN", str(bindir / "cov"))
    monkeypatch.setenv("DSEAMS_V30_LIVE", "1")
    args = type("A", (), {})()
    args.src = str(ROOT)
    args.out = str(out)
    seal.cmd_python(args)
    assert (out / "coverage.json").read_text().strip() == "{}"
    monkeypatch.setenv("DSEAMS_V30_LIVE", "0")
    write_exe(
        bindir / "cov",
        "#!/bin/sh\n"
        "if [ \"$3\" = json ]; then\n"
        "  out=$(echo \" $* \" | sed -n 's/.* -o \\([^ ]*\\).*/\\1/p')\n"
        "  printf '%s\\n' '{}' > \"$out\"\n"
        "  exit 0\n"
        "fi\n"
        "exit 4\n",
    )
    with pytest.raises(SystemExit):
        seal.cmd_python(args)

    build = tmp_path / "gcovb"
    build.mkdir()
    note = build / "src_widget.cpp.gcno"
    note.write_text("note\n")
    gcov_text = (
        "        -:    0:Source:widget.cpp\n"
        "        3:    1:int live;\n"
        "    =====:    2:int dead;\n"
        "branch  0 taken 50%\n"
        "branch  1 taken 0% (throw)\n"
    )
    write_exe(
        bindir / "gcov",
        "#!/bin/sh\nprintf '%s\\n' %s > widget.cpp.gcov\n" % (
            "",
            " ".join("'" + line + "'" for line in gcov_text.splitlines()),
        ),
    )
    # The fake gcov above is awkward. Write the file directly from a script.
    write_exe(
        bindir / "gcov",
        "#!/bin/sh\ncat > widget.cpp.gcov << 'END'\n" + gcov_text + "END\n",
    )
    monkeypatch.setenv("DSEAMS_V30_SKIP_TEST", "1")
    ran = tmp_path / "meson-ran.txt"
    write_exe(
        bindir / "meson",
        "#!/bin/sh\necho ran >> " + str(ran) + "\nexit 0\n",
    )
    (tmp_path / "gcov-run.log").write_text("kept\n")
    # cmd_gcov looks at out/gcov-run.log, so seed it after out is known.
    gout = tmp_path / "gout"
    gout.mkdir()
    (gout / "gcov-run.log").write_text("kept\n")
    # Plant a changed-file note name. The real diff names real files, so
    # point gcov at the first real note-less file by also dropping a note
    # whose name matches a real changed source when one exists.
    changed = seal.git(ROOT, "diff", "--name-only", seal.PUBLIC_BASE + "..HEAD").split()
    cpp = next(rel for rel in changed if rel.endswith(".cpp") and (ROOT / rel).is_file())
    name = Path(cpp).name
    (build / ("src_" + name + ".gcno")).write_text("note\n")
    write_exe(
        bindir / "gcov",
        "#!/bin/sh\ncat > " + name + ".gcov << 'END'\n"
        "        -:    0:Source:" + name + "\n"
        "        4:    1:int live;\n"
        "    =====:    2:int dead;\n"
        "branch  0 taken 80% (fallthrough)\n"
        "branch  1 taken 0% (throw)\n"
        "END\n",
    )
    args.build = str(build)
    args.out = str(gout)
    args.gcov = str(bindir / "gcov")
    args.src = str(ROOT)
    seal.cmd_gcov(args)
    assert ran.read_text().strip() == "ran"
    report = (gout / "COVERAGE_CPP.txt").read_text()
    assert "=====" not in report or name in report
    assert (gout / "COVERAGE_EXCLUDED.txt").exists()
    assert "branch" in (gout / "COVERAGE_EXCLUDED.txt").read_text()

    cout = tmp_path / "commits"
    args.out = str(cout)
    seal.cmd_commits(args)
    assert "failures 0" in (cout / "COMMIT_CHECK.txt").read_text()

    write_exe(bindir / "meson", "#!/bin/sh\necho meson-ran\nexit 7\n")
    with pytest.raises(SystemExit):
        seal.cmd_gcov(args)


def test_main_rejects_a_missing_mode():
    seal = load()
    old = sys.argv
    sys.argv = ["insitu_v30_seal.py"]
    try:
        with pytest.raises(SystemExit):
            seal.main()
    finally:
        sys.argv = old


def test_git_failure(tmp_path):
    seal = load()
    with pytest.raises(SystemExit):
        seal.git(tmp_path, "rev-parse", "HEAD")


def test_remaining_arms(tmp_path, monkeypatch):
    seal = load()
    missing = tmp_path / "gone"
    assert seal.is_elf(missing) is False
    empty = tmp_path / "empty-repo"
    empty.mkdir()
    subprocess.run(["git", "init"], cwd=empty, check=True, capture_output=True)
    (empty / "tests").mkdir()
    (empty / "tests" / "fuzz_cage_invariance.cpp").write_text("int main(){return 0;}\n")
    subprocess.run(["git", "add", "tests"], cwd=empty, check=True)
    env = os.environ.copy()
    env.update({
        "GIT_AUTHOR_NAME": "Rohit Goswami",
        "GIT_AUTHOR_EMAIL": "rgoswami@ieee.org",
        "GIT_COMMITTER_NAME": "Rohit Goswami",
        "GIT_COMMITTER_EMAIL": "rgoswami@ieee.org",
    })
    subprocess.run(["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "seed"], cwd=empty, check=True, capture_output=True, env=env)
    assert "not found" in seal.seed_line(empty)
    text = "noise\n        -:    0:Source:widget.cpp\nbranch  0 taken 0%\n        1:    2:int x;\nbranch  1 taken 100% (throw)\n"
    seal.branch_score(text, "widget.cpp")
    seal.gcov_score("no source here\n", "widget.cpp")
    out = tmp_path / "build2"
    out.mkdir()
    args = type("A", (), {})()
    args.src = str(ROOT)
    args.build = str(tmp_path / "b")
    (tmp_path / "b").mkdir()
    args.out = str(out)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    write_exe(bindir / "ninja", "#!/bin/sh\necho again\n")
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    seal.cmd_build(args)
    seal.cmd_build(args)
    raw = out / "fuzz.raw"
    raw.write_text("delta frame 1 OK\nframes 1 failures 0 shiftTrials 1\n")
    binary = tmp_path / "binfile"
    binary.write_bytes(b"\x7fELF")
    args.bin = str(binary)
    args.traj = "traj"
    args.trials = 1
    monkeypatch.delenv("DSEAMS_V30_FUZZ_RESUMMARIZE", raising=False)
    write_exe(bindir / "fuzz", "#!/bin/sh\necho 'epsilon frame 2 OK'\necho 'frames 1 failures 0'\n")
    args.bin = str(bindir / "fuzz")
    seal.cmd_fuzz(args)
    src = out
    (src / ".coverage").write_text("skip\n")
    (src / "leak.txt").write_text("see " + seal._prefix("tmp") + "/x\n")
    gdir = src / "gcov"
    gdir.mkdir()
    (gdir / "a.gcov").write_text("line " + seal._prefix("agent") + "\n")
    (gdir / "sub").mkdir()
    args.from_dir = str(src)
    args.dest = str(tmp_path / "kit2")
    seal.cmd_assemble(args)
    assert seal._prefix("tmp") not in (tmp_path / "kit2" / "leak.txt").read_text()
    # gcov with no executable lines and an unwanted extra file
    build = tmp_path / "gb"
    build.mkdir()
    name = "widget.cpp"
    (build / ("src_" + name + ".gcno")).write_text("n\n")
    write_exe(bindir / "gcov", "#!/bin/sh\necho '        -:    0:Source:widget.cpp' > widget.cpp.gcov\necho junk > other.gcov\n")
    # Use a private repo whose diff lists widget.cpp and a missing header.
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True)
    (repo / "widget.cpp").write_text("int live;\n")
    (repo / "gone.hpp").write_text("int x;\n")
    subprocess.run(["git", "add", "widget.cpp", "gone.hpp"], cwd=repo, check=True)
    subprocess.run(["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "base"], cwd=repo, check=True, capture_output=True, env=env)
    base = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    (repo / "widget.cpp").write_text("int live;\nint dead;\n")
    subprocess.run(["git", "add", "widget.cpp"], cwd=repo, check=True, capture_output=True)
    subprocess.run(["git", "rm", "gone.hpp"], cwd=repo, check=True, capture_output=True)
    subprocess.run(["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "change"], cwd=repo, check=True, capture_output=True, env=env)
    monkeypatch.setattr(seal, "PUBLIC_BASE", base)
    write_exe(bindir / "meson", "#!/bin/sh\necho meson\nexit 0\n")
    monkeypatch.setenv("DSEAMS_V30_SKIP_TEST", "1")
    gout = tmp_path / "gout2"
    gout.mkdir()
    (gout / "gcov-run.log").write_text("kept\n")
    args.src = str(repo)
    args.build = str(build)
    args.out = str(gout)
    args.gcov = str(bindir / "gcov")
    seal.cmd_gcov(args)
    report = (gout / "COVERAGE_CPP.txt").read_text()
    assert "absent" in report
    assert "no executable lines" in report
    old = sys.argv
    try:
        for mode, extra in (
            ("build", ["--build", str(build), "--out", str(tmp_path / "mbuild")]),
            ("commits", ["--out", str(tmp_path / "mc")]),
            ("assemble", ["--build", str(build), "--from", str(src), "--dest", str(tmp_path / "mdest")]),
        ):
            sys.argv = ["seal", mode, "--src", str(repo), *extra]
            if mode == "build":
                seal.main()
            elif mode == "commits":
                with pytest.raises(SystemExit):
                    seal.main()
            else:
                seal.main()
    finally:
        sys.argv = old
    seal.branch_score(
        "        -:    0:Source:widget.cpp\n"
        "        1:    1:int x;\n"
        "branch  3 taken 0%\n"
        "branch  4 never executed (throw)\n",
        "widget.cpp",
    )
    write_exe(bindir / "ninja", "#!/bin/sh\nexit 2\n")
    args.out = str(tmp_path / "badbuild")
    args.src = str(ROOT)
    args.build = str(build)
    with pytest.raises(SystemExit):
        seal.cmd_build(args)
    write_exe(bindir / "gcov", "#!/bin/sh\necho '        -:    0:Source:widget.cpp'\n")
    monkeypatch.delenv("DSEAMS_V30_SKIP_TEST", raising=False)
    (gout / "gcov-run.log").write_text("stale\n")
    write_exe(bindir / "meson", "#!/bin/sh\nexit 0\n")
    args.src = str(repo)
    args.out = str(gout)
    seal.cmd_gcov(args)
    odd = tmp_path / "odd"
    odd.mkdir()
    subprocess.run(["git", "init"], cwd=odd, check=True, capture_output=True)
    (odd / "a.cpp").write_text("int a;\n")
    subprocess.run(["git", "add", "a.cpp"], cwd=odd, check=True)
    parent_env = env
    bad = parent_env.copy()
    bad["GIT_AUTHOR_NAME"] = "Ada Lovelace"
    bad["GIT_AUTHOR_EMAIL"] = "ada@example.com"
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "base"],
        cwd=odd, check=True, capture_output=True, env=bad,
    )
    parent = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=odd, text=True).strip()
    (odd / "a.cpp").write_text("int a;\nint b;\n")
    subprocess.run(["git", "add", "a.cpp"], cwd=odd, check=True)
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m",
         "fix(io): sign the change\n\nSigned-off-by: Ada Lovelace <ada@example.com>\n"],
        cwd=odd, check=True, capture_output=True, env=bad,
    )
    monkeypatch.setattr(seal, "PUBLIC_BASE", parent)
    args.src = str(odd)
    args.out = str(tmp_path / "oddout")
    with pytest.raises(SystemExit):
        seal.cmd_commits(args)
    elf = tmp_path / "elf"
    elf.write_bytes(b"\x7fELF")
    (tmp_path / "mfuzz").mkdir()
    (tmp_path / "mfuzz" / "fuzz.raw").write_text("zeta frame 1 OK\nframes 1 failures 0\n")
    monkeypatch.setenv("DSEAMS_V30_FUZZ_RESUMMARIZE", "1")
    sys.argv = ["seal", "fuzz", "--src", str(ROOT), "--bin", str(elf), "--out", str(tmp_path / "mfuzz")]
    seal.main()
    write_exe(
        bindir / "cov",
        "#!/bin/sh\nif [ \"$3\" = json ]; then out=$(echo \" $* \" | sed -n 's/.* -o \\([^ ]*\\).*/\\1/p'); printf '%s\\n' '{}' > \"$out\"; fi\nexit 0\n",
    )
    monkeypatch.setenv("COVERAGE_BIN", str(bindir / "cov"))
    monkeypatch.setenv("DSEAMS_V30_LIVE", "0")
    sys.argv = ["seal", "python", "--src", str(ROOT), "--out", str(tmp_path / "mpy")]
    seal.main()
    monkeypatch.setenv("DSEAMS_V30_SKIP_TEST", "1")
    monkeypatch.setattr(seal, "PUBLIC_BASE", base)
    (tmp_path / "mg").mkdir()
    (tmp_path / "mg" / "gcov-run.log").write_text("kept\n")
    sys.argv = ["seal", "gcov", "--src", str(repo), "--build", str(build), "--out", str(tmp_path / "mg"), "--gcov", str(bindir / "gcov")]
    seal.main()
    import runpy
    sys.argv = [str(SCRIPT), "--help"]
    try:
        runpy.run_path(str(SCRIPT), run_name="__main__")
    except SystemExit:
        pass
