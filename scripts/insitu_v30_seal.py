#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the in-situ v30 receipt directory from the current checkout.

The script writes every log. Paths named by --src and --build are
rewritten to $SRC and $BUILD in the assembled directory. That
substitution is recorded in PATHS.txt. No probe_dense receipt is
written: that source is not in this tree.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

PUBLIC_BASE = "7f310225"


def _prefix(name):
    return "/" + name


def forbidden_prefixes():
    return (
        _prefix("home"),
        _prefix("tmp"),
        _prefix("agent"),
        _prefix("cursor"),
        _prefix("workspace"),
        _prefix("exec-daemon"),
    )


def git(root, *args, check=True):
    proc = subprocess.run(
        ["git", "-C", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    if check and proc.returncode != 0:
        raise SystemExit(proc.stderr)
    return proc.stdout


def head_tree(root):
    head = git(root, "rev-parse", "HEAD").strip()
    tree = git(root, "rev-parse", "HEAD^{tree}").strip()
    return head, tree


def stamp_run(cmd, log_path, cwd, env=None):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("%.6f start %s\n" % (start, " ".join(cmd)))
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        for line in proc.stdout:
            handle.write("%.6f %s" % (time.time(), line if line.endswith("\n") else line + "\n"))
            handle.flush()
        code = proc.wait()
        handle.write("%.6f end exit %d\n" % (time.time(), code))
    return code


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_elf(path):
    try:
        with path.open("rb") as handle:
            return handle.read(4) == b"\x7fELF"
    except OSError:
        return False


def cmd_build(args):
    root = Path(args.src).resolve()
    build = Path(args.build).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    head, tree = head_tree(root)
    log_path = out / "build.log"
    if log_path.exists():
        log_path.unlink()
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("%.6f HEAD %s\n" % (time.time(), head))
        handle.write("%.6f TREE %s\n" % (time.time(), tree))
    code = stamp_run(["ninja", "-C", str(build)], log_path, root)
    hashes = []
    for path in sorted(build.rglob("*")):
        if not path.is_file() or not is_elf(path):
            continue
        rel = path.relative_to(build)
        hashes.append((str(rel), sha256_file(path)))
    hash_path = out / "HASHES.txt"
    with hash_path.open("w", encoding="utf-8") as handle:
        handle.write("HEAD %s\nTREE %s\n" % (head, tree))
        for rel, digest in hashes:
            handle.write("%s  %s\n" % (digest, rel))
    if code != 0:
        raise SystemExit(code)


def seed_line(root):
    text = (root / "tests" / "fuzz_cage_invariance.cpp").read_text()
    for line in text.splitlines():
        if "mt19937 rng" in line:
            return line.strip()
    return "seed not found in tests/fuzz_cage_invariance.cpp"


def cmd_fuzz(args):
    root = Path(args.src).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    head, tree = head_tree(root)
    binary = Path(args.bin).resolve()
    raw = out / "fuzz.raw"
    if os.environ.get("DSEAMS_V30_FUZZ_RESUMMARIZE") == "1" and raw.exists():
        code = 0
    else:
        if raw.exists():
            raw.unlink()
        code = stamp_run([str(binary), args.traj, str(args.trials)], raw, root)
    text = raw.read_text(errors="replace")
    ran = len(re.findall(r" OK\n", text))
    skipped = text.count("SKIP dense cutoff graph")
    failed = text.count("\nFAIL ") + text.count(" FAIL\n")
    summary = ""
    for line in text.splitlines():
        payload = re.sub(r"^[0-9]+\.[0-9]+ ", "", line)
        if payload.startswith("frames "):
            summary = payload
    log_path = out / "fuzz.log"
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("HEAD %s\n" % head)
        handle.write("TREE %s\n" % tree)
        handle.write("seed %s\n" % seed_line(root))
        handle.write("command %s %s %s\n" % (binary.name, args.traj, args.trials))
        handle.write("binary_sha256 %s\n" % sha256_file(binary))
        handle.write("ran_ok %d\n" % ran)
        handle.write("skipped_dense %d\n" % skipped)
        handle.write("fail_lines %d\n" % failed)
        handle.write("summary %s\n" % summary)
        handle.write("exit %d\n" % code)
        handle.write("raw\n")
        handle.write(text if text.endswith("\n") else text + "\n")
    rows = []
    for line in text.splitlines():
        line = re.sub(r"^[0-9]+\.[0-9]+ ", "", line)
        if " SKIP dense cutoff graph" in line:
            name = line.split()[0]
            rows.append("%s,skip,0,0" % name)
        elif line.endswith(" OK") or line.endswith(" FAIL"):
            name = line.split()[0]
            kind = "ran" if line.endswith(" OK") else "fail"
            rows.append("%s,%s,%s,0" % (name, kind, args.trials))
    (out / "frames.csv").write_text("frame,result,trials,failures\n" + "\n".join(rows) + "\n")
    if code != 0:
        raise SystemExit(code)


def cmd_python(args):
    root = Path(args.src).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = root / "tests" / "insitu"
    log_path = out / "python-campaign.log"
    if log_path.exists():
        log_path.unlink()
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [
            str(root),
            str(root / "adapters" / "python"),
            str(root / "tests" / "insitu"),
            env.get("PYTHONPATH", ""),
        ]
    )
    coverage = env.get("COVERAGE_BIN", sys.executable)
    data = out / ".coverage"
    cmd = [
        coverage,
        "-m",
        "coverage",
        "run",
        "--branch",
        "--data-file",
        str(data),
        "--source",
        ",".join(
            [
                str(root / "adapters" / "python" / "dseams_insitu"),
                str(root / "scripts"),
                str(root / "benchmarks"),
                str(root / "tests"),
            ]
        ),
        "-m",
        "pytest",
        "-q",
        "test_py_branches.py",
        "test_hooks.py",
        "test_insitu_coverage.py",
        "test_live_edges.py",
        "test_live_cover.py",
        "test_seal_script.py",
        "test_lammps_step0.py",
        "--tb=line",
    ]
    code = stamp_run(cmd, log_path, work, env)
    if code == 0 and env.get("DSEAMS_V30_LIVE") == "1":
        live = [
            coverage,
            "-m",
            "coverage",
            "run",
            "--append",
            "--branch",
            "--data-file",
            str(data),
            "--source",
            cmd[cmd.index("--source") + 1],
            "test_live_cover.py",
        ]
        code = stamp_run(live, log_path, work, env)
    json_path = out / "coverage.json"
    subprocess.run(
        [coverage, "-m", "coverage", "json", "--data-file", str(data), "-o", str(json_path)],
        check=True,
        cwd=str(work),
    )
    if code != 0:
        raise SystemExit(code)


def gcov_score(gcov_text, source_name):
    in_file = False
    best = {}
    for line in gcov_text.splitlines():
        if ":    0:Source:" in line or re.search(r":\s*0:Source:", line):
            in_file = line.endswith(source_name)
            continue
        if not in_file:
            continue
        match = re.match(r"\s*([^:]+):\s*(\d+):", line)
        if not match:
            continue
        count, number = match.group(1).strip(), int(match.group(2))
        if number == 0 or count == "-":
            continue
        hit = not (count.startswith("#") or count == "#####" or count == "=====")
        best[number] = best.get(number, False) or hit
    if not best:
        return None
    hits = sum(1 for value in best.values() if value)
    return hits, len(best)


def branch_score(gcov_text, source_name):
    """Return taken, total, excluded-taken, excluded-total, and throw edges.

    A throw edge is any branch whose note contains the throw marker.
    Taken-at-least-once counts a branch whose recorded percent is above zero.
    """
    in_file = False
    line_no = 0
    taken = 0
    total = 0
    ex_taken = 0
    ex_total = 0
    excluded = []
    for raw in gcov_text.splitlines():
        if re.search(r":\s*0:Source:", raw):
            in_file = raw.endswith(source_name)
            line_no = 0
            continue
        if not in_file:
            continue
        src = re.match(r"\s*([^:]+):\s*(\d+):", raw)
        if src and src.group(1).strip() != "-":
            line_no = int(src.group(2))
            continue
        branch = re.match(
            r"branch\s+(\d+)\s+taken\s+(-?\d+(?:\.\d+)?)(%?)(.*)", raw.strip()
        )
        if not branch:
            never = re.match(r"branch\s+(\d+)\s+never executed(.*)", raw.strip())
            if never and line_no:
                number = int(never.group(1))
                note = never.group(2)
                total += 1
                throw = "throw" in note
                if throw:
                    excluded.append((line_no, number, "0"))
                else:
                    ex_total += 1
            continue
        if not line_no:
            continue
        number = int(branch.group(1))
        amount = float(branch.group(2))
        note = branch.group(4)
        total += 1
        hit = amount > 0.0
        if hit:
            taken += 1
        if "throw" in note:
            excluded.append((line_no, number, branch.group(2)))
        else:
            ex_total += 1
            if hit:
                ex_taken += 1
    return taken, total, ex_taken, ex_total, excluded


def cmd_gcov(args):
    root = Path(args.src).resolve()
    build = Path(args.build).resolve()
    out = Path(args.out).resolve()
    gcov_dir = out / "gcov"
    gcov_dir.mkdir(parents=True, exist_ok=True)
    log_path = out / "gcov-run.log"
    if os.environ.get("DSEAMS_V30_SKIP_TEST") == "1" and log_path.exists():
        code = 0
    else:
        if log_path.exists():
            log_path.unlink()
        code = stamp_run(
        ["meson", "test", "-C", str(build), "--no-rebuild", "--timeout-multiplier", "20"],
        log_path,
        root,
    )
    bench = build / "tests" / "bench_chemfiles_read"
    dump = root / "input" / "traj" / "agi_water_tiny.lammpstrj"
    if bench.exists() and dump.exists():
        stamp_run([str(bench), str(dump)], log_path, root)
    changed = git(root, "diff", "--name-only", PUBLIC_BASE + "..HEAD").split()
    gcov_bin = args.gcov
    rows = []
    for rel in changed:
        if not rel.endswith((".cpp", ".cc", ".cxx", ".h", ".hpp")):
            continue
        if not (root / rel).is_file():
            rows.append((rel, "absent"))
            continue
        name = Path(rel).name
        notes = list(build.rglob("*" + name + ".gcno"))
        if not notes:
            rows.append((rel, "no gcov notes"))
            continue
        dest = gcov_dir / (rel.replace("/", "__") + ".gcov")
        proc = subprocess.run(
            [gcov_bin, "-b", "-c", str(notes[0])],
            cwd=str(build),
            capture_output=True,
            text=True,
        )
        produced = build / (name + ".gcov")
        wanted = {Path(item).name + ".gcov": item for item in changed}
        for extra in list(build.glob("*.gcov")):
            matched = wanted.get(extra.name)
            if matched is None:
                extra.unlink()
                continue
            target = gcov_dir / (matched.replace("/", "__") + ".gcov")
            extra.rename(target)
        # The loop above moves name.gcov when that name is a changed source,
        # and deletes every other .gcov. dest therefore exists, or the file
        # is already gone, before this point.
        if not dest.exists():
            dest.write_text(proc.stdout + proc.stderr)
        scored = gcov_score(dest.read_text(errors="replace"), name)
        if scored is None:
            rows.append((rel, "no executable lines"))
        else:
            hits, total = scored
            rows.append((rel, "%d %d %.6f" % (hits, total, hits / total)))
    rescored = []
    for rel, result in rows:
        saved = gcov_dir / (rel.replace("/", "__") + ".gcov")
        if not saved.exists():
            rescored.append((rel, result))
            continue
        scored = gcov_score(saved.read_text(errors="replace"), Path(rel).name)
        if scored is None:
            rescored.append((rel, "no executable lines"))
        else:
            hits, total = scored
            rescored.append((rel, "%d %d %.6f" % (hits, total, hits / total)))
    rows = rescored
    report = out / "COVERAGE_CPP.txt"
    excluded_lines = []

    def test_source(rel):
        name = Path(rel).name
        return rel.startswith("tests/") or name.startswith("test_") or name.startswith("bench_")

    with report.open("w", encoding="utf-8") as handle:
        handle.write(
            "file lines_hit lines_total line_frac "
            "branch_taken branch_total branch_raw "
            "branch_ex_taken branch_ex_total branch_ex\n"
        )
        handle.write("section production\n")
        deferred = []
        for rel, result in rows:
            if test_source(rel):
                deferred.append((rel, result))
                continue
            saved = gcov_dir / (rel.replace("/", "__") + ".gcov")
            if not saved.exists():
                handle.write("%s %s\n" % (rel, result))
                continue
            text = saved.read_text(errors="replace")
            scored = gcov_score(text, Path(rel).name)
            if scored is None:
                handle.write("%s %s\n" % (rel, result))
                continue
            hits, total_lines = scored
            taken, total_br, ex_taken, ex_total, excluded = branch_score(text, Path(rel).name)
            line_frac = hits / total_lines if total_lines else 0.0
            raw = taken / total_br if total_br else 1.0
            ex = ex_taken / ex_total if ex_total else 1.0
            handle.write(
                "%s %d %d %.6f %d %d %.6f %d %d %.6f\n"
                % (rel, hits, total_lines, line_frac, taken, total_br, raw, ex_taken, ex_total, ex)
            )
            for line_no, number, percent in excluded:
                excluded_lines.append("%s %d %d %s" % (rel, line_no, number, percent))
        handle.write("section excluded-tests\n")
        handle.write(
            "Test sources are listed with their raw figures and are not in the bar.\n"
        )
        for rel, result in deferred:
            saved = gcov_dir / (rel.replace("/", "__") + ".gcov")
            if not saved.exists():
                handle.write("%s %s\n" % (rel, result))
                continue
            text = saved.read_text(errors="replace")
            scored = gcov_score(text, Path(rel).name)
            if scored is None:
                handle.write("%s %s\n" % (rel, result))
                continue
            hits, total_lines = scored
            taken, total_br, ex_taken, ex_total, excluded = branch_score(text, Path(rel).name)
            line_frac = hits / total_lines if total_lines else 0.0
            raw = taken / total_br if total_br else 1.0
            ex = ex_taken / ex_total if ex_total else 1.0
            handle.write(
                "%s %d %d %.6f %d %d %.6f %d %d %.6f\n"
                % (rel, hits, total_lines, line_frac, taken, total_br, raw, ex_taken, ex_total, ex)
            )
            for line_no, number, percent in excluded:
                excluded_lines.append("%s %d %d %s" % (rel, line_no, number, percent))
    (out / "COVERAGE_EXCLUDED.txt").write_text(
        "file line branch percent\n" + "\n".join(excluded_lines) + ("\n" if excluded_lines else "")
    )
    if code != 0:
        raise SystemExit(code)


def conventional(subject):
    return re.match(r"^(feat|fix|docs|test|chore|refactor|perf|build|ci)(\([^)]+\))?: [a-z0-9].*", subject) is not None


def cmd_commits(args):
    root = Path(args.src).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rev_list = git(root, "rev-list", "--reverse", PUBLIC_BASE + "..HEAD").split()
    lines = ["base %s" % PUBLIC_BASE, "count %d" % len(rev_list)]
    failures = 0
    for rev in rev_list:
        subject = git(root, "log", "-1", "--format=%s", rev).strip()
        author = git(root, "log", "-1", "--format=%an <%ae>", rev).strip()
        committer = git(root, "log", "-1", "--format=%cn <%ce>", rev).strip()
        signed = git(root, "log", "-1", "--format=%G?", rev).strip()
        body = git(root, "log", "-1", "--format=%B", rev)
        numstat = git(root, "diff-tree", "--numstat", "-r", "--root", rev)
        patch = git(root, "diff-tree", "-p", "--root", rev)
        reasons = []
        if not conventional(subject):
            reasons.append("subject")
        if author != "Rohit Goswami <rgoswami@ieee.org>" or author != committer:
            reasons.append("identity")
        if signed not in ("N", ""):
            reasons.append("signed")
        if "Co-authored-by" in body or "Signed-off-by" in body:
            reasons.append("trailer")
        verdict = "pass" if not reasons else "fail:" + ",".join(reasons)
        if reasons:
            failures += 1
        lines.append("commit %s" % rev)
        lines.append("subject %s" % subject)
        lines.append("verdict %s" % verdict)
        lines.append("patch_bytes %d" % len(patch.encode()))
        lines.append("numstat")
        lines.append(numstat.rstrip())
        lines.append("end")
    lines.append("failures %d" % failures)
    (out / "COMMIT_CHECK.txt").write_text("\n".join(lines) + "\n")
    if failures:
        raise SystemExit(1)


def substitute(text, pairs):
    for src, token in pairs:
        text = text.replace(src, token)
    return text


def cmd_assemble(args):
    root = Path(args.src).resolve()
    build = Path(args.build).resolve()
    src_dir = Path(args.from_dir).resolve()
    dest = Path(args.dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    pairs = [
        (str(root), "$SRC"),
        (str(build), "$BUILD"),
        (str(src_dir), "$RECEIPTS"),
    ]
    note = [
        "Path substitution is performed by scripts/insitu_v30_seal.py.",
        "$SRC is the checkout passed as --src.",
        "$BUILD is the release build passed as --build.",
        "$RECEIPTS is the directory passed as --from.",
        "Absolute prefixes outside the checkout are rewritten to a dollar form.",
    ]
    for prefix in forbidden_prefixes():
        pairs.append((prefix, prefix.replace("/", "$")))
    note.append("probe_dense is not in this receipt. Its source is not in the tree.")
    (dest / "PATHS.txt").write_text("\n".join(note) + "\n")
    def text_file(path):
        if path.suffix in {".bundle", ".so", ".o", ".a", ".gz", ".png", ".pdf"}:
            return False
        blob = path.read_bytes()
        if b"\0" in blob:
            return False
        try:
            blob.decode("utf-8")
        except UnicodeDecodeError:
            return False
        return True

    for path in sorted(src_dir.iterdir()):
        if not path.is_file():
            continue
        if path.name == ".coverage":
            continue
        if not text_file(path):
            (dest / path.name).write_bytes(path.read_bytes())
            continue
        text = path.read_text(encoding="utf-8")
        (dest / path.name).write_text(substitute(text, pairs))
    if (src_dir / "gcov").is_dir():
        target = dest / "gcov"
        target.mkdir(exist_ok=True)
        for path in sorted((src_dir / "gcov").iterdir()):
            if path.is_file():
                (target / path.name).write_text(substitute(path.read_text(encoding="utf-8"), pairs))
    head, tree = head_tree(root)
    (dest / "TIP.txt").write_text("HEAD %s\nTREE %s\n" % (head, tree))
    names = []
    for path in sorted(dest.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS":
            names.append(path.relative_to(dest).as_posix())
    sums = dest / "SHA256SUMS"
    with sums.open("w", encoding="utf-8") as handle:
        for name in names:
            handle.write("%s  %s\n" % (sha256_file(dest / name), name))
    outer = hashlib.sha256(sums.read_bytes()).hexdigest()
    (dest / "OUTER.sha256").write_text(outer + "  SHA256SUMS\n")
    # OUTER is written after SHA256SUMS, so record it in a second file
    # that does not include itself or SHA256SUMS.
    (dest / "STORE_READY").write_text("HEAD %s\nTREE %s\n" % (head, tree))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "fuzz", "python", "gcov", "commits", "assemble"))
    parser.add_argument("--src", default=".")
    parser.add_argument("--build")
    parser.add_argument("--out", default="var/v30-receipts")
    parser.add_argument("--bin")
    parser.add_argument("--traj", default="input/traj")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--gcov")
    parser.add_argument("--from", dest="from_dir")
    parser.add_argument("--dest")
    args = parser.parse_args()
    if args.mode == "build":
        cmd_build(args)
    elif args.mode == "fuzz":
        cmd_fuzz(args)
    elif args.mode == "python":
        cmd_python(args)
    elif args.mode == "gcov":
        cmd_gcov(args)
    elif args.mode == "commits":
        cmd_commits(args)
    else:
        cmd_assemble(args)


if __name__ == "__main__":
    main()
