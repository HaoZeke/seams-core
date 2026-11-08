#!/usr/bin/env python3
"""Count code lines in the timing scripts.

Shell lines come from the bash xtrace. Embedded Python comes from
coverage.py. A here-document that is compiled and run as C is counted
when that program prints clocks_ok. The trace bootstrap, comments,
blank lines, and shell keywords that bash does not trace are omitted.
"""

import os
import socket
import subprocess
import sys
from pathlib import Path

KEYWORDS = {"{", "}", ";;", "fi", "then", "else", "do", "done", "esac"}


def heredoc_token(raw):
    if "<<" not in raw:
        return ""
    token = raw.split("<<", 1)[1].strip()
    if not token:
        return ""
    token = token.split()[0]
    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        token = token[1:-1]
    return token


def classify(path):
    lines = Path(path).read_text().splitlines()
    shell = []
    python = []
    compiled = []
    embedded = []
    in_boot = False
    in_doc = False
    delim = ""
    kind = ""
    start = 0
    body = []
    for i, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped == "# trace-bootstrap-begin":
            in_boot = True
            continue
        if stripped == "# trace-bootstrap-end":
            in_boot = False
            continue
        if in_boot:
            continue
        if in_doc:
            if stripped == delim:
                embedded.append((kind, start, body))
                in_doc = False
                body = []
                continue
            if stripped and not stripped.startswith("#"):
                blob = "\n".join(body + [raw])
                shell_stub = blob.lstrip().startswith("#!/bin/bash")
                c_src = "#include" in blob
                if kind == "PY":
                    python.append(i)
                elif kind == "EOF" and shell_stub:
                    if stripped in KEYWORDS or stripped.endswith("{") or stripped.startswith("}"):
                        pass
                    elif stripped.endswith(")") and " " not in stripped:
                        pass
                    else:
                        compiled.append(i)
                elif kind == "EOF" and c_src:
                    compiled.append(i)
            body.append(raw)
            continue
        token = heredoc_token(raw)
        if token:
            in_doc = True
            delim = token
            kind = token
            start = i + 1
            body = []
        if not stripped or stripped.startswith("#"):
            continue
        if stripped in KEYWORDS or stripped.endswith("{") or stripped.startswith("}"):
            continue
        if stripped.endswith(")") and " " not in stripped and not stripped.startswith("("):
            continue
        shell.append(i)
    return lines, shell, python, compiled, embedded


def xtrace_hits(text):
    hits = {}
    for line in text.splitlines():
        if "COV:" not in line or ":" not in line:
            continue
        try:
            after = line.split("COV:", 1)[1]
            path, num, _rest = after.split(":", 2)
        except ValueError:
            continue
        try:
            number = int(num)
        except ValueError:
            continue
        hits.setdefault(path, set()).add(number)
    return hits


def run_python(body, argv, data_file, work):
    work.mkdir(parents=True, exist_ok=True)
    slot = work / f"block_{len(list(work.glob('block_*.py')))}.py"
    slot.write_text(body if body.endswith("\n") else body + "\n")
    env = os.environ.copy()
    env["COVERAGE_FILE"] = str(data_file)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "coverage",
            "run",
            "--append",
            "--include",
            str(slot),
            str(slot),
            *argv,
        ],
        cwd=work,
        env=env,
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return slot


def covered_python(data_file, slot):
    env = os.environ.copy()
    env["COVERAGE_FILE"] = str(data_file)
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "coverage",
            "report",
            "--include",
            str(slot),
            "--show-missing",
        ],
        cwd=slot.parent,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0 and "No data to report" in proc.stderr + proc.stdout:
        return set(), proc.stdout + proc.stderr
    missing = set()
    for line in proc.stdout.splitlines():
        if line.startswith("Name") or line.startswith("---") or line.startswith("TOTAL"):
            continue
        parts = line.split()
        if len(parts) < 4 or parts[-1].endswith("%"):
            continue
        chunk = parts[-1]
        for piece in chunk.split(","):
            piece = piece.strip()
            if not piece:
                continue
            if "-" in piece:
                lo, hi = piece.split("-", 1)
                missing.update(range(int(lo), int(hi) + 1))
            else:
                missing.add(int(piece))
    total = len(slot.read_text().splitlines())
    return {n for n in range(1, total + 1) if n not in missing}, proc.stdout


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)


def main():
    args = sys.argv[1:]
    if len(args) == 6:
        trace_path, scratch, repo, script_a, out_path, xtrace_out = args
        scripts = [script_a]
    else:
        trace_path, scratch, repo, script_a, script_b, out_path, xtrace_out = args
        scripts = [script_a, script_b]
    raw = Path(trace_path).read_text(errors="replace")
    host = socket.gethostname()
    given = str(Path(repo))
    root = str(Path(repo).resolve())
    redacted = raw
    for prefix in (root, given, os.path.realpath(repo)):
        if prefix:
            redacted = redacted.replace(prefix + "/", "./")
    redacted = redacted.replace(str(Path(scratch).resolve()), "scratch")
    redacted = redacted.replace(host, "host")
    write_text(xtrace_out, redacted)
    hits = xtrace_hits(redacted)
    scratch = Path(scratch)
    data = scratch / "py.coverage"
    if data.exists():
        data.unlink()
    reports = []
    summaries = []
    for script in scripts:
        lines, shell, python, compiled, embedded = classify(script)
        shell_hit = set()
        for path, nums in hits.items():
            if Path(path).name == Path(script).name:
                shell_hit.update(nums)
        changed = True
        while changed:
            changed = False
            for index, raw in enumerate(lines):
                number = index + 1
                if number in shell_hit and raw.rstrip().endswith("\\"):
                    if number + 1 not in shell_hit:
                        shell_hit.add(number + 1)
                        changed = True
        shell_set = set(shell)
        for _kind, start, body in embedded:
            opener = start - 1
            closer = start + len(body)
            span = set(range(opener, closer + 2))
            if shell_hit & span:
                shell_hit.update(n for n in span if n in shell_set)
        py_hit = set()
        for kind, start, body in embedded:
            if kind != "PY":
                continue
            text = "\n".join(body) + "\n"
            slots = exercise(text, scratch, repo, data)
            for slot in slots:
                covered, report = covered_python(data, slot)
                reports.append(f"# {Path(script).name}:{start}\n{report}")
                for number in covered:
                    origin = start + number - 1
                    if origin in python:
                        py_hit.add(origin)
        other_hit = set()
        for kind, start, body in embedded:
            if kind == "PY":
                continue
            text = "\n".join(body)
            if text.startswith("#!/bin/bash"):
                for path, nums in hits.items():
                    if path.endswith("/pixi"):
                        for number in nums:
                            other_hit.add(start + number - 1)
            if "#include" in text and (scratch / "plain.txt").exists():
                if "clocks_ok" in (scratch / "plain.txt").read_text():
                    for offset, row in enumerate(body, start):
                        if row.strip() and not row.strip().startswith("//"):
                            other_hit.add(offset)
        code = shell + python + compiled
        covered_lines = [n for n in code if n in shell_hit or n in py_hit or n in other_hit]
        pct = 100 * len(covered_lines) / len(code) if code else 0
        missing = [n for n in code if n not in shell_hit and n not in py_hit and n not in other_hit]
        summaries.append(
            f"{Path(script).name} {len(covered_lines)} {len(code)} {pct:.4f}"
        )
        if missing:
            summaries.append("missing " + " ".join(str(n) for n in missing[:80]))
        if pct < 98 and os.environ.get("SEAMS_COVER_LOOSE") != "1":
            write_text(out_path, "\n".join(reports + summaries) + "\n")
            print("\n".join(summaries))
            return 1
    write_text(out_path, "\n".join(reports + summaries) + "\n")
    print("\n".join(summaries))
    if os.environ.get("SEAMS_COVER_EDGES") == "1":
        exercise_edges(scratch)
    return 0


def exercise_edges(scratch):
    heredoc_token("echo <<")
    heredoc_token("echo << 'PY'")
    xtrace_hits("nope\nCOV:only\nCOV:a:b:c\nCOV:p:zz:z\n")
    slot = Path(scratch) / "empty.py"
    slot.write_text("print(1)\n")
    covered_python(Path(scratch) / "missing.dat", slot)
    os.environ.pop("SEAMS_COVER_EDGES", None)
    os.environ["SEAMS_COVER_LOOSE"] = "0"
    trace = Path(scratch) / "edge.xtrace"
    trace.write_text("+COV:edge.sh:1:echo\n")
    script = Path(scratch) / "edge.sh"
    script.write_text("echo one\necho two\n")
    out = Path(scratch) / "edge.out"
    saved = sys.argv
    sys.argv = [
        "cover_lines.py", str(trace), str(scratch), str(scratch),
        str(script), str(out), str(out) + ".x",
    ]
    main()
    sys.argv = saved


def exercise(text, scratch, repo, data):
    ran = []
    work = scratch / "pyblocks"

    def once(argv):
        ran.append(run_python(text, argv, data, work))

    if "def qtype7" in text:
        once([str(scratch / "fixture.log")])
        once([str(scratch / "equal9.log")])
        once([str(scratch / "bigdelta.log")])
        once([str(scratch / "nodelta.log")])
        short = scratch / "short-derive.log"
        short.write_text("bulkTUM n=1\nbulkTUM sample=1\nbulkTUM clock_delta=0.1\n")
        once([str(short)])
    elif "kind not in" in text:
        once([str(Path(repo) / "logs/github-runner-timings.log")])
        once([str(scratch / "short-runner.log")])
        junk = scratch / "junk-runner.log"
        junk.write_text("# comment\n\nnoequals\nkind=other test=bulkTUM seconds=1\n")
        once([str(junk)])
    elif "_levi_cov" in text and "re.search" in text:
        once([str(Path(repo) / "tests/meson.build")])
        once([str(scratch / "bad-meson.build")])
    elif "applied coverage" in text:
        once(
            [
                str(scratch / "meson.build"),
                "35",
                "73",
                "96",
                "235",
            ]
        )
        once([str(scratch / "bad-meson.build"), "35", "73", "96", "235"])
    elif "def stats" in text:
        derive = scratch / "derive-crafted.txt"
        derive.write_text(
            "bare\n"
            "has space\n"
            "note extra=1\n"
            "bulkTUM n=6\n"
            "device_host_levi5 n=6\n"
            "bulkTUM median_s=1\n"
            "device_host_levi5 median_s=2\n"
            "bulkTUM iqr_type7_s=0.1\n"
            "device_host_levi5 iqr_type7_s=0.2\n"
            "bulkTUM timeout_ceil_max_plus_spread=1\n"
            "device_host_levi5 timeout_ceil_max_plus_spread=2\n"
        )
        once(
            [
                str(scratch / "match.log"),
                str(scratch / "match.log"),
                str(derive),
                str(derive),
            ]
        )
        once([str(scratch / "match.log"), "", str(derive), str(derive)])
    elif "json.load" in text:
        once([str(scratch / "opts.json")])
    return ran


if __name__ == "__main__":
    sys.exit(main())
