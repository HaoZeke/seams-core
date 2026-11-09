"""Physical line coverage from gcov text and coverage.py JSON.

Arguments are the source root, then one path per input.
Each input is a gcov file or a coverage.py JSON file.
A directory, a missing path, or any other file is a failed row.
A gcov row counts only when the source text embedded in the gcov
matches the file under the root. The same basename in another
directory is not a match.
A JSON row names each measured file. A missing path is a failed
row, and a path outside the root is a failed row.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

# Public history. A Python file changed since this commit must have a row.
PUBLIC_BASE = "7f310225"


def _python_code(line, quote):
    """Drop a hash comment. A string that contains a C marker stays code."""
    i = 0
    code = False
    while i < len(line):
        if quote:
            if line.startswith(quote, i):
                i += len(quote)
                quote = None
                code = True
                continue
            if len(quote) == 1 and line[i] == "\\":
                i += 2
                continue
            i += 1
            continue
        if line.startswith('"""', i) or line.startswith("'''", i):
            quote = line[i : i + 3]
            i += 3
            code = True
            continue
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            code = True
            continue
        if line[i] == "#":
            break
        if not line[i].isspace():
            code = True
        i += 1
    return code, quote


def _c_code(line, block):
    stripped = line
    if block:
        end = stripped.find("*/")
        if end < 0:
            return False, True
        stripped = stripped[end + 2 :]
        block = False
    while True:
        start = stripped.find("/*")
        if start < 0:
            break
        end = stripped.find("*/", start + 2)
        if end < 0:
            stripped = stripped[:start]
            block = True
            break
        stripped = stripped[:start] + stripped[end + 2 :]
    if "//" in stripped:
        stripped = stripped[: stripped.find("//")]
    return bool(stripped.strip()), block


def physical(text, suffix=".cpp"):
    """Line numbers that contain code after comments are stripped.

    Python uses hash comments. C and C++ use block and line comments.
    A bare call keeps the C rules.
    """
    lines = text.splitlines()
    keep = []
    block = False
    quote = None
    python = suffix == ".py"
    for number, line in enumerate(lines, start=1):
        if python:
            code, quote = _python_code(line, quote)
        else:
            code, block = _c_code(line, block)
        if code:
            keep.append(number)
    return set(keep)


def gcov_counts(path):
    counts = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        match = re.match(r"\s*([-0-9#=]+)\*?:\s*(\d+):", line)
        if not match:
            continue
        mark, number = match.group(1), int(match.group(2))
        if mark == "-":
            continue
        hit = 0 if "#" in mark or mark == "=====" else int(mark)
        counts[number] = max(counts.get(number, 0), hit)
    return counts


def embedded_source(path):
    """Source text stored on each gcov line, keyed by source line number."""
    lines = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        match = re.match(r"\s*[-0-9#=*]+:\s*(\d+):(.*)$", line)
        if not match:
            continue
        number = int(match.group(1))
        if number == 0:
            continue
        lines[number] = match.group(2)
    return lines


def source_matches(candidate, embedded):
    if not embedded:
        return False
    file_lines = candidate.read_text(errors="replace").splitlines()
    for number, text in embedded.items():
        if number < 1 or number > len(file_lines):
            return False
        if file_lines[number - 1] != text:
            return False
    return True


def matching_source(root, gcov_path, source_name):
    embedded = embedded_source(gcov_path)
    for candidate in root.rglob(source_name):
        if not candidate.is_file() or "bbdir" in candidate.parts:
            continue
        if source_matches(candidate, embedded):
            return candidate
    return None


def fail_row(name):
    return (str(name), 0, 0, 0.0)


def changed_python(root, base):
    """Python paths changed between base and HEAD. An empty list means no git history."""
    proc = subprocess.run(
        ["git", "diff", "--name-only", base, "HEAD", "--", "*.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return []
    return [line for line in proc.stdout.splitlines() if line.endswith(".py")]


def main():
    root = Path(sys.argv[1])
    rows = []
    for src in sys.argv[2:]:
        path = Path(src)
        if not path.is_file():
            rows.append(fail_row(path))
            continue
        if path.suffix == ".gcov":
            source = path.name
            if source.endswith(".gcov"):
                source = source[: -len(".gcov")]
            candidate = matching_source(root, path, source)
            if candidate is None:
                rows.append(fail_row(source))
                continue
            text = candidate.read_text(errors="replace")
            rel = candidate.relative_to(root)
            code = physical(text, candidate.suffix)
            counts = gcov_counts(path)
            covered = {n for n in code if counts.get(n, 0) > 0}
            # A physical line gcov marks as non-code stays covered.
            executable = {n for n in code if n in counts}
            missed = executable - covered
            # No executable line means the file has nothing gcov can count.
            if not executable:
                rows.append((str(rel), 0, 0, 1.0))
                continue
            rows.append(
                (
                    str(rel),
                    len(executable),
                    len(executable) - len(missed),
                    (len(executable) - len(missed)) / len(executable),
                )
            )
        elif path.suffix == ".json":
            data = json.loads(path.read_text())
            for name, info in data.get("files", {}).items():
                source = Path(name)
                if not source.is_file():
                    rows.append(fail_row(name))
                    continue
                try:
                    rel = source.resolve().relative_to(root.resolve())
                except ValueError:
                    rows.append(fail_row(name))
                    continue
                executed = set(info.get("executed_lines", []))
                missing = set(info.get("missing_lines", []))
                code = physical(source.read_text(errors="replace"), source.suffix)
                statements = (executed | missing) & code
                if not statements:
                    rows.append((str(rel), 0, 0, 1.0))
                    continue
                hit = statements & executed
                rows.append((str(rel), len(statements), len(hit), len(hit) / len(statements)))
        else:
            rows.append(fail_row(src))
    seen = {name for name, _denom, _num, _frac in rows}
    base = os.environ.get("DSEAMS_COVERAGE_BASE", PUBLIC_BASE)
    for rel in changed_python(root, base):
        if rel not in seen:
            rows.append(fail_row(rel))
    print("file executable covered fraction")
    worst = 1.0
    for name, denom, num, frac in rows:
        print("%s %d %d %.6f" % (name, denom, num, frac))
        worst = min(worst, frac)
    if worst < 0.98:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
