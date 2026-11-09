"""Unique physical line coverage from gcov and coverage.py JSON."""

import json
import re
import sys
from pathlib import Path


def physical(text):
    """Line numbers that contain code after comments are stripped."""
    lines = text.splitlines()
    keep = []
    block = False
    for number, line in enumerate(lines, start=1):
        stripped = line
        if block:
            end = stripped.find("*/")
            if end < 0:
                continue
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
        if stripped.strip():
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


def main():
    root = Path(sys.argv[1])
    rows = []
    for src in sys.argv[2:]:
        path = Path(src)
        if path.suffix == ".gcov":
            source = path.name
            if source.endswith(".gcov"):
                source = source[: -len(".gcov")]
            text = None
            for candidate in root.rglob(source):
                if candidate.is_file() and "bbdir" not in candidate.parts:
                    text = candidate.read_text(errors="replace")
                    rel = candidate.relative_to(root)
                    break
            if text is None:
                rows.append((source, 0, 0, 0.0))
                continue
            code = physical(text)
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
                    continue
                try:
                    rel = source.resolve().relative_to(root.resolve())
                except ValueError:
                    continue
                executed = set(info.get("executed_lines", []))
                missing = set(info.get("missing_lines", []))
                code = physical(source.read_text(errors="replace"))
                statements = (executed | missing) & code
                if not statements:
                    rows.append((str(rel), 0, 0, 1.0))
                    continue
                hit = statements & executed
                rows.append((str(rel), len(statements), len(hit), len(hit) / len(statements)))
        else:
            rows.append((src, 0, 0, 1.0))
    print("file executable covered fraction")
    worst = 1.0
    for name, denom, num, frac in rows:
        print("%s %d %d %.6f" % (name, denom, num, frac))
        worst = min(worst, frac)
    if worst < 0.98:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
