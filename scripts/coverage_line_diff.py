"""Compare gcov branch records between two builds of the same sources.

Each executable source line is identified by the SHA-256 of its text.
A line is listed when the two reports do not carry the same number of
branch records for that hash. Throw notes are kept in the count. Call
records are not branches.

Usage: python3 scripts/coverage_line_diff.py DIR_A DIR_B
DIR_A and DIR_B hold *.gcov files from gcov -b -c.
"""

import hashlib
import re
import sys
from pathlib import Path


def lines_of(path):
    rows = {}
    current = None
    text = None
    branches = 0
    for raw in path.read_text(errors="replace").splitlines():
        if re.search(r":\s*0:Source:", raw):
            if current is not None and text is not None:
                rows.setdefault(current, {"text": text, "branches": 0})
                rows[current]["branches"] += branches
            current = None
            text = None
            branches = 0
            continue
        matched = re.match(r"\s*([^:]+):\s*(\d+):(.*)", raw)
        if matched and matched.group(1).strip() not in ("-",):
            if current is not None and text is not None:
                rows.setdefault(current, {"text": text, "branches": 0})
                rows[current]["branches"] += branches
            number = int(matched.group(2))
            body = matched.group(3)
            if number == 0:
                current = None
                text = None
                branches = 0
                continue
            current = number
            text = body
            branches = 0
            continue
        if current is None:
            continue
        if re.match(r"branch\s+\d+\s+(taken\s+|never executed)", raw.strip()):
            branches += 1
    if current is not None and text is not None:
        rows.setdefault(current, {"text": text, "branches": 0})
        rows[current]["branches"] += branches
    return rows


def load_dir(directory):
    found = {}
    for path in sorted(Path(directory).glob("*.gcov")):
        name = path.name[: -len(".gcov")]
        found[name] = lines_of(path)
    return found


def main(argv):
    if len(argv) != 3:
        print("usage: coverage_line_diff.py DIR_A DIR_B", file=sys.stderr)
        return 2
    left = load_dir(argv[1])
    right = load_dir(argv[2])
    names = sorted(set(left) | set(right))
    diffs = 0
    for name in names:
        a = left.get(name, {})
        b = right.get(name, {})
        keys = set()
        for side in (a, b):
            for number, row in side.items():
                digest = hashlib.sha256(row["text"].encode()).hexdigest()
                keys.add((number, digest, row["text"]))
        for number, digest, text in sorted(keys):
            count_a = a.get(number, {}).get("branches")
            count_b = b.get(number, {}).get("branches")
            if a.get(number, {}).get("text") != text:
                count_a = None
            if b.get(number, {}).get("text") != text:
                count_b = None
            if count_a == count_b:
                continue
            diffs += 1
            print(f"{name}:{number}\t{count_a}\t{count_b}\t{digest[:16]}\t{text.strip()[:80]}")
    print(f"differing_lines {diffs}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
