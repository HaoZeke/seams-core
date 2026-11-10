# Copyright (c) d-SEAMS developers
"""Write a command log one line at a time.

Timestamps come from date -Is. Each output line is rewritten and appended
before the next line is read. Nothing is buffered and written later.
"""

import hashlib
import os
import subprocess
from pathlib import Path


def prefixes():
    pairs = []
    for key, token in (
        ("SRC_SEAMS", "$SRC_SEAMS"),
        ("SRC_LAMMPS", "$SRC_LAMMPS"),
        ("BUILD", "$BUILD"),
        ("MPIRUN", "$MPIRUN"),
        ("HOME", "$HOME"),
        ("CVMFS_PREFIX", "$CVMFS"),
    ):
        value = os.environ.get(key, "")
        if value:
            pairs.append((value, token))
    if "CVMFS_PREFIX" not in os.environ:
        pairs.append(("/" + "cvmfs", "$CVMFS"))
    pairs.sort(key=lambda item: len(item[0]), reverse=True)
    return pairs


def clean(text, pairs=None):
    pairs = prefixes() if pairs is None else pairs
    for src, token in pairs:
        text = text.replace(src, token)
    return text


def date_is():
    return subprocess.check_output(["date", "-Is"], text=True).strip()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def elf_kind(path):
    try:
        with open(path, "rb") as handle:
            head = handle.read(18)
    except OSError:
        return None
    if len(head) < 18 or head[:4] != b"\x7fELF":
        return None
    kind = int.from_bytes(head[16:18], "little" if head[5] == 1 else "big")
    if kind == 2:
        return "exec"
    if kind == 3:
        return "dyn"
    return None


class Log:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("w")
        self.pairs = prefixes()

    def line(self, text):
        self.handle.write(clean(text, self.pairs) + "\n")
        self.handle.flush()

    def close(self):
        self.handle.close()

    def run(self, cmd, cwd, env):
        wrapped = ["stdbuf", "-oL", "-eL", *cmd]
        proc = subprocess.Popen(
            wrapped,
            cwd=cwd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for raw in proc.stdout:
            self.line(raw.rstrip("\n"))
        return proc.wait()

    def hash_tree(self, root):
        root = Path(root)
        paths = []
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix in (".o", ".a", ".gcno", ".gcda"):
                continue
            if elf_kind(path) is None:
                continue
            paths.append(path)
        for path in sorted(paths):
            self.line("binary_sha256: %s %s" % (path, sha256(path)))
