# Copyright (c) d-SEAMS developers
"""Run a command and print each liblammps or plugin path it loaded.

The loader log is the source. A maps poll misses a run that exits first.
"""

import hashlib
import os
import subprocess
import sys
import threading
import time


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parent_of(text):
    mark = text.rfind(")")
    if mark < 0:
        return None
    fields = text[mark + 2 :].split()
    if len(fields) < 2:
        return None
    return int(fields[1])


def ppid(pid):
    try:
        text = open("/proc/%d/stat" % pid, encoding="utf-8", errors="replace").read()
    except OSError:
        return None
    return parent_of(text)


def map_path(line):
    if "liblammps.so" not in line and "libdseams_lammps.so" not in line:
        return ""
    parts = line.split()
    if len(parts) >= 6:
        return parts[-1]
    return ""


def descendants(root):
    found = {root}
    changed = True
    while changed:
        changed = False
        names = os.listdir("/proc")
        for name in names:
            if not name.isdigit():
                continue
            pid = int(name)
            if pid in found:
                continue
            parent = ppid(pid)
            if parent in found:
                found.add(pid)
                changed = True
    return found


def mapped(pid):
    try:
        lines = open("/proc/%d/maps" % pid, encoding="utf-8", errors="replace").read().splitlines()
    except OSError:
        return []
    paths = []
    for line in lines:
        path = map_path(line)
        if path:
            paths.append(path)
    return paths


def paths_from_loader(line):
    found = []
    for marker in ("trying file=", "calling init: ", "file="):
        at = line.find(marker)
        if at < 0:
            continue
        rest = line[at + len(marker):].strip()
        if not rest:
            continue
        path = rest.split()[0]
        if "libdseams_lammps.so" in path or "liblammps.so" in path:
            found.append(path)
    return found


def main():
    env = os.environ.copy()
    env["LD_DEBUG"] = "libs"
    proc = subprocess.Popen(sys.argv[1:], env=env, stderr=subprocess.PIPE)
    seen = set()
    stop = False

    def take(path):
        if not os.path.isabs(path):
            path = os.path.abspath(path)
        if os.path.exists(path):
            seen.add(path)

    def watch():
        while not stop:
            for pid in descendants(proc.pid):
                for path in mapped(pid):
                    take(path)
            time.sleep(0.01)

    thread = threading.Thread(target=watch)
    thread.start()
    assert proc.stderr is not None
    for raw in proc.stderr:
        for path in paths_from_loader(raw.decode("utf-8", "replace")):
            take(path)
    code = proc.wait()
    stop = True
    thread.join()
    for path in sorted(seen):
        print("loaded %s %s" % (path, sha256(path)), flush=True)
    if not seen:
        print("loaded none", flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
