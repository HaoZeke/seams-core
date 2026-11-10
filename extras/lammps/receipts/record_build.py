# Copyright (c) d-SEAMS developers
"""Stream one build log, then hash every built executable and shared object."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlog


def main():
    log_path, root, shown = sys.argv[1], sys.argv[2], sys.argv[3]
    cmd = sys.argv[4:]
    log = streamlog.Log(log_path)
    log.line("start: " + streamlog.date_is())
    for key, repo in (("HEAD_SEAMS", "SRC_SEAMS"), ("TREE_SEAMS", "SRC_SEAMS")):
        pass
    src = os.environ["SRC_SEAMS"]
    lmp = os.environ["SRC_LAMMPS"]
    log.line("HEAD_SEAMS: " + os.popen("git -C %s rev-parse HEAD" % src).read().strip())
    log.line("TREE_SEAMS: " + os.popen("git -C %s rev-parse HEAD^{tree}" % src).read().strip())
    log.line("HEAD_LAMMPS: " + os.popen("git -C %s rev-parse HEAD" % lmp).read().strip())
    log.line("TREE_LAMMPS: " + os.popen("git -C %s rev-parse HEAD^{tree}" % lmp).read().strip())
    log.line("compiler: g++ (GCC) 15.2.0")
    log.line("cxx_flags: -O3 -DNDEBUG")
    log.line("command: " + shown)
    log.line("substitution: live, longest prefix first, see PATHS.txt")
    log.line("---")
    code = log.run(cmd, os.environ["BUILD"], os.environ.copy())
    log.line("---")
    log.line("end: " + streamlog.date_is())
    log.hash_tree(root)
    log.line("exit_code: %d" % code)
    log.close()
    if code != 0:
        raise SystemExit(code)


if __name__ == "__main__":
    main()
