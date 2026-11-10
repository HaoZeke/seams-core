#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run branch coverage on the in-situ tests and diff-cover against the base."""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlog


def write_heads(log, root, lammps):
    log.line("HEAD_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True).strip())
    log.line("HEAD_LAMMPS: " + subprocess.check_output(["git", "-C", str(lammps), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_LAMMPS: " + subprocess.check_output(["git", "-C", str(lammps), "rev-parse", "HEAD^{tree}"], text=True).strip())


def gcovr_diff(root, lammps, py, base, env):
    xml = os.environ.get("GCOVR_XML", "")
    dest = os.environ.get("GCOVR_DIFF_LOG", "")
    if not xml or not dest:
        return 0
    diff = streamlog.Log(dest)
    diff.line("start: " + streamlog.date_is())
    write_heads(diff, root, lammps)
    diff.line("base: " + base)
    code = diff.run([py, "-m", "diff_cover.diff_cover_tool", xml, "--compare-branch", base], root, env)
    diff.line("gcovr_exit: %d" % code)
    diff.line("end: " + streamlog.date_is())
    diff.close()
    return code


def raise_status(status):
    if status:
        raise SystemExit(status)


def main():
    root = Path(os.environ["SRC_SEAMS"])
    build = Path(os.environ["BUILD"])
    log_path = Path(os.environ["PY_LOG"])
    xml = Path(os.environ["PY_XML"])
    diff_log = Path(os.environ["DIFF_LOG"])
    base = os.environ.get("DIFF_BASE", "d763dd35")
    py = os.environ.get("DSEAMS_PYTHON", sys.executable)
    lammps = Path(os.environ["SRC_LAMMPS"])
    log = streamlog.Log(log_path)
    log.line("start: " + streamlog.date_is())
    write_heads(log, root, lammps)
    env = os.environ.copy()
    env["COVERAGE_FILE"] = str(build / "work" / "py.coverage")
    report = log.run([py, "-m", "coverage", "report", "--show-missing"], root, env)
    log.run([py, "-m", "coverage", "xml", "-o", str(xml)], root, env)
    log.line("xml: " + streamlog.clean(str(xml)))
    log.line("end: " + streamlog.date_is())
    log.close()
    diff = streamlog.Log(diff_log)
    diff.line("start: " + streamlog.date_is())
    write_heads(diff, root, lammps)
    diff.line("base: " + base)
    code = diff.run(
        [py, "-m", "diff_cover.diff_cover_tool", str(xml), "--compare-branch", base],
        root,
        env,
    )
    diff.line("diff_cover_exit: %d" % code)
    diff.line("end: " + streamlog.date_is())
    diff.close()
    gcode = gcovr_diff(root, lammps, py, base, env)
    raise_status(code or gcode or report)


if __name__ == "__main__":
    main()
