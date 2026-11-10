#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run branch coverage on the in-situ tests and diff-cover against the base."""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlog


def main():
    root = Path(os.environ["SRC_SEAMS"])
    build = Path(os.environ["BUILD"])
    log_path = Path(os.environ["PY_LOG"])
    xml = Path(os.environ["PY_XML"])
    diff_log = Path(os.environ["DIFF_LOG"])
    base = os.environ.get("DIFF_BASE", "d763dd35")
    py = os.environ.get("DSEAMS_PYTHON", sys.executable)
    log = streamlog.Log(log_path)
    log.line("start: " + streamlog.date_is())
    log.line("HEAD_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True).strip())
    env = os.environ.copy()
    env["COVERAGE_FILE"] = str(build / "work" / "py.coverage")
    report = log.run([py, "-m", "coverage", "report", "--show-missing"], root, env)
    log.run([py, "-m", "coverage", "xml", "-o", str(xml)], root, env)
    log.line("xml: " + streamlog.clean(str(xml)))
    log.line("end: " + streamlog.date_is())
    log.close()
    diff = streamlog.Log(diff_log)
    diff.line("start: " + streamlog.date_is())
    diff.line("HEAD_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip())
    diff.line("TREE_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True).strip())
    diff.line("base: " + base)
    code = diff.run(
        [py, "-m", "diff_cover.diff_cover_tool", str(xml), "--compare-branch", base],
        root,
        env,
    )
    diff.line("diff_cover_exit: %d" % code)
    diff.line("end: " + streamlog.date_is())
    diff.close()
    if report != 0 or code != 0:
        raise SystemExit(code or report)


if __name__ == "__main__":
    main()
