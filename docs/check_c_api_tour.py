#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Compile c_api_tour.c, run it, and require docs/api/c.md to quote every line."""

import importlib.util
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def helpers(root):
    path = root / "examples" / "c" / "check_installed.py"
    spec = importlib.util.spec_from_file_location("quickstart_installed", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    root = Path(sys.argv[1]).resolve()
    build = Path(sys.argv[2]).resolve()
    help_mod = helpers(root)
    page = (root / "docs" / "api" / "c.md").read_text()
    source = root / "docs" / "examples" / "c_api_tour.c"
    header = (root / "include" / "seams.h").read_text()
    functions = re.findall(r"SEAMS_API\s+(?:const\s+char\s+\*|[\w\s\*]+?)\s*(seams_[A-Za-z0-9_]+)\s*\(", header)
    # The tour links an installed copy, so it reads the same headers and
    # library a consumer gets.
    stage = build / "c-api-tour-destdir"
    shutil.rmtree(stage, ignore_errors=True)
    subprocess.check_call(
        ["meson", "install", "-C", str(build), "--no-rebuild", "--quiet", "--destdir", str(stage)]
    )
    libdir = help_mod.staged_libdir(stage)
    headers = list(stage.rglob("seams_version.h"))
    if len(headers) != 1 or not (headers[0].parent / "seams.h").is_file():
        fail("install did not produce seams.h and seams_version.h together")
    binary = build / "c-api-tour-work" / "c_api_tour"
    binary.parent.mkdir(parents=True, exist_ok=True)
    compile_cmd = [
        *shlex.split(os.environ.get("CC", "cc")),
        *help_mod.consumer_flags(),
        "-I" + str(headers[0].parent),
        str(source),
        "-L" + str(libdir),
        "-lseams",
        "-lm",
        "-o",
        str(binary),
    ]
    done = subprocess.run(compile_cmd, text=True, capture_output=True)
    if done.returncode != 0:
        fail(done.stderr[-800:])
    env = help_mod.run_env(os.environ.copy(), libdir)
    help_mod.check_resolved(binary, env, stage)
    env["OMP_NUM_THREADS"] = "1"
    run = subprocess.run([str(binary), str(root)], text=True, capture_output=True, env=env)
    if run.returncode != 0:
        fail(run.stderr[-800:])
    lines = [row for row in run.stdout.splitlines() if row.startswith("capi ")]
    if len(lines) < 60:
        fail("tour printed {} lines".format(len(lines)))
    missing = [row for row in lines if row not in page]
    if missing:
        fail("docs/api/c.md is missing:\n" + "\n".join(missing))
    quiet = [name for name in functions if "capi " + name not in run.stdout]
    if quiet:
        fail("tour did not print " + " ".join(quiet))
    print("c api tour {} lines".format(len(lines)))


if __name__ == "__main__":
    main()
