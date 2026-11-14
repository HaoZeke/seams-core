#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the six family examples out of tree against an installed libseams.

Usage: check_installed.py BUILD SOURCE

Installs BUILD into BUILD/capi-family-destdir, then builds every example in
docs/examples/c twice, once with the flags `pkg-config seams` reports and
once with CMake's find_package(seams), and runs each binary against the
installed library. CC names the C compiler.
"""

import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

FAMILIES = ("neighbours", "rings", "chill_plus", "q4q6", "cages", "errors")
WARN = ["-std=c11", "-Wall", "-Wextra", "-Wpedantic", "-Werror"]


def helpers(source):
    path = source / "examples" / "c" / "check_installed.py"
    spec = importlib.util.spec_from_file_location("quickstart_installed", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_examples(label, binaries, env, source, stage, help_mod):
    for family, binary in binaries:
        args = [str(binary)] if family == "errors" else [str(binary), str(source)]
        help_mod.check_resolved(binary, env, stage)
        listed = subprocess.check_output(["ldd", str(binary)], env=env, text=True)
        resolved = [line for line in listed.splitlines() if line.strip().startswith("libseams.so")]
        if len(resolved) != 1:
            raise SystemExit("{} {} does not load libseams:\n{}".format(label, family, listed))
        done = subprocess.run(args, env=env, capture_output=True, text=True)
        if done.returncode != 0:
            raise SystemExit("{} {} exit {}\n{}{}".format(label, family, done.returncode, done.stdout, done.stderr))
        print("{} {}: {}".format(label, family, done.stdout.strip()))


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    build = Path(argv[1]).resolve()
    source = Path(argv[2]).resolve()
    help_mod = helpers(source)
    stage = build / "capi-family-destdir"
    work = build / "capi-family-work"
    for path in (stage, work):
        shutil.rmtree(path, ignore_errors=True)
    work.mkdir(parents=True)
    subprocess.check_call(
        ["meson", "install", "-C", str(build), "--no-rebuild", "--quiet", "--destdir", str(stage)]
    )

    pc_files = list(stage.rglob("seams.pc"))
    configs = list(stage.rglob("seams-config.cmake"))
    if len(pc_files) != 1 or len(configs) != 1:
        raise SystemExit("install did not produce one seams.pc and seams-config.cmake")
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(pc_files[0].parent)
    env["PKG_CONFIG_SYSROOT_DIR"] = str(stage)
    run_env = help_mod.run_env(env, help_mod.staged_libdir(stage))
    cc = shlex.split(os.environ.get("CC", "cc"))
    extra = help_mod.consumer_flags()

    pkg_config = help_mod.pkg_config_program()
    flags = subprocess.check_output([pkg_config, "--cflags", "--libs", "seams"], env=env, text=True).split()
    examples = source / "docs" / "examples" / "c"
    binaries = []
    for family in FAMILIES:
        binary = work / ("pkg-" + family)
        subprocess.check_call([*cc, *WARN, *extra, str(examples / (family + ".c")), "-o", str(binary), *flags], env=env)
        binaries.append((family, binary))
    run_examples("pkg-config", binaries, run_env, source, stage, help_mod)

    installed = stage / Path(help_mod.prefix_of(build)).relative_to("/")
    cmake_build = work / "cmake"
    subprocess.check_call(
        [
            "cmake", "-S", str(examples), "-B", str(cmake_build),
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DCMAKE_BUILD_TYPE=Release",
            *help_mod.cmake_consumer_flags(cc, extra),
        ],
        env=env,
    )
    subprocess.check_call(["cmake", "--build", str(cmake_build)], env=env)
    run_examples("cmake", [(family, cmake_build / family) for family in FAMILIES], run_env, source, stage, help_mod)
    print("installed family examples ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
