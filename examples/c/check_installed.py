#!/usr/bin/env python3
"""Install libseams and build examples/c against that install."""

import atexit
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def pkg_config_program():
    """The conda pkg-config wrapper applies PKG_CONFIG_SYSROOT_DIR twice."""
    found = shutil.which("pkg-config")
    if found:
        sibling = Path(found).with_name("pkg-config.bin")
        if sibling.is_file():
            return str(sibling)
    return found or "pkg-config"


def run(cmd, env=None, cwd=None):
    subprocess.check_call(cmd, env=env, cwd=cwd)


def prefix_of(build):
    import json

    options = json.loads(
        subprocess.check_output(
            ["meson", "introspect", "--buildoptions", str(build)], text=True
        )
    )
    for option in options:
        if option["name"] == "prefix":
            return option["value"]
    raise SystemExit("meson prefix is missing")


def project_version(build):
    import json

    info = json.loads(
        subprocess.check_output(["meson", "introspect", "--projectinfo", str(build)], text=True)
    )
    return info["version"]


def consumer_flags():
    """Sanitizer flags of an instrumented build, which its consumers share."""
    return shlex.split(os.environ.get("SEAMS_CONSUMER_FLAGS", ""))


def cmake_consumer_flags(cc, extra):
    flags = ["-DCMAKE_C_COMPILER=" + cc[0]]
    if extra:
        joined = " ".join(extra)
        flags += ["-DCMAKE_C_FLAGS=" + joined, "-DCMAKE_EXE_LINKER_FLAGS=" + joined]
    return flags


def scratch(build, prefix):
    """A fresh directory under BUILD that is removed at exit. A --repeat copy
    of the same test runs beside this one and must not share its install."""
    parent = build / "install-stages"
    parent.mkdir(exist_ok=True)
    path = Path(tempfile.mkdtemp(prefix=prefix + "-", dir=parent))
    atexit.register(shutil.rmtree, path, True)
    return path


def staged_libdir(stage):
    libs = list(stage.rglob("libseams.so"))
    if len(libs) != 1:
        raise SystemExit("install did not produce one libseams.so")
    return libs[0].parent


def run_env(env, libdir):
    """Only the staged libdir. The installed libraries find each other
    through their $ORIGIN rpath, so nothing in the build tree is needed."""
    found = dict(env)
    found["LD_LIBRARY_PATH"] = str(libdir)
    return found


def check_resolved(binary, env, stage):
    """Every library the install ships resolves inside the stage."""
    listed = subprocess.check_output(["ldd", str(binary)], text=True, env=env)
    for line in listed.splitlines():
        if "not found" in line:
            raise SystemExit("{} has an unresolved library:\n{}".format(binary, listed))
    shipped = {path.name for path in stage.rglob("*.so*") if path.is_file()}
    for line in listed.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[1] == "=>" and parts[0] in shipped:
            if not parts[2].startswith(str(stage)):
                raise SystemExit("{} loads {} from outside the stage:\n{}".format(binary, parts[0], listed))
    return env


def main():
    build = Path(sys.argv[1]).resolve()
    source = Path(sys.argv[2]).resolve()
    base = scratch(build, "capi")
    stage = base / "destdir"
    work = base / "work"
    work.mkdir()
    run(["meson", "install", "-C", str(build), "--destdir", str(stage)])

    pc_files = list(stage.rglob("seams.pc"))
    configs = list(stage.rglob("seams-config.cmake"))
    if len(pc_files) != 1 or len(configs) != 1:
        raise SystemExit("install did not produce seams.pc and seams-config.cmake")
    pc = pc_files[0]
    libdir = staged_libdir(stage)
    compile_env = os.environ.copy()
    compile_env["PKG_CONFIG_PATH"] = str(pc.parent)
    compile_env["PKG_CONFIG_SYSROOT_DIR"] = str(stage)
    pkg_config = pkg_config_program()
    cflags = subprocess.check_output(
        [pkg_config, "--cflags", "seams"], env=compile_env, text=True
    ).split()
    libs = subprocess.check_output(
        [pkg_config, "--libs", "seams"], env=compile_env, text=True
    ).split()
    binary = work / "quickstart-pkg"
    cc = shlex.split(os.environ.get("CC", "cc"))
    extra = consumer_flags()
    run(
        [*cc, *extra, str(source / "examples" / "c" / "quickstart.c"), "-o", str(binary), *cflags, *libs, "-lm"],
        env=compile_env,
    )
    env = run_env(compile_env, libdir)
    check_resolved(binary, env, stage)
    sample = source / "input" / "traj" / "agi_water_tiny.lammpstrj"
    run([str(binary), str(sample)], env=env)

    version = project_version(build)
    modversion = subprocess.check_output(
        [pkg_config, "--modversion", "seams"], env=env, text=True
    ).strip()
    if modversion != version:
        raise SystemExit("pkg-config says {}, project is {}".format(modversion, version))
    pinned = subprocess.run(
        [pkg_config, "--exists", "seams >= " + version], env=env
    )
    if pinned.returncode != 0:
        raise SystemExit("pkg-config --exists 'seams >= {}' failed".format(version))
    find_check = source / "examples" / "c" / "find_check"
    version_bin = work / "version-pkg"
    run(
        [*cc, *extra, str(find_check / "version.c"), "-o", str(version_bin), *cflags, *libs],
        env=compile_env,
    )
    check_resolved(version_bin, env, stage)
    run([str(version_bin), version], env=env)

    prefix = prefix_of(build)
    installed = stage / Path(prefix).relative_to("/")
    cmake_flags = cmake_consumer_flags(cc, extra)
    cmake_build = work / "cmake"
    run(
        [
            "cmake",
            "-S",
            str(source / "examples" / "c"),
            "-B",
            str(cmake_build),
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DCMAKE_BUILD_TYPE=Release",
            *cmake_flags,
        ],
        env=compile_env,
    )
    run(["cmake", "--build", str(cmake_build)], env=compile_env)
    cmake_bin = cmake_build / "quickstart"
    check_resolved(cmake_bin, env, stage)
    run([str(cmake_bin), str(sample)], env=env)

    find_build = work / "cmake-find"
    run(
        [
            "cmake",
            "-S",
            str(find_check),
            "-B",
            str(find_build),
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DSEAMS_EXPECT_VERSION=" + version,
            *cmake_flags,
        ],
        env=compile_env,
    )
    run(["cmake", "--build", str(find_build)], env=compile_env)
    find_bin = find_build / "version"
    check_resolved(find_bin, env, stage)
    run([str(find_bin), version], env=env)
    print("installed C quickstart ok")


if __name__ == "__main__":
    main()
