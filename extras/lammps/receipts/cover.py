# Copyright (c) d-SEAMS developers
"""Build the adapter with coverage, run it, and record branch counts.

The log is streamed. gcov JSON is rewritten with the same prefix rules.
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlog

SRC = Path(os.environ["SRC_SEAMS"])
LAMMPS = Path(os.environ["SRC_LAMMPS"])
BUILD = Path(os.environ["BUILD"])
WORK = BUILD / "work" / "coverage"
LOG = BUILD / "kit" / "logs" / "coverage.log"


def require_live(code):
    if code != 0:
        raise SystemExit("test_live_cover.py exit %d" % code)


def configure_live(env, build):
    """Point a live run at the engines that are actually installed."""
    build = Path(build)
    for name, key in (
        ("seams-insitu", "SEAMS_INSITU_EXE"),
        ("frame_clock", "DSEAMS_FRAME_CLOCK"),
    ):
        exe = build / "seams" / name
        if exe.is_file():
            env[key] = str(exe)
    mpi_prefix = build / "lmp-mpi"
    (mpi_prefix / "bin").mkdir(parents=True, exist_ok=True)
    (mpi_prefix / "lib").mkdir(parents=True, exist_ok=True)
    mpirun_bin = Path(os.environ.get("MPIRUN", os.environ.get("DSEAMS_MPIRUN", "mpirun")))
    if mpirun_bin.is_file():
        (mpi_prefix / "bin" / "mpirun").unlink(missing_ok=True)
        (mpi_prefix / "bin" / "mpirun").symlink_to(mpirun_bin)
    intree = build / "lammps-intree" / "lmp"
    if intree.is_file():
        (mpi_prefix / "bin" / "lmp_mpi").unlink(missing_ok=True)
        (mpi_prefix / "bin" / "lmp_mpi").symlink_to(intree)
    for so in (build / "lammps-intree").glob("*.so*"):
        dest = mpi_prefix / "lib" / so.name
        if dest.is_symlink() or dest.exists():
            dest.unlink()
        dest.symlink_to(so)
    plugin_so = build / "plugin" / "libdseams_lammps.so"
    if plugin_so.is_file():
        env["DSEAMS_LMP_MPI_PREFIX"] = str(mpi_prefix)
        env["DSEAMS_LMP_MPI_PLUGIN"] = str(plugin_so)
    # Each lmp binary keeps its own liblammps on RUNPATH. Putting one
    # build directory here would make the other binary load the wrong one.
    env["LD_LIBRARY_PATH"] = ":".join(
        part
        for part in (
            str(build / "seams" / "src"),
            str(build / "seams" / "chemfiles-prefix" / "lib"),
            env.get("LD_LIBRARY_PATH", ""),
        )
        if part
    )
    return env


def clear_gcda(root):
    root = Path(root)
    if not root.is_dir():
        return
    for path in root.rglob("*.gcda"):
        path.unlink()


def abort(log, code):
    if code != 0:
        log.line("exit_code: %d" % code)
        log.line("end: " + streamlog.date_is())
        log.close()
        raise SystemExit(code)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    log = streamlog.Log(LOG)
    log.line("start: " + streamlog.date_is())
    head = subprocess.check_output(["git", "-C", str(SRC), "rev-parse", "HEAD"], text=True).strip()
    tree = subprocess.check_output(["git", "-C", str(SRC), "rev-parse", "HEAD^{tree}"], text=True).strip()
    lhead = subprocess.check_output(["git", "-C", str(LAMMPS), "rev-parse", "HEAD"], text=True).strip()
    ltree = subprocess.check_output(["git", "-C", str(LAMMPS), "rev-parse", "HEAD^{tree}"], text=True).strip()
    log.line("HEAD_SEAMS: " + head)
    log.line("TREE_SEAMS: " + tree)
    log.line("HEAD_LAMMPS: " + lhead)
    log.line("TREE_LAMMPS: " + ltree)
    log.line("command: python3 tools/cover.py")
    log.line("substitution: live, longest prefix first, see PATHS.txt")
    log.line("---")
    flags = " ".join([
        "--coverage", "-O0", "-fprofile-update=atomic",
    ])
    plugin = WORK / "plugin"
    code = log.run(
        [
            "cmake", "-S", str(SRC / "src" / "adapters" / "lammps"), "-B", str(plugin),
            "-DLAMMPS_SRC=" + str(LAMMPS / "src"),
            "-DSEAMS_INCLUDE=" + str(SRC / "src" / "include" / "internal"),
            "-DSEAMS_LIBRARY=" + str(BUILD / "seams" / "src" / "libyodaLib.so"),
            "-DDSEAMS_LAMMPS_MPI=ON",
            "-DCMAKE_BUILD_TYPE=Debug",
            "-DCMAKE_CXX_COMPILER=mpicxx",
            "-DCMAKE_CXX_FLAGS=" + flags,
            "-DCMAKE_EXE_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
            "-DCMAKE_SHARED_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
            "-DCMAKE_MODULE_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
        ],
        WORK,
        os.environ.copy(),
    )
    abort(log, code)
    code = log.run(["cmake", "--build", str(plugin), "-j", "4"], WORK, os.environ.copy())
    abort(log, code)
    # Counters from an earlier or a parallel run are not part of this figure.
    clear_gcda(plugin)
    clear_gcda(WORK / "compute")
    clear_gcda(BUILD / "seams" / "seams.p")
    os.environ["OMP_NUM_THREADS"] = "1"
    env = os.environ.copy()
    env["DSEAMS_LAMMPS_SRC"] = str(LAMMPS / "src")
    env["SEAMS_LIB"] = str(BUILD / "seams" / "src" / "libyodaLib.so")
    env["DSEAMS_MPICXX"] = "mpicxx"
    env["DSEAMS_LMP"] = str(BUILD / "lammps-stock" / "lmp")
    env["DSEAMS_LMP_INTREE"] = str(BUILD / "lammps-intree" / "lmp")
    env["DSEAMS_MPIRUN"] = os.environ["MPIRUN"]
    env["DSEAMS_ADAPTER_WORK"] = str(WORK)
    env["SEAMS_CLI"] = str(BUILD / "seams" / "seams")
    env["LD_LIBRARY_PATH"] = str(BUILD / "seams" / "src") + ":" + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    env["SRC_SEAMS"] = str(SRC)
    env["COVERAGE_FILE"] = str(BUILD / "work" / "py.coverage")
    py = os.environ.get("DSEAMS_PYTHON", sys.executable)
    code = log.run(
        [py, "-m", "coverage", "run", "--branch", str(SRC / "tests" / "insitu" / "adapter_cov.py")],
        WORK,
        env,
    )
    log.line("adapter_cov_exit: %d" % code)
    abort(log, code)
    env["DSEAMS_PLUGIN"] = str(BUILD / "plugin" / "libdseams_lammps.so")
    env["DSEAMS_IN_DSEAMS_WORK"] = str(WORK / "shipped")
    for args in (
        [str(SRC / "tests" / "insitu" / "test_shipped_in_dseams.py")],
        [str(SRC / "extras" / "lammps" / "receipts" / "import_checks.py")],
        [str(SRC / "extras" / "lammps" / "receipts" / "receipt_fail.py")],
        [str(SRC / "tests" / "cli_every.py"), str(BUILD / "seams" / "seams"), str(SRC)],
    ):
        run_env = env
        if Path(args[0]).name == "receipt_fail.py":
            # OVITO's Qt libraries need libGL. That directory breaks LAMMPS,
            # so only this process, which does not launch lmp, receives it.
            run_env = dict(env)
            glue = "/lib/x86_64-linux-gnu"
            current = run_env.get("LD_LIBRARY_PATH", "")
            if glue not in current.split(":"):
                run_env["LD_LIBRARY_PATH"] = (current + ":" + glue).strip(":")
        code = log.run([py, "-m", "coverage", "run", "--branch", "--append", *args], WORK, run_env)
        log.line("py %s exit %d" % (Path(args[0]).name, code))
        abort(log, code)
    live_env = dict(env)
    live_env["BUILD"] = str(BUILD)
    live_env["DSEAMS_LIVE_WORK"] = str(BUILD / "work" / "live-cover")
    live_env["SEAMS_CLI"] = str(BUILD / "seams" / "seams")
    live_env["SEAMS_LIB"] = str(BUILD / "seams" / "src" / "libyodaLib.so")
    live_env["SEAMS_LAMMPS_PLUGIN"] = str(BUILD / "plugin" / "libdseams_lammps.so")
    live_env["DSEAMS_MD_PREFIX"] = str(BUILD / "md-prefix")
    live_env["DSEAMS_LMP"] = str(BUILD / "lammps-stock" / "lmp")
    live_env["DSEAMS_LMP_INTREE"] = str(BUILD / "lammps-intree" / "lmp")
    live_env["DSEAMS_PLUGIN"] = str(BUILD / "plugin" / "libdseams_lammps.so")
    live_env["DSEAMS_MPIRUN"] = os.environ.get("MPIRUN", os.environ.get("DSEAMS_MPIRUN", ""))
    configure_live(live_env, BUILD)
    code = log.run(
        [py, "-m", "coverage", "run", "--branch", "--append", str(SRC / "tests" / "insitu" / "test_live_cover.py")],
        WORK,
        live_env,
    )
    log.line("py test_live_cover.py exit %d" % code)
    log.line("live_suite_exit: %d" % code)
    require_live(code)
    log.line("gcov_object_directory: $BUILD/work/coverage")
    fix_notes = next(plugin.rglob("fix_dseams.cpp.gcno"))
    comp_notes = next((WORK / "compute").rglob("compute_dseams.gcno"))
    cli_notes = next(Path(os.environ["BUILD"]).joinpath("seams").rglob("*seams_cli.cpp.gcno"))
    # Meson records seams_cli.cpp as ../../src/seams/... relative to the
    # build root. gcov resolves that path from its working directory.
    gcov_jobs = (
        (fix_notes, WORK),
        (comp_notes, WORK),
        (cli_notes, BUILD / "seams"),
    )
    gcov_code = 0
    for note, cwd in gcov_jobs:
        log.line("command: gcov -b " + str(note))
        log.line("gcov_cwd: " + str(cwd))
        one = log.run(["gcov", "-b", str(note)], cwd, os.environ.copy())
        log.line("gcov exit %d" % one)
        gcov_code = one or gcov_code
    gcovr = os.environ.get("GCOVR", "gcovr")
    targets = (
        ("fix", WORK / "plugin" / "CMakeFiles" / "dseams_lammps.dir", ".*fix_dseams\\.cpp"),
        ("compute", WORK / "compute", ".*compute_dseams\\.cpp"),
        ("cli", BUILD / "seams" / "seams.p", ".*seams_cli\\.cpp"),
    )
    flags = ["--exclude-throw-branches", "--exclude-unreachable-branches"]
    gcovr_code = 0
    for label, folder, filt in targets:
        for name, extra in (("nomarkers", ["--no-markers"]), ("markers", [])):
            cmd = [
                gcovr, "-j", "1", "--root", str(SRC),
                "--object-directory", str(folder),
                "--filter", filt,
                "--txt-metric", "line", "--txt-metric", "branch",
                *flags, *extra,
                "--txt", str(WORK / ("gcovr-%s-%s.txt" % (label, name))),
                "--json", str(WORK / ("gcovr-%s-%s.json" % (label, name))),
            ]
            shown = "gcovr -j 1 --exclude-throw-branches --exclude-unreachable-branches"
            if extra:
                shown += " --no-markers"
            shown += " " + label
            log.line("command: " + shown)
            gcovr_code = log.run(cmd, WORK, os.environ.copy()) or gcovr_code
    log.line("end: " + streamlog.date_is())
    log.hash_tree(WORK)
    log.line("exit_code: %d" % (gcov_code or gcovr_code or code))
    log.close()
    abort(log, gcov_code or gcovr_code)


if __name__ == "__main__":
    main()
