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


def prove_live_report():
    try:
        require_live(1)
    except SystemExit as err:
        if err.code != "test_live_cover.py exit 1":
            raise
    else:
        raise SystemExit("a failing live test was accepted")
    require_live(0)


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
        "--coverage", "-O0", "-fprofile-update=single",
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
        ],
        WORK,
        os.environ.copy(),
    )
    abort(log, code)
    code = log.run(["cmake", "--build", str(plugin), "-j", "4"], WORK, os.environ.copy())
    abort(log, code)
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
        code = log.run([py, "-m", "coverage", "run", "--branch", "--append", *args], WORK, env)
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
    gmx = BUILD / "gromacs-dseams" / "bin" / "gmx_mpi"
    md_bin = BUILD / "md-prefix" / "bin"
    md_bin.mkdir(parents=True, exist_ok=True)
    if gmx.is_file():
        live_env["DSEAMS_GMX"] = str(gmx)
        live_env["DSEAMS_GMX_MPIRUN"] = os.environ.get("MPIRUN", os.environ.get("DSEAMS_MPIRUN", ""))
        data = os.environ.get("GMXDATA", "")
        if data:
            live_env["GMXDATA"] = data
        link = md_bin / "gmx"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(gmx)
        top = Path(data) / "top" if data else Path()
        if top.is_dir():
            share = BUILD / "md-prefix" / "share" / "gromacs"
            share.mkdir(parents=True, exist_ok=True)
            dest = share / "top"
            if dest.is_symlink() or dest.exists():
                dest.unlink()
            dest.symlink_to(top)
    insitu = BUILD / "seams" / "seams-insitu"
    imd = BUILD / "seams" / "seams-gromacs-imd"
    clock = BUILD / "seams" / "frame_clock"
    if insitu.is_file():
        live_env["SEAMS_INSITU_EXE"] = str(insitu)
    if imd.is_file():
        live_env["SEAMS_GRO_IMD_EXE"] = str(imd)
    if clock.is_file():
        live_env["DSEAMS_FRAME_CLOCK"] = str(clock)
    mpi_prefix = BUILD / "lmp-mpi"
    (mpi_prefix / "bin").mkdir(parents=True, exist_ok=True)
    (mpi_prefix / "lib").mkdir(parents=True, exist_ok=True)
    mpirun_bin = Path(os.environ.get("MPIRUN", os.environ.get("DSEAMS_MPIRUN", "mpirun")))
    if mpirun_bin.is_file():
        (mpi_prefix / "bin" / "mpirun").unlink(missing_ok=True)
        (mpi_prefix / "bin" / "mpirun").symlink_to(mpirun_bin)
    intree = BUILD / "lammps-intree" / "lmp"
    if intree.is_file():
        (mpi_prefix / "bin" / "lmp_mpi").unlink(missing_ok=True)
        (mpi_prefix / "bin" / "lmp_mpi").symlink_to(intree)
    for so in (BUILD / "lammps-intree").glob("*.so*"):
        dest = mpi_prefix / "lib" / so.name
        if dest.is_symlink() or dest.exists():
            dest.unlink()
        dest.symlink_to(so)
    plugin_so = BUILD / "plugin" / "libdseams_lammps.so"
    if plugin_so.is_file():
        live_env["DSEAMS_LMP_MPI_PREFIX"] = str(mpi_prefix)
        live_env["DSEAMS_LMP_MPI_PLUGIN"] = str(plugin_so)
    live_env["LD_LIBRARY_PATH"] = ":".join(
        part
        for part in (
            str(BUILD / "gromacs-dseams" / "lib"),
            str(BUILD / "lammps-intree"),
            str(BUILD / "seams" / "src"),
            str(BUILD / "seams" / "chemfiles-prefix" / "lib"),
            live_env.get("LD_LIBRARY_PATH", ""),
        )
        if part
    )
    prove_live_report()
    code = log.run(
        [py, "-m", "coverage", "run", "--branch", "--append", str(SRC / "tests" / "insitu" / "test_live_cover.py")],
        WORK,
        live_env,
    )
    log.line("py test_live_cover.py exit %d" % code)
    require_live(code)
    log.line("gcov_object_directory: $BUILD/work/coverage")
    fix_notes = next(plugin.rglob("fix_dseams.cpp.gcno"))
    comp_notes = next((WORK / "compute").rglob("compute_dseams.gcno"))
    log.line("command: gcov -b " + str(fix_notes))
    gcov_code = log.run(["gcov", "-b", str(fix_notes)], WORK, os.environ.copy())
    log.line("command: gcov -b " + str(comp_notes))
    gcov_code = log.run(["gcov", "-b", str(comp_notes)], WORK, os.environ.copy()) or gcov_code
    cli_notes = next(Path(os.environ["BUILD"]).joinpath("seams").rglob("seams_cli.cpp.gcno"))
    log.line("command: gcov -b " + str(cli_notes))
    gcov_code = log.run(["gcov", "-b", str(cli_notes)], WORK, os.environ.copy()) or gcov_code
    gcovr = os.environ.get("GCOVR", "gcovr")
    common = [
        gcovr, "--root", str(SRC), "--object-directory", str(WORK),
        "--filter", ".*fix_dseams\\.cpp", "--filter", ".*compute_dseams\\.cpp",
        "--filter", ".*seams_cli\\.cpp",
        "--txt-metric", "line",
        "--txt-metric", "branch",
    ]
    line_patterns = [
        r"\s*if \(i < 0\)",
        r"\s*if \(i >= atom->nlocal\)",
        r"\s*if \(j < 0\)",
        r"\s*if \(j >= nall\)",
        r"\s*if \(want\)",
        r"\s*if \(static_cast<int>\(flat.size\(\)\) > begin\)",
        r"\s*\} else if \(total_atoms == 0\)",
        r"\s*if \(comm->me != 0\) \{ // nprocs is already 1",
        r"\s*continue; // builder stores i at or above zero",
        r"\s*continue; // builder stores i below nlocal",
        r"\s*continue; // NEIGHMASK keeps j at or above zero",
        r"\s*continue; // stencil stores j below nall",
        r"\s*return; // the only rank left is rank 0",
    ]
    raw_json = common + ["--json", str(WORK / "gcovr-raw.json"), "--txt", str(WORK / "gcovr-raw.txt")]
    log.line("command: gcovr raw branches")
    gcovr_code = log.run(raw_json, WORK, os.environ.copy())
    throw_json = common + [
        "--exclude-throw-branches",
        "--json", str(WORK / "gcovr-throw.json"),
        "--txt", str(WORK / "gcovr-throw.txt"),
    ]
    log.line("command: gcovr --exclude-throw-branches")
    gcovr_code = log.run(throw_json, WORK, os.environ.copy()) or gcovr_code
    gcovr_cmd = common + ["--exclude-throw-branches"]
    for pattern in line_patterns:
        gcovr_cmd.extend(["--exclude-lines-by-pattern", pattern])
    gcovr_cmd.extend([
        "--json", str(WORK / "gcovr.json"),
        "--xml", str(WORK / "gcovr.xml"),
        "--txt", str(WORK / "gcovr.txt"),
    ])
    log.line("command: gcovr --exclude-throw-branches plus one pattern per unreachable check")
    gcovr_code = log.run(gcovr_cmd, WORK, os.environ.copy()) or gcovr_code
    log.line("end: " + streamlog.date_is())
    log.hash_tree(WORK)
    log.line("exit_code: %d" % (gcov_code or gcovr_code or code))
    log.close()
    abort(log, gcov_code or gcovr_code)


if __name__ == "__main__":
    main()
