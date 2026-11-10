#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Instruction-count ratio of the public steinhardtQl path.

The driver calls chill::steinhardtQl on a chain of N atoms and on a
chain of 2N atoms. valgrind --tool=callgrind counts the instructions.
A linear pass stays near 2. A scan of every stored bond for every atom
stays near 4. The bound is 2.5. Wall time is not the assertion.

Callgrind 3.22 and 3.25 do not decode AVX-512. A sphericart library
built with -march=native dies at the first spherical-harmonic kernel.
Set SPHERICART_LIB to a build configured with SPHERICART_ARCH_NATIVE=OFF.
The steinhardt translation unit is the one in this tree.

usage: bop_insn_ratio.py PARENT_SHA BUILD RECEIPT_DIR
"""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

N = 25000
BOUND = 2.5


def run(cmd, cwd, env=None):
    return subprocess.run(
        cmd, cwd=cwd, text=True, capture_output=True, env=env
    )


def scrub(text):
    kept = []
    for line in text.splitlines():
        if any(
            p in line
            for p in ("/agent", "/tmp", "/cursor", "/cvmfs", "/home/")
        ):
            continue
        kept.append(line)
    return "\n".join(kept) + "\n"


def summary_of(path):
    text = Path(path).read_text(errors="replace")
    match = re.search(r"^summary:\s+(\d+)", text, re.M)
    if not match:
        raise SystemExit("no callgrind summary in %s" % path.name)
    return int(match.group(1))


def compile_entry(build):
    path = build / "compile_commands.json"
    data = json.loads(path.read_text())
    for entry in data:
        if entry.get("file", "").endswith("steinhardt_parallel.cpp"):
            return entry
    raise SystemExit("steinhardt_parallel.cpp is not in the compilation database")


def compile_unit(entry, source, output):
    argv = shlex.split(entry["command"])
    cleaned = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "-MD":
            i += 1
            continue
        if arg in ("-MF", "-MQ") and i + 1 < len(argv):
            i += 2
            continue
        if arg == "-o" and i + 1 < len(argv):
            cleaned.extend(["-o", str(output)])
            i += 2
            continue
        if arg == "-c" and i + 1 < len(argv):
            cleaned.extend(["-c", str(source)])
            i += 2
            continue
        cleaned.append(arg)
        i += 1
    ran = run(cleaned, entry["directory"])
    if ran.returncode != 0:
        sys.stderr.write(ran.stderr[-4000:])
        raise SystemExit(ran.returncode)
    return cleaned[0]


def link_library(compiler, obj, so_path, pieces):
    lib64 = str(Path(compiler).resolve().parent.parent / "lib64")
    cmd = [
        compiler, "-shared", "-fPIC", "-fopenmp", "-o", str(so_path),
        str(obj),
        str(pieces["ozaki"]),
        str(pieces["sphericart_ylm"]),
        "-L", pieces["sphericart_lib"], "-lsphericart",
        "-L", pieces["minimage_lib"], "-lminimage",
        "-lgcov", "-lgomp",
        "-Wl,-rpath," + pieces["sphericart_lib"],
        "-Wl,-rpath," + pieces["minimage_lib"],
        "-Wl,-rpath," + lib64,
        "-Wl,-soname,libyodaLib.so",
    ]
    ran = run(cmd, pieces["root"])
    if ran.returncode != 0:
        sys.stderr.write(ran.stderr[-4000:])
        raise SystemExit(ran.returncode)
    return lib64


def compile_driver(compiler, pieces, so_dir, exe, lib64):
    root = pieces["root"]
    cmd = [
        compiler, "-std=c++20", "-O0", "-g", "-fPIC",
        "-DSEAMS_HAS_OPENMP=1", "-fopenmp",
        "-I", str(root / "src" / "include" / "internal"),
        "-I", str(root / "src" / "include" / "external"),
        "-I", str(root / "build" / "prefix" / "include"),
        str(root / "scripts" / "bop_insn_chain.cpp"),
        "-L", str(so_dir), "-lyodaLib",
        "-Wl,-rpath-link," + pieces["sphericart_lib"],
        "-Wl,-rpath-link," + pieces["minimage_lib"],
        "-Wl,-rpath," + str(so_dir),
        "-Wl,-rpath," + pieces["sphericart_lib"],
        "-Wl,-rpath," + pieces["minimage_lib"],
        "-Wl,-rpath," + lib64,
        "-o", str(exe),
    ]
    ran = run(cmd, root)
    if ran.returncode != 0:
        sys.stderr.write(ran.stderr[-4000:])
        raise SystemExit(ran.returncode)


def measure(exe, work, label, lib_path):
    valgrind = shutil.which("valgrind")
    if not valgrind:
        raise SystemExit("valgrind is not on PATH")
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = lib_path
    env["OMP_NUM_THREADS"] = "1"
    env["GCOV_PREFIX"] = str(work / "gcov")
    env["GCOV_PREFIX_STRIP"] = "20"
    counts = {}
    refs = []
    for n in (N, 2 * N):
        sys.stderr.write("callgrind %s %d\n" % (label, n))
        sys.stderr.flush()
        out = work / ("%s-%d.cg" % (label, n))
        cmd = [
            valgrind, "--tool=callgrind",
            "--callgrind-out-file=%s" % out,
            "--dump-instr=no", "--collect-jumps=no", "--cache-sim=no",
            str(exe), str(n),
        ]
        ran = run(cmd, work, env)
        if ran.returncode != 0:
            sys.stderr.write(scrub(ran.stderr)[-2000:])
            raise SystemExit("%s N=%d exited %d" % (label, n, ran.returncode))
        counts[n] = summary_of(out)
        for line in (ran.stderr + ran.stdout).splitlines():
            if "I   refs:" in line or line.startswith("==") and "Collected" in line:
                refs.append(scrub(line).rstrip())
    ratio = counts[2 * N] / counts[N]
    body = [
        "tool: valgrind --tool=callgrind",
        "event: Ir",
        "chain: %d and %d atoms, one neighbour each, degree 6" % (N, 2 * N),
        "threads: 1",
        "sphericart: SPHERICART_ARCH_NATIVE off",
        "instructions_%d: %d" % (N, counts[N]),
        "instructions_%d: %d" % (2 * N, counts[2 * N]),
        "ratio: %.6f" % ratio,
        "bound: %.1f" % BOUND,
        "result: %s" % ("pass" if ratio <= BOUND else "fail"),
        "",
    ]
    body.extend(line for line in refs if line)
    body.append("")
    return ratio, "\n".join(body)


def pieces_from(root, build):
    sph = os.environ.get("SPHERICART_LIB", "")
    if not sph:
        sph = str(root / "build" / "prefix" / "lib")
    return {
        "root": root,
        "sphericart_lib": sph,
        "minimage_lib": str(build / "subprojects" / "minimage"),
        "ozaki": build / "src" / "libyodaLib.so.p" / "ozaki.cpp.o",
        "sphericart_ylm": build / "src" / "libyodaLib.so.p" / "sphericart_ylm.cpp.o",
    }


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "usage: bop_insn_ratio.py PARENT_SHA BUILD RECEIPT_DIR\n"
        )
        return 2
    parent_sha, build_name, receipt_dir = sys.argv[1:4]
    root = Path(__file__).resolve().parents[1]
    build = Path(build_name)
    if not build.is_absolute():
        build = root / build
    entry = compile_entry(build)
    shown = run(
        ["git", "show", "%s:src/steinhardt_parallel.cpp" % parent_sha], root
    )
    if shown.returncode != 0:
        sys.stderr.write(shown.stderr)
        return shown.returncode
    if "owner[static_cast<size_t>(k)] == i" not in shown.stdout:
        sys.stderr.write("parent source has no bond rescan\n")
        return 2
    tip_src = root / "src" / "steinhardt_parallel.cpp"
    if "used[static_cast<size_t>(i - begin)]" not in tip_src.read_text():
        sys.stderr.write("tip source does not store the bond count\n")
        return 2
    # Keep the objects outside the suite build so its gcda stays untouched.
    work = root / "build" / "bop-insn"
    work.mkdir(parents=True, exist_ok=True)
    parent_cpp = work / "parent.cpp"
    parent_cpp.write_text(shown.stdout)
    pieces = pieces_from(root, build)
    for key in ("ozaki", "sphericart_ylm"):
        if not Path(pieces[key]).is_file():
            sys.stderr.write("missing %s\n" % key)
            return 2
    out_dir = Path(receipt_dir)
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    compiler = compile_unit(entry, parent_cpp, work / "parent.o")
    parent_so = work / "parent" / "libyodaLib.so"
    parent_so.parent.mkdir(parents=True, exist_ok=True)
    lib64 = link_library(compiler, work / "parent.o", parent_so, pieces)
    parent_exe = work / "parent" / "chain"
    compile_driver(compiler, pieces, parent_so.parent, parent_exe, lib64)
    lib_path = ":".join([
        str(parent_so.parent), pieces["sphericart_lib"],
        pieces["minimage_lib"], lib64,
    ])
    parent_ratio, parent_body = measure(parent_exe, work, "parent", lib_path)
    (out_dir / "bop_insn_parent.txt").write_text(parent_body)

    compile_unit(entry, tip_src, work / "tip.o")
    tip_so = work / "tip" / "libyodaLib.so"
    tip_so.parent.mkdir(parents=True, exist_ok=True)
    link_library(compiler, work / "tip.o", tip_so, pieces)
    tip_exe = work / "tip" / "chain"
    compile_driver(compiler, pieces, tip_so.parent, tip_exe, lib64)
    tip_lib_path = ":".join([
        str(tip_so.parent), pieces["sphericart_lib"],
        pieces["minimage_lib"], lib64,
    ])
    tip_ratio, tip_body = measure(tip_exe, work, "tip", tip_lib_path)
    (out_dir / "bop_insn_tip.txt").write_text(tip_body)
    print("parent %.6f tip %.6f" % (parent_ratio, tip_ratio))
    if parent_ratio <= BOUND or tip_ratio > BOUND:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
