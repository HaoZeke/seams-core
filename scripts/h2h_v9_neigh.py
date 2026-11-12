#!/usr/bin/env python3
"""Cutoff neighbour lists, d-SEAMS 064bfef7 against vesin 0.6.2.

Ideal ice Ic, cutoff 3.5 A, both directions, self excluded. Sizes
10648 and 110592. Thread counts 1 through 8. One warm-up, six calls,
median and type-7 IQR. The d-SEAMS library is the release build of
perf/neighbours-threads, passed with --lib.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)
FIELDS = [
    "threads", "code", "version", "operation", "molecules", "warmup", "n",
    "median_ms", "q1_ms", "q3_ms", "iqr_ms", "agreement", "detail",
    "samples_ms", "cpus",
]


def lib_env(lib):
    build = Path(lib).resolve().parents[1]
    return os.pathsep.join([
        str(build / "bindings"), str(build / "src"),
        str(build / "subprojects/readcon-core"),
        str(build / "subprojects/readcon-db"),
        str(build / "subprojects/linkcell"),
        str(build / "subprojects/vesin"),
        str(build / "subprojects/minimage"),
        FLEXI,
    ])


def cpus_now():
    for line in open("/proc/self/status", encoding="utf-8"):
        if line.startswith("Cpus_allowed_list:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def child_env(kind, threads, lib):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["RAYON_NUM_THREADS"] = width
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env.pop("OMP_PROC_BIND", None)
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lib_env(lib)
    else:
        env.pop("LD_LIBRARY_PATH", None)
    return env


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def run_worker(args):
    pinned = cpus_now()
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)
    sys.path.insert(0, str(HERE.parent))
    import ctypes
    import numpy as np
    import h2h_protocol as proto
    data = np.load(args.frame)
    pos = data["positions"]
    side = float(data["side"])
    n = len(pos)

    if args.worker == "dseams":
        lib = ctypes.CDLL(args.lib)
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        xyz = np.ascontiguousarray(pos, dtype=np.float64)
        box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        handle = lib.seams_frame_from_arrays(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if not handle:
            raise RuntimeError(lib.seams_last_error().decode())

        def call(_h):
            status = lib.seams_frame_neighbours(handle, proto.CUTOFF_A, 1, b"cutoff")
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return status

        samples, _ = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": lib.seams_lang_version().decode(),
            "operation": "neighbours", "molecules": n, "samples_ms": samples,
            "rows": n, "cpus": pinned, "threads": args.threads,
            "note": "status 0 means one row per atom; the C entry discards partners",
        })
        return

    import vesin
    from vesin import NeighborList
    points = np.ascontiguousarray(pos, dtype=np.float64)
    boxm = np.eye(3, dtype=np.float64) * side
    nlist = NeighborList(cutoff=proto.CUTOFF_A, full_list=True, n_threads=args.threads)

    def call(_p):
        i, j = nlist.compute(points=points, box=boxm, periodic=True, quantities="ij")
        return int(len(i))

    samples, pairs = proto.time_calls(lambda: points, call)
    emit({
        "code": "vesin", "version": vesin.__version__,
        "operation": "neighbours", "molecules": n, "samples_ms": samples,
        "pairs": pairs, "cpus": pinned, "threads": args.threads,
    })


def launch(python, threads, worker, frame, lib):
    cores = "0" if threads == 1 else f"0-{threads - 1}"
    kind = "dseams" if worker == "dseams" else "other"
    done = subprocess.run(
        ["taskset", "-c", cores, python, str(HERE),
         "--worker", worker, "--threads", str(threads), "--frame", frame,
         "--lib", lib],
        cwd=ROOT, env=child_env(kind, threads, lib), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": worker, "threads": threads}
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        return {"error": f"not pinned: {payload.get('cpus')}", "code": worker,
                "threads": threads}
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/h2h-watch-v9-neigh.csv")
    parser.add_argument("--python")
    parser.add_argument("--lib", default="")
    parser.add_argument("--worker")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--frame")
    args = parser.parse_args()
    if args.worker:
        run_worker(args)
        return
    if not args.lib:
        raise SystemExit("need --lib")
    sys.path.insert(0, str(HERE.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    import numpy as np
    import h2h_protocol as proto
    python = args.python or sys.executable
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    frames = {}
    for reps in (11, 24):
        pos, side = proto.ice_ic(reps)
        path = ROOT / "benchmarks" / f".ice-ic-neigh-{len(pos)}.npz"
        np.savez(path, positions=pos, side=side)
        frames[len(pos)] = str(path)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for molecules in (10648, 110592):
            for threads in range(1, 9):
                for worker in ("dseams", "vesin"):
                    payload = launch(python, threads, worker, frames[molecules], args.lib)
                    if "error" in payload:
                        writer.writerow({
                            "threads": threads, "code": payload.get("code", worker),
                            "version": "", "operation": "neighbours",
                            "molecules": molecules, "warmup": "", "n": 0,
                            "median_ms": "", "q1_ms": "", "q3_ms": "", "iqr_ms": "",
                            "agreement": "", "detail": payload["error"],
                            "samples_ms": "", "cpus": "",
                        })
                    else:
                        summary = proto.summarize(payload["samples_ms"])
                        counted = payload.get("pairs")
                        if counted is None:
                            same = payload.get("rows") == molecules
                            detail = "rows=%s; partners discarded by the C entry" % payload.get("rows")
                        else:
                            same = int(counted) == 4 * molecules
                            detail = f"pairs={counted}"
                        writer.writerow({
                            "threads": threads, "code": payload["code"],
                            "version": payload.get("version", ""),
                            "operation": "neighbours", "molecules": molecules,
                            "warmup": proto.WARMUP, "n": summary["n"],
                            "median_ms": f"{summary['median_ms']:.3f}",
                            "q1_ms": f"{summary['q1_ms']:.3f}",
                            "q3_ms": f"{summary['q3_ms']:.3f}",
                            "iqr_ms": f"{summary['iqr_ms']:.3f}",
                            "agreement": f"{float(same):.6f}",
                            "detail": detail,
                            "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
                            "cpus": payload.get("cpus", ""),
                        })
                    handle.flush()
                    print(f"N={molecules} threads={threads} {worker} "
                          f"{payload.get('samples_ms', payload.get('error'))}", flush=True)
    for path in frames.values():
        Path(path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
