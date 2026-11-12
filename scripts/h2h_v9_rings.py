#!/usr/bin/env python3
"""Ring census on ideal ice Ic, one session.

d-SEAMS cutoff primitive rings, OVITO Find rings, SQQ shortest-path
rings, cycless and fastatomstruct King rings. Sizes 10648 and 110592.
One and eight threads. One warm-up, six calls, median and type-7 IQR.
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
HOST = "/usr/lib/x86_64-linux-gnu"
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)
LIB = str(ROOT / "bbdir-release" / "bindings" / "libseams_lang.so")
FIELDS = [
    "threads", "code", "version", "operation", "molecules", "warmup", "n",
    "median_ms", "q1_ms", "q3_ms", "iqr_ms", "agreement", "detail",
    "samples_ms", "cpus",
]


def lib_env():
    build = Path(LIB).resolve().parents[1]
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


def child_env(kind, threads):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["RAYON_NUM_THREADS"] = width
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env.pop("OMP_PROC_BIND", None)
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lib_env()
    elif kind == "ovito":
        env["LD_LIBRARY_PATH"] = HOST
        env["OVITO_THREAD_COUNT"] = width
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        env.pop("LD_LIBRARY_PATH", None)
        env.pop("OVITO_THREAD_COUNT", None)
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
        lib = ctypes.CDLL(LIB)
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
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
            counts = (ctypes.c_int * (proto.MAX_RING + 1))()
            status = lib.seams_frame_rings(
                handle, proto.CUTOFF_A, 1, b"cutoff", proto.MAX_RING, counts,
                proto.MAX_RING + 1)
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return [int(counts[i]) for i in range(3, proto.MAX_RING + 1)]

        samples, census = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": lib.seams_lang_version().decode(),
            "operation": "rings", "molecules": n, "samples_ms": samples,
            "census": census, "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "ovito":
        os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        import ovito
        from ovito.data import DataCollection
        from ovito.modifiers import CreateBondsModifier, FindRingsModifier

        def prep():
            data = DataCollection()
            data.create_cell(
                matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
                pbc=(True, True, True))
            particles = data.create_particles(count=n)
            particles.create_property("Position", data=np.ascontiguousarray(pos))
            return data

        def call(frame):
            frame.apply(CreateBondsModifier(cutoff=proto.CUTOFF_A))
            frame.apply(FindRingsModifier(
                minimum_ring_size=3, maximum_ring_size=proto.MAX_RING))
            hist = frame.tables["ring-size-histogram"]
            found = dict(zip(
                [int(v) for v in np.asarray(hist["Ring Size"])],
                [int(v) for v in np.asarray(hist["Count"])]))
            return [found.get(size, 0) for size in range(3, proto.MAX_RING + 1)]

        samples, census = proto.time_calls(prep, call)
        emit({
            "code": "OVITO", "version": ovito.version_string,
            "operation": "rings", "molecules": n, "samples_ms": samples,
            "census": census, "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "sqq":
        from sqq.api import analyze_frame, load_config
        from sqq.models import Atom, Frame
        import sqq
        pos_nm = np.asarray(pos, dtype=float) / 10.0
        side_nm = side / 10.0
        atoms = [
            Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
                 xyz=pos_nm[i], molecule_id=i + 1)
            for i in range(n)
        ]
        frame = Frame(name=f"ice-{n}", atoms=atoms, box=np.array([side_nm] * 3))
        cfg = load_config({
            "engine": "py",
            "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
            "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
            "half_cage": {"enabled": False},
            "quasi_cage": {"enabled": False},
            "cage": {"enabled": False},
            "order_parameter": {"enabled": []},
            "parallel": {"workers": args.threads, "math_threads": 1},
        }, engine="py")

        def call(_f):
            result = analyze_frame(frame, cfg, engine="py")
            return [len(result.rings.get(size, [])) for size in (4, 5, 6)]

        samples, census = proto.time_calls(lambda: frame, call)
        emit({
            "code": "SQQ", "version": sqq._version.__version__,
            "operation": "rings", "molecules": n, "samples_ms": samples,
            "census": census, "cpus": pinned, "threads": args.threads,
            "note": "sizes 4-6",
        })
        return

    if args.worker == "cycless":
        import importlib.metadata
        import networkx as nx
        import pairlist
        from cycless.cycles import cycles_iter

        def call(points):
            fractional = points / side
            graph = nx.Graph()
            graph.add_nodes_from(range(n))
            graph.add_edges_from(pairlist.pairs_iter(
                fractional, proto.CUTOFF_A, np.diag([side, side, side]),
                distance=False))
            census = [0, 0, 0, 0]
            for ring in cycles_iter(graph, proto.MAX_RING, pos=fractional):
                size = len(ring)
                if 3 <= size <= proto.MAX_RING:
                    census[size - 3] += 1
            return census

        samples, census = proto.time_calls(lambda: pos, call)
        emit({
            "code": "cycless", "version": importlib.metadata.version("cycless"),
            "operation": "rings", "molecules": n, "samples_ms": samples,
            "census": census, "cpus": pinned, "threads": args.threads,
        })
        return

    import fastatomstruct as fs
    fs.disable_telemetry()
    fs.suppress_telemetry_notice()
    fs.set_num_threads(args.threads)
    elements = [fs.Element.O] * n
    atoms = fs.Atoms(pos, np.eye(3) * side, elements)
    cutoff = fs.analysis.distance.CutoffMode.Fixed(proto.CUTOFF_A)
    method = fs.analysis.rings.RingMethod.King

    def call(_a):
        sizes, counts = fs.analysis.rings.ring_size_distribution(
            atoms, cutoff, proto.MAX_RING, method)
        found = {int(s): int(round(float(c))) for s, c in zip(sizes, counts)}
        return [found.get(size, 0) for size in range(3, proto.MAX_RING + 1)]

    samples, census = proto.time_calls(lambda: atoms, call)
    emit({
        "code": "fastatomstruct", "version": fs.__version__,
        "operation": "rings", "molecules": n, "samples_ms": samples,
        "census": census, "cpus": pinned, "threads": int(fs.get_num_threads()),
    })


def launch(python, threads, worker, frame):
    cores = "0" if threads == 1 else "0-7"
    kind = "dseams" if worker == "dseams" else ("ovito" if worker == "ovito" else "other")
    done = subprocess.run(
        ["taskset", "-c", cores, python, str(HERE),
         "--worker", worker, "--threads", str(threads), "--frame", frame],
        cwd=ROOT, env=child_env(kind, threads), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": worker, "threads": threads, "molecules": 0}
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        return {"error": f"not pinned: {payload.get('cpus')}", "code": worker,
                "threads": threads}
    return payload


def expect_of(payload, molecules):
    census = payload.get("census") or []
    if payload.get("code") == "SQQ":
        return [0, 0, 2 * molecules]
    return [0, 0, 0, 2 * molecules]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/h2h-watch-v9-rings.csv")
    parser.add_argument("--python")
    parser.add_argument("--worker")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--frame")
    args = parser.parse_args()
    if args.worker:
        run_worker(args)
        return
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
        path = ROOT / "benchmarks" / f".ice-ic-{len(pos)}.npz"
        np.savez(path, positions=pos, side=side)
        frames[len(pos)] = str(path)
    workers = ("dseams", "ovito", "sqq", "fastatomstruct", "cycless")
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for molecules in (10648, 110592):
            for threads in (1, 8):
                for worker in workers:
                    payload = launch(python, threads, worker, frames[molecules])
                    payload["threads"] = threads
                    payload["molecules"] = molecules
                    if "error" in payload:
                        writer.writerow({
                            "threads": threads, "code": payload.get("code", worker),
                            "version": "", "operation": "rings", "molecules": molecules,
                            "warmup": "", "n": 0, "median_ms": "", "q1_ms": "",
                            "q3_ms": "", "iqr_ms": "", "agreement": "",
                            "detail": payload["error"], "samples_ms": "", "cpus": "",
                        })
                    else:
                        summary = proto.summarize(payload["samples_ms"])
                        census = payload.get("census")
                        same = census == expect_of(payload, molecules)
                        writer.writerow({
                            "threads": threads, "code": payload["code"],
                            "version": payload.get("version", ""),
                            "operation": "rings", "molecules": molecules,
                            "warmup": proto.WARMUP, "n": summary["n"],
                            "median_ms": f"{summary['median_ms']:.3f}",
                            "q1_ms": f"{summary['q1_ms']:.3f}",
                            "q3_ms": f"{summary['q3_ms']:.3f}",
                            "iqr_ms": f"{summary['iqr_ms']:.3f}",
                            "agreement": f"{float(same):.6f}",
                            "detail": f"census={census}",
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
