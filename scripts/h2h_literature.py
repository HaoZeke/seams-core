#!/usr/bin/env python3
"""Head-to-head of d-SEAMS against mdapy, pyscal3, vesin, IRA and SQQ.

The statistics are h2h_protocol: one warm-up, eight timed calls, the
median of an even sample, and type-7 quartiles. One thread is pinned to
cpu 0. Eight threads are pinned to cpus 0-7. The allowed cpu list is
recorded. Host names are not. Frames are built in memory.

Operations are the cutoff neighbour list, local q4 and q6 on that list,
the diamond common-neighbour labels, metal CNA, primitive rings, and
local ice symmetry against IRA. CHILL and CHILL+ are not timed and
their cuts are not changed.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import gen_ideal_ice as ice
import h2h_protocol as proto

CUTOFF = proto.CUTOFF_A
REPEATS = 8
HAUSDORFF_A = 0.2
FIELDS = [
    "threads", "cores", "cpus", "phase", "molecules", "operation",
    "code", "version", "warmup", "n", "median_ms", "q1_ms", "q3_ms",
    "iqr_ms", "agreement", "detail", "samples_ms",
]


def project_version():
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    match = re.search(r"version:\s*'([0-9.]+)'", text)
    return match.group(1) if match else "unknown"


def core_count(cpus):
    total = 0
    for part in str(cpus).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            total += int(hi) - int(lo) + 1
        else:
            total += 1
    return total


def scrub(text):
    return re.sub(r"(?<![\w.])/[\w.+@=-]+(?:/[\w.+@=-]+)+", "<path>", str(text))


def fmt_samples(samples):
    return " ".join(f"{float(v):.3f}" for v in samples)


def row_of(threads, cpus, frame, operation, code, version, summary, agreement, detail, samples):
    return {
        "threads": threads,
        "cores": core_count(cpus),
        "cpus": cpus,
        "phase": frame["phase"],
        "molecules": frame["molecules"],
        "operation": operation,
        "code": code,
        "version": version,
        "warmup": proto.WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "q1_ms": f"{summary['q1_ms']:.3f}",
        "q3_ms": f"{summary['q3_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "agreement": "" if agreement is None else f"{float(agreement):.6f}",
        "detail": scrub(detail),
        "samples_ms": fmt_samples(samples),
    }


def load_core(path):
    lib = ctypes.CDLL(path)
    c_double_p = ctypes.POINTER(ctypes.c_double)
    c_int_p = ctypes.POINTER(ctypes.c_int)
    lib.seams_cutoff_csr.argtypes = [
        c_double_p, ctypes.c_int, c_double_p, ctypes.c_double,
        c_int_p, c_int_p, ctypes.c_int, c_int_p,
    ]
    lib.seams_cutoff_csr.restype = ctypes.c_int
    lib.seams_ql_cutoff.argtypes = [
        c_double_p, ctypes.c_int, c_double_p, ctypes.c_double,
        ctypes.c_int, c_double_p,
    ]
    lib.seams_ql_cutoff.restype = ctypes.c_int
    lib.seams_diamond_cna.argtypes = [
        c_double_p, ctypes.c_int, c_double_p, ctypes.c_double, c_int_p,
    ]
    lib.seams_diamond_cna.restype = ctypes.c_int
    lib.seams_primitive_rings.argtypes = [
        c_double_p, ctypes.c_int, c_double_p, ctypes.c_double,
        ctypes.c_int, c_int_p, ctypes.c_int,
    ]
    lib.seams_primitive_rings.restype = ctypes.c_int
    lib.seams_symmetry_classify.argtypes = [
        c_double_p, ctypes.c_int, c_double_p, c_int_p,
    ]
    lib.seams_symmetry_classify.restype = ctypes.c_int
    lib.seams_symmetry_template.argtypes = [
        ctypes.c_int, c_double_p, ctypes.c_int, c_int_p,
    ]
    lib.seams_symmetry_template.restype = ctypes.c_int
    lib.seams_symmetry_backend.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.seams_symmetry_backend.restype = ctypes.c_int
    return lib


def pointers(frame):
    xyz = np.ascontiguousarray(frame["positions"], dtype=np.float64)
    box = np.ascontiguousarray(frame["box9"], dtype=np.float64)
    lattice = np.ascontiguousarray(frame["lattice"], dtype=np.float64).reshape(-1)
    return xyz, box, lattice


def rows_from_csr(offsets, cols):
    rows = []
    for i in range(len(offsets) - 1):
        block = cols[int(offsets[i]):int(offsets[i + 1])]
        rows.append(tuple(sorted({int(v) for v in block})))
    return rows


def rows_from_pairs(n, left, right):
    buckets = [[] for _ in range(n)]
    for a, b in zip(left, right):
        buckets[int(a)].append(int(b))
    return [tuple(sorted(set(bucket))) for bucket in buckets]


def rows_from_lists(neigh):
    return [tuple(sorted({int(v) for v in row})) for row in neigh]


def pyscal_neighbor_rows(atoms):
    if "pyscal_neighbors" in atoms.arrays:
        return rows_from_lists(atoms.arrays["pyscal_neighbors"])
    if "pyscal_neighbors" in atoms.info:
        return rows_from_lists(atoms.info["pyscal_neighbors"])
    raise KeyError("pyscal_neighbors")


def set_agreement(reference, other):
    if len(reference) != len(other):
        return 0.0, f"length {len(reference)} vs {len(other)}"
    same = sum(a == b for a, b in zip(reference, other))
    differ = len(reference) - same
    return same / len(reference), f"atoms {differ} differ"


def q_agreement(reference, other):
    ref = np.asarray(reference, dtype=np.float64)
    got = np.asarray(other, dtype=np.float64)
    if ref.shape != got.shape:
        return 0.0, f"shape {ref.shape} vs {got.shape}"
    delta = np.abs(ref - got)
    worst = float(np.max(delta)) if delta.size else 0.0
    frac = float(np.mean(delta <= 1e-6)) if delta.size else 1.0
    return frac, f"max_abs={worst:.3e} within_1e-6={frac:.6f}"


def phase_word(code, source):
    code = int(code)
    if source == "dseams":
        return {1: "cubic", 0: "hexagonal"}.get(code, "other")
    if source == "diamond":
        return {1: "cubic", 4: "hexagonal", 0: "other"}.get(code, "neighbor")
    return {0: "other", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico"}.get(code, "other")


def label_agreement(reference, other, ref_source, other_source):
    if len(reference) != len(other):
        return 0.0, f"length {len(reference)} vs {len(other)}"
    same = 0
    counts = {}
    for a, b in zip(reference, other):
        left = phase_word(a, ref_source)
        right = phase_word(b, other_source)
        if left == right:
            same += 1
        key = f"{left}/{right}"
        counts[key] = counts.get(key, 0) + 1
    top = sorted(counts.items(), key=lambda item: -item[1])[:4]
    shown = " ".join(f"{name}:{count}" for name, count in top)
    return same / len(reference), shown


def time_pair(prep, run):
    samples, last = proto.time_calls(prep, run, warmup=proto.WARMUP, repeats=REPEATS)
    return proto.summarize(samples), samples, last


class Engine:
    def __init__(self, lib_path):
        self.lib = load_core(lib_path)
        self.version = project_version()
        buf = ctypes.create_string_buffer(32)
        self.lib.seams_symmetry_backend(buf, 32)
        self.backend = buf.value.decode()

    def neighbours(self, frame):
        xyz, box, _lattice = pointers(frame)
        n = frame["molecules"]
        offsets = np.empty(n + 1, dtype=np.int32)
        need = ctypes.c_int()
        status = self.lib.seams_cutoff_csr(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CUTOFF,
            offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None, 0, ctypes.byref(need))
        if status != 0:
            raise RuntimeError(f"cutoff list query returned {status}")
        cols = np.empty(int(need.value), dtype=np.int32)
        got = ctypes.c_int()

        def run(_unused):
            code = self.lib.seams_cutoff_csr(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CUTOFF,
                offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                cols.size, ctypes.byref(got))
            if code != 0:
                raise RuntimeError(f"cutoff list returned {code}")
            return int(got.value)

        summary, samples, _last = time_pair(lambda: None, run)
        rows = rows_from_csr(offsets, cols[:int(got.value)])
        return summary, samples, rows

    def ql(self, frame, order):
        xyz, box, _lattice = pointers(frame)
        n = frame["molecules"]
        values = np.empty(n, dtype=np.float64)

        def run(_unused):
            code = self.lib.seams_ql_cutoff(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CUTOFF,
                order, values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            if code != 0:
                raise RuntimeError(f"q{order} returned {code}")
            return None

        summary, samples, _last = time_pair(lambda: None, run)
        return summary, samples, values.copy()

    def diamond(self, frame):
        xyz, box, _lattice = pointers(frame)
        n = frame["molecules"]
        labels = np.empty(n, dtype=np.int32)

        def run(_unused):
            code = self.lib.seams_diamond_cna(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CUTOFF,
                labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
            if code != 0:
                raise RuntimeError(f"diamond returned {code}")
            return None

        summary, samples, _last = time_pair(lambda: None, run)
        return summary, samples, labels.copy()

    def rings(self, frame):
        xyz, box, _lattice = pointers(frame)
        n = frame["molecules"]
        counts = np.zeros(7, dtype=np.int32)

        def run(_unused):
            code = self.lib.seams_primitive_rings(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), CUTOFF,
                6, counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), 7)
            if code != 0:
                raise RuntimeError(f"rings returned {code}")
            return None

        summary, samples, _last = time_pair(lambda: None, run)
        return summary, samples, counts.copy()

    def symmetry(self, frame):
        xyz, _box, lattice = pointers(frame)
        n = frame["molecules"]
        labels = np.empty(n, dtype=np.int32)

        def run(_unused):
            code = self.lib.seams_symmetry_classify(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
            if code != 0:
                raise RuntimeError(f"symmetry returned {code}")
            return None

        summary, samples, _last = time_pair(lambda: None, run)
        return summary, samples, labels.copy()

    def templates(self):
        found = []
        for ident in (0, 1, 2):
            n_out = ctypes.c_int()
            code = self.lib.seams_symmetry_template(ident, None, 0, ctypes.byref(n_out))
            if code != 0 or n_out.value <= 0:
                continue
            buf = np.empty(int(n_out.value) * 3, dtype=np.float64)
            code = self.lib.seams_symmetry_template(
                ident, buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                int(n_out.value), ctypes.byref(n_out))
            if code != 0:
                continue
            found.append((ident, buf.reshape(-1, 3).copy()))
        return found


def mdapy_system(frame):
    from mdapy import System
    return System(positions=frame["positions"], box=frame["lattice"], elements="O")


def time_mdapy(frame, operation, threads):
    from mdapy import set_num_threads
    from mdapy.analysis import RingStatistics, SteinhardtParameters
    from mdapy.identify import CommonNeighborAnalysis, DiamondStructure
    from mdapy.neighbor import NeighborList
    import mdapy
    set_num_threads(threads)
    system = mdapy_system(frame)
    if operation == "neighbours":
        def run(_unused):
            built = NeighborList.build(system, cutoff=CUTOFF)
            return np.asarray(built.offsets), np.asarray(built.j)
        summary, samples, last = time_pair(lambda: system, run)
        rows = rows_from_csr(last[0], last[1])
        return summary, samples, rows, mdapy.__version__
    if operation in ("q4", "q6"):
        order = 4 if operation == "q4" else 6
        tool = SteinhardtParameters(degrees=(order,), cutoff=CUTOFF)

        def run(_unused, tool=tool):
            result = tool.compute(system)
            values = np.asarray(result.q, dtype=np.float64)
            if values.ndim == 2:
                values = values[:, 0]
            return values
        summary, samples, last = time_pair(lambda: system, run)
        return summary, samples, last, mdapy.__version__
    if operation == "diamond":
        tool = DiamondStructure()

        def run(_unused, tool=tool):
            return np.asarray(tool.compute(system).structure)
        summary, samples, last = time_pair(lambda: system, run)
        return summary, samples, last, mdapy.__version__
    if operation == "cna":
        tool = CommonNeighborAnalysis(cutoff=CUTOFF)

        def run(_unused, tool=tool):
            return np.asarray(tool.compute(system).structure)
        summary, samples, last = time_pair(lambda: system, run)
        return summary, samples, last, mdapy.__version__
    if operation == "rings":
        tool = RingStatistics(max_size=6, min_size=3, cutoff=CUTOFF)

        def run(_unused, tool=tool):
            return np.asarray(tool.compute(system).counts)
        summary, samples, last = time_pair(lambda: system, run)
        return summary, samples, last, mdapy.__version__
    raise ValueError(operation)


def pyscal_atoms(frame):
    from ase import Atoms
    return Atoms(
        symbols=["O"] * frame["molecules"],
        positions=frame["positions"],
        cell=frame["lattice"],
        pbc=True,
    )


def time_pyscal(frame, operation, threads):
    import pyscal3
    pyscal3.set_num_threads(threads)
    atoms = pyscal_atoms(frame)
    if operation == "neighbours":
        def run(_unused):
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
            return pyscal_neighbor_rows(atoms)
        summary, samples, last = time_pair(lambda: atoms, run)
        return summary, samples, last, pyscal3.__version__
    if operation in ("q4", "q6"):
        order = 4 if operation == "q4" else 6

        def run(_unused):
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
            return np.asarray(pyscal3.steinhardt_parameter(atoms, order)[0], dtype=np.float64)
        summary, samples, last = time_pair(lambda: atoms, run)
        return summary, samples, last, pyscal3.__version__
    if operation == "diamond":
        def run(_unused):
            pyscal3.diamond_structure(atoms)
            return np.asarray(atoms.arrays["pyscal_structure"])
        summary, samples, last = time_pair(lambda: atoms, run)
        return summary, samples, last, pyscal3.__version__
    if operation == "cna":
        def run(_unused):
            pyscal3.common_neighbor_analysis(atoms)
            return np.asarray(atoms.arrays["pyscal_structure"])
        summary, samples, last = time_pair(lambda: atoms, run)
        return summary, samples, last, pyscal3.__version__
    raise ValueError(operation)


def time_vesin(frame, threads):
    import vesin
    points = frame["positions"]
    box = frame["lattice"]
    finder = vesin.NeighborList(cutoff=CUTOFF, full_list=True, n_threads=threads)

    def run(_unused):
        left, right = finder.compute(
            points=points, box=box, periodic=True, quantities="ij", copy=True)
        return left, right

    summary, samples, last = time_pair(lambda: points, run)
    rows = rows_from_pairs(frame["molecules"], last[0], last[1])
    return summary, samples, rows, vesin.__version__


def time_sqq(frame, threads):
    import sqq
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame
    pos_nm = np.asarray(frame["positions"], dtype=float) / 10.0
    side = frame["box9"]
    # The ideal cells used here are orthorhombic, which is the box SQQ stores.
    box_nm = np.array([side[3], side[4], side[5]], dtype=float) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(len(pos_nm))
    ]
    built = Frame(name=frame["name"], atoms=atoms, box=box_nm)
    cfg = load_config({
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": CUTOFF / 10.0},
        "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": False},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": threads, "math_threads": 1},
    }, engine="py")

    def run(_unused):
        result = analyze_frame(built, cfg, engine="py")
        return [len(result.rings.get(size, [])) for size in (4, 5, 6)]

    summary, samples, last = time_pair(lambda: built, run)
    return summary, samples, last, getattr(sqq, "__version__", "unknown")


def shell_table(frame):
    import vesin
    points = frame["positions"]
    n = frame["molecules"]
    finder = vesin.NeighborList(cutoff=5.0, full_list=True, n_threads=1)
    left, dist, vec = finder.compute(
        points=points, box=frame["lattice"], periodic=True, quantities="idD", copy=True)
    order = np.argsort(left, kind="mergesort")
    left = left[order]
    dist = dist[order]
    vec = vec[order]
    cuts = np.searchsorted(left, np.arange(n + 1))
    shells = np.full((n, 17, 3), np.nan, dtype=np.float64)
    for i in range(n):
        start, stop = int(cuts[i]), int(cuts[i + 1])
        if stop - start < 16:
            continue
        pick = np.argsort(dist[start:stop])[:16]
        if float(dist[start:stop][pick[-1]]) > 5.0:
            continue
        shells[i, 0] = points[i]
        shells[i, 1:] = points[i] + vec[start:stop][pick]
    return shells


def ira_match(frame, templates, threads):
    shells = shell_table(frame)
    usable = [points for ident, points in templates if len(points) == 17]
    idents = [ident for ident, points in templates if len(points) == 17]
    if not usable:
        raise RuntimeError("no 17-atom symmetry template")
    payload = tempfile.NamedTemporaryFile(suffix=".npz", delete=False)
    result = tempfile.NamedTemporaryFile(suffix=".npz", delete=False)
    payload.close()
    result.close()
    try:
        np.savez(
            payload.name, shells=shells,
            templates=np.stack(usable),
            idents=np.asarray(idents, dtype=np.int32),
            threads=np.asarray([threads], dtype=np.int32),
        )
        env = os.environ.copy()
        env["H2H_PINNED"] = "1"
        # pixi runs another Python. A caller module path would load the
        # wrong numpy extension into that interpreter.
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
        env.pop("LD_LIBRARY_PATH", None)
        cmd = [
            "pixi", "exec", "--spec", "ira", "python", str(Path(__file__).resolve()),
            "--worker", "ira", "--payload", payload.name, "--result", result.name,
        ]
        proc = subprocess.run(cmd, env=env, text=True, capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(scrub(proc.stderr[-2000:] or proc.stdout[-2000:]))
        got = np.load(result.name)
        return (
            [float(v) for v in got["samples"]],
            np.asarray(got["labels"]),
            np.asarray(got["hausdorff"]),
            str(got["version"]),
        )
    finally:
        for path in (payload.name, result.name):
            try:
                os.unlink(path)
            except OSError:
                pass


_SHELLS = None
_TEMPLATES = None
_IDENTS = None


def _match_range(bounds):
    import ira_mod
    start, stop = bounds
    ira = ira_mod.IRA()
    labels = np.empty(stop - start, dtype=np.int32)
    haus = np.empty(stop - start, dtype=np.float64)
    nat = 17
    types = np.ones(nat, dtype=np.int32)
    for offset, index in enumerate(range(start, stop)):
        shell = _SHELLS[index]
        if not np.isfinite(shell[0, 0]):
            labels[offset] = -1
            haus[offset] = np.nan
            continue
        best_h = 1e9
        best = -1
        for slot, templ in enumerate(_TEMPLATES):
            _rotation, _translation, _perm, distance = ira.match(
                nat, types, templ, nat, types, shell, 1.5)
            if float(distance) < best_h:
                best_h = float(distance)
                best = int(_IDENTS[slot])
        haus[offset] = best_h
        labels[offset] = best if best_h <= HAUSDORFF_A else -1
    return start, labels, haus


def ira_worker(payload_path, result_path):
    import multiprocessing as mp
    import ira_mod
    global _SHELLS, _TEMPLATES, _IDENTS
    data = np.load(payload_path)
    _SHELLS = data["shells"]
    _TEMPLATES = data["templates"]
    _IDENTS = data["idents"]
    threads = int(data["threads"][0])
    n = len(_SHELLS)
    width = max(threads, 1)
    edges = [int(v) for v in np.linspace(0, n, width + 1)]
    ranges = [(edges[i], edges[i + 1]) for i in range(width) if edges[i] < edges[i + 1]]
    pool = None
    if width > 1:
        pool = mp.get_context("fork").Pool(processes=width)

    def one_call():
        labels = np.full(n, -1, dtype=np.int32)
        haus = np.full(n, np.nan, dtype=np.float64)
        if pool is None:
            parts = [_match_range(ranges[0])]
        else:
            parts = pool.map(_match_range, ranges)
        for start, part_labels, part_haus in parts:
            labels[start:start + len(part_labels)] = part_labels
            haus[start:start + len(part_haus)] = part_haus
        return labels, haus

    try:
        samples, last = proto.time_calls(lambda: None, lambda _unused: one_call(),
                                         warmup=proto.WARMUP, repeats=REPEATS)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    labels, haus = last
    version = ira_mod.version if isinstance(ira_mod.version, str) else str(ira_mod.version)
    np.savez(result_path, labels=labels, hausdorff=haus,
             samples=np.asarray(samples, dtype=np.float64),
             version=np.asarray(version))


def emit(writer, handle, log, record):
    writer.writerow(record)
    handle.flush()
    log.write(
        "threads={threads} cores={cores} cpus={cpus} phase={phase} "
        "molecules={molecules} operation={operation} code={code} "
        "version={version} n={n} median_ms={median_ms} iqr_ms={iqr_ms} "
        "agreement={agreement} samples_ms={samples_ms} detail={detail}\n".format(**record)
    )
    log.flush()


def record_failure(writer, handle, log, threads, cpus, frame, operation, code, detail):
    record = {
        "threads": threads, "cores": core_count(cpus), "cpus": cpus,
        "phase": frame["phase"], "molecules": frame["molecules"],
        "operation": operation, "code": code, "version": "",
        "warmup": proto.WARMUP, "n": 0, "median_ms": "", "q1_ms": "",
        "q3_ms": "", "iqr_ms": "", "agreement": "",
        "detail": scrub(detail), "samples_ms": "",
    }
    emit(writer, handle, log, record)


def finished_rows(path):
    done = set()
    target = Path(path)
    if not target.exists() or target.stat().st_size == 0:
        return done
    with target.open(newline="") as handle:
        for record in csv.DictReader(handle):
            done.add((
                str(record["threads"]), record["phase"], str(record["molecules"]),
                record["operation"], record["code"],
            ))
    return done


def median_degree(frame):
    import vesin
    left, _right = vesin.NeighborList(cutoff=CUTOFF, full_list=True).compute(
        points=frame["positions"], box=frame["lattice"], periodic=True,
        quantities="ij", copy=True)
    counts = np.bincount(left, minlength=frame["molecules"])
    return float(np.median(counts))


def run_session(args, threads):
    cpus = proto.cpus_allowed()
    lib_path = args.lib or os.environ.get("SEAMS_CORE_LIB") or str(
        ROOT / "build" / "src" / "libyodaLib.so")
    engine = Engine(lib_path)
    frames = ice.campaign_frames(tuple(args.sizes), ROOT)
    done = finished_rows(args.out) if args.resume else set()
    if args.phase:
        frames = [frame for frame in frames if frame["phase"] == args.phase]
    if args.molecules:
        frames = [frame for frame in frames if frame["molecules"] == args.molecules]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    new_file = not out.exists() or out.stat().st_size == 0
    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a", newline="") as handle, log_path.open("a", encoding="utf-8") as log:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
            log.write(
                f"cpus={cpus} cores={core_count(cpus)} "
                f"d-SEAMS={engine.version} symmetry_backend={engine.backend} "
                f"repeats={REPEATS} warmup={proto.WARMUP} cutoff_A={CUTOFF}\n"
            )
        for frame in frames:
            if frame["phase"] in ("Ic", "Ih"):
                import vesin
                left, _right = vesin.NeighborList(cutoff=CUTOFF, full_list=True).compute(
                    points=frame["positions"], box=frame["lattice"], periodic=True,
                    quantities="ij", copy=True)
                counts = np.bincount(left, minlength=frame["molecules"])
                if not np.all(counts == 4):
                    raise RuntimeError(
                        f"{frame['name']} neighbour degree {int(counts.min())} to {int(counts.max())}"
                    )
                log.write(f"sanity {frame['name']} neighbours=4\n")
            run_frame(engine, writer, handle, log, threads, cpus, frame, args, done)


def run_frame(engine, writer, handle, log, threads, cpus, frame, args, done):
    wanted = set(args.only.split(",")) if args.only else None

    def keep(name):
        return wanted is None or name in wanted

    def seen(operation, code):
        return (str(threads), frame["phase"], str(frame["molecules"]), operation, code) in done

    ref_diamond = None
    if keep("neighbours") and not seen("neighbours", "d-SEAMS"):
        summary, samples, ref_rows = engine.neighbours(frame)
        emit(writer, handle, log, row_of(
            threads, cpus, frame, "neighbours", "d-SEAMS", engine.version,
            summary, 1.0, "reference cutoff list", samples))
        for code, caller in (
            ("mdapy", lambda: time_mdapy(frame, "neighbours", threads)),
            ("pyscal3", lambda: time_pyscal(frame, "neighbours", threads)),
            ("vesin", lambda: time_vesin(frame, threads)),
        ):
            if seen("neighbours", code):
                continue
            try:
                summary, samples, rows, version = caller()
                agree, detail = set_agreement(ref_rows, rows)
                emit(writer, handle, log, row_of(
                    threads, cpus, frame, "neighbours", code, version,
                    summary, agree, detail, samples))
            except Exception as exc:
                record_failure(writer, handle, log, threads, cpus, frame,
                               "neighbours", code, str(exc))
    for operation, order in (("q4", 4), ("q6", 6)):
        if not keep(operation) or seen(operation, "d-SEAMS"):
            continue
        summary, samples, values = engine.ql(frame, order)
        ref_q = values
        emit(writer, handle, log, row_of(
            threads, cpus, frame, operation, "d-SEAMS", engine.version,
            summary, 1.0, "reference local ql on the cutoff list", samples))
        for code, caller in (
            ("mdapy", lambda operation=operation: time_mdapy(frame, operation, threads)),
            ("pyscal3", lambda operation=operation: time_pyscal(frame, operation, threads)),
        ):
            if seen(operation, code):
                continue
            try:
                summary, samples, other, version = caller()
                agree, detail = q_agreement(ref_q, other)
                emit(writer, handle, log, row_of(
                    threads, cpus, frame, operation, code, version,
                    summary, agree, detail, samples))
            except Exception as exc:
                record_failure(writer, handle, log, threads, cpus, frame,
                               operation, code, str(exc))
    need_diamond = (keep("diamond") and not seen("diamond", "d-SEAMS")) or (
        keep("cna") and not seen("cna", "mdapy"))
    if need_diamond and not seen("diamond", "d-SEAMS"):
        summary, samples, ref_diamond = engine.diamond(frame)
        if keep("diamond"):
            emit(writer, handle, log, row_of(
                threads, cpus, frame, "diamond", "d-SEAMS", engine.version,
                summary, 1.0, "reference diamond common-neighbour labels", samples))
    if keep("diamond") and not seen("diamond", "d-SEAMS"):
        pass
    if keep("diamond") and ref_diamond is not None and not seen("diamond", "mdapy"):
        for code, caller in (
            ("mdapy", lambda: time_mdapy(frame, "diamond", threads)),
            ("pyscal3", lambda: time_pyscal(frame, "diamond", threads)),
        ):
            if seen("diamond", code):
                continue
            try:
                summary, samples, labels, version = caller()
                agree, detail = label_agreement(ref_diamond, labels, "dseams", "diamond")
                emit(writer, handle, log, row_of(
                    threads, cpus, frame, "diamond", code, version,
                    summary, agree, detail, samples))
            except Exception as exc:
                record_failure(writer, handle, log, threads, cpus, frame,
                               "diamond", code, str(exc))
    if keep("cna") and ref_diamond is not None and not seen("cna", "mdapy"):
        for code, caller in (
            ("mdapy", lambda: time_mdapy(frame, "cna", threads)),
            ("pyscal3", lambda: time_pyscal(frame, "cna", threads)),
        ):
            if seen("cna", code):
                continue
            try:
                summary, samples, labels, version = caller()
                agree, detail = label_agreement(ref_diamond, labels, "dseams", "cna")
                detail = (
                    "metal CNA versus diamond labels; fcc/hcp/bcc/ico are not "
                    "ice classes. " + detail
                )
                emit(writer, handle, log, row_of(
                    threads, cpus, frame, "cna", code, version,
                    summary, agree, detail, samples))
            except Exception as exc:
                record_failure(writer, handle, log, threads, cpus, frame,
                               "cna", code, str(exc))
    if keep("rings") and not seen("rings", "mdapy"):
        degree = median_degree(frame)
        # A cutoff list with dozens of neighbours has millions of rings of
        # size at most six. One d-SEAMS census is kept; the other codes are
        # not asked to repeat it.
        crowded = degree > 16.0
        if not seen("rings", "d-SEAMS"):
            summary, samples, ref_rings = engine.rings(frame)
            census = [int(ref_rings[k]) for k in (3, 4, 5, 6)]
            emit(writer, handle, log, row_of(
                threads, cpus, frame, "rings", "d-SEAMS", engine.version,
                summary, 1.0, f"reference primitive sizes 3-6 {census}", samples))
        else:
            census = None
        if crowded:
            note = (
                f"not timed; median coordination inside 3.5 A is {degree:.0f}. "
                "The primitive search returns millions of rings, which is not "
                "a water hydrogen-bond graph"
            )
            for code in ("mdapy", "SQQ"):
                if not seen("rings", code):
                    record_failure(writer, handle, log, threads, cpus, frame,
                                   "rings", code, note)
        else:
            if census is None:
                _summary, _samples, ref_rings = engine.rings(frame)
                census = [int(ref_rings[k]) for k in (3, 4, 5, 6)]
            if not seen("rings", "mdapy"):
                try:
                    summary, samples, counts, version = time_mdapy(frame, "rings", threads)
                    got = [int(counts[k]) for k in (3, 4, 5, 6)] if len(counts) > 6 else []
                    agree = float(got == census)
                    emit(writer, handle, log, row_of(
                        threads, cpus, frame, "rings", "mdapy", version, summary, agree,
                        f"shortest-path sizes 3-6 {got}", samples))
                except Exception as exc:
                    record_failure(writer, handle, log, threads, cpus, frame,
                                   "rings", "mdapy", str(exc))
            if not seen("rings", "SQQ"):
                try:
                    summary, samples, counts, version = time_sqq(frame, threads)
                    expect = census[1:]
                    agree = float(list(counts) == expect)
                    emit(writer, handle, log, row_of(
                        threads, cpus, frame, "rings", "SQQ", version, summary, agree,
                        f"shortest-path chordless 4-6 {list(counts)}; d-SEAMS {expect}; "
                        "SQQ does not search 3-rings", samples))
                except Exception as exc:
                    record_failure(writer, handle, log, threads, cpus, frame,
                                   "rings", "SQQ", str(exc))
    if keep("symmetry") and not seen("symmetry", "d-SEAMS"):
        summary, samples, ref_sym = engine.symmetry(frame)
        emit(writer, handle, log, row_of(
            threads, cpus, frame, "symmetry", "d-SEAMS",
            f"{engine.version}/{engine.backend}", summary, 1.0,
            "reference local symmetry classify", samples))
        try:
            samples, labels, haus, version = ira_match(frame, engine.templates(), threads)
            summary = proto.summarize(samples)
            ref = np.asarray(ref_sym)
            got = np.asarray(labels)
            agree = float(np.mean(ref == got)) if len(ref) else 1.0
            finite = haus[np.isfinite(haus)]
            median_h = float(np.median(finite)) if finite.size else float("nan")
            emit(writer, handle, log, row_of(
                threads, cpus, frame, "symmetry", "IRA", version, summary, agree,
                f"best ice template by Hausdorff; reject above {HAUSDORFF_A} A; "
                f"median_hausdorff={median_h:.4f}", samples))
        except Exception as exc:
            record_failure(writer, handle, log, threads, cpus, frame, "symmetry", "IRA", str(exc))


def launch(args, threads):
    cores = "0" if threads == 1 else f"0-{threads - 1}"
    env = os.environ.copy()
    env["H2H_PINNED"] = "1"
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = str(threads)
    cmd = [
        "taskset", "-c", cores, sys.executable, str(Path(__file__).resolve()),
        "--threads", str(threads),
        "--out", args.out, "--log", args.log,
        "--sizes", ",".join(str(size) for size in args.sizes),
    ]
    if args.lib:
        cmd.extend(["--lib", args.lib])
    if args.only:
        cmd.extend(["--only", args.only])
    if args.resume:
        cmd.append("--resume")
    if args.phase:
        cmd.extend(["--phase", args.phase])
    if args.molecules:
        cmd.extend(["--molecules", str(args.molecules)])
    subprocess.check_call(cmd, env=env)


def run_once(args):
    code, operation = args.once
    threads = args.threads or 1
    if os.environ.get("H2H_PINNED") != "1":
        pinned = argparse.Namespace(**vars(args))
        pinned.threads = threads
        pinned.once = None
        # A single pinned call. The session driver is the timed sample.
        cores = "0" if threads == 1 else f"0-{threads - 1}"
        env = os.environ.copy()
        env["H2H_PINNED"] = "1"
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_DYNAMIC"] = "FALSE"
        env["RAYON_NUM_THREADS"] = str(threads)
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        cmd = [
            "taskset", "-c", cores, sys.executable, str(Path(__file__).resolve()),
            "--once", code, operation, "--threads", str(threads),
            "--molecules", str(args.molecules or 4096),
            "--sizes", str(args.molecules or 4096),
            "--phase", args.phase or "Ic",
        ]
        if args.lib:
            cmd.extend(["--lib", args.lib])
        subprocess.check_call(cmd, env=env)
        return
    lib_path = args.lib or os.environ.get("SEAMS_CORE_LIB") or str(
        ROOT / "build" / "src" / "libyodaLib.so")
    frames = ice.campaign_frames((args.molecules or 4096,), ROOT)
    phase = args.phase or "Ic"
    frame = next(item for item in frames if item["phase"] == phase)
    if code == "d-SEAMS":
        engine = Engine(lib_path)
        if operation == "neighbours":
            engine.neighbours(frame)
        elif operation in ("q4", "q6"):
            engine.ql(frame, int(operation[1]))
        elif operation == "diamond":
            engine.diamond(frame)
        elif operation == "rings":
            engine.rings(frame)
        elif operation == "symmetry":
            engine.symmetry(frame)
        else:
            raise SystemExit(f"unknown operation {operation}")
        return
    if code == "mdapy":
        time_mdapy(frame, operation, threads)
    elif code == "pyscal3":
        time_pyscal(frame, operation, threads)
    elif code == "vesin":
        time_vesin(frame, threads)
    elif code == "SQQ":
        time_sqq(frame, threads)
    elif code == "IRA":
        engine = Engine(lib_path)
        ira_match(frame, engine.templates(), threads)
    else:
        raise SystemExit(f"unknown code {code}")


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="benchmarks/h2h-literature.csv")
    parser.add_argument("--log", default="benchmarks/h2h-literature.txt")
    parser.add_argument("--lib", default="")
    parser.add_argument("--sizes", default="128,4096,110592")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--only", default="")
    parser.add_argument("--phase", default="")
    parser.add_argument("--molecules", type=int, default=0)
    parser.add_argument("--once", nargs=2, metavar=("CODE", "OP"))
    parser.add_argument("--worker", default="")
    parser.add_argument("--payload", default="")
    parser.add_argument("--result", default="")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    args.sizes = tuple(int(part) for part in str(args.sizes).split(",") if part)
    return args


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.worker == "ira":
        ira_worker(args.payload, args.result)
        return
    if args.once:
        run_once(args)
        return
    if args.threads == 0:
        if not args.resume:
            for path in (args.out, args.log):
                target = Path(path)
                if target.exists():
                    target.unlink()
        for threads in (1, 8):
            launch(args, threads)
        return
    if os.environ.get("H2H_PINNED") != "1":
        launch(args, args.threads)
        return
    run_session(args, args.threads)


if __name__ == "__main__":
    main()
