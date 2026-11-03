#!/usr/bin/env python3
"""Repeated, equal-hardware comparisons of water-structure classifiers.

The campaign has four commands:

``prepare``
    Build one-frame oxygen trajectories for ideal ice Ic, ideal ice Ih,
    coexistence liquid water, and a published confined-water nanotube.
``case SYSTEM SCALE``
    Run all available implementations on one prepared frame. Each method is
    isolated in a subprocess and is killed at the stated resource limit.
``worker METHOD NPZ RESULT LABELS``
    Internal Python-method worker used by ``case``.
``summarize``
    Write timing intervals, tests, confusion matrices, and completion status.

The in-memory timed scope includes each code's own neighbour search. Whole
process ring rows use separate drivers and are never mixed with those rows.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import importlib.metadata
import json
import math
import os
import pathlib
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import time
from collections import Counter
from typing import Callable

import numpy as np

CUTOFF_A = 3.5
PTM_RMSD = 0.1
Q6_THRESHOLD = 0.5579227209091187
BOND_A = 2.75
WARMUP = int(os.environ.get("H2H_WARMUP", "5"))
REPEATS = int(os.environ.get("H2H_REPEATS", "30"))
CORE = int(os.environ.get("H2H_CORE", "3"))
WORKER_MEMORY_GIB = 12
NS_PER_MS = 1.0e6
MS_PER_S = 1000.0
# Every driver reads frame DRIVER_FRAME of its trajectory file.
DRIVER_FRAME = 1
# A whole-process ring run of the d-SEAMS driver makes one timed call and no
# warm-up call, so that the process time covers a single census.
PROCESS_WARMUP = 0
PROCESS_REPEATS = 1
# Slow codes run one warm-up call, then as many timed calls as fit in
# EXPENSIVE_TARGET_S seconds, at least EXPENSIVE_MIN_REPEATS and at most REPEATS.
# A warm-up call shorter than EXPENSIVE_FLOOR_MS counts as that long.
EXPENSIVE_TARGET_S = 120.0
EXPENSIVE_MIN_REPEATS = 6
EXPENSIVE_FLOOR_MS = 1.0
# RINGS runs once as warm-up and then RINGS_REPEATS times, one process per run.
RINGS_WARMUP = 1
RINGS_REPEATS = 6
# Every ring census counts the rings of MIN_RING to MAX_RING members.
MIN_RING = 3
MAX_RING = 6
RING_SIZES = range(MIN_RING, MAX_RING + 1)
# The campaign frames come from frame SOURCE_FRAME of each deposited dump.
SOURCE_FRAME = 1
# Atom type of the oxygen atoms in the confined-water dump of figshare 11448768.
CONFINED_OXYGEN_TYPE = 2
# Lines of a LAMMPS dump frame after its "ITEM: TIMESTEP" line: the atom count,
# the first of three box bounds, the atom column names and the first atom.
DUMP_COUNT_LINE = 3
DUMP_BOX_LINE = 5
DUMP_COLUMNS_LINE = 8
DUMP_FIRST_ATOM_LINE = 9
# Every system is prepared at every scale. A scale repeats the ideal unit cell
# the same number of times along each axis, the liquid frame likewise, and the
# confined frame by its own count along each axis.
SYSTEMS = ("ic", "ih", "liquid", "confined")
SCALES = ("1e3", "1e4", "1e5", "1e6")
IDEAL_REPEATS = {"1e3": 5, "1e4": 11, "1e5": 23, "1e6": 50}
LIQUID_REPEATS = {"1e3": 1, "1e4": 2, "1e5": 3, "1e6": 6}
CONFINED_REPEATS = {
    "1e3": (2, 2, 1),
    "1e4": (4, 5, 2),
    "1e5": (10, 10, 4),
    "1e6": (20, 20, 10),
}
LABELS = (
    "other",
    "cubic",
    "hexagonal",
    "interfacial",
    "hydrate",
    "interfacial-hydrate",
    "mixed",
)
LABEL_ID = {name: index for index, name in enumerate(LABELS)}
ICE_LIQUID_LABELS = ("liquid", "ice")
# The CHILL+ types of OVITO and of the Palmer library in the order of their
# numbers, as campaign labels. The Palmer library numbers liquid 4.
OVITO_CHILL_TYPES = (
    "other",
    "hexagonal",
    "cubic",
    "interfacial",
    "hydrate",
    "interfacial-hydrate",
)
PALMER_CHILL_TYPES = (
    "other",
    "cubic",
    "hexagonal",
    "interfacial",
    "other",
    "hydrate",
    "interfacial-hydrate",
)
OVITO_CHILL_LABELS = tuple(LABEL_ID[name] for name in OVITO_CHILL_TYPES)
PALMER_CHILL_LABELS = tuple(LABEL_ID[name] for name in PALMER_CHILL_TYPES)
# OVITO and pyscal number the diamond types alike: cubic diamond with its first
# and second neighbour shells, then the same three for hexagonal diamond.
DIAMOND_CUBIC_TYPES = (1, 2, 3)
DIAMOND_HEXAGONAL_TYPES = (4, 5, 6)
# q6 is the Steinhardt parameter of order Q6_ORDER over the Q6_NEIGHBOURS
# nearest neighbours.
Q6_ORDER = 6
Q6_NEIGHBOURS = 4
# The Palmer library computes CHILL+ from the harmonics of order CHILL_ORDER over
# the CHILL_NEIGHBOURS nearest neighbours. Its buffers hold PALMER_MAX_NEIGHBOURS
# neighbours per molecule, PALMER_HISTOGRAM_BINS bins of the bond correlations
# and PALMER_PLM_COLUMNS Legendre values per order.
CHILL_ORDER = 3
CHILL_NEIGHBOURS = 4
PALMER_MAX_NEIGHBOURS = 128
PALMER_HISTOGRAM_BINS = 20
PALMER_PLM_COLUMNS = 4
# A ratio of medians is resampled BOOTSTRAP_DRAWS times by a generator seeded
# with BOOTSTRAP_SEED, so that its interval is the same on every run.
BOOTSTRAP_SEED = 20261103
BOOTSTRAP_DRAWS = 10000
# Intervals of medians cover at least COVERAGE, and the interval of a ratio of
# medians lies between the INTERVAL_QUANTILES of its bootstrap ratios.
COVERAGE = 0.95
INTERVAL_QUANTILES = (0.025, 0.975)
# A timing difference needs a Holm-adjusted p-value below SIGNIFICANCE.
SIGNIFICANCE = 0.05
# A failed process reports the last STDERR_TAIL characters of its stderr.
STDERR_TAIL = 1000
# Files are hashed HASH_CHUNK_BYTES at a time.
HASH_CHUNK_BYTES = 1 << 20
# Default limits in seconds of a case, a ring case and a reference run, and of
# one run of RINGS; the listings pass their own.
DEFAULT_CASE_TIMEOUT_S = 900.0
DEFAULT_RING_CASE_TIMEOUT_S = 300.0
DEFAULT_REFERENCE_TIMEOUT_S = 900.0
DEFAULT_RINGS_RUN_TIMEOUT_S = 600.0

# Only these pairs implement the same scientific workload. The remaining
# classifiers are still timed and cross-tabulated, but a timing ratio between
# a local classifier and topological cage affiliation is not a winner claim.
EQUIVALENT_PAIRS = {
    ("dseams-chill", "ovito-chill"),
    ("dseams-chill", "palmer-chill"),
    ("dseams-q6", "freud-q6"),
    ("dseams-q6", "pyscal-q6"),
    ("dseams-rings", "matscipy-rings-public"),
    ("dseams-rings", "matscipy-rings-kernel"),
    ("dseams-rings-process", "rings-1.3.5"),
}
METHOD_INTERFACE = {
    "dseams-chill": "d-SEAMS test driver",
    "dseams-seeded": "d-SEAMS test driver",
    "dseams-q6": "d-SEAMS test driver",
    "dseams-rings": "d-SEAMS test driver",
    "ovito-chill": "public Python API",
    "ovito-ptm": "public Python API",
    "ovito-diamond": "public Python API",
    "ovito-cna": "public Python API",
    "palmer-chill": "adapted reference shared library",
    "freud-q6": "public Python API",
    "pyscal-q6": "public Python API",
    "pyscal-diamond": "public Python API",
    "matscipy-rings-public": "public Python API",
    "matscipy-rings-kernel": "documented internal FFI",
    "genice-rings": "pairlist and cycless APIs used by GenIce",
}


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_lammps(
    path: pathlib.Path, frame: int = 1, atom_type: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text().splitlines()
    starts = [
        index
        for index, line in enumerate(lines)
        if line.startswith("ITEM: TIMESTEP")
    ]
    start = starts[frame - 1]
    count = int(lines[start + DUMP_COUNT_LINE])
    low = np.array(
        [float(lines[start + DUMP_BOX_LINE + axis].split()[0]) for axis in range(3)]
    )
    high = np.array(
        [float(lines[start + DUMP_BOX_LINE + axis].split()[1]) for axis in range(3)]
    )
    columns = lines[start + DUMP_COLUMNS_LINE].split()[2:]
    xyz_columns = [columns.index(axis) for axis in ("x", "y", "z")]
    type_column = columns.index("type") if "type" in columns else None
    rows = [
        lines[start + DUMP_FIRST_ATOM_LINE + index].split() for index in range(count)
    ]
    if atom_type is not None and type_column is None:
        raise ValueError(f"{path} has no atom type column")
    positions = []
    for row in rows:
        if atom_type is None or int(row[type_column]) == atom_type:
            positions.append([float(row[index]) for index in xyz_columns])
    box = high - low
    return (np.asarray(positions, dtype=np.float64) - low) % box, box


def write_lammps(path: pathlib.Path, positions: np.ndarray, box: np.ndarray) -> None:
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        handle.write(f"{len(positions)}\nITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write(f"0.0 {length:.12f}\n")
        handle.write("ITEM: ATOMS id type x y z\n")
        for index, point in enumerate(positions, 1):
            handle.write(
                f"{index} 1 {point[0]:.12f} {point[1]:.12f} "
                f"{point[2]:.12f}\n"
            )


def cubic_diamond(repetitions: int) -> tuple[np.ndarray, np.ndarray]:
    lattice = 4.0 * BOND_A / math.sqrt(3.0)
    fcc = (
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    )
    basis = fcc + tuple(
        (x + 0.25, y + 0.25, z + 0.25) for x, y, z in fcc
    )
    cell = np.asarray(
        [
            ((i + x) * lattice, (j + y) * lattice, (k + z) * lattice)
            for i in range(repetitions)
            for j in range(repetitions)
            for k in range(repetitions)
            for x, y, z in basis
        ],
        dtype=np.float64,
    )
    box = np.full(3, repetitions * lattice)
    return cell % box, box


def lonsdaleite(repetitions: int) -> tuple[np.ndarray, np.ndarray]:
    lattice = BOND_A / ((3.0 / 8.0) * math.sqrt(8.0 / 3.0))
    height = lattice * math.sqrt(8.0 / 3.0)
    hexagonal = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    x_length = lattice
    y_length = lattice * math.sqrt(3.0)
    basis = []
    for hx, hy, hz in hexagonal:
        x = (hx - 0.5 * hy) * lattice
        y = hy * math.sqrt(3.0) * lattice / 2.0
        for shift_x, shift_y in (
            (0.0, 0.0),
            (x_length / 2.0, y_length / 2.0),
        ):
            basis.append(
                ((x + shift_x) % x_length, (y + shift_y) % y_length, hz * height)
            )
    cell = np.asarray(
        [
            (i * x_length + x, j * y_length + y, k * height + z)
            for i in range(repetitions)
            for j in range(repetitions)
            for k in range(repetitions)
            for x, y, z in basis
        ],
        dtype=np.float64,
    )
    box = np.array(
        [
            repetitions * x_length,
            repetitions * y_length,
            repetitions * height,
        ]
    )
    return cell % box, box


def replicate(
    positions: np.ndarray, box: np.ndarray, repeats: tuple[int, int, int]
) -> tuple[np.ndarray, np.ndarray]:
    copies = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                copies.append(positions + np.array((i, j, k)) * box)
    return np.concatenate(copies), box * np.asarray(repeats)


def prepare(args: argparse.Namespace) -> int:
    output = args.data
    output.mkdir(parents=True, exist_ok=True)
    liquid, liquid_box = read_lammps(args.liquid, frame=SOURCE_FRAME)
    confined, confined_box = read_lammps(
        args.confined, frame=SOURCE_FRAME, atom_type=CONFINED_OXYGEN_TYPE
    )
    manifest = {
        "created": utc_now(),
        "liquid_source": {
            "name": args.liquid.name,
            "sha256": sha256(args.liquid),
            "frame": SOURCE_FRAME,
        },
        "confined_source": {
            "name": args.confined.name,
            "sha256": sha256(args.confined),
            "doi": "10.6084/m9.figshare.11448768.v1",
            "frame": SOURCE_FRAME,
            "oxygen_type": CONFINED_OXYGEN_TYPE,
        },
        "frames": [],
    }
    for system in SYSTEMS:
        for scale in SCALES:
            if system == "ic":
                positions, box = cubic_diamond(IDEAL_REPEATS[scale])
            elif system == "ih":
                positions, box = lonsdaleite(IDEAL_REPEATS[scale])
            elif system == "liquid":
                rep = LIQUID_REPEATS[scale]
                positions, box = replicate(liquid, liquid_box, (rep, rep, rep))
            else:
                positions, box = replicate(
                    confined, confined_box, CONFINED_REPEATS[scale]
                )
            stem = f"{system}-{scale}"
            npz = output / f"{stem}.npz"
            trajectory = output / f"{stem}.lammpstrj"
            np.savez(npz, positions=positions, box=box)
            write_lammps(trajectory, positions, box)
            manifest["frames"].append(
                {
                    "system": system,
                    "scale": scale,
                    "molecules": len(positions),
                    "box_A": box.tolist(),
                    "npz": npz.name,
                    "trajectory": trajectory.name,
                    "npz_sha256": sha256(npz),
                    "trajectory_sha256": sha256(trajectory),
                }
            )
            print(
                f"prepared system={system} scale={scale} n={len(positions)}",
                flush=True,
            )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return 0


def measure(
    operation: Callable[[], np.ndarray | tuple | list],
    prepare_call: Callable[[], None] | None = None,
) -> tuple[list[float], object]:
    result = None
    samples = []
    for repetition in range(WARMUP + REPEATS):
        if prepare_call is not None:
            prepare_call()
        start = time.perf_counter_ns()
        result = operation()
        elapsed = (time.perf_counter_ns() - start) / NS_PER_MS
        if repetition >= WARMUP:
            samples.append(elapsed)
    return samples, result


def measure_expensive(
    operation: Callable[[], object], target_seconds: float = EXPENSIVE_TARGET_S
) -> tuple[list[float], object]:
    start = time.perf_counter_ns()
    result = operation()
    warmup_ms = (time.perf_counter_ns() - start) / NS_PER_MS
    budget = target_seconds * MS_PER_S / max(warmup_ms, EXPENSIVE_FLOOR_MS)
    repeats = max(EXPENSIVE_MIN_REPEATS, min(REPEATS, int(budget)))
    samples = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        result = operation()
        samples.append((time.perf_counter_ns() - start) / NS_PER_MS)
    return samples, result


def ovito_data(positions: np.ndarray, box: np.ndarray):
    from ovito.data import DataCollection

    data = DataCollection()
    data.create_cell(
        matrix=[
            [box[0], 0.0, 0.0, 0.0],
            [0.0, box[1], 0.0, 0.0],
            [0.0, 0.0, box[2], 0.0],
        ],
        pbc=(True, True, True),
    )
    particles = data.create_particles(count=len(positions))
    particles.create_property("Position", data=positions)
    return data


def worker_ovito(
    method: str, positions: np.ndarray, box: np.ndarray
) -> tuple[list[float], np.ndarray]:
    holder = {"data": None}

    def prep() -> None:
        holder["data"] = ovito_data(positions, box)

    if method == "ovito-chill":
        from ovito.modifiers import ChillPlusModifier

        modifier = ChillPlusModifier(cutoff=CUTOFF_A)

        def operation():
            holder["data"].apply(modifier)
            native = np.asarray(holder["data"].particles["Structure Type"])
            mapping = np.array(OVITO_CHILL_LABELS, dtype=np.int8)
            return mapping[native]

    elif method == "ovito-ptm":
        from ovito.modifiers import PolyhedralTemplateMatchingModifier as PTM

        modifier = PTM(rmsd_cutoff=PTM_RMSD)
        modifier.structures[PTM.Type.CUBIC_DIAMOND].enabled = True
        modifier.structures[PTM.Type.HEX_DIAMOND].enabled = True

        def operation():
            holder["data"].apply(modifier)
            native = np.asarray(holder["data"].particles["Structure Type"])
            labels = np.zeros(len(native), dtype=np.int8)
            labels[native == int(PTM.Type.CUBIC_DIAMOND)] = LABEL_ID["cubic"]
            labels[native == int(PTM.Type.HEX_DIAMOND)] = LABEL_ID["hexagonal"]
            return labels

    elif method == "ovito-diamond":
        from ovito.modifiers import IdentifyDiamondModifier

        modifier = IdentifyDiamondModifier()

        def operation():
            holder["data"].apply(modifier)
            native = np.asarray(holder["data"].particles["Structure Type"])
            labels = np.zeros(len(native), dtype=np.int8)
            labels[np.isin(native, DIAMOND_CUBIC_TYPES)] = LABEL_ID["cubic"]
            labels[np.isin(native, DIAMOND_HEXAGONAL_TYPES)] = LABEL_ID["hexagonal"]
            return labels

    elif method == "ovito-cna":
        from ovito.modifiers import CommonNeighborAnalysisModifier as CNA

        modifier = CNA(mode=CNA.Mode.AdaptiveCutoff)

        def operation():
            holder["data"].apply(modifier)
            # Adaptive CNA has no cubic- or hexagonal-diamond type. Preserve
            # its output as "other" in this ice-polymorph cross-tabulation.
            return np.zeros(len(positions), dtype=np.int8)

    else:
        raise ValueError(method)
    return measure(operation, prep)


def worker_freud(
    positions: np.ndarray, box: np.ndarray
) -> tuple[list[float], np.ndarray]:
    import freud

    freud.parallel.set_num_threads(1)
    fbox = freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2])
    shifted = positions.astype(np.float32) - box.astype(np.float32) / 2.0

    def operation():
        q6 = freud.order.Steinhardt(l=Q6_ORDER, average=False)
        q6.compute((fbox, shifted), neighbors={"num_neighbors": Q6_NEIGHBOURS})
        return (
            np.asarray(q6.particle_order) > Q6_THRESHOLD
        ).astype(np.int8)

    return measure(operation)


def ase_atoms(positions: np.ndarray, box: np.ndarray):
    from ase import Atoms

    return Atoms(
        symbols=["O"] * len(positions),
        positions=positions,
        cell=box,
        pbc=True,
    )


def worker_pyscal(
    method: str, positions: np.ndarray, box: np.ndarray
) -> tuple[list[float], np.ndarray]:
    import pyscal3

    pyscal3.set_num_threads(1)
    holder = {"atoms": None}

    def prep() -> None:
        holder["atoms"] = ase_atoms(positions, box)

    if method == "pyscal-q6":

        def operation():
            atoms = holder["atoms"]
            pyscal3.find_neighbors(
                atoms,
                method="number",
                nmax=Q6_NEIGHBOURS,
                assign_neighbor=True,
                store_rows=False,
            )
            values = pyscal3.steinhardt_parameter(atoms, Q6_ORDER, averaged=False)[0]
            return (np.asarray(values) > Q6_THRESHOLD).astype(np.int8)

    elif method == "pyscal-diamond":

        def operation():
            atoms = holder["atoms"]
            pyscal3.diamond_structure(atoms)
            native = atoms.arrays["pyscal_structure"]
            labels = np.zeros(len(native), dtype=np.int8)
            labels[np.isin(native, DIAMOND_CUBIC_TYPES)] = LABEL_ID["cubic"]
            labels[np.isin(native, DIAMOND_HEXAGONAL_TYPES)] = LABEL_ID["hexagonal"]
            return labels

    else:
        raise ValueError(method)
    return measure(operation, prep)


def worker_matscipy(
    method: str, positions: np.ndarray, box: np.ndarray
) -> tuple[list[float], dict[int, int]]:
    from matscipy.ffi import distances_on_graph, find_sp_rings
    from matscipy.neighbours import neighbour_list
    from matscipy.rings import ring_statistics

    atoms = ase_atoms(positions, box)

    def operation():
        if method == "matscipy-rings-public":
            # matscipy's C kernel treats maxlength as an exclusive bound.
            # Version 1.3.1's public wrapper invokes the kernel twice; this
            # row intentionally measures that released public interface.
            return ring_statistics(atoms, CUTOFF_A, maxlength=MAX_RING + 1)
        first, second, vectors = neighbour_list("ijD", atoms, CUTOFF_A)
        distances = distances_on_graph(first, second)
        return find_sp_rings(first, second, vectors, distances, MAX_RING + 1)

    samples, values = measure_expensive(operation)
    census = {
        size: int(values[size]) if size < len(values) else 0 for size in RING_SIZES
    }
    return samples, census


def worker_genice(
    positions: np.ndarray, box: np.ndarray
) -> tuple[list[float], dict[int, int]]:
    import networkx as nx
    import pairlist
    from cycless.cycles import cycles_iter

    relative = positions / box

    def operation():
        graph = nx.Graph()
        graph.add_nodes_from(range(len(positions)))
        graph.add_edges_from(
            pairlist.pairs_iter(relative, CUTOFF_A, np.diag(box), distance=False)
        )
        return list(cycles_iter(graph, MAX_RING, pos=relative))

    samples, rings = measure_expensive(operation)
    counts = Counter(len(ring) for ring in rings)
    return samples, {size: counts.get(size, 0) for size in RING_SIZES}


def worker_palmer(
    positions: np.ndarray, box: np.ndarray, library: pathlib.Path
) -> tuple[list[float], np.ndarray]:
    library_handle = ctypes.CDLL(str(library))
    order = ctypes.c_int(CHILL_ORDER)
    sph = np.zeros(2 * order.value + 1, dtype=np.float64)
    plm = np.zeros((order.value + 1) * PALMER_PLM_COLUMNS, dtype=np.float64)
    library_handle.initialize_CHILL(
        ctypes.byref(order),
        sph.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        plm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    key_text = b"CHILL+"
    key = (ctypes.c_char * len(key_text))(*key_text)
    key_length = ctypes.c_int(len(key_text))
    atom_count = ctypes.c_int(len(positions))
    maximum_neighbours = ctypes.c_int(PALMER_MAX_NEIGHBOURS)
    nearest_neighbours = ctypes.c_int(CHILL_NEIGHBOURS)
    cutoff_squared = ctypes.c_double(CUTOFF_A * CUTOFF_A)
    xyz = np.ascontiguousarray(positions, dtype=np.float64)
    lengths = np.ascontiguousarray(box, dtype=np.float64)

    def operation():
        correlation = np.zeros(
            nearest_neighbours.value * atom_count.value, dtype=np.float64
        )
        histogram = np.zeros(PALMER_HISTOGRAM_BINS, dtype=np.float64)
        native = np.zeros(atom_count.value, dtype=np.int32)
        ice_count = ctypes.c_int()
        library_handle.call_CHILL(
            sph.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            plm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            key,
            ctypes.byref(key_length),
            ctypes.byref(atom_count),
            ctypes.byref(maximum_neighbours),
            ctypes.byref(nearest_neighbours),
            ctypes.byref(order),
            ctypes.byref(cutoff_squared),
            lengths.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            correlation.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            histogram.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            native.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(ice_count),
        )
        mapping = np.array(PALMER_CHILL_LABELS, dtype=np.int8)
        return mapping[native]

    return measure_expensive(operation)


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unrecorded"


def worker(args: argparse.Namespace) -> int:
    memory_limit = WORKER_MEMORY_GIB * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (memory_limit, memory_limit))
    payload = np.load(args.npz)
    positions = payload["positions"]
    box = payload["box"]
    started = utc_now()
    if args.method.startswith("ovito-"):
        samples, result = worker_ovito(args.method, positions, box)
    elif args.method == "freud-q6":
        samples, result = worker_freud(positions, box)
    elif args.method.startswith("pyscal-"):
        samples, result = worker_pyscal(args.method, positions, box)
    elif args.method.startswith("matscipy-rings-"):
        samples, result = worker_matscipy(args.method, positions, box)
    elif args.method == "genice-rings":
        samples, result = worker_genice(positions, box)
    elif args.method == "palmer-chill":
        samples, result = worker_palmer(positions, box, args.palmer)
    else:
        raise ValueError(args.method)
    record = {
        "method": args.method,
        "warmup": (
            1
            if args.method.startswith("matscipy-rings-")
            or args.method in ("genice-rings", "palmer-chill")
            else WARMUP
        ),
        "repeats": len(samples),
        "samples_ms": samples,
        "started": started,
        "finished": utc_now(),
        "versions": {
            "ovito": package_version("ovito"),
            "freud-analysis": package_version("freud-analysis"),
            "pyscal3": package_version("pyscal3"),
            "matscipy": package_version("matscipy"),
            "GenIce2": package_version("GenIce2"),
            "cycless": package_version("cycless"),
        },
    }
    if isinstance(result, dict):
        record["census"] = {str(key): value for key, value in result.items()}
    else:
        np.save(args.labels, np.asarray(result, dtype=np.int8))
        record["labels"] = args.labels.name
        record["histogram"] = {
            str(key): value
            for key, value in sorted(
                Counter(np.asarray(result).tolist()).items()
            )
        }
    args.result.write_text(json.dumps(record, indent=2) + "\n")
    return 0


def parse_driver_samples(stdout: str) -> list[float]:
    samples = []
    columns_seen = False
    for line in stdout.splitlines():
        if line.startswith("columns "):
            columns_seen = True
            continue
        if columns_seen:
            try:
                samples.append(float(line))
            except ValueError:
                continue
    return samples


def dseams_method(
    method: str,
    trajectory: pathlib.Path,
    bench: pathlib.Path,
    work: pathlib.Path,
    timeout: float,
) -> tuple[dict, pathlib.Path | None]:
    modes = {
        "dseams-chill": "chill",
        "dseams-q6": "q6-full",
        "dseams-seeded": "seeded-adj",
        "dseams-rings": "rings",
    }
    sidecar = work / f"{method}.txt"
    command = [
        "taskset",
        "-c",
        str(CORE),
        str(bench),
        str(trajectory),
        str(DRIVER_FRAME),
        str(WARMUP),
        str(REPEATS),
        modes[method],
        str(sidecar),
    ]
    child_environment = os.environ.copy()
    child_environment["LD_LIBRARY_PATH"] = ":".join(
        entry
        for entry in child_environment.get("LD_LIBRARY_PATH", "").split(":")
        if entry
        and not entry.startswith("/lib/")
        and not entry.startswith("/usr/lib/")
    )
    started = utc_now()
    process = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=child_environment,
    )
    if process.returncode:
        raise RuntimeError(
            f"{method} exited {process.returncode}: {process.stderr[-STDERR_TAIL:]}"
        )
    samples = parse_driver_samples(process.stdout)
    if len(samples) != REPEATS:
        raise RuntimeError(f"{method} returned {len(samples)} samples")
    record = {
        "method": method,
        "warmup": WARMUP,
        "repeats": REPEATS,
        "samples_ms": samples,
        "started": started,
        "finished": utc_now(),
        "command": [pathlib.Path(item).name if "/" in item else item for item in command],
    }
    labels_path = None
    if method in ("dseams-chill", "dseams-seeded"):
        native = sidecar.read_text().splitlines()
        mapping = {
            "cubic": "cubic",
            "reCubic": "cubic",
            "hexagonal": "hexagonal",
            "reHex": "hexagonal",
            "interfacial": "interfacial",
            "clathrate": "hydrate",
            "interClathrate": "interfacial-hydrate",
            "mixed": "mixed",
            "water": "other",
            "unclassified": "other",
        }
        labels = np.asarray([LABEL_ID[mapping[item]] for item in native], dtype=np.int8)
        labels_path = work / f"{method}.npy"
        np.save(labels_path, labels)
        record["labels"] = labels_path.name
        record["histogram"] = {
            str(key): value
            for key, value in sorted(Counter(labels.tolist()).items())
        }
    elif method == "dseams-q6":
        values = []
        value_mode = False
        for line in sidecar.read_text().splitlines():
            if line == "values":
                value_mode = True
                continue
            if value_mode:
                values.append(float(line.split()[0]))
        labels = (np.asarray(values) > Q6_THRESHOLD).astype(np.int8)
        labels_path = work / f"{method}.npy"
        np.save(labels_path, labels)
        record["labels"] = labels_path.name
        record["histogram"] = {
            str(key): value
            for key, value in sorted(Counter(labels.tolist()).items())
        }
    else:
        record["census"] = {
            line.split()[0]: int(line.split()[1])
            for line in sidecar.read_text().splitlines()
        }
    return record, labels_path


def append_jsonl(path: pathlib.Path, record: dict) -> None:
    with path.open("a") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def write_xyz(path: pathlib.Path, positions: np.ndarray, box: np.ndarray) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(positions)}\n")
        handle.write(f"box {box[0]:.12f} {box[1]:.12f} {box[2]:.12f}\n")
        for point in positions:
            handle.write(f"O {point[0]:.12f} {point[1]:.12f} {point[2]:.12f}\n")


def prepare_rings_directory(
    directory: pathlib.Path, positions: np.ndarray, box: np.ndarray
) -> None:
    from timing_protocol import RINGS_OPTIONS

    if directory.exists():
        shutil.rmtree(directory)
    (directory / "data").mkdir(parents=True)
    write_xyz(directory / "data" / "conf.xyz", positions, box)
    (directory / "options").write_text(RINGS_OPTIONS)
    # The grid settings after conf.xyz are those of timing_protocol.py, which
    # explains them; the search depth after them limits rings to twice its value.
    (directory / "input").write_text(
        "#######################################\n"
        "       R.I.N.G.S. input file          #\n"
        "#######################################\n"
        f"h2h\n{len(positions)}\n1\nO\n1\n0\n"
        f"{box[0]:.12f} {box[1]:.12f} {box[2]:.12f}\n"
        "90.0 90.0 90.0\n1.0\nANI\nconf.xyz\n"
        f"200\n500\n25.0\n0.125\n90\n20\n{MAX_RING // 2}\n10\n"
        "#######################################\n"
        f"O O {CUTOFF_A:.3f}\nGrtot {CUTOFF_A:.3f}\n"
        "#######################################\n"
    )


def clean_rings_outputs(directory: pathlib.Path) -> None:
    for name in ("rstat", "tmp", "gr", "sq", "cstat", "bonds", "angles"):
        shutil.rmtree(directory / name, ignore_errors=True)


def rings_census(directory: pathlib.Path) -> dict[str, int]:
    path = directory / "rstat" / "evol-RINGS-5.dat"
    census: dict[str, int] = {}
    seen = False
    for line in path.read_text(errors="replace").splitlines():
        if "Total number of rings" in line:
            seen = True
            continue
        fields = line.split()
        if seen and len(fields) >= 2:
            try:
                census[str(int(float(fields[0])))] = int(float(fields[1]))
            except ValueError:
                continue
    return census


def external_rings(args: argparse.Namespace) -> int:
    os.sched_setaffinity(0, {CORE})
    stem = f"{args.system}-{args.scale}"
    payload = np.load(args.data / f"{stem}.npz")
    positions = payload["positions"]
    box = payload["box"]
    trajectory = args.data / f"{stem}.lammpstrj"
    work = args.output / "work" / f"rings-{stem}"
    prepare_rings_directory(work, positions, box)
    child_environment = os.environ.copy()
    child_environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
        }
    )
    child_environment["LD_LIBRARY_PATH"] = ":".join(
        entry
        for entry in child_environment.get("LD_LIBRARY_PATH", "").split(":")
        if entry
        and not entry.startswith("/lib/")
        and not entry.startswith("/usr/lib/")
    )
    methods = (
        (
            "dseams-rings-process",
            [
                "taskset",
                "-c",
                str(CORE),
                str(args.bench),
                str(trajectory),
                str(DRIVER_FRAME),
                str(PROCESS_WARMUP),
                str(PROCESS_REPEATS),
                "rings",
                str(work / "dseams-rings-process.txt"),
            ],
            None,
            WARMUP,
            REPEATS,
        ),
        (
            "rings-1.3.5",
            ["taskset", "-c", str(CORE), str(args.rings), "input"],
            work,
            RINGS_WARMUP,
            RINGS_REPEATS,
        ),
    )
    for method, command, cwd, warmup, repeats in methods:
        samples = []
        last = None
        try:
            for repetition in range(warmup + repeats):
                if method == "rings-1.3.5":
                    clean_rings_outputs(work)
                start = time.perf_counter_ns()
                last = subprocess.run(
                    command,
                    cwd=cwd,
                    env=child_environment,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                )
                elapsed = (time.perf_counter_ns() - start) / NS_PER_MS
                if last.returncode:
                    raise RuntimeError(
                        f"exit {last.returncode}: {last.stderr[-STDERR_TAIL:]}"
                    )
                if repetition >= warmup:
                    samples.append(elapsed)
            census = (
                rings_census(work)
                if method == "rings-1.3.5"
                else {
                    line.split()[0]: int(line.split()[1])
                    for line in (
                        work / "dseams-rings-process.txt"
                    ).read_text().splitlines()
                }
            )
            append_jsonl(
                args.output / "records.jsonl",
                {
                    "system": args.system,
                    "scale": args.scale,
                    "molecules": len(positions),
                    "operation": "rings-process",
                    "method": method,
                    "warmup": warmup,
                    "repeats": repeats,
                    "samples_ms": samples,
                    "census": census,
                    "scope": "whole-process",
                    "core": CORE,
                    "source_sha256": sha256(args.data / f"{stem}.npz"),
                },
            )
            print(f"complete {stem} {method}", flush=True)
        except subprocess.TimeoutExpired:
            append_jsonl(
                args.output / "errors.jsonl",
                {
                    "system": args.system,
                    "scale": args.scale,
                    "molecules": len(positions),
                    "operation": "rings-process",
                    "method": method,
                    "status": "timeout",
                    "timeout_s": args.timeout,
                },
            )
            print(f"timeout {stem} {method}", flush=True)
        except Exception as error:
            append_jsonl(
                args.output / "errors.jsonl",
                {
                    "system": args.system,
                    "scale": args.scale,
                    "molecules": len(positions),
                    "operation": "rings-process",
                    "method": method,
                    "status": "error",
                    "error": f"{type(error).__name__}: {error}",
                },
            )
            print(f"error {stem} {method}: {error}", flush=True)
    return 0


def truth_label(system: str, operation: str) -> int:
    if operation == "ice-liquid":
        return int(system in ("ic", "ih"))
    if system == "ic":
        return LABEL_ID["cubic"]
    if system == "ih":
        return LABEL_ID["hexagonal"]
    return LABEL_ID["other"]


def confusion(
    left: np.ndarray, right: np.ndarray, size: int
) -> list[list[int]]:
    matrix = np.zeros((size, size), dtype=np.int64)
    np.add.at(matrix, (left, right), 1)
    return matrix.tolist()


def case(args: argparse.Namespace) -> int:
    os.sched_setaffinity(0, {CORE})
    os.environ.update(
        {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "true",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
            "OVITO_THREAD_COUNT": "1",
        }
    )
    stem = f"{args.system}-{args.scale}"
    npz = args.data / f"{stem}.npz"
    trajectory = args.data / f"{stem}.lammpstrj"
    payload = np.load(npz)
    count = len(payload["positions"])
    work = args.output / "work" / stem
    work.mkdir(parents=True, exist_ok=True)
    records_path = args.output / "records.jsonl"
    errors_path = args.output / "errors.jsonl"
    confusions_path = args.output / "confusions.jsonl"
    methods = {
        "chill": ("dseams-chill", "ovito-chill", "palmer-chill"),
        "ice-polymorph": (
            "dseams-seeded",
            "ovito-ptm",
            "ovito-diamond",
            "ovito-cna",
            "pyscal-diamond",
        ),
        "ice-liquid": ("dseams-q6", "freud-q6", "pyscal-q6"),
        "rings": (
            "dseams-rings",
            "matscipy-rings-public",
            "matscipy-rings-kernel",
            "genice-rings",
        ),
    }
    labels: dict[str, np.ndarray] = {}
    operation_of = {
        method: operation
        for operation, operation_methods in methods.items()
        for method in operation_methods
    }
    for method in [item for group in methods.values() for item in group]:
        base = {
            "system": args.system,
            "scale": args.scale,
            "molecules": count,
            "operation": operation_of[method],
            "source_sha256": sha256(npz),
            "core": CORE,
            "scope": "in-memory, including neighbour search",
            "interface": METHOD_INTERFACE[method],
        }
        try:
            if method.startswith("dseams-"):
                record, labels_path = dseams_method(
                    method, trajectory, args.bench, work, args.timeout
                )
                if labels_path is not None:
                    labels[method] = np.load(labels_path)
            else:
                result = work / f"{method}.json"
                label_path = work / f"{method}.npy"
                command = [
                    sys.executable,
                    str(pathlib.Path(__file__).resolve()),
                    "worker",
                    method,
                    str(npz),
                    str(result),
                    str(label_path),
                    "--palmer",
                    str(args.palmer),
                ]
                process = subprocess.run(
                    command,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=(
                        args.ring_timeout
                        if method.startswith("matscipy-rings-")
                        or method == "genice-rings"
                        else args.reference_timeout
                        if method == "palmer-chill"
                        else args.timeout
                    ),
                    env=os.environ,
                )
                if process.returncode:
                    raise RuntimeError(
                        f"{method} exited {process.returncode}: "
                        f"{process.stderr[-STDERR_TAIL:]}"
                    )
                record = json.loads(result.read_text())
                if label_path.exists():
                    labels[method] = np.load(label_path)
            record.update(base)
            append_jsonl(records_path, record)
            print(
                f"complete system={args.system} scale={args.scale} "
                f"method={method}",
                flush=True,
            )
        except subprocess.TimeoutExpired:
            append_jsonl(
                errors_path,
                {
                    **base,
                    "method": method,
                    "status": "timeout",
                    "timeout_s": (
                        args.ring_timeout
                        if method.startswith("matscipy-rings-")
                        or method == "genice-rings"
                        else args.reference_timeout
                        if method == "palmer-chill"
                        else args.timeout
                    ),
                },
            )
            print(
                f"timeout system={args.system} scale={args.scale} method={method}",
                flush=True,
            )
        except Exception as error:
            append_jsonl(
                errors_path,
                {
                    **base,
                    "method": method,
                    "status": "error",
                    "error": f"{type(error).__name__}: {error}",
                },
            )
            print(
                f"error system={args.system} scale={args.scale} method={method}: "
                f"{error}",
                flush=True,
            )
    for operation, operation_methods in methods.items():
        available = [method for method in operation_methods if method in labels]
        if not available:
            continue
        names = ICE_LIQUID_LABELS if operation == "ice-liquid" else LABELS
        if args.system != "confined":
            expected = np.full(
                count, truth_label(args.system, operation), dtype=np.int8
            )
            for method in available:
                append_jsonl(
                    confusions_path,
                    {
                        "system": args.system,
                        "scale": args.scale,
                        "molecules": count,
                        "operation": operation,
                        "labels": names,
                        "left": "truth",
                        "right": method,
                        "matrix": confusion(
                            expected, labels[method], len(names)
                        ),
                    },
                )
        reference = available[0]
        for method in available[1:]:
            append_jsonl(
                confusions_path,
                {
                    "system": args.system,
                    "scale": args.scale,
                    "molecules": count,
                    "operation": operation,
                    "labels": names,
                    "left": reference,
                    "right": method,
                    "matrix": confusion(
                        labels[reference], labels[method], len(names)
                    ),
                },
            )
    return 0


def median_interval(samples: list[float]) -> tuple[float, float]:
    ordered = sorted(samples)
    count = len(ordered)
    selected = 1
    for index in range(1, count // 2 + 1):
        lower_tail = sum(
            math.comb(count, value) for value in range(index)
        ) / (2**count)
        coverage = 1.0 - 2.0 * lower_tail
        if coverage >= COVERAGE:
            selected = index
        else:
            break
    return ordered[selected - 1], ordered[count - selected]


def bootstrap_ratio(
    numerator: list[float],
    denominator: list[float],
    draws: int = BOOTSTRAP_DRAWS,
) -> tuple[float, float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    top = np.asarray(numerator)
    bottom = np.asarray(denominator)
    ratios = np.empty(draws)
    for draw in range(draws):
        ratios[draw] = np.median(
            rng.choice(top, len(top), replace=True)
        ) / np.median(rng.choice(bottom, len(bottom), replace=True))
    low, high = np.quantile(ratios, INTERVAL_QUANTILES)
    return float(low), float(high)


def holm_adjust(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    adjusted = np.empty(len(pvalues))
    running = 0.0
    total = len(pvalues)
    for rank, index in enumerate(order):
        value = min(1.0, (total - rank) * pvalues[index])
        running = max(running, value)
        adjusted[index] = running
    return adjusted.tolist()


def summarize(args: argparse.Namespace) -> int:
    from scipy.stats import fisher_exact, mannwhitneyu

    records = [
        json.loads(line)
        for line in (args.output / "records.jsonl").read_text().splitlines()
        if line
    ]
    confusions = [
        json.loads(line)
        for line in (args.output / "confusions.jsonl").read_text().splitlines()
        if line
    ]
    errors_path = args.output / "errors.jsonl"
    errors = (
        [
            json.loads(line)
            for line in errors_path.read_text().splitlines()
            if line
        ]
        if errors_path.exists()
        else []
    )
    rows = []
    for record in records:
        samples = record["samples_ms"]
        low, high = median_interval(samples)
        rows.append(
            {
                "system": record["system"],
                "scale": record["scale"],
                "molecules": record["molecules"],
                "operation": record["operation"],
                "method": record["method"],
                "repeats": len(samples),
                "median_ms": statistics.median(samples),
                "ci95_low_ms": low,
                "ci95_high_ms": high,
            }
        )
    comparisons = []
    for row in rows:
        if row["method"].startswith("dseams-"):
            continue
        reference = next(
            (
                candidate
                for candidate in rows
                if candidate["system"] == row["system"]
                and candidate["scale"] == row["scale"]
                and candidate["operation"] == row["operation"]
                and candidate["method"].startswith("dseams-")
            ),
            None,
        )
        if reference is None:
            continue
        left = next(
            record["samples_ms"]
            for record in records
            if record["system"] == row["system"]
            and record["scale"] == row["scale"]
            and record["method"] == reference["method"]
        )
        right = next(
            record["samples_ms"]
            for record in records
            if record["system"] == row["system"]
            and record["scale"] == row["scale"]
            and record["method"] == row["method"]
        )
        ratio = statistics.median(right) / statistics.median(left)
        ratio_low, ratio_high = bootstrap_ratio(right, left)
        test = mannwhitneyu(left, right, alternative="two-sided", method="auto")
        comparisons.append(
            {
                "system": row["system"],
                "scale": row["scale"],
                "molecules": row["molecules"],
                "operation": row["operation"],
                "dseams": reference["method"],
                "other": row["method"],
                "other_over_dseams": ratio,
                "ratio_ci95_low": ratio_low,
                "ratio_ci95_high": ratio_high,
                "mann_whitney_u": float(test.statistic),
                "p_raw": float(test.pvalue),
            }
        )
    adjusted = holm_adjust([item["p_raw"] for item in comparisons])
    for item, value in zip(comparisons, adjusted):
        item["p_holm"] = value
        item["timing_difference_supported"] = bool(
            value < SIGNIFICANCE
            and (
                item["ratio_ci95_low"] > 1.0
                or item["ratio_ci95_high"] < 1.0
            )
        )
        item["equivalent_workload"] = (
            item["dseams"], item["other"]
        ) in EQUIVALENT_PAIRS
        item["ordering_supported"] = bool(
            item["equivalent_workload"]
            and item["timing_difference_supported"]
        )
    accuracy = []
    for item in confusions:
        if item["left"] != "truth":
            continue
        matrix = np.asarray(item["matrix"], dtype=np.int64)
        correct = int(np.trace(matrix))
        total = int(matrix.sum())
        accuracy.append(
            {
                "system": item["system"],
                "scale": item["scale"],
                "molecules": item["molecules"],
                "operation": item["operation"],
                "method": item["right"],
                "correct": correct,
                "incorrect": total - correct,
                "accuracy": correct / total,
            }
        )
    fisher_rows = []
    for row in accuracy:
        if row["method"].startswith("dseams-"):
            continue
        reference = next(
            (
                candidate
                for candidate in accuracy
                if candidate["system"] == row["system"]
                and candidate["scale"] == row["scale"]
                and candidate["operation"] == row["operation"]
                and candidate["method"].startswith("dseams-")
            ),
            None,
        )
        if reference is None:
            continue
        table = [
            [reference["correct"], reference["incorrect"]],
            [row["correct"], row["incorrect"]],
        ]
        odds, pvalue = fisher_exact(table)
        fisher_rows.append(
            {
                "system": row["system"],
                "scale": row["scale"],
                "operation": row["operation"],
                "dseams": reference["method"],
                "other": row["method"],
                "odds_ratio": (
                    float(odds) if math.isfinite(float(odds)) else None
                ),
                "p_two_sided": (
                    float(pvalue)
                    if math.isfinite(float(pvalue))
                    else None
                ),
            }
        )
    fisher_pvalues = [
        row["p_two_sided"] for row in fisher_rows
        if row["p_two_sided"] is not None
    ]
    fisher_adjusted = iter(holm_adjust(fisher_pvalues))
    for row in fisher_rows:
        row["p_holm"] = (
            next(fisher_adjusted)
            if row["p_two_sided"] is not None
            else None
        )
    pairwise_agreement = []
    for item in confusions:
        if item["left"] == "truth":
            continue
        matrix = np.asarray(item["matrix"], dtype=np.int64)
        total = int(matrix.sum())
        agreed = int(np.trace(matrix))
        pairwise_agreement.append(
            {
                "system": item["system"],
                "scale": item["scale"],
                "molecules": item["molecules"],
                "operation": item["operation"],
                "left": item["left"],
                "right": item["right"],
                "agreed": agreed,
                "disagreed": total - agreed,
                "agreement": agreed / total,
            }
        )
    ring_census = []
    for record in records:
        if "census" not in record:
            continue
        ring_census.append(
            {
                "system": record["system"],
                "scale": record["scale"],
                "molecules": record["molecules"],
                "operation": record["operation"],
                "method": record["method"],
                **{
                    f"rings_{size}": int(record["census"].get(str(size), 0))
                    for size in RING_SIZES
                },
            }
        )
    census_comparisons = []
    for row in ring_census:
        if row["method"].startswith("dseams-"):
            continue
        reference = next(
            (
                candidate
                for candidate in ring_census
                if candidate["system"] == row["system"]
                and candidate["scale"] == row["scale"]
                and candidate["operation"] == row["operation"]
                and candidate["method"].startswith("dseams-")
            ),
            None,
        )
        if reference is None:
            continue
        exact = all(
            row[f"rings_{size}"] == reference[f"rings_{size}"]
            for size in RING_SIZES
        )
        census_comparisons.append(
            {
                "system": row["system"],
                "scale": row["scale"],
                "operation": row["operation"],
                "dseams": reference["method"],
                "other": row["method"],
                "exact_census_3_to_6": exact,
            }
        )
    payload = {
        "created": utc_now(),
        "timings": rows,
        "timing_comparisons": comparisons,
        "accuracy": accuracy,
        "fisher_exact": fisher_rows,
        "confusions": confusions,
        "pairwise_agreement": pairwise_agreement,
        "ring_census": ring_census,
        "census_comparisons": census_comparisons,
        "errors": errors,
        "statistics": {
            "median_interval": "exact binomial order-statistic, at least 95%",
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "multiple_testing": "Holm correction over timing comparisons",
        },
    }
    (args.output / "summary.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    with (args.output / "timings.tsv").open("w") as handle:
        if rows:
            keys = list(rows[0])
            handle.write("\t".join(keys) + "\n")
            for row in rows:
                handle.write("\t".join(str(row[key]) for key in keys) + "\n")
    with (args.output / "comparisons.tsv").open("w") as handle:
        if comparisons:
            keys = list(comparisons[0])
            handle.write("\t".join(keys) + "\n")
            for row in comparisons:
                handle.write("\t".join(str(row[key]) for key in keys) + "\n")
    with (args.output / "agreements.tsv").open("w") as handle:
        if pairwise_agreement:
            keys = list(pairwise_agreement[0])
            handle.write("\t".join(keys) + "\n")
            for row in pairwise_agreement:
                handle.write("\t".join(str(row[key]) for key in keys) + "\n")
    with (args.output / "ring-census.tsv").open("w") as handle:
        if ring_census:
            keys = list(ring_census[0])
            handle.write("\t".join(keys) + "\n")
            for row in ring_census:
                handle.write("\t".join(str(row[key]) for key in keys) + "\n")
    markdown = [
        "# Head-to-head campaign summary",
        "",
        "Medians and exact order-statistic intervals are in milliseconds.",
        "An ordering is reported only for equivalent workloads.",
        "",
        "| system | scale | molecules | operation | method | n | median | 95% CI |",
        "|---|---:|---:|---|---|---:|---:|---:|",
    ]
    for row in rows:
        markdown.append(
            f"| {row['system']} | {row['scale']} | {row['molecules']} | "
            f"{row['operation']} | {row['method']} | {row['repeats']} | "
            f"{row['median_ms']:.6g} | "
            f"[{row['ci95_low_ms']:.6g}, {row['ci95_high_ms']:.6g}] |"
        )
    markdown.extend(
        [
            "",
            "## Censored or failed rows",
            "",
            "| system | scale | operation | method | status | detail |",
            "|---|---:|---|---|---|---|",
        ]
    )
    for error in errors:
        detail = (
            f"{error.get('timeout_s')} s"
            if error.get("status") == "timeout"
            else str(error.get("error", ""))
        ).replace("|", "\\|")
        markdown.append(
            f"| {error['system']} | {error['scale']} | "
            f"{error['operation']} | {error['method']} | "
            f"{error['status']} | {detail} |"
        )
    (args.output / "summary.md").write_text("\n".join(markdown) + "\n")
    return 0


def environment(args: argparse.Namespace) -> int:
    def version(command: list[str]) -> str:
        return subprocess.run(
            command, check=False, capture_output=True, text=True
        ).stdout.strip()

    details = {
        "created": utc_now(),
        "uname": platform.uname()._asdict(),
        "python": platform.python_version(),
        "affinity": sorted(os.sched_getaffinity(0)),
        "lscpu": subprocess.run(
            ["lscpu"], check=False, capture_output=True, text=True
        ).stdout,
        "memtotal": next(
            line
            for line in pathlib.Path("/proc/meminfo").read_text().splitlines()
            if line.startswith("MemTotal:")
        ),
        "versions": {
            name: package_version(name)
            for name in (
                "numpy",
                "scipy",
                "ovito",
                "freud-analysis",
                "pyscal3",
                "matscipy",
                "GenIce2",
                "cycless",
            )
        },
        "toolchain": {
            "c++": version(["c++", "--version"]),
            "gfortran": version(["gfortran", "--version"]),
            "meson": version(["meson", "--version"]),
        },
        "thread_controls": {
            "core": CORE,
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
            "OVITO_THREAD_COUNT": "1",
            "frequency_pinned": False,
        },
        "bench_sha256": sha256(args.bench),
        "library_sha256": sha256(args.library),
        "source_commit": subprocess.run(
            ["git", "-C", str(args.source), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "source_tree": subprocess.run(
            ["git", "-C", str(args.source), "write-tree"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "source_dirty": bool(
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(args.source),
                    "status",
                    "--porcelain",
                    "--untracked-files=no",
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        ),
        "palmer_commit": subprocess.run(
            ["git", "-C", str(args.palmer_source), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "palmer_library_sha256": sha256(args.palmer_library),
        "rings_archive_sha256": sha256(args.rings_archive),
        "rings_sha256": sha256(args.rings),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "environment.json").write_text(
        json.dumps(details, indent=2) + "\n"
    )
    return 0


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__)
    subcommands = command.add_subparsers(dest="command", required=True)
    prep = subcommands.add_parser("prepare")
    prep.add_argument("--data", type=pathlib.Path, required=True)
    prep.add_argument("--liquid", type=pathlib.Path, required=True)
    prep.add_argument("--confined", type=pathlib.Path, required=True)
    prep.set_defaults(function=prepare)

    run_case = subcommands.add_parser("case")
    run_case.add_argument("system", choices=SYSTEMS)
    run_case.add_argument("scale", choices=SCALES)
    run_case.add_argument("--data", type=pathlib.Path, required=True)
    run_case.add_argument("--output", type=pathlib.Path, required=True)
    run_case.add_argument("--bench", type=pathlib.Path, required=True)
    run_case.add_argument("--palmer", type=pathlib.Path, required=True)
    run_case.add_argument("--timeout", type=float, default=DEFAULT_CASE_TIMEOUT_S)
    run_case.add_argument(
        "--ring-timeout", type=float, default=DEFAULT_RING_CASE_TIMEOUT_S
    )
    run_case.add_argument(
        "--reference-timeout", type=float, default=DEFAULT_REFERENCE_TIMEOUT_S
    )
    run_case.set_defaults(function=case)

    rings = subcommands.add_parser("rings")
    rings.add_argument("system", choices=SYSTEMS)
    rings.add_argument("scale", choices=SCALES)
    rings.add_argument("--data", type=pathlib.Path, required=True)
    rings.add_argument("--output", type=pathlib.Path, required=True)
    rings.add_argument("--bench", type=pathlib.Path, required=True)
    rings.add_argument("--rings", type=pathlib.Path, required=True)
    rings.add_argument("--timeout", type=float, default=DEFAULT_RINGS_RUN_TIMEOUT_S)
    rings.set_defaults(function=external_rings)

    work = subcommands.add_parser("worker")
    work.add_argument("method")
    work.add_argument("npz", type=pathlib.Path)
    work.add_argument("result", type=pathlib.Path)
    work.add_argument("labels", type=pathlib.Path)
    work.add_argument("--palmer", type=pathlib.Path, required=True)
    work.set_defaults(function=worker)

    summary = subcommands.add_parser("summarize")
    summary.add_argument("--output", type=pathlib.Path, required=True)
    summary.set_defaults(function=summarize)

    env = subcommands.add_parser("environment")
    env.add_argument("--output", type=pathlib.Path, required=True)
    env.add_argument("--bench", type=pathlib.Path, required=True)
    env.add_argument("--library", type=pathlib.Path, required=True)
    env.add_argument("--source", type=pathlib.Path, required=True)
    env.add_argument("--palmer-source", type=pathlib.Path, required=True)
    env.add_argument("--palmer-library", type=pathlib.Path, required=True)
    env.add_argument("--rings-archive", type=pathlib.Path, required=True)
    env.add_argument("--rings", type=pathlib.Path, required=True)
    env.set_defaults(function=environment)
    return command


if __name__ == "__main__":
    arguments = parser().parse_args()
    raise SystemExit(arguments.function(arguments))
