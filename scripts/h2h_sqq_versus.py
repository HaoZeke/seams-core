#!/usr/bin/env python3
"""d-SEAMS against SQQ on ring statistics and ice classification.

SQQ is Pang and Sun, J. Phys. Chem. A 130, 8187 (2026),
doi:10.1021/acs.jpca.6c05769, code github.com/pimooni/sqq. The PyPI
wheel 0.5.8 is tag v0.5.8, commit 17cf76f9a2fc18848321025674a0f4d4adfa26e8.
Only the pure-Python engine has rings and the ice stage; the native
engine does cages and F3/F4.

Rings. d-SEAMS `seams_rings_find` on the 3.5 A cutoff graph, sizes up
to 6, against SQQ `ring.definition: shortest_path`, sizes 4 to 6, on
the 0.35 nm O-O graph. Both predicates are Franzblau's shortest-path
test. SQQ does not search three-rings, so d-SEAMS does strictly more
work in the timed call. A chordless cycle is fixed by its vertex set,
so rings are compared as vertex sets.

Ice. d-SEAMS published CHILL+ (Nguyen and Molinero 2015, q3 over the
four nearest neighbours inside 3.5 A, inclusive c_ij edges) against
the SQQ ice stage (`ice.method: chill`): a water is ice-I when it has
four graph neighbours, every neighbour has four, and it sits on at
least two six-rings; four neighbours and one six-ring is interfacial.
SQQ's six-rings are its default chordless cycles. SQQ is asked for
size 6 only, the only size its rule reads.

Like for like. CHILL+ reads bond orientations and no rings, so the ice
row compares two whole tools. bench/ring_ice.c applies SQQ's rule to
d-SEAMS rings through the public C API: on the primitive six-rings of
`seams_rings_find` (ice-primitive), and on the chordless six-cycles of
the same d-SEAMS cutoff graph (ice-chordless), the ring set SQQ counts.
The three d-SEAMS ice calls and one SQQ ice call share each timing
round, so every d-SEAMS ice row is paired with the same SQQ samples.

Correctness writes one row per frame and quantity with the overlap and
every disagreement. Each disagreement is checked against reference
implementations in this file that call neither code: a periodic KD
tree graph, a breadth-first shortest-path test, the SQQ rule, and
CHILL+ from scipy spherical harmonics.

Timing is one process per thread count. One warm-up of each call,
then n paired rounds, the order reversing every round so neither code
always runs first. Frames are built before the clock. One thread is cpu 0; eight
threads are cpus 0-7. The cpu list and core count are logged, host
names are not.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
import platform
import random
import re
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import gen_ideal_ice as ice  # noqa: E402

CUTOFF_A = 3.5
MAX_RING = 6
SQQ_COMMIT = "17cf76f9a2fc18848321025674a0f4d4adfa26e8"
SQQ_DOI = "10.1021/acs.jpca.6c05769"
BOOT = 10000
BOOT_SEED = 20261114
FRAMES = (
    "Ic-4096", "Ih-4096", "mw-liquid-4096", "tip4p2005-4096",
    "Ic-110592", "Ih-110592", "mw-liquid-110592", "tip4p-o-250",
)
OPERATIONS = ("rings", "ice", "ice-primitive", "ice-chordless")
# The d-SEAMS ice calls, in the order of an even round. SQQ's ice call
# is the same for all three.
ICE_CALLS = ("ice", "ice-primitive", "ice-chordless")
RING_ICE_LABEL = {0: "liquid", 1: "interfacial", 2: "ice"}

# seams_label codes.
CUBIC, HEX, WATER, INTERFACIAL, CLATH, INTERCLATH = 0, 1, 2, 3, 4, 5
SEAMS_NAME = {0: "cubic", 1: "hexagonal", 2: "water", 3: "interfacial",
              4: "clathrate", 5: "interClathrate", 6: "unclassified",
              7: "reCubic", 8: "reHex"}
# Both rule sets reduce to ice-I, interfacial, and the rest.
SEAMS_MERGED = {0: "ice", 1: "ice", 3: "interfacial", 7: "ice", 8: "ice"}


def scrub(text):
    return re.sub(r"(?<![\w.])/[\w.+@=-]+(?:/[\w.+@=-]+)+", "<path>", str(text))


def cpus_allowed():
    cpus = sorted(os.sched_getaffinity(0))
    if not cpus:
        return ""
    parts, start, prev = [], cpus[0], cpus[0]
    for cpu in cpus[1:]:
        if cpu == prev + 1:
            prev = cpu
            continue
        parts.append(f"{start}-{prev}" if prev > start else f"{start}")
        start = prev = cpu
    parts.append(f"{start}-{prev}" if prev > start else f"{start}")
    return ",".join(parts)


def code_state():
    """Commit, and whether the code the timings depend on matches it."""
    paths = ["scripts", "bench", "src", "include", "Snakefile", "meson.build"]
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", *paths],
                             capture_output=True, text=True, check=True)
        clean = "clean" if not out.stdout.strip() else "modified"
    except Exception:
        clean = "unknown"
    return f"{git_head()} code {clean}"


def git_head():
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------- frames

def mw_dir():
    return Path(os.environ.get("MW215_DIR", "mw215-frames"))


def tip4p2005_dir():
    return Path(os.environ.get("TIP4P2005_DIR", "tip4p2005-frames"))


def wrapped(frame):
    """Shift the origin to zero and wrap into [0, L). Both codes get this."""
    box9 = np.asarray(frame["box9"], dtype=np.float64)
    if np.any(box9[6:9] != 0.0):
        raise ValueError(f"{frame['name']} is tilted; SQQ reads orthorhombic cells only")
    lengths = box9[3:6].copy()
    pos = np.asarray(frame["positions"], dtype=np.float64) - box9[0:3]
    pos = pos - np.floor(pos / lengths) * lengths
    pos = np.where(pos >= lengths, pos - lengths, pos)
    out = dict(frame)
    out["positions"] = np.ascontiguousarray(pos)
    out["lengths"] = lengths
    out["box9"] = np.array([0, 0, 0, *lengths, 0, 0, 0], dtype=np.float64)
    return out


def load_frame(name):
    if name.startswith(("Ic-", "Ih-")):
        n = int(name.split("-")[1])
        frame = ice.ice_ic(n) if name.startswith("Ic-") else ice.ice_ih(n)
        frame["source"] = "scripts/gen_ideal_ice.py"
    elif name.startswith("mw-liquid-"):
        n = int(name.rsplit("-", 1)[1])
        path = mw_dir() / f"mw-liquid-{n}.lammpstrj"
        if not path.is_file():
            raise FileNotFoundError(
                f"{path.name} is missing; run scripts/h2h_sqq_mw.py with lmp on PATH")
        frame = ice.read_lammps_first(path, 1, name, "liquid-mw")
        frame["source"] = f"{path.name} sha256 {sha256(path)}"
    elif name == "tip4p2005-4096":
        path = tip4p2005_dir() / "tip4p2005-4096.lammpstrj"
        if not path.is_file():
            raise FileNotFoundError(
                f"{path.name} is missing; run input/tip4p2005/regenerate.py with lmp on PATH")
        frame = ice.read_lammps_first(path, 1, name, "liquid-tip4p2005")
        got = sha256(path)
        want = (ROOT / "input" / "tip4p2005" / "tip4p2005-4096.sha256").read_text().split()[0]
        # LAMMPS with OpenMP is not bitwise reproducible across machines, so a
        # regenerated frame is recorded rather than refused.
        frame["source"] = f"{path.name} sha256 {got} " + (
            "committed" if got == want else f"differs from committed {want}")
    elif name == "tip4p-o-250":
        path = ROOT / "input" / "traj" / "exampleTraj.lammpstrj"
        frame = ice.read_lammps_first(path, 2, name, "tip4p-o")
        frame["source"] = "input/traj/exampleTraj.lammpstrj type 2"
    else:
        raise ValueError(name)
    return wrapped(frame)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------- d-SEAMS

class Seams:
    def __init__(self, path):
        self.path = path
        lib = ctypes.CDLL(path)
        dp = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        vp = ctypes.c_void_p
        lib.seams_version_string.restype = ctypes.c_char_p
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_set_num_threads.argtypes = [ctypes.c_int]
        lib.seams_set_num_threads.restype = ctypes.c_int
        lib.seams_get_num_threads.restype = ctypes.c_int
        lib.seams_system_create.restype = vp
        lib.seams_system_destroy.argtypes = [vp]
        lib.seams_system_set_atoms.argtypes = [
            vp, dp, ctypes.c_size_t, dp, ctypes.c_size_t, ip, ip]
        lib.seams_system_set_atoms.restype = ctypes.c_int
        lib.seams_rings_find.argtypes = [
            vp, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(vp)]
        lib.seams_rings_find.restype = ctypes.c_int
        lib.seams_rings_destroy.argtypes = [vp]
        lib.seams_rings_count.argtypes = [vp]
        lib.seams_rings_count.restype = ctypes.c_int
        lib.seams_rings_vertices.argtypes = [
            vp, ctypes.c_int, ip, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        lib.seams_rings_vertices.restype = ctypes.c_int
        lib.seams_chill_plus.argtypes = [
            vp, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int, ip, ip]
        lib.seams_chill_plus.restype = ctypes.c_int
        lib.seams_neighbours_build.argtypes = [
            vp, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.POINTER(vp)]
        lib.seams_neighbours_build.restype = ctypes.c_int
        lib.seams_neighbours_destroy.argtypes = [vp]
        lib.seams_neighbours_copy.argtypes = [
            vp, ip, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
            ip, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        lib.seams_neighbours_copy.restype = ctypes.c_int
        self.lib = lib
        self.version = lib.seams_version_string().decode()

    def check(self, code, what):
        if code != 0:
            raise RuntimeError(f"{what}: status {code}: {self.lib.seams_last_error().decode()}")

    def threads(self, n):
        self.check(self.lib.seams_set_num_threads(int(n)), "set_num_threads")

    def system(self, frame):
        xyz = np.ascontiguousarray(frame["positions"], dtype=np.float64)
        box9 = np.ascontiguousarray(frame["box9"], dtype=np.float64)
        sys_ = self.lib.seams_system_create()
        if not sys_:
            raise RuntimeError("seams_system_create returned NULL")
        dp = ctypes.POINTER(ctypes.c_double)
        self.check(self.lib.seams_system_set_atoms(
            sys_, xyz.ctypes.data_as(dp), xyz.shape[0], box9.ctypes.data_as(dp), 9,
            None, None), "set_atoms")
        return sys_

    def rings_handle(self, sys_):
        out = ctypes.c_void_p()
        self.check(self.lib.seams_rings_find(sys_, CUTOFF_A, 0, 0, MAX_RING,
                                             ctypes.byref(out)), "rings_find")
        return out

    def ring_sets(self, handle):
        total = self.lib.seams_rings_count(handle)
        buf = (ctypes.c_int * 16)()
        got = ctypes.c_size_t()
        by_size = defaultdict(set)
        for k in range(total):
            self.check(self.lib.seams_rings_vertices(handle, k, buf, 16, ctypes.byref(got)),
                       "rings_vertices")
            nodes = tuple(buf[i] for i in range(got.value))
            by_size[len(nodes)].add(frozenset(nodes))
        return by_size

    def chill_plus(self, sys_, labels):
        counts = (ctypes.c_int * 9)()
        self.check(self.lib.seams_chill_plus(
            sys_, CUTOFF_A, 0, 0, 0,
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), counts), "chill_plus")
        return list(counts)

    def edges(self, sys_):
        nb = ctypes.c_void_p()
        self.check(self.lib.seams_neighbours_build(sys_, CUTOFF_A, 0, 0, ctypes.byref(nb)),
                   "neighbours_build")
        try:
            n_off = ctypes.c_size_t()
            nnz = ctypes.c_size_t()
            self.check(self.lib.seams_neighbours_copy(nb, None, 0, ctypes.byref(n_off),
                                                      None, 0, ctypes.byref(nnz)), "copy")
            off = np.zeros(n_off.value, dtype=np.int32)
            cols = np.zeros(max(nnz.value, 1), dtype=np.int32)
            ip = ctypes.POINTER(ctypes.c_int)
            self.check(self.lib.seams_neighbours_copy(
                nb, off.ctypes.data_as(ip), off.size, ctypes.byref(n_off),
                cols.ctypes.data_as(ip), cols.size, ctypes.byref(nnz)), "copy")
        finally:
            self.lib.seams_neighbours_destroy(nb)
        edges = set()
        for i in range(off.size - 1):
            for j in cols[off[i]:off[i + 1]]:
                edges.add((min(i, int(j)), max(i, int(j))))
        return edges


class RingIce:
    """bench/ring_ice.c: SQQ's ice rule on d-SEAMS primitive or chordless six-rings."""

    DEFINITION = {"ice-primitive": 0, "ice-chordless": 1}

    def __init__(self, path, engine):
        lib = ctypes.CDLL(path)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.seams_h2h_ring_ice.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ip, ip,
            ctypes.c_size_t, ip, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        lib.seams_h2h_ring_ice.restype = ctypes.c_int
        self.lib = lib
        self.engine = engine
        maps = Path("/proc/self/maps")
        if maps.is_file():
            seen = {line.split()[-1] for line in maps.read_text().splitlines()
                    if line.rstrip().endswith("/libseams.so")}
            if len(seen) != 1:
                raise RuntimeError(f"{len(seen)} copies of libseams.so are mapped; need one")

    def run(self, sys_, operation, labels, six, rings=False):
        ip = ctypes.POINTER(ctypes.c_int)
        definition = self.DEFINITION[operation]
        got = ctypes.c_size_t()
        n = labels.size
        args = (sys_, CUTOFF_A, definition, 2, labels.ctypes.data_as(ip), six.ctypes.data_as(ip),
                n)
        if not rings:
            self.engine.check(self.lib.seams_h2h_ring_ice(*args, None, 0, None), operation)
            return None
        self.engine.check(self.lib.seams_h2h_ring_ice(*args, None, 0, ctypes.byref(got)),
                          operation)
        buf = np.zeros((max(got.value, 1), 6), dtype=np.int32)
        self.engine.check(self.lib.seams_h2h_ring_ice(
            *args, buf.ctypes.data_as(ip), got.value, ctypes.byref(got)), operation)
        return {frozenset(int(v) for v in row) for row in buf[:got.value]}


# ---------------------------------------------------------------- SQQ

def sqq_version():
    import sqq
    return sqq.__version__


def sqq_frame(frame):
    from sqq.models import Atom, Frame
    pos_nm = np.asarray(frame["positions"], dtype=float) / 10.0
    atoms = [Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
                  xyz=pos_nm[i], molecule_id=i + 1) for i in range(len(pos_nm))]
    return Frame(name=frame["name"], atoms=atoms,
                 box=np.asarray(frame["lengths"], dtype=float) / 10.0)


def sqq_config(operation, threads, sizes=None):
    import warnings
    from sqq.api import load_config
    if operation == "rings":
        ring = {"sizes": sizes or [4, 5, 6], "definition": "shortest_path", "chordless": True}
    else:
        ring = {"sizes": sizes or [6], "definition": "chordless", "chordless": True}
    raw = {
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": CUTOFF_A / 10.0},
        "ring": ring,
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": False},
        "order_parameter": {"parameters": []},
        "ice": {"enabled": operation == "ice", "method": "chill", "min_six_rings": 2,
                "require_four_coord_neighbors": True},
        "parallel": {"workers": threads, "math_threads": threads},
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cfg = load_config(raw, engine="py")
    d = cfg.to_mutable_dict()
    assert d["order"]["parameters"] == [], d["order"]
    assert not d["cage"]["enabled"] and not d["half_cage"]["enabled"]
    assert not d["quasi_cage"]["enabled"]
    assert d["ice"]["enabled"] == (operation == "ice")
    assert d["ring"]["definition"] == ring["definition"]
    return cfg


def sqq_run(built, cfg):
    from sqq.api import analyze_frame
    return analyze_frame(built, cfg, engine="py")


def sqq_ring_sets(result):
    out = defaultdict(set)
    for size, rings in result.rings.items():
        for ring in rings:
            out[int(size)].add(frozenset(int(v) for v in ring.nodes))
    return out


def sqq_ice_labels(result, n):
    labels = np.full(n, "liquid", dtype=object)
    for v in result.interfacial_ice_waters:
        labels[int(v)] = "interfacial"
    for v in result.ice_i_waters:
        labels[int(v)] = "ice"
    return labels


def sqq_edges(result):
    return {(min(int(a), int(b)), max(int(a), int(b))) for a, b in result.graph.edges}


# ---------------------------------------------------------------- references

def reference_graph(frame):
    """Periodic KD tree; pairs strictly inside the cutoff and their distances."""
    from scipy.spatial import cKDTree
    pos = frame["positions"]
    tree = cKDTree(pos, boxsize=frame["lengths"])
    pairs = tree.query_pairs(CUTOFF_A + 0.05, output_type="ndarray")
    d = pos[pairs[:, 0]] - pos[pairs[:, 1]]
    d -= np.round(d / frame["lengths"]) * frame["lengths"]
    r = np.sqrt((d * d).sum(1))
    keep = r < CUTOFF_A
    edges = {(int(min(a, b)), int(max(a, b))) for a, b in pairs[keep]}
    margin = float(np.min(np.abs(r - CUTOFF_A))) if r.size else float("inf")
    adj = [set() for _ in range(len(pos))]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return edges, adj, margin


def ring_cycle(nodes, adj):
    """Cyclic order of a vertex set that induces a cycle, or None."""
    nodes = set(nodes)
    start = min(nodes)
    order = [start]
    prev = None
    cur = start
    while True:
        nxt = sorted((adj[cur] & nodes) - ({prev} if prev is not None else set()))
        if prev is None:
            nxt = nxt[:1]
        if not nxt:
            return None
        step = nxt[0]
        if step == start:
            break
        if step in order:
            return None
        order.append(step)
        prev, cur = cur, step
        if len(order) > len(nodes):
            return None
    return order if len(order) == len(nodes) else None


def is_chordless(nodes, adj):
    nodes = set(nodes)
    return all(len(adj[v] & nodes) == 2 for v in nodes) and ring_cycle(nodes, adj) is not None


def reference_chordless6(adj):
    """Every simple 6-cycle from its lowest vertex, kept when is_chordless accepts it."""
    found = set()
    for s in range(len(adj)):
        for b in adj[s]:
            if b <= s:
                continue
            for c in adj[b]:
                if c <= s or c == b:
                    continue
                for d in adj[c]:
                    if d <= s or d in (b, c):
                        continue
                    for e in adj[d]:
                        if e <= s or e in (b, c, d):
                            continue
                        for f in adj[e]:
                            if f <= s or f in (b, c, d, e) or s not in adj[f]:
                                continue
                            ring = frozenset((s, b, c, d, e, f))
                            if ring not in found and is_chordless(ring, adj):
                                found.add(ring)
    return found


def bfs_dist(adj, src, limit):
    dist = {src: 0}
    queue = deque([src])
    while queue:
        u = queue.popleft()
        if dist[u] >= limit:
            continue
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                queue.append(w)
    return dist


def is_shortest_path_ring(nodes, adj):
    """Franzblau: every pair at graph distance equal to the ring distance."""
    order = ring_cycle(nodes, adj)
    if order is None or not is_chordless(nodes, adj):
        return False
    k = len(order)
    for a in range(k):
        dist = bfs_dist(adj, order[a], k // 2)
        for b in range(k):
            ring_d = min(abs(a - b), k - abs(a - b))
            if dist.get(order[b], k) != ring_d:
                return False
    return True


def shortcut_witness(nodes, adj):
    order = ring_cycle(nodes, adj)
    k = len(order)
    for a in range(k):
        dist = bfs_dist(adj, order[a], k // 2)
        for b in range(a + 1, k):
            ring_d = min(b - a, k - (b - a))
            got = dist.get(order[b], k)
            if got < ring_d:
                return order[a], order[b], ring_d, got
    return None


def sqq_rule(adj, six_rings, n):
    """SQQ ice.method chill: degree, six-ring count, and neighbour degree."""
    count = np.zeros(n, dtype=np.int64)
    for ring in six_rings:
        for v in ring:
            count[v] += 1
    deg = np.array([len(a) for a in adj])
    labels = np.full(n, "liquid", dtype=object)
    for i in range(n):
        if deg[i] != 4 or count[i] == 0:
            continue
        if count[i] >= 2 and all(deg[j] == 4 for j in adj[i]):
            labels[i] = "ice"
        else:
            labels[i] = "interfacial"
    return labels, deg, count


def global_q6(frame, edges):
    """Steinhardt Q6 of the whole bond set. Ideal ice Ic is 0.5. A liquid is near 0."""
    from scipy.special import sph_harm_y
    pairs = np.array(sorted(edges), dtype=np.int64)
    pos, L = frame["positions"], frame["lengths"]
    d = pos[pairs[:, 1]] - pos[pairs[:, 0]]
    d -= np.round(d / L) * L
    r = np.linalg.norm(d, axis=1)
    theta = np.arccos(np.clip(d[:, 2] / r, -1.0, 1.0))
    phi = np.arctan2(d[:, 1], d[:, 0])
    total = 0.0
    for m in range(-6, 7):
        total += abs(np.mean(sph_harm_y(6, m, theta, phi))) ** 2
    return math.sqrt(4.0 * math.pi / 13.0 * total)


def xy_span(frame):
    pos, L = frame["positions"], frame["lengths"]
    d = pos[:, :2] - pos[0, :2]
    d -= np.round(d / L[:2]) * L[:2]
    return d.max(axis=0) - d.min(axis=0)


def reference_chill_plus(frame):
    """CHILL+ from scipy Y_3m, the four nearest inside 3.5 A, inclusive edges."""
    from scipy.spatial import cKDTree
    from scipy.special import sph_harm_y
    pos = frame["positions"]
    L = frame["lengths"]
    n = len(pos)
    tree = cKDTree(pos, boxsize=L)
    dist, idx = tree.query(pos, k=12, distance_upper_bound=CUTOFF_A + 0.5)
    dist, idx = dist[:, 1:], idx[:, 1:]
    shell = np.sum(dist < CUTOFF_A, axis=1)
    keep = np.minimum(shell, 4)
    slot = np.arange(4)[None, :] < keep[:, None]
    rows, cols = np.nonzero(slot)
    nbr = idx[:, :4][rows, cols].astype(np.int64)
    d = pos[nbr] - pos[rows]
    d -= np.round(d / L) * L
    r = np.linalg.norm(d, axis=1)
    theta = np.arccos(np.clip(d[:, 2] / r, -1.0, 1.0))
    phi = np.arctan2(d[:, 1], d[:, 0])
    q = np.zeros((n, 7), dtype=np.complex128)
    for m in range(-3, 4):
        np.add.at(q[:, m + 3], rows, sph_harm_y(3, m, theta, phi))
    chosen = [[] for _ in range(n)]
    for i, j in zip(rows.tolist(), nbr.tolist()):
        chosen[i].append(j)
    norm = np.sqrt(np.real(np.sum(q * np.conj(q), axis=1)))
    cij = [[] for _ in range(n)]
    for i in range(n):
        for j in chosen[i]:
            denom = norm[i] * norm[j]
            c = float(np.real(np.sum(q[i] * np.conj(q[j]))) / denom) if denom > 0 else 0.0
            cij[i].append(c)
    stag = np.array([sum(c <= -0.8 for c in row) for row in cij])
    eclip = np.array([sum((-0.35 <= c <= 0.25) for c in row) for row in cij])
    labels = np.full(n, WATER, dtype=np.int32)
    for i in range(n):
        if shell[i] != 4 or keep[i] != 4:
            continue
        inter = False
        if stag[i] == 2:
            inter = any(stag[j] > 2 for j in chosen[i])
        if stag[i] == 3 and eclip[i] == 0:
            inter = any(stag[j] > 1 for j in chosen[i])
        if eclip[i] == 0 and stag[i] == 4:
            labels[i] = CUBIC
        elif eclip[i] == 1 and stag[i] == 3:
            labels[i] = HEX
        elif inter:
            labels[i] = INTERFACIAL
        elif eclip[i] == 4 and stag[i] == 0:
            labels[i] = CLATH
        elif eclip[i] == 3:
            labels[i] = INTERCLATH
    return labels, shell, stag, eclip, cij


# ---------------------------------------------------------------- correctness

def jaccard(a, b):
    """Two empty sets agree completely."""
    union = len(a | b)
    return len(a & b) / union if union else 1.0


def agree(args):
    engine = Seams(args.lib)
    engine.threads(1)
    census = RingIce(args.ring_ice_lib, engine)
    version = sqq_version()
    rows = []
    log = open(args.log, "w", encoding="utf-8")

    def say(text):
        log.write(scrub(text) + "\n")
        log.flush()
        print(text, flush=True)

    say(f"d-SEAMS {engine.version} commit {code_state()}")
    say(f"SQQ {version} commit {SQQ_COMMIT} doi {SQQ_DOI} engine py")
    say(f"cutoff_A {CUTOFF_A} cores {len(os.sched_getaffinity(0))} cpus {cpus_allowed()}")
    for name in args.frames:
        frame = load_frame(name)
        n = frame["molecules"]
        say(f"frame {name} molecules {n} edges_A {' '.join(f'{x:.6f}' for x in frame['lengths'])} "
            f"source {frame['source']}")
        sys_ = engine.system(frame)
        ref_edges, adj, margin = reference_graph(frame)
        seams_e = engine.edges(sys_)
        built = sqq_frame(frame)
        res_r = sqq_run(built, sqq_config("rings", 1))
        sqq_e = sqq_edges(res_r)
        say(f"graph {name} reference {len(ref_edges)} d-SEAMS {len(seams_e)} SQQ {len(sqq_e)} "
            f"d-SEAMS^ref {len(seams_e ^ ref_edges)} SQQ^ref {len(sqq_e ^ ref_edges)} "
            f"closest_pair_to_cutoff_A {margin:.3e}")
        rows.append(dict(frame=name, molecules=n, quantity="O-O edges", d_seams=len(seams_e),
                         sqq=len(sqq_e), both=len(seams_e & sqq_e),
                         d_seams_only=len(seams_e - sqq_e), sqq_only=len(sqq_e - seams_e),
                         agreement=f"{jaccard(seams_e, sqq_e):.6f}",
                         resolution=f"reference KD tree {len(ref_edges)} edges; nearest pair "
                                    f"{margin:.3e} A from the cutoff"))

        handle = engine.rings_handle(sys_)
        dsets = engine.ring_sets(handle)
        engine.lib.seams_rings_destroy(handle)
        ssets = sqq_ring_sets(res_r)
        say(f"rings {name} d-SEAMS 3-rings {len(dsets.get(3, ()))} (SQQ does not search size 3)")
        for size in (4, 5, 6):
            d, s = dsets.get(size, set()), ssets.get(size, set())
            only_d, only_s = d - s, s - d
            verdict = []
            for tag, group in (("d-SEAMS", only_d), ("SQQ", only_s)):
                good = sum(is_shortest_path_ring(r, adj) for r in group)
                verdict.append(f"{tag}-only {len(group)}, reference shortest-path {good}")
            say(f"rings {name} size {size} d-SEAMS {len(d)} SQQ {len(s)} both {len(d & s)} "
                + "; ".join(verdict))
            rows.append(dict(frame=name, molecules=n, quantity=f"{size}-rings", d_seams=len(d),
                             sqq=len(s), both=len(d & s), d_seams_only=len(only_d),
                             sqq_only=len(only_s),
                             agreement=f"{jaccard(d, s):.6f}",
                             resolution="; ".join(verdict)))
        sample = [r for size in (4, 5, 6) for r in sorted(dsets.get(size, ()), key=sorted)[:200]]
        bad = sum(not is_shortest_path_ring(r, adj) for r in sample)
        say(f"rings {name} reference check of {len(sample)} shared rings: {bad} fail")

        labels = np.zeros(n, dtype=np.int32)
        engine.chill_plus(sys_, labels)
        res_i = sqq_run(built, sqq_config("ice", 1))
        s_lab = sqq_ice_labels(res_i, n)
        res_i456 = sqq_run(built, sqq_config("ice", 1, sizes=[4, 5, 6]))
        same456 = bool(np.all(sqq_ice_labels(res_i456, n) == s_lab))
        chord6 = sqq_ring_sets(res_i).get(6, set())
        ref_lab, deg, six = sqq_rule(adj, chord6, n)
        prim_lab, _deg, six_p = sqq_rule(adj, dsets.get(6, set()), n)
        r_lab, shell, stag, eclip, cij = reference_chill_plus(frame)
        say(f"ice {name} SQQ labels sizes [6] == sizes [4,5,6]: {same456}")
        say(f"ice {name} SQQ rule re-derived from the reference graph and SQQ chordless 6-rings: "
            f"{int(np.sum(ref_lab == s_lab))}/{n} equal")
        say(f"ice {name} chordless 6-rings {len(chord6)} primitive 6-rings {len(dsets.get(6, ()))}; "
            f"SQQ rule on primitive 6-rings changes {int(np.sum(prim_lab != s_lab))} labels")
        say(f"ice {name} CHILL+ re-derived with scipy Y_3m: {int(np.sum(r_lab == labels))}/{n} equal")
        seams_m = np.array([SEAMS_MERGED.get(int(v), "liquid") for v in labels], dtype=object)
        conf = Counter(zip((SEAMS_NAME[int(v)] for v in labels), s_lab))
        say(f"ice {name} CHILL+ counts "
            + " ".join(f"{SEAMS_NAME[k]}={int(np.sum(labels == k))}" for k in range(6)))
        say(f"ice {name} SQQ counts "
            + " ".join(f"{k}={int(np.sum(s_lab == k))}" for k in ("ice", "interfacial", "liquid")))
        say(f"ice {name} confusion (CHILL+ label, SQQ label): "
            + "; ".join(f"{a}/{b} {c}" for (a, b), c in sorted(conf.items())))
        mismatch = np.flatnonzero(seams_m != s_lab)
        say(f"ice {name} merged agreement {n - len(mismatch)}/{n}")

        def sqq_clause(i):
            if deg[i] != 4:
                return f"degree {int(deg[i])}"
            if six[i] == 0:
                return "no six-ring"
            if six[i] == 1:
                return "one six-ring"
            if not all(deg[j] == 4 for j in adj[i]):
                return "two+ six-rings, a neighbour not 4-coordinated"
            return "two+ six-rings, neighbours 4-coordinated"

        def chill_clause(i):
            if shell[i] != 4:
                return f"{int(shell[i])} neighbours inside 3.5 A"
            return f"{int(stag[i])} staggered {int(eclip[i])} eclipsed"

        causes = Counter((SEAMS_NAME[int(labels[i])], s_lab[i], chill_clause(i), sqq_clause(i))
                         for i in mismatch)
        say(f"ice {name} disagreements by rule clause (CHILL+ label | SQQ label | CHILL+ clause "
            f"| SQQ clause): count")
        shown = 0
        ranked = sorted(causes.items(), key=lambda kv: (-kv[1], kv[0]))
        for (a, b, c, d), count in ranked[:args.causes]:
            say(f"  {a} | {b} | {c} | {d}: {count}")
            shown += count
        if shown < len(mismatch):
            say(f"  remaining {len(mismatch) - shown} in smaller groups of the same clauses")

        # Bulk ice I is the diamond net: every site sits on 12 primitive
        # six-rings and on no 3-, 4- or 5-ring. Neither rule reads that.
        q6g, span = global_q6(frame, ref_edges), xy_span(frame)
        say(f"ice {name} global Steinhardt Q6 over {len(ref_edges)} bonds {q6g:.4f}; "
            f"minimum-image xy span of the oxygens {span[0]:.2f} by {span[1]:.2f} A")
        tri = sum(1 for a in range(n) for b in adj[a] if b > a for c in adj[a] & adj[b] if c > b)
        say(f"ice {name} reference triangles {tri} d-SEAMS 3-rings {len(dsets.get(3, ()))}")
        site = np.zeros((n, 7), dtype=np.int64)
        for size, group in dsets.items():
            if size < 7:
                for r in group:
                    for v in r:
                        site[v, size] += 1
        bulk = (site[:, 3] == 0) & (site[:, 4] == 0) & (site[:, 5] == 0) & (site[:, 6] == 12)
        small = (site[:, 3] + site[:, 4] + site[:, 5]) > 0
        for tag, mask in (("CHILL+ ice", seams_m == "ice"), ("SQQ ice", s_lab == "ice"),
                          ("CHILL+ interfacial", seams_m == "interfacial"),
                          ("SQQ interfacial", s_lab == "interfacial")):
            k = int(mask.sum())
            if k == 0:
                say(f"ice {name} {tag}: 0 molecules")
                continue
            med6 = float(np.median(site[mask, 6]))
            say(f"ice {name} {tag}: {k} molecules; bulk ice-I ring signature "
                f"(12 six-rings, no 3/4/5-ring) {int((bulk & mask).sum())}; on a 3/4/5-ring "
                f"{int((small & mask).sum())}; median primitive six-rings {med6:g}")
        for i in mismatch[:args.examples]:
            say(f"  example {name} molecule {int(i)} CHILL+ {SEAMS_NAME[int(labels[i])]} "
                f"SQQ {s_lab[i]} degree {int(deg[i])} six-rings {int(six[i])} "
                f"neighbour degrees {sorted(int(deg[j]) for j in adj[i])} "
                f"c_ij {' '.join(f'{c:.3f}' for c in cij[i])}")
        rows.append(dict(frame=name, molecules=n, quantity="ice label (merged)",
                         d_seams=int(np.sum(seams_m != "liquid")),
                         sqq=int(np.sum(s_lab != "liquid")),
                         both=int(n - len(mismatch)), d_seams_only="", sqq_only="",
                         agreement=f"{(n - len(mismatch)) / n:.6f}",
                         resolution=(f"CHILL+ reproduced {int(np.sum(r_lab == labels))}/{n} by "
                                     f"scipy; SQQ rule reproduced {int(np.sum(ref_lab == s_lab))}/{n}"
                                     f" from the reference graph; mismatches {len(mismatch)}")))

        # SQQ's rule on d-SEAMS rings. The rule itself is checked against
        # sqq_rule on the reference graph; the ring sets against SQQ's.
        missing = {v for e in (ref_edges ^ sqq_e) for v in e}
        near_missing = {i for i in range(n) if i in missing or adj[i] & missing}
        primitive6 = dsets.get(6, set())
        c_lab_raw = np.zeros(n, dtype=np.int32)
        c_six = np.zeros(n, dtype=np.int32)
        for operation in ("ice-primitive", "ice-chordless"):
            got = census.run(sys_, operation, c_lab_raw, c_six, rings=True)
            c_lab = np.array([RING_ICE_LABEL[int(v)] for v in c_lab_raw], dtype=object)
            rule_lab, _deg, rule_six = sqq_rule(adj, got, n)
            rule_ok = int(np.sum(rule_lab == c_lab))
            six_ok = int(np.sum(rule_six == c_six))
            off = np.flatnonzero(c_lab != s_lab)
            if operation == "ice-primitive":
                ring_tag = "primitive"
                say(f"ice {name} {operation}: six-rings {len(got)} equal to seams_rings_find "
                    f"size 6: {got == primitive6}")
                chord_count = np.zeros(n, dtype=np.int64)
                for r in chord6:
                    for v in r:
                        chord_count[v] += 1
                differ = int(sum(1 for i in off if chord_count[i] != c_six[i]
                                 or any(chord_count[j] != c_six[j] for j in adj[i])))
                why = (f"each mismatch has a primitive and a chordless six-ring count that differ "
                       f"on the molecule or a neighbour: {differ}/{len(off)}")
            else:
                ring_tag = "chordless"
                only_d, only_s = got - chord6, chord6 - got
                ok_d = sum(is_chordless(r, adj) for r in only_d)
                ok_s = sum(is_chordless(r, adj) for r in only_s)
                on_d = sum(1 for r in only_d if any(v in missing for v in r))
                on_s = sum(1 for r in only_s if any(v in missing for v in r))
                if n <= args.enumerate_max:
                    ref6 = reference_chordless6(adj)
                    complete = (f"reference enumeration {len(ref6)}, equal to d-SEAMS "
                                f"{ref6 == got}, to SQQ {ref6 == chord6}")
                else:
                    complete = f"reference enumeration skipped above {args.enumerate_max} molecules"
                say(f"rings {name} chordless 6-rings {complete}")
                say(f"rings {name} chordless 6-rings d-SEAMS {len(got)} SQQ ice stage {len(chord6)} "
                    f"both {len(got & chord6)}; d-SEAMS-only {len(only_d)} (reference chordless "
                    f"{ok_d}, on an edge SQQ lacks {on_d}); SQQ-only {len(only_s)} (reference "
                    f"chordless {ok_s}, on an edge SQQ lacks {on_s}); primitive inside "
                    f"{len(primitive6 & got)}/{len(primitive6)}")
                rows.append(dict(frame=name, molecules=n, quantity="chordless 6-rings",
                                 d_seams=len(got), sqq=len(chord6), both=len(got & chord6),
                                 d_seams_only=len(only_d), sqq_only=len(only_s),
                                 agreement=f"{jaccard(got, chord6):.6f}",
                                 resolution=(f"d-SEAMS-only {len(only_d)}, reference chordless "
                                             f"{ok_d}; SQQ-only {len(only_s)}, reference "
                                             f"chordless {ok_s}; every primitive 6-ring inside: "
                                             f"{primitive6 <= got}; {complete}")))
                beside = sum(1 for i in off if i in near_missing)
                why = f"mismatches on or next to an edge SQQ lacks: {beside}/{len(off)}"
            say(f"ice {name} {operation}: SQQ rule re-derived on the reference graph and these "
                f"rings {rule_ok}/{n} labels, {six_ok}/{n} six-ring counts; against SQQ "
                f"{n - len(off)}/{n}; {why}")
            say(f"ice {name} {operation} counts "
                + " ".join(f"{k}={int(np.sum(c_lab == k))}" for k in ("ice", "interfacial", "liquid")))
            rows.append(dict(frame=name, molecules=n,
                             quantity=f"ice label, SQQ rule on d-SEAMS {ring_tag} 6-rings",
                             d_seams=int(np.sum(c_lab != "liquid")),
                             sqq=int(np.sum(s_lab != "liquid")),
                             both=int(n - len(off)), d_seams_only="", sqq_only="",
                             agreement=f"{(n - len(off)) / n:.6f}",
                             resolution=(f"rule reproduced {rule_ok}/{n} from the reference "
                                         f"graph; {why}")))
        engine.lib.seams_system_destroy(sys_)
    fields = ["frame", "molecules", "quantity", "d_seams", "sqq", "both", "d_seams_only",
              "sqq_only", "agreement", "resolution"]
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    log.close()


# ---------------------------------------------------------------- timing

def time_session(args):
    threads = args.threads
    os.sched_setaffinity(0, set(range(threads)))
    engine = Seams(args.lib)
    engine.threads(threads)
    census = RingIce(args.ring_ice_lib, engine)
    version = sqq_version()
    log = open(args.log, "w", encoding="utf-8")

    def say(text):
        log.write(scrub(text) + "\n")
        log.flush()
        print(text, flush=True)

    say(f"d-SEAMS {engine.version} commit {code_state()} build release -O3")
    say(f"SQQ {version} commit {SQQ_COMMIT} engine py")
    say(f"threads {threads} cores {len(os.sched_getaffinity(0))} cpus {cpus_allowed()} "
        f"python {platform.python_version()} warmup 1 n {args.n}")
    groups = []
    if "rings" in args.operations:
        groups.append(("rings", ["rings"]))
    ice_ops = [op for op in ICE_CALLS if op in args.operations]
    if ice_ops:
        groups.append(("ice", ice_ops))
    for name in args.frames:
        frame = load_frame(name)
        n = frame["molecules"]
        sys_ = engine.system(frame)
        built = sqq_frame(frame)
        labels = np.zeros(n, dtype=np.int32)
        six = np.zeros(n, dtype=np.int32)
        for sqq_op, seams_ops in groups:
            cfg = sqq_config(sqq_op, threads)

            def seams_call(op):
                if op == "rings":
                    h = engine.rings_handle(sys_)
                    engine.lib.seams_rings_destroy(h)
                elif op == "ice":
                    engine.chill_plus(sys_, labels)
                else:
                    census.run(sys_, op, labels, six)

            calls = [((op, "d-SEAMS"), lambda op=op: seams_call(op)) for op in seams_ops]
            calls.append(((sqq_op, "SQQ"), lambda: sqq_run(built, cfg)))
            for _key, call in calls:
                call()
            samples = {key: [] for key, _call in calls}
            for k in range(args.n):
                for key, call in (calls if k % 2 == 0 else calls[::-1]):
                    t0 = time.perf_counter()
                    call()
                    samples[key].append((time.perf_counter() - t0) * 1e3)
            for (op, code), vals in samples.items():
                say(f"sample {name} molecules {n} operation {op} threads {threads} "
                    f"code {code} ms {' '.join(f'{v:.4f}' for v in vals)}")
        engine.lib.seams_system_destroy(sys_)
    log.close()


# ---------------------------------------------------------------- provenance

def provenance(args):
    """Every .py file of the installed sqq against the same path at the tag."""
    import sqq
    src = Path(args.sqq_src)
    head = subprocess.run(["git", "-C", str(src), "rev-parse", f"{args.tag}^{{commit}}"],
                          capture_output=True, text=True, check=True).stdout.strip()
    site = Path(sqq.__file__).resolve().parent
    files = sorted(site.rglob("*.py"))
    differ = []
    for path in files:
        rel = path.relative_to(site.parent).as_posix()
        shown = subprocess.run(["git", "-C", str(src), "show", f"{args.tag}:{rel}"],
                               capture_output=True, check=False)
        if shown.returncode != 0 or shown.stdout != path.read_bytes():
            differ.append(rel)
    lines = [
        f"sqq {sqq.__version__} installed from PyPI",
        f"github.com/pimooni/sqq tag {args.tag} commit {head}",
        f"paper Pang and Sun, J. Phys. Chem. A 130 (41), 8187-8199 (2026), doi {SQQ_DOI}",
        f"installed .py files {len(files)} identical to the tag {len(files) - len(differ)}",
    ] + [f"differs {rel}" for rel in differ]
    if head != SQQ_COMMIT:
        lines.append(f"tag commit {head} is not the recorded {SQQ_COMMIT}")
    Path(args.out).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


# ---------------------------------------------------------------- summary

def paired_ratio_ci(dseams, sqq, seed=BOOT_SEED, boot=BOOT):
    rng = random.Random(seed)
    n = len(dseams)
    ratios = []
    for _ in range(boot):
        pick = [rng.randrange(n) for _ in range(n)]
        ratios.append(sum(sqq[i] for i in pick) / sum(dseams[i] for i in pick))
    ratios.sort()
    lo = ratios[int(0.025 * boot)]
    hi = ratios[int(math.ceil(0.975 * boot)) - 1]
    return lo, hi


def parse_log(path):
    meta, samples = {}, {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if line.startswith("threads "):
            meta["threads"] = int(parts[1])
            meta["cores"] = int(parts[3])
            meta["cpus"] = parts[5]
        elif line.startswith("d-SEAMS "):
            meta["dseams"] = parts[1]
            meta["commit"] = parts[3]
            meta["code"] = parts[5] if len(parts) > 5 and parts[4] == "code" else ""
        elif line.startswith("SQQ "):
            meta["sqq"] = parts[1]
        elif line.startswith("sample "):
            key = (parts[1], int(parts[3]), parts[5], int(parts[7]), parts[9])
            samples[key] = [float(v) for v in parts[11:]]
    return meta, samples


def summarize(args):
    agreement = {}
    with open(args.agree, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            agreement[(row["frame"], row["quantity"])] = row["agreement"]
    fields = ["frame", "molecules", "operation", "threads", "cores", "cpus", "n",
              "dseams_mean_ms", "dseams_sd_ms", "sqq_mean_ms", "sqq_sd_ms",
              "ratio_sqq_over_dseams", "ratio_ci95_lo", "ratio_ci95_hi", "bootstrap",
              "agreement", "dseams_version", "dseams_commit", "dseams_code", "sqq_version",
              "sqq_commit", "sqq_engine", "note"]
    out_rows = []
    for path in args.logs:
        meta, samples = parse_log(path)
        for (name, n, op, threads, code), vals in samples.items():
            if code != "d-SEAMS":
                continue
            other = samples[(name, n, "rings" if op == "rings" else "ice", threads, "SQQ")]
            lo, hi = paired_ratio_ci(vals, other)
            if op == "rings":
                keys = [(name, f"{k}-rings") for k in (4, 5, 6)]
                note = "d-SEAMS primitive 3-6 vs SQQ shortest_path 4-6"
            elif op == "ice":
                keys = [(name, "ice label (merged)")]
                note = ("whole tools: d-SEAMS CHILL+ (published; no rings) vs SQQ ice chill "
                        "on chordless 6-rings")
            elif op == "ice-primitive":
                keys = [(name, "ice label, SQQ rule on d-SEAMS primitive 6-rings")]
                note = ("SQQ rule on d-SEAMS primitive 6-rings (bench/ring_ice.c) vs SQQ ice; "
                        "SQQ samples shared by the three ice rows")
            else:
                keys = [(name, "ice label, SQQ rule on d-SEAMS chordless 6-rings"),
                        (name, "chordless 6-rings")]
                note = ("like for like: SQQ rule on chordless 6-cycles of the d-SEAMS graph "
                        "(bench/ring_ice.c) vs SQQ ice; SQQ samples shared by the three ice rows")
            agree_vals = [float(agreement[k]) for k in keys if k in agreement]
            out_rows.append(dict(
                frame=name, molecules=n, operation=op, threads=threads, cores=meta["cores"],
                cpus=meta["cpus"], n=len(vals),
                dseams_mean_ms=f"{statistics.fmean(vals):.4f}",
                dseams_sd_ms=f"{statistics.stdev(vals):.4f}",
                sqq_mean_ms=f"{statistics.fmean(other):.4f}",
                sqq_sd_ms=f"{statistics.stdev(other):.4f}",
                ratio_sqq_over_dseams=f"{sum(other) / sum(vals):.3f}",
                ratio_ci95_lo=f"{lo:.3f}", ratio_ci95_hi=f"{hi:.3f}",
                bootstrap=f"paired percentile B={BOOT} seed={BOOT_SEED}",
                agreement=f"{min(agree_vals):.6f}" if agree_vals else "",
                dseams_version=meta["dseams"], dseams_commit=meta["commit"],
                dseams_code=meta.get("code", ""),
                sqq_version=meta["sqq"], sqq_commit=SQQ_COMMIT, sqq_engine="py",
                note=note))
    order = {name: i for i, name in enumerate(FRAMES)}
    out_rows.sort(key=lambda r: (OPERATIONS.index(r["operation"]), order.get(r["frame"], 99),
                                 int(r["threads"])))
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    lib_default = os.environ.get("SEAMS_C_LIB", str(ROOT / "bbdir-release" / "src" / "libseams.so"))
    ring_ice_default = os.environ.get(
        "SEAMS_RING_ICE_LIB", str(ROOT / "bbdir-release" / "bench" / "libseams_ring_ice.so"))
    p = sub.add_parser("agree")
    p.add_argument("--lib", default=lib_default)
    p.add_argument("--ring-ice-lib", default=ring_ice_default)
    p.add_argument("--frames", nargs="+", default=list(FRAMES))
    p.add_argument("--examples", type=int, default=5)
    p.add_argument("--causes", type=int, default=12)
    p.add_argument("--enumerate-max", type=int, default=110592,
                   help="largest frame for the brute-force chordless 6-cycle reference")
    p.add_argument("--out", default="benchmarks/h2h-sqq-agree.csv")
    p.add_argument("--log", default="benchmarks/h2h-sqq-agree.txt")
    p = sub.add_parser("time")
    p.add_argument("--lib", default=lib_default)
    p.add_argument("--ring-ice-lib", default=ring_ice_default)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--n", type=int, default=8)
    p.add_argument("--frames", nargs="+", default=list(FRAMES))
    p.add_argument("--operations", nargs="+", default=list(OPERATIONS))
    p.add_argument("--log", required=True)
    p = sub.add_parser("provenance")
    p.add_argument("--sqq-src", required=True, help="a clone of github.com/pimooni/sqq")
    p.add_argument("--tag", default="v0.5.8")
    p.add_argument("--out", default="benchmarks/h2h-sqq-provenance.txt")
    p = sub.add_parser("summarize")
    p.add_argument("--agree", default="benchmarks/h2h-sqq-agree.csv")
    p.add_argument("--logs", nargs="+", required=True)
    p.add_argument("--out", default="benchmarks/h2h-sqq.csv")
    args = parser.parse_args()
    if args.cmd == "agree":
        agree(args)
    elif args.cmd == "time":
        time_session(args)
    elif args.cmd == "provenance":
        provenance(args)
    else:
        summarize(args)


if __name__ == "__main__":
    main()
