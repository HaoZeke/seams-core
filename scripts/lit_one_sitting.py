#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Literature sitting for ice Ic.

Each code runs in its own process. The parent keeps every process except
the one being timed stopped, so a rival's OpenMP, Rayon, or Qt pool cannot
spin on the cores the sample is using. OMP_PROC_BIND is set only for
d-SEAMS, and only together with OMP_PLACES, the same way scripts/omp_scaling.py
pins its own rows. Rival pools are not given that bind: under taskset it
pins the main thread to one core and every pool created afterwards inherits
the mask.

OpenBLAS, MKL, and NumExpr stay at one thread in every process. The code's
own pool (OpenMP, mdapy, pyscal3, freud, vesin, OVITO) uses the row width.
OVITO_THREAD_COUNT is in the environment before the interpreter starts, so
it is set before ovito is imported.

d-SEAMS is timed on the same boundary as bench/omp-scaling/scaling.csv.
Arrays and the system handle are built before the timer. neighbours, q4,
q6, and rings are the flat C API (seams_cutoff_csr, seams_ql_cutoff,
seams_primitive_rings). CHILL+ is seams_chill_plus with the OVITO edge
rule. The Lechner-Dellago average is not in the q_l sample. The team uses
an active wait so a gap between samples does not park the workers.

An earlier sitting imported every library into one process, or left every
process runnable, with OMP_WAIT_POLICY unset and OMP_PROC_BIND=close. Idle
teams then stayed inside libgomp's spin window and the 8-thread kernel that
scaling.csv records at a few milliseconds was measured at 50-70 ms.
"""

from __future__ import annotations

import csv
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paper_stats as stats  # noqa: E402

CUTOFF = 3.5
MAX_RING = 6
WARMUP = 3
REPEATS = 30
SIZES = (1000, 10648, 110592)
WIDTHS = (1, 8)
PHASE = "Ic"
DEFAULT_CSV = ROOT / "benchmarks" / "literature" / "h2h-lit.csv"

# code -> operations it actually implements. A missing operation is not a
# skipped test; it is not part of the plan.
PLAN = (
    ("d-SEAMS", ("neighbours", "q4", "q6", "chill_plus", "rings")),
    ("freud", ("neighbours", "q4", "q6")),
    ("OVITO", ("neighbours", "chill_plus", "rings")),
    ("mdapy", ("neighbours", "q4", "q6", "chill_plus", "rings")),
    ("pyscal3", ("neighbours", "q4", "q6")),
    ("vesin", ("neighbours",)),
)

ENV_KEYS = (
    "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
    "OMP_WAIT_POLICY", "GOMP_SPINCOUNT", "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "RAYON_NUM_THREADS",
    "OVITO_THREAD_COUNT", "MDAPY_NUM_THREADS", "PYSCAL_NUM_THREADS",
    "SEAMS_FAIR_WIDTH",
)

# scaling.csv "after", idle, 110592 molecules, 8 threads. The gate accepts
# a median within a factor of two of these, which still rejects the 50-70 ms
# collapse.
# Uncontended q6 at 8 threads on 110592 molecules, the window of the
# sitting that still left other pools runnable. OpenBLAS, MKL and NumExpr
# stay at one thread here so they do not form a second team. The gate is
# this window; a median outside it is reported, not widened.
Q6_8_MS = (15.0, 22.0)
SCALING_AFTER_8 = {
    "neighbours": 2.510,
    "q4": 2.367,
    "q6": 3.452,
    "rings": 18.510,
    "chill_plus": 6.949,
}


def ice_ic(molecules: int):
    import numpy as np
    ncell = round((molecules / 8) ** (1.0 / 3.0))
    if 8 * ncell ** 3 != molecules:
        raise SystemExit(f"{molecules} is not an 8-molecule cubic supercell")
    bond = 2.75
    lattice = 4.0 * bond / np.sqrt(3.0)
    side = ncell * lattice
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    pos = np.empty((molecules, 3), dtype=np.float64)
    cursor = 0
    for i in range(ncell):
        for j in range(ncell):
            for k in range(ncell):
                for base in range(8):
                    f = base % 4
                    shift = 0.0 if base < 4 else 0.25
                    pos[cursor, 0] = (i + fcc[f][0] + shift) * lattice
                    pos[cursor, 1] = (j + fcc[f][1] + shift) * lattice
                    pos[cursor, 2] = (k + fcc[f][2] + shift) * lattice
                    cursor += 1
    box9 = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
    cell = np.diag(np.array([side, side, side], dtype=np.float64))
    return pos, box9, cell, side


def write_frames(directory: Path, sizes: tuple[int, ...]) -> None:
    import numpy as np
    directory.mkdir(parents=True, exist_ok=True)
    for molecules in sizes:
        pos, box9, cell, side = ice_ic(molecules)
        np.savez(directory / f"ic-{molecules}.npz", pos=pos, box9=box9, cell=cell, side=np.array(side))


def child_env(code: str, threads: int, lib_dir: Path, lang_dir: Path) -> dict[str, str]:
    """Environment for one code.

    d-SEAMS matches the idle rows of scripts/omp_scaling.py, plus an active
    OpenMP wait. Rivals do not inherit OMP_PROC_BIND. BLAS stays serial so
    it does not build a second pool beside the code's own.
    """
    env = dict(os.environ)
    for key in ENV_KEYS:
        env.pop(key, None)
    width = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_NUM_THREADS"] = width
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = width
    env["MDAPY_NUM_THREADS"] = width
    env["PYSCAL_NUM_THREADS"] = width
    if code == "d-SEAMS":
        env["OMP_PROC_BIND"] = "close"
        env["OMP_PLACES"] = "cores"
        env["OMP_WAIT_POLICY"] = "ACTIVE"
    else:
        # Leave the wait policy unset. Forcing PASSIVE would change the
        # rival's own samples; forcing ACTIVE while another code ran was
        # the collapse. The parent stops this process whenever it is not
        # the one being timed.
        env.pop("OMP_PROC_BIND", None)
        env.pop("OMP_PLACES", None)
    if code == "OVITO":
        env["OVITO_THREAD_COUNT"] = width
        env["QT_QPA_PLATFORM"] = "offscreen"
    lib_path = os.pathsep.join(
        part for part in (str(lib_dir), str(lang_dir), env.get("LD_LIBRARY_PATH", "")) if part)
    env["LD_LIBRARY_PATH"] = lib_path
    env["LIT_CODE_THREADS"] = width
    return env


def env_record(env: dict) -> str:
    return ";".join(
        f"{key}={env[key]}" if key in env else f"{key}=unset" for key in ENV_KEYS)


def taskset_cores(threads: int) -> str:
    if threads < 1:
        raise ValueError("threads")
    return "0" if threads == 1 else f"0-{threads - 1}"


def jobs(sizes: tuple[int, ...]) -> list[dict]:
    rows = []
    for molecules in sizes:
        for threads in WIDTHS:
            for code, operations in PLAN:
                for operation in operations:
                    rows.append({
                        "code": code, "operation": operation,
                        "threads": threads, "molecules": molecules,
                    })
    rows.sort(key=lambda row: (row["molecules"], row["operation"], row["threads"], row["code"]))
    return rows


class Worker:
    def __init__(self, code: str, threads: int, cpus: str, env: dict, script: Path):
        self.code = code
        self.threads = threads
        cmd = [
            "taskset", "-c", cpus, sys.executable, "-u", str(script),
            "--worker", code, str(threads),
        ]
        self.proc = subprocess.Popen(
            cmd, env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, bufsize=1,
        )
        ready = self._read(180)
        if not ready.get("ready"):
            raise RuntimeError(f"{code} t{threads} did not start: {ready}")
        self.version = ready.get("version", "")
        self.timer = ready.get("timer", "")
        self.omp_max = ready.get("omp_max")
        self.affinity = ready.get("affinity", "")
        self.env = ready.get("env", "")
        self._stop()

    def begin(self) -> None:
        """Let this process run. The caller stops every other worker first."""
        self._cont()

    def end(self) -> None:
        """Freeze the process, including a pool that is still spinning."""
        self._stop()

    def call(self, operation: str, molecules: int) -> dict:
        assert self.proc.stdin is not None
        self.proc.stdin.write(json.dumps({
            "operation": operation, "molecules": molecules,
        }) + "\n")
        self.proc.stdin.flush()
        payload = self._read(600)
        if "ms" not in payload:
            raise RuntimeError(f"{self.code} {operation}: {payload}")
        return payload

    def close(self) -> None:
        if self.proc.poll() is None:
            self._cont()
            if self.proc.stdin is not None:
                self.proc.stdin.close()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()

    def _stop(self) -> None:
        if self.proc.poll() is None:
            os.kill(self.proc.pid, signal.SIGSTOP)

    def _cont(self) -> None:
        if self.proc.poll() is None:
            os.kill(self.proc.pid, signal.SIGCONT)

    def _read(self, timeout: float) -> dict:
        assert self.proc.stdout is not None
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.proc.poll() is not None:
                err = ""
                if self.proc.stderr is not None:
                    err = self.proc.stderr.read()[-800:]
                raise RuntimeError(f"{self.code} exited {self.proc.returncode}: {err}")
            line = self.proc.stdout.readline()
            if not line:
                time.sleep(0.02)
                continue
            if line.startswith("{"):
                return json.loads(line)
        raise RuntimeError(f"{self.code} timed out")


def omp_max_threads() -> int:
    import ctypes
    try:
        lib = ctypes.CDLL("libgomp.so.1")
    except OSError:
        return -1
    lib.omp_get_max_threads.restype = ctypes.c_int
    return int(lib.omp_get_max_threads())


def affinity_now() -> str:
    cpus = sorted(os.sched_getaffinity(0))
    text = ",".join(str(cpu) for cpu in cpus)
    ntask = len(os.listdir("/proc/self/task"))
    return f"sched={text};ntask={ntask}"


def load_frame(molecules: int):
    import numpy as np
    path = Path(os.environ["LIT_FRAMES"]) / f"ic-{molecules}.npz"
    data = np.load(path)
    return (
        np.ascontiguousarray(data["pos"]),
        np.ascontiguousarray(data["box9"]),
        np.ascontiguousarray(data["cell"]),
        float(data["side"]),
    )


def degree_check(degrees) -> str:
    import numpy as np
    degrees = np.asarray(degrees, dtype=np.int64)
    if degrees.size == 0:
        return "deg=empty"
    counts = np.bincount(degrees)
    parts = [f"{degree}:{int(count)}" for degree, count in enumerate(counts) if count]
    return "deg=" + ",".join(parts)


def q_check(values) -> str:
    import numpy as np
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    return f"mean={float(values.mean()):.10f};span={float(values.max() - values.min()):.3e}"


def hist_check(prefix: str, labels, names: dict[int, str]) -> str:
    import numpy as np
    labels = np.asarray(labels, dtype=np.int64).reshape(-1)
    counts: dict[str, int] = {}
    for value, count in zip(*np.unique(labels, return_counts=True)):
        name = names.get(int(value), f"id{int(value)}")
        counts[name] = counts.get(name, 0) + int(count)
    body = ",".join(f"{name}:{counts[name]}" for name in sorted(counts))
    return f"{prefix}{body}"


DSEAMS_CHILL = {
    0: "cubic", 1: "hexagonal", 2: "water", 3: "interfacial",
    4: "clathrate", 5: "interclathrate", 6: "unclassified",
    7: "recubic", 8: "rehex",
}
# OVITO ChillPlusModifier.Type integer values used by OVITO 3.17.
OVITO_CHILL = {
    0: "water", 1: "hexagonal", 2: "cubic", 3: "interfacial",
    4: "clathrate", 5: "interclathrate",
}
# mdapy ChillPlus.NAMES order.
MDAPY_CHILL = {
    0: "water", 1: "hexagonal", 2: "cubic", 3: "interfacial",
    4: "clathrate", 5: "interclathrate",
}


def worker_main(code: str, threads: int) -> None:
    import numpy as np
    prepared = prepare(code, threads)
    sys.stdout.write(json.dumps({
        "ready": True,
        "version": prepared["version"],
        "timer": prepared["timer"],
        "omp_max": omp_max_threads(),
        "affinity": affinity_now(),
        "env": env_record(os.environ),
    }) + "\n")
    sys.stdout.flush()
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        request = json.loads(line)
        operation = request["operation"]
        molecules = int(request["molecules"])
        t0 = time.perf_counter()
        prepared["run"](operation, molecules)
        ms = (time.perf_counter() - t0) * 1e3
        check = prepared["digest"](operation, molecules)
        sys.stdout.write(json.dumps({
            "ms": ms, "check": check, "ntask": len(os.listdir("/proc/self/task")),
        }) + "\n")
        sys.stdout.flush()


def prepare(code: str, threads: int) -> dict:
    if code == "d-SEAMS":
        return prepare_dseams(threads)
    if code == "freud":
        return prepare_freud(threads)
    if code == "OVITO":
        return prepare_ovito(threads)
    if code == "mdapy":
        return prepare_mdapy(threads)
    if code == "pyscal3":
        return prepare_pyscal(threads)
    if code == "vesin":
        return prepare_vesin(threads)
    raise SystemExit(f"unknown code {code}")


def prepare_dseams(threads: int) -> dict:
    import ctypes
    import numpy as np
    core = ctypes.CDLL(os.environ["SEAMS_CORE_LIB"])
    api = ctypes.CDLL(os.environ["SEAMS_API_LIB"])
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    vp = ctypes.c_void_p
    core.seams_cutoff_csr.restype = ctypes.c_int
    core.seams_cutoff_csr.argtypes = [
        dp, ctypes.c_int, dp, ctypes.c_double, ip, ip, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    core.seams_ql_cutoff.restype = ctypes.c_int
    core.seams_ql_cutoff.argtypes = [dp, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, dp]
    core.seams_primitive_rings.restype = ctypes.c_int
    core.seams_primitive_rings.argtypes = [
        dp, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, ip, ctypes.c_int]
    api.seams_system_create.restype = vp
    api.seams_system_set_ortho.argtypes = [vp, dp, ctypes.c_int] + [ctypes.c_double] * 6 + [ip, ip]
    api.seams_chill_plus.argtypes = [vp, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int, ip, ip]
    cache: dict = {}

    def slot(molecules: int) -> dict:
        if molecules in cache:
            return cache[molecules]
        pos, box9, _cell, side = load_frame(molecules)
        n = len(pos)
        offsets = np.empty(n + 1, dtype=np.int32)
        need = ctypes.c_int(0)
        rc = core.seams_cutoff_csr(
            pos.ctypes.data_as(dp), n, box9.ctypes.data_as(dp), CUTOFF,
            offsets.ctypes.data_as(ip), None, 0, ctypes.byref(need))
        if rc != 0:
            raise RuntimeError(f"cutoff query returned {rc}")
        handle = api.seams_system_create()
        rc = api.seams_system_set_ortho(
            handle, pos.ctypes.data_as(dp), n, 0.0, 0.0, 0.0, side, side, side, None, None)
        if rc != 0:
            raise RuntimeError(f"seams_system_set_ortho returned {rc}")
        cache[molecules] = {
            "pos": pos, "box9": box9, "n": n,
            "offsets": offsets, "cols": np.empty(max(int(need.value), 1), dtype=np.int32),
            "need": ctypes.c_int(0), "ql": np.empty(n, dtype=np.float64),
            "rings": np.zeros(MAX_RING + 1, dtype=np.int32),
            "labels": np.zeros(n, dtype=np.int32),
            "counts": np.zeros(9, dtype=np.int32),
            "handle": handle,
        }
        return cache[molecules]

    def run(operation: str, molecules: int) -> None:
        s = slot(molecules)
        pos, n, box9 = s["pos"], s["n"], s["box9"]
        if operation == "neighbours":
            rc = core.seams_cutoff_csr(
                pos.ctypes.data_as(dp), n, box9.ctypes.data_as(dp), CUTOFF,
                s["offsets"].ctypes.data_as(ip), s["cols"].ctypes.data_as(ip),
                len(s["cols"]), ctypes.byref(s["need"]))
            if rc != 0:
                raise RuntimeError(f"seams_cutoff_csr returned {rc}")
        elif operation in ("q4", "q6"):
            order = 4 if operation == "q4" else 6
            rc = core.seams_ql_cutoff(
                pos.ctypes.data_as(dp), n, box9.ctypes.data_as(dp), CUTOFF, order,
                s["ql"].ctypes.data_as(dp))
            if rc != 0:
                raise RuntimeError(f"seams_ql_cutoff returned {rc}")
        elif operation == "rings":
            rc = core.seams_primitive_rings(
                pos.ctypes.data_as(dp), n, box9.ctypes.data_as(dp), CUTOFF, MAX_RING,
                s["rings"].ctypes.data_as(ip), len(s["rings"]))
            if rc != 0:
                raise RuntimeError(f"seams_primitive_rings returned {rc}")
        elif operation == "chill_plus":
            # rule 1 is the OVITO edge. graph 0 is the cutoff shell.
            rc = api.seams_chill_plus(
                s["handle"], CUTOFF, 0, 0, 1, s["labels"].ctypes.data_as(ip),
                s["counts"].ctypes.data_as(ip))
            if rc != 0:
                raise RuntimeError(f"seams_chill_plus returned {rc}")
        else:
            raise RuntimeError(operation)

    def digest(operation: str, molecules: int) -> str:
        s = slot(molecules)
        if operation == "neighbours":
            degrees = np.diff(s["offsets"][:s["n"] + 1])
            return degree_check(degrees)
        if operation in ("q4", "q6"):
            return q_check(s["ql"])
        if operation == "rings":
            return f"six={int(s['rings'][6])}"
        if operation == "chill_plus":
            return hist_check("", s["labels"], DSEAMS_CHILL)
        raise RuntimeError(operation)

    return {
        "run": run, "digest": digest, "version": "2.14.0",
        "timer": "arrays and the system handle are built before the timer; "
                 "the flat C call (or seams_chill_plus) is inside; q_l has no q_l-bar",
    }


def prepare_freud(threads: int) -> dict:
    import freud
    import numpy as np
    freud.parallel.set_num_threads(threads)
    cache: dict = {}

    def slot(molecules: int):
        if molecules not in cache:
            pos, _box9, _cell, side = load_frame(molecules)
            pts = np.ascontiguousarray(pos - 0.5 * side)
            box = freud.box.Box.cube(side)
            cache[molecules] = {"pts": pts, "box": box, "n": len(pos), "ql": None, "deg": None}
        return cache[molecules]

    def run(operation: str, molecules: int) -> None:
        s = slot(molecules)
        if operation == "neighbours":
            query = freud.locality.AABBQuery(s["box"], s["pts"])
            found = query.query(s["pts"], {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
            s["deg"] = np.bincount(found.query_point_indices, minlength=s["n"])
        elif operation in ("q4", "q6"):
            order = 4 if operation == "q4" else 6
            op = freud.order.Steinhardt(l=order, average=False)
            op.compute((s["box"], s["pts"]), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            s["ql"] = np.array(op.particle_order, copy=True)
        else:
            raise RuntimeError(operation)

    def digest(operation: str, molecules: int) -> str:
        s = slot(molecules)
        if operation == "neighbours":
            return degree_check(s["deg"])
        return q_check(s["ql"])

    return {
        "run": run, "digest": digest, "version": freud.__version__,
        "timer": "box built before the timer; the query or Steinhardt.compute is inside",
    }


def prepare_ovito(threads: int) -> dict:
    # OVITO_THREAD_COUNT was set in the process environment before Python started.
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier, CreateBondsModifier, FindRingsModifier
    import ovito
    cache: dict = {}

    def fresh(molecules: int):
        pos, _box9, cell, _side = load_frame(molecules)
        data = DataCollection()
        matrix = np.zeros((3, 4), dtype=np.float64)
        matrix[:, :3] = cell
        data.create_cell(matrix=matrix, pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=np.ascontiguousarray(pos))
        return data

    def run(operation: str, molecules: int) -> None:
        data = fresh(molecules)
        if operation == "neighbours":
            data.apply(CreateBondsModifier(cutoff=CUTOFF))
            bonds = np.asarray(data.particles.bonds.topology)
            degree = np.bincount(bonds.reshape(-1), minlength=data.particles.count)
            cache[molecules] = {"kind": "neighbours", "deg": degree}
        elif operation == "chill_plus":
            data.apply(ChillPlusModifier(cutoff=CUTOFF))
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int64)
            cache[molecules] = {"kind": "chill_plus", "labels": labels}
        elif operation == "rings":
            data.apply(CreateBondsModifier(cutoff=CUTOFF))
            data.apply(FindRingsModifier(minimum_ring_size=3, maximum_ring_size=6))
            table = data.tables["ring-size-histogram"]
            found = {
                int(size): int(count)
                for size, count in zip(np.asarray(table["Ring Size"]), np.asarray(table["Count"]))
            }
            cache[molecules] = {"kind": "rings", "six": found.get(6, 0)}
        else:
            raise RuntimeError(operation)

    def digest(operation: str, molecules: int) -> str:
        got = cache[molecules]
        if operation == "neighbours":
            return degree_check(got["deg"])
        if operation == "chill_plus":
            return hist_check("", got["labels"], OVITO_CHILL)
        return f"six={int(got['six'])}"

    version = getattr(ovito, "version_string", None) or ".".join(str(v) for v in tuple(ovito.version)[:3])
    return {
        "run": run, "digest": digest, "version": version,
        "timer": "OVITO_THREAD_COUNT was set before import; the modifier pipeline is inside the timer",
    }


def prepare_mdapy(threads: int) -> dict:
    import mdapy
    import numpy as np
    from mdapy import System, set_num_threads
    from mdapy.analysis import RingStatistics, SteinhardtParameters
    from mdapy.identify import ChillPlus
    from mdapy.neighbor import NeighborList
    set_num_threads(threads)
    cache: dict = {}

    def system(molecules: int):
        if molecules not in cache:
            pos, _box9, cell, _side = load_frame(molecules)
            cache[molecules] = {
                "system": System(positions=pos, box=cell),
                "q4": SteinhardtParameters(degrees=(4,), cutoff=CUTOFF),
                "q6": SteinhardtParameters(degrees=(6,), cutoff=CUTOFF),
                "rings": RingStatistics(max_size=6, min_size=3, cutoff=CUTOFF),
                "chill": ChillPlus(cutoff=CUTOFF),
                "last": None,
            }
        return cache[molecules]

    def run(operation: str, molecules: int) -> None:
        built = system(molecules)
        set_num_threads(threads)
        if operation == "neighbours":
            listed = NeighborList.build(built["system"], cutoff=CUTOFF)
            degree = np.diff(np.asarray(listed.offsets))
            built["last"] = ("neighbours", degree)
        elif operation in ("q4", "q6"):
            result = built[operation].compute(built["system"])
            built["last"] = ("q", np.asarray(result.q).reshape(-1))
        elif operation == "chill_plus":
            result = built["chill"].compute(built["system"])
            built["last"] = ("chill", np.asarray(result.structure))
        elif operation == "rings":
            result = built["rings"].compute(built["system"])
            sizes = np.asarray(result.sizes)
            six = int(np.count_nonzero(sizes == 6))
            built["last"] = ("rings", six)
        else:
            raise RuntimeError(operation)

    def digest(operation: str, molecules: int) -> str:
        kind, value = system(molecules)["last"]
        if kind == "neighbours":
            return degree_check(value)
        if kind == "q":
            return q_check(value)
        if kind == "chill":
            return hist_check("", value, MDAPY_CHILL)
        return f"six={int(value)}"

    return {
        "run": run, "digest": digest, "version": mdapy.__version__,
        "timer": "System built before the timer; the neighbour, Steinhardt, CHILL+, or ring compute is inside",
    }


def prepare_pyscal(threads: int) -> dict:
    import numpy as np
    import pyscal3
    from ase import Atoms
    pyscal3.set_num_threads(threads)
    cache: dict = {}

    def atoms(molecules: int):
        if molecules not in cache:
            pos, _box9, cell, _side = load_frame(molecules)
            cache[molecules] = Atoms(
                symbols=["O"] * len(pos), positions=np.ascontiguousarray(pos),
                cell=cell, pbc=True)
        return cache[molecules]

    last: dict = {}

    def run(operation: str, molecules: int) -> None:
        built = atoms(molecules)
        pyscal3.set_num_threads(threads)
        if operation == "neighbours":
            pyscal3.find_neighbors(built, method="cutoff", cutoff=CUTOFF)
            # The per-atom row list is omitted once the bond CSR is stored.
            degree = np.diff(np.asarray(built.info["pyscal_bond_offsets"]))
            last[molecules] = ("neighbours", degree)
        elif operation in ("q4", "q6"):
            order = 4 if operation == "q4" else 6
            pyscal3.find_neighbors(built, method="cutoff", cutoff=CUTOFF)
            values = np.asarray(pyscal3.steinhardt_parameter(built, order)[0]).reshape(-1)
            last[molecules] = ("q", values)
        else:
            raise RuntimeError(operation)

    def digest(operation: str, molecules: int) -> str:
        kind, value = last[molecules]
        if kind == "neighbours":
            return degree_check(value)
        return q_check(value)

    return {
        "run": run, "digest": digest, "version": pyscal3.__version__,
        "timer": "Atoms built before the timer; find_neighbors and, for q_l, the Steinhardt call are inside",
    }


def prepare_vesin(threads: int) -> dict:
    import numpy as np
    import vesin
    cache: dict = {}

    def slot(molecules: int):
        if molecules not in cache:
            pos, _box9, cell, _side = load_frame(molecules)
            cache[molecules] = {
                "finder": vesin.NeighborList(cutoff=CUTOFF, full_list=True, n_threads=threads),
                "pos": pos, "cell": cell, "deg": None,
            }
        return cache[molecules]

    def run(operation: str, molecules: int) -> None:
        if operation != "neighbours":
            raise RuntimeError(operation)
        s = slot(molecules)
        left, _right = s["finder"].compute(
            points=s["pos"], box=s["cell"], periodic=True, quantities="ij", copy=True)
        s["deg"] = np.bincount(np.asarray(left, dtype=np.int64), minlength=len(s["pos"]))

    def digest(operation: str, molecules: int) -> str:
        return degree_check(slot(molecules)["deg"])

    return {
        "run": run, "digest": digest, "version": vesin.__version__,
        "timer": "NeighborList built before the timer; compute is inside",
    }


def parse_q(check: str) -> tuple[float, float]:
    parts = dict(item.split("=", 1) for item in check.split(";"))
    return float(parts["mean"]), float(parts["span"])


def agree(operation: str, reference: str, other: str) -> tuple[bool, str]:
    """Whether two check strings describe the same result, and why if not."""
    if operation == "neighbours":
        if reference == other:
            return True, ""
        return False, (
            f"per-atom coordination {other} against d-SEAMS {reference}. "
            "Both counts are the number of neighbours inside 3.5 A, self excluded."
        )
    if operation in ("q4", "q6"):
        mean_a, span_a = parse_q(reference)
        mean_b, span_b = parse_q(other)
        if abs(mean_a - mean_b) <= 1e-4 and span_a <= 1e-3 and span_b <= 1e-3:
            return True, ""
        if abs(mean_a - mean_b) <= 1e-3:
            return True, (
                f"means differ by {abs(mean_a - mean_b):.3e}, inside 1e-3, "
                "on the cutoff Steinhardt parameter without the Lechner-Dellago average"
            )
        return False, (
            f"cutoff Steinhardt means {mean_a:.6f} and {mean_b:.6f} "
            f"(spans {span_a:.3e}, {span_b:.3e})"
        )
    if operation == "chill_plus":
        if reference == other:
            return True, ""
        return False, (
            f"CHILL+ classes {other} against d-SEAMS OVITO-rule {reference}. "
            "Names are aligned (cubic, hexagonal, water, interfacial, clathrate); "
            "the integer codes of the libraries are not the same numbers."
        )
    if operation == "rings":
        if reference == other:
            return True, ""
        return False, (
            f"six-ring count {other} against d-SEAMS primitive rings {reference}. "
            "d-SEAMS counts Franzblau primitive rings. OVITO FindRings counts rings "
            "on the bond network and is a different census. mdapy counts shortest-path "
            "rings. A different six-ring count is that definition, not a failed run."
        )
    raise ValueError(operation)


def scaling_reference(path: Path) -> dict[str, float]:
    found = dict(SCALING_AFTER_8)
    if not path.is_file():
        return found
    with path.open(newline="") as handle:
        # The file's first line is a comment, then a header.
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        if row.get("code") != "dseams" or row.get("build") != "after":
            continue
        if row.get("load") != "idle" or row.get("threads") != "8":
            continue
        if row.get("molecules") != "110592":
            continue
        operation = row["operation"]
        if operation == "chill_plus_ovito":
            found["chill_plus"] = float(row["median_ms"])
        elif operation in found:
            found[operation] = float(row["median_ms"])
    return found


def within_scaling(median: float, reference: float) -> bool:
    if reference <= 0:
        return False
    ratio = median / reference
    return 0.5 <= ratio <= 2.0


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker_main(sys.argv[2], int(sys.argv[3]))
        return 0
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks" / "literature")
    parser.add_argument("--lib", type=Path, default=ROOT / "bbdir-release" / "src" / "libyodaLib.so")
    parser.add_argument("--api", type=Path, default=ROOT / "bbdir-release" / "src" / "libseams.so")
    parser.add_argument("--repeats", type=int, default=int(os.environ.get("LIT_N", str(REPEATS))))
    parser.add_argument("--warmup", type=int, default=int(os.environ.get("LIT_WARMUP", str(WARMUP))))
    parser.add_argument("--sizes", default=",".join(str(s) for s in SIZES))
    parser.add_argument(
        "--handoff", type=Path,
        default=Path(os.environ.get("LIT_HANDOFF", str(DEFAULT_CSV))),
        help="second copy of h2h-lit.csv; skipped when it is the --out table itself")
    args = parser.parse_args()
    if args.repeats < REPEATS and os.environ.get("LIT_ALLOW_SHORT") != "1":
        parser.error(f"--repeats must be at least {REPEATS}")
    sizes = tuple(int(part) for part in args.sizes.split(",") if part)
    for molecules in sizes:
        ice_ic(molecules)
    only = os.environ.get("LIT_ONLY", "").strip()
    partial = bool(only)

    import release_build
    build_dir = release_build.build_dir_of(args.lib.resolve())
    if build_dir is None:
        raise SystemExit(f"refusing to time {args.lib}: not inside a meson build")
    build_label = release_build.assert_release(build_dir)
    if not args.api.is_file():
        raise SystemExit(f"missing {args.api}")

    frames = args.out / "frames"
    write_frames(frames, sizes)
    plan = jobs(sizes)
    if partial:
        code, operation, width = only.split(",")
        width = int(width)
        plan = [
            row for row in plan
            if row["code"] == code and row["operation"] == operation and row["threads"] == width
        ]
        if not plan:
            raise SystemExit(f"LIT_ONLY matches nothing: {only}")
    model = stats.machine(ROOT)["cpu_model"]
    cores, affinity = stats.cores_now()
    commit = stats.source_commit(ROOT)
    tree = stats.source_tree(ROOT)
    reference = scaling_reference(ROOT / "bench" / "omp-scaling" / "scaling.csv")

    keys = []
    for row in plan:
        key = (row["code"], row["threads"])
        if key not in keys:
            keys.append(key)
    workers: dict[tuple, Worker] = {}
    samples: dict[tuple, list[float]] = {}
    checks: dict[tuple, str] = {}
    versions: dict[tuple, str] = {}
    timers: dict[tuple, str] = {}
    omp_max: dict[tuple, object] = {}
    affinities: dict[tuple, str] = {}
    envs: dict[tuple, str] = {}
    script = Path(__file__).resolve()
    try:
        for code, threads in keys:
            print(f"start {code} t{threads}", flush=True)
            env = child_env(code, threads, args.lib.parent, args.lib.parent.parent / "bindings")
            env["LIT_FRAMES"] = str(frames)
            env["SEAMS_CORE_LIB"] = str(args.lib.resolve())
            env["SEAMS_API_LIB"] = str(args.api.resolve())
            worker = Worker(code, threads, taskset_cores(threads), env, script)
            workers[(code, threads)] = worker
            versions[(code, threads)] = worker.version
            timers[(code, threads)] = worker.timer
            omp_max[(code, threads)] = worker.omp_max
            affinities[(code, threads)] = worker.affinity
            envs[(code, threads)] = worker.env
        for row in plan:
            samples[(row["code"], row["operation"], row["threads"], row["molecules"])] = []
        args.out.mkdir(parents=True, exist_ok=True)
        log_path = args.out / "raw-lit.txt"
        with log_path.open("w", encoding="utf-8") as handle:
            header = [
                f"cpu_model={model}",
                f"cores={cores}",
                f"affinity={affinity}",
                f"seams_commit={commit}",
                f"source_tree={tree}",
                f"build={build_label}",
                f"warmup={args.warmup}",
                f"n={args.repeats}",
                f"cutoff_A={CUTOFF}",
                "layout=one process per code and thread width; every other process is SIGSTOP'd during a sample",
                "d-SEAMS OMP_PROC_BIND=close OMP_PLACES=cores OMP_WAIT_POLICY=ACTIVE OPENBLAS_NUM_THREADS=1",
                "rivals OMP_PROC_BIND=unset OMP_PLACES=unset OPENBLAS_NUM_THREADS=1; their own pool uses the row width",
                "OVITO_THREAD_COUNT is set before the OVITO process starts",
                "timer=flat C API for neighbours, q4, q6, rings; seams_chill_plus OVITO rule for CHILL+; no q_l-bar",
                "cause=a previous sitting left rival thread pools runnable, so idle OpenMP teams spun on the same cores",
            ]
            handle.write("\n".join(header) + "\n")
            for code, threads in keys:
                handle.write(
                    f"worker code={code} threads={threads} version={versions[(code, threads)]} "
                    f"omp_max={omp_max[(code, threads)]} affinity={affinities[(code, threads)]} "
                    f"env={envs[(code, threads)]}\n")
            total = args.warmup + args.repeats
            # One code's samples run back to back. Stopping that process
            # between samples freezes its OpenMP team, and the next sample
            # pays the wake-up. Other codes stay stopped for the whole block.
            for row in plan:
                print(
                    f"{row['code']} {row['molecules']} {row['operation']} t{row['threads']}",
                    flush=True)
                worker = workers[(row["code"], row["threads"])]
                key = (row["code"], row["operation"], row["threads"], row["molecules"])
                worker.begin()
                try:
                    for round_index in range(total):
                        payload = worker.call(row["operation"], row["molecules"])
                        if round_index >= args.warmup:
                            samples[key].append(float(payload["ms"]))
                            checks[key] = payload["check"]
                        handle.write(
                            f"round={round_index} molecules={row['molecules']} operation={row['operation']} "
                            f"code={row['code']} threads={row['threads']} ms={payload['ms']:.6f} "
                            f"check={payload['check']} ntask={payload.get('ntask', '')} "
                            f"warmup={int(round_index < args.warmup)}\n")
                finally:
                    worker.end()
                handle.flush()
    finally:
        for worker in workers.values():
            worker.close()

    fields = [
        "entry", "phase", "molecules", "threads", "operation", "code", "version",
        "median_ms", "min_ms", "max_ms", "q1_ms", "q3_ms", "iqr_ms", "samples_ms",
        "timer", "check", "agreement", "note", "seams_commit", "source_tree",
        "build_type", "warmup", "n", "cpus", "cores", "omp_proc_bind", "omp_places",
        "omp_num_threads", "omp_wait_policy", "omp_max", "affinity", "env", "cpu_model",
    ]
    out_csv = args.out / "h2h-lit.csv"
    notes = []
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in plan:
            key = (row["code"], row["operation"], row["threads"], row["molecules"])
            summary = stats.summarize(samples[key])
            ref_key = ("d-SEAMS", row["operation"], row["threads"], row["molecules"])
            if row["code"] == "d-SEAMS" or ref_key not in checks:
                agreement = "reference" if row["code"] == "d-SEAMS" else "uncompared"
                note = "" if ref_key in checks or row["code"] == "d-SEAMS" else "partial plan has no d-SEAMS row"
                if row["molecules"] < 2048 and row["threads"] > 1 and row["code"] == "d-SEAMS":
                    note = (
                        "below 2048 atoms the OpenMP team is one thread; "
                        "eight threads are slower than the serial walk on this call"
                    )
                elif row["molecules"] < 8192 and row["threads"] > 4:
                    note = "below 8192 atoms the OpenMP team is at most four threads"
            else:
                ok, note = agree(row["operation"], checks[ref_key], checks[key])
                agreement = "agree" if ok else "explained" if note else "differ"
                if not ok:
                    notes.append(
                        f"{row['molecules']} {row['operation']} t{row['threads']} "
                        f"{row['code']}: {note}")
            env = envs[(row["code"], row["threads"])]
            proc_bind = "close" if row["code"] == "d-SEAMS" else "unset"
            places = "cores" if row["code"] == "d-SEAMS" else "unset"
            wait = "ACTIVE" if row["code"] == "d-SEAMS" else "unset"
            writer.writerow({
                "entry": "capi" if row["code"] == "d-SEAMS" and row["operation"] != "chill_plus" else (
                    "system" if row["code"] == "d-SEAMS" else "native"),
                "phase": PHASE,
                "molecules": row["molecules"],
                "threads": row["threads"],
                "operation": row["operation"],
                "code": row["code"],
                "version": versions[(row["code"], row["threads"])],
                "median_ms": f"{summary['median_ms']:.6f}",
                "min_ms": f"{min(samples[key]):.6f}",
                "max_ms": f"{max(samples[key]):.6f}",
                "q1_ms": f"{summary['q1_ms']:.6f}",
                "q3_ms": f"{summary['q3_ms']:.6f}",
                "iqr_ms": f"{summary['iqr_ms']:.6f}",
                "samples_ms": stats.fmt_samples(samples[key]),
                "timer": timers[(row["code"], row["threads"])],
                "check": checks[key],
                "agreement": agreement,
                "note": note,
                "seams_commit": commit,
                "source_tree": tree,
                "build_type": build_label,
                "warmup": args.warmup,
                "n": summary["n"],
                "cpus": taskset_cores(row["threads"]),
                "cores": cores,
                "omp_proc_bind": proc_bind,
                "omp_places": places,
                "omp_num_threads": row["threads"],
                "omp_wait_policy": wait,
                "omp_max": omp_max[(row["code"], row["threads"])],
                "affinity": affinities[(row["code"], row["threads"])],
                "env": env,
                "cpu_model": model,
            })
            print(
                f"{row['code']:8} {row['molecules']:7} {row['operation']:12} t{row['threads']} "
                f"{summary['median_ms']:10.3f} ms  {agreement} {checks[key]}",
                flush=True)

    failures = []
    full = (not partial) and set(SIZES) <= set(sizes) and args.repeats >= REPEATS
    if not full:
        print(f"wrote {out_csv}", flush=True)
        return 0
    for operation, ref in reference.items():
        key = ("d-SEAMS", operation, 8, 110592)
        if key not in samples or len(samples[key]) < 2:
            continue
        median = stats.summarize(samples[key])["median_ms"]
        if not within_scaling(median, ref):
            failures.append(
                f"d-SEAMS {operation} 8 threads median {median:.3f} ms is outside "
                f"half to twice scaling.csv after {ref:.3f} ms")
    low, high = Q6_8_MS
    for code in ("mdapy", "pyscal3"):
        key = (code, "q6", 8, 110592)
        if 110592 not in sizes:
            continue
        if key not in samples:
            failures.append(f"{code} q6 at 8 threads was not measured")
            continue
        median = stats.summarize(samples[key])["median_ms"]
        if not low <= median <= high:
            failures.append(
                f"{code} q6 at 8 threads median {median:.3f} ms is outside {low:.0f}-{high:.0f} ms")
    for note in notes:
        if "differ" in note and "definition" not in note and "integer codes" not in note:
            failures.append(note)
    # A neighbour or q_l disagreement is a failure. Ring-definition notes are recorded.
    for row in plan:
        if row["code"] == "d-SEAMS":
            continue
        key = (row["code"], row["operation"], row["threads"], row["molecules"])
        ref_key = ("d-SEAMS", row["operation"], row["threads"], row["molecules"])
        ok, note = agree(row["operation"], checks[ref_key], checks[key])
        if not ok and row["operation"] != "rings":
            failures.append(
                f"{row['code']} {row['operation']} {row['molecules']} t{row['threads']}: {note}")
    handoff = args.handoff.resolve()
    if handoff != out_csv.resolve():
        handoff.parent.mkdir(parents=True, exist_ok=True)
        handoff.write_text(out_csv.read_text(encoding="utf-8"), encoding="utf-8")
        print(f"wrote {handoff}", flush=True)
    print(f"wrote {out_csv}", flush=True)
    if failures:
        print("GATES FAILED", flush=True)
        for line in failures:
            print(line, flush=True)
        return 1
    print("gates passed", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
