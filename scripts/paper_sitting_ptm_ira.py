#!/usr/bin/env python3
"""PTM and IRA structure identification at 884736 molecules.

The frames are the ones scripts/paper_sitting_large.py times. Ice is
ideal ice Ic, 48^3 conventional cells. The liquid is the 4096-molecule
TIP4P/ice frame in benchmarks/paper-cells/tip4pice-liquid-4096.txt,
replicated 6 by 6 by 6.

d-SEAMS is seams_symmetry_classify. OVITO 3.17.2 is
PolyhedralTemplateMatchingModifier: the default template set, and the
set with only cubic and hexagonal diamond, both at RMSD cutoff 0.1.
IRA is ira_mod.IRA.match, the Python interface. It scores the 17-atom
ice Ih, ice Ic and ice Isd templates by Hausdorff distance and rejects
a distance above 0.2 A. The search-radius factor is 1.5. Those gates
are the ones already used for the published rows. Nothing here moves
them to raise the agreement.

One warm-up, then eight calls, serial and eight threads, in this
process. Each width is its own process, pinned with taskset. The core
count is the lscpu CPU count. The host name is not read.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_sitting_large as large
import paper_stats as stats

MOLECULES = large.MOLECULES
WARMUP = 1
REPEATS = 8
WIDTHS = (1, 8)
HAUSDORFF_A = 0.2
KMAX_FACTOR = 1.5
# A Hausdorff distance below this is a congruent template. The search
# stops there because no other template can score lower. It is not the
# 0.2 A classification gate.
CONGRUENT_A = 1.0e-8
PTM_RMSD = 0.1
SHELL_CUTOFF = 5.0
BOND_CUTOFF = 3.5
BOND_IGNORE = 2.15
BOND_LO = 2.30
BOND_HI = 3.35

# d-SEAMS kind, then the OVITO Structure Type of the same crystal.
# ice Isd has no PTM template.
DSEAMS_NAME = {-1: "other", 0: "Ih", 1: "Ic", 2: "Isd"}
PTM_NAME = {
    0: "Other",
    1: "FCC",
    2: "HCP",
    3: "BCC",
    4: "ICO",
    5: "SC",
    6: "cubic-diamond",
    7: "hex-diamond",
    8: "graphene",
}
# Classes both the classifier and ice-only PTM can name.
PTM_OF = {-1: 0, 0: 7, 1: 6}
IRA_NAME = {-1: "other", 0: "Ih", 1: "Ic", 2: "Isd"}

FIELDS = [
    "source_commit", "cpu_model", "cpu_cores", "affinity",
    "phase", "operation", "code", "version",
    "molecules", "threads", "warmup", "n",
    "median_ms", "iqr_ms", "q1_ms", "q3_ms", "samples_ms", "detail",
]
AGREE_FIELDS = [
    "phase", "threads", "left", "left_label", "right", "right_label",
    "atoms", "molecules", "fraction",
]


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload) + "\n")
    sys.stdout.flush()


def scrub(text: str) -> str:
    for prefix in ("home", "tmp", "Users"):
        text = re.sub("/" + prefix + r"/\S+", "<path>", text)
    return text


def project_version() -> str:
    return large.project_version()


def hold_affinity(threads: int) -> None:
    large.hold_affinity(threads)


def frame_of(phase: str):
    if phase == "liquid":
        return large.liquid()
    if phase == "ice":
        return large.ice()
    raise RuntimeError(f"unknown phase {phase}")


def open_core(path: str):
    import ctypes

    lib = ctypes.CDLL(path)
    lib.seams_symmetry_classify.restype = ctypes.c_int
    lib.seams_symmetry_classify.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_template.restype = ctypes.c_int
    lib.seams_symmetry_template.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    return lib


def ice_templates(lib):
    import ctypes
    import numpy as np

    found = []
    for ident in range(8):
        n_out = ctypes.c_int()
        code = lib.seams_symmetry_template(ident, None, 0, ctypes.byref(n_out))
        if code != 0:
            continue
        count = int(n_out.value)
        if count != 17:
            continue
        buf = np.empty(count * 3, dtype=np.float64)
        code = lib.seams_symmetry_template(
            ident,
            buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            count,
            ctypes.byref(n_out),
        )
        if code != 0:
            raise RuntimeError(f"template {ident} was not written")
        found.append((ident, np.ascontiguousarray(buf.reshape(count, 3))))
    if len(found) != 3:
        raise RuntimeError(f"expected three 17-atom ice templates, got {len(found)}")
    idents = np.asarray([ident for ident, _pts in found], dtype=np.int32)
    templates = np.ascontiguousarray(np.stack([pts for _ident, pts in found]))
    return idents, templates


def classify_into(lib, xyz, side, labels) -> None:
    import ctypes
    import numpy as np

    lat = np.zeros(9, dtype=np.float64)
    lat[0] = lat[4] = lat[8] = side
    code = lib.seams_symmetry_classify(
        np.ascontiguousarray(xyz).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(labels.shape[0]),
        lat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if code != 0:
        raise RuntimeError(f"seams_symmetry_classify returned {code}")


def save_array(directory: str, threads: int, phase: str, code: str, operation: str, values) -> None:
    import numpy as np

    path = Path(directory) / f"t{threads}-{phase}-{code}-{operation}.npy"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.asarray(values))


def load_array(directory: Path, threads: int, phase: str, code: str, operation: str):
    import numpy as np

    path = directory / f"t{threads}-{phase}-{code}-{operation}.npy"
    if not path.is_file():
        return None
    return np.load(path)


def time_calls(call, warmup: int, repeats: int):
    samples = []
    last = None
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        last = call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index >= warmup:
            samples.append(elapsed)
    if len(samples) < 8:
        raise RuntimeError(f"n={len(samples)} is below 8")
    return samples, last


def _done(operation: str, code: str, version: str, samples, detail: str) -> None:
    emit({
        "event": "done",
        "operation": operation,
        "code": code,
        "version": version,
        "molecules": MOLECULES,
        "samples_ms": samples,
        "detail": detail,
    })


def worker_dseams(phase: str, threads: int, lib_path: str, agree: str) -> None:
    import numpy as np

    hold_affinity(threads)
    xyz, _box, side = frame_of(phase)
    if len(xyz) != MOLECULES:
        raise RuntimeError(f"{phase} has {len(xyz)} molecules")
    lib = open_core(lib_path)
    labels = np.empty(MOLECULES, dtype=np.int32)

    def run():
        hold_affinity(threads)
        classify_into(lib, xyz, side, labels)
        return int(labels[0])

    samples, _probe = time_calls(run, WARMUP, REPEATS)
    save_array(agree, threads, phase, "d-SEAMS", "classifier", labels)
    counted = {name: int(np.count_nonzero(labels == code)) for code, name in DSEAMS_NAME.items()}
    _done(
        "classifier", "d-SEAMS", project_version(), samples,
        "seams_symmetry_classify "
        + " ".join(f"{name}={counted[name]}" for name in ("other", "Ih", "Ic", "Isd")),
    )


def ovito_labels(xyz, side: float, kind: str, with_rmsd: bool):
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
    parts = data.create_particles(count=len(xyz))
    parts.create_property("Position", data=np.ascontiguousarray(xyz))
    pipe = Pipeline(source=StaticSource(data=data))
    mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=PTM_RMSD)
    if kind == "ice":
        for structure in mod.structures:
            structure.enabled = structure.id in (6, 7)
    elif kind != "default":
        raise RuntimeError(kind)
    mod.output_rmsd = with_rmsd
    pipe.modifiers.append(mod)
    result = pipe.compute()
    labels = np.asarray(result.particles["Structure Type"], dtype=np.int32).copy()
    rmsd = None
    if with_rmsd:
        rmsd = np.asarray(result.particles["RMSD"], dtype=np.float64).copy()
    if labels.shape[0] != len(xyz):
        raise RuntimeError("PTM dropped atoms")
    return labels, rmsd


def worker_ovito(phase: str, threads: int, agree: str) -> None:
    import ovito

    hold_affinity(threads)
    version = ".".join(str(part) for part in tuple(ovito.version)[:3])
    if version != "3.17.2":
        raise RuntimeError(f"OVITO {version} is not 3.17.2")
    xyz, _box, side = frame_of(phase)
    jobs = (
        ("ptm-default", "default", "templates=FCC+HCP+BCC"),
        ("ptm-ice", "ice", "templates=cubic+hexagonal-diamond"),
    )
    for operation, kind, note in jobs:
        def run(kind=kind):
            hold_affinity(threads)
            labels, _rmsd = ovito_labels(xyz, side, kind, False)
            return int(labels.shape[0])

        samples, _probe = time_calls(run, WARMUP, REPEATS)
        labels, rmsd = ovito_labels(xyz, side, kind, True)
        # The untimed RMSD read must not change the structure types.
        check, _ignored = ovito_labels(xyz, side, kind, False)
        if not (labels == check).all():
            raise RuntimeError(f"{operation} labels changed when RMSD was requested")
        save_array(agree, threads, phase, "OVITO", operation, labels)
        save_array(agree, threads, phase, "OVITO", operation + "-rmsd", rmsd)
        _done(
            operation, "OVITO", version, samples,
            f"rows={len(labels)} OVITO_THREAD_COUNT={threads} rmsd={PTM_RMSD} {note}",
        )


_SHELLS = None
_TEMPLATES = None
_IDENTS = None
_IRA = None


def _init_ira() -> None:
    global _IRA
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    import ira_mod

    _IRA = ira_mod.IRA()


def _match_range(bounds):
    import numpy as np

    start, stop = bounds
    ira = _IRA
    labels = np.empty(stop - start, dtype=np.int32)
    haus = np.empty(stop - start, dtype=np.float64)
    types = np.ones(17, dtype=np.int32)
    templates = _TEMPLATES
    idents = _IDENTS
    shells = _SHELLS
    ntempl = int(templates.shape[0])
    for offset, index in enumerate(range(start, stop)):
        shell = shells[index]
        if not np.isfinite(shell[0, 0]):
            labels[offset] = -1
            haus[offset] = np.nan
            continue
        best_h = 1.0e9
        best = -1
        for slot in range(ntempl):
            _rotation, _translation, _perm, distance = ira.match(
                17, types, templates[slot], 17, types, shell, KMAX_FACTOR)
            distance = float(distance)
            if distance < best_h:
                best_h = distance
                best = int(idents[slot])
            if best_h <= CONGRUENT_A:
                break
        haus[offset] = best_h
        labels[offset] = best if best_h <= HAUSDORFF_A else -1
    return start, labels, haus


def shell_table(points, side: float, threads: int):
    """Centre plus the 16 nearest periodic neighbours inside 5 A."""
    import numpy as np
    import vesin

    n = len(points)
    box = np.diag([side, side, side])
    finder = vesin.NeighborList(cutoff=SHELL_CUTOFF, full_list=True, n_threads=threads)
    left, dist, vec = finder.compute(
        points=points, box=box, periodic=True, quantities="idD", copy=True)
    order = np.argsort(left, kind="mergesort")
    left = left[order]
    dist = dist[order]
    vec = vec[order]
    cuts = np.searchsorted(left, np.arange(n + 1))
    shells = np.full((n, 17, 3), np.nan, dtype=np.float64)
    for index in range(n):
        start = int(cuts[index])
        stop = int(cuts[index + 1])
        if stop - start < 16:
            continue
        block = dist[start:stop]
        pick = np.argpartition(block, 15)[:16]
        if float(block[pick].max()) > SHELL_CUTOFF:
            continue
        pick = pick[np.argsort(block[pick])]
        shells[index, 0] = points[index]
        shells[index, 1:] = points[index] + vec[start:stop][pick]
    return shells


def ira_serial():
    """One process. The library is loaded once, in this process."""
    if _IRA is None:
        _init_ira()
    _start, labels, haus = _match_range((0, len(_SHELLS)))
    return labels, haus


def ira_parallel(threads: int, work: Path):
    """One fresh interpreter per core.

    A fork after vesin or FlexiBLAS has already started threads deadlocks.
    Each child reads the shells through a memory map and writes its slice.
    """
    import numpy as np

    n = len(_SHELLS)
    edges = [int(v) for v in np.linspace(0, n, threads + 1)]
    ranges = [(edges[i], edges[i + 1]) for i in range(threads) if edges[i] < edges[i + 1]]
    procs = []
    for slot, (start, stop) in enumerate(ranges):
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OMP_DYNAMIC"] = "FALSE"
        # close binding would pin every child to the first core.
        env["OMP_PROC_BIND"] = "false"
        env.pop("OMP_PLACES", None)
        cmd = [
            sys.executable, str(Path(__file__).resolve()),
            "--worker", "ira-slice",
            "--work", str(work),
            "--start", str(start),
            "--stop", str(stop),
            "--slot", str(slot),
            "--threads", str(threads),
        ]
        err = open(work / f"slot-{slot}.err", "w", encoding="utf-8")
        procs.append((slot, subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=err), err))
    failed = []
    for slot, proc, err in procs:
        proc.wait()
        err.close()
        if proc.returncode != 0:
            text = (work / f"slot-{slot}.err").read_text(encoding="utf-8", errors="replace")
            failed.append(scrub(text[-500:]))
    if failed:
        raise RuntimeError("IRA slice failed: " + " ".join(failed))
    labels = np.full(n, -1, dtype=np.int32)
    haus = np.full(n, np.nan, dtype=np.float64)
    for slot, (start, stop) in enumerate(ranges):
        labels[start:stop] = np.load(work / f"slot-{slot}-labels.npy")
        haus[start:stop] = np.load(work / f"slot-{slot}-haus.npy")
    return labels, haus


def _require_cores(width: int) -> None:
    if width <= 1:
        return
    try:
        got = os.sched_getaffinity(0)
    except AttributeError:
        return
    if len(got) < width:
        raise RuntimeError(f"affinity has {len(got)} cores, need {width}")


def worker_ira_slice(work: Path, start: int, stop: int, slot: int, threads: int) -> None:
    global _SHELLS, _TEMPLATES, _IDENTS
    import numpy as np

    _require_cores(threads)
    _SHELLS = np.load(work / "shells.npy", mmap_mode="r")
    _TEMPLATES = np.load(work / "templates.npy", mmap_mode="r")
    _IDENTS = np.load(work / "idents.npy")
    _init_ira()
    _require_cores(threads)
    _start, labels, haus = _match_range((start, stop))
    np.save(work / f"slot-{slot}-labels.npy", labels)
    np.save(work / f"slot-{slot}-haus.npy", haus)


def worker_ira(phase: str, threads: int, lib_path: str, agree: str) -> None:
    global _SHELLS, _TEMPLATES, _IDENTS
    import ira_mod
    import numpy as np

    hold_affinity(threads)
    version = str(ira_mod.version)
    xyz, _box, side = frame_of(phase)
    lib = open_core(lib_path)
    _IDENTS, _TEMPLATES = ice_templates(lib)
    _SHELLS = shell_table(xyz, side, threads)
    built = int(np.count_nonzero(np.isfinite(_SHELLS[:, 0, 0])))
    work = Path(agree) / f"ira-{phase}-t{threads}"
    work.mkdir(parents=True, exist_ok=True)
    np.save(work / "shells.npy", _SHELLS)
    np.save(work / "templates.npy", _TEMPLATES)
    np.save(work / "idents.npy", _IDENTS)
    kept = {}

    def run():
        hold_affinity(threads)
        if threads == 1:
            labels, haus = ira_serial()
        else:
            labels, haus = ira_parallel(threads, work)
        kept["labels"] = labels
        kept["haus"] = haus
        return int(labels.shape[0])

    samples, _probe = time_calls(run, WARMUP, REPEATS)
    labels, haus = kept["labels"], kept["haus"]
    save_array(agree, threads, phase, "IRA", "ira", labels)
    save_array(agree, threads, phase, "IRA", "ira-haus", haus)
    import shutil
    shutil.rmtree(work, ignore_errors=True)
    counted = {name: int(np.count_nonzero(labels == code)) for code, name in IRA_NAME.items()}
    finite = haus[np.isfinite(haus)]
    median_h = float(np.median(finite)) if finite.size else float("nan")
    _done(
        "ira", "IRA", version, samples,
        f"ira_mod.IRA.match shells={built}/{MOLECULES} "
        f"hausdorff_A={HAUSDORFF_A} kmax_factor={KMAX_FACTOR} "
        f"congruent_A={CONGRUENT_A} "
        f"median_hausdorff={median_h:.4f} "
        + " ".join(f"{name}={counted[name]}" for name in ("other", "Ih", "Ic", "Isd"))
        + f"; shells built once, match timed, processes={threads}",
    )


def bond_audit(points, side: float, labels) -> dict[str, int]:
    """Four nearest bonds of every d-SEAMS ice label, against the gates."""
    import numpy as np
    import vesin

    ice = np.flatnonzero(labels >= 0)
    if ice.size == 0:
        return {"ice": 0, "window": 0, "isolated": 0}
    box = np.diag([side, side, side])
    finder = vesin.NeighborList(cutoff=BOND_CUTOFF, full_list=True, n_threads=8)
    left, right, dist = finder.compute(
        points=points, box=box, periodic=True, quantities="ijd", copy=True)
    order = np.argsort(left, kind="mergesort")
    left = left[order]
    right = right[order]
    dist = dist[order]
    cuts = np.searchsorted(left, np.arange(len(points) + 1))
    window = 0
    isolated = 0
    for index in ice:
        start = int(cuts[index])
        stop = int(cuts[index + 1])
        block_d = dist[start:stop]
        block_j = right[start:stop]
        keep = block_d > BOND_IGNORE
        block_d = block_d[keep]
        block_j = block_j[keep]
        if block_d.size < 4:
            window += 1
            continue
        pick = np.argpartition(block_d, 3)[:4]
        nearest = block_d[pick]
        if float(nearest.min()) < BOND_LO or float(nearest.max()) > BOND_HI:
            window += 1
            continue
        if not np.any(labels[block_j[pick]] >= 0):
            isolated += 1
    return {"ice": int(ice.size), "window": window, "isolated": isolated}


def crosstab(left, right, left_names: dict, right_names: dict):
    import numpy as np

    rows = []
    for left_code in sorted(left_names):
        mask = left == left_code
        total = int(np.count_nonzero(mask))
        if total == 0:
            continue
        seen = right[mask]
        values, counts = np.unique(seen, return_counts=True)
        for value, count in zip(values, counts):
            rows.append((left_names[left_code], right_names.get(int(value), str(int(value))), int(count)))
    return rows


def _fmt(summary: dict, samples) -> dict:
    return {
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.6f}",
        "iqr_ms": f"{summary['iqr_ms']:.6f}",
        "q1_ms": f"{summary['q1_ms']:.6f}",
        "q3_ms": f"{summary['q3_ms']:.6f}",
        "samples_ms": stats.fmt_samples(samples),
    }


def row_of(msg: dict, threads: int, affinity: str, info: dict, phase: str) -> dict:
    summary = stats.summarize(msg["samples_ms"])
    formed = _fmt(summary, msg["samples_ms"])
    return {
        "source_commit": info["source_commit"],
        "cpu_model": info["cpu_model"],
        "cpu_cores": info["cpu_cores"],
        "affinity": affinity,
        "phase": phase,
        "operation": msg["operation"],
        "code": msg["code"],
        "version": msg["version"],
        "molecules": msg["molecules"],
        "threads": threads,
        "warmup": WARMUP,
        "n": formed["n"],
        "median_ms": formed["median_ms"],
        "iqr_ms": formed["iqr_ms"],
        "q1_ms": formed["q1_ms"],
        "q3_ms": formed["q3_ms"],
        "samples_ms": formed["samples_ms"],
        "detail": msg["detail"],
    }


def launch(code: str, phase: str, threads: int, lib: str, agree: str) -> list[dict]:
    env = stats.pin_env(threads)
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["PYTHONUNBUFFERED"] = "1"
    env["OMP_NUM_THREADS"] = str(threads)
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        sys.executable, str(Path(__file__).resolve()),
        "--worker", code,
        "--phase", phase,
        "--threads", str(threads),
        "--lib", lib,
        "--agree", agree,
    ]
    proc = subprocess.Popen(cmd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    messages = []
    for line in proc.stdout:
        line = line.rstrip("\n")
        if line.startswith("{"):
            messages.append(json.loads(line))
        elif line:
            print(scrub(line), flush=True)
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    code_rc = proc.wait()
    if code_rc != 0:
        raise RuntimeError(scrub(stderr[-2000:] or f"{code} exited {code_rc}"))
    done = [msg for msg in messages if msg.get("event") == "done"]
    if not done:
        raise RuntimeError(f"{code} {phase} threads {threads} returned no row")
    return done


def _note(row: dict, text: str) -> None:
    row["detail"] = f"{row['detail']} {text}".strip()


def paired_fraction(left, right, phase: str) -> str:
    import numpy as np

    if phase == "ice":
        matched = int(np.count_nonzero((left == 1) & (right == 6)))
    else:
        same = (left < 0) & (right == 0)
        same |= (left == 0) & (right == 7)
        same |= (left == 1) & (right == 6)
        matched = int(np.count_nonzero(same))
    return f"{matched}/{len(left)}"


def ira_fraction(labels, ira) -> str:
    import numpy as np

    return f"{int(np.count_nonzero(labels == ira))}/{len(labels)}"


def median_of(values, mask) -> str:
    import numpy as np

    picked = np.asarray(values)[mask]
    picked = picked[np.isfinite(picked)]
    if picked.size == 0:
        return "nan"
    return f"{float(np.median(picked)):.4f}"


def explain(lines: list[str], phase: str, threads: int, labels, ptm, default, ira, haus, ice_rmsd, default_rmsd, rules: bool) -> None:
    import numpy as np

    if rules:
        lines.append(
            "rules: d-SEAMS names a centre from the four nearest neighbours inside 3.5 A "
        "whose lengths sit in [2.30, 3.35] A. Unequal bonds need tetrahedral "
        "order at least 0.85. The second shell is the twelve atoms reached "
        "through those bonds. An equal ideal shell is named by the staggered "
        "and eclipsed fingerprint. Any other shell is named only when the "
        "scale-free residual is within 0.15 and a bonded neighbour is already "
        "crystal; an ice label with no crystal neighbour is removed. "
            "ice Isd is that third template."
        )
        lines.append(
            "OVITO PTM matches the template neighbour star and keeps the lowest "
            f"RMSD at or below {PTM_RMSD}. The default search enables FCC, HCP and "
            "BCC. Cubic and hexagonal diamond stay disabled, so a diamond "
            "neighbourhood is Other in that search. An RMSD of zero on an Other "
            "atom is the unset value when no template was accepted, not a perfect "
            "match. The ice search enables only "
            "those two diamond templates. Diamond is an extension of Larsen, "
            "Schmidt and Schiotz 2016. There is no stacking-disordered template."
        )
        lines.append(
            "IRA matches the centre plus the 16 nearest periodic neighbours, when "
            f"the 16th lies inside {SHELL_CUTOFF:.0f} A, to the three 17-atom ice "
            f"templates. The label is the lowest Hausdorff distance. Above {HAUSDORFF_A} A "
            "the atom is other. kmax_factor is "
            f"{KMAX_FACTOR}. A search stops below {CONGRUENT_A} A because that "
            "template is congruent and no other score can be lower. That stop is "
            "not the classification gate. The Hausdorff distance is not the "
            "scale-free residual, and the 16 nearest atoms are not the bonded "
            "second shell."
        )
    pairs = (
        ("d-SEAMS", DSEAMS_NAME, labels, "PTM-ice", PTM_NAME, ptm),
        ("d-SEAMS", DSEAMS_NAME, labels, "PTM-default", PTM_NAME, default),
        ("d-SEAMS", DSEAMS_NAME, labels, "IRA", IRA_NAME, ira),
        ("IRA", IRA_NAME, ira, "PTM-ice", PTM_NAME, ptm),
    )
    for left_name, left_names, left, right_name, right_names, right in pairs:
        for our, their, count in crosstab(left, right, left_names, right_names):
            lines.append(
                f"crosstab {phase} threads {threads} {left_name} {our} "
                f"{right_name} {their} {count}"
            )
    # Off-diagonal cells of the paired crystal classes, with the scores.
    for our_code, our_name in DSEAMS_NAME.items():
        mask_our = labels == our_code
        if not np.any(mask_our):
            continue
        for ptm_code, ptm_name in PTM_NAME.items():
            both = mask_our & (ptm == ptm_code)
            count = int(np.count_nonzero(both))
            if count == 0:
                continue
            expected = PTM_OF.get(our_code)
            if expected == ptm_code or (our_code == 2 and ptm_code in (0, 6, 7)):
                kind = "paired" if expected == ptm_code else "unpaired-Isd"
            else:
                kind = "differ"
            lines.append(
                f"score {phase} threads {threads} d-SEAMS {our_name} PTM-ice {ptm_name} "
                f"{kind} n={count} median_rmsd={median_of(ice_rmsd, both)} "
                f"median_hausdorff={median_of(haus, both)}"
            )
        for ira_code, ira_name in IRA_NAME.items():
            both = mask_our & (ira == ira_code)
            count = int(np.count_nonzero(both))
            if count == 0 or ira_code == our_code:
                continue
            lines.append(
                f"score {phase} threads {threads} d-SEAMS {our_name} IRA {ira_name} "
                f"differ n={count} median_hausdorff={median_of(haus, both)} "
                f"median_ptm_rmsd={median_of(ice_rmsd, both)}"
            )
    for ptm_code, ptm_name in ((1, "FCC"), (2, "HCP"), (3, "BCC"), (0, "Other")):
        both = default == ptm_code
        count = int(np.count_nonzero(both))
        if count == 0:
            continue
        lines.append(
            f"score {phase} threads {threads} PTM-default {ptm_name} "
            f"n={count} median_rmsd={median_of(default_rmsd, both)}"
        )
    n = len(labels)
    finite = haus[np.isfinite(haus)]
    median_h = float(np.median(finite)) if finite.size else float("nan")
    if (
        int(np.count_nonzero(labels == 1)) == n
        and int(np.count_nonzero(ptm == 6)) == n
        and int(np.count_nonzero(ira == 1)) == n
        and int(np.count_nonzero(default == 0)) == n
    ):
        lines.append(
            f"reading {phase} threads {threads}: every atom is ice Ic, "
            "cubic diamond, and the IRA ice Ic template. "
            f"agree_paired={n}/{n} agree_ira={n}/{n}. "
            "The default PTM search names every atom Other because cubic "
            "diamond is disabled. That is the template set."
        )
    elif (
        int(np.count_nonzero(labels < 0)) == n
        and int(np.count_nonzero(ptm == 0)) == n
        and int(np.count_nonzero(ira < 0)) == n
        and int(np.count_nonzero(default == 0)) == n
    ):
        lines.append(
            f"reading {phase} threads {threads}: every atom is other. "
            f"agree_paired={n}/{n} agree_ira={n}/{n}. "
            f"shells={int(finite.size)}/{n} median_hausdorff={median_h:.4f}. "
            "No centre passes the bond gates, no diamond template is at or "
            f"below RMSD {PTM_RMSD}, and no shell is at or below Hausdorff "
            f"{HAUSDORFF_A} A. The three rules agree because this liquid fails "
            "each test on its own terms."
        )


def write_table(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_sitting(args) -> int:
    import numpy as np

    lib = args.lib or str(ROOT / "bbdir-release" / "src" / "libyodaLib.so")
    if not Path(lib).is_file():
        raise SystemExit(f"missing library {lib}")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    agree = out / ".ptm-ira-agree"
    agree.mkdir(parents=True, exist_ok=True)
    info = stats.machine(ROOT)
    log = [
        *stats.machine_lines(ROOT),
        "lscpu",
        *large.lscpu_details(),
        "sitting one",
        f"warmup {WARMUP}",
        f"n {REPEATS}",
        f"molecules {MOLECULES}",
        "ice Ic 48^3 conventional cells",
        "liquid TIP4P/ice 4096 molecules replicated 6^3",
        "widths 1 8",
        "one sitting; each width is one process",
        "host name is not recorded",
        f"ptm_rmsd {PTM_RMSD}",
        f"ira_hausdorff_A {HAUSDORFF_A}",
        f"ira_kmax_factor {KMAX_FACTOR}",
        "gates are not adjusted",
    ]
    collected = []
    for phase in large.PHASES:
        for threads in WIDTHS:
            affinity = stats.taskset_cores(threads)
            for code in ("dseams", "ovito", "ira"):
                print(f"phase {phase} threads {threads} code {code}", flush=True)
                messages = launch(code, phase, threads, lib, str(agree))
                for msg in messages:
                    collected.append(row_of(msg, threads, affinity, info, phase))
                    print(
                        f"done {phase} {msg['operation']} {msg['code']} "
                        f"threads {threads} median_ms "
                        f"{stats.summarize(msg['samples_ms'])['median_ms']:.3f}",
                        flush=True,
                    )
                    log.append(
                        f"{phase} {msg['operation']} {msg['code']} threads {threads} "
                        f"n {len(msg['samples_ms'])} "
                        f"median_ms {stats.summarize(msg['samples_ms'])['median_ms']:.6f} "
                        f"affinity {affinity} {msg['detail']}"
                    )
    by = {(row["phase"], str(row["threads"]), row["code"], row["operation"]): row for row in collected}
    agree_rows = []
    problems = []
    for phase in large.PHASES:
        xyz, _box, side = frame_of(phase)
        reference = None
        for threads in WIDTHS:
            labels = load_array(agree, threads, phase, "d-SEAMS", "classifier")
            ptm = load_array(agree, threads, phase, "OVITO", "ptm-ice")
            default = load_array(agree, threads, phase, "OVITO", "ptm-default")
            ira = load_array(agree, threads, phase, "IRA", "ira")
            haus = load_array(agree, threads, phase, "IRA", "ira-haus")
            ice_rmsd = load_array(agree, threads, phase, "OVITO", "ptm-ice-rmsd")
            default_rmsd = load_array(agree, threads, phase, "OVITO", "ptm-default-rmsd")
            if any(item is None for item in (labels, ptm, default, ira, haus, ice_rmsd, default_rmsd)):
                raise SystemExit(f"missing labels for {phase} threads {threads}")
            if reference is None:
                reference = (labels, ptm, default, ira)
                audit = bond_audit(xyz, side, labels)
                log.append(
                    f"audit {phase} d-SEAMS ice={audit['ice']} "
                    f"bond_window_violations={audit['window']} "
                    f"isolated_violations={audit['isolated']}"
                )
                if audit["window"] or audit["isolated"]:
                    problems.append(
                        f"{phase} d-SEAMS ice labels break the bond gates "
                        f"window={audit['window']} isolated={audit['isolated']}"
                    )
            else:
                names = ("d-SEAMS", "PTM-ice", "PTM-default", "IRA")
                for name, got, ref in zip(names, (labels, ptm, default, ira), reference):
                    if not np.array_equal(got, ref):
                        problems.append(f"{phase} {name} labels differ between 1 and 8 threads")
            text_ptm = "agree_paired=" + paired_fraction(labels, ptm, phase)
            text_ira = "agree_ira=" + ira_fraction(labels, ira)
            if phase == "ice":
                ic = int(np.count_nonzero(labels == 1))
                cubic = int(np.count_nonzero(ptm == 6))
                other = int(np.count_nonzero(default == 0))
                ira_ic = int(np.count_nonzero(ira == 1))
                text_ice = (
                    f"ice_Ic={ic}/{MOLECULES} ptm_cubic_diamond={cubic}/{MOLECULES} "
                    f"ptm_default_Other={other}/{MOLECULES} ira_Ic={ira_ic}/{MOLECULES}"
                )
                if ic != MOLECULES or cubic != MOLECULES or ira_ic != MOLECULES:
                    problems.append(f"ice threads {threads} {text_ice}")
                if other != MOLECULES:
                    problems.append(
                        f"ice default PTM threads {threads} Other {other}/{MOLECULES}"
                    )
                log.append(f"identity {phase} threads {threads} {text_ice}")
            else:
                crystal_d = int(np.count_nonzero(labels >= 0))
                crystal_p = int(np.count_nonzero((ptm == 6) | (ptm == 7)))
                crystal_i = int(np.count_nonzero(ira >= 0))
                log.append(
                    f"identity {phase} threads {threads} "
                    f"crystal_dseams={crystal_d} crystal_ptm={crystal_p} crystal_ira={crystal_i}"
                )
            for code, operation, text in (
                ("d-SEAMS", "classifier", text_ptm + " " + text_ira),
                ("OVITO", "ptm-ice", text_ptm),
                ("IRA", "ira", text_ira),
            ):
                row = by[(phase, str(threads), code, operation)]
                _note(row, text)
            blocks = (
                ("d-SEAMS", DSEAMS_NAME, labels, "PTM-ice", PTM_NAME, ptm),
                ("d-SEAMS", DSEAMS_NAME, labels, "PTM-default", PTM_NAME, default),
                ("d-SEAMS", DSEAMS_NAME, labels, "IRA", IRA_NAME, ira),
                ("IRA", IRA_NAME, ira, "PTM-ice", PTM_NAME, ptm),
            )
            for left_name, left_names, left, right_name, right_names, right in blocks:
                for our, their, count in crosstab(left, right, left_names, right_names):
                    agree_rows.append({
                        "phase": phase,
                        "threads": threads,
                        "left": left_name,
                        "left_label": our,
                        "right": right_name,
                        "right_label": their,
                        "atoms": count,
                        "molecules": MOLECULES,
                        "fraction": f"{count / MOLECULES:.8f}",
                    })
            explain(
                log, phase, threads, labels, ptm, default, ira, haus,
                ice_rmsd, default_rmsd, rules=(threads == 1),
            )
    if problems:
        log.append("failures")
        log.extend(problems)
    (out / "raw-ptm-ira.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    write_table(out / "sitting-ptm-ira.csv", FIELDS, collected)
    write_table(out / "ptm-ira-agreement.csv", AGREE_FIELDS, agree_rows)
    # Label arrays are working files. The tables are the record.
    for path in agree.glob("*.npy"):
        path.unlink()
    if problems:
        raise SystemExit("structure labels break a published gate\n" + "\n".join(problems))
    print(f"wrote {out / 'sitting-ptm-ira.csv'}", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="benchmarks/paper-cells")
    parser.add_argument("--lib", default="")
    parser.add_argument("--worker", choices=("dseams", "ovito", "ira", "ira-slice"))
    parser.add_argument("--phase", choices=large.PHASES)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--agree", default="")
    parser.add_argument("--work", default="")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=0)
    parser.add_argument("--slot", type=int, default=0)
    args = parser.parse_args()
    if args.worker == "ira-slice":
        if not args.work or args.stop <= args.start:
            raise SystemExit("ira-slice needs --work, --start and --stop")
        worker_ira_slice(Path(args.work), args.start, args.stop, args.slot, args.threads)
        return 0
    if args.worker:
        if not args.phase or not args.agree:
            raise SystemExit("a worker needs --phase and --agree")
        lib = args.lib or str(ROOT / "bbdir-release" / "src" / "libyodaLib.so")
        if args.worker == "dseams":
            worker_dseams(args.phase, args.threads, lib, args.agree)
        elif args.worker == "ovito":
            worker_ovito(args.phase, args.threads, args.agree)
        else:
            worker_ira(args.phase, args.threads, lib, args.agree)
        return 0
    return run_sitting(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        sys.stderr.write(scrub(str(exc)) + "\n")
        raise
