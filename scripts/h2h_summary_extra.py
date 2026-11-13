#!/usr/bin/env python3
"""Competitors named by the paper, in the summary sitting.

Each timed row asserts the outputs match. An operation the library does
not provide is absent. A sample that does not return inside the wall-clock
budget is DNF. Adaptive CNA is a label comparison, not a timing.
"""

from __future__ import annotations

import ctypes
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BUDGET_S = 600
STAMP = ROOT / "bbdir-release" / "h2h-budget.json"
PTM_REV = "162cc834"

PLUMED_KERNEL = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/PLUMED/2.10.0-foss-2026.1/lib/libplumedKernel.so"
)
Q6_SCALE = (4.0 * 3.14159265358979323846 / 13.0) ** 0.5


def absent(row_of, competitor, operation, molecules, threads, cores, model, version, why):
    made = row_of(competitor, operation, molecules, threads, cores, [1.0], [1.0], why, version, model)
    made["n"] = 0
    made["dseams_median_ms"] = ""
    made["dseams_min_ms"] = ""
    made["other_median_ms"] = ""
    made["other_min_ms"] = ""
    made["ratio"] = ""
    made["verdict"] = "absent"
    made["budget_s"] = ""
    made["largest_finished"] = ""
    return made


def blank_times(made, verdict):
    made["n"] = 0
    made["dseams_median_ms"] = ""
    made["dseams_min_ms"] = ""
    made["other_median_ms"] = ""
    made["other_min_ms"] = ""
    made["ratio"] = ""
    made["verdict"] = verdict
    return made


def dnf(row_of, competitor, operation, molecules, threads, cores, model, version, largest, why):
    made = row_of(competitor, operation, molecules, threads, cores, [1.0], [1.0], why, version, model)
    blank_times(made, "DNF")
    made["budget_s"] = str(BUDGET_S)
    made["largest_finished"] = "" if not largest else str(largest)
    return made


def correctness(row_of, competitor, operation, molecules, threads, cores, model, version, text):
    made = row_of(competitor, operation, molecules, threads, cores, [1.0], [1.0], text, version, model)
    blank_times(made, "correctness")
    made["n"] = 1
    made["budget_s"] = ""
    made["largest_finished"] = ""
    return made


def _stamp():
    if STAMP.is_file():
        return json.loads(STAMP.read_text(encoding="ascii"))
    return {}


def _save_stamp(data):
    STAMP.parent.mkdir(parents=True, exist_ok=True)
    STAMP.write_text(json.dumps(data), encoding="ascii")


def note_finished(kind, molecules):
    data = _stamp()
    key = f"finished:{kind}"
    data[key] = max(int(data.get(key, 0)), int(molecules))
    _save_stamp(data)


def largest_finished(kind):
    return int(_stamp().get(f"finished:{kind}", 0))


def smaller_timeout(kind, molecules, threads):
    serial = kind in {"cycless", "ira", "plumed"}
    width = 1 if serial else threads
    found = 0
    for key, rec in _stamp().items():
        parts = key.split(":")
        if len(parts) != 3 or parts[0] != kind or parts[2] != str(width):
            continue
        if int(parts[1]) < molecules and not rec["ok"]:
            found = max(found, int(parts[1]))
    return found


def budget_probe(kind, molecules, threads):
    """One sample. Serial kinds share the threads=1 stamp."""
    serial = kind in {"cycless", "ira", "plumed"}
    width = 1 if serial else threads
    key = f"{kind}:{molecules}:{width}"
    data = _stamp()
    if key in data:
        return data[key]
    print(f"probe {kind} molecules {molecules} threads {width} budget {BUDGET_S}", file=sys.stderr, flush=True)
    cmd = [sys.executable, "-u", str(ROOT / "scripts" / "h2h_budget.py"), kind, str(molecules), str(width)]
    t0 = time.perf_counter()
    try:
        subprocess.run(cmd, cwd=ROOT, env=os.environ.copy(), timeout=BUDGET_S, check=True)
        rec = {"ok": True, "seconds": round(time.perf_counter() - t0, 3)}
    except subprocess.TimeoutExpired:
        rec = {"ok": False, "seconds": BUDGET_S}
    data = _stamp()
    data[key] = rec
    _save_stamp(data)
    print(f"probe {kind} molecules {molecules} ok {rec['ok']} seconds {rec['seconds']}", file=sys.stderr, flush=True)
    return rec


def bind(lib):
    if not hasattr(lib, "seams_system_symmetry") or lib.seams_system_symmetry.argtypes is None:
        lib.seams_system_symmetry.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ]
        lib.seams_system_symmetry.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int


def ring_hist(lib, held, h):
    rings = ctypes.c_void_p()
    h.require(lib, lib.seams_rings_find(held, h.CUTOFF, 1, 0, 6, ctypes.byref(rings)), "rings")
    hist = (ctypes.c_int * 7)()
    h.require(lib, lib.seams_rings_histogram(rings, hist, 7), "histogram")
    lib.seams_rings_destroy(rings)
    return [int(hist[k]) for k in range(7)]


def load_ptm():
    lib = ctypes.CDLL(str(ROOT / "bbdir-release" / "libptm_ice.so"))
    lib.ptm_ice_classify.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.ptm_ice_classify.restype = ctypes.c_int
    return lib


def ptm_codes(ptm, pos, lengths, threads):
    n = int(pos.shape[0])
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    box = (ctypes.c_double * 3)(*lengths)
    types = (ctypes.c_int * n)()
    rmsd = (ctypes.c_double * n)()
    code = ptm.ptm_ice_classify(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n, box, 0.1, threads, types, rmsd,
    )
    if code != 0:
        raise SystemExit(f"ptm_ice_classify returned {code}")
    return np.ctypeslib.as_array(types).copy(), float(np.max(np.ctypeslib.as_array(rmsd)))


def ptm_rows(lib, held, pos, lengths, threads, cores, model, version, row_of, h, pin, abab):
    n = int(pos.shape[0])
    if n > 110592:
        rec = budget_probe("ptm-lib", n, threads)
        if not rec["ok"]:
            return [dnf(
                row_of, "PTM", "ptm", n, threads, cores, model, f"{version}/{PTM_REV}",
                largest_finished("ptm"),
                f"one sample exceeded {BUDGET_S} s",
            )]
    ptm = load_ptm()
    codes, rmsd_max = ptm_codes(ptm, pos, lengths, threads)
    sym = (ctypes.c_int * n)()
    h.require(lib, lib.seams_set_num_threads(threads), "threads")
    h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")
    sym_np = np.ctypeslib.as_array(sym)
    # d-SEAMS symmetry 0 is hexagonal diamond. The library's 7 is the same label.
    if not np.all(sym_np == 0) or not np.all(codes == 7) or not (rmsd_max <= 0.1):
        raise SystemExit(
            f"PTM labels differ at {n}: symmetry {int((sym_np == 0).sum())} "
            f"library {int((codes == 7).sum())} rmsd_max {rmsd_max:.3e}"
        )

    def ours():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")

    def other():
        pin(threads)
        ptm_codes(ptm, pos, lengths, threads)

    left, right = abab(ours, other)
    note_finished("ptm", n)
    return [row_of(
        "PTM", "ptm", n, threads, cores, left, right,
        f"hexagonal {n} equal rmsd_max {rmsd_max:.3e}",
        f"{version}/{PTM_REV}", model,
    )]


def cna_codes(pos, lengths):
    from ovito.data import DataCollection
    from ovito.modifiers import CommonNeighborAnalysisModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(np.diag(lengths), pbc=(True, True, True))
    parts = data.create_particles(count=len(pos))
    parts.create_property("Position", data=np.ascontiguousarray(pos))
    mod = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(mod)
    return np.asarray(pipe.compute().particles["Structure Type"])


def chill_names(lib, h, pos, lengths):
    n = int(pos.shape[0])
    held = h.frame_system(lib, pos, lengths)
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    h.dseams_chill_plus(lib, held, labels, counts, h.SEAMS_CHILL_OVITO)
    lib.seams_system_destroy(held)
    return [h.DSEAMS_AS_OVITO[int(v)] for v in labels]


def mismatch_text(phase, names, codes):
    crystal = np.array([name != "other" for name in names])
    other = np.asarray(codes) == 0
    rate = float(np.mean(crystal == other))
    # crystal==other means CHILL+ says ice and CNA says Other, or the reverse
    # of that equality is agreement. The rate above is the mismatch.
    chill = {}
    for name in names:
        chill[name] = chill.get(name, 0) + 1
    cna_name = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}
    cna = {}
    for value in np.asarray(codes).tolist():
        label = cna_name.get(int(value), str(int(value)))
        cna[label] = cna.get(label, 0) + 1
    chill_s = " ".join(f"{key}={chill[key]}" for key in sorted(chill))
    cna_s = " ".join(f"{key}={cna[key]}" for key in sorted(cna))
    return f"phase {phase} mismatch_rate {rate:.6f} chill {chill_s} cna {cna_s}"


def liquid_positions(h):
    pos, lengths = h.ice_positions(8)
    rng = np.random.Generator(np.random.PCG64(20261113))
    moved = np.ascontiguousarray(pos + rng.normal(scale=1.5, size=pos.shape))
    box = np.asarray(lengths, dtype=np.float64)
    return np.ascontiguousarray(moved % box), lengths


def cna_rows(lib, h, threads, cores, model, version, row_of, ih_pos, ih_lengths):
    from gen_ice_ic import ice_ic

    ic_pos, side = ice_ic(8)
    frames = (
        ("Ih", ih_pos, ih_lengths),
        ("Ic", np.ascontiguousarray(ic_pos), (side, side, side)),
        ("liquid", *liquid_positions(h)),
    )
    rows = []
    for phase, pos, lengths in frames:
        names = chill_names(lib, h, pos, lengths)
        codes = cna_codes(pos, lengths)
        rows.append(correctness(
            row_of, "OVITO", "cna", int(len(pos)), threads, cores, model, version,
            mismatch_text(phase, names, codes),
        ))
    return rows


def sqq_once(pos, lengths, threads):
    import numpy as np
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame

    pos_nm = np.asarray(pos, dtype=np.float64) / 10.0
    box = np.asarray(lengths, dtype=np.float64) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(len(pos_nm))
    ]
    frame = Frame(name=f"ice-{len(atoms)}", atoms=atoms, box=box)
    cfg = load_config({
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
        "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": True},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": threads, "math_threads": 1},
    }, engine="py")
    result = analyze_frame(frame, cfg, engine="py")
    return len(result.cages)


def dseams_cage_total(lib, held, h):
    if lib.seams_cages_find.argtypes is None:
        lib.seams_cages_find.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p),
        ]
        lib.seams_cages_find.restype = ctypes.c_int
        lib.seams_cages_destroy.argtypes = [ctypes.c_void_p]
        lib.seams_cages_counts.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
        ]
        lib.seams_cages_counts.restype = ctypes.c_int
    cages = ctypes.c_void_p()
    h.require(lib, lib.seams_cages_find(held, h.CUTOFF, 1, 0, 0, 0, ctypes.byref(cages)), "cages")
    hex_n = ctypes.c_int()
    cubic_n = ctypes.c_int()
    water_n = ctypes.c_int()
    both_n = ctypes.c_int()
    h.require(lib, lib.seams_cages_counts(
        cages, ctypes.byref(hex_n), ctypes.byref(cubic_n),
        ctypes.byref(water_n), ctypes.byref(both_n),
    ), "cage counts")
    lib.seams_cages_destroy(cages)
    return int(hex_n.value + cubic_n.value + water_n.value + both_n.value)


def sqq_rows(lib, held, pos, lengths, threads, cores, model, version, row_of, h, pin, abab):
    n = int(pos.shape[0])
    blocked = smaller_timeout("sqq", n, threads)
    if blocked:
        return [dnf(
            row_of, "SQQ", "cages", n, threads, cores, model, "0.5.8",
            largest_finished("sqq"),
            f"{blocked} molecules exceeded {BUDGET_S} s, so this size was not started",
        )]
    rec = budget_probe("sqq", n, threads)
    if not rec["ok"]:
        if largest_finished("sqq") == 0:
            small = budget_probe("sqq", 512, threads)
            if small["ok"]:
                note_finished("sqq", 512)
        return [dnf(
            row_of, "SQQ", "cages", n, threads, cores, model, "0.5.8",
            largest_finished("sqq"),
            f"one sample exceeded {BUDGET_S} s",
        )]
    ours_n = dseams_cage_total(lib, held, h)
    other_n = sqq_once(pos, lengths, threads)
    if ours_n != 0 or other_n != 0:
        raise SystemExit(f"ice cages d-SEAMS {ours_n} SQQ {other_n} at {n}")

    def ours():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        dseams_cage_total(lib, held, h)

    def other():
        pin(threads)
        sqq_once(pos, lengths, threads)

    left, right = abab(ours, other)
    note_finished("sqq", n)
    return [row_of(
        "SQQ", "cages", n, threads, cores, left, right,
        f"cages 0 equal", f"{version}/0.5.8", model,
    )]


def extra_rows(lib, held, pos, lengths, threads, cores, model, seams, pin, abab, row_of, h):
    n = int(pos.shape[0])
    rows = []
    version = seams

    if threads == 1 and n == 4096:
        rows.extend(cna_rows(lib, h, threads, cores, model, version, row_of, pos, lengths))
    rows.extend(ptm_rows(
        lib, held, pos, lengths, threads, cores, model, version, row_of, h, pin, abab,
    ))
    rows.append(absent(
        row_of, "ClavisCage", "chill_plus", n, threads, cores, model, version,
        "ClavisCage is not installed; the clavis package here is unrelated",
    ))
    rows.append(absent(
        row_of, "freud", "rings", n, threads, cores, model, version,
        "freud has no primitive-ring finder",
    ))
    from ovito.data import DataCollection
    from ovito.modifiers import (
        CreateBondsModifier,
        FindRingsModifier,
        IdentifyDiamondModifier,
        PolyhedralTemplateMatchingModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().setMaxThreadCount(threads)
    data = DataCollection()
    data.create_cell(np.diag(lengths), pbc=(True, True, True))
    parts = data.create_particles(count=n)
    parts.create_property("Position", data=np.ascontiguousarray(pos))

    sym = (ctypes.c_int * n)()
    h.require(lib, lib.seams_set_num_threads(threads), "threads")
    h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")
    sym_np = np.ctypeslib.as_array(sym).copy()

    def fresh_pipe(*mods):
        pipe = Pipeline(source=StaticSource(data=data))
        for mod in mods:
            pipe.modifiers.append(mod)
        return pipe

    # 884736 is attempted under the wall-clock budget. A smaller timeout
    # means this size is not started again.
    do_rings = n <= 110592
    do_ptm = n <= 110592
    if n > 110592:
        if smaller_timeout("ovito-rings", n, threads):
            do_rings = False
        else:
            do_rings = budget_probe("ovito-rings", n, threads)["ok"]
        if smaller_timeout("ovito-ptm", n, threads):
            do_ptm = False
        else:
            do_ptm = budget_probe("ovito-ptm", n, threads)["ok"]
        if not do_rings:
            why = (
                f"{smaller_timeout('ovito-rings', n, threads)} molecules exceeded {BUDGET_S} s, so this size was not started"
                if smaller_timeout("ovito-rings", n, threads)
                else f"one sample exceeded {BUDGET_S} s"
            )
            rows.append(dnf(
                row_of, "OVITO", "rings", n, threads, cores, model, version,
                largest_finished("ovito-rings") or 110592, why,
            ))
        if not do_ptm:
            why = (
                f"{smaller_timeout('ovito-ptm', n, threads)} molecules exceeded {BUDGET_S} s, so this size was not started"
                if smaller_timeout("ovito-ptm", n, threads)
                else f"one sample exceeded {BUDGET_S} s"
            )
            rows.append(dnf(
                row_of, "OVITO", "ptm", n, threads, cores, model, version,
                largest_finished("ovito-ptm") or 110592, why,
            ))
    ours_hist = ring_hist(lib, held, h)
    ovito_hex = ours_hist[6]
    ring_pipe = None
    ptm_pipe = None
    if do_rings:
        bonds = CreateBondsModifier(cutoff=h.CUTOFF)
        rings_mod = FindRingsModifier(maximum_ring_size=6)
        ring_pipe = fresh_pipe(bonds, rings_mod)
        table = ring_pipe.compute().tables["ring-size-histogram"]
        ovito_hex = 0
        for size, count in zip(table["Ring Size"][...], table["Count"][...]):
            if int(size) == 6:
                ovito_hex = int(count)
        if ours_hist[6] != ovito_hex:
            raise SystemExit(f"OVITO rings {ovito_hex} differ from d-SEAMS {ours_hist[6]} at {n}")
    if do_ptm:
        ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
        ptm_pipe = fresh_pipe(ptm)
        ptm_lab = np.asarray(ptm_pipe.compute().particles["Structure Type"])
        # d-SEAMS symmetry 0 is hexagonal. OVITO PTM 7 is hexagonal diamond.
        if not np.all(sym_np == 0) or not np.all(ptm_lab == 7):
            raise SystemExit(f"PTM labels do not match hexagonal ice at {n}")

    dia = IdentifyDiamondModifier()
    dia_pipe = fresh_pipe(dia)
    dia_lab = np.asarray(dia_pipe.compute().particles["Structure Type"])
    # Identify diamond: 4 hexagonal, 5 first shell, 6 second shell.
    if not np.all(np.isin(dia_lab, (4, 5, 6))):
        raise SystemExit(f"identify-diamond is not hexagonal at {n}")

    def dseams_rings():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        ring_hist(lib, held, h)

    def ovito_rings():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        for mod in ring_pipe.modifiers:
            mod.notify_object_changed()
        ring_pipe.compute()

    def dseams_sym():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")

    def ovito_ptm():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        ptm_pipe.modifiers[0].notify_object_changed()
        ptm_pipe.compute()

    def ovito_dia():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        dia_pipe.modifiers[0].notify_object_changed()
        dia_pipe.compute()

    import ovito
    ov = ".".join(str(v) for v in ovito.version)
    jobs = [("OVITO", "identify_diamond", dseams_sym, ovito_dia, f"hexagonal {n} equal", "")]
    if do_rings:
        jobs.append(("OVITO", "rings", dseams_rings, ovito_rings, f"hexagons {ovito_hex} equal", "ovito-rings"))
    if do_ptm:
        jobs.append(("OVITO", "ptm", dseams_sym, ovito_ptm, f"hexagonal {n} equal", "ovito-ptm"))
    for name, op, ours, other, note, finished in jobs:
        left, right = abab(ours, other)
        rows.append(row_of(
            name, op, n, threads, cores, left, right, note, f"{version}/{ov}", model,
        ))
        if finished:
            note_finished(finished, n)

    # PLUMED Q6, scaled to the local Steinhardt parameter.
    run_plumed = os.path.isfile(PLUMED_KERNEL)
    if run_plumed and n > 110592:
        blocked = smaller_timeout("plumed", n, threads)
        run_plumed = not blocked and budget_probe("plumed", n, threads)["ok"]
        if not run_plumed:
            rows.append(dnf(
                row_of, "PLUMED", "q6", n, threads, cores, model, version,
                largest_finished("plumed") or 110592,
                f"{blocked} molecules exceeded {BUDGET_S} s, so this size was not started"
                if blocked else f"one sample exceeded {BUDGET_S} s",
            ))
    if run_plumed:
        os.environ.setdefault("PLUMED_KERNEL", PLUMED_KERNEL)
        import plumed

        q = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()
        h.dseams_q6(lib, held, q, bar)
        local = np.ctypeslib.as_array(q).copy()
        line = f"obs: Q6 SPECIES=1-{n} SWITCH={{CUSTOM FUNC=step(3.5-x) R_0=1 D_MAX=4.2}}"
        calc = plumed.Plumed()
        calc.cmd("setMDEngine", "python")
        calc.cmd("setNatoms", n)
        calc.cmd("setTimestep", 1.0)
        calc.cmd("setKbT", 1.0)
        calc.cmd("setLogFile", "/dev/null")
        calc.cmd("init")
        calc.cmd("readInputLine", line)
        rank = np.zeros(1, dtype=np.int64)
        calc.cmd("getDataRank obs", rank)
        shape = np.zeros(int(rank[0]), dtype=np.int64)
        calc.cmd("getDataShape obs", shape)
        data = np.zeros(tuple(int(v) for v in shape), dtype=np.float64)
        calc.cmd("setMemoryForData obs", data)
        box = np.diag(np.asarray(lengths, dtype=np.float64))
        mass = np.ones(n, dtype=np.float64)
        charge = np.zeros(n, dtype=np.float64)
        forces = np.zeros((n, 3), dtype=np.float64)
        virial = np.zeros((3, 3), dtype=np.float64)
        step = {"n": 0}

        def plumed_q():
            step["n"] += 1
            calc.cmd("setStep", step["n"])
            calc.cmd("setBox", box)
            calc.cmd("setPositions", pos)
            calc.cmd("setMasses", mass)
            calc.cmd("setCharges", charge)
            calc.cmd("setForces", forces)
            calc.cmd("setVirial", virial)
            calc.cmd("calc")

        plumed_q()
        got = np.array(data, copy=True).reshape(-1) * Q6_SCALE
        diff = float(np.max(np.abs(local - got)))
        if diff > 1e-8:
            raise SystemExit(f"PLUMED q6 max_abs {diff:.3e} at {n}")

        def dseams_local():
            pin(threads)
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            h.dseams_q6(lib, held, q, bar)

        left, right = abab(dseams_local, plumed_q)
        rows.append(row_of(
            "PLUMED", "q6", n, threads, cores, left, right,
            f"q6_max_abs {diff:.3e}", f"{version}/2.10.0", model,
        ))
        note_finished("plumed", n)
    elif not os.path.isfile(PLUMED_KERNEL):
        rows.append(absent(
            row_of, "PLUMED", "q6", n, threads, cores, model, version,
            "libplumedKernel.so was not found",
        ))

    # rings-statistics-matter, the Rust shortest-path census.
    import rings_statistics_matter
    from ase import Atoms

    atoms = Atoms(symbols=["O"] * n, positions=pos, cell=list(lengths), pbc=True)

    def rsm_counts():
        _counts, blocks = rings_statistics_matter.ring_statistics(
            atoms, cutoff=h.CUTOFF, maxlvl=8, supercell=False,
            deduplicate=True, backend="rust", num_threads=threads,
        )
        found = [0] * 7
        for block in blocks:
            arr = np.asarray(block)
            if arr.size == 0 or arr.ndim != 2:
                continue
            size = int(arr.shape[1])
            if 0 <= size <= 6:
                found[size] = int(arr.shape[0])
        return found

    rsm = rsm_counts()
    if rsm[6] != ours_hist[6]:
        raise SystemExit(f"rings-statistics-matter {rsm[6]} hexagons vs d-SEAMS {ours_hist[6]}")

    def dseams_r():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        ring_hist(lib, held, h)

    def rsm_call():
        pin(threads)
        rsm_counts()

    left, right = abab(dseams_r, rsm_call)
    rows.append(row_of(
        "rings-statistics-matter", "rings", n, threads, cores, left, right,
        f"hexagons {rsm[6]} equal", f"{version}/0.3.0", model,
    ))

    run_cycless = n <= 4096 or (
        not smaller_timeout("cycless", n, threads)
        and budget_probe("cycless", n, threads)["ok"]
    )
    run_ira = n <= 4096 or (
        not smaller_timeout("ira", n, threads)
        and budget_probe("ira", n, threads)["ok"]
    )
    if not run_cycless:
        blocked = smaller_timeout("cycless", n, threads)
        rows.append(dnf(
            row_of, "cycless", "rings", n, threads, cores, model, "0.7",
            largest_finished("cycless") or 4096,
            f"{blocked} molecules exceeded {BUDGET_S} s, so this size was not started"
            if blocked else f"one sample exceeded {BUDGET_S} s",
        ))
    if not run_ira:
        blocked = smaller_timeout("ira", n, threads)
        rows.append(dnf(
            row_of, "IRA", "shape", n, threads, cores, model, "2.2.0",
            largest_finished("ira") or 4096,
            f"{blocked} molecules exceeded {BUDGET_S} s, so this size was not started"
            if blocked else f"one sample exceeded {BUDGET_S} s",
        ))
    if run_cycless or run_ira:
        import networkx as nx
        from cycless.cycles import cycles_iter
        import freud

        box = freud.box.Box(*lengths)
        query = freud.locality.AABBQuery(box, pos).query(
            pos, {"r_max": h.CUTOFF, "exclude_ii": True}
        ).toNeighborList()
        graph = nx.Graph()
        graph.add_nodes_from(range(n))
        ii = np.asarray(query[:, 0])
        jj = np.asarray(query[:, 1])
        graph.add_edges_from(
            (int(a), int(b)) for a, b in zip(ii, jj) if int(a) < int(b)
        )
        if run_cycless:
            cyc = sum(1 for _ in cycles_iter(graph, 6))
            if cyc != ours_hist[6]:
                raise SystemExit(f"cycless {cyc} rings vs d-SEAMS {ours_hist[6]}")

            def cyc_call():
                pin(threads)
                sum(1 for _ in cycles_iter(graph, 6))

            left, right = abab(dseams_r, cyc_call)
            rows.append(row_of(
                "cycless", "rings", n, threads, cores, left, right,
                f"hexagons {cyc} equal; cycles_iter has no thread argument",
                f"{version}/0.7", model,
            ))
            note_finished("cycless", n)

        if not run_ira:
            rows.extend(sqq_rows(
                lib, held, pos, lengths, threads, cores, model, version, row_of, h, pin, abab,
            ))
            return rows

        import ira_mod

        ira = ira_mod.IRA()
        from h2h_ira import classify_serial, cluster

        ih = cluster(pos, np.asarray(lengths), 0)
        ic_pos, side = __import__("gen_ice_ic", fromlist=["ice_ic"]).ice_ic(4)
        ic = cluster(
            np.ascontiguousarray(ic_pos),
            np.asarray((side, side, side)),
            0,
        )
        labels_ira = classify_serial(ira, ih, ic, pos, np.asarray(lengths))
        # Ideal Ih is hexagonal under symmetry. IRA's nearer template is Ih.
        if int(np.count_nonzero(labels_ira == 0)) != n:
            raise SystemExit(f"IRA labelled {int(np.count_nonzero(labels_ira == 0))} of {n} hexagonal")

        def ira_call():
            pin(threads)
            classify_serial(ira, ih, ic, pos, np.asarray(lengths))

        left, right = abab(dseams_sym, ira_call)
        ver = ira.get_version()
        rows.append(row_of(
            "IRA", "shape", n, threads, cores, left, right,
            f"hexagonal {n} equal", f"{version}/{ver[0]}", model,
        ))
        note_finished("ira", n)
    rows.extend(sqq_rows(
        lib, held, pos, lengths, threads, cores, model, version, row_of, h, pin, abab,
    ))
    return rows
