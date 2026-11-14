#!/usr/bin/env python3
"""One sample of a competitor, for the 600 s summary budget.

Exit 0 only when that sample returns. The caller kills the process
when the wall clock exceeds the budget.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

CUTOFF = 3.5


def cell(molecules: int):
    from h2h_mdapy_pyscal3_vesin import ice_positions

    for count in (2, 4, 8, 24, 48):
        if count ** 3 * 8 == molecules:
            return ice_positions(count)
    raise SystemExit(f"no ice cell has {molecules} molecules")


def cycless(molecules: int) -> None:
    import networkx as nx
    import numpy as np
    from cycless.cycles import cycles_iter
    import freud

    pos, lengths = cell(molecules)
    pos = np.ascontiguousarray(pos)
    box = freud.box.Box(*lengths)
    query = freud.locality.AABBQuery(box, pos).query(
        pos, {"r_max": CUTOFF, "exclude_ii": True}
    ).toNeighborList()
    graph = nx.Graph()
    graph.add_nodes_from(range(molecules))
    ii = np.asarray(query[:, 0])
    jj = np.asarray(query[:, 1])
    graph.add_edges_from((int(a), int(b)) for a, b in zip(ii, jj) if int(a) < int(b))
    found = sum(1 for _ in cycles_iter(graph, 6))
    if found < 1:
        raise SystemExit("cycless returned no rings")
    print(f"cycless {molecules} rings {found}", file=sys.stderr, flush=True)


def ira(molecules: int) -> None:
    import numpy as np
    import ira_mod
    from gen_ice_ic import ice_ic
    from h2h_ira import classify_serial, cluster

    pos, lengths = cell(molecules)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.asarray(lengths, dtype=np.float64)
    engine = ira_mod.IRA()
    ih = cluster(pos, box, 0)
    ic_pos, side = ice_ic(4)
    ic = cluster(np.ascontiguousarray(ic_pos), np.asarray((side, side, side)), 0)
    labels = classify_serial(engine, ih, ic, pos, box)
    if int(labels.shape[0]) != molecules:
        raise SystemExit("IRA short")
    print(f"ira {molecules} hexagonal {int((labels == 0).sum())}", file=sys.stderr, flush=True)


def sqq(molecules: int, threads: int) -> None:
    import numpy as np
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame

    pos, lengths = cell(molecules)
    pos_nm = np.asarray(pos, dtype=np.float64) / 10.0
    box = np.asarray(lengths, dtype=np.float64) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(molecules)
    ]
    frame = Frame(name=f"ice-{molecules}", atoms=atoms, box=box)
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
    print(f"sqq {molecules} cages {len(result.cages)}", file=sys.stderr, flush=True)


def _ovito_frame(molecules: int):
    import numpy as np
    from ovito.data import DataCollection

    pos, lengths = cell(molecules)
    data = DataCollection()
    data.create_cell(np.diag(lengths), pbc=(True, True, True))
    parts = data.create_particles(count=molecules)
    parts.create_property("Position", data=np.ascontiguousarray(pos))
    return data


def ovito_rings(molecules: int, threads: int) -> None:
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().setMaxThreadCount(threads)
    data = _ovito_frame(molecules)
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=6))
    table = pipe.compute().tables["ring-size-histogram"]
    hexagons = 0
    for size, count in zip(table["Ring Size"][...], table["Count"][...]):
        if int(size) == 6:
            hexagons = int(count)
    if hexagons < 1:
        raise SystemExit("OVITO rings returned no hexagons")
    print(f"ovito-rings {molecules} hexagons {hexagons}", file=sys.stderr, flush=True)


def ovito_ptm(molecules: int, threads: int) -> None:
    import numpy as np
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().setMaxThreadCount(threads)
    data = _ovito_frame(molecules)
    mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
    mod.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
    mod.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(mod)
    labels = np.asarray(pipe.compute().particles["Structure Type"])
    if int(labels.shape[0]) != molecules:
        raise SystemExit("OVITO PTM short")
    print(f"ovito-ptm {molecules} hexagonal {int((labels == 7).sum())}", file=sys.stderr, flush=True)


def plumed(molecules: int) -> None:
    import numpy as np
    import plumed

    kernel = (
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
        "sapphirerapids/software/PLUMED/2.10.0-foss-2026.1/lib/libplumedKernel.so"
    )
    os.environ.setdefault("PLUMED_KERNEL", kernel)
    pos, lengths = cell(molecules)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    line = f"obs: Q6 SPECIES=1-{molecules} SWITCH={{CUSTOM FUNC=step(3.5-x) R_0=1 D_MAX=4.2}}"
    calc = plumed.Plumed()
    calc.cmd("setMDEngine", "python")
    calc.cmd("setNatoms", molecules)
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
    calc.cmd("setStep", 1)
    calc.cmd("setBox", box)
    calc.cmd("setPositions", pos)
    calc.cmd("setMasses", np.ones(molecules, dtype=np.float64))
    calc.cmd("setCharges", np.zeros(molecules, dtype=np.float64))
    calc.cmd("setForces", np.zeros((molecules, 3), dtype=np.float64))
    calc.cmd("setVirial", np.zeros((3, 3), dtype=np.float64))
    calc.cmd("calc")
    if data.size < molecules:
        raise SystemExit("PLUMED q6 short")
    print(f"plumed {molecules} q6 {float(np.max(np.abs(data))):.6e}", file=sys.stderr, flush=True)


def ptm_lib(molecules: int, threads: int) -> None:
    import ctypes
    import numpy as np

    pos, lengths = cell(molecules)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    lib = ctypes.CDLL(str(ROOT / "bbdir-release" / "libptm_ice.so"))
    lib.ptm_ice_classify.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.ptm_ice_classify.restype = ctypes.c_int
    box = (ctypes.c_double * 3)(*lengths)
    types = (ctypes.c_int * molecules)()
    rmsd = (ctypes.c_double * molecules)()
    code = lib.ptm_ice_classify(
        pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        molecules, box, 0.1, threads, types, rmsd,
    )
    if code != 0:
        raise SystemExit(f"ptm_ice_classify {code}")
    print(f"ptm-lib {molecules} hexagonal {sum(1 for v in types if v == 7)}", file=sys.stderr, flush=True)


def main() -> None:
    kind = sys.argv[1]
    molecules = int(sys.argv[2])
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    os.sched_setaffinity(0, set(range(8)) if threads > 1 else {0})
    if kind == "cycless":
        cycless(molecules)
    elif kind == "ira":
        ira(molecules)
    elif kind == "sqq":
        sqq(molecules, threads)
    elif kind == "ovito-rings":
        ovito_rings(molecules, threads)
    elif kind == "ovito-ptm":
        ovito_ptm(molecules, threads)
    elif kind == "plumed":
        plumed(molecules)
    elif kind == "ptm-lib":
        ptm_lib(molecules, threads)
    else:
        raise SystemExit(f"unknown probe {kind}")


if __name__ == "__main__":
    main()
