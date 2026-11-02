#!/usr/bin/env python3
"""Fresh Table 2 rows for OVITO, freud and GenIce2, plus the d-SEAMS side.

Five warmup repetitions and 30 recorded repetitions, pinned to core 3.
Samples are appended to timings.json. Medians are left to render_timings.py.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import types
from collections import Counter
from pathlib import Path

import numpy as np

from timing_protocol import CORE, CUT, REPS, WARMUP

# Documented OVITO default, not a per-system fit.
PTM_RMSD = 0.1

_ROOT = Path(__file__).resolve().parents[1]
_H2H = Path(os.environ.get("H2H_ROOT", str(_ROOT / "h2h")))
PY = os.environ.get("PYH2H", str(_H2H / "sotaenv" / "bin" / "python"))
BENCH = os.environ.get("BENCH_PAPER", str(_ROOT / "build" / "tests" / "bench_paper"))
PROTOCOL = os.environ.get("BENCH_PROTOCOL", str(_ROOT / "build" / "tests" / "bench_protocol"))
SEAMS = os.environ.get("SEAMS", str(_ROOT / "build" / "seams"))
RINGS = os.environ.get("RINGS", str(_H2H / "rings" / "bin" / "rings"))
JSON_PATH = Path(os.environ.get("TIMING_JSON", str(_H2H / "out" / "timings.json")))
DATA = Path(os.environ.get("H2H_DATA", str(_H2H / "data")))

DS_CLASS = {
    "cubic": "cubic",
    "reCubic": "cubic",
    "hexagonal": "hexagonal",
    "reHex": "hexagonal",
    "interfacial": "interfacial",
    "clathrate": "hydrate",
    "interClathrate": "interfacial-hydrate",
    "water": "other",
    "unclassified": "other",
    "unknown": "other",
}
OV_CLASS = {0: "other", 1: "hexagonal", 2: "cubic", 3: "interfacial", 4: "hydrate", 5: "interfacial-hydrate"}
DIAMOND_SHELL = {1: "cubic", 2: "cubic", 3: "cubic", 4: "hexagonal", 5: "hexagonal", 6: "hexagonal"}


def env():
    """C++ binaries are linked against the EESSI libc. Do not prefer /lib."""
    e = os.environ.copy()
    e["OMP_NUM_THREADS"] = "1"
    e["OMP_PROC_BIND"] = "true"
    parts = [p for p in e.get("LD_LIBRARY_PATH", "").split(":") if p and not p.startswith("/lib/") and not p.startswith("/usr/lib/")]
    if parts:
        e["LD_LIBRARY_PATH"] = ":".join(parts)
    else:
        e.pop("LD_LIBRARY_PATH", None)
    return e


def samples_of(fn):
    last = None
    for _ in range(WARMUP):
        last = fn()
    out = []
    for _ in range(REPS):
        t0 = time.perf_counter()
        last = fn()
        out.append([(time.perf_counter() - t0) * 1e3])
    return out, last


def session(sid, comparison, frame, side, code, command, samples, counts, agreement, scope="kernel"):
    return {
        "id": sid,
        "role": "head-to-head",
        "comparison": comparison,
        "scope": scope,
        "frame": frame,
        "side": side,
        "code_name": code,
        "command": command,
        "cwd": "",
        "warmup": WARMUP,
        "columns": ["kernel_ms"] if scope == "kernel" else ["wall_ms"],
        "samples": samples,
        "counts": counts,
        "agreement": agreement,
        "load_start": Path("/proc/loadavg").read_text().strip(),
        "load_end": Path("/proc/loadavg").read_text().strip(),
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def run_bench(mode, traj, sidecar):
    cmd = ["taskset", "-c", str(CORE), BENCH, str(traj), "1", str(WARMUP), str(REPS), mode, str(sidecar)]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env())
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-500:])
    samples = []
    for line in proc.stdout.splitlines():
        try:
            samples.append([float(line)])
        except ValueError:
            continue
    if len(samples) != REPS:
        raise RuntimeError(f"{mode} produced {len(samples)} samples\n{proc.stdout[-400:]}")
    return samples, cmd


def read_lmp(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[3])
    box = []
    for line in lines[5:8]:
        a, b = map(float, line.split()[:2])
        box.append(b - a)
    pos = []
    for line in lines[9 : 9 + n]:
        p = line.split()
        pos.append((float(p[2]), float(p[3]), float(p[4])))
    return np.asarray(pos), np.asarray(box)


def ovito_data(pos, box):
    from ovito.data import DataCollection

    data = DataCollection()
    data.create_cell(matrix=[[box[0], 0, 0, 0], [0, box[1], 0, 0], [0, 0, box[2], 0]], pbc=(True, True, True))
    particles = data.create_particles(count=len(pos))
    particles.create_property("Position", data=np.asarray(pos, dtype=np.float64))
    return data


def chill_pair(frame_name, traj):
    side = DATA / f"{frame_name}-chill.txt"
    samples, cmd = run_bench("chill", traj, side)
    ds = [DS_CLASS[line.strip()] for line in side.read_text().splitlines() if line.strip()]
    pos, box = read_lmp(traj)

    def once():
        data = ovito_data(pos, box)
        from ovito.modifiers import ChillPlusModifier

        data.apply(ChillPlusModifier(cutoff=CUT))
        return np.asarray(data.particles["Structure Type"]).copy()

    ov_samples, ov = samples_of(once)
    ov_m = [OV_CLASS[int(v)] for v in ov]
    same = sum(a == b for a, b in zip(ds, ov_m))
    agree = f"{same} of {len(ds)}"
    return [
        session(f"chill-{frame_name}-dseams", f"CHILL+ cutoff {CUT} A", frame_name, "dseams", "d-SEAMS 2.0", cmd, samples, {"labels": dict(Counter(ds))}, agree),
        session(f"chill-{frame_name}-ovito", f"CHILL+ cutoff {CUT} A", frame_name, "ovito", "OVITO 3.15.5 ChillPlus", [f"ChillPlusModifier(cutoff={CUT})"], ov_samples, {"labels": dict(Counter(ov_m))}, agree),
    ]


def q6_pair(frame_name, traj):
    side = DATA / f"{frame_name}-q6.txt"
    samples, cmd = run_bench("q6", traj, side)
    text = side.read_text().splitlines()
    bonds, values = [], []
    mode = None
    for line in text:
        if line == "bonds":
            mode = "bonds"
            continue
        if line == "values":
            mode = "values"
            continue
        if mode == "bonds":
            a, b = map(int, line.split())
            bonds.append((a, b))
        elif mode == "values":
            q, qb = map(float, line.split())
            values.append((q, qb))
    q_ds = np.array([v[0] for v in values])
    qb_ds = np.array([v[1] for v in values])
    pos, box = read_lmp(traj)
    import freud

    freud.parallel.set_num_threads(1)
    fbox = freud.box.Box(Lx=float(box[0]), Ly=float(box[1]), Lz=float(box[2]))
    pts = np.asarray(pos, dtype=np.float32) - box / 2.0
    qi = np.array([a for a, _ in bonds], dtype=np.uint32)
    pj = np.array([b for _, b in bonds], dtype=np.uint32)
    vec = fbox.wrap(pts[pj] - pts[qi])
    nlist = freud.locality.NeighborList.from_arrays(len(pos), len(pos), qi, pj, vec)

    def once():
        op = freud.order.Steinhardt(l=6, average=False)
        op.compute((fbox, pts), neighbors=nlist)
        opb = freud.order.Steinhardt(l=6, average=True)
        opb.compute((fbox, pts), neighbors=nlist)
        return np.asarray(op.particle_order).copy(), np.asarray(opb.particle_order).copy()

    fr_samples, (q_fr, qb_fr) = samples_of(once)
    d6 = float(np.nanmax(np.abs(q_ds - q_fr)))
    d6b = float(np.nanmax(np.abs(qb_ds - qb_fr)))
    agree = f"max |q6| {d6:.3e}, |q6bar| {d6b:.3e}"
    return [
        session(f"q6-{frame_name}-dseams", "q6 shared mutual 4-NN", frame_name, "dseams", "d-SEAMS 2.0", cmd, samples, {}, agree),
        session(f"q6-{frame_name}-freud", "q6 shared mutual 4-NN", frame_name, "freud", "freud 3.5.0", ["Steinhardt(l=6) on the d-SEAMS mutual list"], fr_samples, {}, agree),
    ]


def cna_pair(frame_name, traj):
    side = DATA / f"{frame_name}-seeded.txt"
    samples, cmd = run_bench("seeded", traj, side)
    ds = [line.strip() for line in side.read_text().splitlines() if line.strip()]
    pos, box = read_lmp(traj)

    def once():
        data = ovito_data(pos, box)
        from ovito.modifiers import IdentifyDiamondModifier

        data.apply(IdentifyDiamondModifier())
        return np.asarray(data.particles["Structure Type"]).copy()

    ov_samples, dia = samples_of(once)
    a = np.array([x if x in ("cubic", "hexagonal") else "other" for x in ds])
    shell = np.array([DIAMOND_SHELL.get(int(t), "other") for t in dia])
    same = int((a == shell).sum())
    agree = f"{same} of {len(ds)}"
    return [
        session(f"cna-{frame_name}-dseams", "seeded cages vs identify diamond", frame_name, "dseams", "d-SEAMS 2.0", cmd, samples, {}, agree),
        session(f"cna-{frame_name}-ovito", "seeded cages vs identify diamond", frame_name, "ovito", "OVITO 3.15.5 IdentifyDiamond", ["IdentifyDiamondModifier"], ov_samples, {}, agree),
    ]


def ptm_pair(traj):
    side = DATA / "ic512-seeded-adj.txt"
    samples, cmd = run_bench("seeded-adj", traj, side)
    ds = [line.strip() for line in side.read_text().splitlines() if line.strip()]
    pos, box = read_lmp(traj)

    def once():
        from ovito.data import DataCollection
        from ovito.modifiers import PolyhedralTemplateMatchingModifier
        from ovito.pipeline import Pipeline, StaticSource

        data = DataCollection()
        data.create_cell(matrix=[[box[0], 0, 0, 0], [0, box[1], 0, 0], [0, 0, box[2], 0]], pbc=(True, True, True))
        particles = data.create_particles(count=len(pos))
        particles.create_property("Position", data=np.asarray(pos, dtype=np.float64))
        pipeline = Pipeline(source=StaticSource(data=data))
        # OVITO documents 0.1 as a reasonable RMSD cutoff:
        # https://www.ovito.org/docs/current/reference/pipelines/modifiers/polyhedral_template_matching.html
        mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=PTM_RMSD)
        mod.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
        mod.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
        pipeline.modifiers.append(mod)
        out = pipeline.compute()
        return np.asarray(out.particles["Structure Type"]).copy()

    ov_samples, types = samples_of(once)
    from ovito.modifiers import PolyhedralTemplateMatchingModifier as P

    n_ds = sum(x == "cubic" for x in ds)
    n_ptm = int((types == int(P.Type.CUBIC_DIAMOND)).sum())
    agree = f"{n_ds} and {n_ptm} Ic"
    return [
        session("ptm-ic-dseams", "completed seeded cages vs PTM", "jittered Ic, 512", "dseams", "d-SEAMS 2.0", cmd, samples, {}, agree),
        session("ptm-ic-ovito", "completed seeded cages vs PTM", "jittered Ic, 512", "ovito", "OVITO 3.15.5 PTM", [f"PolyhedralTemplateMatchingModifier(rmsd_cutoff={PTM_RMSD})"], ov_samples, {}, agree),
    ]


def genice_rings(frame_name, traj):
    # d-SEAMS kernel is bench_protocol rings mode.
    cmd = ["taskset", "-c", str(CORE), PROTOCOL, str(traj), "1", str(WARMUP), str(REPS), "rings"]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env())
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-400:])
    samples, census = [], ""
    for line in proc.stdout.splitlines():
        if line.startswith("ring_census "):
            census = line.split(" ", 1)[1]
        else:
            try:
                samples.append([float(line)])
            except ValueError:
                pass
    pos, box = read_lmp(traj)

    def once():
        import networkx as nx
        import pairlist as pl
        from cycless.cycles import cycles_iter

        rel = pos / box
        graph = nx.Graph()
        graph.add_nodes_from(range(len(pos)))
        graph.add_edges_from(pl.pairs_iter(rel, CUT, np.diag(box), distance=False))
        return [tuple(int(x) for x in ring) for ring in cycles_iter(graph, 6, pos=rel)]

    gi_samples, rings = samples_of(once)
    c_gi = Counter(len(r) for r in rings)
    gi_census = " ".join(f"{s}:{c_gi.get(s, 0)}" for s in range(3, 7))
    return [
        session(f"gring-{frame_name}-dseams", "shortcut-free rings D=6", frame_name, "dseams", "d-SEAMS 2.0", cmd, samples, {"ring_census": census}, ""),
        session(f"gring-{frame_name}-genice", "shortcut-free rings D=6", frame_name, "genice", "GenIce2 2.2.13.3 cycless 0.7", ["cycless.cycles.cycles_iter"], gi_samples, {"ring_census": gi_census}, ""),
    ]


def genice_cages(traj, frame_name):
    cmd = ["taskset", "-c", str(CORE), PROTOCOL, str(traj), "1", str(WARMUP), str(REPS), "hydrate"]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env())
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-400:])
    samples = []
    sigs = {}
    for line in proc.stdout.splitlines():
        if line.startswith("sig "):
            body = line[4:]
            name, count = body.rsplit(" ", 1)
            sigs[name] = int(count)
        else:
            try:
                samples.append([float(line)])
            except ValueError:
                pass
    pos, box = read_lmp(traj)

    def once():
        import networkx as nx
        import pairlist as pl
        from cycless.cycles import cycles_iter
        from cycless.polyhed import polyhedra_iter

        rel = pos / box
        graph = nx.Graph()
        graph.add_nodes_from(range(len(pos)))
        graph.add_edges_from(pl.pairs_iter(rel, CUT, np.diag(box), distance=False))
        ringlist = [tuple(int(x) for x in ring) for ring in cycles_iter(graph, 6, pos=rel)]
        counts = Counter()
        for cage in polyhedra_iter(ringlist, 14):
            if len(cage) in (12, 14) and all(len(ringlist[i]) in (5, 6) for i in cage):
                faces = Counter(len(ringlist[i]) for i in cage)
                if faces.get(5) == 12 and faces.get(6, 0) == 0:
                    counts["5:12"] += 1
                elif faces.get(5) == 12 and faces.get(6) == 2:
                    counts["5:12,6:2"] += 1
        return counts

    gi_samples, counts = samples_of(once)
    cages = {"5:12": sigs.get("5:12", 0), "5:12,6:2": sigs.get("5:12,6:2", 0), "5:12,6:4": sigs.get("5:12,6:4", 0)}
    other = {"5:12": int(counts.get("5:12", 0)), "5:12,6:2": int(counts.get("5:12,6:2", 0)), "5:12,6:4": 0}
    return [
        session(f"gcage-{frame_name}-dseams", "5^12 and 5^12 6^2 cages", frame_name, "dseams", "d-SEAMS 2.0", cmd, samples, {"cages": cages}, ""),
        session(f"gcage-{frame_name}-genice", "5^12 and 5^12 6^2 cages", frame_name, "genice", "GenIce2 cycless polyhedra", ["cycless.polyhed.polyhedra_iter"], gi_samples, {"cages": other}, ""),
    ]


def ensure_frames():
    if (DATA / "niu_final.lammpstrj").exists() and (DATA / "crn1.lammpstrj").exists() and (DATA / "ic512.lammpstrj").exists():
        return
    raise SystemExit("frames missing; generate them before paper_rows.py")


def merge(new_sessions):
    payload = json.loads(JSON_PATH.read_text())
    by_id = {s["id"]: s for s in payload["sessions"]}
    for session_ in new_sessions:
        by_id[session_["id"]] = session_
    # Preserve order: old sessions, then any brand new ids.
    order = [s["id"] for s in payload["sessions"]]
    for session_ in new_sessions:
        if session_["id"] not in order:
            order.append(session_["id"])
    payload["sessions"] = [by_id[i] for i in order]
    payload["paper_rows_updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print("merged", len(new_sessions), "sessions into", JSON_PATH)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "short"
    os.sched_setaffinity(0, {CORE})
    ensure_frames()
    made = []
    mw = str(_ROOT / "input" / "traj" / "mW_cubic.lammpstrj")
    niu = DATA / "niu_final.lammpstrj"
    crn = DATA / "crn1.lammpstrj"
    ic = DATA / "ic512.lammpstrj"
    si = Path(os.environ.get("H2H_WORK", str(_H2H / "work"))) / "timing-frames" / "si2.lammpstrj"
    if which in ("short", "all", "dseams"):
        made += chill_pair("mW cubic, 4096", mw)
        made += chill_pair("Niu final, 2880", niu)
        made += q6_pair("mW cubic, 4096", mw)
        made += q6_pair("Niu final, 2880", niu)
        made += cna_pair("mW cubic, 4096", mw)
        made += cna_pair("Niu final, 2880", niu)
        made += cna_pair("jittered Ic, 512", ic)
        made += ptm_pair(ic)
        made += genice_rings("mW cubic, 4096", mw)
        made += genice_rings("CRN1, 1000", crn)
        made += genice_cages(si, "sI 2x2x2, 368")
    if which in ("rings", "all"):
        pass  # the long R.I.N.G.S. rows are a separate driver
    merge(made)
    print("done", which, "n", len(made))


if __name__ == "__main__":
    main()
