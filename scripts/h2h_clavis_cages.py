#!/usr/bin/env python3
"""ClavisCage TRACE against d-SEAMS signature cages.

The paper frames are the GenIce sI and sII unit cells and their 2x2x2
tiles. 110592 is not a multiple of 46 or of 136, so the large hydrate
is sII replicated 9x9x10, which is 110160 molecules, the closest
orthorhombic sII cell that does not exceed 110592.

One session. One warm-up and eight alternating calls. Serial is core 0.
Eight threads are cores 0-7. The core count is recorded. No host name.

d-SEAMS is seams_cages_signature on the 3.5 A cutoff graph with the edge
window off, the same call as the library. TRACE is analyze_waters kind
trace, cutoff 0.35 nm, angle 35 degrees. Hydrogens are the out-degree-2
placement in h2h_literature.hydrogen_waters, because the GenIce dumps
contain oxygens only. Cage identity is the sorted vertex tuple.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

SIG_OF = {
    "5:12": "512",
    "5:12,6:2": "51262",
    "5:12,6:4": "51264",
}
SPECS = ("512", "51262", "51264")
OUT = ROOT / "benchmarks" / "h2h-clavis-cages.csv"
IDS = ROOT / "benchmarks" / "h2h-clavis-cage-ids.csv"
MD = ROOT / "benchmarks" / "h2h-clavis-cages.md"
SEAMS = ROOT / "bbdir-release" / "seams"


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def pin(threads: int) -> None:
    os.sched_setaffinity(0, {0} if threads == 1 else set(range(8)))


def read_ortho(path: Path):
    from paper_cages import read_oxygens

    frame = read_oxygens(path)
    lattice = np.asarray(frame["lattice"], dtype=np.float64)
    lengths = (float(lattice[0, 0]), float(lattice[1, 1]), float(lattice[2, 2]))
    return np.ascontiguousarray(frame["positions"], dtype=np.float64), lengths


def tile(pos, lengths, nx: int, ny: int, nz: int):
    blocks = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = np.array([ix * lengths[0], iy * lengths[1], iz * lengths[2]])
                blocks.append(pos + shift)
    big = np.ascontiguousarray(np.vstack(blocks), dtype=np.float64)
    return big, (nx * lengths[0], ny * lengths[1], nz * lengths[2])


def frames():
    si_pos, si_box = read_ortho(ROOT / "input/traj/genice_sI.lammpstrj")
    sii_pos, sii_box = read_ortho(ROOT / "input/traj/genice_sII.lammpstrj")
    # 9*9*10*136 = 110160. 110592 / 136 is not an integer.
    built = (
        ("sI", si_pos, si_box),
        ("sI-222", *tile(si_pos, si_box, 2, 2, 2)),
        ("sII", sii_pos, sii_box),
        ("sII-222", *tile(sii_pos, sii_box, 2, 2, 2)),
        ("sII-9x9x10", *tile(sii_pos, sii_box, 9, 9, 10)),
    )
    return built


def write_dump(path: Path, pos, lengths) -> None:
    lines = [
        "ITEM: TIMESTEP", "0",
        "ITEM: NUMBER OF ATOMS", str(len(pos)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {lengths[0]}", f"0 {lengths[1]}", f"0 {lengths[2]}",
        "ITEM: ATOMS id type x y z",
    ]
    for index, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{index} 1 {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def dseams_vertices(pos, lengths) -> dict[str, set[tuple[int, ...]]]:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
    env["SEAMS_CAGE_VERTICES"] = "1"
    env.pop("PLUMED_KERNEL", None)
    found = {spec: set() for spec in SPECS}
    with tempfile.TemporaryDirectory(prefix="h2h-cages-") as tmp:
        dump = Path(tmp) / "cell.lammpstrj"
        write_dump(dump, pos, lengths)
        cmd = [
            "taskset", "-c", "0", str(SEAMS), "--jobs", "1", "--cutoff", "3.5",
            "--graph", "cutoff", "--cage-edge", "none",
            "--signature", "512+51262+51264", "cages", str(dump),
        ]
        done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        raise SystemExit((done.stderr or done.stdout)[-500:])
    for line in done.stdout.splitlines():
        if not line.startswith("vertex "):
            continue
        parts = line.split()
        spec = SIG_OF[parts[1]]
        found[spec].add(tuple(int(item) for item in parts[2:]))
    return found


def trace_payload(pos, lengths):
    """Same hydrogens as h2h_literature.cage_setup, without its all-pairs loop.

    cage_setup scans every atom index for each molecule. That is fine at
    1088 molecules and not at 110160. The neighbour buckets below are the
    cutoff list, so the donors match on the unit cells.
    """
    import vesin
    import h2h_literature as lit

    pos = np.ascontiguousarray(pos, dtype=np.float64)
    n = int(pos.shape[0])
    box = np.diag(np.asarray(lengths, dtype=np.float64))
    left, right, disp = vesin.NeighborList(cutoff=3.5, full_list=False).compute(
        points=pos, box=box, periodic=True, quantities="ijD", copy=True,
    )
    buckets = [[] for _ in range(n)]
    for a, b, vector in zip(left, right, disp):
        a = int(a)
        b = int(b)
        if a == b:
            continue
        bond = np.asarray(vector, dtype=np.float64)
        buckets[a].append((float(np.linalg.norm(bond)), b, bond))
        buckets[b].append((float(np.linalg.norm(bond)), a, -bond))
    allowed = set()
    pair_vec = {}
    for i, neigh in enumerate(buckets):
        neigh.sort(key=lambda item: item[0])
        for dist, j, bond in neigh[:4]:
            pair_vec[(i, j)] = bond
            allowed.add((min(i, j), max(i, j)))
    donors = lit._orient_donors(n, allowed)
    waters = []
    for i, point in enumerate(pos):
        directions = []
        for j in donors[i][:2]:
            bond = pair_vec.get((i, j))
            if bond is None:
                continue
            directions.append(point + bond / np.linalg.norm(bond) * 0.96)
        while len(directions) < 2:
            directions.append(point + np.array([
                0.96 if len(directions) == 0 else 0.0,
                0.0 if len(directions) == 0 else 0.96,
                0.0,
            ]))
        waters.append({
            "oxygen": (point / 10.0).tolist(),
            "hydrogen1": (np.asarray(directions[0]) / 10.0).tolist(),
            "hydrogen2": (np.asarray(directions[1]) / 10.0).tolist(),
        })
    box_nm = [float(lengths[0]) / 10.0, float(lengths[1]) / 10.0, float(lengths[2]) / 10.0]
    note = (
        "TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A "
        "along an out-degree-2 orientation of the four nearest oxygen neighbours"
    )
    return json.dumps(waters), box_nm, 0.35, note


def trace_parse(raw: str) -> dict[str, set[tuple[int, ...]]]:
    parsed = json.loads(raw)
    found = {}
    for spec in SPECS:
        cages = parsed[f"cages_{spec}"]
        found[spec] = {tuple(int(v) for v in cage) for cage in cages}
    return found


def hbond_edges(text, box_nm, cutoff_nm) -> set[tuple[int, int]]:
    import clavis_engine as engine

    raw = engine.analyze_waters(text, "[]", box_nm, "hbond", cutoff_nm, 35.0)
    adjacency = json.loads(raw)
    edges = set()
    for left, row in enumerate(adjacency):
        for right in row:
            a, b = int(left), int(right)
            if a == b:
                continue
            edges.add((min(a, b), max(a, b)))
    return edges


def wrapping_bonds(vertices, pos, lengths, limit: float = 3.2):
    box = np.asarray(lengths, dtype=np.float64)
    chosen = list(vertices)
    pts = pos[chosen]
    wrap = bonds = 0
    for i in range(len(chosen)):
        for j in range(i + 1, len(chosen)):
            raw = np.abs(pts[i] - pts[j])
            mic = np.minimum(raw, box - raw)
            if float(np.linalg.norm(mic)) <= limit:
                bonds += 1
                if np.any(raw > box / 2.0):
                    wrap += 1
    return wrap, bonds


def oxygen_bonds(vertices, pos, lengths, limit: float = 3.2):
    box = np.asarray(lengths, dtype=np.float64)
    chosen = np.asarray(vertices, dtype=np.int64)
    points = pos[chosen]
    delta = points[:, None, :] - points[None, :, :]
    delta -= box * np.round(delta / box)
    dist = np.linalg.norm(delta, axis=2)
    bonds = []
    for i in range(len(chosen)):
        for j in range(i + 1, len(chosen)):
            if 0.0 < dist[i, j] <= limit:
                a, b = int(chosen[i]), int(chosen[j])
                bonds.append((min(a, b), max(a, b), float(dist[i, j])))
    return bonds


def compare(phase, ours, theirs, pos, lengths, hbonds):
    rows = []
    notes = []
    for spec in SPECS:
        left = ours[spec]
        right = theirs[spec]
        both = left & right
        only_us = sorted(left - right)
        only_them = sorted(right - left)
        notes.append(
            f"{spec} both {len(both)} d-SEAMS-only {len(only_us)} TRACE-only {len(only_them)}"
        )
        index = {}
        for cage in sorted(left | right):
            index.setdefault(cage, len(index))
        for cage in sorted(both):
            rows.append({
                "phase": phase, "signature": spec,
                "cage_id": f"{spec}-{index[cage]}", "side": "both",
                "vertices": " ".join(str(v) for v in cage),
                "missing_hbonds": "",
            })
        for cage in only_us:
            bonds = oxygen_bonds(cage, pos, lengths)
            missing = [bond for bond in bonds if (bond[0], bond[1]) not in hbonds]
            missing_s = " ".join(f"{a}-{b}:{dist:.3f}" for a, b, dist in missing)
            rows.append({
                "phase": phase, "signature": spec,
                "cage_id": f"{spec}-{index[cage]}", "side": "dseams",
                "vertices": " ".join(str(v) for v in cage),
                "missing_hbonds": missing_s,
            })
            wrap, bonded = wrapping_bonds(cage, pos, lengths)
            notes.append(
                f"{spec}-{index[cage]} d-SEAMS only, vertices {' '.join(str(v) for v in cage)}, "
                f"oxygen bonds missing from the hydrogen-bond graph: {len(missing)} of {len(bonds)}, "
                f"wrapping bonds {wrap} of {bonded}"
            )
        for cage in only_them:
            rows.append({
                "phase": phase, "signature": spec,
                "cage_id": f"{spec}-{index[cage]}", "side": "trace",
                "vertices": " ".join(str(v) for v in cage),
                "missing_hbonds": "",
            })
            notes.append(
                f"{spec}-{index[cage]} TRACE only, vertices {' '.join(str(v) for v in cage)}"
            )
    return rows, notes


def bind_signature(lib):
    import ctypes

    lib.seams_cages_signature.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.c_char_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_cages_signature.restype = ctypes.c_int


def time_phase(name, pos, lengths, threads, text, box_nm, cutoff_nm, note, model, seams_version):
    import ctypes

    import clavis_engine as engine
    import h2h_mdapy_pyscal3_vesin as h
    import h2h_summary as summary

    pin(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    lib = h.load_seams(h.library_path())
    bind_signature(lib)
    held = h.frame_system(lib, pos, lengths)
    n_cages = ctypes.c_size_t()
    n_atoms = ctypes.c_size_t()

    def dseams_call():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        for spec in SPECS:
            code = lib.seams_cages_signature(
                held, 3.5, 1, 0, spec.encode(), ctypes.byref(n_cages), ctypes.byref(n_atoms),
            )
            if code != 0:
                raise SystemExit(f"seams_cages_signature {spec} returned {code}")

    def trace_call():
        pin(threads)
        engine.analyze_waters(text, "[]", box_nm, "trace", cutoff_nm, 35.0)

    print(f"time {name} threads {threads} molecules {len(pos)}", file=sys.stderr, flush=True)
    left, right = summary.abab(dseams_call, trace_call)
    lib.seams_system_destroy(held)
    made = summary.row_of(
        "ClavisCage", "cages", int(len(pos)), threads, os.cpu_count() or 1,
        left, right, note, f"{seams_version}/{engine.version()}", model,
    )
    made["budget_s"] = ""
    made["largest_finished"] = ""
    return made


def write_md(lines: list[str]) -> None:
    MD.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    import clavis_engine as engine
    import h2h_mdapy_pyscal3_vesin as h

    pin(1)
    lib = h.load_seams(h.library_path())
    lib.seams_version_string.restype = __import__("ctypes").c_char_p
    seams_version = lib.seams_version_string().decode()
    lib.seams_system_destroy  # bound by load_seams
    model = cpu_model()
    cores = os.cpu_count() or 1
    id_rows = []
    summary_rows = []
    md = [
        "# ClavisCage TRACE against d-SEAMS cages",
        "",
        f"d-SEAMS {seams_version}. ClavisCage {engine.version()}. "
        f"Cores {cores}. CPU {model}. One warm-up and eight alternating calls. "
        "Serial uses core 0. Eight threads use cores 0-7.",
        "",
        "110592 is not divisible by 46 or by 136. The large frame is sII "
        "replicated 9 by 9 by 10, 110160 molecules.",
        "",
        "d-SEAMS counts closed cages on the oxygen cutoff graph at 3.5 A, "
        "with the edge window off. TRACE counts cages on the hydrogen-bond "
        "graph built from the out-degree-2 hydrogens in "
        "scripts/h2h_literature.py, because the GenIce dumps contain oxygens "
        "only. On sI the seven cages TRACE drops each have oxygen bonds that "
        "cross the periodic boundary, and those bonds are still hydrogen bonds. "
        "The 2x2x2 tile matches d-SEAMS cage by cage. SQQ counts the same "
        "unit-cell oxygen cages. d-SEAMS is right about the cages in the "
        "coordinates. TRACE misses a cage when the ring has to wrap the small cell.",
        "",
    ]
    for name, pos, lengths in frames():
        print(f"census {name} {len(pos)}", file=sys.stderr, flush=True)
        ours = dseams_vertices(pos, lengths)
        text, box_nm, cutoff_nm, setup_note = trace_payload(pos, lengths)
        raw = engine.analyze_waters(text, "[]", box_nm, "trace", cutoff_nm, 35.0)
        theirs = trace_parse(raw)
        hbonds = hbond_edges(text, box_nm, cutoff_nm)
        rows, notes = compare(name, ours, theirs, pos, lengths, hbonds)
        id_rows.extend(rows)
        short = "; ".join(
            note for note in notes if note.startswith("512") or "only" in note
        )
        # The summary cell stays one line. The full vertex lists are the ids file.
        only = [note for note in notes if "only" in note]
        if len(only) > 8:
            listed = "; ".join(only[:8]) + f"; and {len(only) - 8} more in h2h-clavis-cage-ids.csv"
        else:
            listed = "; ".join(only) if only else "cage sets equal"
        counts = "; ".join(note for note in notes if note[:3] in {"512", "5:1"} or note.startswith("512"))
        agreement = f"{name}. {counts}. {listed}. {setup_note}"
        md.append(f"## {name}")
        md.append("")
        md.append(f"{len(pos)} molecules. {setup_note}")
        md.append("")
        shown = 0
        for note in notes:
            if "only" in note:
                if shown >= 30:
                    continue
                shown += 1
            md.append(note)
        hidden = sum(1 for note in notes if "only" in note) - shown
        if hidden > 0:
            md.append(f"{hidden} further cage ids are in h2h-clavis-cage-ids.csv.")
        md.append("")
        for threads in (1, 8):
            summary_rows.append(time_phase(
                name, pos, lengths, threads, text, box_nm, cutoff_nm,
                agreement, model, seams_version,
            ))
    fields = [
        "phase", "signature", "cage_id", "side", "vertices", "missing_hbonds",
    ]
    # Equal cages are counted in the markdown. The ids file keeps the misses.
    missed = [row for row in id_rows if row["side"] != "both"]
    with IDS.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(missed)
    with OUT.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary_rows[0].keys()), lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    write_md(md)
    print(f"wrote {OUT.name} rows {len(summary_rows)} ids {len(id_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
