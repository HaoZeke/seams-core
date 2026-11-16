#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Per-atom CHILL+ labels of two libseams builds, frame by frame.

Each build runs in its own process at each OpenMP width. A case is one
frame, one rule (published or OVITO), one graph and one species filter;
the record is the SHA-256 of the int32 label array and the nine class
counts. The two builds agree when every record is the same. counts.csv
holds the class counts of every case; labels-*.jsonl hold the records.

Frames: every frame of every LAMMPS dump in input/traj; each 4096-atom
frame again tiled 3x3x3 (110592) and 6x6x6 (884736); ice Ic and ice Ih
at 4096, 110592 and 884736, perfect and with seeded Gaussian noise; a
random liquid-density frame at each size; ice Ic in a tilted cell; and a
two-species frame classified on one species.

With --a-lang and --b-lang the language front (seams_frame_chill_plus_rule)
runs the same cases, and one more frame whose atom ids repeat.

    python3 scripts/chill_plus_equivalence.py \\
        --a bbdir-before/src/libseams.so --b bbdir-release/src/libseams.so \\
        --a-lang bbdir-before/bindings/libseams_lang.so \\
        --b-lang bbdir-release/bindings/libseams_lang.so \\
        --a-tree 00129f1e --b-tree "$(git rev-parse HEAD^{tree})" \\
        --out benchmarks/chill-plus-1t/equivalence
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTOFF = 3.5
TILES = {4096: (1, 3, 6)}
SIZES = (4096, 110592, 884736)
NOISE = (0.0, 0.15, 0.35)
RULES = {"published": 0, "ovito": 1}
GRAPHS = {"cutoff": 0, "mutual": 1}
WIDTHS = (1, 8)


def read_dump(path: Path):
    """Frames of a LAMMPS dump: (positions, types, ids, box9)."""
    import numpy as np
    frames = []
    with path.open(encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    k = 0
    while k < len(lines):
        if not lines[k].startswith("ITEM: TIMESTEP"):
            k += 1
            continue
        n = int(lines[k + 3].split()[0])
        header = lines[k + 4].split()
        bounds = [[float(v) for v in lines[k + 5 + d].split()] for d in range(3)]
        if "xy" in header:
            xy, xz, yz = bounds[0][2], bounds[1][2], bounds[2][2]
            xlo = bounds[0][0] - min(0.0, xy, xz, xy + xz)
            xhi = bounds[0][1] - max(0.0, xy, xz, xy + xz)
            ylo = bounds[1][0] - min(0.0, yz)
            yhi = bounds[1][1] - max(0.0, yz)
        else:
            xy = xz = yz = 0.0
            xlo, xhi = bounds[0][:2]
            ylo, yhi = bounds[1][:2]
        zlo, zhi = bounds[2][:2]
        box9 = np.array([xlo, ylo, zlo, xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz])
        cols = lines[k + 8].split()[2:]
        rows = np.array([line.split()[:len(cols)] for line in lines[k + 9:k + 9 + n]],
                        dtype=np.float64)
        at = {name: c for c, name in enumerate(cols)}
        ids = rows[:, at["id"]].astype(np.int32)
        types = rows[:, at["type"]].astype(np.int32)
        if "x" in at:
            pos = rows[:, [at["x"], at["y"], at["z"]]]
        else:
            frac = rows[:, [at["xs"], at["ys"], at["zs"]]]
            h = np.array([[box9[3], 0, 0], [xy, box9[4], 0], [xz, yz, box9[5]]])
            pos = frac @ h + box9[:3]
        order = np.argsort(ids, kind="stable")
        frames.append((np.ascontiguousarray(pos[order]), types[order], ids[order], box9))
        k += 9 + n
    return frames


def tile(pos, types, box9, m: int):
    import numpy as np
    if m == 1:
        return pos, types, box9
    if any(box9[6:]):
        raise ValueError("tiling a tilted cell")
    shifts = np.array([(i, j, k) for i in range(m) for j in range(m) for k in range(m)],
                      dtype=np.float64) * box9[3:6]
    big = (pos[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    out = box9.copy()
    out[3:6] *= m
    return np.ascontiguousarray(big), np.tile(types, m ** 3), out


def ice(kind: str, molecules: int, sigma: float, seed: int):
    """Ice Ic (cubic cell) or ice Ih (orthohexagonal cell), 8 molecules a cell."""
    import numpy as np
    ncell = round((molecules / 8) ** (1.0 / 3.0))
    assert 8 * ncell ** 3 == molecules
    bond = 2.75
    if kind == "Ic":
        a = 4.0 * bond / np.sqrt(3.0)
        cell = np.array([a, a, a])
        base = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        frac = [b for b in base] + [(x + 0.25, y + 0.25, z + 0.25) for x, y, z in base]
    else:
        a = bond * np.sqrt(8.0 / 3.0)
        c = 8.0 * bond / 3.0
        cell = np.array([a, a * np.sqrt(3.0), c])
        u = 3.0 / 8.0
        frac = [(x, y, z + dz)
                for x, y, z in ((0.0, 1 / 3, 0.0), (0.5, 5 / 6, 0.0),
                                (0.5, 1 / 6, 0.5), (0.0, 2 / 3, 0.5))
                for dz in (0.0, u)]
    grid = np.array([(i, j, k) for i in range(ncell) for j in range(ncell) for k in range(ncell)],
                    dtype=np.float64)
    pos = ((grid[:, None, :] + np.array(frac)[None, :, :]) * cell).reshape(-1, 3)
    if sigma > 0:
        pos = pos + np.random.default_rng(seed).normal(0.0, sigma, pos.shape)
    box9 = np.array([0.0, 0.0, 0.0, *(cell * ncell), 0.0, 0.0, 0.0])
    return np.ascontiguousarray(pos), np.ones(molecules, dtype=np.int32), box9


def cases():
    """(name, positions, types, box9, species filters, ids or None)."""
    import numpy as np
    for path in sorted((ROOT / "input" / "traj").glob("*.lammpstrj")):
        for f, (pos, types, _ids, box9) in enumerate(read_dump(path), start=1):
            common = int(np.bincount(types).argmax())
            filters = (0,) if np.all(types == types[0]) else (0, common)
            for m in TILES.get(len(pos), (1,)):
                if m > 1 and any(box9[6:]):
                    continue
                p, t, b = tile(pos, types, box9, m)
                yield f"{path.stem}/f{f}/x{m}", p, t, b, filters, None
    for kind in ("Ic", "Ih"):
        for molecules in SIZES:
            for s, sigma in enumerate(NOISE):
                pos, types, box9 = ice(kind, molecules, sigma, 1000 * s + molecules)
                yield f"{kind}/{molecules}/sigma{sigma}", pos, types, box9, (0,), None
    for molecules in SIZES:
        side = (molecules / 0.0334) ** (1.0 / 3.0)
        pos = np.random.default_rng(molecules).uniform(0.0, side, (molecules, 3))
        box9 = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
        yield f"random/{molecules}", pos, np.ones(molecules, dtype=np.int32), box9, (0,), None
    pos, types, box9 = ice("Ic", 4096, 0.15, 7)
    box9 = box9.copy()
    box9[6:] = (1.1, -0.7, 0.4)
    yield "Ic/4096/tilted", pos, types, box9, (0,), None
    pos, types, box9 = ice("Ih", 110592, 0.25, 11)
    types = types.copy()
    types[::5] = 2
    yield "Ih/110592/two-species", pos, types, box9, (1, 2), None
    pos, types, box9 = ice("Ic", 4096, 0.15, 13)
    ids = np.arange(1, 4097, dtype=np.int32)
    ids[1::7] = ids[0::7][:len(ids[1::7])]
    yield "Ic/4096/repeated-ids", pos, types, box9, (0,), ids


LANG_RULES = {"published": b"chill+", "ovito": b"ovito-chill+"}


def worker(api_path: str, lang_path: str | None) -> None:
    import ctypes
    import numpy as np
    api = ctypes.CDLL(api_path)
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    vp = ctypes.c_void_p
    api.seams_system_create.restype = vp
    api.seams_system_destroy.argtypes = [vp]
    api.seams_system_set_frame.argtypes = [vp, dp, ctypes.c_int, dp, ip, ip]
    api.seams_chill_plus.argtypes = [vp, ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                     ctypes.c_int, ip, ip]
    lang = None
    if lang_path:
        lang = ctypes.CDLL(lang_path)
        lang.seams_frame_from_arrays.restype = vp
        lang.seams_frame_from_arrays.argtypes = [dp, ctypes.c_int, dp, ip, ip]
        lang.seams_frame_free.argtypes = [vp]
        lang.seams_frame_chill_plus_rule.argtypes = [
            vp, ctypes.c_double, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ip, ip]

    def emit(front, name, n, species, graph, rule, rc, labels, counts):
        sys.stdout.write(json.dumps({
            "front": front, "case": name, "n": n, "type": species, "graph": graph,
            "rule": rule, "rc": rc, "sha256": hashlib.sha256(labels.tobytes()).hexdigest(),
            "counts": counts.tolist(),
        }) + "\n")
        sys.stdout.flush()

    for name, pos, types, box9, filters, ids in cases():
        n = len(pos)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        box9 = np.ascontiguousarray(box9, dtype=np.float64)
        types = np.ascontiguousarray(types, dtype=np.int32)
        if ids is None:
            handle = api.seams_system_create()
            rc = api.seams_system_set_frame(handle, pos.ctypes.data_as(dp), n,
                                            box9.ctypes.data_as(dp), types.ctypes.data_as(ip),
                                            None)
            if rc != 0:
                raise SystemExit(f"{name}: seams_system_set_frame returned {rc}")
            for species in filters:
                for graph, g in GRAPHS.items():
                    for rule, r in RULES.items():
                        labels = np.zeros(n, dtype=np.int32)
                        counts = np.zeros(9, dtype=np.int32)
                        rc = api.seams_chill_plus(handle, CUTOFF, species, g, r,
                                                  labels.ctypes.data_as(ip),
                                                  counts.ctypes.data_as(ip))
                        emit("c", name, n, species, graph, rule, rc, labels, counts)
            api.seams_system_destroy(handle)
        if lang is None:
            continue
        idarr = None if ids is None else np.ascontiguousarray(ids, dtype=np.int32)
        frame = lang.seams_frame_from_arrays(
            pos.ctypes.data_as(dp), n, box9.ctypes.data_as(dp), types.ctypes.data_as(ip),
            None if idarr is None else idarr.ctypes.data_as(ip))
        if not frame:
            raise SystemExit(f"{name}: seams_frame_from_arrays failed")
        for species in filters:
            for graph in GRAPHS:
                for rule, text in LANG_RULES.items():
                    labels = np.zeros(n, dtype=np.int32)
                    counts = np.zeros(9, dtype=np.int32)
                    rc = lang.seams_frame_chill_plus_rule(
                        frame, CUTOFF, species, graph.encode(), text,
                        labels.ctypes.data_as(ip), counts.ctypes.data_as(ip))
                    emit("lang", name, n, species, graph, rule, rc, labels, counts)
        lang.seams_frame_free(frame)


def run(api: Path, lang: Path | None, threads: int, out: Path) -> list[dict]:
    env = dict(os.environ)
    width = str(threads)
    env.update({"OMP_NUM_THREADS": width, "RAYON_NUM_THREADS": width,
                "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"})
    cpus = "0" if threads == 1 else f"0-{threads - 1}"
    with out.open("w", encoding="utf-8") as handle:
        subprocess.run(["taskset", "-c", cpus, sys.executable, str(Path(__file__).resolve()),
                        "--worker", str(api.resolve()),
                        *([str(lang.resolve())] if lang else [])],
                       env=env, stdout=handle, check=True)
    return [json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]


def main() -> int:
    if len(sys.argv) > 2 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
        return 0
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--a", type=Path, required=True, help="libseams.so of the first build")
    parser.add_argument("--b", type=Path, required=True, help="libseams.so of the second build")
    parser.add_argument("--a-lang", type=Path, help="libseams_lang.so of the first build")
    parser.add_argument("--b-lang", type=Path, help="libseams_lang.so of the second build")
    parser.add_argument("--a-tree", default="", help="source tree of the first build")
    parser.add_argument("--b-tree", default="", help="source tree of the second build")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if (args.a_lang is None) != (args.b_lang is None):
        parser.error("give both --a-lang and --b-lang, or neither")
    args.out.mkdir(parents=True, exist_ok=True)
    records = {}
    for tag, api, lang in (("a", args.a, args.a_lang), ("b", args.b, args.b_lang)):
        for threads in WIDTHS:
            records[(tag, threads)] = run(api, lang, threads,
                                          args.out / f"labels-{tag}-t{threads}.jsonl")
    key = lambda row: (row["front"], row["case"], row["type"], row["graph"],  # noqa: E731
                       row["rule"])
    reference = {key(row): row for row in records[("a", 1)]}
    failures = []
    for (tag, threads), rows in records.items():
        got = {key(row): row for row in rows}
        if got.keys() != reference.keys():
            failures.append(f"{tag} t{threads}: different case list")
        for k, row in got.items():
            ref = reference.get(k)
            if ref is None:
                continue
            if row["rc"] != 0 or (row["sha256"], row["counts"]) != (ref["sha256"], ref["counts"]):
                failures.append(f"{tag} t{threads} {k}: {row['counts']} against {ref['counts']}")
    names = ("cubic", "hexagonal", "water", "interfacial", "clathrate", "interclathrate",
             "unclassified", "recubic", "rehex")
    with (args.out / "counts.csv").open("w", encoding="utf-8") as handle:
        handle.write(",".join(["front", "case", "n", "type", "graph", "rule", *names]) + "\n")
        for row in records[("a", 1)]:
            handle.write(",".join(str(v) for v in (
                row["front"], row["case"], row["n"], row["type"], row["graph"], row["rule"],
                *row["counts"])) + "\n")
    frames = len({row["case"] for row in records[("a", 1)]})
    summary = (f"{len(reference)} cases over {frames} frames, builds a and b, "
               f"widths {', '.join(map(str, WIDTHS))}: "
               + ("labels identical" if not failures else f"{len(failures)} differ"))
    (args.out / "summary.txt").write_text(
        "\n".join([f"a_source_tree={args.a_tree}", f"b_source_tree={args.b_tree}",
                   summary, *failures]) + "\n", encoding="utf-8")
    print(summary)
    for line in failures[:20]:
        print(line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
