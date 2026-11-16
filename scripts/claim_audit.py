#!/usr/bin/env python3
"""Re-derive the CHILL+ counts and the open-face padding claim.

The published rule and the OVITO rule are the library's. A pre-gate
rule, the one that scored four nearest bonds when a fifth oxygen sat
inside 3.5 A, is evaluated here from the same coordinates. It is not
a looser threshold: the staggered and eclipsed windows are unchanged,
and a shell that is not four is liquid under the rule the library
ships. OVITO 3.17.2 labels the same oxygen frames.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)
CODE = {name: i for i, name in enumerate(NAMES)}
CUTOFF = 3.5
CUT2 = CUTOFF * CUTOFF
CANCELLED = 1e-12
PI = math.pi
NORM = (
    0.25 * math.sqrt(7.0 / PI),
    0.125 * math.sqrt(21.0 / PI),
    0.25 * math.sqrt(105.0 / (2.0 * PI)),
    0.125 * math.sqrt(35.0 / PI),
)
# OVITO ChillPlusModifier structure codes.
OVITO_TO_SEAMS = {
    0: CODE["water"],
    1: CODE["hexagonal"],
    2: CODE["cubic"],
    3: CODE["interfacial"],
    4: CODE["clathrate"],
    5: CODE["interClathrate"],
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_dump(path: Path):
    frames = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i] != "ITEM: TIMESTEP":
            i += 1
            continue
        step = int(lines[i + 1])
        n = int(lines[i + 3])
        box_header = lines[i + 4]
        bounds = []
        tilt = [0.0, 0.0, 0.0]
        for axis in range(3):
            parts = lines[i + 5 + axis].split()
            lo = float(parts[0])
            hi = float(parts[1])
            bounds.append((lo, hi))
            if len(parts) > 2:
                tilt[axis] = float(parts[2])
        atom_at = i + 8
        header = lines[atom_at].split()
        cols = {name: k for k, name in enumerate(header[2:])}
        need = ("id", "type", "x", "y", "z")
        if any(name not in cols for name in need):
            raise SystemExit(f"{path} frame {step} has no id type x y z")
        ids = np.empty(n, dtype=np.int64)
        types = np.empty(n, dtype=np.int32)
        xyz = np.empty((n, 3), dtype=np.float64)
        for row in range(n):
            parts = lines[atom_at + 1 + row].split()
            ids[row] = int(parts[cols["id"]])
            types[row] = int(parts[cols["type"]])
            xyz[row, 0] = float(parts[cols["x"]])
            xyz[row, 1] = float(parts[cols["y"]])
            xyz[row, 2] = float(parts[cols["z"]])
        lengths = np.array([hi - lo for lo, hi in bounds], dtype=np.float64)
        origin = np.array([lo for lo, _hi in bounds], dtype=np.float64)
        frames.append(
            {
                "step": step,
                "ids": ids,
                "types": types,
                "xyz": xyz,
                "lengths": lengths,
                "origin": origin,
                "tilt": np.array(tilt, dtype=np.float64),
                "triclinic": "xy" in box_header,
            }
        )
        i = atom_at + 1 + n
    if not frames:
        raise SystemExit(f"no frames in {path}")
    return frames


def minimum_image(delta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    return delta - lengths * np.round(delta / lengths)


class Library:
    def __init__(self, path: Path):
        self.lib = ctypes.CDLL(str(path))
        self.lib.seams_frame_read.restype = ctypes.c_void_p
        self.lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
        self.lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        self.lib.seams_frame_nop.restype = ctypes.c_int
        self.lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
        self.lib.seams_frame_ids.restype = ctypes.POINTER(ctypes.c_int)
        self.lib.seams_frame_ids.argtypes = [ctypes.c_void_p]
        self.lib.seams_frame_chill_plus_rule.restype = ctypes.c_int
        self.lib.seams_frame_chill_plus_rule.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        self.lib.seams_last_error.restype = ctypes.c_char_p

    def labels(self, path: Path, frame: int, type_filter: int, rule: str) -> tuple[np.ndarray, np.ndarray]:
        held = self.lib.seams_frame_read(str(path).encode(), frame, type_filter)
        if not held:
            err = self.lib.seams_last_error()
            raise SystemExit(err.decode() if err else f"cannot read {path} frame {frame}")
        try:
            n = self.lib.seams_frame_nop(held)
            ids_p = self.lib.seams_frame_ids(held)
            ids = np.array([ids_p[i] for i in range(n)], dtype=np.int64)
            out = (ctypes.c_int * n)()
            counts = (ctypes.c_int * 9)()
            code = self.lib.seams_frame_chill_plus_rule(
                held, CUTOFF, 0, b"cutoff", rule.encode(), out, counts
            )
            if code != 0:
                err = self.lib.seams_last_error()
                raise SystemExit(err.decode() if err else f"{rule} failed")
            labels = np.array([out[i] for i in range(n)], dtype=np.int32)
            return ids, labels
        finally:
            self.lib.seams_frame_free(held)


def ylm3(dx: float, dy: float, dz: float) -> np.ndarray:
    acc = np.zeros(7, dtype=np.float64)
    r2 = dx * dx + dy * dy + dz * dz
    if not (r2 > 0.0):
        acc[0] = NORM[0] * 2.0
        return acc
    inv = 1.0 / (r2 * math.sqrt(r2))
    z2 = dz * dz
    t1 = 5.0 * z2 - r2
    acc[0] = NORM[0] * (dz * (5.0 * z2 - 3.0 * r2)) * inv
    acc[1] = NORM[1] * (dy * t1) * inv
    acc[2] = NORM[1] * (dx * t1) * inv
    acc[3] = NORM[2] * (dz * (dy * dy - dx * dx)) * inv
    acc[4] = NORM[2] * (2.0 * dx * dy * dz) * inv
    acc[5] = NORM[3] * (dy * (dy * dy - 3.0 * dx * dx)) * inv
    acc[6] = NORM[3] * (dx * (3.0 * dy * dy - dx * dx)) * inv
    return acc


def correlation(qi: np.ndarray, qj: np.ndarray) -> float:
    ni = qi[0] * qi[0] + 2.0 * float(np.dot(qi[1:], qi[1:]))
    nj = qj[0] * qj[0] + 2.0 * float(np.dot(qj[1:], qj[1:]))
    if not (ni > CANCELLED and nj > CANCELLED):
        return 0.0
    dot = qi[0] * qj[0] + 2.0 * float(np.dot(qi[1:], qj[1:]))
    return dot / math.sqrt(ni * nj)


def shells(xyz: np.ndarray, lengths: np.ndarray):
    n = xyz.shape[0]
    neigh = []
    for i in range(n):
        delta = minimum_image(xyz - xyz[i], lengths)
        d2 = np.einsum("ij,ij->i", delta, delta)
        d2[i] = np.inf
        which = np.flatnonzero(d2 < CUT2)
        order = sorted(which.tolist(), key=lambda j: (d2[j], j))
        rows = []
        for j in order:
            rows.append((j, float(delta[j, 0]), float(delta[j, 1]), float(delta[j, 2]), float(d2[j])))
        neigh.append(rows)
    return neigh


def bond_kind(value: float, closed: bool) -> str:
    if closed:
        if value <= -0.8:
            return "staggered"
        if -0.35 <= value <= 0.25:
            return "eclipsed"
    else:
        if value < -0.8:
            return "staggered"
        if -0.35 < value < 0.25:
            return "eclipsed"
    return "out"


def prepare(xyz: np.ndarray, lengths: np.ndarray):
    neigh = shells(xyz, lengths)
    q4 = []
    qall = []
    c4 = []
    call = []
    for i, rows in enumerate(neigh):
        four = rows[:4]
        acc4 = np.zeros(7)
        for _j, dx, dy, dz, _d2 in four:
            acc4 += ylm3(dx, dy, dz)
        acca = np.zeros(7)
        for _j, dx, dy, dz, _d2 in rows:
            acca += ylm3(dx, dy, dz)
        q4.append(acc4)
        qall.append(acca)
        c4.append([j for j, *_rest in four])
        call.append([j for j, *_rest in rows])
    return neigh, q4, qall, c4, call


def classify_one(i, neigh, q_used, chosen, closed: bool, neighbour_test: bool, gate: bool) -> tuple[int, list[float]]:
    shell = len(neigh[i])
    bonds = chosen[i]
    keep = len(bonds)
    values = [correlation(q_used[i], q_used[j]) for j in bonds]
    kinds = [bond_kind(value, closed) for value in values]
    nstag = kinds.count("staggered")
    neclip = kinds.count("eclipsed")
    four = shell == 4 and keep == 4
    if gate and not four:
        return CODE["water"], values
    if keep != 4:
        return CODE["water"], values

    def staggered_of(j: int) -> int:
        return sum(bond_kind(correlation(q_used[j], q_used[k]), closed) == "staggered" for k in chosen[j])

    inter = False
    if neighbour_test:
        if nstag == 2 and any(staggered_of(j) > 2 for j in bonds):
            inter = True
        if nstag == 3 and neclip == 0 and any(staggered_of(j) > 1 for j in bonds):
            inter = True
    elif nstag == 2 or (nstag == 3 and neclip == 0):
        inter = True
    if neclip == 0 and nstag == 4:
        return CODE["cubic"], values
    if neclip == 1 and nstag == 3:
        return CODE["hexagonal"], values
    if inter:
        return CODE["interfacial"], values
    if neclip == 4 and nstag == 0:
        return CODE["clathrate"], values
    if neclip == 3:
        return CODE["interClathrate"], values
    return CODE["water"], values


def classify_cloud(xyz, lengths, mode: str) -> tuple[np.ndarray, list[list[float]]]:
    neigh, q4, qall, c4, call = prepare(xyz, lengths)
    labels = np.empty(xyz.shape[0], dtype=np.int32)
    values = []
    if mode == "published":
        q_used, chosen, closed, test, gate = q4, c4, True, True, True
    elif mode == "ovito":
        q_used, chosen, closed, test, gate = qall, call, False, False, True
    elif mode == "pregate":
        # The rule before the shell-length gate. Four nearest bonds, closed
        # edges, neighbour test. A fifth neighbour inside the cutoff does
        # not by itself make the molecule liquid.
        q_used, chosen, closed, test, gate = q4, c4, True, True, False
    elif mode == "notest":
        q_used, chosen, closed, test, gate = q4, c4, True, False, True
    elif mode == "shellclosed":
        q_used, chosen, closed, test, gate = qall, call, True, False, True
    else:
        raise SystemExit(mode)
    for i in range(xyz.shape[0]):
        label, bond_values = classify_one(i, neigh, q_used, chosen, closed, test, gate)
        labels[i] = label
        values.append(bond_values)
    return labels, values


def cause_of(published: int, notest: int, shell_closed: int, ovito: int) -> str:
    if published == ovito:
        return "agree"
    neighbour = published != notest
    shell = notest != shell_closed
    edge = shell_closed != ovito
    flips = int(neighbour) + int(shell) + int(edge)
    if flips == 1 and neighbour:
        return "neighbour_test"
    if flips == 1 and shell:
        return "q3_shell"
    if flips == 1 and edge:
        return "open_edge"
    return "mixed"


def histogram(labels: np.ndarray) -> dict[str, int]:
    counts = {name: 0 for name in NAMES}
    for label in labels:
        counts[NAMES[int(label)]] += 1
    return counts


def require_match(tag: str, got: np.ndarray, ref: np.ndarray, ids: np.ndarray) -> None:
    if got.shape != ref.shape or np.any(got != ref):
        bad = np.flatnonzero(got != ref)
        sample = []
        for index in bad[:8]:
            sample.append(f"id {int(ids[index])} python {NAMES[int(got[index])]} library {NAMES[int(ref[index])]}")
        raise SystemExit(f"{tag} disagrees on {len(bad)} atoms: " + "; ".join(sample))


def oxygen_frame(frame: dict, type_filter: int) -> dict:
    if type_filter < 0:
        keep = np.ones(frame["types"].shape[0], dtype=bool)
    else:
        keep = frame["types"] == type_filter
    return {
        "ids": frame["ids"][keep],
        "xyz": frame["xyz"][keep],
        "lengths": frame["lengths"],
        "step": frame["step"],
    }


def ovito_labels(path: Path, type_filter: int) -> list[dict[int, int]]:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from ovito.io import import_file
    from ovito.modifiers import ChillPlusModifier, ExpressionSelectionModifier, DeleteSelectedModifier

    pipeline = import_file(str(path))
    if type_filter > 0:
        pipeline.modifiers.append(ExpressionSelectionModifier(expression=f"ParticleType != {type_filter}"))
        pipeline.modifiers.append(DeleteSelectedModifier())
    modifier = ChillPlusModifier()
    modifier.cutoff = CUTOFF
    pipeline.modifiers.append(modifier)
    out = []
    for index in range(pipeline.source.num_frames):
        data = pipeline.compute(index)
        ids = np.asarray(data.particles["Particle Identifier"])
        kind = np.asarray(data.particles["Structure Type"])
        mapped = {}
        for atom_id, label in zip(ids, kind):
            if int(label) not in OVITO_TO_SEAMS:
                raise SystemExit(f"OVITO structure code {int(label)} is not in the Chill+ table")
            mapped[int(atom_id)] = OVITO_TO_SEAMS[int(label)]
        out.append(mapped)
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def padding_frames(probe: Path, path: Path, lane: str, temperature: str) -> list[dict]:
    nframes = sum(1 for line in path.read_text().splitlines() if line == "ITEM: TIMESTEP")
    done = subprocess.run(
        [
            str(probe),
            "layers",
            str(path),
            lane,
            temperature,
            "1",
            str(nframes),
            "1",
            "2",
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        raise SystemExit(done.stderr)
    rows = []
    for line in done.stdout.splitlines():
        if not line.startswith("frame,"):
            continue
        fields = line.split(",")
        fp_off = int(fields[10])
        fp_on = int(fields[11])
        rows.append(
            {
                "lane": lane,
                "temperature_K": temperature,
                "frame": fields[3],
                "molecules": fields[6],
                "chill_ice": fields[7],
                "dseams_off": fields[8],
                "dseams_on": fields[9],
                "fp_off": str(fp_off),
                "fp_on": str(fp_on),
                "kind_changed": fields[12],
                "holds": "yes" if fp_on <= fp_off else "no",
                "source": str(path),
            }
        )
    if len(rows) != nframes:
        raise SystemExit(f"{lane} returned {len(rows)} frame rows for {nframes} frames")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lib", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--nucleation", type=Path, default=ROOT / "input" / "traj" / "nucleation_f1.lammpstrj")
    parser.add_argument("--example", type=Path, default=ROOT / "input" / "traj" / "exampleTraj.lammpstrj")
    parser.add_argument("--tip4p", type=Path, required=True)
    parser.add_argument("--mw", type=Path, required=True)
    args = parser.parse_args()
    library = Library(args.lib)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    # --- nucleation, frame 1, 4096 mW molecules ---
    nuc_ids, nuc_pub = library.labels(args.nucleation, 1, 1, "chill+")
    nuc_ids_o, nuc_ovi = library.labels(args.nucleation, 1, 1, "ovito-chill+")
    if not np.array_equal(nuc_ids, nuc_ids_o):
        raise SystemExit("nucleation ids moved between rules")
    nuc_dump = oxygen_frame(load_dump(args.nucleation)[0], 1)
    if not np.array_equal(nuc_dump["ids"], nuc_ids):
        raise SystemExit("nucleation id order differs from the library")
    pub_py, _ = classify_cloud(nuc_dump["xyz"], nuc_dump["lengths"], "published")
    ovi_py, _ = classify_cloud(nuc_dump["xyz"], nuc_dump["lengths"], "ovito")
    none_py, _ = classify_cloud(nuc_dump["xyz"], nuc_dump["lengths"], "notest")
    shell_py, _ = classify_cloud(nuc_dump["xyz"], nuc_dump["lengths"], "shellclosed")
    require_match("nucleation published", pub_py, nuc_pub, nuc_ids)
    require_match("nucleation ovito-chill+", ovi_py, nuc_ovi, nuc_ids)
    ovito_sw = ovito_labels(args.nucleation, 1)[0]
    if set(ovito_sw) != set(int(i) for i in nuc_ids):
        raise SystemExit("OVITO nucleation ids do not match the dump")
    sw = np.array([ovito_sw[int(i)] for i in nuc_ids], dtype=np.int32)
    disagree_rows = []
    cause_counts = {"neighbour_test": 0, "q3_shell": 0, "open_edge": 0, "mixed": 0}
    software_extra = []
    for i, atom_id in enumerate(nuc_ids):
        cause = cause_of(int(pub_py[i]), int(none_py[i]), int(shell_py[i]), int(ovi_py[i]))
        if int(sw[i]) != int(nuc_ovi[i]):
            software_extra.append(i)
        if cause == "agree" and int(sw[i]) == int(nuc_pub[i]):
            continue
        if cause != "agree":
            cause_counts[cause] += 1
        disagree_rows.append(
            {
                "id": int(atom_id),
                "index": i,
                "published": NAMES[int(nuc_pub[i])],
                "ovito_rule": NAMES[int(nuc_ovi[i])],
                "ovito_3_17_2": NAMES[int(sw[i])],
                "cause": cause if int(sw[i]) == int(nuc_ovi[i]) else "ovito_software",
                "shell": "",
            }
        )
    # Shell length of each disagreeing atom, from the full neighbour list.
    full = shells(nuc_dump["xyz"], nuc_dump["lengths"])
    for row in disagree_rows:
        row["shell"] = str(len(full[int(row["index"])]))
    n_dis = sum(cause_counts.values())
    write_csv(
        out / "nucleation-disagreements.csv",
        disagree_rows,
        ["id", "index", "published", "ovito_rule", "ovito_3_17_2", "cause", "shell"],
    )

    # --- exampleTraj oxygen recount ---
    ex_ids, ex_pub = library.labels(args.example, 1, 2, "chill+")
    ex_ids_o, ex_ovi = library.labels(args.example, 1, 2, "ovito-chill+")
    if not np.array_equal(ex_ids, ex_ids_o):
        raise SystemExit("example ids moved between rules")
    ex_dump = oxygen_frame(load_dump(args.example)[0], 2)
    if not np.array_equal(ex_dump["ids"], ex_ids):
        raise SystemExit("example id order differs from the library")
    ex_pub_py, ex_pub_c = classify_cloud(ex_dump["xyz"], ex_dump["lengths"], "published")
    ex_ovi_py, ex_ovi_c = classify_cloud(ex_dump["xyz"], ex_dump["lengths"], "ovito")
    ex_old, ex_old_c = classify_cloud(ex_dump["xyz"], ex_dump["lengths"], "pregate")
    ex_none, _ = classify_cloud(ex_dump["xyz"], ex_dump["lengths"], "notest")
    ex_shell, _ = classify_cloud(ex_dump["xyz"], ex_dump["lengths"], "shellclosed")
    require_match("example published", ex_pub_py, ex_pub, ex_ids)
    require_match("example ovito-chill+", ex_ovi_py, ex_ovi, ex_ids)
    ex_sw = ovito_labels(args.example, 2)[0]
    ex_sw_arr = np.array([ex_sw[int(i)] for i in ex_ids], dtype=np.int32)
    pub_hist = histogram(ex_pub)
    old_hist = histogram(ex_old)
    moved = []
    for i, atom_id in enumerate(ex_ids):
        if int(ex_old[i]) == int(ex_pub[i]) and int(ex_sw_arr[i]) == int(ex_ovi[i]) and int(ex_pub[i]) == int(ex_ovi[i]):
            continue
        moved.append(
            {
                "id": int(atom_id),
                "index": i,
                "shell": str(len(shells(ex_dump["xyz"], ex_dump["lengths"])[i])),
                "pregate": NAMES[int(ex_old[i])],
                "published": NAMES[int(ex_pub[i])],
                "ovito_rule": NAMES[int(ex_ovi[i])],
                "ovito_3_17_2": NAMES[int(ex_sw_arr[i])],
                "pregate_c": " ".join(f"{v:.6f}" for v in ex_old_c[i]),
                "published_c": " ".join(f"{v:.6f}" for v in ex_pub_c[i]),
                "ovito_c": " ".join(f"{v:.6f}" for v in ex_ovi_c[i]),
                "cause_vs_ovito_rule": cause_of(
                    int(ex_pub[i]), int(ex_none[i]), int(ex_shell[i]), int(ex_ovi[i])
                ),
            }
        )
    # The shell column above recomputes the neighbour list per atom. Replace
    # it with one pass so the recorded shell is the same list as the labels.
    ex_full = shells(ex_dump["xyz"], ex_dump["lengths"])
    for row in moved:
        row["shell"] = str(len(ex_full[int(row["index"])]))
    write_csv(
        out / "example-recount.csv",
        moved,
        [
            "id",
            "index",
            "shell",
            "pregate",
            "published",
            "ovito_rule",
            "ovito_3_17_2",
            "pregate_c",
            "published_c",
            "ovito_c",
            "cause_vs_ovito_rule",
        ],
    )

    # --- padding on the generated interface trajectories ---
    pad_rows = []
    pad_rows.extend(padding_frames(args.probe, args.tip4p, "tip4p-ice-basal", "230"))
    pad_rows.extend(padding_frames(args.probe, args.mw, "mw-interface", "270"))
    write_csv(
        out / "padding-per-frame.csv",
        pad_rows,
        [
            "lane",
            "temperature_K",
            "frame",
            "molecules",
            "chill_ice",
            "dseams_off",
            "dseams_on",
            "fp_off",
            "fp_on",
            "kind_changed",
            "holds",
            "source",
        ],
    )

    # Interface frames: per-atom published vs OVITO software.
    iface_rows = []
    iface_causes = {"neighbour_test": 0, "q3_shell": 0, "open_edge": 0, "mixed": 0, "ovito_software": 0}
    for lane, path, type_filter in (
        ("tip4p-ice-basal", args.tip4p, 1),
        ("mw-interface", args.mw, 1),
    ):
        software = ovito_labels(path, type_filter)
        dumps = load_dump(path)
        for frame, dump in enumerate(dumps, start=1):
            ids, pub = library.labels(path, frame, type_filter, "chill+")
            ids2, ovi = library.labels(path, frame, type_filter, "ovito-chill+")
            if not np.array_equal(ids, ids2):
                raise SystemExit(f"{lane} frame {frame} ids moved")
            oxy = oxygen_frame(dump, type_filter)
            if not np.array_equal(oxy["ids"], ids):
                raise SystemExit(f"{lane} frame {frame} id order differs")
            pub_py, _ = classify_cloud(oxy["xyz"], oxy["lengths"], "published")
            ovi_py, _ = classify_cloud(oxy["xyz"], oxy["lengths"], "ovito")
            none_py, _ = classify_cloud(oxy["xyz"], oxy["lengths"], "notest")
            shell_py, _ = classify_cloud(oxy["xyz"], oxy["lengths"], "shellclosed")
            require_match(f"{lane} frame {frame} published", pub_py, pub, ids)
            require_match(f"{lane} frame {frame} ovito-chill+", ovi_py, ovi, ids)
            neigh = shells(oxy["xyz"], oxy["lengths"])
            swmap = software[frame - 1]
            for i, atom_id in enumerate(ids):
                soft = swmap[int(atom_id)]
                cause = cause_of(int(pub_py[i]), int(none_py[i]), int(shell_py[i]), int(ovi_py[i]))
                if soft != int(ovi[i]):
                    cause = "ovito_software"
                if cause == "agree":
                    continue
                iface_causes[cause] += 1
                iface_rows.append(
                    {
                        "lane": lane,
                        "frame": frame,
                        "id": int(atom_id),
                        "shell": str(len(neigh[i])),
                        "published": NAMES[int(pub[i])],
                        "ovito_rule": NAMES[int(ovi[i])],
                        "ovito_3_17_2": NAMES[soft],
                        "cause": cause,
                    }
                )
    write_csv(
        out / "interface-label-disagreements.csv",
        iface_rows,
        ["lane", "frame", "id", "shell", "published", "ovito_rule", "ovito_3_17_2", "cause"],
    )

    pregate_only = [row for row in moved if row["pregate"] != row["published"]]
    software_mismatch = [row for row in moved if row["ovito_rule"] != row["ovito_3_17_2"]]
    nuc_software_mismatch = int(np.sum(sw != nuc_ovi))
    summary = {
        "nucleation_n": str(len(nuc_ids)),
        "nucleation_disagree_published_vs_ovito_rule": str(n_dis),
        "nucleation_neighbour_test": str(cause_counts["neighbour_test"]),
        "nucleation_q3_shell": str(cause_counts["q3_shell"]),
        "nucleation_open_edge": str(cause_counts["open_edge"]),
        "nucleation_mixed": str(cause_counts["mixed"]),
        "nucleation_ovito_software_vs_ovito_rule": str(nuc_software_mismatch),
        "example_n": str(len(ex_ids)),
        "example_published_interClathrate": str(pub_hist["interClathrate"]),
        "example_published_water": str(pub_hist["water"]),
        "example_pregate_interClathrate": str(old_hist["interClathrate"]),
        "example_pregate_water": str(old_hist["water"]),
        "example_pregate_ids": " ".join(str(row["id"]) for row in pregate_only),
        "example_ovito_software_mismatches": str(len(software_mismatch)),
        "padding_frames": str(len(pad_rows)),
        "padding_frames_fail": str(sum(row["holds"] == "no" for row in pad_rows)),
        "interface_label_rows": str(len(iface_rows)),
        "interface_neighbour_test": str(iface_causes["neighbour_test"]),
        "interface_q3_shell": str(iface_causes["q3_shell"]),
        "interface_open_edge": str(iface_causes["open_edge"]),
        "interface_mixed": str(iface_causes["mixed"]),
        "interface_ovito_software": str(iface_causes["ovito_software"]),
    }
    write_csv(out / "summary.csv", [summary], list(summary.keys()))
    print(" ".join(f"{key}={value}" for key, value in summary.items()))


if __name__ == "__main__":
    main()
