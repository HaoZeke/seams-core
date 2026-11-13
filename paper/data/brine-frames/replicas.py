#!/usr/bin/env python3
"""Per-replica composition, CHILL+ census and ice front of the 230 K brine frames.

The frames are the labelled extended-XYZ files of the Hugging Face
dataset d-SEAMS/brine-front-230K, read at one pinned revision. fetch
downloads them by revision and checks them against the dataset's
SHA256SUMS. Nothing is uploaded.

table writes one row per replica x loading x frame for the expected grid
(replicas r0 and r1; 0, 10, 20 and 40 NaCl pairs; t0, t1, t2 and tend).
A frame that the revision does not hold is a row with status
not_in_dataset and empty numbers. range writes, for every quantity, each
replica's value and the range over the replicas that have one; the range
is left empty when fewer than two replicas exist.

The ice front is the d-SEAMS locateFronts rule (src/brine_front.cpp):
20 bins along z, a bin is solid when at least half of its waters are
crystal (cubic, hexagonal, reCubic or reHex), and the two fronts are the
linear 0.5 crossings at the ends of the longest solid run of occupied
bins. z_front_hi is where that run leaves ice going up z, folded into
the cell. z_front_lo is where it enters ice, shifted down by one cell
length when the slab wraps through z = 0, so z_front_lo < z_front_hi and
the ice thickness is their difference.

    python3 paper/data/brine-frames/replicas.py fetch \
        --dest d-SEAMS/brine-front-230K
    python3 paper/data/brine-frames/replicas.py table \
        --frames d-SEAMS/brine-front-230K/data/frames \
        --csv paper/data/brine-frames/replicas.csv
    python3 paper/data/brine-frames/replicas.py range \
        --csv paper/data/brine-frames/replicas.csv \
        --out paper/data/brine-frames/replicas_range.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import regenerate as deck  # noqa: E402

HF_DATASET = "d-SEAMS/brine-front-230K"
HF_REVISION = "4828e33ecf419c70643b2d466bfa4f283b4138dc"
REPLICAS = ("r0", "r1")
PAIRS = (0, 10, 20, 40)
FRAMES = (("t0", 0), ("t1", 830_000), ("t2", 1_670_000), ("tend", 2_500_000))
TIMESTEP_FS = 2.0
N_BINS = 20
FRONT_CUT = 0.5
CRYSTAL = {"cubic", "hexagonal", "reCubic", "reHex"}

COLUMNS = (
    "replica",
    "pairs",
    "frame",
    "step",
    "t_ns",
    "status",
    "atoms",
    "waters",
    "na",
    "cl",
    "hexagonal",
    "cubic",
    "interfacial",
    "liquid",
    "clathrate",
    "interClathrate",
    "unclassified",
    "crystal",
    "lx_A",
    "ly_A",
    "lz_A",
    "z_front_lo_A",
    "z_front_hi_A",
    "ice_thickness_A",
    "cage_ice",
    "nmax",
    "ion_ice",
    "ion_front",
    "ion_liquid",
    "velocity_seed",
    "ion_seed",
)
QUANTITIES = (
    "waters",
    "na",
    "cl",
    "hexagonal",
    "cubic",
    "interfacial",
    "liquid",
    "clathrate",
    "interClathrate",
    "crystal",
    "lz_A",
    "z_front_lo_A",
    "z_front_hi_A",
    "ice_thickness_A",
    "cage_ice",
    "nmax",
    "ion_ice",
    "ion_front",
    "ion_liquid",
)


def frame_name(replica: str, pairs: int, frame: str) -> str:
    return f"T230_m{pairs}_{replica}_{frame}.xyz"


def _get(url: str, attempts: int = 4) -> bytes:
    delay = 4.0
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(delay)
            delay *= 2
    raise AssertionError("unreachable")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(dest: Path, dataset: str, revision: str, manifest: Path) -> None:
    """Mirror every file that SHA256SUMS lists at the revision, checked against it.

    SHA256SUMS itself is not written; rule brine_sha256sums writes the
    local one once the new runs are in the same tree.
    """
    base = f"https://huggingface.co/datasets/{dataset}/resolve/{revision}"
    sums = {}
    for line in _get(f"{base}/SHA256SUMS").decode().splitlines():
        digest, _, path = line.strip().partition("  ")
        if path:
            sums[path.lstrip("./")] = digest
    for rel in sorted(sums):
        path = dest / rel
        if not path.exists() or _sha256(path) != sums[rel]:
            blob = _get(f"{base}/{rel}")
            if hashlib.sha256(blob).hexdigest() != sums[rel]:
                raise SystemExit(f"{rel} does not match SHA256SUMS at {revision}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(blob)
        print(f"{rel} ok", flush=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(f"# {dataset} at {revision}\n" + "".join(f"{sums[r]}  {r}\n" for r in sorted(sums)))
    print(f"{len(sums)} files at {revision}")


def sums(dataset: Path) -> None:
    """Write SHA256SUMS in the dataset's format: every file but itself and .gitattributes."""
    skip = {"SHA256SUMS", ".gitattributes"}
    files = sorted(
        p.relative_to(dataset).as_posix()
        for p in dataset.rglob("*")
        if p.is_file() and p.relative_to(dataset).as_posix() not in skip
    )
    lines = [f"{_sha256(dataset / rel)}  {rel}\n" for rel in files]
    (dataset / "SHA256SUMS").write_text("".join(lines))
    print(f"wrote {dataset / 'SHA256SUMS'}: {len(lines)} files")


def parse_xyz(path: Path):
    lines = path.read_text().splitlines()
    header = lines[1]
    step = int(header.split("step=")[1].split()[0])
    origin = [float(v) for v in header.split('origin="')[1].split('"')[0].split()]
    lattice = [float(v) for v in header.split('Lattice="')[1].split('"')[0].split()]
    length = (lattice[0], lattice[4], lattice[8])
    rows = []
    for line in lines[2 : 2 + int(lines[0])]:
        name, _x, _y, z, label = line.split()
        rows.append((name, float(z), label))
    return step, origin, length, rows


def locate_fronts(z_water, ice, z0: float, lz: float, n_bins: int = N_BINS, cut: float = FRONT_CUT):
    """Port of site::concentrationProfile and site::locateFronts along z."""
    width = lz / n_bins
    n_water = [0] * n_bins
    n_ice = [0] * n_bins
    for z, flag in zip(z_water, ice):
        z = z0 + ((z - z0) % lz)
        b = min(max(int(math.floor((z - z0) / width)), 0), n_bins - 1)
        n_water[b] += 1
        n_ice[b] += 1 if flag else 0
    mid = [z0 + width * (i + 0.5) for i in range(n_bins)]
    frac = [n_ice[i] / n_water[i] if n_water[i] else 0.0 for i in range(n_bins)]
    occupied = [i for i in range(n_bins) if n_water[i] > 0]
    m = len(occupied)
    solid = [frac[occupied[i]] >= cut for i in range(m)]
    if m == 0 or not any(solid):
        return None
    if all(solid):
        return z0, z0, lz
    best_start, best_len = 0, -1
    for i in range(m):
        if not solid[i] or solid[(i + m - 1) % m]:
            continue
        run = 0
        while run < m and solid[(i + run) % m]:
            run += 1
        if run > best_len:
            best_start, best_len = i, run

    def fold(z):
        return z0 + ((z - z0) % lz)

    def crossing(left, right):
        lb, rb = occupied[left], occupied[right]
        zl, zr = mid[lb], mid[rb]
        if zr < zl:
            zr += lz
        denom = frac[rb] - frac[lb]
        f = 0.5 if abs(denom) < 1e-12 else (cut - frac[lb]) / denom
        f = min(max(f, 0.0), 1.0)
        return fold(zl + f * (zr - zl))

    enter_solid = best_start
    enter_liquid = (best_start + m - 1) % m
    exit_solid = (best_start + best_len - 1) % m
    exit_liquid = (exit_solid + 1) % m
    z_lo = crossing(enter_liquid, enter_solid)
    z_hi = crossing(exit_solid, exit_liquid)
    if z_lo > z_hi:
        z_lo -= lz
    return z_lo, z_hi, z_hi - z_lo


def cage_columns(path: Path, seams: str) -> dict:
    """Campaign observables on one frame: the largest seeded and completed cage
    cluster (K=4, candidate 5.5 A) and the ions with an all-ice, mixed or
    liquid 3.5 A shell, as regenerate.series computes them."""
    lines = path.read_text().splitlines()
    header = lines[1]
    step = int(header.split("step=")[1].split()[0])
    origin = [float(v) for v in header.split('origin="')[1].split('"')[0].split()]
    lattice = [float(v) for v in header.split('Lattice="')[1].split('"')[0].split()]
    length = [lattice[0], lattice[4], lattice[8]]
    species = {"O": 1, "H": 2, "Na": 3, "Cl": 4}
    atoms = []
    for index, line in enumerate(lines[2 : 2 + int(lines[0])], start=1):
        name, x, y, z, _label = line.split()
        atoms.append((index, species[name], float(x), float(y), float(z)))
    with tempfile.TemporaryDirectory() as tmp:
        one = Path(tmp) / "frame.lammpstrj"
        deck.write_one_frame_dump(one, step, origin, length, atoms)
        flags, _text = deck.cage_flags(seams, one, Path(tmp) / "cage.lammpstrj")
    oxy, ice, ions = [], [], []
    for aid, typ, x, y, z in atoms:
        if typ == 1:
            oxy.append((x, y, z))
            ice.append(flags.get(aid, 0) != 0)
        elif typ in (3, 4):
            ions.append((x, y, z))
    n_ice, n_front, n_liq = deck.ions_in_liquid(oxy, ice, ions, length)
    return {
        "cage_ice": sum(ice),
        "nmax": deck.largest_cluster(oxy, length, ice, k=4),
        "ion_ice": n_ice,
        "ion_front": n_front,
        "ion_liquid": n_liq,
    }


def analyse(path: Path, replica: str, pairs: int, frame: str, seams: str | None = None) -> dict:
    step, origin, length, rows = parse_xyz(path)
    labels = [label for name, _z, label in rows if name == "O"]
    z_water = [z for name, z, _label in rows if name == "O"]
    count = {label: labels.count(label) for label in set(labels)}
    ice = [label in CRYSTAL for label in labels]
    fronts = locate_fronts(z_water, ice, origin[2], length[2])
    row = {
        "replica": replica,
        "pairs": pairs,
        "frame": frame,
        "step": step,
        "t_ns": f"{step * TIMESTEP_FS * 1e-6:.4f}",
        "status": "ok",
        "atoms": len(rows),
        "waters": len(labels),
        "na": sum(1 for name, _z, _l in rows if name == "Na"),
        "cl": sum(1 for name, _z, _l in rows if name == "Cl"),
        "hexagonal": count.get("hexagonal", 0) + count.get("reHex", 0),
        "cubic": count.get("cubic", 0) + count.get("reCubic", 0),
        "interfacial": count.get("interfacial", 0),
        "liquid": count.get("water", 0),
        "clathrate": count.get("clathrate", 0),
        "interClathrate": count.get("interClathrate", 0),
        "unclassified": count.get("unclassified", 0),
        "crystal": sum(ice),
        "lx_A": f"{length[0]:.4f}",
        "ly_A": f"{length[1]:.4f}",
        "lz_A": f"{length[2]:.4f}",
    }
    known = sum(row[k] for k in ("hexagonal", "cubic", "interfacial", "liquid", "clathrate", "interClathrate", "unclassified"))
    if known != row["waters"]:
        raise SystemExit(f"{path.name}: {row['waters'] - known} oxygens carry an unknown label")
    if row["na"] != pairs or row["cl"] != pairs:
        raise SystemExit(f"{path.name}: {row['na']} Na and {row['cl']} Cl, expected {pairs} each")
    if fronts is not None:
        row["z_front_lo_A"] = f"{fronts[0]:.3f}"
        row["z_front_hi_A"] = f"{fronts[1]:.3f}"
        row["ice_thickness_A"] = f"{fronts[2]:.3f}"
    if seams:
        row.update(cage_columns(path, seams))
    return row


def seed_columns(replica: str) -> dict:
    number = int(replica[1:])
    return {"velocity_seed": int(deck.velocity_seed(number)), "ion_seed": number}


def table(frames: Path, out: Path, seams: str | None = None) -> None:
    rows = []
    missing = []
    for replica in REPLICAS:
        for pairs in PAIRS:
            for frame, step in FRAMES:
                path = frames / frame_name(replica, pairs, frame)
                if path.exists():
                    row = analyse(path, replica, pairs, frame, seams)
                    row.update(seed_columns(replica))
                    if row["step"] != step:
                        raise SystemExit(f"{path.name} is step {row['step']}, expected {step}")
                else:
                    row = {"replica": replica, "pairs": pairs, "frame": frame, "step": step,
                           "t_ns": f"{step * TIMESTEP_FS * 1e-6:.4f}", "status": "not_in_dataset"}
                    missing.append(path.name)
                rows.append(row)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        handle.write(f"# {HF_DATASET} at {HF_REVISION}, data/frames/T230_m{{pairs}}_{{replica}}_{{frame}}.xyz\n")
        handle.write(f"# front: {N_BINS} z bins, solid when crystal/waters >= {FRONT_CUT}; crystal is cubic+hexagonal\n")
        handle.write("# crystal is the published CHILL+ rule (cubic+hexagonal); cage_ice and nmax are the seeded,\n")
        handle.write("# completed cage labels; ion_* are 3.5 A shells against the cage labels; velocity_seed is\n")
        handle.write("# the in.brine.lmp seed, ion_seed the numpy seed of the ion draw\n")
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}: {len(rows) - len(missing)} frames analysed, {len(missing)} not in the dataset")
    for name in missing:
        print(f"  missing {name}")


def read_rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def value_range(path: Path, out: Path) -> None:
    rows = read_rows(path)
    index = {(r["replica"], int(r["pairs"]), r["frame"]): r for r in rows}
    lines = []
    for pairs in PAIRS:
        for frame, _step in FRAMES:
            for quantity in QUANTITIES:
                values = {}
                for replica in REPLICAS:
                    row = index.get((replica, pairs, frame))
                    if row and row["status"] == "ok" and row.get(quantity, "") != "":
                        values[replica] = float(row[quantity])
                record = {"pairs": pairs, "frame": frame, "quantity": quantity, "n": len(values)}
                for replica in REPLICAS:
                    record[replica] = _fmt(values[replica]) if replica in values else ""
                if len(values) >= 2:
                    record["min"] = _fmt(min(values.values()))
                    record["max"] = _fmt(max(values.values()))
                    record["span"] = _fmt(max(values.values()) - min(values.values()))
                lines.append(record)
        # Table 17 slope: (nmax(5 ns) - nmax(0)) / 5000 ps, per replica.
        values = {}
        for replica in REPLICAS:
            start = index.get((replica, pairs, "t0"))
            end = index.get((replica, pairs, "tend"))
            if start and end and start["status"] == end["status"] == "ok" and start.get("nmax") and end.get("nmax"):
                values[replica] = (float(end["nmax"]) - float(start["nmax"])) / 5000.0
        record = {"pairs": pairs, "frame": "t0-tend", "quantity": "nmax_slope_per_ps", "n": len(values)}
        for replica in REPLICAS:
            record[replica] = f"{values[replica]:.4f}" if replica in values else ""
        if len(values) >= 2:
            record["min"] = f"{min(values.values()):.4f}"
            record["max"] = f"{max(values.values()):.4f}"
            record["span"] = f"{max(values.values()) - min(values.values()):.4f}"
        lines.append(record)
    fields = ("pairs", "frame", "quantity", *REPLICAS, "n", "min", "max", "span")
    with out.open("w", newline="") as handle:
        handle.write(f"# {HF_DATASET} at {HF_REVISION}. Empty replica cell: frame not in the dataset.\n")
        handle.write("# min, max and span are left empty when fewer than two replicas have the frame.\n")
        writer = csv.DictWriter(handle, fieldnames=fields, restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(lines)
    print(f"wrote {out}: {len(lines)} rows")


def _fmt(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--dataset", default=HF_DATASET)
    p.add_argument("--revision", default=HF_REVISION)
    p.add_argument("--manifest", type=Path, required=True)
    p = sub.add_parser("sums")
    p.add_argument("--dataset", type=Path, required=True)
    p = sub.add_parser("table")
    p.add_argument("--frames", type=Path, required=True)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--seams", default=None, help="seams CLI for the cage columns; omitted leaves them empty")
    p = sub.add_parser("range")
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.cmd == "fetch":
        fetch(args.dest, args.dataset, args.revision, args.manifest)
    elif args.cmd == "sums":
        sums(args.dataset)
    elif args.cmd == "table":
        table(args.frames, args.csv, args.seams)
    else:
        value_range(args.csv, args.out)


if __name__ == "__main__":
    sys.exit(main())
