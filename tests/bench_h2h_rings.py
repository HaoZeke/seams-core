#!/usr/bin/env python3
"""Head-to-head ring-statistics timings for d-SEAMS vs OVITO.

Protocol mirrors neighbour timing:
- fresh process per sample
- interleaved order with per-repetition rotation
- n >= 6 from confidence target
- median as mean of the middle pair
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

MIN_REPEATS = 6
WARMUP_CALLS = 2
HEXAGONAL_RING_SIZE = 6


@dataclass(frozen=True)
class Frame:
    positions: np.ndarray
    H: np.ndarray
    origin: np.ndarray
    is_triclinic: bool
    source: str


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    path: str
    frame: int
    atom_type: int
    rep_xyz: tuple[int, int, int]


def quantile_linear(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    count = len(ordered)
    if count < MIN_REPEATS:
        raise ValueError(f"need n >= {MIN_REPEATS}, got {count}")
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def summarize_samples(samples_ms: list[float]) -> dict[str, float]:
    if len(samples_ms) < 2:
        raise ValueError("need at least two samples")
    s = sorted(samples_ms)
    mid = len(s) // 2
    if len(s) % 2 == 0:
        median = 0.5 * (s[mid - 1] + s[mid])
    else:
        median = s[mid]
    q1 = quantile_linear(s, 0.25)
    q3 = quantile_linear(s, 0.75)
    return {
        "median_ms": float(median),
        "q1_ms": float(q1),
        "q3_ms": float(q3),
        "iqr_ms": float(q3 - q1),
    }


def derive_repeats(confidence_target: float) -> int:
    base = int(math.ceil(1.0 - math.log2(1.0 - confidence_target)))
    n = max(MIN_REPEATS, base)
    if n % 2 != 0:
        n += 1
    return n


def derive_warmup(repeats: int) -> int:
    if repeats < WARMUP_CALLS:
        raise ValueError(f"repeats={repeats} must be >= warmup calls {WARMUP_CALLS}")
    # One cold-start sink plus one stabilization call.
    return WARMUP_CALLS


def self_check_summary_rules() -> None:
    even = summarize_samples([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    if abs(even["median_ms"] - 3.5) > 1e-12:
        raise AssertionError(f"even-n median regression: {even['median_ms']}")
    odd = summarize_samples([1.0, 2.0, 3.0, 4.0, 9.0])
    if abs(odd["median_ms"] - 3.0) > 1e-12:
        raise AssertionError(f"odd-n median regression: {odd['median_ms']}")


def rotate(values: list[str], shift: int) -> list[str]:
    if not values:
        return []
    k = shift % len(values)
    return values[k:] + values[:k]


def display_cmd(cmd: list[str]) -> str:
    cwd = os.getcwd()
    shown = []
    for part in cmd:
        if part.startswith(cwd + os.sep):
            shown.append(os.path.relpath(part, cwd))
        else:
            shown.append(part)
    return " ".join(shown)


def pinned(cmd: list[str]) -> list[str]:
    return ["taskset", "-c", "0", *cmd]


def git_sha() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()


def git_dirty() -> str:
    unstaged = subprocess.call(
        ["git", "diff", "--quiet", "--ignore-submodules", "HEAD"]
    )
    staged = subprocess.call(
        ["git", "diff", "--cached", "--quiet", "--ignore-submodules", "HEAD"]
    )
    return "dirty" if unstaged != 0 or staged != 0 else "clean"


def published_cpu_model(name: str) -> str:
    """Drop the register mark so the record reads 'Intel Xeon Processor'."""
    return " ".join(name.replace("(R)", "").split())


def cpu_model_name() -> str:
    try:
        with Path("/proc/cpuinfo").open("r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("model name"):
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        return published_cpu_model(parts[1].strip())
    except OSError:
        pass
    return "unknown"


def runtime_machine_header() -> dict[str, object]:
    return {
        "uname_m": platform.machine(),
        "nproc": int(os.cpu_count() or 1),
        "cpu_model": cpu_model_name(),
    }


def parse_box(bounds_header: str, bounds_lines: list[str]) -> tuple[np.ndarray, np.ndarray, bool]:
    tric = "xy" in bounds_header.lower()
    if not tric:
        xlo, xhi = map(float, bounds_lines[0].split()[:2])
        ylo, yhi = map(float, bounds_lines[1].split()[:2])
        zlo, zhi = map(float, bounds_lines[2].split()[:2])
        H = np.array(
            [[xhi - xlo, 0.0, 0.0], [0.0, yhi - ylo, 0.0], [0.0, 0.0, zhi - zlo]],
            dtype=np.float64,
        )
        return H, np.array([xlo, ylo, zlo], dtype=np.float64), False
    xlo_b, xhi_b, xy = map(float, bounds_lines[0].split()[:3])
    ylo_b, yhi_b, xz = map(float, bounds_lines[1].split()[:3])
    zlo_b, zhi_b, yz = map(float, bounds_lines[2].split()[:3])
    xlo = xlo_b - min(0.0, xy, xz, xy + xz)
    xhi = xhi_b - max(0.0, xy, xz, xy + xz)
    ylo = ylo_b - min(0.0, yz)
    yhi = yhi_b - max(0.0, yz)
    H = np.array(
        [[xhi - xlo, xy, xz], [0.0, yhi - ylo, yz], [0.0, 0.0, zhi_b - zlo_b]],
        dtype=np.float64,
    )
    return H, np.array([xlo, ylo, zlo_b], dtype=np.float64), True


def read_lammps_frame(path: str, target_frame: int) -> Frame:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(path)
    with p.open("r", encoding="utf-8") as fh:
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"frame {target_frame} not found in {path}")
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(fh.readline().strip())
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"malformed dump: {path}")
            n_atoms = int(fh.readline().strip())
            header = fh.readline().strip()
            bounds = [fh.readline().strip(), fh.readline().strip(), fh.readline().strip()]
            atoms_header = fh.readline().strip()
            cols = atoms_header.split()[2:]
            colmap = {name: idx for idx, name in enumerate(cols)}
            if timestep != target_frame and target_frame != 1:
                for _ in range(n_atoms):
                    fh.readline()
                continue
            pts = np.empty((n_atoms, 3), dtype=np.float64)
            ids = np.empty(n_atoms, dtype=np.int64)
            for i in range(n_atoms):
                row = fh.readline().split()
                ids[i] = int(row[colmap["id"]])
                pts[i, 0] = float(row[colmap["x"]])
                pts[i, 1] = float(row[colmap["y"]])
                pts[i, 2] = float(row[colmap["z"]])
            order = np.argsort(ids, kind="mergesort")
            pts = pts[order]
            H, origin, tric = parse_box(header, bounds)
            return Frame(
                positions=pts,
                H=H,
                origin=origin,
                is_triclinic=tric,
                source=f"{path}@{target_frame}",
            )


def replicate_frame(frame: Frame, nx: int, ny: int, nz: int) -> Frame:
    if nx == 1 and ny == 1 and nz == 1:
        return frame
    a = frame.H[:, 0]
    b = frame.H[:, 1]
    c = frame.H[:, 2]
    offsets = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                offsets.append(ix * a + iy * b + iz * c)
    pts = np.concatenate([frame.positions + off for off in offsets], axis=0)
    H = frame.H.copy()
    H[:, 0] *= nx
    H[:, 1] *= ny
    H[:, 2] *= nz
    return Frame(
        positions=pts,
        H=H,
        origin=frame.origin.copy(),
        is_triclinic=frame.is_triclinic,
        source=f"{frame.source} rep={nx}x{ny}x{nz}",
    )


def write_lammps_dump(path: str, frame: Frame, atom_type: int) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    H = frame.H
    lx = float(H[0, 0])
    xy = float(H[0, 1])
    xz = float(H[0, 2])
    ly = float(H[1, 1])
    yz = float(H[1, 2])
    lz = float(H[2, 2])
    xlo = float(frame.origin[0])
    ylo = float(frame.origin[1])
    zlo = float(frame.origin[2])
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    xlo_b = xlo + xmin
    xhi_b = xlo + lx + xmax
    ylo_b = ylo + ymin
    yhi_b = ylo + ly + ymax
    zlo_b = zlo
    zhi_b = zlo + lz
    with p.open("w", encoding="utf-8") as fh:
        fh.write("ITEM: TIMESTEP\n1\n")
        fh.write("ITEM: NUMBER OF ATOMS\n")
        fh.write(f"{frame.positions.shape[0]}\n")
        fh.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        fh.write(f"{xlo_b:.12f} {xhi_b:.12f} {xy:.12f}\n")
        fh.write(f"{ylo_b:.12f} {yhi_b:.12f} {xz:.12f}\n")
        fh.write(f"{zlo_b:.12f} {zhi_b:.12f} {yz:.12f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, pos in enumerate(frame.positions, start=1):
            fh.write(f"{i:d} {atom_type:d} {pos[0]:.12f} {pos[1]:.12f} {pos[2]:.12f}\n")


def save_frame_npz(frame: Frame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        positions=frame.positions,
        H=frame.H,
        origin=frame.origin,
        triclinic=np.array([1 if frame.is_triclinic else 0], dtype=np.int8),
        source=np.array([frame.source], dtype=object),
    )


def load_frame_npz(path: str) -> Frame:
    with np.load(path, allow_pickle=True) as data:
        return Frame(
            positions=data["positions"],
            H=data["H"],
            origin=data["origin"],
            is_triclinic=bool(int(data["triclinic"][0])),
            source=str(data["source"][0]),
        )


def parse_bench_stages(out: str) -> tuple[int, int, float, float]:
    rings_pat = re.compile(r"^# .* rings ([0-9]+) six ([0-9]+) ddc ([0-9]+).*$", re.MULTILINE)
    neigh_pat = re.compile(r"^neigh_ms ([0-9.eE+-]+)$", re.MULTILINE)
    ring_pat = re.compile(r"^rings_ms ([0-9.eE+-]+)$", re.MULTILINE)
    m_rings = rings_pat.search(out)
    m_neigh = neigh_pat.search(out)
    m_ring = ring_pat.search(out)
    if not (m_rings and m_neigh and m_ring):
        raise RuntimeError(f"cannot parse bench_stages output\n{out}")
    return int(m_rings.group(1)), int(m_rings.group(2)), float(m_neigh.group(1)), float(m_ring.group(1))


def run_dseams_once(
    bench_stages_bin: str,
    traj: str,
    frame: int,
    atom_type: int,
    threads: int,
    git_sha_value: str,
    dirty_value: str,
) -> dict[str, object]:
    cmd = pinned([bench_stages_bin, traj, str(frame), str(atom_type), "1"])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["SEAMS_GIT_SHA"] = git_sha_value
    env["SEAMS_GIT_DIRTY"] = dirty_value
    out = subprocess.check_output(cmd, text=True, env=env)
    rings_total, rings_six, neigh_ms, ring_ms = parse_bench_stages(out)
    cmd_text = (
        f"OMP_NUM_THREADS={threads} OMP_DYNAMIC=FALSE "
        f"SEAMS_GIT_SHA={git_sha_value} "
        f"SEAMS_GIT_DIRTY={dirty_value} "
        + display_cmd(cmd)
    )
    return {
        "command": cmd_text,
        "rings_total": rings_total,
        "rings_size6": rings_six,
        "ring_ms": ring_ms,
        "total_ms": neigh_ms + ring_ms,
    }


def run_ovito_worker(
    script_path: str,
    frame_npz: str,
    cutoff: float,
    max_ring_size: int,
    ovito_threads: int,
) -> dict[str, object]:
    cmd = [
        sys.executable,
        script_path,
        "--worker-ovito",
        "--frame-npz",
        frame_npz,
        "--cutoff",
        f"{cutoff:.8f}",
        "--max-ring-size",
        str(max_ring_size),
        "--ovito-threads",
        str(ovito_threads),
    ]
    pinned_cmd = pinned(cmd)
    out = json.loads(subprocess.check_output(pinned_cmd, text=True))
    out["command"] = display_cmd(pinned_cmd)
    return out


def run_interleaved(
    tools: list[str],
    warmup_runs: int,
    repeats: int,
    run_tool: Callable[[str], dict[str, object]],
    label: str,
) -> dict[str, object]:
    records: dict[str, dict[str, object]] = {}
    active: list[str] = []

    for tool in tools:
        res = run_tool(tool)
        if "unsupported" in res:
            continue
        active.append(tool)
        records[tool] = {
            "first_call_ms": float(res["first_ms"]),
            "warmup_runs": warmup_runs,
            "repeats": repeats,
            "samples_ms": [],
            "count": int(res["count"]),
            "command": res.get("command"),
        }

    if not active:
        return {"unsupported": "no active tools"}

    for warm in range(warmup_runs):
        for tool in rotate(active, warm):
            res = run_tool(tool)
            if "unsupported" in res:
                raise RuntimeError(f"{label} warm-up became unsupported for {tool}")
            if int(res["count"]) != int(records[tool]["count"]):
                raise RuntimeError(f"{label} warm-up count drift for {tool}")

    for rep in range(repeats):
        for tool in rotate(active, rep):
            res = run_tool(tool)
            if "unsupported" in res:
                raise RuntimeError(f"{label} sample became unsupported for {tool}")
            if int(res["count"]) != int(records[tool]["count"]):
                raise RuntimeError(f"{label} sample count drift for {tool}")
            records[tool]["samples_ms"].append(float(res["sample_ms"]))

    for tool in active:
        records[tool].update(summarize_samples(records[tool]["samples_ms"]))

    return records


def worker_ovito(args: argparse.Namespace) -> None:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    os.environ["OVITO_THREAD_COUNT"] = str(args.ovito_threads)
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    frame = load_frame_npz(args.frame_npz)
    data = DataCollection()
    data.create_cell(frame.H, pbc=(True, True, True))
    parts = data.create_particles(count=frame.positions.shape[0])
    parts.create_property("Position", data=frame.positions - frame.origin)

    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=args.cutoff))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=args.max_ring_size))

    t0 = time.perf_counter()
    out = pipe.compute()
    t1 = time.perf_counter()
    hist = out.tables["ring-size-histogram"]
    sizes = hist["Ring Size"][...]
    counts = hist["Count"][...]
    total = int(np.sum(counts))
    six = (
        int(counts[np.where(sizes == HEXAGONAL_RING_SIZE)[0][0]])
        if np.any(sizes == HEXAGONAL_RING_SIZE)
        else 0
    )
    print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": total, "rings_size6": six}))


def main() -> None:
    self_check_summary_rules()

    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-ovito", action="store_true")
    ap.add_argument("--frame-npz")
    ap.add_argument("--bench-stages-bin")
    ap.add_argument("--output-json")
    ap.add_argument("--cutoff", type=float, default=3.5)
    ap.add_argument("--max-ring-size", type=int, default=HEXAGONAL_RING_SIZE)
    ap.add_argument("--confidence", type=float, default=0.95)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--dseams-threads", type=int, default=1)
    ap.add_argument("--ovito-threads", type=int, default=1)
    args = ap.parse_args()

    if args.worker_ovito:
        worker_ovito(args)
        return

    if not args.bench_stages_bin or not args.output_json:
        raise SystemExit("--bench-stages-bin and --output-json are required")

    sha = git_sha()
    dirty = git_dirty()
    repeats = derive_repeats(args.confidence)
    warmup = derive_warmup(repeats)
    machine = runtime_machine_header()

    print(
        f"# header git_sha {sha} dirty {dirty} seed {args.seed} repeats {repeats} warmup {warmup}",
        flush=True,
    )
    print(
        "# machine "
        f"uname_m {machine['uname_m']} nproc {machine['nproc']} cpu_model {machine['cpu_model']}",
        flush=True,
    )

    generated_specs: list[DatasetSpec] = []
    for name, path, rep in [
        ("orth-1e5", "input/traj/mW_cubic.lammpstrj", (3, 3, 3)),
        ("tri-1e5", "input/traj/genice_sH.lammpstrj", (10, 10, 10)),
    ]:
        base = read_lammps_frame(path, 1)
        rep_frame = replicate_frame(base, *rep)
        out_path = f"input/traj/{name}-ring-scale.lammpstrj"
        write_lammps_dump(out_path, rep_frame, atom_type=1)
        generated_specs.append(DatasetSpec(name, out_path, 1, 1, (1, 1, 1)))

    datasets = [
        DatasetSpec("orth-4k", "input/traj/mW_cubic.lammpstrj", 1, 1, (1, 1, 1)),
        DatasetSpec("tri-136", "input/traj/genice_sH.lammpstrj", 1, 1, (1, 1, 1)),
        *generated_specs,
    ]

    out: dict[str, object] = {
        "git_sha": sha,
        "git_dirty": dirty,
        "seed": args.seed,
        "cutoff": args.cutoff,
        "max_ring_size": args.max_ring_size,
        "machine": machine,
        "protocol": {
            "confidence_target": args.confidence,
            "repeat_lower_bound": "n >= ceil(1 - log2(1 - confidence_target))",
            "median_rule": "type 7 median and quartiles",
            "repeats": repeats,
            "warmup_runs": warmup,
            "warmup_rule": (
                "fixed 2 calls (one cold-start sink plus one stabilization call)"
            ),
            "fresh_process_per_sample": True,
            "interleaving": "per repetition, rotate tool order",
            "first_call_reported_separately": True,
        },
        "thread_policy": {
            "dseams_threads_reported": [1, args.dseams_threads],
            "ovito_threads": args.ovito_threads,
        },
        "datasets": {},
    }

    frame_store = Path("input/.bench_frames")
    frame_store.mkdir(parents=True, exist_ok=True)

    for spec in datasets:
        frame = read_lammps_frame(spec.path, spec.frame)
        n_atoms = int(frame.positions.shape[0])
        frame_npz = str(frame_store / f"ring-{spec.name}.npz")
        save_frame_npz(frame, frame_npz)

        print(f"[ring] {spec.name} n={n_atoms}", flush=True)

        def run_ring_tool(tool: str) -> dict[str, object]:
            if tool.startswith("dseams"):
                threads = 1 if tool.endswith("omp1") else args.dseams_threads
                res = run_dseams_once(
                    args.bench_stages_bin,
                    spec.path,
                    spec.frame,
                    spec.atom_type,
                    threads,
                    sha,
                    dirty,
                )
                return {
                    "first_ms": float(res["ring_ms"]),
                    "sample_ms": float(res["ring_ms"]),
                    "count": int(res["rings_total"]),
                    "command": res["command"],
                }

            ov = run_ovito_worker(
                script_path=str(Path(__file__)),
                frame_npz=frame_npz,
                cutoff=args.cutoff,
                max_ring_size=args.max_ring_size,
                ovito_threads=args.ovito_threads,
            )
            return {
                "first_ms": float(ov["ms"]),
                "sample_ms": float(ov["ms"]),
                "count": int(ov["count"]),
                "command": ov["command"],
            }

        ring_stage = run_interleaved(
            ["dseams_omp1", "ovito"],
            warmup,
            repeats,
            run_ring_tool,
            f"{spec.name}/ring",
        )
        dseams_count = int(ring_stage["dseams_omp1"]["count"])
        ovito_count = int(ring_stage["ovito"]["count"])
        if dseams_count != ovito_count:
            raise SystemExit(
                f"{spec.name} ring counts differ: dseams {dseams_count} ovito {ovito_count}"
            )
        print(f"[ring] {spec.name} counts equal {dseams_count}", flush=True)
        for tool in ("dseams_omp1", "ovito"):
            samples = [float(value) for value in ring_stage[tool]["samples_ms"]]
            median = type7(samples, 0.5)
            shown = " ".join(f"{value:.6f}" for value in samples)
            print(
                f"[ring] {spec.name} {tool} n={len(samples)} "
                f"median_ms {median:.6f} samples_ms {shown}",
                flush=True,
            )

        def run_total_tool(tool: str) -> dict[str, object]:
            threads = 1 if tool.endswith("omp1") else args.dseams_threads
            res = run_dseams_once(
                args.bench_stages_bin,
                spec.path,
                spec.frame,
                spec.atom_type,
                threads,
                sha,
                dirty,
            )
            return {
                "first_ms": float(res["total_ms"]),
                "sample_ms": float(res["total_ms"]),
                "count": int(res["rings_total"]),
                "command": res["command"],
            }

        total_stage = run_interleaved(
            ["dseams_omp1"],
            warmup,
            repeats,
            run_total_tool,
            f"{spec.name}/total",
        )

        counts_ref = run_dseams_once(
            args.bench_stages_bin,
            spec.path,
            spec.frame,
            spec.atom_type,
            1,
            sha,
            dirty,
        )

        out["datasets"][spec.name] = {
            "source": frame.source,
            "n_atoms": n_atoms,
            "triclinic": bool(frame.is_triclinic),
            "rings_total": int(counts_ref["rings_total"]),
            "rings_size6": int(counts_ref["rings_size6"]),
            "ring_stage_ms": ring_stage,
            "neigh_plus_ring_ms": total_stage,
        }

    Path(args.output_json).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {args.output_json}")


if __name__ == "__main__":
    main()
