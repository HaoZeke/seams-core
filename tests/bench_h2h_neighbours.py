#!/usr/bin/env python3
"""Head-to-head neighbour-stage benchmarks for d-SEAMS and external tools.

Protocol highlights:
- fresh process per sample for every measured tool path
- interleaved timing order with per-repetition rotation
- n >= 6 samples from 1 - 2^(1-n) >= confidence target
- median reported as mean of middle pair of sorted samples
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
from mw_gr_window import (
    MW_K_CANDIDATE_ANGSTROM,
    MW_K_CANDIDATE_DERIVATION,
)

MIN_REPEATS = 6
WARMUP_CALLS = 2
LINKCELL_FIT_POINTS_PER_GEOMETRY_STAGE = 3
LINKCELL_FIT_SELECTION_RULE = (
    "fit each geometry/stage on the largest three measured sizes; "
    "treat smaller sizes as overhead-dominated"
)

SKEW_FRAME_GRID_NX = 8
SKEW_FRAME_GRID_NY = 8
SKEW_FRAME_GRID_NZ = 8
SKEW_FRAME_LX_ANGSTROM = 24.0
SKEW_FRAME_LY_ANGSTROM = 24.0
SKEW_FRAME_LZ_ANGSTROM = 3.0
SKEW_FRAME_XY_ANGSTROM = 7.2
SKEW_FRAME_XZ_ANGSTROM = 4.8
SKEW_FRAME_YZ_ANGSTROM = 1.5
SKEW_FRAME_JITTER_FRACTION = 0.08
SKEW_FRAME_REFERENCE_CUTOFF_ANGSTROM = 3.5


@dataclass(frozen=True)
class Frame:
    positions: np.ndarray
    types: np.ndarray
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
    check_cutoff_set_equality: bool = False


@dataclass
class Sample:
    first_ms: float
    sample_ms: float
    count: int
    sample_path: str | None = None
    command: str | None = None


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
    # Keep the minimum confidence bound and an even sample count so median
    # stays the mean of the two center samples.
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


def repo_relative(path: Path) -> str:
    repo_root = Path(__file__).resolve().parents[1]
    try:
        return str(path.resolve().relative_to(repo_root))
    except ValueError:
        return str(path)


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
            if not header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"missing BOX BOUNDS in {path}")
            bounds = [fh.readline().strip(), fh.readline().strip(), fh.readline().strip()]
            atoms_header = fh.readline().strip()
            if not atoms_header.startswith("ITEM: ATOMS"):
                raise ValueError(f"missing ATOMS header in {path}")
            cols = atoms_header.split()[2:]
            colmap = {name: idx for idx, name in enumerate(cols)}
            required = {"id", "type", "x", "y", "z"}
            if not required.issubset(colmap):
                raise ValueError(f"frame lacks required columns {required} in {path}")
            if timestep != target_frame and target_frame != 1:
                for _ in range(n_atoms):
                    fh.readline()
                continue
            pts = np.empty((n_atoms, 3), dtype=np.float64)
            types = np.empty(n_atoms, dtype=np.int32)
            ids = np.empty(n_atoms, dtype=np.int64)
            for i in range(n_atoms):
                row = fh.readline().split()
                ids[i] = int(row[colmap["id"]])
                types[i] = int(row[colmap["type"]])
                pts[i, 0] = float(row[colmap["x"]])
                pts[i, 1] = float(row[colmap["y"]])
                pts[i, 2] = float(row[colmap["z"]])
            order = np.argsort(ids, kind="mergesort")
            pts = pts[order]
            types = types[order]
            H, origin, tric = parse_box(header, bounds)
            return Frame(
                positions=pts,
                types=types,
                H=H,
                origin=origin,
                is_triclinic=tric,
                source=f"{path}@{target_frame}",
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


def make_skew_narrow_frame(seed: int) -> Frame:
    nx = SKEW_FRAME_GRID_NX
    ny = SKEW_FRAME_GRID_NY
    nz = SKEW_FRAME_GRID_NZ
    lx = SKEW_FRAME_LX_ANGSTROM
    ly = SKEW_FRAME_LY_ANGSTROM
    lz = SKEW_FRAME_LZ_ANGSTROM
    xy = SKEW_FRAME_XY_ANGSTROM
    xz = SKEW_FRAME_XZ_ANGSTROM
    yz = SKEW_FRAME_YZ_ANGSTROM
    rng = np.random.default_rng(seed)
    pts = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                f = np.array(
                    [(i + 0.5) / nx, (j + 0.5) / ny, (k + 0.5) / nz],
                    dtype=np.float64,
                )
                f += (rng.random(3) - 0.5) * np.array(
                    [
                        SKEW_FRAME_JITTER_FRACTION / nx,
                        SKEW_FRAME_JITTER_FRACTION / ny,
                        SKEW_FRAME_JITTER_FRACTION / nz,
                    ]
                )
                x = f[0] * lx + f[1] * xy + f[2] * xz
                y = f[1] * ly + f[2] * yz
                z = f[2] * lz
                pts.append((x, y, z))
    H = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=np.float64)
    a = H[:, 0]
    b = H[:, 1]
    c = H[:, 2]
    volume = abs(float(np.linalg.det(H)))
    widths = np.array(
        [
            volume / np.linalg.norm(np.cross(b, c)),
            volume / np.linalg.norm(np.cross(c, a)),
            volume / np.linalg.norm(np.cross(a, b)),
        ],
        dtype=np.float64,
    )
    if not float(np.min(widths)) < SKEW_FRAME_REFERENCE_CUTOFF_ANGSTROM:
        raise AssertionError(
            "skew frame no longer has perpendicular width below cutoff reference"
        )
    return Frame(
        positions=np.array(pts, dtype=np.float64),
        types=np.full((nx * ny * nz,), 1, dtype=np.int32),
        H=H,
        origin=np.array([0.0, 0.0, 0.0], dtype=np.float64),
        is_triclinic=True,
        source=f"generated-skew-narrow-seed-{seed}",
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
    types = np.tile(frame.types, nx * ny * nz)
    H = frame.H.copy()
    H[:, 0] *= nx
    H[:, 1] *= ny
    H[:, 2] *= nz
    return Frame(
        positions=pts,
        types=types,
        H=H,
        origin=frame.origin.copy(),
        is_triclinic=frame.is_triclinic,
        source=f"{frame.source} rep={nx}x{ny}x{nz}",
    )


def save_frame_npz(frame: Frame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        positions=frame.positions,
        types=frame.types,
        H=frame.H,
        origin=frame.origin,
        triclinic=np.array([1 if frame.is_triclinic else 0], dtype=np.int8),
        source=np.array([frame.source], dtype=object),
    )


def load_frame_npz(path: str) -> Frame:
    with np.load(path, allow_pickle=True) as data:
        return Frame(
            positions=data["positions"],
            types=data["types"],
            H=data["H"],
            origin=data["origin"],
            is_triclinic=bool(int(data["triclinic"][0])),
            source=str(data["source"][0]),
        )


def parse_dseams_bench(text: str) -> dict[str, object]:
    out: dict[str, object] = {}
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    for line in lines:
        if line.startswith("meta "):
            out["meta"] = line
        elif line.startswith("cutoff_first_ms "):
            out["cutoff_first_ms"] = float(line.split()[1])
        elif line.startswith("cutoff_samples_ms "):
            body = line.split(" ", 1)[1].replace("cutoff_samples_ms ", "", 1)
            out["cutoff_samples_ms"] = [float(x) for x in body.split(",") if x]
        elif line.startswith("knn_first_ms "):
            out["knn_first_ms"] = float(line.split()[1])
        elif line.startswith("knn_samples_ms "):
            body = line.split(" ", 1)[1].replace("knn_samples_ms ", "", 1)
            out["knn_samples_ms"] = [float(x) for x in body.split(",") if x]
        elif line.startswith("cutoff_directed_edges "):
            out["cutoff_directed_edges"] = int(line.split()[1])
        elif line.startswith("knn_directed_edges "):
            out["knn_directed_edges"] = int(line.split()[1])
        elif line.startswith("knn_path_first "):
            out["knn_path_first"] = line.split()[1]
        elif line.startswith("knn_path_samples "):
            body = line.split(" ", 1)[1].replace("knn_path_samples ", "", 1)
            out["knn_path_samples"] = [x for x in body.split(",") if x]
    required = [
        "cutoff_first_ms",
        "cutoff_samples_ms",
        "knn_first_ms",
        "knn_samples_ms",
        "cutoff_directed_edges",
        "knn_directed_edges",
        "knn_path_first",
        "knn_path_samples",
    ]
    for key in required:
        if key not in out:
            raise RuntimeError(f"missing key {key} in d-SEAMS bench output\n{text}")
    return out


def run_dseams_once(
    bench_bin: str,
    spec: DatasetSpec,
    cutoff: float,
    k: int,
    cand: float,
    threads: int,
    git_sha_value: str,
    dirty_value: str,
    env_extra: dict[str, str] | None = None,
    dump_cutoff_edges: str | None = None,
) -> dict[str, object]:
    cmd = [
        bench_bin,
        "--traj",
        spec.path,
        "--frame",
        str(spec.frame),
        "--type",
        str(spec.atom_type),
        "--reps",
        "1",
        "--warmup",
        "0",
        "--cutoff",
        f"{cutoff:.8f}",
        "--k",
        str(k),
        "--k-cand",
        f"{cand:.8f}",
        "--nx",
        str(spec.rep_xyz[0]),
        "--ny",
        str(spec.rep_xyz[1]),
        "--nz",
        str(spec.rep_xyz[2]),
        "--threads",
        str(threads),
    ]
    if dump_cutoff_edges is not None:
        cmd.extend(["--dump-cutoff-edges", dump_cutoff_edges])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["SEAMS_GIT_SHA"] = git_sha_value
    env["SEAMS_GIT_DIRTY"] = dirty_value
    if env_extra:
        env.update(env_extra)
    pinned_cmd = pinned(cmd)
    text = subprocess.check_output(pinned_cmd, text=True, env=env)
    parsed = parse_dseams_bench(text)
    extra_env = ""
    if env_extra:
        extra_env = " ".join(f"{k}={v}" for k, v in sorted(env_extra.items())) + " "
    parsed["command"] = (
        f"OMP_NUM_THREADS={threads} OMP_DYNAMIC=FALSE "
        f"SEAMS_GIT_SHA={git_sha_value} SEAMS_GIT_DIRTY={dirty_value} "
        + extra_env
        + display_cmd(pinned_cmd)
    )
    parsed["raw_output"] = text
    return parsed


def triclinic_ckdtree_directed_edges(frame: Frame, cutoff: float) -> set[tuple[int, int]]:
    a = frame.H[:, 0]
    b = frame.H[:, 1]
    c = frame.H[:, 2]
    volume = abs(float(np.linalg.det(frame.H)))
    bxc = np.linalg.norm(np.cross(b, c))
    cxa = np.linalg.norm(np.cross(c, a))
    axb = np.linalg.norm(np.cross(a, b))
    widths = np.array([volume / bxc, volume / cxa, volume / axb], dtype=np.float64)
    reach = int(math.ceil(cutoff / float(np.min(widths)))) + 1

    image_points = []
    image_owner = []
    image_zero = []
    for na in range(-reach, reach + 1):
        for nb in range(-reach, reach + 1):
            for nc in range(-reach, reach + 1):
                shift = na * a + nb * b + nc * c
                zero = na == 0 and nb == 0 and nc == 0
                shifted = frame.positions + shift
                image_points.append(shifted)
                image_owner.append(np.arange(frame.positions.shape[0], dtype=np.int32))
                image_zero.append(np.full(frame.positions.shape[0], zero, dtype=bool))

    all_points = np.concatenate(image_points, axis=0)
    all_owner = np.concatenate(image_owner, axis=0)
    all_zero = np.concatenate(image_zero, axis=0)

    from scipy.spatial import cKDTree

    tree = cKDTree(all_points)
    directed: set[tuple[int, int]] = set()
    cutoff_sq = cutoff * cutoff
    for i, p in enumerate(frame.positions):
        cand = tree.query_ball_point(p, cutoff)
        if not cand:
            continue
        ids = np.array(cand, dtype=np.int64)
        owners = all_owner[ids]
        zeros = all_zero[ids]
        delta = all_points[ids] - p
        d2 = np.einsum("ij,ij->i", delta, delta, dtype=np.float64)
        keep = d2 <= cutoff_sq
        for owner, zero, ok in zip(owners, zeros, keep):
            if not ok:
                continue
            j = int(owner)
            if i == j and bool(zero):
                continue
            directed.add((i, j))
    return directed


def read_directed_edges(path: str) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            i_str, j_str = line.split()
            edges.add((int(i_str), int(j_str)))
    return edges


def ovito_directed_and_unique_edges(
    frame: Frame, cutoff: float, ovito_threads: int
) -> dict[str, int]:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    os.environ["OVITO_THREAD_COUNT"] = str(ovito_threads)
    from ovito.data import CutoffNeighborFinder, DataCollection

    data = DataCollection()
    data.create_cell(frame.H, pbc=(True, True, True))
    parts = data.create_particles(count=frame.positions.shape[0])
    parts.create_property("Position", data=frame.positions - frame.origin)
    pair_index = CutoffNeighborFinder(cutoff, data).find_all()[0]
    directed_edges = int(pair_index.shape[0])
    unique_directed_edges = int(np.unique(pair_index, axis=0).shape[0])
    return {
        "directed_edges": directed_edges,
        "unique_directed_edges": unique_directed_edges,
        "duplicate_image_edges": directed_edges - unique_directed_edges,
    }


def fit_exponent(series: list[dict[str, object]]) -> dict[str, object]:
    n = np.array([float(x["n_atoms"]) for x in series], dtype=np.float64)
    t = np.array([float(x["median_ms"]) for x in series], dtype=np.float64)
    slope, intercept = np.polyfit(np.log(n), np.log(t), deg=1)
    fit = np.exp(intercept) * n ** slope
    ss_res = float(np.sum((t - fit) ** 2))
    ss_tot = float(np.sum((t - np.mean(t)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    rmse = float(np.sqrt(np.mean((t - fit) ** 2)))
    residuals: list[dict[str, object]] = []
    for row, fitted in zip(series, fit):
        observed = float(row["median_ms"])
        residual_ms = observed - float(fitted)
        residuals.append(
            {
                "dataset": str(row.get("dataset", "unknown")),
                "n_atoms": int(row["n_atoms"]),
                "observed_median_ms": observed,
                "fitted_median_ms": float(fitted),
                "residual_ms": residual_ms,
                "residual_percent": (
                    100.0 * residual_ms / observed if observed != 0.0 else 0.0
                ),
            }
        )
    return {
        "slope_alpha": float(slope),
        "intercept_logC": float(intercept),
        "r2": float(r2),
        "rmse_ms": rmse,
        "residuals": residuals,
    }


def project_from_fit(fit: dict[str, float], n_atoms: int) -> float:
    return float(math.exp(fit["intercept_logC"]) * (float(n_atoms) ** fit["slope_alpha"]))


def run_external_worker(
    script_path: str,
    frame_npz: str,
    tool: str,
    cutoff: float,
    k: int,
    freud_threads: int,
    ovito_threads: int,
    scipy_threads: int,
) -> dict[str, object]:
    cmd = [
        sys.executable,
        script_path,
        "--worker-external",
        "--frame-npz",
        frame_npz,
        "--tool",
        tool,
        "--cutoff",
        f"{cutoff:.8f}",
        "--k",
        str(k),
        "--freud-threads",
        str(freud_threads),
        "--ovito-threads",
        str(ovito_threads),
        "--scipy-threads",
        str(scipy_threads),
    ]
    pinned_cmd = pinned(cmd)
    text = subprocess.check_output(pinned_cmd, text=True)
    out = json.loads(text)
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
    unsupported: dict[str, str] = {}

    for tool in tools:
        res = run_tool(tool)
        if "unsupported" in res:
            unsupported[tool] = str(res["unsupported"])
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
        if "sample_path" in res:
            records[tool]["path_samples"] = []
            records[tool]["path_first"] = str(res["sample_path"])

    if not active:
        return {"unsupported": unsupported}

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
            if "sample_path" in res:
                records[tool]["path_samples"].append(str(res["sample_path"]))

    for tool in active:
        samples = records[tool]["samples_ms"]
        records[tool].update(summarize_samples(samples))
        if "path_samples" in records[tool]:
            counts: dict[str, int] = {}
            for name in records[tool]["path_samples"]:
                counts[name] = counts.get(name, 0) + 1
            records[tool]["path_counts"] = counts

    if unsupported:
        records["unsupported"] = unsupported
    return records


def worker_external(args: argparse.Namespace) -> None:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    frame = load_frame_npz(args.frame_npz)
    tool = args.tool

    if tool.startswith("freud"):
        import freud

        freud.parallel.set_num_threads(args.freud_threads)
        box = freud.box.Box.from_matrix(frame.H)
        points = box.wrap(frame.positions - frame.origin)

        if tool == "freud_aabb_cutoff":
            try:
                t0 = time.perf_counter()
                count = int(
                    len(
                        freud.locality.AABBQuery(box, points)
                        .query(points, {"r_max": args.cutoff, "exclude_ii": True})
                        .toNeighborList()
                    )
                )
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            except RuntimeError as exc:
                print(json.dumps({"unsupported": str(exc), "tool": tool}))
            return

        if tool == "freud_linkcell_cutoff":
            try:
                t0 = time.perf_counter()
                count = int(
                    len(
                        freud.locality.LinkCell(
                            box, points, cell_width=args.cutoff
                        )
                        .query(points, {"r_max": args.cutoff, "exclude_ii": True})
                        .toNeighborList()
                    )
                )
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            except RuntimeError as exc:
                print(json.dumps({"unsupported": str(exc), "tool": tool}))
            return

        if tool == "freud_aabb_knn":
            try:
                t0 = time.perf_counter()
                count = int(
                    len(
                        freud.locality.AABBQuery(box, points)
                        .query(
                            points,
                            {"mode": "nearest", "num_neighbors": args.k, "exclude_ii": True},
                        )
                        .toNeighborList()
                    )
                )
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            except RuntimeError as exc:
                print(json.dumps({"unsupported": str(exc), "tool": tool}))
            return

        if tool == "freud_linkcell_knn":
            try:
                t0 = time.perf_counter()
                count = int(
                    len(
                        freud.locality.LinkCell(
                            box, points, cell_width=args.cutoff
                        )
                        .query(
                            points,
                            {
                                "mode": "nearest",
                                "num_neighbors": args.k,
                                "exclude_ii": True,
                            },
                        )
                        .toNeighborList()
                    )
                )
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            except RuntimeError as exc:
                print(json.dumps({"unsupported": str(exc), "tool": tool}))
            return

    if tool.startswith("ovito"):
        os.environ["OVITO_THREAD_COUNT"] = str(args.ovito_threads)
        from ovito.data import CutoffNeighborFinder, DataCollection, NearestNeighborFinder

        data = DataCollection()
        data.create_cell(frame.H, pbc=(True, True, True))
        parts = data.create_particles(count=frame.positions.shape[0])
        parts.create_property("Position", data=frame.positions - frame.origin)

        if tool == "ovito_cutoff":
            t0 = time.perf_counter()
            count = int(CutoffNeighborFinder(args.cutoff, data).find_all()[0].shape[0])
            t1 = time.perf_counter()
            print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            return

        if tool == "ovito_knn":
            t0 = time.perf_counter()
            count = int(NearestNeighborFinder(args.k, data).find_all()[0].shape[0])
            t1 = time.perf_counter()
            print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
            return

    if tool.startswith("scipy"):
        from scipy.spatial import cKDTree

        if tool == "scipy_ckdtree_knn" and frame.is_triclinic:
            print(json.dumps({"unsupported": "triclinic exact kNN with scipy cKDTree boxsize", "tool": tool}))
            return
        if tool == "scipy_ckdtree_cutoff" and frame.is_triclinic:
            print(
                json.dumps(
                    {
                        "unsupported": (
                            "triclinic periodic cutoff unsupported in scipy cKDTree boxsize; "
                            "fractional workaround undercounts tri-skew-512"
                        ),
                        "tool": tool,
                    }
                )
            )
            return

        if not frame.is_triclinic:
            box_lengths = np.array([frame.H[0, 0], frame.H[1, 1], frame.H[2, 2]], dtype=np.float64)
            tree = cKDTree(frame.positions, boxsize=box_lengths)
            if tool == "scipy_ckdtree_cutoff":
                t0 = time.perf_counter()
                count = int(2 * tree.query_pairs(args.cutoff, output_type="ndarray").shape[0])
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
                return
            if tool == "scipy_ckdtree_knn":
                t0 = time.perf_counter()
                dists = tree.query(
                    frame.positions, k=args.k + 1, workers=args.scipy_threads
                )[0][:, 1:]
                count = int(np.isfinite(dists).sum())
                t1 = time.perf_counter()
                print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "tool": tool}))
                return

    print(json.dumps({"unsupported": f"unknown tool {tool}", "tool": tool}))


def main() -> None:
    self_check_summary_rules()

    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-external", action="store_true")
    ap.add_argument("--frame-npz")
    ap.add_argument("--tool")
    ap.add_argument("--bench-bin")
    ap.add_argument("--output-json")
    ap.add_argument("--cutoff", type=float, default=3.5)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument(
        "--k-candidate",
        type=float,
        default=MW_K_CANDIDATE_ANGSTROM,
    )
    ap.add_argument("--confidence", type=float, default=0.95)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--dseams-threads", type=int, default=1)
    ap.add_argument("--freud-threads", type=int, default=1)
    ap.add_argument("--ovito-threads", type=int, default=1)
    ap.add_argument("--scipy-threads", type=int, default=1)
    ap.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Restrict to one or more dataset names (repeatable).",
    )
    ap.add_argument(
        "--skip-fallback-scaling",
        action="store_true",
        help="Skip forced fallback scaling series (debug/smoke only).",
    )
    args = ap.parse_args()

    if args.worker_external:
        worker_external(args)
        return

    if not args.bench_bin or not args.output_json:
        raise SystemExit("--bench-bin and --output-json are required")

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

    skew_path = "input/traj/skew_triclinic_width_lt_cutoff.lammpstrj"
    write_lammps_dump(skew_path, make_skew_narrow_frame(args.seed), atom_type=1)

    datasets = [
        DatasetSpec("orth-4k", "input/traj/mW_cubic.lammpstrj", 1, 1, (1, 1, 1)),
        DatasetSpec("orth-32k", "input/traj/mW_cubic.lammpstrj", 1, 1, (2, 2, 2)),
        DatasetSpec("orth-1e5", "input/traj/mW_cubic.lammpstrj", 1, 1, (3, 3, 3)),
        DatasetSpec("orth-1e6", "input/traj/mW_cubic.lammpstrj", 1, 1, (6, 6, 7)),
        DatasetSpec(
            "tri-136",
            "input/traj/genice_sH.lammpstrj",
            1,
            1,
            (1, 1, 1),
            check_cutoff_set_equality=True,
        ),
        DatasetSpec("tri-1088", "input/traj/genice_sH.lammpstrj", 1, 1, (2, 2, 2)),
        DatasetSpec("tri-17k", "input/traj/genice_sH.lammpstrj", 1, 1, (5, 5, 5)),
        DatasetSpec("tri-46k", "input/traj/genice_sH.lammpstrj", 1, 1, (7, 7, 7)),
        DatasetSpec(
            "tri-skew-512",
            skew_path,
            1,
            1,
            (1, 1, 1),
            check_cutoff_set_equality=True,
        ),
        DatasetSpec("tri-1e5", "input/traj/genice_sH.lammpstrj", 1, 1, (10, 10, 10)),
        DatasetSpec("tri-1e6", "input/traj/genice_sH.lammpstrj", 1, 1, (20, 20, 20)),
    ]

    selected = set(args.dataset)
    out: dict[str, object] = {
        "git_sha": sha,
        "git_dirty": dirty,
        "seed": args.seed,
        "cutoff": args.cutoff,
        "k": args.k,
        "k_candidate": args.k_candidate,
        "k_candidate_reference": {
            "name": "mW radial distribution second-shell window midpoint",
            "source_file": repo_relative(MW_K_CANDIDATE_DERIVATION.source_path),
            "source_script": "tests/derive_mw_gr_window_midpoint.py",
            "source_frame": MW_K_CANDIDATE_DERIVATION.frame,
            "source_atom_type": MW_K_CANDIDATE_DERIVATION.atom_type,
            "method": MW_K_CANDIDATE_DERIVATION.method,
            "window_left_angstrom": (
                MW_K_CANDIDATE_DERIVATION.shell_window_left_angstrom
            ),
            "window_right_angstrom": (
                MW_K_CANDIDATE_DERIVATION.shell_window_right_angstrom
            ),
            "window_midpoint_angstrom": MW_K_CANDIDATE_DERIVATION.candidate_angstrom,
        },
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
            "freud_threads": args.freud_threads,
            "ovito_threads": args.ovito_threads,
            "scipy_threads": args.scipy_threads,
        },
        "linkcell_policy": {
            "cell_width": "r_max",
            "projection_trigger_rule": (
                "project when n_atoms exceeds the largest fitted n_atoms for each "
                "geometry/stage pair"
            ),
            "fit_selection_rule": LINKCELL_FIT_SELECTION_RULE,
            "minimum_fit_points": LINKCELL_FIT_POINTS_PER_GEOMETRY_STAGE,
            "fits": {},
            "projections": {},
        },
        "datasets": {},
        "fallback_scaling": {},
    }

    frame_store = Path("input/.bench_frames")
    frame_store.mkdir(parents=True, exist_ok=True)
    linkcell_fit_points: dict[tuple[str, str], list[dict[str, object]]] = {}
    linkcell_fit_models: dict[tuple[str, str], dict[str, object]] = {}

    for spec in datasets:
        if selected and spec.name not in selected:
            continue
        base = read_lammps_frame(spec.path, spec.frame)
        frame = replicate_frame(base, *spec.rep_xyz)
        n_atoms = int(frame.positions.shape[0])

        frame_npz = str(frame_store / f"{spec.name}.npz")
        save_frame_npz(frame, frame_npz)

        print(f"[dataset] {spec.name} n={n_atoms} triclinic={frame.is_triclinic}", flush=True)

        geometry = "triclinic" if frame.is_triclinic else "orthorhombic"

        def projected_linkcell_skip(stage: str) -> dict[str, object] | None:
            fit_key = (stage, geometry)
            fit = linkcell_fit_models.get(fit_key)
            if fit is None:
                return None
            trigger_n_atoms = int(fit["projection_trigger_n_atoms"])
            if n_atoms <= trigger_n_atoms:
                return None
            projected_sample_ms = project_from_fit(fit, n_atoms)
            sample_calls = 1 + warmup + repeats
            projected_campaign_hours = projected_sample_ms * sample_calls / 3_600_000.0
            projection = {
                "fit_key": f"{geometry}/{stage}",
                "n_atoms": n_atoms,
                "projection_trigger_n_atoms": trigger_n_atoms,
                "projected_sample_ms": projected_sample_ms,
                "projected_calls": sample_calls,
                "projected_campaign_hours": projected_campaign_hours,
                "fit_datasets": fit["datasets"],
                "fit": {
                    "slope_alpha": fit["slope_alpha"],
                    "intercept_logC": fit["intercept_logC"],
                    "r2": fit["r2"],
                    "rmse_ms": fit["rmse_ms"],
                },
            }
            out["linkcell_policy"]["projections"][f"{spec.name}/{stage}"] = projection
            return {
                "unsupported": (
                    "not run (projected "
                    f"{projected_sample_ms:.3f} ms/sample, "
                    f"{projected_campaign_hours:.3f} h for {sample_calls} calls)"
                )
            }

        def run_stage_tool(stage: str, tool: str) -> dict[str, object]:
            if tool.startswith("dseams"):
                threads = 1 if tool.endswith("omp1") else args.dseams_threads
                parsed = run_dseams_once(
                    args.bench_bin,
                    spec,
                    cutoff=args.cutoff,
                    k=args.k,
                    cand=args.k_candidate,
                    threads=threads,
                    git_sha_value=sha,
                    dirty_value=dirty,
                )
                key = "cutoff" if stage == "cutoff" else "knn"
                res: dict[str, object] = {
                    "first_ms": float(parsed[f"{key}_first_ms"]),
                    "sample_ms": float(parsed[f"{key}_samples_ms"][0]),
                    "count": int(parsed[f"{key}_directed_edges"]),
                    "command": str(parsed["command"]),
                }
                if stage == "knn":
                    path_samples = parsed.get("knn_path_samples", [])
                    res["sample_path"] = str(path_samples[0]) if path_samples else "unknown"
                return res

            if tool.startswith("freud_linkcell"):
                projected = projected_linkcell_skip(stage)
                if projected is not None:
                    return projected

            ext = run_external_worker(
                script_path=str(Path(__file__)),
                frame_npz=frame_npz,
                tool=tool,
                cutoff=args.cutoff,
                k=args.k,
                freud_threads=args.freud_threads,
                ovito_threads=args.ovito_threads,
                scipy_threads=args.scipy_threads,
            )
            if "unsupported" in ext:
                return {"unsupported": ext["unsupported"], "command": ext.get("command")}
            return {
                "first_ms": float(ext["ms"]),
                "sample_ms": float(ext["ms"]),
                "count": int(ext["count"]),
                "command": ext.get("command"),
            }

        cutoff_tools = [
            "dseams_omp1",
            "freud_aabb_cutoff",
            "freud_linkcell_cutoff",
            "ovito_cutoff",
            "scipy_ckdtree_cutoff",
        ]
        knn_tools = [
            "dseams_omp1",
            "freud_aabb_knn",
            "freud_linkcell_knn",
            "ovito_knn",
            "scipy_ckdtree_knn",
        ]

        cutoff_result = run_interleaved(
            cutoff_tools,
            warmup,
            repeats,
            lambda tool: run_stage_tool("cutoff", tool),
            f"{spec.name}/cutoff",
        )
        knn_result = run_interleaved(
            knn_tools,
            warmup,
            repeats,
            lambda tool: run_stage_tool("knn", tool),
            f"{spec.name}/knn",
        )

        dataset_out: dict[str, object] = {
            "source": frame.source,
            "replication": {"nx": spec.rep_xyz[0], "ny": spec.rep_xyz[1], "nz": spec.rep_xyz[2]},
            "n_atoms": n_atoms,
            "triclinic": bool(frame.is_triclinic),
            "count_semantics": {
                "cutoff_directed_edges": (
                    "directed neighbour-image records returned by each backend"
                ),
                "dseams_knn_directed_edges": (
                    "directed edges in the mutual-k graph; this can be below n_atoms * k"
                ),
                "freud_knn_directed_edges": (
                    "directed neighbours from fixed-k query; equals n_atoms * k when supported"
                ),
            },
            "cutoff_stage": cutoff_result,
            "knn_stage": knn_result,
        }

        expected_directed_knn = n_atoms * args.k
        dseams_knn = None
        freud_knn = None
        if "dseams_omp1" in knn_result:
            dseams_knn = int(knn_result["dseams_omp1"]["count"])
        if "freud_aabb_knn" in knn_result:
            freud_knn = int(knn_result["freud_aabb_knn"]["count"])
        if dseams_knn is not None or freud_knn is not None:
            dataset_out["knn_count_explanation"] = {
                "expected_directed_edges_n_times_k": expected_directed_knn,
                "dseams_mutual_graph_directed_edges": dseams_knn,
                "freud_fixed_k_directed_edges": freud_knn,
                "dseams_minus_expected": (
                    None if dseams_knn is None else dseams_knn - expected_directed_knn
                ),
                "freud_minus_expected": (
                    None if freud_knn is None else freud_knn - expected_directed_knn
                ),
            }

        if spec.check_cutoff_set_equality:
            dump_path = f"input/{spec.name}-dseams-cutoff-edges.txt"
            parity = run_dseams_once(
                args.bench_bin,
                spec,
                cutoff=args.cutoff,
                k=args.k,
                cand=args.k_candidate,
                threads=1,
                git_sha_value=sha,
                dirty_value=dirty,
                dump_cutoff_edges=dump_path,
            )
            dseams_edges = read_directed_edges(dump_path)
            scipy_edges = triclinic_ckdtree_directed_edges(frame, args.cutoff)
            only_dseams = len(dseams_edges - scipy_edges)
            only_scipy = len(scipy_edges - dseams_edges)
            equal = only_dseams == 0 and only_scipy == 0
            if not equal:
                raise AssertionError(
                    f"{spec.name}: cutoff-set mismatch (dseams_only={only_dseams}, scipy_only={only_scipy})"
                )
            dataset_out["ckdtree_cutoff_set_equality"] = {
                "asserted": True,
                "equal": True,
                "edge_count_directed": len(dseams_edges),
                "dseams_only_edges": 0,
                "scipy_only_edges": 0,
                "dseams_dump_command": parity["command"],
            }
            ovito_counts = ovito_directed_and_unique_edges(
                frame, args.cutoff, args.ovito_threads
            )
            dataset_out["ovito_cutoff_edge_accounting"] = ovito_counts

        for stage, tool_key in (
            ("cutoff", "freud_linkcell_cutoff"),
            ("knn", "freud_linkcell_knn"),
        ):
            stage_map = dataset_out[f"{stage}_stage"]
            if tool_key not in stage_map:
                continue
            stage_row = stage_map[tool_key]
            if "median_ms" not in stage_row:
                continue
            fit_key = (stage, geometry)
            points = linkcell_fit_points.setdefault(fit_key, [])
            points.append(
                {
                    "dataset": spec.name,
                    "n_atoms": n_atoms,
                    "median_ms": float(stage_row["median_ms"]),
                }
            )
            points.sort(key=lambda row: int(row["n_atoms"]), reverse=True)
            if len(points) > LINKCELL_FIT_POINTS_PER_GEOMETRY_STAGE:
                del points[LINKCELL_FIT_POINTS_PER_GEOMETRY_STAGE:]
            if len(points) >= LINKCELL_FIT_POINTS_PER_GEOMETRY_STAGE:
                fit_points = sorted(points, key=lambda row: int(row["n_atoms"]))
                fit = fit_exponent(fit_points)
                max_fit_n_atoms = max(int(p["n_atoms"]) for p in fit_points)
                fit_with_meta: dict[str, object] = {
                    "slope_alpha": fit["slope_alpha"],
                    "intercept_logC": fit["intercept_logC"],
                    "r2": fit["r2"],
                    "rmse_ms": fit["rmse_ms"],
                    "datasets": [str(p["dataset"]) for p in fit_points],
                    "projection_trigger_n_atoms": max_fit_n_atoms,
                    "residuals": fit["residuals"],
                }
                linkcell_fit_models[fit_key] = fit_with_meta
                out["linkcell_policy"]["fits"][f"{geometry}/{stage}"] = fit_with_meta

        out["datasets"][spec.name] = dataset_out

    # O(N^2) removal evidence: forced fallback cell-list vs brute across 4 sizes.
    if args.skip_fallback_scaling:
        out["fallback_scaling"] = {"skipped": True}
        out["fallback_scaling_fit"] = {"skipped": True}
    else:
        fallback_specs = [
            DatasetSpec("tri-fallback-136", "input/traj/genice_sH.lammpstrj", 1, 1, (1, 1, 1)),
            DatasetSpec("tri-fallback-1088", "input/traj/genice_sH.lammpstrj", 1, 1, (2, 2, 2)),
            DatasetSpec("tri-fallback-3672", "input/traj/genice_sH.lammpstrj", 1, 1, (3, 3, 3)),
            DatasetSpec("tri-fallback-8704", "input/traj/genice_sH.lammpstrj", 1, 1, (4, 4, 4)),
        ]

        fit_cell: list[dict[str, object]] = []
        fit_brute: list[dict[str, object]] = []
        for spec in fallback_specs:
            frame = replicate_frame(read_lammps_frame(spec.path, spec.frame), *spec.rep_xyz)
            n_atoms = int(frame.positions.shape[0])
            print(f"[fallback] {spec.name} n={n_atoms}", flush=True)

            def run_fallback(mode: str) -> dict[str, object]:
                parsed = run_dseams_once(
                    args.bench_bin,
                    spec,
                    cutoff=args.cutoff,
                    k=args.k,
                    cand=args.k_candidate,
                    threads=1,
                    git_sha_value=sha,
                    dirty_value=dirty,
                    env_extra={"SEAMS_KNN_FORCE_FALLBACK": "1", "SEAMS_KNN_FALLBACK": mode},
                )
                return {
                    "first_ms": float(parsed["knn_first_ms"]),
                    "sample_ms": float(parsed["knn_samples_ms"][0]),
                    "count": int(parsed["knn_directed_edges"]),
                    "sample_path": str(parsed.get("knn_path_samples", ["unknown"])[0]),
                    "command": str(parsed["command"]),
                }

            fallback_result = run_interleaved(
                ["cell", "brute"],
                warmup,
                repeats,
                lambda tool: run_fallback("cell" if tool == "cell" else "brute"),
                f"{spec.name}/fallback",
            )

            cell = fallback_result["cell"]
            brute = fallback_result["brute"]
            speedup = (
                float(brute["median_ms"]) / float(cell["median_ms"])
                if float(cell["median_ms"]) > 0.0
                else 0.0
            )

            out["fallback_scaling"][spec.name] = {
                "n_atoms": n_atoms,
                "cell": cell,
                "brute": brute,
                "speedup_brute_over_cell": speedup,
            }
            fit_cell.append({"n_atoms": n_atoms, "median_ms": float(cell["median_ms"])})
            fit_brute.append({"n_atoms": n_atoms, "median_ms": float(brute["median_ms"])})

        out["fallback_scaling_fit"] = {
            "cell": fit_exponent(fit_cell),
            "brute": fit_exponent(fit_brute),
        }

    Path(args.output_json).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {args.output_json}")


if __name__ == "__main__":
    main()
