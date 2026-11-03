#!/usr/bin/env python3
"""Derive the mW k-candidate from first-frame neighbour-shell geometry."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).resolve().parents[1]
MW_SOURCE_TRAJECTORY = REPO_ROOT / "input" / "traj" / "mW_cubic.lammpstrj"
SHELL_SCAN_MAX_RADIUS_ANGSTROM = 7.0
SHELL_DISTANCE_DECIMALS = 3
FIRST_SHELL_CUTOFF_ANGSTROM = 3.5
MIN_GAP_FOR_SHELL_BREAK_ANGSTROM = 0.1
WINDOW_MIDPOINT_DECIMALS = 3


@dataclass(frozen=True)
class KCandidateDerivation:
    candidate_angstrom: float
    shell_window_left_angstrom: float
    shell_window_right_angstrom: float
    source_path: Path
    frame: int
    atom_type: int
    method: str


def _read_orthorhombic_frame_positions(
    traj_path: Path, frame: int, atom_type: int
) -> tuple[np.ndarray, np.ndarray]:
    if frame < 1:
        raise ValueError(f"frame index must be >= 1, got {frame}")
    frame_index = 0
    with traj_path.open("r", encoding="utf-8") as fh:
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"frame {frame} not found in {traj_path}")
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            _ = int(fh.readline().strip())
            frame_index += 1
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"malformed dump in {traj_path}")
            n_atoms = int(fh.readline().strip())
            header = fh.readline().strip()
            if not header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"missing box bounds in {traj_path}")
            if "xy" in header.lower():
                raise ValueError(
                    f"{traj_path} frame {frame} is triclinic; expected orthorhombic mW source"
                )
            xlo, xhi = map(float, fh.readline().split()[:2])
            ylo, yhi = map(float, fh.readline().split()[:2])
            zlo, zhi = map(float, fh.readline().split()[:2])
            atoms_header = fh.readline().strip()
            if not atoms_header.startswith("ITEM: ATOMS"):
                raise ValueError(f"missing atom header in {traj_path}")
            cols = atoms_header.split()[2:]
            colmap = {name: idx for idx, name in enumerate(cols)}
            required = {"id", "type", "x", "y", "z"}
            if not required.issubset(colmap):
                raise ValueError(f"missing columns {required} in {traj_path}")

            rows = []
            for _ in range(n_atoms):
                tok = fh.readline().split()
                if int(tok[colmap["type"]]) != atom_type:
                    continue
                rows.append(
                    (
                        int(tok[colmap["id"]]),
                        float(tok[colmap["x"]]),
                        float(tok[colmap["y"]]),
                        float(tok[colmap["z"]]),
                    )
                )
            if frame_index != frame:
                continue
            rows.sort(key=lambda r: r[0])
            pts = np.array([(x, y, z) for _, x, y, z in rows], dtype=np.float64)
            box_lengths = np.array([xhi - xlo, yhi - ylo, zhi - zlo], dtype=np.float64)
            return pts, box_lengths


def _second_shell_window_from_pairs(
    positions: np.ndarray, box_lengths: np.ndarray
) -> tuple[float, float]:
    tree = cKDTree(positions, boxsize=box_lengths)
    pairs = tree.query_pairs(SHELL_SCAN_MAX_RADIUS_ANGSTROM, output_type="ndarray")
    if pairs.size == 0:
        raise ValueError("no pairs in requested shell-scan radius")
    delta = positions[pairs[:, 0]] - positions[pairs[:, 1]]
    delta -= box_lengths * np.round(delta / box_lengths)
    dist = np.linalg.norm(delta, axis=1)
    shell_distances = np.unique(np.round(dist, SHELL_DISTANCE_DECIMALS))
    shell_distances = shell_distances[shell_distances > FIRST_SHELL_CUTOFF_ANGSTROM]
    if shell_distances.size < 3:
        raise ValueError("insufficient shell distances for second-shell window derivation")
    gaps = np.diff(shell_distances)
    shell_break_idx = np.where(gaps >= MIN_GAP_FOR_SHELL_BREAK_ANGSTROM)[0]
    if shell_break_idx.size < 2:
        raise ValueError(
            "need at least two shell breaks; adjust scan radius or shell-gap threshold"
        )
    top_breaks = sorted(shell_break_idx, key=lambda i: gaps[i], reverse=True)[:2]
    mids = sorted(
        0.5 * (shell_distances[i] + shell_distances[i + 1]) for i in top_breaks
    )
    return float(mids[0]), float(mids[1])


def derive_k_candidate_from_data(
    source_path: Path = MW_SOURCE_TRAJECTORY, frame: int = 1, atom_type: int = 1
) -> KCandidateDerivation:
    pts, box_lengths = _read_orthorhombic_frame_positions(source_path, frame, atom_type)
    shell_left, shell_right = _second_shell_window_from_pairs(pts, box_lengths)
    midpoint_quant = Decimal(f"1e-{WINDOW_MIDPOINT_DECIMALS}")
    candidate = float(
        Decimal(str(0.5 * (shell_left + shell_right))).quantize(
            midpoint_quant, rounding=ROUND_HALF_UP
        )
    )
    return KCandidateDerivation(
        candidate_angstrom=candidate,
        shell_window_left_angstrom=shell_left,
        shell_window_right_angstrom=shell_right,
        source_path=source_path,
        frame=frame,
        atom_type=atom_type,
        method=(
            "midpoint between the two largest pair-distance shell breaks above "
            f"{FIRST_SHELL_CUTOFF_ANGSTROM} A (scan radius "
            f"{SHELL_SCAN_MAX_RADIUS_ANGSTROM} A)"
        ),
    )


MW_K_CANDIDATE_DERIVATION = derive_k_candidate_from_data()
MW_K_CANDIDATE_ANGSTROM = MW_K_CANDIDATE_DERIVATION.candidate_angstrom
