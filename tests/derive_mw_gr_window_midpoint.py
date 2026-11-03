#!/usr/bin/env python3
"""Derive the mW k-candidate from the mW trajectory radial-distribution data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mw_gr_window import derive_k_candidate_from_data


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source-traj",
        default="input/traj/mW_cubic.lammpstrj",
        help="Path to the orthorhombic mW trajectory used for derivation.",
    )
    ap.add_argument("--frame", type=int, default=1)
    ap.add_argument("--atom-type", type=int, default=1)
    args = ap.parse_args()

    source = Path(args.source_traj)
    if not source.is_absolute():
        source = (Path(__file__).resolve().parents[1] / source).resolve()
    deriv = derive_k_candidate_from_data(
        source_path=source, frame=args.frame, atom_type=args.atom_type
    )
    out = {
        "source_traj": str(deriv.source_path),
        "frame": deriv.frame,
        "atom_type": deriv.atom_type,
        "candidate_angstrom": deriv.candidate_angstrom,
        "shell_window_left_angstrom": deriv.shell_window_left_angstrom,
        "shell_window_right_angstrom": deriv.shell_window_right_angstrom,
        "method": deriv.method,
    }
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
