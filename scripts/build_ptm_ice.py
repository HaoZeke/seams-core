#!/usr/bin/env python3
"""Build libptm_ice.so from the standalone PTM library.

Clones github.com/pmla/polyhedral-template-matching at a pinned commit,
applies scripts/ptm_index_rmsd.patch, and links scripts/ptm_ice.cpp.
The library and the clone stay in the build tree.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "bbdir-release" / "ptm-src"
SO = ROOT / "bbdir-release" / "libptm_ice.so"
PATCH = ROOT / "scripts" / "ptm_index_rmsd.patch"
WRAPPER = ROOT / "scripts" / "ptm_ice.cpp"
URL = "https://github.com/pmla/polyhedral-template-matching"
REV = "162cc83401bee802cebe987bdbf3db2a3c1b1e71"
UNITS = (
    "ptm_alloy_types.cpp",
    "ptm_canonical_coloured.cpp",
    "ptm_convex_hull_incremental.cpp",
    "ptm_deformation_gradient.cpp",
    "ptm_graph_data.cpp",
    "ptm_graph_tools.cpp",
    "ptm_index.cpp",
    "ptm_initialize_data.cpp",
    "ptm_multishell.cpp",
    "ptm_neighbour_ordering.cpp",
    "ptm_normalize_vertices.cpp",
    "ptm_polar.cpp",
    "ptm_quat.cpp",
    "ptm_structure_matcher.cpp",
    "ptm_voronoi_cell.cpp",
)


def rev_of() -> str:
    return subprocess.check_output(
        ["git", "-C", str(SRC), "rev-parse", "HEAD"], text=True
    ).strip()


def ensure_source() -> str:
    SRC.parent.mkdir(parents=True, exist_ok=True)
    if not (SRC / "ptm_functions.h").is_file():
        subprocess.check_call(
            ["git", "clone", "--depth", "1", URL, str(SRC)]
        )
    if not rev_of().startswith(REV[:12]):
        subprocess.check_call(["git", "-C", str(SRC), "fetch", "--depth", "1", "origin", REV])
        subprocess.check_call(["git", "-C", str(SRC), "checkout", REV])
    # git apply is a no-op once the hunk is present.
    check = subprocess.run(
        ["git", "-C", str(SRC), "apply", "--check", str(PATCH)],
        capture_output=True,
        text=True,
    )
    if check.returncode == 0:
        subprocess.check_call(["git", "-C", str(SRC), "apply", str(PATCH)])
    return rev_of()


def build() -> tuple[Path, str]:
    rev = ensure_source()
    newest = max(WRAPPER.stat().st_mtime, PATCH.stat().st_mtime)
    newest = max(newest, *( (SRC / name).stat().st_mtime for name in UNITS ))
    if SO.is_file() and SO.stat().st_mtime >= newest:
        return SO, rev
    cmd = [
        "g++", "-O3", "-fPIC", "-std=c++17", "-fopenmp",
        f"-I{SRC}", "-shared", str(WRAPPER),
        *[str(SRC / name) for name in UNITS],
        "-o", str(SO),
    ]
    subprocess.check_call(cmd)
    return SO, rev


def main() -> None:
    path, rev = build()
    print(f"{path} {rev[:12]}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode or 1) from exc
