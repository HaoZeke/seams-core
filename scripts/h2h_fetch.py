#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Download the head-to-head frames from d-SEAMS/h2h-systems.

Each file is an xz-compressed LAMMPS dump. huggingface_hub fetches it
without a token. A missing package or a failed download falls back to
scripts/h2h_cpc.py prepare, which rebuilds the frame from the ice
lattice and input/traj/dump-l.lammpstrj.
"""

from __future__ import annotations

import lzma
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = ROOT / "h2h-cpc" / "systems"
REPO = "d-SEAMS/h2h-systems"
NAMES = (
    "ice-1k",
    "ice-10k",
    "ice-50k",
    "ice-110k",
    "liq-1k",
    "liq-10k",
    "liq-50k",
    "liq-110k",
    "interface-50k",
)


def download(name: str, dest: Path) -> None:
    from huggingface_hub import hf_hub_download

    cached = hf_hub_download(
        repo_id=REPO,
        filename=f"{name}.lammpstrj.xz",
        repo_type="dataset",
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".lammpstrj.part")
    with lzma.open(cached, "rb") as src, tmp.open("wb") as out:
        shutil.copyfileobj(src, out)
    tmp.replace(dest)


def generate(name: str) -> None:
    subprocess.run(
        ["python3.12", "scripts/h2h_cpc.py", "prepare", name],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def fetch_one(name: str) -> Path:
    dest = SYSTEMS / f"{name}.lammpstrj"
    if dest.is_file():
        print(dest.relative_to(ROOT))
        return dest
    try:
        download(name, dest)
    except Exception as exc:
        print(f"{name} download failed ({exc.__class__.__name__}); generating", file=sys.stderr)
        generate(name)
    if not dest.is_file():
        raise SystemExit(f"missing {dest.name}")
    print(dest.relative_to(ROOT))
    return dest


def main() -> None:
    names = sys.argv[1:] or list(NAMES)
    for name in names:
        if name not in NAMES:
            raise SystemExit(f"unknown system {name}")
        fetch_one(name)


if __name__ == "__main__":
    main()
