"""Count published CHILL+ labels on a LAMMPS dump using chill_independent.

The reference does not import d-SEAMS. Oxygen rows are the chosen type.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chill_independent import classify


def load(path, atom_type):
    lines = Path(path).read_text().splitlines()
    cols = []
    pos = []
    bounds = []
    atoms = False
    frame = 0
    for line in lines:
        if line.startswith("ITEM: TIMESTEP"):
            frame += 1
            if frame > 1:
                break
            atoms = False
            continue
        if line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            atoms = True
            continue
        if atoms:
            parts = line.split()
            if len(parts) < len(cols):
                continue
            rec = dict(zip(cols, parts))
            if int(float(rec["type"])) != atom_type:
                continue
            pos.append([float(rec["x"]), float(rec["y"]), float(rec["z"])])
        elif not line.startswith("ITEM:") and len(bounds) < 3 and line.strip():
            parts = line.split()
            if len(parts) < 2:
                continue
            lo, hi = parts[:2]
            bounds.append(float(hi) - float(lo))
    return np.asarray(pos, dtype=float), np.asarray(bounds, dtype=float)


def main():
    path = sys.argv[1]
    atom_type = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    pos, box = load(path, atom_type)
    result = classify(pos, box, "shell", "nearest", True, "plus")
    labels = result["label"]
    names, counts = np.unique(labels, return_counts=True)
    print("n", len(labels))
    for name, count in zip(names, counts):
        print(name, int(count))


if __name__ == "__main__":
    main()
