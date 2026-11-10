"""Write input/traj/ice1h.lammpstrj from GenIce ice Ih.

    python scripts/gen_ice1h.py [output]

Runs `genice2 ice1h --rep 2 2 2` (GenIce 2.2.13).
Oxygen rows are nanometres and are stored in angstroms.
The cell on stderr is the orthorhombic box, also in angstroms.
"""

import pathlib
import re
import subprocess
import sys


def parse(stdout, stderr):
    rows = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) < 6 or parts[1] != "O":
            continue
        rows.append(
            (float(parts[3]) * 10.0, float(parts[4]) * 10.0, float(parts[5]) * 10.0)
        )
    cell = {}
    # The last cell before the molecules are placed is in nanometres.
    for line in stderr.splitlines():
        if line.startswith("INFO:Stage2:"):
            break
        match = re.search(r"([abc]) = ([0-9.]+)", line)
        if match:
            cell[match.group(1)] = float(match.group(2))
    if set(cell) != {"a", "b", "c"} or not rows:
        raise SystemExit("genice2 did not return a cell and oxygen rows")
    if cell["a"] < 5.0:
        cell = {key: value * 10.0 for key, value in cell.items()}
    return (cell["a"], cell["b"], cell["c"]), rows


def dump(box, rows):
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(rows)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0.0 {box[0]:.10f}",
        f"0.0 {box[1]:.10f}",
        f"0.0 {box[2]:.10f}",
        "ITEM: ATOMS id type x y z",
    ]
    for index, (x, y, z) in enumerate(rows, start=1):
        lines.append(f"{index} 1 {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines) + "\n"


def main():
    root = pathlib.Path(__file__).resolve().parents[1]
    dest = (
        pathlib.Path(sys.argv[1])
        if len(sys.argv) > 1
        else root / "input" / "traj" / "ice1h.lammpstrj"
    )
    run = subprocess.run(
        ["genice2", "ice1h", "--rep", "2", "2", "2"],
        check=True,
        text=True,
        capture_output=True,
    )
    box, rows = parse(run.stdout, run.stderr)
    dest.write_text(dump(box, rows))
    print(f"{len(rows)} oxygens box {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}")


if __name__ == "__main__":
    main()
