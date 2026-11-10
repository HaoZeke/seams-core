"""Write per-atom labels and print every agreement figure.

Environment:
  FRAME_DIR LIQUID OUT
  BASE_BIN TIP_BIN BASE_LIB TIP_LIB TIP_LANG
The liquid dump has no type column.
The reader stores an atom only when a type column is present.
This run inserts type 1 and keeps every coordinate.
"""
import ctypes
import os
import subprocess
from pathlib import Path

import numpy as np

FRAME_DIR = Path(os.environ["FRAME_DIR"])
LIQUID = Path(os.environ["LIQUID"])
OUT = Path(os.environ["OUT"])
OUT.mkdir(parents=True, exist_ok=True)


def write_frame(path, pos, box):
    n = pos.shape[0]
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n" % n)
        handle.write("ITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write("0 %s\n" % repr(float(length)))
        handle.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos, start=1):
            handle.write("%d 1 %.10f %.10f %.10f\n" % (i, x, y, z))


def run_bin(binary, lib, mode, path):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    proc = subprocess.run(
        [binary, mode, str(path), "1", "1"],
        check=False, capture_output=True, text=True, env=env,
    )
    if proc.returncode != 0:
        raise RuntimeError(mode + " failed on " + path.name)
    return [int(line) for line in proc.stdout.splitlines() if line.strip()]


def save(path, values):
    path.write_text("".join("%d\n" % v for v in values))


def ovito_of(path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from ovito.io import import_file
    from ovito.modifiers import ChillPlusModifier

    pipeline = import_file(str(path))
    pipeline.modifiers.append(ChillPlusModifier(cutoff=3.5))
    data = pipeline.compute()
    structure = np.array(data.particles["Structure Type"])
    mapping = {0: 2, 1: 1, 2: 0, 3: 3, 4: 4, 5: 5}
    return np.array([mapping[int(v)] for v in structure], dtype=np.int32)


def library_rule(pos, box, rule):
    lib = ctypes.CDLL(os.environ["TIP_LANG"])
    lib.seams_chill_plus_rule.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_double,
        ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus_rule.restype = ctypes.c_int
    n = pos.shape[0]
    xyz = np.ascontiguousarray(pos, dtype=np.float64).reshape(-1)
    cell = np.array([0, 0, 0, box[0], box[1], box[2]], dtype=np.float64)
    labels = np.empty(n, dtype=np.int32)
    rc = lib.seams_chill_plus_rule(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        cell.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 3.5,
        b"cutoff", rule.encode(),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if rc != 0:
        raise RuntimeError("library " + rule + " returned " + str(rc))
    return labels


def frames():
    jobs = []
    for name in ["Ih", "Ic", "sI", "sII", "liquid", "interface"]:
        data = np.load(FRAME_DIR / (name + ".npz"))
        for index in range(data["pos"].shape[0]):
            path = OUT / ("%s-%d.lammpstrj" % (name, index))
            write_frame(path, data["pos"][index], data["box"][index])
            jobs.append((path, data["pos"][index], data["box"][index]))
    for name in ["counterexample-five", "counterexample-interface", "counterexample-q3"]:
        data = np.load(FRAME_DIR / (name + ".npz"))
        path = OUT / (name + ".lammpstrj")
        write_frame(path, data["pos"], data["box"])
        jobs.append((path, data["pos"], data["box"]))
    return jobs


def typed_liquid():
    text = LIQUID.read_text().splitlines()
    typed = OUT / "liquid4096-typed.lammpstrj"
    lines = []
    atoms = False
    bounds = []
    pos = []
    timesteps = 0
    for line in text:
        if line.startswith("ITEM: TIMESTEP"):
            timesteps += 1
            if timesteps > 1:
                break
        if line.startswith("ITEM: ATOMS"):
            lines.append("ITEM: ATOMS id type x y z")
            atoms = True
            continue
        if atoms and line and not line.startswith("ITEM:"):
            parts = line.split()
            lines.append("%s 1 %s %s %s" % (parts[0], parts[1], parts[2], parts[3]))
            pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
        else:
            if line.startswith("ITEM:"):
                atoms = False
            lines.append(line)
            if len(bounds) < 3 and line and not line.startswith("ITEM:") and " " in line and not atoms:
                lo, hi = line.split()[:2]
                try:
                    bounds.append(float(hi) - float(lo))
                except ValueError:
                    pass
    typed.write_text("\n".join(lines) + "\n")
    return typed, np.asarray(pos, dtype=np.float64), np.asarray(bounds[:3], dtype=np.float64)


def main():
    jobs = frames()
    typed, liquid_pos, liquid_box = typed_liquid()
    changed = 0
    atoms = 0
    published_miss = [0]
    ovito_miss = [0]
    for path, pos, box in jobs:
        base = run_bin(os.environ["BASE_BIN"], os.environ["BASE_LIB"], "chill", path)
        tip = run_bin(os.environ["TIP_BIN"], os.environ["TIP_LIB"], "chill", path)
        save(Path(str(path) + ".base-chill.txt"), base)
        save(Path(str(path) + ".tip-chill.txt"), tip)
        published = library_rule(pos, box, "chill+")
        ovito_rule = library_rule(pos, box, "ovito-chill+")
        save(Path(str(path) + ".tip-published.txt"), published.tolist())
        save(Path(str(path) + ".tip-ovito.txt"), ovito_rule.tolist())
        reference = ovito_of(path)
        save(Path(str(path) + ".ovito.txt"), reference.tolist())
        if base != tip:
            changed += sum(a != b for a, b in zip(base, tip))
        atoms += len(base)
        published_miss[0] += int(np.sum(published != reference))
        ovito_miss[0] += int(np.sum(ovito_rule != reference))
    base4096 = run_bin(os.environ["BASE_BIN"], os.environ["BASE_LIB"], "chill", typed)
    tip4096 = run_bin(os.environ["TIP_BIN"], os.environ["TIP_LIB"], "chill", typed)
    save(OUT / "liquid4096.base-chill.txt", base4096)
    save(OUT / "liquid4096.tip-chill.txt", tip4096)
    if base4096 != tip4096:
        changed += sum(a != b for a, b in zip(base4096, tip4096))
    atoms += len(base4096)
    print("plain_chill frames", len(jobs) + 1, "atoms", atoms, "label_changes", changed)
    print("frames", len(jobs), "published_mismatch_atoms", published_miss[0], "ovito_rule_mismatch_atoms", ovito_miss[0])
    base_plus = run_bin(os.environ["BASE_BIN"], os.environ["BASE_LIB"], "chillplus", typed)
    save(OUT / "liquid4096.base-chillplus.txt", base_plus)
    published = library_rule(liquid_pos, liquid_box, "chill+")
    ovito_rule = library_rule(liquid_pos, liquid_box, "ovito-chill+")
    save(OUT / "liquid4096.tip-published.txt", published.tolist())
    save(OUT / "liquid4096.tip-ovito.txt", ovito_rule.tolist())
    ovito = ovito_of(typed)
    save(OUT / "liquid4096.ovito.txt", ovito.tolist())
    base_arr = np.asarray(base_plus, dtype=np.int32)
    print("liquid4096 n", len(base_arr))
    print("base_chillplus_ovito_agree", int(np.sum(base_arr == ovito)))
    print("tip_published_ovito_agree", int(np.sum(published == ovito)))
    print("tip_ovito_rule_agree", int(np.sum(ovito_rule == ovito)))


if __name__ == "__main__":
    main()
