#!/usr/bin/env python3
"""Time the seams commands that sit on the hot kernels.

One process starts each sample. The first call is a warm-up and is not
recorded. Eight samples follow. The report is the median and the
interquartile range. taskset pins the process to the first N cores and
OMP_NUM_THREADS is that same N.

Structures are written with tempfile. Reference times from mdapy, OVITO
and rings-statistics-matter are the published comparison points; this
script records the seams side of that comparison.
"""

import math
import os
import statistics
import subprocess
import tempfile
from pathlib import Path


def ice_ic(ncell, path):
    bond = 2.75
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = ncell * lattice
    n = ncell * ncell * ncell * 8
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    with path.open("w") as out:
        out.write("ITEM: TIMESTEP\n0\n")
        out.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
        out.write("ITEM: BOX BOUNDS pp pp pp\n")
        out.write(f"0 {side}\n0 {side}\n0 {side}\n")
        out.write("ITEM: ATOMS id type x y z\n")
        aid = 1
        for i in range(ncell):
            for j in range(ncell):
                for k in range(ncell):
                    for b in range(8):
                        fx, fy, fz = fcc[b % 4]
                        shift = 0.0 if b < 4 else 0.25
                        x = ((i + fx + shift) * lattice) % side
                        y = ((j + fy + shift) * lattice) % side
                        z = ((k + fz + shift) * lattice) % side
                        out.write(f"{aid} 1 {x:.8f} {y:.8f} {z:.8f}\n")
                        aid += 1
    return n


def percentile(xs, p):
    k = (len(xs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] * (1 - (k - lo)) + xs[hi] * (k - lo)


def samples(seams, args, dump, cores, n=8):
    cpus = ",".join(str(i) for i in range(cores))
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_WAIT_POLICY"] = "active"
    cmd = ["taskset", "-c", cpus, str(seams), *args, str(dump), "--type", "1"]
    subprocess.run(cmd, env=env, check=True, stdout=subprocess.DEVNULL)
    import time

    got = []
    for _ in range(n):
        t0 = time.perf_counter()
        subprocess.run(
            cmd,
            env=env,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        got.append((time.perf_counter() - t0) * 1000.0)
    got.sort()
    med = statistics.median(got)
    iqr = percentile(got, 0.75) - percentile(got, 0.25)
    return med, iqr


def main():
    root = Path(__file__).resolve().parents[1]
    seams = root / "build" / "seams"
    if not seams.is_file():
        raise SystemExit("build/seams is missing")
    cores_available = os.cpu_count() or 1
    with tempfile.TemporaryDirectory(prefix="seams-cli-") as tmp:
        big = Path(tmp) / "iceIc.lammpstrj"
        small = Path(tmp) / "ice1000.lammpstrj"
        n_big = ice_ic(24, big)
        n_small = ice_ic(5, small)
        rows = []
        plan = [
            ("read", ["read"], big, n_big),
            ("chill-plus", ["chill-plus"], big, n_big),
            ("symmetry", ["symmetry"], big, n_big),
            ("steinhardt --l 4,6", ["steinhardt", "--l", "4,6"], big, n_big),
            ("cages", ["cages"], big, n_big),
            ("rings cutoff depth 6", None, small, n_small),
        ]
        lib = root / "build" / "bindings" / "libseams_lang.so"
        for cores in (1, 8):
            if cores > cores_available:
                continue
            for name, args, dump, n in plan:
                if args is None:
                    med, iqr = ring_samples(lib, dump, cores)
                else:
                    med, iqr = samples(seams, args, dump, cores)
                rows.append((name, n, cores, med, iqr))
                print(
                    f"{name} n={n} cores={cores} median_ms={med:.2f} iqr_ms={iqr:.2f}",
                    flush=True,
                )
    out = root / "markdown" / "cli-fast-paths.md"
    out.write_text(render(rows, cores_available))
    plot(rows, root / "markdown" / "cli-fast-plots")
    print(f"wrote {out.relative_to(root)}")


def ring_samples(lib, dump, cores, n=8):
    """Primitive rings through the C ABI, neighbour list included."""
    import textwrap

    code = textwrap.dedent(
        """
        import ctypes, sys, time
        lib = ctypes.CDLL(sys.argv[1])
        lib.seams_frame_read.restype = ctypes.c_void_p
        lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_rings.restype = ctypes.c_int
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        counts = (ctypes.c_int * 9)()
        frame = lib.seams_frame_read(sys.argv[2].encode(), 1, 1)
        t0 = time.perf_counter()
        rc = lib.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, counts, 9)
        dt = (time.perf_counter() - t0) * 1000.0
        lib.seams_frame_free(frame)
        if rc != 0:
            raise SystemExit(rc)
        print(f"{dt:.6f}")
        """
    )
    cpus = ",".join(str(i) for i in range(cores))
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_WAIT_POLICY"] = "active"
    cmd = [
        "taskset",
        "-c",
        cpus,
        "python3",
        "-c",
        code,
        str(lib),
        str(dump),
    ]
    subprocess.run(cmd, env=env, check=True, stdout=subprocess.DEVNULL)
    import time

    got = []
    for _ in range(n):
        proc = subprocess.run(
            cmd, env=env, check=True, capture_output=True, text=True
        )
        got.append(float(proc.stdout.strip()))
    got.sort()
    return statistics.median(got), percentile(got, 0.75) - percentile(got, 0.25)


def render(rows, cores_available):
    lines = [
        "# CLI times",
        "",
        "Ice Ic, orthorhombic, one oxygen per molecule. "
        f"The machine reported {cores_available} cores. "
        "Each row is one warm-up and eight samples, median and interquartile range, "
        "with the process pinned to that many cores and `OMP_NUM_THREADS` set to the same number.",
        "",
        "Cited reference points, not remeasured in this run:",
        "",
        "- mdapy 2.0.0 cutoff CHILL+, 65.6 ms at 1 thread, on the 110592-molecule cell.",
        "- OVITO 3.17.2 polyhedral template matching, 208.7 ms at 8 threads, same cell.",
        "- rings-statistics-matter 0.2.0, 3.56 ms at 8 threads on 1000 molecules. "
        "The previous d-SEAMS figure on that cell was 4.02 ms.",
        "",
        "| command | molecules | cores | median (ms) | IQR (ms) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for name, n, cores, med, iqr in rows:
        lines.append(
            f"| `{name}` | {n} | {cores} | {med:.2f} | {iqr:.2f} |"
        )
    lines.append("")
    return "\n".join(lines)


def plot(rows, directory):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    directory.mkdir(parents=True, exist_ok=True)
    names = []
    for name, _n, cores, _med, _iqr in rows:
        if cores == 1 and name not in names:
            names.append(name)
    for cores, marker in ((1, "o"), (8, "s")):
        xs = []
        ys = []
        for name in names:
            for row_name, _n, row_cores, med, _iqr in rows:
                if row_name == name and row_cores == cores:
                    xs.append(name)
                    ys.append(med)
        if ys:
            plt.plot(xs, ys, marker=marker, label=f"{cores} cores")
    plt.ylabel("median (ms)")
    plt.xticks(rotation=20, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(directory / "cli-fast.png", dpi=120)
    plt.close()


if __name__ == "__main__":
    main()
