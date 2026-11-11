#!/usr/bin/env python3
"""One session of the block cut and the linkcell cell-grid cut.

Strong scaling uses 110592 and 884736 molecules at 1, 2, 4 and 8 ranks.
Weak scaling keeps 110592 molecules on each rank: 24x24x24, 24x24x48,
24x48x48 and 48x48x48 cells. One thread per rank. One warm-up, six
timed calls. The block-cut library runs, then the cell-grid library,
before the next configuration.

The host is one Intel Xeon socket in a KVM guest: 8 cores, 1 thread
per core, AVX-512, 320 MiB of L3, and MemTotal 49342816 kB.

    python3 paper/tex/repro/partition_session.py

The benchmark line is:

    OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE \\
    OMP_DYNAMIC=FALSE OMPI_MCA_btl=^openib OMPI_MCA_pml=ob1 \\
    mpirun --bind-to core --map-by node:PE=1 -n N \\
        /tmp/rank_bench 24 6 rings
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "paper" / "tex" / "figures" / "data"
OUT = ROOT / "paper" / "tex" / "figures"
CSV_PATH = DATA / "partition-session.csv"
TABLE_PATH = DATA / "tab-partition.tex"
PROSE_PATH = DATA / "partition-prose.tex"
BENCH = Path("/tmp/rank_bench")
BUILD_SO = Path("/tmp/builds/chill-mpi/src/libyodaLib.so")
SLAB = Path("/tmp/lib_slab.so")
PART = Path("/tmp/lib_part.so")
LIBS = ":".join(
    [
        "/tmp/builds/chill-mpi/src",
        "/tmp/builds/chill-mpi/subprojects/readcon-core",
        "/tmp/builds/chill-mpi/subprojects/readcon-db",
        "/tmp/builds/chill-mpi/subprojects/linkcell",
        "/tmp/builds/chill-mpi/subprojects/vesin",
        "/tmp/builds/chill-mpi/subprojects/minimage",
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
        "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib",
    ]
)
OPS = ("neigh", "chill", "q6", "rings")
OP_LABEL = {"neigh": "neighbours", "chill": "CHILL+", "q6": "$q_6$", "rings": "rings"}
STRONG = ((24, 24, 24, 110592), (48, 48, 48, 884736))
WEAK = ((1, 24, 24, 24), (2, 24, 24, 48), (4, 24, 48, 48), (8, 48, 48, 48))
REPS = 6
COLOURS = {"neigh": "#0072B2", "chill": "#D55E00", "q6": "#009E73", "rings": "#CC79A7"}


def quantile_linear(ordered, q):
    if len(ordered) == 1:
        return ordered[0]
    pos = q * (len(ordered) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(ordered):
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac


def summarize(values):
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return median, q1, q3, q3 - q1


def use(lib):
    shutil.copyfile(lib, BUILD_SO)


def run_bench(mx, my, mz, ranks, op):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env.pop("SEAMS_RING_PROFILE", None)
    env["LD_LIBRARY_PATH"] = LIBS
    env["OMP_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = "1"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMPI_MCA_btl"] = "^openib"
    env["OMPI_MCA_pml"] = "ob1"
    cmd = [
        "mpirun",
        "--bind-to",
        "core",
        "--map-by",
        "node:PE=1",
        "-n",
        str(ranks),
        str(BENCH),
        str(mx),
        str(my),
        str(mz),
        str(REPS),
        op,
    ]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, cwd=ROOT / "input")
    with (DATA / "partition-raw.log").open("a") as raw:
        raw.write(f"# {mx} {my} {mz} ranks {ranks} {op}\n")
        raw.write(proc.stdout)
        if proc.stderr:
            raw.write(proc.stderr)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"benchmark failed ({proc.returncode}): {' '.join(cmd)}")
    walls, halos, imbalances = [], [], []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        walls.append(float(parts[0]))
        halos.append(float(parts[3]))
        imbalances.append(float(parts[4]) if len(parts) >= 5 else 1.0)
    if len(walls) < REPS:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"expected {REPS} samples, got {len(walls)}")
    return walls, halos, imbalances


def row_of(kind, scale, molecules, ranks, op, walls, halos, imbalances):
    median, q1, q3, iqr = summarize(walls)
    halo, _, _, _ = summarize(halos)
    imbalance, _, _, _ = summarize(imbalances)
    return {
        "cut": kind,
        "scale": scale,
        "molecules": molecules,
        "ranks": ranks,
        "operation": op,
        "median_ms": f"{median:.6f}",
        "q1_ms": f"{q1:.6f}",
        "q3_ms": f"{q3:.6f}",
        "iqr_ms": f"{iqr:.6f}",
        "halo_fraction": f"{halo:.6f}",
        "imbalance": f"{imbalance:.6f}",
        "speedup": "",
        "efficiency": "",
    }


def fill_efficiency(rows):
    base = {}
    for row in rows:
        if int(row["ranks"]) == 1:
            base[(row["cut"], row["scale"], row["operation"], int(row["molecules"]))] = float(
                row["median_ms"]
            )
    for row in rows:
        if row["scale"] == "strong":
            key = (row["cut"], row["scale"], row["operation"], int(row["molecules"]))
            atoms = int(row["molecules"])
        else:
            # Weak efficiency is against the one-rank time of the same cut and operation.
            key = (row["cut"], row["scale"], row["operation"], 110592)
            atoms = int(row["molecules"])
        origin = base[key]
        med = float(row["median_ms"])
        n = int(row["ranks"])
        if row["scale"] == "strong":
            row["speedup"] = f"{origin / med:.6f}"
            row["efficiency"] = f"{origin / (med * n):.6f}"
        else:
            row["speedup"] = f"{origin / med:.6f}"
            row["efficiency"] = f"{origin / med:.6f}"
        row["molecules"] = atoms


def measure():
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "partition-raw.log").write_text("")
    rows = []
    for mx, my, mz, molecules in STRONG:
        for ranks in (1, 2, 4, 8):
            for op in OPS:
                for kind, lib in (("block", SLAB), ("cell", PART)):
                    use(lib)
                    print(f"strong {kind} {molecules} n={ranks} {op}", flush=True)
                    walls, halos, imbalances = run_bench(mx, my, mz, ranks, op)
                    rows.append(row_of(kind, "strong", molecules, ranks, op, walls, halos, imbalances))
                    dump(rows)
    for ranks, mx, my, mz in WEAK:
        molecules = mx * my * mz * 8
        for op in OPS:
            for kind, lib in (("block", SLAB), ("cell", PART)):
                use(lib)
                print(f"weak {kind} {molecules} n={ranks} {op}", flush=True)
                walls, halos, imbalances = run_bench(mx, my, mz, ranks, op)
                rows.append(row_of(kind, "weak", molecules, ranks, op, walls, halos, imbalances))
                dump(rows)
    fill_efficiency(rows)
    dump(rows)
    return rows


def dump(rows):
    DATA.mkdir(parents=True, exist_ok=True)
    fields = [
        "cut",
        "scale",
        "molecules",
        "ranks",
        "operation",
        "median_ms",
        "q1_ms",
        "q3_ms",
        "iqr_ms",
        "speedup",
        "efficiency",
        "halo_fraction",
        "imbalance",
    ]
    with CSV_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def pick(rows, scale, molecules, op, cut, ranks):
    for row in rows:
        if (
            row["scale"] == scale
            and int(row["molecules"]) == molecules
            and row["operation"] == op
            and row["cut"] == cut
            and int(row["ranks"]) == ranks
        ):
            return row
    raise KeyError((scale, molecules, op, cut, ranks))


def write_table(rows):
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{Strong scaling of the block cut and the cell-grid cut.}",
        "Both run in the same session, one thread per rank.",
        "Efficiency is the one-rank time divided by the eight-rank time and by eight.",
        "The imbalance is the eight-rank ratio of the cell-grid cut.}",
        "\\label{tab:partition-strong}",
        "\\begin{tabular}{llrrrrrrr}",
        "\\toprule",
        "Molecules & Operation & Block 1 & Block 8 & Block eff. & Cell 1 & Cell 8 & Cell eff. & Imbalance \\\\",
        "\\midrule",
    ]
    last = None
    for molecules in (110592, 884736):
        for op in OPS:
            b1 = float(pick(rows, "strong", molecules, op, "block", 1)["median_ms"])
            b8 = float(pick(rows, "strong", molecules, op, "block", 8)["median_ms"])
            c1 = float(pick(rows, "strong", molecules, op, "cell", 1)["median_ms"])
            c8 = float(pick(rows, "strong", molecules, op, "cell", 8)["median_ms"])
            imb = float(pick(rows, "strong", molecules, op, "cell", 8)["imbalance"])
            shown = str(molecules) if molecules != last else ""
            last = molecules
            lines.append(
                f"{shown} & {OP_LABEL[op]} & {b1:.3f} & {b8:.3f} & {b1 / (b8 * 8):.3f} "
                f"& {c1:.3f} & {c8:.3f} & {c1 / (c8 * 8):.3f} & {imb:.3f} \\\\"
            )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    lines.extend(
        [
            "\\begin{table}[tbp]",
            "\\centering",
            "\\caption{\\textbf{Weak scaling at 110592 molecules per rank.}",
            "Efficiency is the one-rank time divided by the eight-rank time.",
            "The imbalance is the eight-rank ratio of the cell-grid cut.}",
            "\\label{tab:partition-weak}",
            "\\begin{tabular}{lrrrrr}",
            "\\toprule",
            "Operation & Block 8 (ms) & Block eff. & Cell 8 (ms) & Cell eff. & Imbalance \\\\",
            "\\midrule",
        ]
    )
    for op in OPS:
        b1 = float(pick(rows, "weak", 110592, op, "block", 1)["median_ms"])
        b8 = float(pick(rows, "weak", 884736, op, "block", 8)["median_ms"])
        c1 = float(pick(rows, "weak", 110592, op, "cell", 1)["median_ms"])
        c8 = float(pick(rows, "weak", 884736, op, "cell", 8)["median_ms"])
        imb = float(pick(rows, "weak", 884736, op, "cell", 8)["imbalance"])
        lines.append(
            f"{OP_LABEL[op]} & {b8:.3f} & {b1 / b8:.3f} & {c8:.3f} & {c1 / c8:.3f} & {imb:.3f} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    TABLE_PATH.write_text("\n".join(lines))


def write_prose(rows):
    bits = []
    for molecules in (110592, 884736):
        for op in OPS:
            cell = pick(rows, "strong", molecules, op, "cell", 8)
            block = pick(rows, "strong", molecules, op, "block", 8)
            bits.append(
                f"At {molecules} molecules the eight-rank {OP_LABEL[op]} wall is "
                f"{float(block['median_ms']):.3f} ms on the block cut "
                f"(efficiency {float(block['efficiency']):.3f}, imbalance {float(block['imbalance']):.3f}) "
                f"and {float(cell['median_ms']):.3f} ms on the cell-grid cut "
                f"(efficiency {float(cell['efficiency']):.3f}, imbalance {float(cell['imbalance']):.3f}, "
                f"halo {float(cell['halo_fraction']):.3f})."
            )
    for op in OPS:
        cell = pick(rows, "weak", 884736, op, "cell", 8)
        block = pick(rows, "weak", 884736, op, "block", 8)
        bits.append(
            f"Weak scaling at eight ranks, {OP_LABEL[op]}: the block cut is "
            f"{float(block['median_ms']):.3f} ms (efficiency {float(block['efficiency']):.3f}) "
            f"and the cell-grid cut is {float(cell['median_ms']):.3f} ms "
            f"(efficiency {float(cell['efficiency']):.3f}, imbalance {float(cell['imbalance']):.3f})."
        )
    PROSE_PATH.write_text("\n\n".join(bits) + "\n")


def write_figure(rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    for ax, molecules, title in (
        (axes[0], 110592, "110592 molecules"),
        (axes[1], 884736, "884736 molecules"),
    ):
        ax.plot([1, 8], [1, 8], color="#000000", lw=0.8, ls="--", label="ideal")
        for op in OPS:
            for cut, style, marker in (("cell", "-", "o"), ("block", "--", "s")):
                series = [
                    row
                    for row in rows
                    if row["scale"] == "strong"
                    and int(row["molecules"]) == molecules
                    and row["operation"] == op
                    and row["cut"] == cut
                ]
                series.sort(key=lambda row: int(row["ranks"]))
                label = None
                if molecules == 110592:
                    label = f"{OP_LABEL[op].replace('$', '')} {cut}"
                ax.plot(
                    [int(row["ranks"]) for row in series],
                    [float(row["speedup"]) for row in series],
                    marker=marker,
                    ms=3.5,
                    ls=style,
                    color=COLOURS[op],
                    lw=1.2,
                    label=label,
                )
        ax.set_title(title)
        ax.set_xlabel("Ranks (1 thread each)")
        ax.set_ylabel("Speedup vs 1 rank")
        ax.set_xticks([1, 2, 4, 8])
        ax.set_ylim(0, 9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, ncol=2)
    fig.savefig(OUT / "partition-scale.pdf")
    fig.savefig(OUT / "partition-scale.svg")
    plt.close(fig)


def main():
    for path in (BENCH, SLAB, PART, BUILD_SO):
        if not path.exists():
            raise SystemExit(f"missing {path}")
    if "--plot-only" in sys.argv:
        with CSV_PATH.open() as handle:
            rows = list(csv.DictReader(handle))
    else:
        try:
            rows = measure()
        finally:
            if PART.exists():
                use(PART)
    write_table(rows)
    write_prose(rows)
    write_figure(rows)
    print(f"wrote {TABLE_PATH}")
    print(f"wrote {OUT / 'partition-scale.pdf'}")


if __name__ == "__main__":
    main()
