# Copyright (c) d-SEAMS developers
"""Serial remeasurement of the results-section claims.

Subcommands write raw logs under benchmarks/cpc/records. Paths in those
logs stay relative to the repository.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (
    BIN,
    MERGED,
    NAMES,
    OVITO_MERGED,
    PY,
    RECORDS,
    REPEATS,
    ROOT,
    SEAMS,
    WARMUP,
    loadavg,
    parse_counts,
    parse_repeats,
    rel,
    run,
    scrub,
    summary,
)

FRAMES = ROOT / "benchmarks" / "cpc" / "frames"
TIP_LIB = ROOT.parent / "seams-chill-published" / "build" / "src" / "libyodaLib.so"


def log_path(measurement: str, name: str) -> Path:
    return RECORDS / measurement / name


def engine(mode: str, path: Path, frame: int = 1, typ: int = 1, kind: str = "time") -> str:
    cmd = [str(BIN), kind, mode, rel(path), str(frame), str(typ)]
    if kind == "time":
        cmd.extend([str(WARMUP), str(REPEATS)])
    proc = run(cmd, log_path(mode if kind == "time" else "labels", f"{path.stem}-{mode}-{kind}.log"))
    return proc.stdout


def read_dump(path: Path, frame: int = 1):
    lines = path.read_text().splitlines()
    seen = 0
    i = 0
    while i < len(lines):
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        while not lines[i].startswith("ITEM: NUMBER"):
            i += 1
        n = int(lines[i + 1].split()[0])
        i += 2
        while not lines[i].startswith("ITEM: BOX"):
            i += 1
        i += 1
        lo = []
        hi = []
        while not lines[i].startswith("ITEM: ATOMS"):
            parts = lines[i].split()
            lo.append(float(parts[0]))
            hi.append(float(parts[1]))
            i += 1
        header = lines[i].split()[2:]
        cols = {name: idx for idx, name in enumerate(header)}
        i += 1
        seen += 1
        if seen == frame:
            pos = np.empty((n, 3), dtype=np.float64)
            for row in range(n):
                parts = lines[i + row].split()
                pos[row] = [float(parts[cols["x"]]), float(parts[cols["y"]]), float(parts[cols["z"]])]
            origin = np.asarray(lo[:3])
            lengths = np.asarray(hi[:3]) - origin
            return pos, origin, lengths
        i += n
    raise RuntimeError(f"frame {frame} missing in {path.name}")


def n_frames(path: Path) -> int:
    return sum(1 for line in path.read_text().splitlines() if line.startswith("ITEM: TIMESTEP"))


def per_class(left: list[str], right: list[str], classes: list[str]) -> list[str]:
    rows = ["class both dseams ovito recall_given_ovito jaccard"]
    n = len(left)
    same = sum(a == b for a, b in zip(left, right))
    rows.append(f"overall {same} {n} {n} {same / n if n else 0:.6f} {same / n if n else 0:.6f}")
    for name in classes:
        both = sum(a == name and b == name for a, b in zip(left, right))
        nl = sum(a == name for a in left)
        nr = sum(b == name for b in right)
        union = nl + nr - both
        recall = both / nr if nr else float("nan")
        jaccard = both / union if union else float("nan")
        rows.append(f"{name} {both} {nl} {nr} {recall:.6f} {jaccard:.6f}")
    return rows


def ovito_chill(path: Path):
    import os
    import sys

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("OVITO_THREAD_COUNT", "1")
    sys.path.insert(0, str(Path(PY).parent))
    # Import from the venv by executing that interpreter when this process
    # is already the venv. Callers run this file with the venv python.
    from ovito.io import import_file
    from ovito.modifiers import ChillPlusModifier

    pipeline = import_file(str(path))
    samples = []
    data = None
    for index in range(WARMUP + REPEATS):
        pipeline.modifiers.clear()
        pipeline.modifiers.append(ChillPlusModifier(cutoff=3.5))
        t0 = time.perf_counter()
        data = pipeline.compute()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index >= WARMUP:
            samples.append(elapsed)
    labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
    return labels, samples


def tip_rule(pos, origin, lengths, rule: str):
    import ctypes

    lib = ctypes.CDLL(str(TIP_LIB))
    lib.seams_chill_plus_rule.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_double,
        ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus_rule.restype = ctypes.c_int
    n = pos.shape[0]
    xyz = np.ascontiguousarray(pos, dtype=np.float64).reshape(-1)
    box = np.array([origin[0], origin[1], origin[2], lengths[0], lengths[1], lengths[2]], dtype=np.float64)
    labels = np.empty(n, dtype=np.int32)
    rc = lib.seams_chill_plus_rule(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 3.5,
        b"cutoff", rule.encode(),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if rc != 0:
        raise RuntimeError("chill-plus rule failed")
    return labels


def names_from_ints(values) -> list[str]:
    return [NAMES[int(v)] if 0 <= int(v) < len(NAMES) else "unclassified" for v in values]


def write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(scrub("\n".join(lines) + "\n"))


def measure_a() -> None:
    systems = [
        ("ice-Ih", FRAMES / "ih-1000.lammpstrj"),
        ("ice-Ic", FRAMES / "ic-1000.lammpstrj"),
        ("sI", ROOT / "input" / "traj" / "genice_sI.lammpstrj"),
        ("mW-liquid-T215", FRAMES / "mw-liquid-512.lammpstrj"),
        ("mW-cubic", ROOT / "input" / "traj" / "mW_cubic.lammpstrj"),
    ]
    out = RECORDS / "a-chill"
    out.mkdir(parents=True, exist_ok=True)
    classes = ["cubic", "hexagonal", "water", "interfacial", "clathrate", "interClathrate",
               "unclassified", "reCubic", "reHex"]
    merged_classes = ["cubic", "hexagonal", "other", "interfacial", "hydrate", "interfacial-hydrate"]
    for name, path in systems:
        if not path.is_file():
            write_lines(out / f"{name}.status", [f"absent {path.name}"])
            continue
        for mode in ("chill", "chillplus"):
            labels = run(
                [str(BIN), "labels", mode, rel(path), "1", "1"],
                out / f"{name}-{mode}-labels.log",
            ).stdout
            ints = [int(line) for line in labels.splitlines() if line.strip().lstrip("-").isdigit()]
            (out / f"{name}-{mode}.txt").write_text("\n".join(str(v) for v in ints) + "\n")
            timed = run(
                [str(BIN), "time", mode, rel(path), "1", "1", str(WARMUP), str(REPEATS)],
                out / f"{name}-{mode}-time.log",
            ).stdout
            (out / f"{name}-{mode}-time.csv").write_text(timed)
        ov, samples = ovito_chill(path)
        (out / f"{name}-ovito.txt").write_text("\n".join(str(int(v)) for v in ov) + "\n")
        write_lines(out / f"{name}-ovito-time.log",
                    [f"loadavg {loadavg()}", "command ovito ChillPlusModifier cutoff=3.5"]
                    + [f"repeat {i} path_ms {ms:.6f}" for i, ms in enumerate(samples, 1)])
        dseams = names_from_ints(np.loadtxt(out / f"{name}-chillplus.txt", dtype=int))
        raw_ovito = [NAMES[{"0": 2, "1": 1, "2": 0, "3": 3, "4": 4, "5": 5}[str(int(v))]] for v in ov]
        # Map OVITO integers onto the d-SEAMS enum before the merged comparison.
        mapped = []
        for value in ov:
            mapped.append({0: "water", 1: "hexagonal", 2: "cubic", 3: "interfacial",
                           4: "clathrate", 5: "interClathrate"}[int(value)])
        write_lines(out / f"{name}-chillplus-vs-ovito-raw.tsv", per_class(dseams, mapped, classes))
        left = [MERGED.get(label, label) for label in dseams]
        right = [OVITO_MERGED[int(v)] for v in ov]
        write_lines(out / f"{name}-chillplus-vs-ovito-merged.tsv", per_class(left, right, merged_classes))
        chill = names_from_ints(np.loadtxt(out / f"{name}-chill.txt", dtype=int))
        write_lines(out / f"{name}-chill-vs-ovito-merged.tsv",
                    per_class([MERGED.get(label, label) for label in chill], right, merged_classes))
        if TIP_LIB.is_file():
            pos, origin, lengths = read_dump(path, 1)
            for rule in ("chill+", "ovito-chill+"):
                got = tip_rule(pos, origin, lengths, rule)
                (out / f"{name}-tip-{rule.replace('+','plus')}.txt").write_text(
                    "\n".join(str(int(v)) for v in got) + "\n")
                tip_names = names_from_ints(got)
                write_lines(
                    out / f"{name}-tip-{rule.replace('+','plus')}-vs-ovito-merged.tsv",
                    per_class([MERGED.get(label, label) for label in tip_names], right, merged_classes),
                )
        else:
            write_lines(out / f"{name}-tip.status", ["chill-plus library absent"])
    write_lines(out / "liquid-campaign.status", [
        "The comparator-campaign liquid dump is not in this checkout.",
        "mW-liquid-T215 is dseams2_repro/data/mW_liquid_T215_512.data, 512 mW sites.",
    ])
    print("a-chill")


def measure_b() -> None:
    out = RECORDS / "b-cages"
    out.mkdir(parents=True, exist_ok=True)
    frames = sorted(FRAMES.glob("ic-*.lammpstrj")) + sorted(FRAMES.glob("ih-*.lammpstrj"))
    frames.append(ROOT / "input" / "traj" / "genice_sI.lammpstrj")
    frames.append(ROOT / "input" / "traj" / "mW_cubic.lammpstrj")
    for path in frames:
        if not path.is_file():
            continue
        timed = run(
            [str(BIN), "time", "cages", rel(path), "1", "1", str(WARMUP), str(REPEATS)],
            out / f"{path.stem}.log",
        ).stdout
        (out / f"{path.stem}.csv").write_text(timed)
    print("b-cages")


def measure_c() -> None:
    out = RECORDS / "c-cluster"
    out.mkdir(parents=True, exist_ok=True)
    for path in sorted((ROOT / "input" / "traj").glob("*.lammpstrj")):
        frames = n_frames(path)
        typ = "2" if path.name == "exampleTraj.lammpstrj" else "1"
        for mode in ("cluster-q6", "cluster-chill"):
            timed = run(
                [str(BIN), "time", mode, rel(path), "1", typ, str(WARMUP), str(REPEATS)],
                out / f"{path.stem}-{mode}.log",
            ).stdout
            (out / f"{path.stem}-{mode}.csv").write_text(timed)
            largest = 0
            best_frame = 1
            scan = [f"file {path.name}", f"frames {frames}", f"mode {mode}"]
            for frame in range(1, frames + 1):
                text = run(
                    [str(BIN), "time", mode, rel(path), str(frame), typ, "0", "1"],
                    out / f"{path.stem}-{mode}-f{frame}.log",
                ).stdout
                size = 0
                for line in text.splitlines():
                    parts = line.split()
                    if len(parts) == 4 and parts[0] == "repeat" and parts[2] == "largest":
                        size = int(float(parts[3]))
                scan.append(f"frame {frame} largest {size}")
                if size >= largest:
                    largest = size
                    best_frame = frame
            scan.append(f"max_largest {largest}")
            scan.append(f"max_frame {best_frame}")
            write_lines(out / f"{path.stem}-{mode}-scan.txt", scan)
    print("c-cluster")


def measure_d() -> None:
    out = RECORDS / "d-qlm"
    out.mkdir(parents=True, exist_ok=True)
    frames = [
        FRAMES / "ic-1000.lammpstrj",
        FRAMES / "ih-1000.lammpstrj",
        FRAMES / "mw-liquid-512.lammpstrj",
        ROOT / "input" / "traj" / "genice_sI.lammpstrj",
        ROOT / "input" / "traj" / "mW_cubic.lammpstrj",
    ]
    for path in frames:
        if not path.is_file():
            continue
        for mode in ("q3", "q6", "q6avg"):
            timed = run(
                [str(BIN), "time", mode, rel(path), "1", "1", str(WARMUP), str(REPEATS)],
                out / f"{path.stem}-{mode}.log",
            ).stdout
            (out / f"{path.stem}-{mode}.csv").write_text(timed)
    print("d-qlm")


def seams_line(args: list[str], log: Path) -> str:
    proc = run([str(SEAMS), *args], log)
    return proc.stdout


def measure_f4() -> None:
    out = RECORDS / "f-f4"
    out.mkdir(parents=True, exist_ok=True)
    jobs = [
        ("exampleTraj", ROOT / "input" / "traj" / "exampleTraj.lammpstrj", ["--type", "2", "--htype", "1"]),
        ("genice-sI-oxygen", ROOT / "input" / "traj" / "genice_sI.lammpstrj", ["--type", "1", "--htype", "2"]),
        ("clathrate-thf", ROOT / "input" / "traj" / "clathrate-thf.lammpstrj", ["--type", "1", "--htype", "2"]),
    ]
    for name, path, flags in jobs:
        samples = []
        text = ""
        for index in range(WARMUP + REPEATS):
            t0 = time.perf_counter()
            proc = run([str(SEAMS), "f4", rel(path), *flags], out / f"{name}-r{index}.log")
            elapsed = (time.perf_counter() - t0) * 1000.0
            text = proc.stdout.strip()
            if index >= WARMUP:
                samples.append(elapsed)
        write_lines(out / f"{name}.csv",
                    [f"loadavg {loadavg()}", f"stdout {text}"]
                    + [f"repeat {i} process_ms {ms:.6f}" for i, ms in enumerate(samples, 1)])
    write_lines(out / "definition.txt", [
        "seams f4 is topoparam::rodgerF4 on the mutual four-nearest oxygen graph.",
        "exampleTraj keeps oxygen type 2 and hydrogen type 1, the pair used by tests/test_order_parameter.cpp.",
        "Oxygen-only genice sI has no H-O-O-H dihedral, so the mean is not a hydrate F4.",
        "tests/sota_f4_zeron.cpp reports Rodger F4 undefined on the oxygen-only mW lattices.",
    ])
    print("f-f4")


def write_points(path: Path, pos: np.ndarray, box: np.ndarray, types=None) -> None:
    n = len(pos)
    if types is None:
        types = np.ones(n, dtype=int)
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n" % n)
        handle.write("ITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write("0.0 %.12f\n" % float(length))
        handle.write("ITEM: ATOMS id type x y z\n")
        for i, (xyz, typ) in enumerate(zip(pos, types), start=1):
            handle.write("%d %d %.10f %.10f %.10f\n" % (i, int(typ), xyz[0], xyz[1], xyz[2]))


def cubic_diamond(reps: int):
    bond = 2.75
    a = 4.0 * bond / np.sqrt(3.0)
    fcc = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    basis = fcc + [(x + 0.25, y + 0.25, z + 0.25) for x, y, z in fcc]
    pos = [((i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a)
           for i in range(reps) for j in range(reps) for k in range(reps) for b in basis]
    box = np.array([reps * a] * 3)
    return np.asarray(pos) % box, box


def lonsdaleite(rx, ry, rz):
    bond = 2.75
    a = bond / ((3.0 / 8.0) * np.sqrt(8.0 / 3.0))
    c = a * np.sqrt(8.0 / 3.0)
    hexagonal = [(1 / 3, 2 / 3, 0), (2 / 3, 1 / 3, 0.5), (1 / 3, 2 / 3, 3 / 8), (2 / 3, 1 / 3, 7 / 8)]
    x_len, y_len = a, a * np.sqrt(3.0)
    basis = []
    for hx, hy, hz in hexagonal:
        x = (hx - 0.5 * hy) * a
        y = hy * np.sqrt(3.0) * a / 2.0
        for shift in ((0.0, 0.0), (x_len / 2.0, y_len / 2.0)):
            basis.append(((x + shift[0]) % x_len, (y + shift[1]) % y_len, hz * c))
    box = np.array([rx * x_len, ry * y_len, rz * c])
    pos = []
    for i in range(rx):
        for j in range(ry):
            for k in range(rz):
                for x, y, z in basis:
                    pos.append(((i * x_len + x) % box[0], (j * y_len + y) % box[1], (k * c + z) % box[2]))
    return np.asarray(pos), box


def dense_null(n, box, seed):
    rng = np.random.default_rng(seed)
    cell = 2.3
    pos = []
    while len(pos) < n:
        cand = rng.random(3) * box
        ok = True
        pool = pos[-400:] if len(pos) > 400 else pos
        for prev in pool:
            d = np.abs(cand - prev)
            d = np.minimum(d, box - d)
            if (d * d).sum() < cell * cell:
                ok = False
                break
        if ok:
            pos.append(cand)
    return np.asarray(pos)


def parse_cages(text: str) -> dict[str, int]:
    # nop N graph KIND hexagonal H cubic C both B water W
    parts = text.split()
    out = {}
    for key in ("nop", "hexagonal", "cubic", "both", "water"):
        if key in parts:
            out[key] = int(parts[parts.index(key) + 1])
    return out


def score_counts(counts: dict[str, int], truth: str) -> str:
    n = counts.get("nop", 0)
    cubic = counts.get("cubic", 0) / n if n else float("nan")
    hexa = counts.get("hexagonal", 0) / n if n else float("nan")
    both = counts.get("both", 0) / n if n else float("nan")
    # truth matches only the disjoint class, as score() in sota_compare.py does.
    acc = {"cubic": cubic, "hexagonal": hexa, "water": counts.get("water", 0) / n if n else float("nan")}[truth]
    crystal = cubic + hexa
    return f"acc={acc:.6f} cubic={cubic:.6f} hex={hexa:.6f} both={both:.6f} crystal={crystal:.6f}"


def measure_sota() -> None:
    out = RECORDS / "e-sota"
    frames = out / "frames"
    frames.mkdir(parents=True, exist_ok=True)
    ic, ic_box = cubic_diamond(4)
    ih, ih_box = lonsdaleite(5, 3, 3)
    null = dense_null(360, ic_box, 88172645463325252)
    base = {"ic": (ic, ic_box, "cubic"), "ih": (ih, ih_box, "hexagonal"), "null": (null, ic_box, "water")}
    sigmas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    rows = ["method system sigma n " + "acc cubic hex both crystal median_ms q1_ms q3_ms iqr_ms"]
    for system, (pos0, box, truth) in base.items():
        for sigma in sigmas:
            rng = np.random.default_rng(88172645463325252 + int(sigma * 1000))
            pos = pos0.copy()
            if sigma:
                pos = (pos + rng.normal(0.0, sigma, pos.shape)) % box
            path = frames / f"{system}-{sigma:.2f}.lammpstrj"
            write_points(path, pos, box)
            samples = []
            text = ""
            for index in range(WARMUP + REPEATS):
                t0 = time.perf_counter()
                proc = run(
                    [str(SEAMS), "cages", "--graph", "cutoff", rel(path)],
                    out / f"{system}-{sigma:.2f}-r{index}.log",
                )
                if index >= WARMUP:
                    samples.append((time.perf_counter() - t0) * 1000.0)
                text = proc.stdout
            counts = parse_cages(text)
            stat = summary(samples)
            rows.append(
                f"dseams-topo {system} {sigma:.2f} {counts.get('nop', 0)} {score_counts(counts, truth)} "
                f"median_ms={stat['median']:.6f} q1_ms={stat['q1']:.6f} q3_ms={stat['q3']:.6f} "
                f"iqr_ms={stat['iqr']:.6f}"
            )
            # CHILL+ labels, one classification plus the in-process repeats.
            label_text = run(
                [str(BIN), "labels", "chillplus", rel(path), "1", "1"],
                out / f"{system}-{sigma:.2f}-chillplus-labels.log",
            ).stdout
            ints = [int(line) for line in label_text.splitlines() if line.strip().lstrip("-").isdigit()]
            chill_names = [MERGED.get(NAMES[v], "other") for v in ints]
            n = len(chill_names)
            acc = sum(label == truth or (truth == "water" and label == "other") for label in chill_names) / n
            cubic = sum(label == "cubic" for label in chill_names) / n
            hexa = sum(label == "hexagonal" for label in chill_names) / n
            timed = run(
                [str(BIN), "time", "chillplus", rel(path), "1", "1", str(WARMUP), str(REPEATS)],
                out / f"{system}-{sigma:.2f}-chillplus-time.log",
            ).stdout
            series = parse_repeats(timed)["path_ms"]
            stat = summary(series)
            rows.append(
                f"chill+ {system} {sigma:.2f} {n} acc={acc:.6f} cubic={cubic:.6f} hex={hexa:.6f} "
                f"both=nan crystal={cubic + hexa:.6f} median_ms={stat['median']:.6f} "
                f"q1_ms={stat['q1']:.6f} q3_ms={stat['q3']:.6f} iqr_ms={stat['iqr']:.6f}"
            )
    write_lines(out / "sota.tsv", rows)
    print("e-sota")


def measure_ions() -> None:
    out = RECORDS / "e-ions"
    out.mkdir(parents=True, exist_ok=True)
    hosts = {
        "Ic": (cubic_diamond(3), "cubic"),
        "Ih": (lonsdaleite(5, 3, 3), "hexagonal"),
    }
    rows = ["host pairs draw n_water n_ion ice in_ice front liquid shell"]
    seed = 20260902
    for host, ((sites, box), _truth) in hosts.items():
        for pairs in (0, 1, 2, 4):
            draws = 1 if pairs == 0 else 3
            for draw in range(draws):
                rng = np.random.default_rng(seed + 1000 * pairs + draw)
                types = np.ones(len(sites), dtype=int)
                if pairs:
                    chosen = rng.choice(len(sites), size=2 * pairs, replace=False)
                    types[chosen[:pairs]] = 3
                    types[chosen[pairs:]] = 4
                path = out / "frames" / f"{host}-p{pairs}-d{draw}.lammpstrj"
                path.parent.mkdir(parents=True, exist_ok=True)
                write_points(path, sites, box, types)
                samples = []
                text = ""
                for index in range(WARMUP + REPEATS):
                    proc = run(
                        [str(SEAMS), "ions", rel(path), "--ion-types", "3,4",
                         "--type", "1", "--complete"],
                        out / f"{host}-p{pairs}-d{draw}-r{index}.log",
                    )
                    text = proc.stdout.strip()
                    if index >= WARMUP:
                        samples.append(text)
                write_lines(out / f"{host}-p{pairs}-d{draw}.csv",
                            [f"stdout {text}"] + [f"repeat {i} {line}" for i, line in enumerate(samples, 1)])
                rows.append(f"{host} {pairs} {draw} {text}")
    write_lines(out / "ions.tsv", rows)
    print("e-ions")


def measure_brine() -> None:
    out = RECORDS / "e-brine"
    out.mkdir(parents=True, exist_ok=True)
    root = ROOT.parent / "dseams2_repro" / "results" / "brine"
    script = ROOT.parent / "dseams2_repro" / "scripts" / "brine_summary.py"
    if not root.is_dir() or not script.is_file():
        write_lines(out / "status.txt", [
            "BRINE_INPUTS_ABSENT",
            "Coordinate trajectories of the ice/brine campaign are not in the checkout.",
            "results/brine holds the stored DSEAMS_CAGES columns, not a fresh classification.",
        ])
        print("e-brine absent")
        return
    proc = run(["benchmarks/cpc/.venv/bin/python", rel(script), rel(root)], out / "brine_summary.log")
    (out / "summary.tsv").write_text(proc.stdout)
    write_lines(out / "status.txt", [
        "The atomic brine trajectories are not in dseams2_repro.",
        "summary.tsv is brine_summary.py reread of the stored BRINE columns.",
        "It is not a new engine classification. Reference conditions name GenIce Ih TIP4P/2005, 1152 molecules, engine v2.9.2.",
    ])
    print("e-brine")


def measure_gpu() -> None:
    out = RECORDS / "e-gpu"
    out.mkdir(parents=True, exist_ok=True)
    write_lines(out / "host-device-ladder.txt", [
        "GPU_REQUIRED",
        "Host-device ladder n=1000 to 32000 was not run.",
        "nvidia-smi is absent and /dev/nvidia* is absent.",
        "No device times were invented.",
    ])
    print("e-gpu")


def measure_e() -> None:
    out = RECORDS / "e-other"
    out.mkdir(parents=True, exist_ok=True)
    # Clathrate cage census on the shipped sI and sII unit cells.
    for name in ("genice_sI.lammpstrj", "genice_sII.lammpstrj"):
        path = ROOT / "input" / "traj" / name
        samples = []
        text = ""
        for index in range(WARMUP + REPEATS):
            t0 = time.perf_counter()
            proc = run(
                [str(SEAMS), "cages", "--graph", "cutoff", "--signature", "512+51262+51264", rel(path)],
                out / f"{path.stem}-sig-r{index}.log",
            )
            text = proc.stdout
            if index >= WARMUP:
                samples.append((time.perf_counter() - t0) * 1000.0)
        stat = summary(samples)
        write_lines(out / f"{path.stem}-signature.txt",
                    [text.strip(),
                     f"median_ms {stat['median']:.6f}",
                     f"q1_ms {stat['q1']:.6f}",
                     f"q3_ms {stat['q3']:.6f}",
                     f"iqr_ms {stat['iqr']:.6f}"])
    write_lines(out / "not-measured.txt", [
        "Neighbour-list and ring head-to-heads are outside this run.",
        "d-SEAMS 1 commit 6e51956f is not in this checkout, so the version-comparison ratios were not recomputed.",
        "The 148-frame nucleation trajectory is not shipped; only input/traj/nucleation_f1.lammpstrj is present.",
        "Parallel strong-scaling speedups were not run; the protocol is one thread.",
        "submission/cpc-v3 and REMOVED.tsv are not on a reachable remote, so the 57 v1 decimals were not listed.",
    ])
    print("e-other")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("which", choices=["a", "b", "c", "d", "e", "f4", "sota", "ions", "brine", "gpu"])
    args = parser.parse_args()
    {
        "a": measure_a,
        "b": measure_b,
        "c": measure_c,
        "d": measure_d,
        "e": measure_e,
        "f4": measure_f4,
        "sota": measure_sota,
        "ions": measure_ions,
        "brine": measure_brine,
        "gpu": measure_gpu,
    }[args.which]()


if __name__ == "__main__":
    main()
