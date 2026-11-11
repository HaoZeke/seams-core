#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Serial CHILL+ timings against OVITO ChillPlusModifier.

One warm-up, then six fresh processes, interleaved. The frame is read
before the timer. Neighbour and ring rows are copied from summary.tsv.
"""

from __future__ import annotations

import math
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "h2h-cpc" / "raw" / "chill-plus"
BENCH = ROOT / "bbdir" / "h2h-release" / "bench_h2h_cpc"
SYSTEMS = ROOT / "h2h-cpc" / "systems"
SUMMARY = ROOT / "h2h-cpc" / "summary.tsv"
OUT = ROOT / "h2h-cpc" / "timings.csv"
CUTOFF = 3.5
WARMUP = 1
REPEATS = 6
FRAMES = (
    "ice-1k", "ice-10k", "ice-50k", "ice-110k",
    "liq-1k", "liq-10k", "liq-50k", "liq-110k",
    "interface-50k",
)
TOOLS = ("dseams", "ovito")
ENV = (
    "OMP_NUM_THREADS=1",
    "OMP_DYNAMIC=FALSE",
    "OVITO_THREAD_COUNT=1",
    "OPENBLAS_NUM_THREADS=1",
    "MKL_NUM_THREADS=1",
    "NUMEXPR_NUM_THREADS=1",
    "PYTHONPATH=.h2h-py",
)


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    count = len(ordered)
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def loadavg() -> str:
    return Path("/proc/loadavg").read_text(encoding="utf-8").strip()


def cpu_model() -> str:
    for line in subprocess.check_output(["lscpu"], text=True).splitlines():
        if line.startswith("Model name:"):
            return line
    return "Model name:"


def parse_kv(text: str) -> dict[str, str]:
    record = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition(" ")
        record[key] = value
    return record


def run(argv: list[str]) -> dict[str, str]:
    before = loadavg()
    started = time.perf_counter()
    done = subprocess.run(argv, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall = (time.perf_counter() - started) * 1000.0
    record = parse_kv(done.stdout)
    record["loadavg_before"] = before
    record["loadavg_after"] = loadavg()
    record["wall_ms"] = f"{wall:.3f}"
    record["exit"] = str(done.returncode)
    record["command"] = " ".join(argv)
    if done.returncode != 0:
        raise SystemExit(done.stderr[-2000:])
    return record


def dseams_argv(dump: Path) -> list[str]:
    return [
        "taskset", "-c", "0", "env", *ENV,
        str(BENCH.relative_to(ROOT)), "chill-plus",
        str(dump.relative_to(ROOT)), f"{CUTOFF:.8f}", "6",
        "h2h-cpc/raw/chill-plus/last-dseams.csv",
    ]


def ovito_argv(dump: Path) -> list[str]:
    return [
        "taskset", "-c", "0", "env", *ENV,
        "python3.12", "scripts/h2h_chill.py", "worker",
        str(dump.relative_to(ROOT)), f"{CUTOFF:.8f}",
    ]


def worker(dump: str, cutoff: float) -> None:
    os.environ["OVITO_THREAD_COUNT"] = "1"
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource

    lines = Path(dump).read_text(encoding="utf-8").splitlines()
    n = int(lines[3])
    box = [float(lines[5 + k].split()[1]) - float(lines[5 + k].split()[0]) for k in range(3)]
    pos = np.empty((n, 3), dtype=np.float64)
    for i, line in enumerate(lines[9 : 9 + n]):
        parts = line.split()
        pos[i] = (float(parts[2]), float(parts[3]), float(parts[4]))
    data = DataCollection()
    data.create_cell(
        [[box[0], 0, 0], [0, box[1], 0], [0, 0, box[2]]],
        pbc=(True, True, True),
    )
    parts = data.create_particles(count=n)
    parts.create_property("Position", data=pos)
    modifier = ChillPlusModifier(cutoff=cutoff)
    names = {int(s.id): s.name for s in modifier.structures}
    canon = {
        "Other": "other",
        "Hexagonal ice": "hexagonal",
        "Cubic ice": "cubic",
        "Interfacial ice": "interfacial",
        "Hydrate": "hydrate",
        "Interfacial hydrate": "interfacial-hydrate",
    }
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(modifier)
    started = time.perf_counter()
    out = pipe.compute()
    elapsed = (time.perf_counter() - started) * 1000.0
    types = np.asarray(out.particles["Structure Type"])
    counts: dict[str, int] = {}
    for value in types:
        label = canon.get(names.get(int(value), ""), "other")
        counts[label] = counts.get(label, 0) + 1
    hist = " ".join(f"{k}:{counts[k]}" for k in sorted(counts))
    print("status ok")
    print(f"molecules {n}")
    print(f"time_ms {elapsed:.6f}")
    print(f"histogram {hist}")


def write_log(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{k} {v}\n" for k, v in rows), encoding="utf-8")


def summary_rows() -> list[str]:
    if not SUMMARY.is_file():
        return []
    lines = SUMMARY.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    out = []
    for line in lines[1:]:
        cols = line.split("\t")
        if cols[idx["median_ms"]] == "":
            continue
        out.append(",".join([
            cols[idx["comparison"]],
            cols[idx["system"]],
            cols[idx["molecules"]],
            cols[idx["tool"]],
            cols[idx["repeats"]],
            cols[idx["median_ms"]],
            cols[idx["iqr_ms"]],
            '"' + cols[idx["histogram"]].replace('"', "") + '"',
        ]))
    return out


def ensure_frame(name: str) -> Path | None:
    done = subprocess.run(
        ["python3.12", "scripts/h2h_cpc.py", "prepare", name],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if done.stdout.startswith("skip"):
        return None
    if done.returncode != 0:
        raise SystemExit(done.stderr[-2000:] or done.stdout)
    return SYSTEMS / f"{name}.lammpstrj"


def main_worker(argv: list[str]) -> None:
    worker(argv[0], float(argv[1]))


def main() -> None:
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        main_worker(sys.argv[2:])
        return
    if not BENCH.is_file():
        raise SystemExit("missing bbdir/h2h-release/bench_h2h_cpc")
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "environment.txt").write_text(
        "\n".join([
            cpu_model(),
            Path("/proc/meminfo").read_text(encoding="utf-8").splitlines()[0],
            "compiler " + subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
            "threads taskset -c 0 OMP_NUM_THREADS=1 OVITO_THREAD_COUNT=1",
            "warmup 1",
            "repeats 6",
            "cutoff 3.5",
            "rule d-SEAMS ChillPlusRule::Ovito and OVITO ChillPlusModifier",
            "plain CHILL is only in d-SEAMS; OVITO publishes ChillPlusModifier",
            "",
        ]),
        encoding="utf-8",
    )
    chill_lines = []
    slower = []
    for name in FRAMES:
        dump = ensure_frame(name)
        if dump is None:
            print(f"[chill-plus] {name} skipped, liquid seed absent", flush=True)
            continue
        print(f"[chill-plus] {name}", flush=True)
        samples: dict[str, list[float]] = {tool: [] for tool in TOOLS}
        hist: dict[str, str] = {}
        phases: dict[str, list[float]] = {"neigh": [], "class": []}
        calls = [("warmup", i) for i in range(WARMUP)] + [("timed", i) for i in range(REPEATS)]
        for phase, index in calls:
            order = TOOLS if index % 2 == 0 else tuple(reversed(TOOLS))
            for tool in order:
                argv = dseams_argv(dump) if tool == "dseams" else ovito_argv(dump)
                record = run(argv)
                hist[tool] = record.get("histogram", "")
                if phase == "timed":
                    samples[tool].append(float(record["time_ms"]))
                    if tool == "dseams":
                        phases["neigh"].append(float(record.get("neigh_ms", "nan")))
                        phases["class"].append(float(record.get("class_ms", "nan")))
                write_log(RAW / name / tool / f"{phase}-{index:02d}.txt", [
                    ("comparison", "chill-plus"),
                    ("system", name),
                    ("tool", tool),
                    ("phase", phase),
                    ("repeat", str(index)),
                    ("molecules", record.get("molecules", "")),
                    ("time_ms", record.get("time_ms", "")),
                    ("neigh_ms", record.get("neigh_ms", "")),
                    ("class_ms", record.get("class_ms", "")),
                    ("histogram", record.get("histogram", "")),
                    ("command", record["command"]),
                    ("loadavg_before", record["loadavg_before"]),
                    ("loadavg_after", record["loadavg_after"]),
                    ("wall_ms", record["wall_ms"]),
                ])
        if hist["dseams"] != hist["ovito"]:
            print(f"[chill-plus] {name} histogram differs", flush=True)
            print(" dseams", hist["dseams"], flush=True)
            print(" ovito", hist["ovito"], flush=True)
        else:
            print(f"[chill-plus] {name} histogram match", flush=True)
        molecules = Path(RAW / name / "dseams" / "timed-00.txt").read_text(encoding="utf-8")
        mol = parse_kv(molecules)["molecules"]
        medians = {}
        for tool in TOOLS:
            median = type7(samples[tool], 0.5)
            iqr = type7(samples[tool], 0.75) - type7(samples[tool], 0.25)
            medians[tool] = median
            chill_lines.append(",".join([
                "chill-plus", name, mol, tool, str(REPEATS),
                f"{median:.6f}", f"{iqr:.6f}", '"' + hist[tool] + '"',
            ]))
            print(f"[chill-plus] {name} {tool} median_ms {median:.3f} iqr {iqr:.3f}", flush=True)
        if medians["dseams"] > medians["ovito"]:
            slower.append((name, medians, phases))
    rows = ["comparison,system,molecules,tool,n,median_ms,iqr_ms,histogram"]
    rows.extend(summary_rows())
    rows.extend(chill_lines)
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}", flush=True)
    if slower:
        lines = [
            "d-SEAMS CHILL+ is slower than OVITO ChillPlusModifier on the frames below.",
            "time_ms is the cutoff neighbour list plus applyChillPlus(Ovito).",
            "neigh_ms is neighListO. class_ms is applyChillPlus.",
            "",
        ]
        for name, medians, phases in slower:
            lines.append(
                f"{name}: d-SEAMS {medians['dseams']:.3f} ms, OVITO {medians['ovito']:.3f} ms, "
                f"neigh median {type7(phases['neigh'], 0.5):.3f} ms, "
                f"classify median {type7(phases['class'], 0.5):.3f} ms"
            )
        note = ROOT / "h2h-cpc" / "chill-note.txt"
        note.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"wrote {note.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
