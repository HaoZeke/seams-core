#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Identity of neighbours, six-rings and CHILL+ across one mW quench window.

The window dump is read once. readcon writes a multi-frame .con file that
keeps the LAMMPS atom id, the element and the cell on every frame. The same
frames are written as XYZ with a different atom order on each frame. seams
then builds the 3.5 A cutoff list, the primitive rings through seven members,
and the published CHILL+ label. Persistence is the fraction of that record
at frame t that is still there at t+tau when the two frames are joined by
atom id. For the XYZ file that id is the order of the lines.

Usage:
  python3 scripts/readcon_connectivity.py DUMP --seams-build BUILD
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCAN_SRC = HERE / "readcon_connectivity_scan.cpp"
LAGS = (1, 5, 10, 25)
REPEATS = 8


def features_ok(seams: Path) -> bool:
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    text = subprocess.run(
        [str(seams), "--features"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout
    return "readcon-core: enabled" in text and "readcon-db: enabled" in text


def compile_flags(build: Path) -> list[str]:
    commands = json.loads((build / "compile_commands.json").read_text())
    raw = ""
    for entry in commands:
        if str(entry.get("file", "")).endswith("seams_cli.cpp"):
            raw = entry["command"]
            break
    if not raw:
        raise SystemExit("compile_commands.json has no seams_cli.cpp")
    args = shlex.split(raw)
    flags: list[str] = []
    i = 1
    while i < len(args):
        arg = args[i]
        if arg in ("-c", "-MD"):
            i += 1
            continue
        if arg in ("-o", "-MF", "-MQ"):
            i += 2
            continue
        if arg.endswith(".cpp"):
            i += 1
            continue
        flags.append(arg)
        i += 1
    return [args[0], *flags]


def link_flags(build: Path) -> list[str]:
    text = (build / "build.ninja").read_text()
    start = text.find("\nbuild seams: cpp_LINKER")
    if start < 0:
        raise SystemExit("build.ninja has no seams link line")
    end = text.find("\nbuild ", start + 1)
    block = text[start:end if end > 0 else None]
    for line in block.splitlines():
        if "LINK_ARGS" in line and "=" in line:
            value = line.split("=", 1)[1].strip().replace("$$", "$")
            return shlex.split(value)
    raise SystemExit("seams LINK_ARGS missing")


def build_scan(build: Path, work: Path) -> Path:
    binary = work / "scan"
    compiler = compile_flags(build)
    obj = work / "scan.o"
    subprocess.run(
        [*compiler, "-c", str(SCAN_SRC), "-o", str(obj)],
        cwd=build,
        check=True,
    )
    rpath = [
        build / "src",
        build / "subprojects" / "readcon-core",
        build / "subprojects" / "readcon-db",
    ]
    subprocess.run(
        [
            compiler[0],
            str(obj),
            "-o",
            str(binary),
            *[f"-Wl,-rpath,{path}" for path in rpath],
            *link_flags(build),
        ],
        cwd=build,
        check=True,
    )
    return binary


def run_text(cmd: list[str], env: dict[str, str]) -> str:
    done = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if done.stderr:
        sys.stderr.write(done.stderr)
    if done.returncode != 0:
        raise SystemExit(done.returncode)
    return done.stdout


def parse_emit(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        key, value = line.split(None, 1)
        out[key] = value
    return out


def load_jsonl(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        if line:
            rows.append(json.loads(line))
    return rows


def mean_at(frames: list[dict], lag: int, kind: str) -> float:
    hits = 0.0
    used = 0
    for t in range(len(frames) - lag):
        left = frames[t]
        right = frames[t + lag]
        if kind == "neighbour":
            pairs = left["pairs"]
            now = {(pairs[i], pairs[i + 1]) for i in range(0, len(pairs), 2)}
            later_pairs = right["pairs"]
            later = {
                (later_pairs[i], later_pairs[i + 1]) for i in range(0, len(later_pairs), 2)
            }
            if not now:
                continue
            hits += len(now & later) / len(now)
        elif kind == "ring6":
            now = set(left["ring6"])
            later = set(right["ring6"])
            if not now:
                continue
            hits += len(now & later) / len(now)
        else:
            later = dict(zip(right["ice_id"], right["ice_code"]))
            same = 0
            total = 0
            for atom, code in zip(left["ice_id"], left["ice_code"]):
                if atom in later:
                    total += 1
                    same += int(later[atom] == code)
            if total == 0:
                continue
            hits += same / total
        used += 1
    if used == 0:
        raise SystemExit(f"no frame pairs at lag {lag}")
    return hits / used


def aggregate_mismatch(con: list[dict], xyz: list[dict]) -> tuple[int, int]:
    if len(con) != len(xyz):
        return len(con) + len(xyz), len(con) + len(xyz)
    rings = 0
    ice = 0
    for a, b in zip(con, xyz):
        if a["nop"] != b["nop"]:
            rings += abs(a["nop"] - b["nop"])
        for x, y in zip(a["rings"], b["rings"]):
            rings += abs(x - y)
        for x, y in zip(a["ice"], b["ice"]):
            ice += abs(x - y)
    return rings, ice


def draw(con: list[dict], xyz: list[dict], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    last = min(40, len(con) - 1)
    lags = list(range(1, last + 1))
    panels = (
        ("neighbour", "Neighbour pairs"),
        ("ring6", "Atoms in a 6-ring"),
        ("ice", "CHILL+ ice type"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), sharey=True)
    for ax, (kind, title) in zip(axes, panels):
        ax.plot(
            lags,
            [mean_at(con, lag, kind) for lag in lags],
            color="#1f4e79",
            lw=1.8,
            label=".con",
        )
        ax.plot(
            lags,
            [mean_at(xyz, lag, kind) for lag in lags],
            color="#c46b2c",
            lw=1.8,
            ls="--",
            label="shuffled XYZ",
        )
        ax.set_title(title)
        ax.set_xlabel("lag (frames)")
        ax.set_xlim(1, last)
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, lw=0.4, alpha=0.5)
    axes[0].set_ylabel("fraction still true at t+tau")
    axes[0].legend(frameon=False, loc="center")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def series(frames: list[dict], kind: str) -> list[float]:
    return [mean_at(frames, lag, kind) for lag in LAGS]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump", type=Path)
    parser.add_argument("--seams-build", type=Path, default=None)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "bench" / "readcon-connectivity",
    )
    args = parser.parse_args()
    build = args.seams_build
    if build is None:
        env_build = os.environ.get("SEAMS_BUILD", "")
        if not env_build:
            raise SystemExit("pass --seams-build or set SEAMS_BUILD")
        build = Path(env_build)
    seams = build / "seams"
    if not features_ok(seams):
        raise SystemExit("seams --features must show readcon-core and readcon-db enabled")

    out_dir = args.out_dir
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    cores = os.cpu_count() or 1
    env["OMP_NUM_THREADS"] = str(cores)

    scan = build_scan(build, work)
    con_path = work / "window.con"
    xyz_path = work / "window.xyz"
    emitted = parse_emit(
        run_text([str(scan), "emit", str(args.dump), str(con_path), str(xyz_path)], env)
    )
    con_rows = work / "con.jsonl"
    xyz_rows = work / "xyz.jsonl"
    run_text([str(scan), "scan-con", str(con_path), str(con_rows)], env)
    run_text([str(scan), "scan-xyz", str(xyz_path), str(xyz_rows)], env)
    con = load_jsonl(con_rows)
    xyz = load_jsonl(xyz_rows)
    ring_mismatch, ice_mismatch = aggregate_mismatch(con, xyz)
    trip = json.loads(
        run_text(
            [str(scan), "roundtrip", str(con_path), str(work / "corpus"), str(REPEATS)],
            env,
        )
    )
    melt = int(emitted["melt_steps"])
    step_fs = float(emitted["step_fs"])
    first = int(emitted["first_step"])
    last = int(emitted["last_step"])
    ps = lambda step: (step - melt) * step_fs / 1000.0
    result = {
        "seed": 20261118,
        "rng_seed": int(emitted["rng_seed"]),
        "nop": int(emitted["nop"]),
        "frames": int(emitted["frames"]),
        "first_step": first,
        "last_step": last,
        "ps_first": ps(first),
        "ps_last": ps(last),
        "cutoff_A": float(emitted["cutoff"]),
        "ring_depth": int(emitted["ring_depth"]),
        "n_cores": int(trip["n_cores"]),
        "dump_regenerated": 1,
        "lags": list(LAGS),
        "con_neighbour": series(con, "neighbour"),
        "con_ring6": series(con, "ring6"),
        "con_ice": series(con, "ice"),
        "xyz_neighbour": series(xyz, "neighbour"),
        "xyz_ring6": series(xyz, "ring6"),
        "xyz_ice": series(xyz, "ice"),
        "aggregate_ring_mismatch": ring_mismatch,
        "aggregate_ice_mismatch": ice_mismatch,
        "roundtrip_frames": int(trip["frames"]),
        "roundtrip_exact": int(trip["exact"]),
        "roundtrip_repeats": int(trip["repeats"]),
        "roundtrip_max_abs_position": float(trip["max_abs_position"]),
        "roundtrip_max_abs_cell": float(trip["max_abs_cell"]),
        "roundtrip_seconds": [float(v) for v in trip["seconds"]],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    draw(con, xyz, out_dir / "persistence.png")
    print(json.dumps(result, indent=2))
    return 0 if ring_mismatch == 0 and ice_mismatch == 0 and trip["exact"] == 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
