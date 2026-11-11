# Copyright (c) d-SEAMS developers
"""Shared helpers for the serial results measurements.

Paths written to logs are relative to the repository root. Absolute
locations under the agent, home, workspace, cursor or tmp trees are stripped.
"""

from __future__ import annotations

import math
import os
import re
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECORDS = ROOT / "benchmarks" / "cpc" / "records"
BIN = ROOT / "build" / "benchmarks" / "cpc" / "measure_cpc"
SEAMS = ROOT / "build" / "seams"
PY = ROOT / "benchmarks" / "cpc" / ".venv" / "bin" / "python"
WARMUP = 1
REPEATS = 6
CORE = "0"
CUTOFF = 3.5
FORBIDDEN = re.compile(r"/(?:home|workspace|cursor|tmp|agent)(?:/[^\s\"']*)?")

ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "OVITO_THREAD_COUNT": "1",
    "QT_QPA_PLATFORM": "offscreen",
})
ENV["NO_COLOR"] = "1"
ENV["TERM"] = "dumb"
ENV.pop("SEAMS_ROOT", None)


def scrub(text: str) -> str:
    return FORBIDDEN.sub("[path]", text)


def loadavg() -> str:
    return Path("/proc/loadavg").read_text().strip()


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    count = len(ordered)
    if count == 0:
        raise ValueError("empty sample")
    if count == 1:
        return ordered[0]
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def summary(values: list[float]) -> dict[str, float]:
    med = type7(values, 0.5)
    q1 = type7(values, 0.25)
    q3 = type7(values, 0.75)
    return {
        "n": float(len(values)),
        "median": med,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
    }


def run(cmd: list[str], log: Path, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    log.parent.mkdir(parents=True, exist_ok=True)
    full = ["taskset", "-c", CORE, *[str(part) for part in cmd]]
    shown = " ".join(full)
    proc = subprocess.run(
        full,
        cwd=cwd or ROOT,
        env=ENV,
        text=True,
        capture_output=True,
        check=False,
    )
    body = [
        f"command {shown}",
        f"loadavg {loadavg()}",
        f"returncode {proc.returncode}",
        "stdout",
        proc.stdout,
        "stderr",
        proc.stderr,
    ]
    log.write_text(scrub("\n".join(body)))
    if proc.returncode != 0:
        raise RuntimeError(f"command failed, see {log.name}")
    return proc


def parse_repeats(stdout: str) -> dict[str, list[float]]:
    series: dict[str, list[float]] = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] == "repeat":
            key = parts[2]
            series.setdefault(key, []).append(float(parts[3]))
    return series


def parse_counts(stdout: str) -> dict[str, str]:
    counts = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == "count":
            counts[parts[1]] = parts[2]
    return counts


def rel(path: Path) -> str:
    return os.path.relpath(path, ROOT)


NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)

MERGED = {
    "cubic": "cubic",
    "reCubic": "cubic",
    "hexagonal": "hexagonal",
    "reHex": "hexagonal",
    "water": "other",
    "unclassified": "other",
    "interfacial": "interfacial",
    "clathrate": "hydrate",
    "interClathrate": "interfacial-hydrate",
}

OVITO_MERGED = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial-hydrate",
}
