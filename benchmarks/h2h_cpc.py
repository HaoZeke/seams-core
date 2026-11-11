# Copyright (c) d-SEAMS developers
"""Airspeed Velocity arms for the neighbour and primitive-ring comparison.

Each method is one serial call of the same worker the Snakemake rules use.
The interleaved medians are the campaign logs, not these single calls.
freud publishes no primitive-ring count, so the ring class has no freud arm.
"""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OVITO_THREAD_COUNT": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONPATH": ".h2h-py",
})
SYSTEMS = [
    "ice-1k", "ice-10k", "ice-50k", "ice-110k",
    "liq-1k", "liq-10k", "liq-50k", "liq-110k",
    "interface-50k",
]


def _prepare(system: str) -> None:
    done = subprocess.run(
        ["python3.12", "scripts/h2h_cpc.py", "prepare", system],
        cwd=ROOT, env=ENV, check=False, text=True, stdout=subprocess.PIPE,
    )
    if done.stdout.startswith("skip"):
        raise NotImplementedError(done.stdout.strip())
    if done.returncode != 0:
        raise RuntimeError(done.stderr or done.stdout)


def _run(comparison: str, tool: str, system: str) -> None:
    cmd = [
        "taskset", "-c", "0",
        "python3.12", "scripts/h2h_cpc.py", "asv-call",
        comparison, tool, system,
    ]
    subprocess.run(cmd, cwd=ROOT, env=ENV, check=True)


class NeighbourLists:
    params = [SYSTEMS]
    param_names = ["system"]

    def setup(self, system):
        _prepare(system)

    def time_dseams(self, system):
        _run("neighbours", "dseams", system)

    def time_freud(self, system):
        _run("neighbours", "freud", system)

    def time_ovito(self, system):
        _run("neighbours", "ovito", system)


class PrimitiveRings:
    params = [SYSTEMS]
    param_names = ["system"]

    def setup(self, system):
        _prepare(system)

    def time_dseams(self, system):
        _run("rings", "dseams", system)

    def time_ovito(self, system):
        _run("rings", "ovito", system)


for _method in (
    NeighbourLists.time_dseams,
    NeighbourLists.time_freud,
    NeighbourLists.time_ovito,
    PrimitiveRings.time_dseams,
    PrimitiveRings.time_ovito,
):
    _method.repeat = 1
    _method.number = 1
    _method.warmup_time = 0
