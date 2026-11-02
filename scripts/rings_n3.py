#!/usr/bin/env python3
"""The mW rows of the paper: three recorded repetitions, one per core.

No warmup is discarded. One earlier warmup on core 3 is not a sample.
d-SEAMS uses cores 3, 4 and 5 and finishes before R.I.N.G.S. starts, so
the two sides do not share the machine. scripts/rings_rows.py is the
5+30 protocol and does not produce these rows.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import threading
import time
from pathlib import Path

from timing_protocol import CUT

# The recorded n=3 set. This is not CORE from timing_protocol.
CORES = [3, 4, 5]
_ROOT = Path(__file__).resolve().parents[1]
H2H = Path(os.environ.get("H2H_ROOT", str(_ROOT / "h2h")))
RINGS = os.environ.get("RINGS", str(H2H / "rings" / "bin" / "rings"))
SEAMS = os.environ.get("SEAMS", str(_ROOT / "build" / "seams"))
TRAJ = os.environ.get("H2H_TRAJ", str(_ROOT / "input" / "traj" / "mW_cubic.lammpstrj"))
OUT = _ROOT / "logs" / "rings-n3-results.json"
LOG = _ROOT / "logs" / "rings-n3.log"


def load_tp():
    spec = importlib.util.spec_from_file_location(
        "tp", Path(__file__).with_name("timing_protocol.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def child_env():
    e = os.environ.copy()
    e["OMP_NUM_THREADS"] = "1"
    e["OMP_PROC_BIND"] = "true"
    e.pop("LD_LIBRARY_PATH", None)
    return e


def log(msg: str) -> None:
    # The stamp is the guest wall clock. The durations are perf_counter.
    LOG.parent.mkdir(parents=True, exist_ok=True)
    line = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + msg
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def save(payload) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(OUT)


def main() -> int:
    os.sched_setaffinity(0, set(range(8)))
    tp = load_tp()
    pos, box = tp.read_lammps(Path(TRAJ))
    payload = {"cores": CORES, "dseams": [], "rings": [], "errors": []}
    lock = threading.Lock()

    def guard(fn, core: int) -> None:
        try:
            fn(core)
        except Exception as err:
            with lock:
                payload["errors"].append(f"core {core}: {err}")
                save(payload)
            log(f"ERROR core {core} {err}")

    def dseams(core: int) -> None:
        cmd = [
            "taskset", "-c", str(core), SEAMS, "fingerprint", TRAJ,
            "--cutoff", str(CUT), "--type", "1",
        ]
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=child_env())
        dt = (time.perf_counter() - t0) * 1e3
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-400:])
        row = {"core": core, "wall_ms": dt, "stdout": proc.stdout.strip(), "command": cmd}
        with lock:
            payload["dseams"].append(row)
            save(payload)
        log(f"dseams core {core} wall {dt:.3f}")

    threads = [threading.Thread(target=guard, args=(dseams, core)) for core in CORES]
    log("dseams start")
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if payload["errors"] or len(payload["dseams"]) != 3:
        raise SystemExit(f"dseams failed {payload['errors']}")

    works = {}
    for core in CORES:
        work = H2H / "work" / "rings-n3" / f"core{core}"
        tp.prepare_rings(work, pos, box)
        works[core] = work

    def rings(core: int) -> None:
        work = works[core]
        cmd = ["taskset", "-c", str(core), RINGS, "input"]
        stdout_path = work / "stdout.txt"
        stderr_path = work / "stderr.txt"
        tp.clean_rings(work)
        with stdout_path.open("w") as out, stderr_path.open("w") as err:
            t0 = time.perf_counter()
            proc = subprocess.Popen(cmd, cwd=work, stdout=out, stderr=err, env=child_env())
            rc = proc.wait()
            dt = (time.perf_counter() - t0) * 1e3
        if rc != 0:
            raise RuntimeError(stderr_path.read_text(errors="replace")[-500:])
        text = stdout_path.read_text(errors="replace")
        stage = None
        for line in text.splitlines():
            if line.strip().startswith("Primitive Rings"):
                stage = float(line.split("|")[1]) * 1000.0
        if stage is None:
            raise RuntimeError(f"core {core} has no Primitive Rings line")
        census = " ".join(f"{k}:{v}" for k, v in sorted(tp.parse_rings_census(work).items()))
        row = {
            "core": core,
            "wall_ms": dt,
            "stage_ms": stage,
            "ring_census": census,
            "command": cmd,
        }
        with lock:
            payload["rings"].append(row)
            save(payload)
        log(f"rings core {core} wall {dt:.3f} stage {stage:.3f} census {census}")

    log("rings start")
    threads = [threading.Thread(target=guard, args=(rings, core)) for core in CORES]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if len(payload["rings"]) != 3 or len(payload["dseams"]) != 3:
        raise SystemExit(f"incomplete {payload}")
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
