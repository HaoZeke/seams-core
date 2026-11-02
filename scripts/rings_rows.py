#!/usr/bin/env python3
"""R.I.N.G.S. rows for the 5+30 protocol.

The Primitive Rings line is read from each completed run. Both series are
one session: WARMUP unrecorded repetitions, then REPS recorded ones, pinned
to CORE. The mW rows in the paper are not this session. They are three
recorded repetitions and no further warmup, one on each of cores 3, 4 and
5, from scripts/rings_n3.py.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import importlib.util

from timing_protocol import CORE, CUT, REPS, WARMUP

_ROOT = Path(__file__).resolve().parents[1]
_H2H = Path(os.environ.get("H2H_ROOT", str(_ROOT / "h2h")))
RINGS = os.environ.get("RINGS", str(_H2H / "rings" / "bin" / "rings"))
SEAMS = os.environ.get("SEAMS", str(_ROOT / "build" / "seams"))
PROTOCOL = os.environ.get("BENCH_PROTOCOL", str(_ROOT / "build" / "tests" / "bench_protocol"))
JSON_PATH = Path(os.environ.get("TIMING_JSON", str(_H2H / "out" / "timings.json")))
DATA = Path(os.environ.get("H2H_DATA", str(_H2H / "data")))


def load_tp():
    spec = importlib.util.spec_from_file_location(
        "tp", Path(__file__).with_name("timing_protocol.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cpp_env():
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_PROC_BIND"] = "true"
    parts = [p for p in env.get("LD_LIBRARY_PATH", "").split(":") if p and not p.startswith("/lib/") and not p.startswith("/usr/lib/")]
    if parts:
        env["LD_LIBRARY_PATH"] = ":".join(parts)
    else:
        env.pop("LD_LIBRARY_PATH", None)
    return env


def kernel(tp_mod, traj, sid, frame):
    cmd = ["taskset", "-c", str(CORE), PROTOCOL, str(traj), "1", str(WARMUP), str(REPS), "rings"]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=cpp_env())
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-400:])
    samples, census = [], ""
    for line in proc.stdout.splitlines():
        if line.startswith("ring_census "):
            census = line.split(" ", 1)[1]
        else:
            try:
                samples.append([float(line)])
            except ValueError:
                pass
    return {
        "id": sid,
        "role": "head-to-head",
        "comparison": "primitive rings D=6, ring stage",
        "scope": "ring-stage",
        "frame": frame,
        "side": "dseams",
        "code_name": "d-SEAMS 2.0",
        "command": cmd,
        "warmup": WARMUP,
        "columns": ["kernel_ms"],
        "samples": samples,
        "counts": {"ring_census": census},
        "agreement": "",
    }


def fingerprint(traj, sid, frame):
    cmd = ["taskset", "-c", str(CORE), SEAMS, "fingerprint", str(traj), "--cutoff", str(CUT), "--type", "1"]
    env = cpp_env()
    samples = []
    text = ""
    for i in range(WARMUP + REPS):
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
        dt = (time.perf_counter() - t0) * 1e3
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-400:])
        text = proc.stdout.strip()
        if i >= WARMUP:
            samples.append([dt])
    return {
        "id": sid,
        "role": "head-to-head",
        "comparison": "primitive rings D=6, whole process",
        "scope": "whole-process",
        "frame": frame,
        "side": "dseams",
        "code_name": "d-SEAMS 2.0 seams fingerprint",
        "command": cmd,
        "warmup": WARMUP,
        "columns": ["wall_ms"],
        "samples": samples,
        "counts": {"stdout": text, "ring_census": text},
        "agreement": "",
    }


def rings_session(tp, work, pos, box, frame):
    tp.prepare_rings(work, pos, box)
    cmd = ["taskset", "-c", str(CORE), RINGS, "input"]
    env = cpp_env()
    walls, stages, census = [], [], ""
    log = Path(os.environ.get("H2H_WORK", str(_H2H / "work"))) / "rings-rows.log"
    for i in range(WARMUP + REPS):
        tp.clean_rings(work)
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, cwd=work, check=False, capture_output=True, text=True, env=env)
        dt = (time.perf_counter() - t0) * 1e3
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-500:])
        stage = None
        for line in proc.stdout.splitlines():
            if line.strip().startswith("Primitive Rings"):
                stage = float(line.split("|")[1]) * 1000.0
        if stage is None:
            raise RuntimeError("no Primitive Rings line")
        census = " ".join(f"{k}:{v}" for k, v in sorted(tp.parse_rings_census(work).items()))
        if i >= WARMUP:
            walls.append([dt])
            stages.append([stage])
        with log.open("a") as fh:
            fh.write(f"{frame} rep {i} wall {dt:.3f} stage {stage:.3f}\n")
    common = {
        "role": "head-to-head",
        "frame": frame,
        "side": "rings",
        "code_name": "R.I.N.G.S. 1.3.5",
        "command": cmd,
        "warmup": WARMUP,
        "counts": {"ring_census": census},
        "agreement": "",
    }
    return [
        {**common, "id": f"rings-stage-{frame}", "comparison": "primitive rings D=6, ring stage", "scope": "ring-stage", "columns": ["report_ms"], "samples": stages},
        {**common, "id": f"rings-wall-{frame}", "comparison": "primitive rings D=6, whole process", "scope": "whole-process", "columns": ["wall_ms"], "samples": walls},
    ]


def merge(new):
    payload = json.loads(JSON_PATH.read_text())
    by_id = {s["id"]: s for s in payload["sessions"]}
    order = [s["id"] for s in payload["sessions"]]
    for session in new:
        if session["id"] not in by_id:
            order.append(session["id"])
        by_id[session["id"]] = session
    payload["sessions"] = [by_id[i] for i in order]
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n")


def main():
    os.sched_setaffinity(0, {CORE})
    tp = load_tp()
    made = []
    frames = [
        ("mW cubic, 4096", _ROOT / "input" / "traj" / "mW_cubic.lammpstrj"),
        ("CRN1, 1000", DATA / "crn1.lammpstrj"),
    ]
    # Fast d-SEAMS sides first, then the long R.I.N.G.S. processes, CRN1 before mW.
    for frame, traj in frames:
        made.append(kernel(tp, traj, f"rings-kernel-{frame}", frame))
        made.append(fingerprint(traj, f"rings-fp-{frame}", frame))
        merge(made)
        print("dseams", frame, flush=True)
    for frame, traj in (("CRN1, 1000", DATA / "crn1.lammpstrj"), ("mW cubic, 4096", _ROOT / "input" / "traj" / "mW_cubic.lammpstrj")):
        pos, box = tp.read_lammps(traj)
        work = Path(os.environ.get("H2H_WORK", str(_H2H / "work"))) / "rings" / frame.split(",")[0].replace(" ", "")
        print("RINGS start", frame, flush=True)
        made.extend(rings_session(tp, work, pos, box, frame))
        merge(made)
        print("RINGS done", frame, flush=True)


if __name__ == "__main__":
    main()
