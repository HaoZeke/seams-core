#!/usr/bin/env python3
"""Turn the 600 s DNF rows into measured outcomes.

Each comparator call is allowed 3600 s. A size that returns is repeated
eight times after one warm-up, in one session, alternating with d-SEAMS.
A call that does not return is not repeated. Peak RSS is the process
high-water mark. Serial tools have no thread argument; the eight-thread
row is the same call pinned to cores 0-7.
"""

from __future__ import annotations

import argparse
import os
import resource
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
LOG = Path(os.environ.get("H2H_DNF_LOG", str(ROOT / "benchmarks" / "logs" / "h2h-dnf-hour.txt")))
BUDGET_S = int(os.environ.get("H2H_BUDGET_S", "3600"))
SAMPLES = 8


class CallTimeout(Exception):
    pass


def _alarm(_signum, _frame):
    raise CallTimeout(f"exceeded {BUDGET_S} s")


def cpu_model() -> str:
    text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def pin(threads: int) -> str:
    forced = os.environ.get("H2H_PIN_CPU", "")
    if forced != "":
        cpu = int(forced)
        os.sched_setaffinity(0, {cpu})
        return str(cpu)
    group = os.environ.get("H2H_PIN_CPUS", "")
    if group != "":
        cpus = {int(part) for part in group.split(",") if part != ""}
        os.sched_setaffinity(0, cpus)
        ordered = sorted(cpus)
        if len(ordered) == 1:
            return str(ordered[0])
        return ",".join(str(cpu) for cpu in ordered)
    os.sched_setaffinity(0, {0} if threads == 1 else set(range(8)))
    return "0" if threads == 1 else "0-7"


def rss_kb() -> int:
    # Linux ru_maxrss is kilobytes.
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def append(line: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="ascii") as handle:
        handle.write(line.rstrip("\n") + "\n")
        handle.flush()


def _deadline_child() -> str:
    # A C extension that keeps the interpreter cannot see SIGALRM.
    # This child sleeps, then records the stop and kills the parent.
    return (
        "import os, sys, time, signal\n"
        "pid = int(sys.argv[1]); seconds = float(sys.argv[2]); path = sys.argv[3]; line = sys.argv[4]\n"
        "time.sleep(seconds)\n"
        "try:\n"
        "    os.kill(pid, 0)\n"
        "except OSError:\n"
        "    raise SystemExit(0)\n"
        "rss = '0'\n"
        "try:\n"
        "    for row in open(f'/proc/{pid}/status', encoding='ascii', errors='replace'):\n"
        "        if row.startswith('VmHWM:'):\n"
        "            rss = row.split()[1]\n"
        "except OSError:\n"
        "    raise SystemExit(0)\n"
        "with open(path, 'a', encoding='ascii') as handle:\n"
        "    handle.write(line.replace('RSS', rss) + '\\n')\n"
        "    handle.flush()\n"
        "os.kill(pid, signal.SIGKILL)\n"
    )


def arm_deadline(line: str) -> subprocess.Popen:
    return subprocess.Popen(
        [sys.executable, "-c", _deadline_child(), str(os.getpid()), str(BUDGET_S), str(LOG), line],
        start_new_session=True,
    )


def disarm(proc: subprocess.Popen | None) -> None:
    if proc is None or proc.poll() is not None:
        return
    proc.kill()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()


def guarded(fn, stop_line: str):
    killer = arm_deadline(stop_line)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(BUDGET_S)
    t0 = time.perf_counter()
    try:
        value = fn()
        elapsed = time.perf_counter() - t0
        return {"ok": True, "seconds": elapsed, "value": value, "rss_kb": rss_kb(), "error": ""}
    except CallTimeout as exc:
        return {"ok": False, "seconds": float(BUDGET_S), "value": "", "rss_kb": rss_kb(), "error": str(exc)}
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        text = f"{type(exc).__name__}: {exc}"
        return {"ok": False, "seconds": elapsed, "value": "", "rss_kb": rss_kb(), "error": text[:500]}
    finally:
        signal.alarm(0)
        disarm(killer)


def ice(molecules: int):
    from h2h_mdapy_pyscal3_vesin import ice_positions

    for count in (2, 4, 8, 24, 48):
        if count ** 3 * 8 == molecules:
            return ice_positions(count)
    raise SystemExit(f"no ice cell has {molecules} molecules")


def prepare_ira(molecules: int):
    import numpy as np
    import ira_mod
    from gen_ice_ic import ice_ic
    from h2h_ira import classify_serial, cluster

    pos, lengths = ice(molecules)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.asarray(lengths, dtype=np.float64)
    engine = ira_mod.IRA()
    ih = cluster(pos, box, 0)
    ic_pos, side = ice_ic(4)
    ic = cluster(np.ascontiguousarray(ic_pos), np.asarray((side, side, side)), 0)

    def call():
        labels = classify_serial(engine, ih, ic, pos, box)
        return int((labels == 0).sum())

    return call


def prepare_cycless(molecules: int):
    import networkx as nx
    import numpy as np
    import freud
    from cycless.cycles import cycles_iter

    pos, lengths = ice(molecules)
    pos = np.ascontiguousarray(pos)
    box = freud.box.Box(*lengths)
    t0 = time.perf_counter()
    query = freud.locality.AABBQuery(box, pos).query(
        pos, {"r_max": 3.5, "exclude_ii": True}
    ).toNeighborList()
    graph = nx.Graph()
    graph.add_nodes_from(range(molecules))
    ii = np.asarray(query[:, 0])
    jj = np.asarray(query[:, 1])
    graph.add_edges_from((int(a), int(b)) for a, b in zip(ii, jj) if int(a) < int(b))
    built = time.perf_counter() - t0

    def call():
        return sum(1 for _ in cycles_iter(graph, 6))

    call.graph_s = built
    return call


def prepare_sqq(molecules: int, threads: int):
    import numpy as np
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame

    pos, lengths = ice(molecules)
    pos_nm = np.asarray(pos, dtype=np.float64) / 10.0
    box = np.asarray(lengths, dtype=np.float64) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(molecules)
    ]
    frame = Frame(name=f"ice-{molecules}", atoms=atoms, box=box)
    # The python shortest-path engine does not return on these cells.
    # SQQ 0.5.8's cpp engine accepts only the chordless definition.
    cfg = load_config({
        "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
        "ring": {"sizes": [4, 5, 6], "definition": "chordless", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": True},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": threads, "math_threads": 1},
    }, engine="cpp")

    def call():
        result = analyze_frame(frame, cfg, engine="cpp")
        return len(result.cages)

    return call


def prepare_dseams(kind: str, molecules: int, threads: int):
    import ctypes
    import h2h_mdapy_pyscal3_vesin as h

    lib = h.load_seams(h.library_path())
    pos, lengths = ice(molecules)
    held = h.frame_system(lib, pos, lengths)
    if kind == "ira":
        # A float argument without a prototype raises ArgumentError and
        # the call never enters the library. The 4096-molecule summary
        # sets this prototype in bind().
        from h2h_summary_extra import bind

        bind(lib)
        n = molecules
        sym = (ctypes.c_int * n)()

        def call():
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")
            return int(sum(1 for value in sym if int(value) == 0))

        return call
    if kind == "cycless":
        from h2h_summary_extra import bind, ring_hist

        bind(lib)

        def call():
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            hist = ring_hist(lib, held, h)
            return int(hist[6])

        return call
    from h2h_summary_extra import dseams_cage_total

    def call():
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        return dseams_cage_total(lib, held, h)

    return call


def session(kind: str, molecules: int, threads: int) -> None:
    cpus = pin(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    cores = os.cpu_count() or 1
    model = cpu_model()
    print(f"session {kind} {molecules} threads {threads}", file=sys.stderr, flush=True)
    if kind == "ira":
        other = prepare_ira(molecules)
        note = "classify_serial has no thread argument"
    elif kind == "cycless":
        other = prepare_cycless(molecules)
        note = f"cycles_iter has no thread argument; graph_build_s {other.graph_s:.3f}"
    elif kind == "sqq":
        other = prepare_sqq(molecules, threads)
        note = "SQQ 0.5.8 engine cpp; ring.definition chordless. The python shortest_path engine does not return on these cells"
    else:
        raise SystemExit(kind)
    ours = prepare_dseams(kind, molecules, threads)
    # One warm-up. A warm-up that misses the budget ends the session.
    # RSS is filled by the outside killer if the call never returns to Python.
    warm_stop = (
        f"stop {kind} molecules {molecules} threads {threads} n 0 "
        f"reason exceeded {BUDGET_S} s rss_kb RSS {note}"
    )
    warm = guarded(other, warm_stop)
    append(
        f"warm {kind} molecules {molecules} threads {threads} cores {cores} cpus {cpus} "
        f"ok {warm['ok']} seconds {warm['seconds']:.3f} rss_kb {warm['rss_kb']} "
        f"value {warm['value']} error {warm['error']} cpu {model}"
    )
    if not warm["ok"]:
        append(
            f"stop {kind} molecules {molecules} threads {threads} n 0 "
            f"reason {warm['error'] or 'time'} rss_kb {warm['rss_kb']} {note}"
        )
        return
    left, right = [], []
    for index in range(SAMPLES):
        order = (("d", ours), ("o", other)) if index % 2 == 0 else (("o", other), ("d", ours))
        got = {}
        for name, fn in order:
            sample_stop = (
                f"stop {kind} molecules {molecules} threads {threads} n {index} "
                f"reason exceeded {BUDGET_S} s rss_kb RSS {note}"
            )
            got[name] = guarded(fn, sample_stop)
            if not got[name]["ok"]:
                append(
                    f"call {kind} molecules {molecules} threads {threads} index {index} "
                    f"side {name} ok 0 seconds {got[name]['seconds']:.3f} "
                    f"rss_kb {got[name]['rss_kb']} error {got[name]['error']}"
                )
                append(
                    f"stop {kind} molecules {molecules} threads {threads} n {index} "
                    f"reason {got[name]['error'] or 'time'} rss_kb {rss_kb()} {note}"
                )
                return
        left.append(got["d"]["seconds"] * 1000.0)
        right.append(got["o"]["seconds"] * 1000.0)
        append(
            f"call {kind} molecules {molecules} threads {threads} index {index} "
            f"dseams_ms {left[-1]:.4f} other_ms {right[-1]:.4f} "
            f"dseams_value {got['d']['value']} other_value {got['o']['value']} "
            f"rss_kb {rss_kb()}"
        )
        print(
            f"call {kind} {molecules} t {threads} {index} other_ms {right[-1]:.1f}",
            file=sys.stderr, flush=True,
        )
    append(
        f"done {kind} molecules {molecules} threads {threads} cores {cores} cpus {cpus} "
        f"n {len(right)} rss_kb {rss_kb()} {note} cpu {model} "
        f"dseams_ms {' '.join(f'{v:.4f}' for v in left)} "
        f"other_ms {' '.join(f'{v:.4f}' for v in right)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("ira", "cycless", "sqq"))
    parser.add_argument("molecules", type=int)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    session(args.kind, args.molecules, args.threads)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        append("crash " + traceback.format_exc().replace("\n", " | ")[:800])
        raise
