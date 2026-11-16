#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run one benchmark command on fixed cpus at a fixed thread width.

    python3 scripts/pin_run.py --cpus 0 --width 1 --record OUT.cpus -- CMD ...

The command starts with OMP_NUM_THREADS and RAYON_NUM_THREADS at --width
and its cpu mask set to --cpus before it executes, as taskset -c does.
While it runs, this process watches the command and every process it
starts from cpus outside --cpus when there are any: the cpus each thread
may run on, the cpus they ran on, and the largest thread count of one
process. When the command exits one line is appended to --record (and
written to stderr):

    cpus=0 width=1 cpus_effective=1 allowed=0 ran=0 peak_threads=1 ...

cpus_effective is the size of the union of the masks seen. A mask that
leaves --cpus makes the exit status 3 unless the command already failed:
the run did not measure the width it names.
"""

from __future__ import annotations

import argparse
import ctypes
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import paper_stats as stats  # noqa: E402

ESCAPED = 3
PR_SET_PDEATHSIG = 1
LIBC = ctypes.CDLL(None, use_errno=True)


def text_of(cpus: set[int]) -> str:
    return ",".join(str(cpu) for cpu in sorted(cpus)) or "none"


def children_of(root: int) -> set[int]:
    parent = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
        except OSError:
            continue
        parent[int(entry.name)] = int(stat[stat.rfind(")") + 2:].split()[1])
    found = {root}
    grew = True
    while grew:
        grew = False
        for pid, ppid in parent.items():
            if ppid in found and pid not in found:
                found.add(pid)
                grew = True
    return found


class Watch:
    def __init__(self, root: int, period: float):
        self.root = root
        self.period = period
        self.allowed: set[int] = set()
        self.ran: set[int] = set()
        self.peak = 0
        self.processes: set[int] = set()
        self.done = threading.Event()
        self.thread = threading.Thread(target=self.run, daemon=True)

    def sample(self) -> None:
        for pid in children_of(self.root):
            try:
                tids = os.listdir(f"/proc/{pid}/task")
            except OSError:
                continue
            self.processes.add(pid)
            self.peak = max(self.peak, len(tids))
            for tid in tids:
                base = f"/proc/{pid}/task/{tid}"
                try:
                    status = Path(base, "status").read_text()
                    stat = Path(base, "stat").read_text()
                except OSError:
                    continue
                for line in status.splitlines():
                    if line.startswith("Cpus_allowed_list:"):
                        self.allowed |= stats.cpu_set(line.split(":", 1)[1])
                fields = stat[stat.rfind(")") + 2:].split()
                if len(fields) > 36:
                    self.ran.add(int(fields[36]))

    def run(self) -> None:
        while not self.done.is_set():
            self.sample()
            self.done.wait(self.period)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpus", required=True, help="taskset cpu list, such as 0 or 0-7")
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--record", type=Path, help="file the record line is appended to")
    parser.add_argument("--period", type=float, default=0.05, help="seconds between samples")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.width < 1:
        parser.error("need --width of at least 1 and a command after --")
    cpus = stats.cpu_set(args.cpus)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.width)
    env["RAYON_NUM_THREADS"] = str(args.width)
    own = os.sched_getaffinity(0)
    if not cpus <= own:
        parser.error(f"--cpus {args.cpus} is not inside this process's cpus {text_of(own)}")
    if own - cpus:
        os.sched_setaffinity(0, own - cpus)

    def start() -> None:
        os.sched_setaffinity(0, cpus)
        # A caller that kills this process on a timeout must not leave the
        # command running on cpus the next run will time.
        LIBC.prctl(PR_SET_PDEATHSIG, signal.SIGKILL)

    child = subprocess.Popen(command, env=env, preexec_fn=start)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, lambda number, _frame: child.send_signal(number))
    watch = Watch(child.pid, args.period)
    watch.thread.start()
    status = child.wait()
    watch.done.set()
    watch.thread.join()
    allowed = watch.allowed or cpus
    escaped = not allowed <= cpus or not watch.ran <= cpus
    record = (
        f"cpus={args.cpus} width={args.width} OMP_NUM_THREADS={args.width} "
        f"RAYON_NUM_THREADS={args.width} cpus_effective={len(allowed)} "
        f"allowed={text_of(allowed)} ran={text_of(watch.ran)} peak_threads={watch.peak} "
        f"processes={len(watch.processes)} watcher_cpus={text_of(os.sched_getaffinity(0))} "
        f"escaped={int(escaped)} exit={status} "
        f"time={time.strftime('%Y-%m-%dT%H:%M:%S')} cmd={shlex.join(command)}"
    )
    sys.stderr.write(record + "\n")
    if args.record is not None:
        with args.record.open("a", encoding="utf-8") as handle:
            handle.write(record + "\n")
    if status < 0:
        return 128 - status
    if status != 0:
        return status
    return ESCAPED if escaped else 0


if __name__ == "__main__":
    sys.exit(main())
