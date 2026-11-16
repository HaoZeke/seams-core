#!/usr/bin/env python3
"""scripts/pin_run.py pins the width and the cpus and records what ran."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIN_RUN = ROOT / "scripts" / "pin_run.py"


def pinned(cpus: str, width: int, code: str, record: Path) -> subprocess.CompletedProcess:
    env = {key: value for key, value in os.environ.items()
           if key not in ("OMP_NUM_THREADS", "RAYON_NUM_THREADS")}
    return subprocess.run(
        [sys.executable, str(PIN_RUN), "--cpus", cpus, "--width", str(width),
         "--record", str(record), "--period", "0.01", "--", sys.executable, "-c", code],
        env=env, capture_output=True, text=True, timeout=60)


class PinRunTest(unittest.TestCase):
    def setUp(self):
        self.cpu = str(min(os.sched_getaffinity(0)))
        self.record = Path(tempfile.mkdtemp()) / "run.cpus"

    def test_width_and_mask_reach_the_command(self):
        done = pinned(self.cpu, 1, (
            "import os, time\n"
            "print(os.environ['OMP_NUM_THREADS'], os.environ['RAYON_NUM_THREADS'],\n"
            "      ','.join(map(str, sorted(os.sched_getaffinity(0)))))\n"
            "time.sleep(0.2)\n"), self.record)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stdout.split(), ["1", "1", self.cpu])
        line = self.record.read_text().strip()
        self.assertIn(f"cpus={self.cpu} width=1 ", line)
        self.assertIn("cpus_effective=1 ", line)
        self.assertIn("peak_threads=1 ", line)
        self.assertIn("escaped=0 ", line)

    @unittest.skipIf(len(os.sched_getaffinity(0)) < 2, "widening needs a second cpu")
    def test_a_command_that_widens_its_mask_fails(self):
        done = pinned(self.cpu, 1, (
            "import os, time\n"
            f"os.sched_setaffinity(0, {sorted(os.sched_getaffinity(0))!r})\n"
            "time.sleep(0.2)\n"), self.record)
        self.assertEqual(done.returncode, 3, done.stderr)
        self.assertIn("escaped=1 ", self.record.read_text())

    def test_the_command_status_is_kept(self):
        done = pinned(self.cpu, 1, "raise SystemExit(5)", self.record)
        self.assertEqual(done.returncode, 5)

    def test_killing_pin_run_stops_the_command(self):
        pid_file = self.record.with_name("child.pid")
        runner = subprocess.Popen(
            [sys.executable, str(PIN_RUN), "--cpus", self.cpu, "--width", "1", "--",
             sys.executable, "-c",
             f"import os, time\nopen({str(pid_file)!r}, 'w').write(str(os.getpid()))\n"
             "time.sleep(60)\n"],
            stderr=subprocess.DEVNULL)
        deadline = time.monotonic() + 30
        while not pid_file.is_file() or not pid_file.read_text():
            self.assertLess(time.monotonic(), deadline, "the command did not start")
            time.sleep(0.02)
        child = int(pid_file.read_text())
        runner.kill()
        runner.wait()
        while Path(f"/proc/{child}").exists():
            state = Path(f"/proc/{child}/stat").read_text().rsplit(")", 1)[-1].split()[0]
            if state == "Z":
                break
            self.assertLess(time.monotonic(), deadline, "the command outlived pin_run")
            time.sleep(0.02)


if __name__ == "__main__":
    unittest.main()
