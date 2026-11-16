#!/usr/bin/env python3
"""The literature sitting keeps one code runnable and does not pin rivals."""

from __future__ import annotations

import csv
import os
import signal
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import lit_one_sitting as lit  # noqa: E402


class LiteratureSittingTest(unittest.TestCase):
    def test_rival_env_does_not_bind_or_grow_a_blas_pool(self):
        env = lit.child_env("OVITO", 8, Path("/tmp/lib"), Path("/tmp/lang"))
        self.assertNotIn("OMP_PROC_BIND", env)
        self.assertNotIn("OMP_PLACES", env)
        self.assertEqual(env["OVITO_THREAD_COUNT"], "8")
        self.assertEqual(env["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual(env["MKL_NUM_THREADS"], "1")
        self.assertEqual(env["OMP_NUM_THREADS"], "8")
        self.assertEqual(env["QT_QPA_PLATFORM"], "offscreen")

    def test_dseams_env_matches_the_scaling_pin_and_stays_active(self):
        env = lit.child_env("d-SEAMS", 8, Path("/tmp/lib"), Path("/tmp/lang"))
        self.assertEqual(env["OMP_PROC_BIND"], "close")
        self.assertEqual(env["OMP_PLACES"], "cores")
        self.assertEqual(env["OMP_WAIT_POLICY"], "ACTIVE")
        self.assertEqual(env["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual(env["OMP_NUM_THREADS"], "8")
        self.assertNotIn("OVITO_THREAD_COUNT", env)

    def test_plan_covers_every_requested_code_in_its_own_key(self):
        rows = lit.jobs((1000, 10648, 110592))
        seen = {(row["code"], row["operation"], row["threads"], row["molecules"]) for row in rows}
        for molecules in (1000, 10648, 110592):
            for threads in (1, 8):
                self.assertIn(("d-SEAMS", "q6", threads, molecules), seen)
                self.assertIn(("freud", "q6", threads, molecules), seen)
                self.assertIn(("OVITO", "chill_plus", threads, molecules), seen)
                self.assertIn(("OVITO", "rings", threads, molecules), seen)
                self.assertIn(("mdapy", "q6", threads, molecules), seen)
                self.assertIn(("pyscal3", "q6", threads, molecules), seen)
                self.assertIn(("vesin", "neighbours", threads, molecules), seen)
        codes = {row["code"] for row in rows}
        self.assertEqual(codes, {"d-SEAMS", "freud", "OVITO", "mdapy", "pyscal3", "vesin"})

    def test_q_agreement_uses_the_mean_and_ring_definitions_are_explained(self):
        ok, note = lit.agree("q6", "mean=0.6285393611;span=1.000e-16", "mean=0.6285393611;span=2.000e-16")
        self.assertTrue(ok)
        self.assertEqual(note, "")
        ok, note = lit.agree("rings", "six=221184", "six=110592")
        self.assertFalse(ok)
        self.assertIn("Franzblau", note)
        self.assertIn("FindRings", note)

    def test_stopped_process_is_not_runnable(self):
        proc = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            stdout=subprocess.DEVNULL)
        try:
            os.kill(proc.pid, signal.SIGSTOP)
            state = Path(f"/proc/{proc.pid}/status").read_text(encoding="utf-8")
            self.assertIn("State:\tT", state)
            os.kill(proc.pid, signal.SIGCONT)
            state = Path(f"/proc/{proc.pid}/status").read_text(encoding="utf-8")
            self.assertNotIn("State:\tT", state)
        finally:
            proc.kill()
            proc.wait(timeout=5)

    def test_scaling_band_rejects_the_collapsed_median(self):
        self.assertTrue(lit.within_scaling(3.6, 3.452))
        self.assertFalse(lit.within_scaling(70.4, 3.452))
        self.assertFalse(lit.within_scaling(51.9, 2.510))

    def test_fallback_references_are_the_scaling_csv_idle_rows(self):
        self.assertTrue(lit.SCALING_PATH.is_file())
        self.assertEqual(lit.rival_reference(lit.SCALING_PATH), lit.RIVAL_Q6_IDLE_8)
        self.assertEqual(lit.scaling_reference(lit.SCALING_PATH)["q6"], 3.452)
        self.assertEqual(lit.RIVAL_Q6_IDLE_8["mdapy"][0], 14.741)
        self.assertEqual(lit.RIVAL_Q6_IDLE_8["pyscal3"][0], 21.875)

    def test_rival_band_is_the_wider_of_the_fence_and_the_relative_bound(self):
        low, high = lit.rival_band(*lit.RIVAL_Q6_IDLE_8["mdapy"])
        self.assertAlmostEqual(low, 14.741 * 0.88)
        self.assertAlmostEqual(high, 14.741 * 1.12)
        low, high = lit.rival_band(*lit.RIVAL_Q6_IDLE_8["pyscal3"])
        self.assertAlmostEqual(low, 20.542 - 3 * (22.497 - 20.542))
        self.assertAlmostEqual(high, 22.497 + 3 * (22.497 - 20.542))

    def test_gate_passes_the_committed_v93_sitting(self):
        self.assertEqual(lit.check_csv(lit.DEFAULT_CSV), [])

    def test_gate_fails_the_v92_sitting_on_dseams_and_mdapy(self):
        failures = lit.check_csv(V92)
        self.assertTrue(any(line.startswith("d-SEAMS q6 8 threads median 70.433") for line in failures))
        self.assertTrue(any(line.startswith("d-SEAMS q6 8 threads median 75.892") for line in failures))
        self.assertTrue(any(line.startswith("mdapy q6 at 8 threads median 17.046") for line in failures))

    def test_each_gate_fails_alone_on_a_v92_median(self):
        rows = v93_rows()
        for code, v92_median in (("d-SEAMS", "70.433"), ("mdapy", "17.046")):
            changed = [
                dict(row, median_ms=v92_median) if gate_row(row, code) else row
                for row in rows
            ]
            failures = lit.idle_gate(changed)
            self.assertEqual(len(failures), 1, failures)
            self.assertTrue(failures[0].startswith(code), failures)

    def test_gate_requires_the_q6_rows(self):
        rows = [row for row in v93_rows() if not gate_row(row, "d-SEAMS")]
        self.assertIn("d-SEAMS q6 at 8 threads was not measured", lit.idle_gate(rows))

    def test_check_only_times_nothing_and_exits_on_the_gate(self):
        script = ROOT / "scripts" / "lit_one_sitting.py"
        done = subprocess.run(
            [sys.executable, str(script), "--check-only"],
            capture_output=True, text=True, timeout=60)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("gates passed", done.stdout)
        done = subprocess.run(
            [sys.executable, str(script), "--check-only", str(V92)],
            capture_output=True, text=True, timeout=60)
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("GATES FAILED", done.stdout)

    def test_default_table_is_in_the_repository(self):
        self.assertEqual(lit.DEFAULT_CSV, ROOT / "benchmarks" / "literature" / "h2h-lit.csv")


# The v92 sitting left rival pools runnable on the same cores. Rows are
# copied unchanged from its h2h-lit.csv (seams_commit 8d381700).
V92 = ROOT / "tests" / "data" / "lit" / "h2h-lit-v92-q6.csv"


def v93_rows() -> list[dict]:
    with lit.DEFAULT_CSV.open(newline="") as handle:
        return list(csv.DictReader(handle))


def gate_row(row: dict, code: str) -> bool:
    return (row["code"] == code and row["operation"] == "q6"
            and row["threads"] == "8" and row["molecules"] == "110592")


if __name__ == "__main__":
    unittest.main()
