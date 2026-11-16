#!/usr/bin/env python3
"""The literature sitting keeps one code runnable and does not pin rivals."""

from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
