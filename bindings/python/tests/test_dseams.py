# Copyright (c) d-SEAMS developers
"""The C API, reached through nanobind, on one cubic cell."""

import math
import os
import unittest

import numpy as np

try:
    from ase import Atoms
except ImportError:
    Atoms = None

import dseams
from dseams import CUTOFF, OVITO, PUBLISHED


def diamond(nrep=2):
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    points = [
        (
            (i + shift[0] + site[0]) * a,
            (j + shift[1] + site[1]) * a,
            (k + shift[2] + site[2]) * a,
        )
        for i in range(nrep)
        for j in range(nrep)
        for k in range(nrep)
        for shift in fcc
        for site in basis
    ]
    length = nrep * a
    return np.asarray(points, dtype=np.float64), length


def chordless_reference():
    """Rows and per-atom counts from chordless_ref, or None outside meson."""
    path = os.environ.get("SEAMS_CHORDLESS_REF")
    traj = os.environ.get("SEAMS_CHORDLESS_TRAJ")
    if not path or not traj:
        return None
    with open(path, encoding="ascii") as handle:
        values = np.array(handle.read().split(), dtype=np.int64)
    n, count, size = (int(v) for v in values[:3])
    rows = values[3 : 3 + count * size].reshape(count, size)
    per = values[3 + count * size :]
    assert per.shape == (n,)
    return traj, size, rows, per


class ApiTest(unittest.TestCase):
    def test_both_rules_and_order(self):
        pos, length = diamond()
        dseams.set_num_threads(1)
        self.assertEqual(dseams.num_threads(), 1)
        self.assertTrue(dseams.version())
        frame = dseams.Frame()
        frame.set_frame(pos, (length, length, length))
        self.assertEqual(frame.count(), pos.shape[0])
        self.assertTrue(np.all(frame.chill_plus()[0] == 0))
        fq4, fq6 = frame.q4q6()
        self.assertEqual(fq4.shape, (pos.shape[0],))
        self.assertGreater(float(fq6.mean()), 0.3)
        foff, fcols = frame.neighbours()
        self.assertEqual(foff.shape[0], pos.shape[0] + 1)
        self.assertEqual(fcols.shape[0], int(np.diff(foff).sum()))
        fsizes, fsym, fhist, ffrac = frame.rings()
        self.assertGreater(int(fhist[6]), 0)
        self.assertEqual(fsym.shape, fsizes.shape)
        self.assertGreater(ffrac, 0.5)
        for rule in (PUBLISHED, OVITO, "paper", "ovito"):
            labels, counts = dseams.chill_plus(pos, (length, length, length), rule=rule)
            self.assertEqual(labels.shape, (pos.shape[0],))
            self.assertTrue(np.all(labels == 0))
            self.assertEqual(int(counts[0]), pos.shape[0])
        dseams.set_num_threads(1)
        paper1, _ = dseams.chill_plus(pos, (length, length, length), rule="paper")
        ovito1, _ = dseams.chill_plus(pos, (length, length, length), rule="ovito")
        dseams.set_num_threads(8)
        paper8, _ = dseams.chill_plus(pos, (length, length, length), rule="paper")
        ovito8, _ = dseams.chill_plus(pos, (length, length, length), rule="ovito")
        self.assertTrue(np.array_equal(paper1, paper8))
        self.assertTrue(np.array_equal(ovito1, ovito8))
        q4, q6 = dseams.q4q6(pos, (length, length, length))
        self.assertGreater(float(q6.mean()), 0.3)
        self.assertTrue(np.all(np.isfinite(q4)))
        offsets, cols = dseams.neighbours(pos, (length, length, length), graph=CUTOFF)
        degree = np.diff(offsets)
        self.assertTrue(np.all(degree == 4))
        self.assertEqual(cols.shape[0], int(degree.sum()))
        sizes, sym, hist, fraction = dseams.rings(pos, (length, length, length))
        self.assertGreater(int(hist[6]), 0)
        self.assertEqual(int(np.sum(sizes == 6)), int(hist[6]))
        self.assertGreater(fraction, 0.5)
        self.assertEqual(sym.shape, sizes.shape)

    def test_rings_chordless_diamond(self):
        pos, length = diamond()
        rows, per = dseams.rings_chordless(pos, (length, length, length), size=6)
        # Every diamond atom sits on twelve six-rings and no shorter cycle.
        self.assertEqual(rows.shape, (pos.shape[0] * 12 // 6, 6))
        self.assertTrue(np.all(per == 12))
        self.assertTrue(np.all(rows[:, 0] == rows.min(axis=1)))
        self.assertTrue(np.all(rows[:, 1] < rows[:, -1]))
        flat = [tuple(r) for r in rows.tolist()]
        self.assertEqual(flat, sorted(flat))
        four, per4 = dseams.rings_chordless(pos, (length, length, length), size=4)
        self.assertEqual(four.shape, (0, 4))
        self.assertTrue(np.all(per4 == 0))
        with self.assertRaises(RuntimeError):
            dseams.rings_chordless(pos, (length, length, length), size=13)

    def test_rings_chordless_matches_c(self):
        ref = chordless_reference()
        if ref is None:
            self.skipTest("SEAMS_CHORDLESS_REF is not set")
        traj, size, rows, per = ref
        frame = dseams.read(traj, frame=1, type_filter=2)
        self.assertEqual(frame.count(), per.shape[0])
        dseams.set_num_threads(1)
        got_rows, got_per = frame.rings_chordless(size, 3.5)
        dseams.set_num_threads(4)
        again, _ = frame.rings_chordless(size=size, cutoff=3.5, type=0, graph=CUTOFF)
        dseams.set_num_threads(0)
        self.assertGreater(rows.shape[0], 0)
        self.assertTrue(np.array_equal(got_rows, rows))
        self.assertTrue(np.array_equal(got_per, per))
        self.assertTrue(np.array_equal(again, rows))

    def test_ase_and_mdanalysis(self):
        pos, length = diamond()
        if Atoms is None:
            self.skipTest("ASE is not installed")
        atoms = Atoms(symbols=["O"] * len(pos), positions=pos, cell=[length, length, length], pbc=True)
        apos, box = dseams.from_ase(atoms)
        labels, _ = dseams.chill_plus(apos, box)
        self.assertTrue(np.all(labels == 0))
        try:
            import MDAnalysis as mda
        except ImportError:
            self.skipTest("MDAnalysis is not installed")

        u = mda.Universe.empty(len(pos), trajectory=True)
        u.atoms.positions = pos
        u.dimensions = [length, length, length, 90.0, 90.0, 90.0]
        mpos, mbox = dseams.from_mdanalysis(u)
        labels, _ = dseams.chill_plus(mpos, mbox)
        self.assertTrue(np.all(labels == 0))

    def test_triclinic_box_is_accepted(self):
        pos, length = diamond()
        box = np.array([0, 0, 0, length, length, length, 0.3, 0.0, 0.0])
        labels, _ = dseams.chill_plus(pos, box, rule=PUBLISHED)
        self.assertEqual(labels.shape, (len(pos),))


if __name__ == "__main__":
    unittest.main()
