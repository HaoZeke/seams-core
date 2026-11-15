"""Rings, cages, signature cages, and symmetry agree on one ideal lattice."""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bindings", "python"))

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

from frames import cages, rings, rings_chordless, signature_cages, symmetry


def diamond(nrep=2):
    oo = 2.75
    a = oo * 4.0 / math.sqrt(3.0)
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    points = []
    for i in range(nrep):
        for j in range(nrep):
            for k in range(nrep):
                for shift in fcc:
                    for site in basis:
                        points.append(
                            (
                                (i + shift[0] + site[0]) * a,
                                (j + shift[1] + site[1]) * a,
                                (k + shift[2] + site[2]) * a,
                            )
                        )
    return np.asarray(points, dtype=np.float64), (nrep * a, nrep * a, nrep * a)


def main():
    pos, box = diamond()
    try:
        hist = rings(pos, box)
    except OSError as exc:
        # The interpreter and libseams_lang were not built for the same
        # runtime (a conda libstdc++ or a different libpython). Skip.
        sys.stderr.write(f"python/library mismatch: {exc}\n")
        raise SystemExit(77) from exc
    if hist.shape != (7,) or int(hist[6]) <= 0:
        raise SystemExit("ideal lattice has no six-rings")
    per, counts = cages(pos, box)
    if per.shape != (pos.shape[0],) or sum(counts) != pos.shape[0]:
        raise SystemExit("cage populations do not cover the lattice")
    n_cages, n_atoms = signature_cages(pos, box, spec="512")
    if n_cages < 0 or n_atoms < 0:
        raise SystemExit("signature census is negative")
    labels = symmetry(pos, box)
    if labels.shape != (pos.shape[0],):
        raise SystemExit("symmetry labels do not match the lattice")
    rows, per = rings_chordless(pos, box, size=6)
    if rows.shape != (pos.shape[0] * 2, 6) or not np.all(per == 12):
        raise SystemExit("diamond chordless six-cycles are not twelve per atom")
    check_chordless_reference()


def check_chordless_reference():
    """The Frame API against chordless_ref on the same frame."""
    ref = os.environ.get("SEAMS_CHORDLESS_REF")
    traj = os.environ.get("SEAMS_CHORDLESS_TRAJ")
    if not ref or not traj:
        return
    with open(ref, encoding="ascii") as handle:
        values = np.array(handle.read().split(), dtype=np.int64)
    n, count, size = (int(v) for v in values[:3])
    want = values[3 : 3 + count * size].reshape(count, size)
    want_per = values[3 + count * size :]
    rows, per = rings_chordless(traj, size=size, cutoff=3.5, type_filter=2)
    if not (np.array_equal(rows, want) and np.array_equal(per, want_per) and per.shape == (n,)):
        raise SystemExit("Frame API chordless cycles differ from the C reference")


if __name__ == "__main__":
    try:
        main()
    except OSError as exc:
        sys.stderr.write(f"python/library mismatch: {exc}\n")
        raise SystemExit(77) from exc
