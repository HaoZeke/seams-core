#!/usr/bin/env python3
"""The twelve-neighbour waters are in the shifted frame.

Ice Ih plus a 1.5 A Gaussian draw (PCG64 seed 20261113) puts some
oxygens inside 3.5 A of twelve others. The count is an independent
minimum-image loop, the same count freud returns, and no neighbour id
is repeated. A neighbour-list artefact would fail one of those three.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from h2h_summary_v2 import SHIFTED_SHELL_NOTE

EXPECTED = {
    690: (4.490, 15.864, 50.345),
    1673: (16.865, 18.399, 11.795),
    1675: (15.210, 17.036, 12.464),
    2304: (16.902, 34.179, 1.355),
}


def minimum_image(pos, lengths):
    box = np.asarray(lengths, dtype=np.float64)
    delta = pos[:, None, :] - pos[None, :, :]
    delta -= box * np.round(delta / box)
    dist = np.linalg.norm(delta, axis=2)
    np.fill_diagonal(dist, np.inf)
    return dist


def test_shell():
    import freud
    import h2h_mdapy_pyscal3_vesin as h
    import h2h_summary_extra as extra

    pos, lengths = extra.liquid_positions(h)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    dist = minimum_image(pos, lengths)
    nearest = float(dist.min())
    assert abs(nearest - 0.183) < 0.001
    query = freud.locality.AABBQuery(freud.box.Box(*lengths), pos).query(
        pos, {"r_max": 3.5, "exclude_ii": True}
    ).toNeighborList()
    counts = np.asarray(query.neighbor_counts)
    ii = np.asarray(query.point_indices)
    jj = np.asarray(query.query_point_indices)
    for atom, rounded in EXPECTED.items():
        manual = np.flatnonzero(dist[atom] <= 3.5)
        assert len(manual) == 12
        assert len(set(manual.tolist())) == 12
        assert int(counts[atom]) == 12
        mask = ii == atom
        assert len(set(jj[mask].tolist())) == int(mask.sum()) == 12
        got = tuple(round(float(v), 3) for v in pos[atom])
        assert got == rounded
        assert str(atom) in SHIFTED_SHELL_NOTE
    assert "0.183 A" in SHIFTED_SHELL_NOTE
    assert "repeated periodic image" in SHIFTED_SHELL_NOTE


if __name__ == "__main__":
    test_shell()
    print("ok")
