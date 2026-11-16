#!/usr/bin/env python3
"""IRA shape matching against d-SEAMS ice labels.

ira_mod 2.2.0 is commit 819586078663488a8544797bae47061f9add2554 of
github.com/mammasmias/IterativeRotationsAssignments.

Each atom is the centre of its 16 nearest neighbours. IRA matches that
cluster to the same cluster taken from ideal ice Ih and from ideal ice
Ic, and the label is the template with the smaller Hausdorff distance
when that distance is below 0.5 A. The cross-template distance on the
two ideal clusters is several angstroms, so 0.5 A sits in that gap.
Pinned ABAB, 20 samples. ira_mod is the Python interface of
github.com/mammasmias/IterativeRotationsAssignments.
"""

from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gen_ice_ic import ice_ic  # noqa: E402
from h2h_freud_ovito import (  # noqa: E402
    SYM,
    core_count,
    ice_positions,
    library_path,
    load_seams,
    mixed_frames,
    require,
)
from h2h_ptm import confusion, dseams_sym, pin, stats  # noqa: E402

N = 20
K = 16
CUT = 0.5
KMAX = 1.5


def cluster(xyz, lengths, index: int) -> np.ndarray:
    delta = xyz - xyz[index]
    delta -= lengths * np.round(delta / lengths)
    distance = np.linalg.norm(delta, axis=1)
    keep = np.argsort(distance)[: K + 1]
    return np.ascontiguousarray(delta[keep])


def hausdorff(ira, template, cloud) -> float:
    typ = np.ones(len(template), dtype=np.int32)
    _rot, _tr, _perm, hd = ira.match(len(template), typ, template, len(cloud), typ, cloud, KMAX)
    return float(hd)


def label_one(ira, ih, ic, xyz, lengths, index: int) -> int:
    cloud = cluster(xyz, lengths, index)
    left = hausdorff(ira, ih, cloud)
    right = hausdorff(ira, ic, cloud)
    if left <= right and left <= CUT:
        return 0
    if right < left and right <= CUT:
        return 1
    return -1


def classify_serial(ira, ih, ic, xyz, lengths) -> np.ndarray:
    n = len(xyz)
    return np.fromiter(
        (label_one(ira, ih, ic, xyz, lengths, i) for i in range(n)),
        dtype=np.int32,
        count=n,
    )


def classify_pool(pool, width, ira, ih, ic, xyz, lengths) -> np.ndarray:
    n = len(xyz)
    chunks = np.array_split(np.arange(n), width)

    def run(indices):
        values = [label_one(ira, ih, ic, xyz, lengths, int(index)) for index in indices]
        return np.asarray(indices, dtype=np.int32), np.asarray(values, dtype=np.int32)

    out = np.empty(n, dtype=np.int32)
    for indices, values in pool.map(run, chunks):
        out[indices] = values
    return out


def report(name, width, n, a, b) -> bool:
    am, ai, amin = stats(a)
    bm, bi, bmin = stats(b)
    win = am < bm and amin < bmin
    print(
        f"robust op {name} shape threads {width} molecules {n} "
        f"dseams_median {am:.4f} dseams_iqr {ai:.4f} dseams_min {amin:.4f} "
        f"ira_median {bm:.4f} ira_iqr {bi:.4f} ira_min {bmin:.4f} "
        f"{'WIN' if win else 'LOSE'}",
        flush=True,
    )
    return win


def main():
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    pin(width)
    import ira_mod

    pin(width)
    ira = ira_mod.IRA()
    version = ira.get_version()
    seams = load_seams(library_path())
    require(seams, seams.seams_set_num_threads(width), "threads")
    pin(width)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print(f"samples {N}", flush=True)
    print(f"ira {version[0]}", flush=True)
    print(f"neighbours {K} hausdorff_cut_A {CUT} kmax_factor {KMAX}", flush=True)

    ih_pos, ih_len = ice_positions(4)
    ih_pos = np.ascontiguousarray(ih_pos, dtype=np.float64)
    ih_box = np.asarray(ih_len, dtype=np.float64)
    ic_pts, ic_side = ice_ic(4)
    ic_pos = np.ascontiguousarray(ic_pts, dtype=np.float64)
    ic_box = np.asarray((ic_side, ic_side, ic_side), dtype=np.float64)
    ih = cluster(ih_pos, ih_box, 0)
    ic = cluster(ic_pos, ic_box, 0)
    cross = hausdorff(ira, ih, ic)
    print(f"cross_template_hausdorff_A {cross:.6f}", flush=True)

    frames = []
    for ncell, label in ((8, "Ih4096"), (16, "Ih32768")):
        pos, lengths = ice_positions(ncell)
        frames.append((label, np.ascontiguousarray(pos, dtype=np.float64), np.asarray(lengths)))
    for rep, label in ((8, "Ic4096"), (16, "Ic32768")):
        pts, side = ice_ic(rep)
        frames.append((label, np.ascontiguousarray(pts, dtype=np.float64), np.asarray((side, side, side))))
    for name, (pos, lengths) in mixed_frames():
        frames.append((name, np.ascontiguousarray(pos, dtype=np.float64), np.asarray(lengths)))

    wins = True
    pool = ThreadPoolExecutor(max_workers=width) if width > 1 else None
    try:
        if pool is not None:
            probe = np.ascontiguousarray(ice_positions(2)[0], dtype=np.float64)
            box = np.asarray(ice_positions(2)[1], dtype=np.float64)
            serial = np.array([label_one(ira, ih, ic, probe, box, i) for i in range(len(probe))])
            parallel = classify_pool(pool, width, ira, ih, ic, probe, box)
            if not np.array_equal(serial, parallel):
                raise SystemExit("threaded IRA labels differ")
            print("thread_check 64-molecule ice matches the serial labels", flush=True)
        for name, pos, lengths in frames:
            n = len(pos)
            d_sym = dseams_sym(seams, pos, tuple(lengths))
            found = (
                classify_serial(ira, ih, ic, pos, lengths)
                if pool is None
                else classify_pool(pool, width, ira, ih, ic, pos, lengths)
            )
            d_name = [SYM.get(int(v), "other") for v in d_sym]
            i_name = [SYM.get(int(v), "other") for v in found]
            same = sum(a == b for a, b in zip(d_name, i_name))
            print(
                f"check {name} molecules {n} agree {same} of {n} " + confusion(d_name, i_name),
                flush=True,
            )

            def ours(pos=pos, lengths=lengths):
                pin(width)
                dseams_sym(seams, pos, tuple(lengths))

            def other(pos=pos, lengths=lengths):
                pin(width)
                if pool is None:
                    classify_serial(ira, ih, ic, pos, lengths)
                else:
                    classify_pool(pool, width, ira, ih, ic, pos, lengths)

            ours()
            other()
            a, b = [], []
            for _ in range(N):
                t0 = time.perf_counter_ns()
                ours()
                a.append((time.perf_counter_ns() - t0) / 1e6)
                t0 = time.perf_counter_ns()
                other()
                b.append((time.perf_counter_ns() - t0) / 1e6)
            if not report(name, width, n, a, b):
                wins = False
    finally:
        if pool is not None:
            pool.shutdown()
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
