#!/usr/bin/env python3
"""Every O-O edge on which SQQ and d-SEAMS differ, and the rings it moves.

For each edge in one graph and not the other the log gives the pair
distance from the float64 coordinates, the distance after a float32
round trip (MDAnalysis distance kernels are float32), and whether
MDAnalysis `self_capped_distance` with SQQ's arguments returns the
pair. Every ring found by one code only is then tested for the edge,
and so is every chordless six-ring and ice label on which SQQ's rule
applied to d-SEAMS rings (bench/ring_ice.c) and SQQ itself differ.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h2h_sqq_versus as hv  # noqa: E402


def pair_distance(pos, lengths, a, b, dtype):
    p = np.asarray(pos, dtype=dtype)
    L = np.asarray(lengths, dtype=dtype)
    d = p[b] - p[a]
    d = d - np.round(d / L) * L
    return float(np.sqrt(np.sum(d * d, dtype=dtype)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lib", default=str(hv.ROOT / "bbdir-release" / "src" / "libseams.so"))
    parser.add_argument("--ring-ice-lib",
                        default=str(hv.ROOT / "bbdir-release" / "bench" / "libseams_ring_ice.so"))
    parser.add_argument("--frames", nargs="+", default=list(hv.FRAMES))
    parser.add_argument("--log", default="benchmarks/h2h-sqq-edge.txt")
    args = parser.parse_args()
    from MDAnalysis.lib.distances import self_capped_distance
    import MDAnalysis

    engine = hv.Seams(args.lib)
    engine.threads(1)
    census = hv.RingIce(args.ring_ice_lib, engine)
    lines = [f"d-SEAMS {engine.version} commit {hv.code_state()}",
             f"SQQ {hv.sqq_version()} commit {hv.SQQ_COMMIT}; MDAnalysis {MDAnalysis.__version__}",
             f"cutoff {hv.CUTOFF_A} A; SQQ asks self_capped_distance for pairs within "
             f"{hv.CUTOFF_A / 10.0} + 1e-7 nm and keeps a pair when its float64 distance <= 0.35 nm"]
    for name in args.frames:
        frame = hv.load_frame(name)
        sys_ = engine.system(frame)
        ref, adj, margin = hv.reference_graph(frame)
        seams_e = engine.edges(sys_)
        built = hv.sqq_frame(frame)
        res = hv.sqq_run(built, hv.sqq_config("rings", 1))
        sqq_e = hv.sqq_edges(res)
        diff = sorted((seams_e ^ sqq_e) | (seams_e ^ ref) | (sqq_e ^ ref))
        lines.append(f"frame {name} edges reference {len(ref)} d-SEAMS {len(seams_e)} "
                     f"SQQ {len(sqq_e)} differing {len(diff)}")
        if not diff:
            engine.lib.seams_system_destroy(sys_)
            continue
        nm = np.asarray(frame["positions"], dtype=float) / 10.0
        box = np.array([*(frame["lengths"] / 10.0), 90.0, 90.0, 90.0])
        found = self_capped_distance(nm, max_cutoff=hv.CUTOFF_A / 10.0 + 1.0e-7, box=box,
                                     return_distances=True)
        cand = {tuple(sorted(map(int, p))): float(d) for p, d in zip(found[0], found[1])}
        for a, b in diff:
            d64 = pair_distance(frame["positions"], frame["lengths"], a, b, np.float64)
            d32 = pair_distance(frame["positions"], frame["lengths"], a, b, np.float32)
            got = cand.get((a, b))
            lines.append(
                f"edge {name} {a}-{b} reference {(a, b) in ref} d-SEAMS {(a, b) in seams_e} "
                f"SQQ {(a, b) in sqq_e}; float64 {d64:.9f} A (cutoff - d = {hv.CUTOFF_A - d64:.3e}); "
                f"float32 {d32:.9f} A; self_capped_distance "
                + ("returns it at {:.9f} nm".format(got) if got is not None else "does not return it"))
        n = frame["molecules"]
        sqq_adj = [set() for _ in range(n)]
        for x, y in sqq_e:
            sqq_adj[x].add(y)
            sqq_adj[y].add(x)
        handle = engine.rings_handle(sys_)
        dsets = engine.ring_sets(handle)
        engine.lib.seams_rings_destroy(handle)
        ssets = hv.sqq_ring_sets(res)
        ice = hv.sqq_run(built, hv.sqq_config("ice", 1))
        labels = hv.sqq_ice_labels(ice, n)
        six = hv.sqq_ring_sets(ice).get(6, set())
        on_ref, deg_ref, _ = hv.sqq_rule(adj, six, n)
        on_sqq, deg_sqq, _ = hv.sqq_rule(sqq_adj, six, n)
        for i in np.flatnonzero(on_ref != labels):
            ends = [f"{a}-{b}" for a, b in diff if i in (a, b) or i in adj[a] | adj[b]]
            lines.append(
                f"ice {name} molecule {i} SQQ {labels[i]}; SQQ rule on the reference graph "
                f"{on_ref[i]} (degree {deg_ref[i]}), on SQQ's graph {on_sqq[i]} "
                f"(degree {deg_sqq[i]}); next to {', '.join(ends) or 'no differing edge'}")
        lines.append(f"ice {name} SQQ rule on SQQ's own graph reproduces "
                     f"{int(np.sum(on_sqq == labels))}/{n} SQQ labels")
        raw = np.zeros(n, dtype=np.int32)
        count = np.zeros(n, dtype=np.int32)
        chordless = census.run(sys_, "ice-chordless", raw, count, rings=True)
        mine = np.array([hv.RING_ICE_LABEL[int(v)] for v in raw], dtype=object)
        for tag, group in (("d-SEAMS only", chordless - six), ("SQQ only", six - chordless)):
            for ring in sorted(group, key=sorted):
                touches = [f"{a}-{b}" for a, b in diff if a in ring and b in ring]
                why = ("chordless in the reference graph" if hv.is_chordless(ring, adj)
                       else "not chordless in the reference graph")
                if not hv.is_chordless(ring, adj) and hv.is_chordless(ring, sqq_adj):
                    why = (f"chordless only in SQQ's graph; in the reference graph "
                           f"{', '.join(touches)} is a chord of it")
                lines.append(f"chordless {name} {tag} {sorted(ring)} contains "
                             f"{', '.join(touches) or 'no differing edge'}; {why}")
        for i in np.flatnonzero(mine != labels):
            ends = [f"{a}-{b}" for a, b in diff if i in (a, b) or i in adj[a] | adj[b]]
            lines.append(
                f"ice {name} molecule {i} SQQ rule on d-SEAMS chordless 6-rings {mine[i]}, SQQ "
                f"{labels[i]}; SQQ rule on the reference graph {on_ref[i]}; next to "
                f"{', '.join(ends) or 'no differing edge'}")
        for size in (4, 5, 6):
            for tag, group in (("d-SEAMS only", dsets.get(size, set()) - ssets.get(size, set())),
                               ("SQQ only", ssets.get(size, set()) - dsets.get(size, set()))):
                for ring in sorted(group, key=sorted):
                    touches = [f"{a}-{b}" for a, b in diff if a in ring and b in ring]
                    order = hv.ring_cycle(ring, adj) if hv.is_chordless(ring, adj) else None
                    witness = None
                    if order is None:
                        why = "not a chordless cycle of the reference graph"
                        if hv.is_chordless(ring, sqq_adj):
                            why = (f"a chordless cycle of SQQ's graph; in the reference graph "
                                   f"{', '.join(touches)} is a chord of it")
                    elif hv.is_shortest_path_ring(ring, adj):
                        why = "shortest-path ring of the reference graph"
                    else:
                        witness = hv.shortcut_witness(ring, adj)
                        why = "chordless but not shortest-path in the reference graph"
                    if witness:
                        why += (f"; reference shortcut {witness[0]}-{witness[1]} graph distance "
                                f"{witness[3]} < ring distance {witness[2]}")
                    lines.append(f"ring {name} size {size} {tag} {sorted(ring)} contains "
                                 f"{', '.join(touches) or 'no differing edge'}; {why}")
        engine.lib.seams_system_destroy(sys_)
    Path(args.log).write_text("\n".join(hv.scrub(x) for x in lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
