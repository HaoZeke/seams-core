# One-thread OVITO-rule CHILL+

The timed call is `seams_chill_plus` with the OVITO rule on the cutoff
graph, the d-SEAMS CHILL+ row of the literature table. The base build
(source tree `00129f1e`) is compared with the build that puts that rule
on the in-tree cutoff CSR (source tree `b36fb066`), with mdapy 2.0.0
`ChillPlus` in the same sitting. Ice Ic, cutoff 3.5 A, 3 warm-up calls,
n=30, one worker process per build and width, every other worker stopped
during a block. The machine is an 8-core Intel Xeon cloud VM.

    CHILL_BEFORE=bbdir-before/src CHILL_BEFORE_TREE=00129f1e... \
        snakemake chill_plus_one_thread chill_plus_profile

| file | contents |
|---|---|
| `chill-sitting.csv`, `raw-chill.txt` | base, current and mdapy, 1000 to 884736 molecules, 1 and 8 threads |
| `no-linkcell/` | the same, with both d-SEAMS builds configured `-Dwith_linkcell=disabled` |
| `literature/` | the full `lit_one_sitting.py` sitting on the current build |
| `profile/` | perf self time at one thread on 110592 molecules, by stage |
| `asv/` | `CutoffIce.time_{published,ovito}_chill_plus` for both builds |
| `equivalence/` | per-atom labels of both builds, every frame, 1 and 8 threads |
| `affinity/` | the 110592-molecule one-thread row on `taskset -c 0`, with each round's cpu mask |

Each file names the tree that measured it (`source_tree`) and, where
two builds appear, the tree each library was built from (`build_tree`).

## Medians (ms), with the interquartile range

| molecules | threads | base | current | mdapy |
|---:|---:|---:|---:|---:|
| 1000 | 1 | 0.474 (0.460-0.519) | 0.223 (0.216-0.231) | 0.893 (0.889-0.903) |
| 1000 | 8 | 0.357 (0.348-0.365) | 0.223 (0.221-0.230) | 0.264 (0.258-0.269) |
| 10648 | 1 | 4.653 (4.603-4.714) | 2.334 (2.325-2.381) | 6.912 (6.852-6.946) |
| 10648 | 8 | 2.485 (2.413-2.507) | 0.587 (0.559-0.651) | 1.585 (1.531-1.639) |
| 110592 | 1 | 51.984 (51.287-53.163) | 23.679 (23.581-24.074) | 67.095 (63.555-69.080) |
| 110592 | 8 | 9.549 (8.517-9.945) | 4.495 (3.988-4.986) | 12.625 (12.394-13.325) |
| 884736 | 1 | 505.668 (493.252-515.664) | 297.759 (293.571-303.480) | 756.426 (736.516-768.908) |
| 884736 | 8 | 94.895 (90.293-100.774) | 42.082 (40.905-44.218) | 128.396 (121.548-133.181) |

Every row of every build gives the same classes as the base
(`cubic:<molecules>`), and mdapy agrees.

Without linkcell (`no-linkcell/`), the base takes 83.302 ms at one
thread on 110592 molecules and the current build 24.824 ms; mdapy takes
66.076 ms. Every base sample there has `ntask=1`, as every d-SEAMS sample
of the v93 literature sitting does: that sitting's 98.9 ms is a build
without linkcell, where the one-thread rows come from vesin.

## Where the time went

`profile/` at one thread on 110592 molecules, scaled to the medians above:

| stage | base | current |
|---|---:|---:|
| neighbour build | 22.1 | 9.4 |
| shell sort, deltas, c_ij and classification | 16.2 | 6.6 |
| q3 harmonics | 5.7 | 5.3 |
| allocation and free | 3.8 | 0.3 |
| API copy and labels | 1.5 | 1.1 |
| other | 2.6 | 0.9 |

The published rule already read the in-tree cutoff CSR. The OVITO rule
built id-keyed `std::vector` rows instead: linkcell's pair list (vesin
without linkcell) converted to id rows, then per row a heap sort of the
shell, an `idIndexMap` lookup per bond, a second shell count per atom and
a release of every row. The current build reads the same CSR as the
published rule, sorts each shell once into a stack buffer, sums q3 into
flat arrays and classifies in a second pass. The q3 work is unchanged.

mdapy at one thread spends about 46 ms building a neighbour list that
stores each pair's distance, vector and image, and about 17 ms in its
native CHILL+ pass over that list.

## Literature sitting

`literature/` is the full sitting, n=30, on the current build. Its
CHILL+ rows match the table above (110592 molecules: 23.680 ms at one
thread against mdapy's 62.776 ms; 3.793 ms against 11.647 ms at eight).
The gate fails on a rival row: mdapy q6 at eight threads takes 12.031 ms,
under the 12.972-16.510 ms band built from `bench/omp-scaling/scaling.csv`.
The base sitting on this machine gives 12.647 ms for the same row and
fails the same gate, so the band reflects the machine that wrote
`scaling.csv`. The committed `benchmarks/literature/` table is unchanged.

## asv

One thread on core 0, through `seams_frame_chill_plus_rule`: the OVITO
rule on 110592 molecules takes 54.5 ms on the base and 30.9 ms now; the
published rule takes 27.3 ms and 26.3 ms.

## Labels

`equivalence/summary.txt`: 3932 cases over 485 frames, the C API and the
language front, both rules, the cutoff and mutual graphs, at one and
eight threads, give the same label array and counts on both builds.
`equivalence/counts.csv` lists the class counts of each case.
