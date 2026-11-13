# Ice Ic rings, matched timers

`seams_version_string` is 2.14.0. OVITO is 3.17.2. cycless is 0.7. Eight cores, Intel Xeon Processor. Eight samples after one warm-up. Serial calls are pinned to core 0. Eight-thread calls use cores 0-7. The raw rows are `benchmarks/h2h-ice-ic-rings.csv`. `snakemake ice_ic_rings` regenerates that file.

Two variants, the same work on both sides.

`rings` builds the cutoff graph before the timer. d-SEAMS then calls `ringNetwork`. OVITO calls `FindRingsModifier` on bonds from `CreateBondsModifier`.

`bonds_and_rings` puts that build inside the timer. d-SEAMS rebuilds the neighbour list and searches. OVITO applies `CreateBondsModifier` and `FindRingsModifier`.

The census is 2 hexagons per molecule (2000, 21296, 221184).

## 110592 molecules, milliseconds

Median, then minimum.

| Variant | Threads | d-SEAMS | OVITO |
|---|---:|---|---|
| rings | 1 | 101.912 / 100.673 | 499.459 / 486.623 |
| bonds and rings | 1 | 126.616 / 121.957 | 573.221 / 558.723 |
| rings | 8 | 28.482 / 25.928 | 316.943 / 282.208 |
| bonds and rings | 8 | 42.002 / 37.065 | 297.174 / 286.205 |

On one thread the neighbour list is 24.7 ms of the d-SEAMS call (126.616 minus 101.912) and bond creation is 73.8 ms of the OVITO call (573.221 minus 499.459).

## Where 82.3 vs 474 and 129.3 vs 612.5 come from

`markdown/rings-one-thread.md` prints 82.343 ms against 474.104 ms for 110592 molecules at one thread. The samples stored beside that table, `markdown/rings-one-thread.json`, are 68.773 ms against 415.508 ms. That driver's timing line says the cutoff graph was built before both timers, and the OVITO worker builds bonds once and then times `FindRingsModifier`. The table cell and the stored samples are not the same run. Neither number is a version string. The logs that print d-SEAMS 2.13.0 time `seams_frame_rings`, which builds the cutoff graph inside the call.

The later paper-review row, 129.317 ms against 612.537 ms, is a different operation. `scripts/h2h_paper_review.py` calls `seams_rings_find`, which builds the cutoff graph inside the timer, and it appends `CreateBondsModifier` and `FindRingsModifier` to one pipeline and notifies both on every sample. The refresh of that log in this sitting is 128.638 ms against 628.092 ms. The d-SEAMS side matches the bonds-and-rings row above (126.616 ms). The OVITO side is slower than the fresh-collection bonds-and-rings row (573.221 ms) because that pipeline recomputes both modifiers after the cache is dropped.

So 82.3 against 474 is not the bonds-and-rings comparison, and 129.3 against 612.5 is. The gap is the neighbour list on the d-SEAMS side and `CreateBondsModifier` plus the pipeline restart on the OVITO side. This file is the pair of boundaries on 2.14.0.

## cycless

`cycles_iter` takes a graph and a maximum size. It has no thread argument, so the walk is serial and is pinned to core 0. There is no eight-thread cycless row. On the 1000-molecule cell the serial medians are 3109.4 ms with the graph already built and 3194.1 ms when the neighbour list is inside the timer. Both counts are 2000 hexagons. The same rule is in `scripts/h2h_rings.py`: the eight-thread cycless lines in the older log were this serial walk while the process was allowed cores 0-7.
