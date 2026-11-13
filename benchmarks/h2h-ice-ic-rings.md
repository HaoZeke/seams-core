# Ice Ic rings, matched timers

`seams_version_string` is 2.14.0. OVITO is 3.17.2. cycless is 0.7. Eight cores, Intel Xeon Processor. Eight samples after one warm-up. Serial calls are pinned to core 0. Eight-thread calls use cores 0-7. The raw rows are `benchmarks/h2h-ice-ic-rings.csv`. `snakemake ice_ic_rings` regenerates that file.

Two variants, the same work on both sides.

`rings` builds the cutoff graph before the timer. d-SEAMS then calls `ringNetwork`. OVITO calls `FindRingsModifier` on bonds from `CreateBondsModifier`.

`bonds_and_rings` puts that build inside the timer. d-SEAMS rebuilds the neighbour list and searches. OVITO applies `CreateBondsModifier` and `FindRingsModifier`.

The census is 2 hexagons per molecule (2000, 21296, 221184).

## One sitting

Every row below is from one process. Both OVITO clocks use one pipeline built before the timer. Rings-only notifies `FindRingsModifier`. Bonds and rings notifies `CreateBondsModifier` and `FindRingsModifier`. The two variants alternate which one runs first, the first pass is discarded, and eight samples remain. One thread is core 0. Eight threads are cores 0-7, and the OVITO pool reported that count.

Median, then minimum, milliseconds.

| Molecules | Variant | Threads | d-SEAMS | OVITO |
|---:|---|---:|---|---|
| 1000 | rings | 1 | 1.151 / 1.141 | 3.849 / 3.825 |
| 1000 | bonds and rings | 1 | 1.857 / 1.802 | 4.346 / 4.308 |
| 1000 | rings | 8 | 0.729 / 0.494 | 4.096 / 3.969 |
| 1000 | bonds and rings | 8 | 1.192 / 0.923 | 4.598 / 4.489 |
| 10648 | rings | 1 | 9.864 / 9.697 | 37.804 / 37.669 |
| 10648 | bonds and rings | 1 | 13.141 / 12.394 | 41.872 / 41.514 |
| 10648 | rings | 8 | 4.044 / 3.872 | 26.907 / 26.166 |
| 10648 | bonds and rings | 8 | 6.018 / 5.321 | 33.498 / 29.608 |
| 110592 | rings | 1 | 99.170 / 98.658 | 466.652 / 457.794 |
| 110592 | bonds and rings | 1 | 127.116 / 125.899 | 500.230 / 487.436 |
| 110592 | rings | 8 | 36.371 / 32.442 | 287.152 / 252.499 |
| 110592 | bonds and rings | 8 | 50.806 / 48.501 | 277.778 / 256.658 |

At one thread on 110592 molecules the neighbour list is 27.9 ms of the d-SEAMS call and bond creation is 33.6 ms of the OVITO call. Bonds and rings is slower than rings-only at one thread, and at eight threads on 1000 and 10648 molecules.

At eight threads on 110592 molecules the OVITO medians go the other way by 9.4 ms. The rings-only samples run from 252.5 to 311.8 ms and the bonds-and-rings samples run from 256.7 to 295.6 ms, so the ranges overlap and the bond build is not resolved. An earlier eight-thread sitting built a new pipeline inside the rings-only clock and ran that clock to completion before bonds, and reported 316.943 ms against 297.174 ms. That gap was the pipeline construction and the warm core. It is not a property of FindRings.

The eight-thread rings-only ratio on this sitting is 287.152 / 36.371 = 7.90. The bonds-and-rings ratio is 277.778 / 50.806 = 5.47. A seven-times claim is the rings-only ratio. The bonds-and-rings median is not a lower bound on the ring search.

## Where 82.3 vs 474 and 129.3 vs 612.5 come from

`markdown/rings-one-thread.md` prints 82.343 ms against 474.104 ms for 110592 molecules at one thread. The samples stored beside that table, `markdown/rings-one-thread.json`, are 68.773 ms against 415.508 ms. That driver's timing line says the cutoff graph was built before both timers, and the OVITO worker builds bonds once and then times `FindRingsModifier`. The table cell and the stored samples are not the same run. Neither number is a version string. The logs that print d-SEAMS 2.13.0 time `seams_frame_rings`, which builds the cutoff graph inside the call.

The later paper-review row, 129.317 ms against 612.537 ms, is a different operation. `scripts/h2h_paper_review.py` calls `seams_rings_find`, which builds the cutoff graph inside the timer, and it appends `CreateBondsModifier` and `FindRingsModifier` to one pipeline and notifies both on every sample. This file's one-thread bonds-and-rings row is 127.116 ms against 500.230 ms. The d-SEAMS side matches that review. The OVITO side of the review is slower because that pipeline recomputes both modifiers after the cache is dropped.

So 82.3 against 474 is not the bonds-and-rings comparison, and 129.3 against 612.5 is. This file is the pair of boundaries on 2.14.0.

## cycless

`cycles_iter` takes a graph and a maximum size. It has no thread argument, so the walk is serial and is pinned to core 0. There is no eight-thread cycless row. On the 1000-molecule cell in this sitting the serial medians are 3065.5 ms with the graph already built and 3121.0 ms when the neighbour list is inside the timer. Both counts are 2000 hexagons. The same rule is in `scripts/h2h_rings.py`: the eight-thread cycless lines in the older log were this serial walk while the process was allowed cores 0-7.
