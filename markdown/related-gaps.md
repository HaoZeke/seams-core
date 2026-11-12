# Where OVITO is faster, and what to change

These notes are from the v7 and v8 timings and from the OVITO 3.17.2 manual plus the `b0aeaf41` shell classifier. The kernels are not changed on this branch.

The session is ideal ice Ic unless a row names a 200 ps frame. One warm-up, six calls, median and type-7 IQR. One core is `taskset -c 0`. Eight cores are `taskset -c 0-7`.

## q4 and q6

At 110592 molecules, one thread, OVITO `SteinhardtParametersModifier` is 58.710 ms and `seams_frame_steinhardt` is 72.020 ms. At eight threads the figures are 13.113 ms and 29.327 ms. Both return the same per-atom \(q_l\) (maximum absolute difference 1.4e-7 for freud and under 1e-15 for OVITO and pyscal3). The mean values are \(q_4 = 0.509175\) and \(q_6 = 0.628539\).

What each call does:

- d-SEAMS, `seams_frame_steinhardt` in `bindings/cxx/seams_lang.cpp`. The neighbour argument is not the 3.5 Å cutoff list. It calls `nneigh::kNearestNeighbourList` with \(k = 4\) and a candidate window of `cutoff + 1.5` (5.0 Å), then `chill::steinhardtQl`. On this crystal the four nearest neighbours are the cutoff shell, which is why the values match a cutoff Steinhardt. The search is still a wider k-nearest, returned as `vector<vector<int>>`. `flatten` in `src/steinhardt_parallel.cpp` walks those rows, drops duplicate indices with a per-atom `std::vector`, and copies the bonds into a CSR (`offsets`, `cols`, `dr`).
- `steinhardtQl` always fills both `ql` and `qlBar`. Pass 1 accumulates \(q_{lm}\). Pass 2 (`runPass2`) forms \(q_l\) and then the Lechner–Dellago average over the atom and its neighbours. The head-to-head reads `ql` only. The second pass still runs.
- The default sum is `SEAMS_SUM` unset, which is `Reduction::Left` (`src/ozaki.cpp`). For \(l = 6\) that selects `qlmOneAtomL6Batch`, the Cartesian harmonic, eight bonds at a time. \(l = 4\) uses `qlmOneAtomDr`. sphericart is not linked in the release build this watch timed (`with_sphericart` was off), so the Cartesian kernel is the one that ran.
- OVITO 3.17.2, manual page for the modifier added in 3.17.0. `neighbor_mode=Cutoff`, `cutoff=3.5`, `compute_averaged=False`, `max_l` equal to the order being timed. One cutoff list, uniform weights, and \(q_l\) for every \(l = 2, \ldots, l_{\max}\) in that same call. With averaging off there is no second walk over neighbours. A q6 call therefore also fills q2, q3, q4 and q5, and it is still faster than a single-order d-SEAMS call.

The cutoff neighbour list by itself, `seams_frame_neighbours`, is 25.778 ms at one thread and 23.709 ms at eight on this frame. It barely moves when the thread count goes from 1 to 8. q6 at eight threads is 2.5 times the one-thread time (72.0 to 29.3). OVITO's q6 is 4.5 times faster (58.7 to 13.1). The list build is the part that does not scale. The harmonic-plus-average sits on top of a k-nearest that is at least as expensive as that 26 ms cutoff list, because the window is 5.0 Å and the result is an array of vectors.

Fixes, in the order they pay:

1. Skip pass 2 when the caller does not want \(\bar q_l\). The C entry always writes both buffers. A flag on `seams_frame_steinhardt`, defaulting to the current both-values behaviour, lets the benchmark and the CHILL path ask for \(q_l\) alone.
2. For a cutoff Steinhardt, build the CSR from the cutoff bond list that `seams_frame_neighbours` already uses. Do not go through `kNearestNeighbourList` and `vector<vector<int>>` when the crystal's four neighbours are the cutoff shell. Write `offsets`, `cols` and `dr` in that list's own pass.
3. Evaluate q4 and q6 from one neighbour list. The l=6 Cartesian batch already has the bond vectors. q4 can consume the same `dr`. OVITO's one call covers q2 through q6 for less than our q6 alone.
4. Make the cutoff list scale before expecting q6 to scale. The eight-thread list is the section below. Until that list scales, q6 at eight threads cannot land near OVITO's 13 ms.

## Identify diamond on the 200 ps frames

Seed 91001, 13553 atoms, one thread: OVITO identify diamond 18.717 ms, `classifyIce` on `b0aeaf41` 79.158 ms. The other two seeds are the same shape (about 19 ms against about 80 ms). Eight threads: about 6 ms against about 13 ms. The labels are not the same decision. OVITO names 288 cubic and 39 hexagonal diamond on that seed. The kernel names 57 hexagonal, 278 cubic and 1 stacking-disordered.

OVITO's manual for 3.17.2 states the method, which is Appendix A of Maras, Trushin, Stukowski, Ala-Nissila and Jónsson, Computer Physics Communications 205, 13 (2016), doi:10.1016/j.cpc.2016.04.001. The four nearest neighbours are found. Each of those contributes its own neighbours, giving twelve second neighbours. A common-neighbour fingerprint of those twelve plus the central atom decides the label: an FCC shell is cubic diamond, an HCP shell is hexagonal diamond. The signatures are the integer CNA triples, 12 times (4,2,1) for cubic and 6 times (4,2,1) plus 6 times (4,2,2) for hexagonal. Atoms that fail are then marked as a first or second neighbour of a diamond atom. There is no RMSD and no radius gate.

`diamondSecondShell` in `src/shell_classify.cpp` on `b0aeaf41` is the same staggered-versus-eclipsed split, but only after a radius test: the twelve second neighbours must satisfy `rmax <= rmin * 1.0404` (about two percent). A wider shell returns -1 and `shellFit` runs, which is the permutation overlay of the seventeen-atom shell. The q gate at 0.88 does skip the second shell for a bad tetrahedron, so liquid atoms do not enter Kabsch. A 220 K crystal atom often fails the two-percent test and does enter it. That is the 79 ms. The bond rows are still stored for every atom (`kMaxBond` slots) before the gate.

Fixes:

1. On a failed radius test, compute the CNA fingerprint of the twelve second neighbours before `shellFit`. FCC and HCP are the OVITO labels and they do not need the two-percent shell. Keep `shellFit` for stacking-disordered ice and for cages, which that fingerprint does not name.
2. Let the staggered-versus-eclipsed count run without `rmax/rmin < 1.0404`. A bond that is neither eclipsed (`cosine >= 0.95`) nor staggered (`cosine <= 0.65`) already returns -1 from `diamondSecondShell`. The radius gate is what sends a warm crystal into Kabsch.
3. Build the four-neighbour rows once and read second neighbours from those rows. The fast path currently re-walks each neighbour's bond row and rejects images with a 0.25 Å² test inside the atom loop.

## Why the cutoff list does not scale

vesin 0.6.2 at eight threads is 20.808 ms at 110592 molecules and 2.415 ms at 10648. d-SEAMS on the same session is 27.520 ms and 2.813 ms. At one thread the same list is 25.653 ms and 2.441 ms, so eight threads are slower than one. The source compared here is vesin tag v0.6.2, `vesin/src/cpu_cell_list.cpp` and `vesin/src/threadpool.hpp`. The kernels are not changed on this branch.

`seams_frame_neighbours` calls `bond_list` with the cutoff graph, which is `nneigh::neighListO`. That function returns `cellListRowsThreaded` as soon as `omp_get_max_threads() > 1` and the cloud has at least `kThreadedCellRowsMinAtoms` (2048) atoms. The eight-thread timing is that function. The one-thread timing never enters it, because the OpenMP team size is 1, and falls through to `cutoffRows`. With `threads == 1`, `cutoffRows` places pairs straight into the row arrays and sorts each row (`place`, around lines 212–215).

`cellListRowsThreaded` (from line 528) does four things the one-thread placement does not:

- The cell histogram is serial: fractional coordinates, then a prefix sum, then a fill of `cellAtoms`.
- The parallel region is `#pragma omp parallel` over atoms, created and joined on every call. Each atom walks the 27 cells around it, pushes every candidate into a fresh vector, resizes `qs` and `d2`, and calls `mi_dist2_many` for that atom alone. The row is then `std::sort`ed. On this crystal a row has four partners, so the sort and the library call dominate the distance tests.
- After the parallel region returns, `neighListO` copies every row into a new `vector<vector<int>>` and prepends the atom id. That copy is serial (around lines 1174–1178).
- The orthorhombic cell walk `orthoFullList`, which packs coordinates and skips the periodic wrap on an interior cell, sits later in `neighListO` and is never reached once `cellListRowsThreaded` returns true.

vesin splits cells, not atoms, across `ThreadPool::global()`. Workers stay up between calls. Thread 0 is the caller and writes the main pair buffer. Threads 1 through n−1 append to their own `GrowableNeighborList` with no lock on the push. One later pass copies those buffers into the main pair array. `neighbors.sort()` runs only when `options.sorted` is set. The timed call left it unset, and asked only for the pair indices, so there is no sort and no distance square root. The cutoff test is an inline dot product of the already shifted vector.

The other parallel path in `cutoffRows` has the same shape of overhead and is not what the eight-thread row measured. If `cellListRowsThreaded` returned false, `threads > 1` would allocate one `std::pair` per neighbour, fill that array, `parallelStableSort` it by `(row, partner)`, and scatter (from line 217). At 110592 molecules that is 442368 pairs. Sorting them is why a parallel cutoff list can lose to the one-thread placement. Both parallel strategies lose on this frame. The one-thread placement is the faster of the two d-SEAMS paths.

Fixes, in the order they pay:

1. On an orthorhombic box, do not take `cellListRowsThreaded` ahead of the cutoff placement. Walk cells the way vesin does. Give each thread a private pair buffer, concatenate the buffers, and build the row offsets from the per-thread counts. Use the packed orthorhombic difference that `orthoFullList` already uses for an interior cell. Do not call `mi_dist2_many` once per atom.
2. Do not sort a four-integer row inside the atom loop, and do not `parallelStableSort` the full pair array, unless that caller asked for ordered partners. `seams_frame_neighbours` discards the list, so both sorts are overhead on this timer. A CHILL or ring caller can sort the row it reads.
3. Return the CSR `cutoffRows` already holds. The serial `insert` into a second `vector<vector<int>>` copies every partner again after the list exists.
4. A persistent pool, as in `ThreadPool::global()`, removes the fork on every `seams_frame_neighbours` call. That is the smaller term. The per-atom library call, the row sort, and the serial copy are the 27.5 ms.

## What 064bfef7 changed

`perf/neighbours-threads` at `064bfef72eb67c48cf23e77ed328243cdebefcbf` is the orthorhombic CSR this note asked for. `neighListO` calls `orthoFullCsr` (`src/cutoff_csr.cpp`) when the box has three lengths. That function bins atoms with a per-thread histogram and writes a CSR. The per-atom `cellListRowsThreaded` walk runs only if that CSR returns false. The watch branch does not contain the change.

A release build of that commit, timed against vesin 0.6.2 in one later session, is ahead at every width from 1 to 8. At 110592 molecules the medians are 30.820 ms at one thread and 6.986 ms at eight. vesin is 52.161 ms and 19.892 ms. The full table is in `references/related_work.md`.

`idRowsFromCsr` still copies the CSR into `vector<vector<int>>`, and `seams_frame_neighbours` still discards those rows. The copy is inside the 6.986 ms. Skipping it is the remaining step on that branch. The one-thread median, 30.820 ms, sits above the watch-tree one-thread median from the v8 session, 25.653 ms. Those are different libraries and different sessions. The eight-thread median is far below the watch-tree 27.520 ms.
