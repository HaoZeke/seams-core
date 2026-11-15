# d-SEAMS and SQQ

SQQ 0.5.8 is the PyPI wheel. The code is github.com/pimooni/sqq, tag v0.5.8, commit 17cf76f9a2fc18848321025674a0f4d4adfa26e8. The paper is J. Pang and Q. Sun, J. Phys. Chem. A 130 (41), 8187-8199 (2026), doi 10.1021/acs.jpca.6c05769. That citation matches `sqq/presentation/citation.py`. All 128 `.py` files in the wheel are byte-identical to the tag (`benchmarks/h2h-sqq-provenance.txt`). Rings and the ice stage exist only in SQQ's Python engine (`py`), and that is the engine used here. d-SEAMS is 2.14.0, built in release mode at -O3 and called through `libseams.so`.

The operations SQQ supports, and that are compared here, are primitive ring statistics and the ice-I label. Cages, half cages, quasi cages and SQQ order parameters are turned off. The O-O graph has a 3.5 A cutoff, and CHILL+ uses the published thresholds. No threshold on either side was changed.

The ice label is compared two ways. The whole-tool comparison is d-SEAMS CHILL+ against SQQ's ice stage, which are different rules. The like-for-like comparison runs SQQ's own rule on d-SEAMS rings (`bench/ring_ice.c`), so both sides do the same work: find six-rings, then apply the two-six-ring rule.

## Frames

| Frame | Molecules | Source |
|---|---:|---|
| Ic-4096, Ic-110592, Ic-884736 | 4096, 110592, 884736 | `scripts/gen_ideal_ice.py` |
| Ih-4096, Ih-110592, Ih-884736 | 4096, 110592, 884736 | `scripts/gen_ideal_ice.py` |
| mw-liquid-4096 | 4096 | `input/mw215/in.equilibrate`, 16 cells, sha256 1b18b513... |
| mw-liquid-110592 | 110592 | `input/mw215/in.equilibrate`, 48 cells, sha256 a576324b... |
| mw-liquid-884736 | 884736 | `input/mw215/in.equilibrate`, 96 cells, 8 OpenMP threads, sha256 f7fecdf2... |
| tip4p2005-4096 | 4096 | `input/tip4p2005/in.bulk`, sha256 4d019021... |
| tip4p-o-250 | 250 | `input/traj/exampleTraj.lammpstrj`, type 2 (oxygen) |

The 215 K mW liquid frames come from LAMMPS 22 Jul 2025 Update 4 running the deck in the repository. The banner, the full hashes and each frame's OpenMP thread count are in `benchmarks/h2h-sqq-mw.txt`. The two smaller frames ran on one thread. The 96-cell frame ran on eight with `pair sw/omp` and took 41 minutes; it is fixed by the seed only at that thread count. The published Hugging Face copy (`d-SEAMS/mw-liquid-215K`) returned HTTP 401 "Repository Not Found". `input/traj/dump-l.lammpstrj` cannot stand in for it: it is not liquid water, at 0.311 A^-3 with 56 neighbours inside 3.5 A.

The TIP4P/2005 frame is bulk water at 300 K and 1 atm from the same LAMMPS build.
- **Model.** eps 0.1852 kcal/mol and sigma 3.1589 A on oxygen, qH 0.5564, d_OM 0.1546 A, rigid by Shake.
- **Electrostatics.** `pppm/tip4p` at 1e-4, with a 10 A LJ cutoff and tail correction.
- **Run.** 4096 molecules start on a grid with random orientations. After a minimisation, they melt for 20 ps in NVT and equilibrate for 100 ps in NPT on a 2 fs step.
- **Result.** The last thermo rows give 300-302 K and 0.992-0.999 g/cm^3, against the model's 0.997. The oxygens of the last step are the frame: a cubic box of 49.70 A, global Q6 0.0121.
- **Provenance.** `input/tip4p2005/regenerate.py` writes the frame outside git. The deck, the banner, the thermo rows and the hash are in `benchmarks/h2h-sqq-tip4p2005.txt`, and the hash is also in `input/tip4p2005/tip4p2005-4096.sha256`. The harness records whether the frame it read matches that hash, and the agreement log says "committed".

The 250-molecule TIP4P frame is water confined in a nanotube. It is an extra row, not a bulk liquid. The oxygens' minimum-image xy span is 4.93 by 5.00 A in a box 180 A long.

## Definitions

- **SQQ rings.** `ring.definition: shortest_path` over sizes 4-6. These are chordless cycles that also pass Franzblau's shortest-path test, which is the d-SEAMS primitive-ring predicate. SQQ does not search 3-rings; d-SEAMS does, and its 3-rings match a reference triangle count.
- **SQQ ice** (`ice.method: chill`, `sqq/core/sqq_py/ice.py`).
  - `ice_i`: the molecule has degree 4, lies on at least 2 chordless six-rings, and every neighbour has degree 4.
  - `interfacial`: degree 4 and at least one six-ring, but not `ice_i`.
  - Otherwise `liquid_like`.
  - The rule reads no bond orientation.
- **d-SEAMS CHILL+** (published, `src/bop.cpp`). q3 is computed over the four nearest neighbours inside 3.5 A. A bond is staggered when c_ij <= -0.8 and eclipsed when -0.35 <= c_ij <= 0.25. Cubic means 4 staggered; hexagonal means 3 staggered and 1 eclipsed. The interfacial, clathrate and interClathrate clauses follow the published rule.
- **d-SEAMS ring census** (`bench/ring_ice.c`, built as `libseams_ring_ice.so`). It is SQQ's rule above, on top of the public C API (`seams.h`) and nothing else. It starts from coordinates, like the other calls. It has two ring sets:
  - `ice-primitive`: the primitive six-rings of `seams_rings_find` with a maximum size of 6.
  - `ice-chordless`: every chordless six-cycle of the `seams_neighbours_build` graph, from `seams_rings_chordless` (see the census profile below). This is the ring set SQQ's ice stage counts, so this is the like-for-like call.
- **Merged comparison.** CHILL+ cubic and hexagonal map to ice, CHILL+ interfacial maps to interfacial, and everything else maps to liquid. SQQ `ice_i`, `interfacial` and `liquid_like` map to ice, interfacial and liquid. In the agreement CSV's ice rows, the `d_seams` and `sqq` columns count non-liquid molecules and `both` counts molecules with the same merged label.
- **Arbiter.** A bulk ice-I site lies on exactly 12 primitive six-rings and on no 3-, 4- or 5-ring. Neither rule uses this test, so it can judge both.

Every quantity is also re-derived by code that calls neither tool:
- a SciPy KD-tree graph with a strict `<` cutoff;
- a BFS shortest-path ring test;
- a brute-force enumeration of chordless six-cycles on that graph, on every frame;
- the SQQ rule on that graph;
- CHILL+ with SciPy `sph_harm_y`.

## Agreement

The full table is `benchmarks/h2h-sqq-agree.csv`; the log is `benchmarks/h2h-sqq-agree.txt`. Agreement is shared over union for edges and rings, and the matching fraction for labels. Two empty sets count as full agreement.

| Frame | Edges | 4-rings | 5-rings | 6-rings | Ice label (merged) |
|---|---:|---:|---:|---:|---:|
| Ic-4096 | 8192 = 8192 | 0 = 0 | 0 = 0 | 8192 = 8192 | 4096/4096 |
| Ih-4096 | 8192 = 8192 | 0 = 0 | 0 = 0 | 8192 = 8192 | 4096/4096 |
| mw-liquid-4096 | 9579 = 9579 | 1315 = 1315 | 3190 = 3190 | 3207 = 3207 | 2301/4096 |
| tip4p2005-4096 | 10495 = 10495 | 1758 = 1758 | 2763 = 2763 | 2485 = 2485 | 2969/4096 |
| Ic-110592 | 221184 = 221184 | 0 = 0 | 0 = 0 | 221184 = 221184 | 110592/110592 |
| Ih-110592 | 221184 = 221184 | 0 = 0 | 0 = 0 | 221184 = 221184 | 110592/110592 |
| mw-liquid-110592 | 259187 / 259186 | 34097 / 34098 | 86218 / 86217 | 85092 / 85094 | 62148/110592 |
| Ic-884736 | 1769472 = 1769472 | 0 = 0 | 0 = 0 | 1769472 = 1769472 | 884736/884736 |
| Ih-884736 | 1769472 = 1769472 | 0 = 0 | 0 = 0 | 1769472 = 1769472 | 884736/884736 |
| mw-liquid-884736 | 2075074 / 2075073 | 273457 / 273454 | 692145 / 692147 | 680918 = 680918 | 498791/884736 |
| tip4p-o-250 | 509 = 509 | 273 = 273 | 34 = 34 | 167 = 167 | 44/250 |

The pairs give d-SEAMS / SQQ. Rings are compared as vertex sets, ring by ring, not as counts.

Both independent re-derivations reproduce their tool's labels on every molecule of every frame:
- CHILL+ matches `libseams` labels everywhere.
- The SQQ rule on the reference graph reproduces SQQ's labels everywhere except one molecule of mw-liquid-110592, explained below. On SQQ's own graph it reproduces all of them.

### Like for like

| Frame | Chordless 6-rings | Ice label, chordless census | Ice label, primitive census |
|---|---:|---:|---:|
| Ic-4096 | 8192 = 8192 | 4096/4096 | 4096/4096 |
| Ih-4096 | 8192 = 8192 | 4096/4096 | 4096/4096 |
| mw-liquid-4096 | 4592 = 4592 | 4096/4096 | 4071/4096 |
| tip4p2005-4096 | 4045 = 4045 | 4096/4096 | 4042/4096 |
| Ic-110592 | 221184 = 221184 | 110592/110592 | 110592/110592 |
| Ih-110592 | 221184 = 221184 | 110592/110592 | 110592/110592 |
| mw-liquid-110592 | 123915 / 123917 | 110591/110592 | 109800/110592 |
| Ic-884736 | 1769472 = 1769472 | 884736/884736 | 884736/884736 |
| Ih-884736 | 1769472 = 1769472 | 884736/884736 | 884736/884736 |
| mw-liquid-884736 | 992073 / 992075 | 884736/884736 | 878476/884736 |
| tip4p-o-250 | 203 = 203 | 250/250 | 234/250 |

Each label column is the census's per-molecule SQQ label against SQQ's own.

- **Chordless census.** Its rings equal the brute-force enumeration on every frame, and SQQ's ice-stage rings on every frame except the two mW liquids at 110592 and 884736. Its labels equal SQQ's on every molecule but one. All the exceptions come from one edge per frame that SQQ's graph lacks, described in the next section. In mw-liquid-110592, SQQ's two extra chordless rings and its label for molecule 21148 hold only in SQQ's graph. In mw-liquid-884736 the missing edge changes four chordless rings but no label. The census's labels match the SQQ rule on the reference graph on every molecule of every frame.
- **Bytes.** `seams_rings_chordless` returns canonical rows (lowest vertex first, then its smaller ring neighbour) in sorted order. The agreement log hashes those rows, the census labels and per-molecule six-ring counts, and the same quantities built from the brute-force enumeration and from SQQ's rings. The library and brute-force sha256 digests are equal on all 11 frames, and so are the SQQ digests wherever SQQ's graph has every edge. `tests/test_rings_chordless.c` checks the same thing in C, for every ring size from 3 to 8 and for 1 and 4 threads.
- **Primitive census.** It differs only where SQQ counts a chordless six-ring that is not primitive. Every mismatch (25, 54, 792, 6260 and 16 labels on the five non-crystal frames) has a primitive and a chordless six-ring count that differ on the molecule or on one of its neighbours. On the crystals every chordless six-ring is primitive, so the two censuses agree.

**Chordless or primitive.** The two censuses apply the same two-six-ring rule and differ only in which six-rings count. A chordless six-cycle has no chord, but two of its vertices can still be joined by a path through the network that is shorter than either way around the ring. Franzblau's shortest-path test rejects such a cycle, because it is the sum of two smaller rings that are already in the graph. Counting it credits a molecule with six-rings that exist only as combinations of the 4- and 5-rings around it, which is the liquid topology the rule is meant to reject. The double count is why we hold the primitive set to be the correct one, and why d-SEAMS counts it. The switch moves 25, 54, 792 and 16 labels on mw-liquid-4096, tip4p2005-4096, mw-liquid-110592 and tip4p-o-250, and none on the crystals, where every chordless six-ring is primitive. The chordless rule only inflates SQQ's ice and interfacial sets in the liquids.

### One edge at the cutoff (mw-liquid-110592 and mw-liquid-884736)

These are the only ring and edge disagreements in the study (`benchmarks/h2h-sqq-edge.txt`).

The pair 21148-30374 is 3.499998773 A apart in float64, which is 1.2e-6 A inside the cutoff. SQQ gets candidate pairs from MDAnalysis `self_capped_distance` at 0.35 + 1e-7 nm. That routine computes in float32, where the distance rounds to 3.500004292 A, so the pair is never returned. SQQ's float64 re-check therefore never sees it. The reference KD tree and d-SEAMS both keep the edge, so SQQ's graph has exactly this one edge fewer.

Each of the four one-sided primitive rings contains that edge:

- **d-SEAMS-only 5-ring** [414, 11936, 18660, 21148, 30374]. It passes the shortest-path test in the reference graph, so it is a primitive ring and SQQ misses it.
- **SQQ-only 4-ring** [2917, 4879, 21148, 30374] and **6-rings** [414, 2917, 11936, 18660, 21148, 30374] and [414, 4879, 11936, 18660, 21148, 30374]. Each is chordless only in SQQ's graph. With the missing edge restored, 21148-30374 is a chord of each, so none is a ring. The same two 6-rings are SQQ's only extra chordless six-rings.
- **Molecule 21148.** With the edge, it has degree 5, and SQQ's own rule calls it liquid. Without the edge it has degree 4, and SQQ labels it interfacial. The chordless census calls it liquid.

On all four rings and on that molecule, d-SEAMS is right and SQQ's result is a float32 rounding artefact.

mw-liquid-884736 has one such edge too, 577713-605742 at 3.499991286 A, 8.7e-6 A inside the cutoff. Here the float32 distance, 3.499984264 A, is also inside, yet `self_capped_distance` still does not return the pair, so the loss happens earlier in its float32 search. Every one-sided ring in that frame contains the edge. Three 4-rings and one 5-ring are found only by d-SEAMS and are primitive in the reference graph, and three 5-rings are found only by SQQ. Among chordless six-rings, SQQ has three extra and d-SEAMS one. No label changes. SQQ's rule on SQQ's own graph reproduces all 884736 SQQ labels. The closest pair in mw-liquid-4096 is also 1.25e-6 A inside the cutoff, but float32 rounds it inside, and the graphs agree there. In tip4p2005-4096 the closest pair is 6.2e-5 A from the cutoff, and the graphs agree.

### Ice labels

The ideal crystals agree on every molecule. Every site has the bulk ice-I ring signature, and both tools call all of them ice (CHILL+ cubic for Ic and hexagonal for Ih).

On the liquids and in the nanotube, the tools answer different questions. Per-molecule disagreement there is a difference in definitions, not a bug: each tool's labels are reproduced exactly from its own published rule.

Counts on the liquid and nanotube frames:

| Frame | CHILL+ ice | CHILL+ interfacial | SQQ ice | SQQ interfacial | Global Q6 |
|---|---:|---:|---:|---:|---:|
| mw-liquid-4096 | 11 | 70 | 192 | 1659 | 0.0097 |
| tip4p2005-4096 | 2 | 2 | 15 | 1114 | 0.0121 |
| mw-liquid-110592 | 492 | 2017 | 5379 | 44742 | 0.0016 |
| mw-liquid-884736 | 3950 | 15105 | 40357 | 358126 | 0.0004 |
| tip4p-o-250 | 0 | 0 | 129 | 77 | 0.2234 |
| ideal Ic / Ih | all | 0 | all | 0 | 0.6285 / 0.5303 |

These frames contain essentially no ice I. Global Steinhardt Q6 is near zero for the bulk liquids, against 0.53-0.63 for the crystals.

Few molecules labelled ice by either tool have the ice-I ring signature:
- mw-liquid-4096, tip4p2005-4096 and tip4p-o-250: none of either tool's ice molecules has it.
- mw-liquid-110592: 3 of 492 CHILL+ ice molecules and 13 of 5379 SQQ ice molecules have it.
- mw-liquid-884736: 24 of 3950 CHILL+ ice molecules and 106 of 40357 SQQ ice molecules have it.
- Most molecules labelled ice by either tool lie on a 3-, 4- or 5-ring: 182 of 192 SQQ ice in mw-liquid-4096, 14 of 15 in tip4p2005-4096, 5181 of 5379 in mw-liquid-110592, 38789 of 40357 in mw-liquid-884736, and all 129 in the nanotube.

So on these frames neither tool's ice set is ice I. SQQ's set is larger because its threshold is low: 2 six-rings, against the 12 of a real site, with no orientation test. The nanotube shows the same effect. The 129 SQQ ice molecules are 4-coordinated with 4-coordinated neighbours on a 5 A wide tube, but their bond correlations are not staggered (e.g. c_ij = 0.85, 0.03, -0.07, 0.65), so CHILL+ calls them water.

Nearly every disagreement falls in one cell: CHILL+ water against SQQ interfacial (degree 4 with at least one six-ring). In tip4p2005-4096 that cell holds 903 of the 1127 disagreements, and CHILL+ interClathrate against SQQ interfacial holds another 167. The clause-by-clause counts are in the agreement log.

SQQ's ice stage counts chordless six-rings, not primitive ones. There are 4592 chordless against 3207 primitive in mw-liquid-4096, 4045 against 2485 in tip4p2005-4096, 123917 against 85092 in mw-liquid-110592, and 992075 against 680918 in mw-liquid-884736. This is the gap the primitive census measures, and it changes neither conclusion.

## Timing

Each of the 22 frame and thread-count pairs runs two groups of calls in one session, serially:
- the rings group: d-SEAMS rings and SQQ rings;
- the ice group: d-SEAMS CHILL+, the primitive census, the chordless census and SQQ ice.

Each group has one warm-up round and then n = 8 timed rounds. Even rounds call in the order above, and odd rounds in reverse. All three d-SEAMS ice rows are paired with the same SQQ ice samples.

- **CPUs.** The machine has 8 logical CPUs. The 1-thread runs pin CPU 0, and the 8-thread runs pin CPUs 0-7, i.e. the whole machine. The pinned set and core count are recorded in each log header. Nothing else ran during the session. The machine is a virtual machine that its host paused at times, and the timer does not advance while it is paused. Two Ic-884736 one-thread rows each have one slow sample, 2.6 s for rings and 1.65 s for the chordless census, which widens their sd.
- **Threads.** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` match the thread count. SQQ's `parallel.workers` and `math_threads` are set to the same count. The census runs on `seams_get_num_threads()` threads.
- **What is timed.** Every call starts from coordinates, so the time includes each tool's own neighbour search.
  - The d-SEAMS ring call also finds 3-rings.
  - The SQQ ice call includes the chordless six-ring search its rule needs; CHILL+ needs no rings.
  - The census calls include their six-ring search and the rule, which is the same work as the SQQ ice call.
- **Ratio.** The ratio is sum(SQQ) / sum(d-SEAMS). The interval is a 95% paired-bootstrap percentile interval with B = 10000 and seed 20261114.
- **Raw data.** The raw samples are `benchmarks/h2h-sqq-time-1.txt` and `-8.txt`; the table is `benchmarks/h2h-sqq.csv`. Its `agreement` column is the lowest agreement among the quantities the row's call produces.

All timing was measured in one Snakemake session at commit 46c308ee ("feat(bench): run the H2HSqq asv class as a Snakemake rule"; the subject identifies it if the hash is re-stamped). The session ran agree, time, summarize, asv and edge in that order, and every log header records that the code tree was clean. That commit's source_tree is 721a44619adba20b2e29e44ad3dbb292139b8590. Later commits change only benchmark outputs and this file, with one exception: `benchmarks/h2h-sqq-asv.txt` was rerun on its own at commit c4ac6491 ("fix(bench): resolve relative SQQ frame directories from the checkout"), source_tree fe6d53f3ddff334f40c244d0899a23386c455d0a, with the interpreter in `.venv-sqq`. Between the two commits no file under `src/` that goes into `libseams.so` changed; the library-side differences are the header-only `include/seams.hpp` and test targets in `src/meson.build`.

Times are milliseconds, mean ± sd. In the ice-primitive and ice-chordless rows, the SQQ column is the SQQ ice call.

| Frame | Operation | Threads | d-SEAMS | SQQ | SQQ / d-SEAMS [95% CI] |
|---|---|---:|---:|---:|---:|
| Ic-4096 | rings | 1 | 5.69 ± 0.86 | 667.5 ± 46.2 | 117 [106, 127] |
| Ic-4096 | rings | 8 | 3.28 ± 0.88 | 607.7 ± 10.7 | 185 [159, 224] |
| Ih-4096 | rings | 1 | 5.29 ± 0.39 | 705.8 ± 27.3 | 133 [127, 138] |
| Ih-4096 | rings | 8 | 3.32 ± 0.82 | 654.5 ± 11.3 | 197 [169, 231] |
| mw-liquid-4096 | rings | 1 | 19.97 ± 1.72 | 923.4 ± 46.4 | 46 [44, 48] |
| mw-liquid-4096 | rings | 8 | 7.31 ± 0.73 | 858.3 ± 15.0 | 117 [109, 126] |
| tip4p2005-4096 | rings | 1 | 23.25 ± 0.83 | 1082.7 ± 31.3 | 47 [45, 48] |
| tip4p2005-4096 | rings | 8 | 8.06 ± 0.80 | 1050.9 ± 12.9 | 130 [122, 139] |
| Ic-110592 | rings | 1 | 131.20 ± 6.93 | 19020.2 ± 586.2 | 145 [138, 153] |
| Ic-110592 | rings | 8 | 62.44 ± 7.27 | 18033.4 ± 213.0 | 289 [268, 316] |
| Ih-110592 | rings | 1 | 139.29 ± 17.35 | 20391.9 ± 379.6 | 146 [134, 157] |
| Ih-110592 | rings | 8 | 64.63 ± 6.51 | 19642.9 ± 223.8 | 304 [288, 325] |
| mw-liquid-110592 | rings | 1 | 546.02 ± 17.89 | 26961.8 ± 997.4 | 49 [49, 50] |
| mw-liquid-110592 | rings | 8 | 149.47 ± 10.80 | 26373.5 ± 544.1 | 176 [168, 184] |
| Ic-884736 | rings | 1 | 1365.92 ± 511.96 | 169902.8 ± 5810.9 | 124 [99, 150] |
| Ic-884736 | rings | 8 | 512.20 ± 56.78 | 159923.8 ± 1606.3 | 312 [291, 336] |
| Ih-884736 | rings | 1 | 1053.77 ± 65.53 | 170932.4 ± 3728.5 | 162 [157, 166] |
| Ih-884736 | rings | 8 | 490.86 ± 62.53 | 171530.8 ± 2813.4 | 349 [323, 381] |
| mw-liquid-884736 | rings | 1 | 4778.89 ± 109.11 | 247112.0 ± 5639.7 | 52 [51, 53] |
| mw-liquid-884736 | rings | 8 | 1158.80 ± 107.05 | 239335.4 ± 3231.7 | 207 [195, 218] |
| tip4p-o-250 | rings | 1 | 0.52 ± 0.05 | 24.2 ± 0.2 | 47 [44, 50] |
| tip4p-o-250 | rings | 8 | 0.53 ± 0.11 | 24.1 ± 0.6 | 46 [40, 51] |
| Ic-4096 | ice | 1 | 1.18 ± 0.25 | 622.8 ± 42.3 | 528 [451, 598] |
| Ic-4096 | ice | 8 | 0.71 ± 0.32 | 572.6 ± 14.4 | 804 [605, 1037] |
| Ih-4096 | ice | 1 | 0.95 ± 0.14 | 647.8 ± 29.3 | 682 [612, 755] |
| Ih-4096 | ice | 8 | 0.71 ± 0.42 | 602.6 ± 12.7 | 845 [589, 1166] |
| mw-liquid-4096 | ice | 1 | 1.71 ± 0.18 | 743.2 ± 18.8 | 436 [409, 462] |
| mw-liquid-4096 | ice | 8 | 1.35 ± 0.44 | 739.0 ± 22.3 | 547 [442, 696] |
| tip4p2005-4096 | ice | 1 | 1.82 ± 0.15 | 930.9 ± 14.0 | 512 [484, 545] |
| tip4p2005-4096 | ice | 8 | 1.06 ± 0.47 | 925.4 ± 17.4 | 870 [656, 1156] |
| Ic-110592 | ice | 1 | 25.46 ± 4.50 | 17002.0 ± 318.3 | 668 [589, 737] |
| Ic-110592 | ice | 8 | 10.49 ± 1.93 | 16836.3 ± 250.2 | 1605 [1425, 1826] |
| Ih-110592 | ice | 1 | 26.19 ± 6.49 | 18295.9 ± 394.8 | 699 [592, 813] |
| Ih-110592 | ice | 8 | 11.46 ± 1.94 | 18372.6 ± 561.6 | 1603 [1453, 1795] |
| mw-liquid-110592 | ice | 1 | 56.38 ± 6.22 | 22253.5 ± 276.5 | 395 [368, 425] |
| mw-liquid-110592 | ice | 8 | 17.43 ± 2.26 | 22240.0 ± 348.8 | 1276 [1190, 1387] |
| Ic-884736 | ice | 1 | 282.54 ± 43.76 | 153517.0 ± 5549.3 | 543 [500, 591] |
| Ic-884736 | ice | 8 | 94.39 ± 12.27 | 148056.4 ± 2271.2 | 1569 [1438, 1711] |
| Ih-884736 | ice | 1 | 233.55 ± 11.22 | 157410.6 ± 3532.7 | 674 [650, 697] |
| Ih-884736 | ice | 8 | 66.71 ± 9.90 | 153168.2 ± 1530.0 | 2296 [2074, 2531] |
| mw-liquid-884736 | ice | 1 | 609.22 ± 16.63 | 196608.2 ± 1477.7 | 323 [317, 329] |
| mw-liquid-884736 | ice | 8 | 156.09 ± 10.48 | 195904.7 ± 3123.2 | 1255 [1196, 1312] |
| tip4p-o-250 | ice | 1 | 0.11 ± 0.02 | 19.1 ± 0.1 | 179 [155, 206] |
| tip4p-o-250 | ice | 8 | 0.14 ± 0.04 | 19.6 ± 0.4 | 143 [122, 176] |
| Ic-4096 | ice-primitive | 1 | 6.31 ± 0.73 | 622.8 ± 42.3 | 99 [92, 104] |
| Ic-4096 | ice-primitive | 8 | 3.28 ± 0.34 | 572.6 ± 14.4 | 175 [165, 186] |
| Ih-4096 | ice-primitive | 1 | 5.66 ± 0.48 | 647.8 ± 29.3 | 114 [110, 120] |
| Ih-4096 | ice-primitive | 8 | 3.66 ± 0.56 | 602.6 ± 12.7 | 165 [149, 178] |
| mw-liquid-4096 | ice-primitive | 1 | 19.55 ± 0.64 | 743.2 ± 18.8 | 38 [37, 39] |
| mw-liquid-4096 | ice-primitive | 8 | 7.60 ± 0.83 | 739.0 ± 22.3 | 97 [90, 105] |
| tip4p2005-4096 | ice-primitive | 1 | 23.59 ± 0.25 | 930.9 ± 14.0 | 39 [39, 40] |
| tip4p2005-4096 | ice-primitive | 8 | 8.41 ± 0.82 | 925.4 ± 17.4 | 110 [103, 118] |
| Ic-110592 | ice-primitive | 1 | 137.90 ± 3.72 | 17002.0 ± 318.3 | 123 [121, 126] |
| Ic-110592 | ice-primitive | 8 | 68.00 ± 4.62 | 16836.3 ± 250.2 | 248 [236, 257] |
| Ih-110592 | ice-primitive | 1 | 139.20 ± 7.64 | 18295.9 ± 394.8 | 131 [126, 136] |
| Ih-110592 | ice-primitive | 8 | 74.03 ± 4.87 | 18372.6 ± 561.6 | 248 [239, 259] |
| mw-liquid-110592 | ice-primitive | 1 | 564.01 ± 21.48 | 22253.5 ± 276.5 | 39 [38, 40] |
| mw-liquid-110592 | ice-primitive | 8 | 164.32 ± 20.30 | 22240.0 ± 348.8 | 135 [125, 145] |
| Ic-884736 | ice-primitive | 1 | 1287.13 ± 178.17 | 153517.0 ± 5549.3 | 119 [110, 127] |
| Ic-884736 | ice-primitive | 8 | 619.87 ± 41.96 | 148056.4 ± 2271.2 | 239 [229, 248] |
| Ih-884736 | ice-primitive | 1 | 1202.05 ± 26.02 | 157410.6 ± 3532.7 | 131 [129, 133] |
| Ih-884736 | ice-primitive | 8 | 520.93 ± 48.25 | 153168.2 ± 1530.0 | 294 [277, 311] |
| mw-liquid-884736 | ice-primitive | 1 | 4814.61 ± 20.80 | 196608.2 ± 1477.7 | 41 [41, 41] |
| mw-liquid-884736 | ice-primitive | 8 | 1208.53 ± 36.70 | 195904.7 ± 3123.2 | 162 [159, 166] |
| tip4p-o-250 | ice-primitive | 1 | 0.53 ± 0.04 | 19.1 ± 0.1 | 36 [34, 38] |
| tip4p-o-250 | ice-primitive | 8 | 1.07 ± 0.34 | 19.6 ± 0.4 | 18 [15, 23] |
| Ic-4096 | ice-chordless | 1 | 1.64 ± 0.22 | 622.8 ± 42.3 | 379 [351, 415] |
| Ic-4096 | ice-chordless | 8 | 1.15 ± 0.72 | 572.6 ± 14.4 | 497 [352, 775] |
| Ih-4096 | ice-chordless | 1 | 1.39 ± 0.16 | 647.8 ± 29.3 | 466 [438, 501] |
| Ih-4096 | ice-chordless | 8 | 1.08 ± 0.70 | 602.6 ± 12.7 | 557 [382, 917] |
| mw-liquid-4096 | ice-chordless | 1 | 3.40 ± 0.17 | 743.2 ± 18.8 | 219 [211, 228] |
| mw-liquid-4096 | ice-chordless | 8 | 1.87 ± 0.93 | 739.0 ± 22.3 | 395 [297, 578] |
| tip4p2005-4096 | ice-chordless | 1 | 4.30 ± 0.30 | 930.9 ± 14.0 | 216 [206, 228] |
| tip4p2005-4096 | ice-chordless | 8 | 1.45 ± 0.50 | 925.4 ± 17.4 | 640 [517, 806] |
| Ic-110592 | ice-chordless | 1 | 36.18 ± 5.75 | 17002.0 ± 318.3 | 470 [426, 526] |
| Ic-110592 | ice-chordless | 8 | 10.88 ± 1.69 | 16836.3 ± 250.2 | 1547 [1416, 1711] |
| Ih-110592 | ice-chordless | 1 | 33.33 ± 3.40 | 18295.9 ± 394.8 | 549 [515, 592] |
| Ih-110592 | ice-chordless | 8 | 11.63 ± 2.42 | 18372.6 ± 561.6 | 1580 [1398, 1812] |
| mw-liquid-110592 | ice-chordless | 1 | 96.63 ± 3.28 | 22253.5 ± 276.5 | 230 [224, 236] |
| mw-liquid-110592 | ice-chordless | 8 | 21.64 ± 1.55 | 22240.0 ± 348.8 | 1028 [981, 1082] |
| Ic-884736 | ice-chordless | 1 | 551.26 ± 463.61 | 153517.0 ± 5549.3 | 278 [175, 462] |
| Ic-884736 | ice-chordless | 8 | 102.06 ± 24.11 | 148056.4 ± 2271.2 | 1451 [1249, 1692] |
| Ih-884736 | ice-chordless | 1 | 307.08 ± 19.95 | 157410.6 ± 3532.7 | 513 [491, 536] |
| Ih-884736 | ice-chordless | 8 | 87.73 ± 22.29 | 153168.2 ± 1530.0 | 1746 [1478, 2021] |
| mw-liquid-884736 | ice-chordless | 1 | 843.19 ± 18.94 | 196608.2 ± 1477.7 | 233 [230, 237] |
| mw-liquid-884736 | ice-chordless | 8 | 202.83 ± 15.25 | 195904.7 ± 3123.2 | 966 [913, 1018] |
| tip4p-o-250 | ice-chordless | 1 | 0.20 ± 0.02 | 19.1 ± 0.1 | 97 [91, 105] |
| tip4p-o-250 | ice-chordless | 8 | 0.45 ± 0.20 | 19.6 ± 0.4 | 43 [33, 59] |

**Like for like.** The chordless census does SQQ's ice work and reproduces SQQ's labels. It is 97-549 times faster than SQQ on one thread and 43-1746 times faster on eight. The whole-tool ratio, CHILL+ against SQQ, is 179-699 on one thread. On the bulk liquids the census is 216-233 times faster like for like against 323-512 whole-tool, so the gap that remains is mostly CHILL+ needing no ring search. Both are far above the ring-statistics ratios (46-52 on the bulk liquids), because the census needs only six-rings and only their counts per molecule.

**Threads.** SQQ gains at most 10% from 1 to 8 threads on any frame. d-SEAMS gains up to 4.1 times for rings and 3.9 times for CHILL+ (both on mw-liquid-884736), and the chordless census gains 4.2-4.5 times on the mW liquids at 110592 and 884736. So the ratios grow with threads. On the 250-molecule nanotube, eight threads are slower than one for CHILL+ and both censuses: the work per call is under a millisecond, and thread start-up dominates.

**884736 molecules.** No SQQ call exceeded the 600 s budget. The longest, SQQ rings on mw-liquid-884736, took 247 s on one thread. From 110592 to 884736 molecules, 8 times as many, SQQ's time grew 8.3-9.2 times and d-SEAMS's 5.8-11.1 times (15.2 for the noisy Ic-884736 census row), so both stay close to linear.

### Census profile

Callgrind on one thread over the mw-liquid-4096 census, from coordinates to labels, counted 80.1 million instructions in the v2 census. 60.1 million (75%) were in the census's own depth-first search, and the d-SEAMS neighbour pass was 6.7 million (8%). The search grew one path per start vertex and re-scanned the path for chords at every step. `seams_rings_chordless` (`src/chordless.cpp`) replaces it. It splits each cycle at its lowest vertex into two induced half-paths that meet at one vertex, walks every half-path once per start vertex, buckets the long halves by end vertex, and keeps a pair when the interiors are disjoint and no edge joins them. A per-vertex count of path neighbours replaces the chord scan. At commit 46c308ee the same call counts 32.5 million instructions: 13.1 million in the ring search and the unchanged 6.7 million in the neighbour pass. Timed in the session above, the one-thread census went from 20.61 to 3.40 ms on mw-liquid-4096, from 28.56 to 4.30 ms on tip4p2005-4096 and from 645.01 to 96.63 ms on mw-liquid-110592. Its labels and six-ring counts are byte-identical to the old search's.

## Reproduce

```bash
meson setup bbdir-release --buildtype=release -Doptimization=3 && meson compile -C bbdir-release -j4
python3 -m venv .venv-sqq && .venv-sqq/bin/python -m pip install sqq==0.5.8 MDAnalysis scipy snakemake asv
git clone --branch v0.5.8 https://github.com/pimooni/sqq sqq-src
export H2H_PYTHON=.venv-sqq/bin/python MW215_DIR=mw215-frames TIP4P2005_DIR=tip4p2005-frames
snakemake -j1 h2h_sqq_mw_liquid h2h_sqq_tip4p2005 h2h_sqq_provenance h2h_sqq_agree h2h_sqq_edge h2h_sqq_time h2h_sqq h2h_sqq_asv
```

- `h2h_sqq_mw_liquid` needs `lmp` on `PATH`, with the OPENMP package for the 96-cell frame, and only runs when a frame in `MW215_DIR` is missing.
- `h2h_sqq_tip4p2005` needs `lmp` with the OPENMP package on `PATH`, and only runs when `benchmarks/h2h-sqq-tip4p2005.txt` is missing. It took 31 minutes on 8 threads. A regenerated frame need not match the committed hash, because the OpenMP run is not bitwise reproducible across machines. The agreement log records whether it does.
- `meson compile` builds the census as `bbdir-release/bench/libseams_ring_ice.so`; `SEAMS_RING_ICE_LIB` overrides the path.
- Run the rules with `-j1`: the timing rule pins CPU 0 and CPUs 0-7.
- The session took 17 hours, most of it SQQ calls at 884736 molecules and the brute-force enumeration, which the agreement rule now runs on every frame.
- `h2h_sqq_asv` runs `asv run --python=same --machine anon --bench H2HSqq` under `H2H_PYTHON` and keeps its output in `benchmarks/h2h-sqq-asv.txt`, after the commit, the clean-tree check, the source tree, the CPU set and the interpreter's path in the checkout. asv names the environment after the interpreter's absolute path; the log writes the checkout's part of that name as `<root>`. The asv class covers the same operations and thread counts with eight repeats, on every frame except the three at 884736, which would add about 5 hours of SQQ calls. It skips a configuration when `libseams.so`, the census library, sqq or a frame is missing, and runs `time_sqq` only for `rings` and `ice`. The tilted-cell check is `SqqTiltedBox`, so `--bench H2HSqq` selects only this class. `asv.conf.json` benchmarks `HEAD`, whichever branch is checked out.
