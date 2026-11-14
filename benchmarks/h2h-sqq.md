# d-SEAMS and SQQ

SQQ 0.5.8 is the PyPI wheel. The code is github.com/pimooni/sqq, tag v0.5.8, commit 17cf76f9a2fc18848321025674a0f4d4adfa26e8. The paper is J. Pang and Q. Sun, J. Phys. Chem. A 130 (41), 8187-8199 (2026), doi 10.1021/acs.jpca.6c05769. That citation matches `sqq/presentation/citation.py`. All 128 `.py` files in the wheel are byte-identical to the tag (`benchmarks/h2h-sqq-provenance.txt`). Rings and the ice stage exist only in SQQ's Python engine (`py`), and that is the engine used here. d-SEAMS is 2.14.0, built in release mode at -O3 and called through `libseams.so`.

The operations SQQ supports, and that are compared here, are primitive ring statistics and the ice-I label. Cages, half cages, quasi cages and SQQ order parameters are turned off. The O-O graph has a 3.5 A cutoff, and CHILL+ uses the published thresholds. No threshold on either side was changed.

The ice label is compared two ways. The whole-tool comparison is d-SEAMS CHILL+ against SQQ's ice stage, which are different rules. The like-for-like comparison runs SQQ's own rule on d-SEAMS rings (`bench/ring_ice.c`), so both sides do the same work: find six-rings, then apply the two-six-ring rule.

## Frames

| Frame | Molecules | Source |
|---|---:|---|
| Ic-4096, Ic-110592 | 4096, 110592 | `scripts/gen_ideal_ice.py` |
| Ih-4096, Ih-110592 | 4096, 110592 | `scripts/gen_ideal_ice.py` |
| mw-liquid-4096 | 4096 | `input/mw215/in.equilibrate`, 16 cells, sha256 1b18b513... |
| mw-liquid-110592 | 110592 | `input/mw215/in.equilibrate`, 48 cells, sha256 a576324b... |
| tip4p2005-4096 | 4096 | `input/tip4p2005/in.bulk`, sha256 4d019021... |
| tip4p-o-250 | 250 | `input/traj/exampleTraj.lammpstrj`, type 2 (oxygen) |

The 215 K mW liquid frames come from LAMMPS 22 Jul 2025 Update 4 running the deck in the repository. The banner and the full hashes are in `benchmarks/h2h-sqq-mw.txt`. The published Hugging Face copy (`d-SEAMS/mw-liquid-215K`) returned HTTP 401 "Repository Not Found". `input/traj/dump-l.lammpstrj` cannot stand in for it: it is not liquid water, at 0.311 A^-3 with 56 neighbours inside 3.5 A.

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
  - `ice-chordless`: every chordless six-cycle of the `seams_neighbours_build` graph, found by a depth-first search over OpenMP threads. This is the ring set SQQ's ice stage counts, so this is the like-for-like call.
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
| tip4p-o-250 | 509 = 509 | 273 = 273 | 34 = 34 | 167 = 167 | 44/250 |

The pairs give d-SEAMS / SQQ. Rings are compared as vertex sets, ring by ring, not as counts.

Both independent re-derivations reproduce their tool's labels on every molecule of every frame:
- CHILL+ matches `libseams` labels everywhere.
- The SQQ rule on the reference graph reproduces SQQ's labels everywhere except one molecule, explained below. On SQQ's own graph it reproduces all 110592.

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
| tip4p-o-250 | 203 = 203 | 250/250 | 234/250 |

Each label column is the census's per-molecule SQQ label against SQQ's own.

- **Chordless census.** Its rings equal the brute-force enumeration on every frame, and SQQ's ice-stage rings on every frame but one. Its labels equal SQQ's on every molecule but one. Both exceptions are in mw-liquid-110592 and come from the edge described in the next section: SQQ's two extra chordless rings and its label for molecule 21148 hold only in SQQ's graph, which lacks that edge. The census's labels match the SQQ rule on the reference graph on all 110592 molecules.
- **Primitive census.** It differs only where SQQ counts a chordless six-ring that is not primitive. Every mismatch (25, 54, 792 and 16 labels on the four non-crystal frames) has a primitive and a chordless six-ring count that differ on the molecule or on one of its neighbours. On the crystals every chordless six-ring is primitive, so the two censuses agree.

### One edge at the cutoff (mw-liquid-110592)

These are the only ring and edge disagreements in the study (`benchmarks/h2h-sqq-edge.txt`).

The pair 21148-30374 is 3.499998773 A apart in float64, which is 1.2e-6 A inside the cutoff. SQQ gets candidate pairs from MDAnalysis `self_capped_distance` at 0.35 + 1e-7 nm. That routine computes in float32, where the distance rounds to 3.500004292 A, so the pair is never returned. SQQ's float64 re-check therefore never sees it. The reference KD tree and d-SEAMS both keep the edge, so SQQ's graph has exactly this one edge fewer.

Each of the four one-sided primitive rings contains that edge:

- **d-SEAMS-only 5-ring** [414, 11936, 18660, 21148, 30374]. It passes the shortest-path test in the reference graph, so it is a primitive ring and SQQ misses it.
- **SQQ-only 4-ring** [2917, 4879, 21148, 30374] and **6-rings** [414, 2917, 11936, 18660, 21148, 30374] and [414, 4879, 11936, 18660, 21148, 30374]. Each is chordless only in SQQ's graph. With the missing edge restored, 21148-30374 is a chord of each, so none is a ring. The same two 6-rings are SQQ's only extra chordless six-rings.
- **Molecule 21148.** With the edge, it has degree 5, and SQQ's own rule calls it liquid. Without the edge it has degree 4, and SQQ labels it interfacial. The chordless census calls it liquid.

On all four rings and on that molecule, d-SEAMS is right and SQQ's result is a float32 rounding artefact. The closest pair in mw-liquid-4096 is also 1.25e-6 A inside the cutoff, but float32 rounds it inside, and the graphs agree there. In tip4p2005-4096 the closest pair is 6.2e-5 A from the cutoff, and the graphs agree.

### Ice labels

The ideal crystals agree on every molecule. Every site has the bulk ice-I ring signature, and both tools call all of them ice (CHILL+ cubic for Ic and hexagonal for Ih).

On the liquids and in the nanotube, the tools answer different questions. Per-molecule disagreement there is a difference in definitions, not a bug: each tool's labels are reproduced exactly from its own published rule.

Counts on the liquid and nanotube frames:

| Frame | CHILL+ ice | CHILL+ interfacial | SQQ ice | SQQ interfacial | Global Q6 |
|---|---:|---:|---:|---:|---:|
| mw-liquid-4096 | 11 | 70 | 192 | 1659 | 0.0097 |
| tip4p2005-4096 | 2 | 2 | 15 | 1114 | 0.0121 |
| mw-liquid-110592 | 492 | 2017 | 5379 | 44742 | 0.0016 |
| tip4p-o-250 | 0 | 0 | 129 | 77 | 0.2234 |
| ideal Ic / Ih | all | 0 | all | 0 | 0.6285 / 0.5303 |

These frames contain essentially no ice I. Global Steinhardt Q6 is near zero for the bulk liquids, against 0.53-0.63 for the crystals.

Few molecules labelled ice by either tool have the ice-I ring signature:
- mw-liquid-4096, tip4p2005-4096 and tip4p-o-250: none of either tool's ice molecules has it.
- mw-liquid-110592: 3 of 492 CHILL+ ice molecules and 13 of 5379 SQQ ice molecules have it.
- Most molecules labelled ice by either tool lie on a 3-, 4- or 5-ring: 182 of 192 SQQ ice in mw-liquid-4096, 14 of 15 in tip4p2005-4096, 5181 of 5379 in mw-liquid-110592, and all 129 in the nanotube.

So on these frames neither tool's ice set is ice I. SQQ's set is larger because its threshold is low: 2 six-rings, against the 12 of a real site, with no orientation test. The nanotube shows the same effect. The 129 SQQ ice molecules are 4-coordinated with 4-coordinated neighbours on a 5 A wide tube, but their bond correlations are not staggered (e.g. c_ij = 0.85, 0.03, -0.07, 0.65), so CHILL+ calls them water.

Nearly every disagreement falls in one cell: CHILL+ water against SQQ interfacial (degree 4 with at least one six-ring). In tip4p2005-4096 that cell holds 903 of the 1127 disagreements, and CHILL+ interClathrate against SQQ interfacial holds another 167. The clause-by-clause counts are in the agreement log.

SQQ's ice stage counts chordless six-rings, not primitive ones. There are 4592 chordless against 3207 primitive in mw-liquid-4096, 4045 against 2485 in tip4p2005-4096, and 123917 against 85092 in mw-liquid-110592. This is the gap the primitive census measures, and it changes neither conclusion.

## Timing

Each of the 16 frame and thread-count pairs runs two groups of calls in one session, serially:
- the rings group: d-SEAMS rings and SQQ rings;
- the ice group: d-SEAMS CHILL+, the primitive census, the chordless census and SQQ ice.

Each group has one warm-up round and then n = 8 timed rounds. Even rounds call in the order above, and odd rounds in reverse. All three d-SEAMS ice rows are paired with the same SQQ ice samples.

- **CPUs.** The machine has 8 logical CPUs. The 1-thread runs pin CPU 0, and the 8-thread runs pin CPUs 0-7, i.e. the whole machine. The pinned set is recorded in each log header. Nothing else ran during the session.
- **Threads.** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` match the thread count. SQQ's `parallel.workers` and `math_threads` are set to the same count. The census runs on `seams_get_num_threads()` threads.
- **What is timed.** Every call starts from coordinates, so the time includes each tool's own neighbour search.
  - The d-SEAMS ring call also finds 3-rings.
  - The SQQ ice call includes the chordless six-ring search its rule needs; CHILL+ needs no rings.
  - The census calls include their six-ring search and the rule, which is the same work as the SQQ ice call.
- **Ratio.** The ratio is sum(SQQ) / sum(d-SEAMS). The interval is a 95% paired-bootstrap percentile interval with B = 10000 and seed 20261114.
- **Raw data.** The raw samples are `benchmarks/h2h-sqq-time-1.txt` and `-8.txt`; the table is `benchmarks/h2h-sqq.csv`. Its `agreement` column is the lowest agreement among the quantities the row's call produces.

All timing was measured in one Snakemake session at commit 2ef12775 ("feat(bench): time SQQ's ice rule on d-SEAMS rings against SQQ"; the subject identifies it if the hash is re-stamped), and both log headers record that the code tree was clean. Later commits change only benchmark outputs and this file.

Times are milliseconds, mean ± sd. In the ice-primitive and ice-chordless rows, the SQQ column is the SQQ ice call.

| Frame | Operation | Threads | d-SEAMS | SQQ | SQQ / d-SEAMS [95% CI] |
|---|---|---:|---:|---:|---:|
| Ic-4096 | rings | 1 | 4.89 ± 0.27 | 608.3 ± 12.6 | 124 [121, 128] |
| Ic-4096 | rings | 8 | 3.81 ± 1.00 | 624.0 ± 18.9 | 164 [142, 196] |
| Ih-4096 | rings | 1 | 4.76 ± 0.18 | 650.1 ± 9.9 | 137 [134, 140] |
| Ih-4096 | rings | 8 | 3.40 ± 0.64 | 653.5 ± 10.7 | 192 [173, 219] |
| mw-liquid-4096 | rings | 1 | 18.50 ± 0.33 | 859.8 ± 11.3 | 46 [46, 47] |
| mw-liquid-4096 | rings | 8 | 7.42 ± 0.95 | 872.3 ± 25.3 | 117 [108, 128] |
| tip4p2005-4096 | rings | 1 | 22.60 ± 0.23 | 1054.4 ± 15.5 | 47 [46, 47] |
| tip4p2005-4096 | rings | 8 | 8.22 ± 1.11 | 1056.1 ± 23.6 | 129 [118, 141] |
| Ic-110592 | rings | 1 | 124.79 ± 5.87 | 18027.8 ± 156.6 | 144 [140, 150] |
| Ic-110592 | rings | 8 | 60.29 ± 3.51 | 18238.0 ± 265.2 | 302 [291, 316] |
| Ih-110592 | rings | 1 | 132.69 ± 5.14 | 19406.4 ± 323.6 | 146 [143, 150] |
| Ih-110592 | rings | 8 | 66.25 ± 8.82 | 19951.4 ± 434.2 | 301 [278, 328] |
| mw-liquid-110592 | rings | 1 | 556.89 ± 26.46 | 27272.0 ± 633.7 | 49 [47, 50] |
| mw-liquid-110592 | rings | 8 | 143.11 ± 5.46 | 26157.1 ± 753.5 | 183 [180, 185] |
| tip4p-o-250 | rings | 1 | 0.52 ± 0.07 | 24.6 ± 1.3 | 47 [42, 52] |
| tip4p-o-250 | rings | 8 | 0.54 ± 0.11 | 23.9 ± 0.2 | 44 [38, 50] |
| Ic-4096 | ice | 1 | 1.01 ± 0.13 | 555.7 ± 11.6 | 548 [502, 589] |
| Ic-4096 | ice | 8 | 0.79 ± 0.45 | 562.5 ± 18.6 | 708 [501, 1026] |
| Ih-4096 | ice | 1 | 0.89 ± 0.14 | 596.4 ± 9.9 | 670 [605, 741] |
| Ih-4096 | ice | 8 | 0.60 ± 0.22 | 611.6 ± 9.4 | 1023 [804, 1284] |
| mw-liquid-4096 | ice | 1 | 1.69 ± 0.10 | 731.5 ± 10.1 | 432 [413, 450] |
| mw-liquid-4096 | ice | 8 | 0.95 ± 0.31 | 753.1 ± 34.4 | 790 [652, 964] |
| tip4p2005-4096 | ice | 1 | 1.79 ± 0.12 | 914.3 ± 11.2 | 512 [486, 535] |
| tip4p2005-4096 | ice | 8 | 1.58 ± 0.82 | 916.3 ± 13.2 | 580 [425, 822] |
| Ic-110592 | ice | 1 | 27.39 ± 4.78 | 16937.5 ± 316.0 | 618 [557, 687] |
| Ic-110592 | ice | 8 | 11.80 ± 2.28 | 16973.4 ± 331.4 | 1439 [1277, 1632] |
| Ih-110592 | ice | 1 | 23.24 ± 2.84 | 17753.9 ± 122.4 | 764 [701, 821] |
| Ih-110592 | ice | 8 | 11.18 ± 2.08 | 18366.3 ± 341.0 | 1643 [1456, 1851] |
| mw-liquid-110592 | ice | 1 | 68.95 ± 13.70 | 24022.3 ± 1884.0 | 348 [323, 382] |
| mw-liquid-110592 | ice | 8 | 18.76 ± 2.36 | 22649.2 ± 473.2 | 1207 [1130, 1313] |
| tip4p-o-250 | ice | 1 | 0.10 ± 0.02 | 19.3 ± 0.4 | 188 [165, 212] |
| tip4p-o-250 | ice | 8 | 0.18 ± 0.07 | 19.3 ± 0.2 | 109 [86, 149] |
| Ic-4096 | ice-primitive | 1 | 5.35 ± 0.17 | 555.7 ± 11.6 | 104 [102, 105] |
| Ic-4096 | ice-primitive | 8 | 3.91 ± 0.99 | 562.5 ± 18.6 | 144 [124, 169] |
| Ih-4096 | ice-primitive | 1 | 5.04 ± 0.14 | 596.4 ± 9.9 | 118 [115, 121] |
| Ih-4096 | ice-primitive | 8 | 3.42 ± 0.53 | 611.6 ± 9.4 | 179 [160, 195] |
| mw-liquid-4096 | ice-primitive | 1 | 19.18 ± 0.19 | 731.5 ± 10.1 | 38 [38, 39] |
| mw-liquid-4096 | ice-primitive | 8 | 7.70 ± 0.47 | 753.1 ± 34.4 | 98 [93, 103] |
| tip4p2005-4096 | ice-primitive | 1 | 23.30 ± 0.12 | 914.3 ± 11.2 | 39 [39, 40] |
| tip4p2005-4096 | ice-primitive | 8 | 8.41 ± 0.83 | 916.3 ± 13.2 | 109 [102, 116] |
| Ic-110592 | ice-primitive | 1 | 139.22 ± 5.71 | 16937.5 ± 316.0 | 122 [119, 124] |
| Ic-110592 | ice-primitive | 8 | 72.44 ± 6.41 | 16973.4 ± 331.4 | 234 [223, 247] |
| Ih-110592 | ice-primitive | 1 | 133.18 ± 3.74 | 17753.9 ± 122.4 | 133 [131, 136] |
| Ih-110592 | ice-primitive | 8 | 80.87 ± 15.60 | 18366.3 ± 341.0 | 227 [199, 254] |
| mw-liquid-110592 | ice-primitive | 1 | 641.14 ± 80.38 | 24022.3 ± 1884.0 | 37 [36, 39] |
| mw-liquid-110592 | ice-primitive | 8 | 155.64 ± 6.04 | 22649.2 ± 473.2 | 146 [142, 149] |
| tip4p-o-250 | ice-primitive | 1 | 0.52 ± 0.03 | 19.3 ± 0.4 | 37 [36, 39] |
| tip4p-o-250 | ice-primitive | 8 | 1.04 ± 0.34 | 19.3 ± 0.2 | 19 [15, 23] |
| Ic-4096 | ice-chordless | 1 | 4.38 ± 0.17 | 555.7 ± 11.6 | 127 [125, 129] |
| Ic-4096 | ice-chordless | 8 | 1.39 ± 0.57 | 562.5 ± 18.6 | 405 [320, 521] |
| Ih-4096 | ice-chordless | 1 | 4.42 ± 0.15 | 596.4 ± 9.9 | 135 [131, 139] |
| Ih-4096 | ice-chordless | 8 | 1.67 ± 0.90 | 611.6 ± 9.4 | 367 [268, 546] |
| mw-liquid-4096 | ice-chordless | 1 | 20.61 ± 0.13 | 731.5 ± 10.1 | 35 [35, 36] |
| mw-liquid-4096 | ice-chordless | 8 | 6.91 ± 0.55 | 753.1 ± 34.4 | 109 [102, 116] |
| tip4p2005-4096 | ice-chordless | 1 | 28.56 ± 0.13 | 914.3 ± 11.2 | 32 [32, 32] |
| tip4p2005-4096 | ice-chordless | 8 | 7.97 ± 0.44 | 916.3 ± 13.2 | 115 [111, 119] |
| Ic-110592 | ice-chordless | 1 | 98.96 ± 3.62 | 16937.5 ± 316.0 | 171 [168, 175] |
| Ic-110592 | ice-chordless | 8 | 21.01 ± 3.29 | 16973.4 ± 331.4 | 808 [731, 882] |
| Ih-110592 | ice-chordless | 1 | 99.22 ± 3.90 | 17753.9 ± 122.4 | 179 [175, 183] |
| Ih-110592 | ice-chordless | 8 | 21.87 ± 3.10 | 18366.3 ± 341.0 | 840 [777, 923] |
| mw-liquid-110592 | ice-chordless | 1 | 645.01 ± 38.97 | 24022.3 ± 1884.0 | 37 [36, 38] |
| mw-liquid-110592 | ice-chordless | 8 | 86.70 ± 2.86 | 22649.2 ± 473.2 | 261 [256, 267] |
| tip4p-o-250 | ice-chordless | 1 | 0.50 ± 0.04 | 19.3 ± 0.4 | 39 [36, 40] |
| tip4p-o-250 | ice-chordless | 8 | 0.77 ± 0.22 | 19.3 ± 0.2 | 25 [21, 29] |

**Like for like.** The chordless census does SQQ's ice work and reproduces SQQ's labels. It is 32-179 times faster than SQQ on one thread and 25-840 times faster on eight. The whole-tool ratio, CHILL+ against SQQ, is 188-764 on one thread. Most of that ratio comes from CHILL+ needing no ring search, not from faster code. The bulk liquids show it most clearly, at 32-37 times like for like against 348-512 whole-tool. The census ratios sit close to the ring-statistics ratios (46-49 on the bulk liquids), because both are dominated by the six-ring search.

**Threads.** SQQ gains at most 6% from 1 to 8 threads on any frame. d-SEAMS gains up to 3.9 times for rings and 3.7 times for CHILL+ (both on mw-liquid-110592), and the chordless census gains 7.4 times there. So the ratios grow with threads. On the 250-molecule nanotube, eight threads are slower than one for CHILL+ and both censuses: the work per call is under a millisecond, and thread start-up dominates.

## Reproduce

```bash
meson setup bbdir-release --buildtype=release -Doptimization=3 && meson compile -C bbdir-release -j4
pip install sqq==0.5.8 MDAnalysis scipy snakemake asv
git clone --branch v0.5.8 https://github.com/pimooni/sqq sqq-src
export H2H_PYTHON=python3 MW215_DIR=mw215-frames TIP4P2005_DIR=tip4p2005-frames
snakemake -j1 h2h_sqq_mw_liquid h2h_sqq_tip4p2005 h2h_sqq_provenance h2h_sqq_agree h2h_sqq_edge h2h_sqq_time h2h_sqq
asv run --python=same --bench H2HSqq
```

- `h2h_sqq_mw_liquid` needs `lmp` on `PATH`, and only runs when a frame in `MW215_DIR` is missing.
- `h2h_sqq_tip4p2005` needs `lmp` with the OPENMP package on `PATH`, and only runs when `benchmarks/h2h-sqq-tip4p2005.txt` is missing. It took 31 minutes on 8 threads. A regenerated frame need not match the committed hash, because the OpenMP run is not bitwise reproducible across machines. The agreement log records whether it does.
- `meson compile` builds the census as `bbdir-release/bench/libseams_ring_ice.so`; `SEAMS_RING_ICE_LIB` overrides the path.
- Run the rules with `-j1`: the timing rule pins CPU 0 and CPUs 0-7.
- The asv class `H2HSqq` covers the same frames, operations and thread counts with eight repeats. It skips a configuration when `libseams.so`, the census library, sqq or a frame is missing, and runs `time_sqq` only for `rings` and `ice`. The tilted-cell check is `SqqTiltedBox`, so `--bench H2HSqq` selects only this class. `asv.conf.json` benchmarks `HEAD`, whichever branch is checked out.
