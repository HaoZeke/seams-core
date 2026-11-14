# d-SEAMS and SQQ

SQQ 0.5.8 is the PyPI wheel. The code is github.com/pimooni/sqq, tag v0.5.8, commit 17cf76f9a2fc18848321025674a0f4d4adfa26e8. The paper is J. Pang and Q. Sun, J. Phys. Chem. A 130 (41), 8187-8199 (2026), doi 10.1021/acs.jpca.6c05769. That citation matches `sqq/presentation/citation.py`. All 128 `.py` files in the wheel are byte-identical to the tag (`benchmarks/h2h-sqq-provenance.txt`). Rings and the ice stage exist only in SQQ's Python engine (`py`), and that is the engine used here. d-SEAMS is 2.14.0, built in release mode at -O3 and called through `libseams.so`.

The operations SQQ supports, and that are compared here, are primitive ring statistics and the ice-I label. Cages, half cages, quasi cages and SQQ order parameters are turned off. The O-O graph has a 3.5 A cutoff, and CHILL+ uses the published thresholds. No threshold on either side was changed.

## Frames

| Frame | Molecules | Source |
|---|---:|---|
| Ic-4096, Ic-110592 | 4096, 110592 | `scripts/gen_ideal_ice.py` |
| Ih-4096, Ih-110592 | 4096, 110592 | `scripts/gen_ideal_ice.py` |
| mw-liquid-4096 | 4096 | `input/mw215/in.equilibrate`, 16 cells, sha256 1b18b513... |
| mw-liquid-110592 | 110592 | `input/mw215/in.equilibrate`, 48 cells, sha256 a576324b... |
| tip4p-o-250 | 250 | `input/traj/exampleTraj.lammpstrj`, type 2 (oxygen) |

The 215 K mW liquid frames come from LAMMPS 22 Jul 2025 Update 4 running the deck in the repository. The banner and the full hashes are in `benchmarks/h2h-sqq-mw.txt`. The published Hugging Face copy (`d-SEAMS/mw-liquid-215K`) returned HTTP 401 "Repository Not Found". `input/traj/dump-l.lammpstrj` cannot stand in for it: it is not liquid water, at 0.311 A^-3 with 56 neighbours inside 3.5 A. The TIP4P frame is water confined in a nanotube. The oxygens' minimum-image xy span is 4.93 by 5.00 A in a box 180 A long.

## Definitions

- **SQQ rings.** `ring.definition: shortest_path` over sizes 4-6. These are chordless cycles that also pass Franzblau's shortest-path test, which is the d-SEAMS primitive-ring predicate. SQQ does not search 3-rings; d-SEAMS does, and its 3-rings match a reference triangle count.
- **SQQ ice** (`ice.method: chill`, `sqq/core/sqq_py/ice.py`).
  - `ice_i`: the molecule has degree 4, lies on at least 2 chordless six-rings, and every neighbour has degree 4.
  - `interfacial`: degree 4 and at least one six-ring, but not `ice_i`.
  - Otherwise `liquid_like`.
  - The rule reads no bond orientation.
- **d-SEAMS CHILL+** (published, `src/bop.cpp`). q3 is computed over the four nearest neighbours inside 3.5 A. A bond is staggered when c_ij <= -0.8 and eclipsed when -0.35 <= c_ij <= 0.25. Cubic means 4 staggered; hexagonal means 3 staggered and 1 eclipsed. The interfacial, clathrate and interClathrate clauses follow the published rule.
- **Merged comparison.** CHILL+ cubic and hexagonal map to ice, CHILL+ interfacial maps to interfacial, and everything else maps to liquid. SQQ `ice_i`, `interfacial` and `liquid_like` map to ice, interfacial and liquid. In the agreement CSV's ice rows, the `d_seams` and `sqq` columns count non-liquid molecules and `both` counts molecules with the same merged label.
- **Arbiter.** A bulk ice-I site lies on exactly 12 primitive six-rings and on no 3-, 4- or 5-ring. Neither rule uses this test, so it can judge both.

Every quantity is also re-derived by code that calls neither tool:
- a SciPy KD-tree graph with a strict `<` cutoff;
- a BFS shortest-path ring test;
- the SQQ rule on that graph;
- CHILL+ with SciPy `sph_harm_y`.

## Agreement

The full table is `benchmarks/h2h-sqq-agree.csv`; the log is `benchmarks/h2h-sqq-agree.txt`. Agreement is shared over union for edges and rings, and the matching fraction for labels. Two empty sets count as full agreement.

| Frame | Edges | 4-rings | 5-rings | 6-rings | Ice label (merged) |
|---|---:|---:|---:|---:|---:|
| Ic-4096 | 8192 = 8192 | 0 = 0 | 0 = 0 | 8192 = 8192 | 4096/4096 |
| Ih-4096 | 8192 = 8192 | 0 = 0 | 0 = 0 | 8192 = 8192 | 4096/4096 |
| mw-liquid-4096 | 9579 = 9579 | 1315 = 1315 | 3190 = 3190 | 3207 = 3207 | 2301/4096 |
| tip4p-o-250 | 509 = 509 | 273 = 273 | 34 = 34 | 167 = 167 | 44/250 |
| Ic-110592 | 221184 = 221184 | 0 = 0 | 0 = 0 | 221184 = 221184 | 110592/110592 |
| Ih-110592 | 221184 = 221184 | 0 = 0 | 0 = 0 | 221184 = 221184 | 110592/110592 |
| mw-liquid-110592 | 259187 / 259186 | 34097 / 34098 | 86218 / 86217 | 85092 / 85094 | 62148/110592 |

The pairs give d-SEAMS / SQQ. Rings are compared as vertex sets, ring by ring, not as counts.

Both independent re-derivations reproduce their tool's labels on every molecule of every frame:
- CHILL+ matches `libseams` labels everywhere.
- The SQQ rule on the reference graph reproduces SQQ's labels everywhere except one molecule, explained below. On SQQ's own graph it reproduces all 110592.

### One edge at the cutoff (mw-liquid-110592)

These are the only ring and edge disagreements in the study (`benchmarks/h2h-sqq-edge.txt`).

The pair 21148-30374 is 3.499998773 A apart in float64, which is 1.2e-6 A inside the cutoff. SQQ gets candidate pairs from MDAnalysis `self_capped_distance` at 0.35 + 1e-7 nm. That routine computes in float32, where the distance rounds to 3.500004292 A, so the pair is never returned. SQQ's float64 re-check therefore never sees it. The reference KD tree and d-SEAMS both keep the edge, so SQQ's graph has exactly this one edge fewer.

Each of the four one-sided rings contains that edge:

- **d-SEAMS-only 5-ring** [414, 11936, 18660, 21148, 30374]. It passes the shortest-path test in the reference graph, so it is a primitive ring and SQQ misses it.
- **SQQ-only 4-ring** [2917, 4879, 21148, 30374] and **6-rings** [414, 2917, 11936, 18660, 21148, 30374] and [414, 4879, 11936, 18660, 21148, 30374]. Each is chordless only in SQQ's graph. With the missing edge restored, 21148-30374 is a chord of each, so none is a ring.
- **Molecule 21148.** With the edge, it has degree 5, and SQQ's own rule calls it liquid. Without the edge it has degree 4, and SQQ labels it interfacial.

On all four rings and on that molecule, d-SEAMS is right and SQQ's result is a float32 rounding artefact. The closest pair in mw-liquid-4096 is also 1.25e-6 A inside the cutoff, but float32 rounds it inside, and the graphs agree there.

### Ice labels

The ideal crystals agree on every molecule. Every site has the bulk ice-I ring signature, and both tools call all of them ice (CHILL+ cubic for Ic and hexagonal for Ih).

On the liquids and in the nanotube, the tools answer different questions. Per-molecule disagreement there is a difference in definitions, not a bug: each tool's labels are reproduced exactly from its own published rule.

Counts on the liquid and nanotube frames:

| Frame | CHILL+ ice | CHILL+ interfacial | SQQ ice | SQQ interfacial | Global Q6 |
|---|---:|---:|---:|---:|---:|
| mw-liquid-4096 | 11 | 70 | 192 | 1659 | 0.0097 |
| mw-liquid-110592 | 492 | 2017 | 5379 | 44742 | 0.0016 |
| tip4p-o-250 | 0 | 0 | 129 | 77 | 0.2234 |
| ideal Ic / Ih | all | 0 | all | 0 | 0.6285 / 0.5303 |

These frames contain essentially no ice I. Global Steinhardt Q6 is near zero for the mW liquid, against 0.53-0.63 for the crystals.

Few molecules labelled ice by either tool have the ice-I ring signature:
- mw-liquid-4096 and tip4p-o-250: none of either tool's ice molecules has it.
- mw-liquid-110592: 3 of 492 CHILL+ ice molecules and 13 of 5379 SQQ ice molecules have it.
- Most molecules labelled ice by either tool lie on a 3-, 4- or 5-ring: 182 of 192 SQQ ice in mw-liquid-4096, 5181 of 5379 in mw-liquid-110592, and all 129 in the nanotube.

So on these frames neither tool's ice set is ice I. SQQ's set is ten times larger because its threshold is low: 2 six-rings, against the 12 of a real site, with no orientation test. The nanotube shows the same effect. The 129 SQQ ice molecules are 4-coordinated with 4-coordinated neighbours on a 5 A wide tube, but their bond correlations are not staggered (e.g. c_ij = 0.85, 0.03, -0.07, 0.65), so CHILL+ calls them water.

Nearly every disagreement falls in one cell: CHILL+ water with 0-2 staggered bonds, against SQQ interfacial (at least two six-rings and a neighbour that is not 4-coordinated). The clause-by-clause counts are in the agreement log.

SQQ's ice stage counts chordless six-rings, not primitive ones. There are 4592 chordless against 3207 primitive in mw-liquid-4096, and 123917 against 85092 in mw-liquid-110592. Switching SQQ's rule to primitive six-rings changes 25, 792 and 16 labels (mw-liquid-4096, mw-liquid-110592, tip4p-o-250). That changes neither conclusion.

## Timing

Each of the 28 configurations has n = 8 paired calls in one session, serial, after one warm-up call per code. The order alternates between d-SEAMS first and SQQ first.

- **CPUs.** The machine has 8 logical CPUs. The 1-thread runs pin CPU 0, and the 8-thread runs pin CPUs 0-7, i.e. the whole machine. The pinned set is recorded in each log header.
- **Threads.** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` match the thread count. SQQ's `parallel.workers` and `math_threads` are set to the same count.
- **What is timed.** Both calls start from coordinates, so the time includes each tool's own neighbour search.
  - The d-SEAMS ring call also finds 3-rings.
  - The SQQ ice call includes the chordless six-ring search its rule needs; CHILL+ needs no rings.
- **Ratio.** The ratio is sum(SQQ) / sum(d-SEAMS). The interval is a 95% paired-bootstrap percentile interval with B = 10000 and seed 20261114.
- **Raw data.** The raw samples are `benchmarks/h2h-sqq-time-1.txt` and `-8.txt`; the table is `benchmarks/h2h-sqq.csv`.

The timing was measured at commit 66e81515. The only later script change fixed the agreement value for two empty ring sets, and does not touch the timed calls.

Times are milliseconds, mean ± sd.

| Frame | Operation | Threads | d-SEAMS | SQQ | SQQ / d-SEAMS [95% CI] |
|---|---|---:|---:|---:|---:|
| Ic-4096 | rings | 1 | 5.57 ± 0.26 | 678.8 ± 16.8 | 122 [119, 124] |
| Ic-4096 | rings | 8 | 4.14 ± 1.01 | 620.6 ± 20.9 | 150 [131, 180] |
| Ih-4096 | rings | 1 | 4.95 ± 0.24 | 680.9 ± 35.3 | 138 [132, 145] |
| Ih-4096 | rings | 8 | 3.26 ± 0.72 | 655.3 ± 13.8 | 201 [177, 235] |
| mw-liquid-4096 | rings | 1 | 18.90 ± 0.72 | 880.7 ± 20.7 | 47 [45, 48] |
| mw-liquid-4096 | rings | 8 | 7.41 ± 0.79 | 861.7 ± 10.2 | 116 [109, 126] |
| tip4p-o-250 | rings | 1 | 0.54 ± 0.09 | 24.0 ± 0.3 | 44 [40, 49] |
| tip4p-o-250 | rings | 8 | 0.53 ± 0.09 | 23.8 ± 0.4 | 45 [40, 49] |
| Ic-110592 | rings | 1 | 144.03 ± 11.25 | 20618.1 ± 1008.9 | 143 [137, 152] |
| Ic-110592 | rings | 8 | 67.00 ± 8.50 | 18495.8 ± 293.5 | 276 [256, 300] |
| Ih-110592 | rings | 1 | 150.27 ± 9.40 | 22306.0 ± 702.3 | 148 [143, 154] |
| Ih-110592 | rings | 8 | 74.80 ± 31.68 | 20462.0 ± 1219.6 | 274 [217, 338] |
| mw-liquid-110592 | rings | 1 | 587.37 ± 39.32 | 30583.9 ± 2043.7 | 52 [50, 54] |
| mw-liquid-110592 | rings | 8 | 161.73 ± 9.93 | 28292.2 ± 1286.9 | 175 [167, 184] |
| Ic-4096 | ice | 1 | 1.32 ± 0.28 | 615.9 ± 20.2 | 468 [415, 537] |
| Ic-4096 | ice | 8 | 1.27 ± 0.69 | 557.6 ± 12.2 | 441 [321, 644] |
| Ih-4096 | ice | 1 | 1.18 ± 0.26 | 644.0 ± 44.8 | 545 [470, 635] |
| Ih-4096 | ice | 8 | 1.48 ± 0.95 | 638.8 ± 74.7 | 431 [297, 733] |
| mw-liquid-4096 | ice | 1 | 1.91 ± 0.30 | 750.7 ± 16.7 | 393 [352, 443] |
| mw-liquid-4096 | ice | 8 | 1.51 ± 0.72 | 732.4 ± 12.7 | 485 [363, 693] |
| tip4p-o-250 | ice | 1 | 0.12 ± 0.03 | 19.6 ± 0.9 | 157 [137, 183] |
| tip4p-o-250 | ice | 8 | 0.12 ± 0.02 | 19.0 ± 0.1 | 161 [142, 189] |
| Ic-110592 | ice | 1 | 41.46 ± 10.79 | 19569.6 ± 1212.2 | 472 [401, 573] |
| Ic-110592 | ice | 8 | 13.38 ± 2.49 | 17059.8 ± 135.8 | 1275 [1149, 1459] |
| Ih-110592 | ice | 1 | 33.17 ± 6.33 | 20977.1 ± 854.8 | 632 [574, 712] |
| Ih-110592 | ice | 8 | 14.74 ± 4.05 | 18555.0 ± 245.3 | 1259 [1084, 1538] |
| mw-liquid-110592 | ice | 1 | 63.90 ± 11.99 | 22980.3 ± 532.0 | 360 [320, 414] |
| mw-liquid-110592 | ice | 8 | 19.46 ± 3.27 | 23151.6 ± 666.1 | 1190 [1058, 1360] |

SQQ gains at most 15% from 1 to 8 threads, and none on the 110592-molecule liquid for ice. d-SEAMS gains up to 3.6 times (mw-liquid-110592 rings), so the ratios grow with threads.

## Reproduce

```bash
meson setup bbdir-release --buildtype=release -Doptimization=3 && meson compile -C bbdir-release -j4
pip install sqq==0.5.8 MDAnalysis scipy snakemake asv
git clone --branch v0.5.8 https://github.com/pimooni/sqq sqq-src
export H2H_PYTHON=python3 MW215_DIR=mw215-frames
snakemake -j1 h2h_sqq_mw_liquid h2h_sqq_provenance h2h_sqq_agree h2h_sqq_edge h2h_sqq_time h2h_sqq
asv run --python=same --bench 'H2HSqq\.'
```

- `h2h_sqq_mw_liquid` needs `lmp` on `PATH`, and only runs when a frame in `MW215_DIR` is missing.
- Run the rules with `-j1`: the timing rule pins CPU 0 and CPUs 0-7.
- The asv class `H2HSqq` (not the older `H2HSqqBox`) covers the same frames, operations and thread counts with eight repeats. It skips a configuration when `libseams.so`, sqq or a frame is missing.
