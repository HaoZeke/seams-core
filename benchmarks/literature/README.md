# Literature sitting

`h2h-lit.csv` and `raw-lit.txt` are written by `scripts/lit_one_sitting.py`
(Snakemake rule `literature_sitting`): ice Ic at 1000, 10648 and 110592
molecules, 1 and 8 threads, n=30, one process per code, every other process
stopped during a sample.

`python3 scripts/lit_one_sitting.py --check-only [CSV]` applies the gates to
an existing table (default: this `h2h-lit.csv`) and times nothing. A full
sitting applies the same gates before it exits.

## Idle gate

A contaminated sitting is one where other codes' OpenMP pools spin on the
cores being timed. In the v92 sitting that moved d-SEAMS q6 at 8 threads,
110592 molecules, to 70.4 ms (C API 75.9 ms) and mdapy q6 to 17.05 ms. The
gate compares the 8-thread, 110592-molecule medians with idle references
measured independently of the literature sitting.

References come from `bench/omp-scaling/scaling.csv`, written by
`scripts/omp_scaling.py` (raw samples in `bench/omp-scaling/raw-scaling.txt`,
source tree `daea7c72`, the same 8-core Xeon, n=10, one process per row, no
other load):

| row                      | median ms | Q1     | Q3     |
|--------------------------|-----------|--------|--------|
| dseams after q6 idle     | 3.452     | 3.338  | 3.549  |
| mdapy 2.0.0 q6 idle      | 14.741    | 14.656 | 14.958 |
| pyscal3 4.1.0 q6 idle    | 21.875    | 20.542 | 22.497 |

Tolerances:

- d-SEAMS, every operation in the "after" idle rows (neighbours, q4, q6,
  rings, CHILL+): the median must be within half to twice the reference.
  d-SEAMS changes between sittings. scaling.csv puts q6 at 2.360 (before),
  3.452 (after) and 3.879 ms (after-fixed), 0.68x to 1.12x of after. The
  isolated sitting committed in f1861267 measured 2.485 ms (0.72x). The
  contaminated v92 value is 20x.
- mdapy and pyscal3 q6 are pinned versions, so only machine noise moves them.
  The band is `[min(Q1 - 3 IQR, 0.88 median), max(Q3 + 3 IQR, 1.12 median)]`
  of the idle row. 3 IQR is Tukey's outer fence on the reference samples. The
  12% is the largest clean-to-clean shift of a pinned rival q6 recorded before
  v93 (mdapy 13.236 ms in f1861267 against 14.741 ms, 10.2%) plus about one
  reference IQR (0.30 ms, 2.1%). This gives mdapy 12.97-16.51 ms and pyscal3
  14.68-28.36 ms.

These constants come only from scaling.csv and the f1861267 sitting. They
were set before v93 was checked against them. The d-SEAMS q6 and both rival
q6 rows must be present.

On the recorded sittings:

| sitting           | d-SEAMS q6 | mdapy q6 | pyscal3 q6 | gate |
|-------------------|-----------:|---------:|-----------:|------|
| f1861267 isolated | 2.485      | 13.236   | 21.725     | pass |
| v92 contaminated  | 70.433     | 17.046   | 23.044     | fail |
| v93 (this table)  | 3.474      | 14.052   | 22.993     | pass |

pyscal3's v92 shift (+5%) is inside its own idle spread, so pyscal3 alone does
not separate the two sittings. The d-SEAMS and mdapy checks do.
`tests/data/lit/h2h-lit-v92-q6.csv` keeps the v92 q6 rows unchanged for
`tests/test_lit_one_sitting.py`.

`--handoff PATH` (or `LIT_HANDOFF`) writes a second copy of the table. The
default is this `h2h-lit.csv`, so a plain run writes one file.
