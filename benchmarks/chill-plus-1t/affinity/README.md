# One-thread OVITO-rule CHILL+ on `taskset -c 0`

The 110592-molecule, one-thread row of the parent sitting, measured again
after the library stopped widening its cpu mask past the launch mask.
`before` is the build of the tree this work started from (`98e5f421`),
`after` is the build that keeps the launch mask, and mdapy 2.0.0
`ChillPlus` runs in the same sitting. Ice Ic, cutoff 3.5 A, 3 warm-up
calls, n=30, every worker on `taskset -c 0` with `OMP_NUM_THREADS=1` and
`RAYON_NUM_THREADS=1`, every other worker stopped during a block.

    python3 scripts/chill_plus_sitting.py \
        --before bbdir-before/src --before-tree 98e5f421800a2c8b0e79d4bf026ef76f2212925f \
        --after bbdir-release/src --sizes 110592 --widths 1 --repeats 30 \
        --out benchmarks/chill-plus-1t/affinity
    python3 scripts/chill_plus_rayon_unset.py \
        --before bbdir-before/src --before-tree 98e5f421800a2c8b0e79d4bf026ef76f2212925f \
        --after bbdir-release/src --repeats 30 --out benchmarks/chill-plus-1t/affinity

| file | contents |
|---|---|
| `chill-sitting.csv`, `raw-chill.txt` | before, after and mdapy, every round with its cpu mask |
| `rayon-unset.txt` | before and after with `RAYON_NUM_THREADS` unset, with the worker's thread count |

## Medians (ms), with the interquartile range

| build | ms | cpus_effective |
|---|---:|---:|
| before | 24.326 (24.072-24.498) | 8 |
| after | 23.783 (23.367-24.146) | 1 |
| mdapy | 62.310 (61.513-63.965) | 1 |

All three give `cubic:110592`. The parent sitting had 23.679 ms against
67.095 ms for mdapy; the d-SEAMS row is unchanged and mdapy is within its
spread across sittings.

Every `before` round carries `escaped=1` and the mask `0-7`: the library
moved the process off cpu 0 although it was started on it. Every `after`
round stays on cpu 0.

## Rayon width unset

With `RAYON_NUM_THREADS` unset, `before` takes 23.485 ms (23.417-23.617)
and `after` 23.516 ms (23.446-23.570). Both workers have one thread after
every call; `before` again ends on `0-7`, `after` on `0`. On this host
the CHILL+ OVITO path does not start a rayon pool at one thread, so the
widened mask does not change its time; the eight-worker case is the ring
search, which `tests/test_serial_affinity` covers.
