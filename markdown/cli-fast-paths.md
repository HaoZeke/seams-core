# CLI times

Ice Ic, orthorhombic, one oxygen per molecule. The machine reported 8 cores. Each row is one warm-up and eight samples, median and interquartile range, with the process pinned to that many cores and `OMP_NUM_THREADS` set to the same number.

Cited reference points, not remeasured in this run:

- mdapy 2.0.0 cutoff CHILL+, 65.6 ms at 1 thread, on the 110592-molecule cell.
- OVITO 3.17.2 polyhedral template matching, 208.7 ms at 8 threads, same cell.
- rings-statistics-matter 0.2.0, 3.56 ms at 8 threads on 1000 molecules. The previous d-SEAMS figure on that cell was 4.02 ms.

| command | molecules | cores | median (ms) | IQR (ms) |
| --- | ---: | ---: | ---: | ---: |
| `read` | 110592 | 1 | 25.74 | 1.01 |
| `chill-plus` | 110592 | 1 | 63.38 | 1.27 |
| `symmetry` | 110592 | 1 | 62.66 | 0.87 |
| `steinhardt --l 4,6` | 110592 | 1 | 92.13 | 1.67 |
| `cages` | 110592 | 1 | 562.49 | 1.22 |
| `rings cutoff depth 6` | 1000 | 1 | 1.45 | 0.01 |
| `read` | 110592 | 8 | 25.78 | 0.31 |
| `chill-plus` | 110592 | 8 | 50.31 | 1.78 |
| `symmetry` | 110592 | 8 | 39.57 | 1.00 |
| `steinhardt --l 4,6` | 110592 | 8 | 57.18 | 1.46 |
| `cages` | 110592 | 8 | 263.24 | 2.97 |
| `rings cutoff depth 6` | 1000 | 8 | 1.98 | 0.18 |
