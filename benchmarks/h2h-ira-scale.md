# IRA, pinned, beside d-SEAMS

`ira_mod` 2.2.0 is commit `819586078663488a8544797bae47061f9add2554` of IterativeRotationsAssignments. The driver is `scripts/h2h_ira.py`. `scripts/h2h_ira_scale.py` writes `benchmarks/h2h-ira-scale.csv`. `snakemake h2h_ira_scale` regenerates it. Every row carries the d-SEAMS commit and the IRA commit.

On the 4096-molecule ice Ih cell, one thread, eight interleaved calls after one warm-up, d-SEAMS symmetry is 0.756 ms and IRA is 1203.8 ms. Every atom is hexagonal in both.

`cluster()` scans every molecule for each atom. Sixteen atoms take 7.513 ms at 110592 molecules and 48.203 ms at 884736 molecules, 6.42 times, while the cell is 8 times larger. The full-frame estimates are 831 s and 42647 s. A published 27700 ms at 110592 molecules is a different driver, with the 17-atom shells built once. Scaling that linearly is about 222 s at 884736 molecules. 24992 ms is below both estimates, so that row did not classify every molecule.
