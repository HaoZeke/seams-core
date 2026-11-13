# IRA, pinned, beside d-SEAMS

`ira_mod` 2.2.0 is commit `819586078663488a8544797bae47061f9add2554` of IterativeRotationsAssignments. The driver is `scripts/h2h_ira.py`. `scripts/h2h_ira_scale.py` writes `benchmarks/h2h-ira-scale.csv`. `snakemake h2h_ira_scale` regenerates it. Every row carries the d-SEAMS commit and the IRA commit.

On the 4096-molecule ice Ih cell, one thread, eight interleaved calls after one warm-up, d-SEAMS symmetry is 0.756 ms and IRA is 1203.8 ms. Every atom is hexagonal in both.

`cluster()` scans every molecule for each atom. Sixteen atoms take 7.638 ms at 110592 molecules and 47.905 ms at 884736 molecules, 6.27 times, while the cell is 8 times larger. Those rates are a probe, not a full-frame time. The extrapolated full frames are 845 s and 42383 s. The 884736 shape row is DNF: the full cloud was not classified, the budget is 600 s, and there is no measured median. A published 27700 ms at 110592 molecules is a different driver, with the 17-atom shells built once. Scaling that linearly is about 222 s at 884736 molecules. 24992 ms is below both estimates and is not kept.
