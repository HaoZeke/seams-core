# Freud cutoff list, thread count set

`scripts/bench_neighbours_threads.py` calls `freud.parallel.set_num_threads` with the same width as the pin. One warm-up, eight calls, even median. Commit `ef38adb99d5f8236b9d6dcae0d701421fcafed99` is on every line. `snakemake h2h_freud_threads` regenerates `benchmarks/h2h-freud-threads.txt`.

Ice Ic, 110592 molecules, 442368 pairs, the same count as d-SEAMS. Freud's median is 196.0 ms at one thread and 195.7 ms at eight. The eight samples at one thread run from 193.6 to 202.2 ms. A one-run gap of 102 against 103 ms sits inside a spread this wide, and raising the thread count does not move the median. matscipy is not installed in this environment, so those lines say `pairs=absent`.
