# Carried rings, measured crossover

Eight cores. Eight samples. The frame is 4096 molecules of a jittered tetrahedral cell, cutoff 3.5 A, rings of size at most 6. The raw lines are `benchmarks/h2h-carry-crossover.txt`.

The partial search is the rows whose neighbour list changed. On one thread it is faster while that set is at most 0.55 of the vertices (2255 of 4096 won, 2742 of 4096 lost). With a team the share is 0.33 (1180 of 4096 won, 1555 of 4096 lost). A partial search that would not fill the team also loses, so that case calls `ringNetwork` as well. The shares are the constants in `RingUpdater::update`.

Below the share, the partial search wins. A local move of the first 1200 atoms changes 850 rows and 0.141 bonds per molecule, which is the old full-pass rule, and the partial search is 8.283 ms against 18.103 ms at one thread and 3.001 ms against 5.063 ms at eight threads.

Above the share the update calls the whole-frame search. On the four-of-five frame that call was 4.8771 ms against a bare `ringNetwork` of 4.6538 ms at eight threads (`c38ed3c592caaf48493fc78bb66644d21e2b96bb`), because the assignment freed the previous frame's rings inside the timer. The search now writes into the buffers it already holds, and it stops the row count once the moved rows exceed the share. The same frame at `41ccf8a6afb2b0422f03b0e7cbd46b6b6e262338` is 4.3006 ms against 5.5584 ms. `benchmarks/h2h-carry-fraction.csv` keeps both rows, with the commit in the row, and the modulus-20 sweep from fraction 0 through 1.
