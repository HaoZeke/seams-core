# ClavisCage TRACE against d-SEAMS cages

d-SEAMS 2.14.0. ClavisCage 0.3.18. Cores 8. CPU Intel Xeon Processor. One warm-up and eight alternating calls. Serial uses core 0. Eight threads use cores 0-7.

110592 is not divisible by 46 or by 136. The large frame is sII replicated 9 by 9 by 10, 110160 molecules.

d-SEAMS counts closed cages on the oxygen cutoff graph at 3.5 A, with the edge window off. TRACE counts cages on the hydrogen-bond graph built from the out-degree-2 hydrogens in scripts/h2h_literature.py, because the GenIce dumps contain oxygens only. On sI the seven cages TRACE drops each have 10 to 18 oxygen bonds that cross the periodic boundary, and every one of those bonds is still a hydrogen bond under that placement (0 missing). The 5^12 TRACE does return has no wrapping bond. The 2x2x2 tile, where those cages no longer have to wrap a 12 A cell, matches d-SEAMS cage by cage: 16 of 5^12 and 48 of 5^12 6^2, eight copies of the unit-cell census. SQQ counts the same 2 and 6 oxygen cages on the unit cell. d-SEAMS is right about the cages in the coordinates. TRACE misses a cage when the ring has to wrap the small cell.

## sI

46 molecules. TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A along an out-degree-2 orientation of the four nearest oxygen neighbours

512 both 1 d-SEAMS-only 1 TRACE-only 0
512-0 d-SEAMS only, vertices 0 1 4 6 9 14 16 19 22 23 24 25 28 31 35 36 37 38 42 45, oxygen bonds missing from the hydrogen-bond graph: 0 of 30
51262 both 0 d-SEAMS-only 6 TRACE-only 0
51262-0 d-SEAMS only, vertices 0 1 3 4 8 9 10 13 14 15 19 20 25 27 30 33 34 35 38 39 41 43 44 45, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51262-1 d-SEAMS only, vertices 0 2 3 5 8 9 11 13 18 21 22 23 26 27 28 29 31 32 34 35 37 38 42 44, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51262-2 d-SEAMS only, vertices 0 2 4 5 6 10 11 12 16 18 22 24 25 26 27 32 33 36 37 38 39 40 41 44, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51262-3 d-SEAMS only, vertices 1 2 4 5 6 7 10 11 12 13 15 17 19 20 21 24 28 30 31 32 35 36 42 43, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51262-4 d-SEAMS only, vertices 1 6 7 8 12 14 16 17 19 21 22 23 27 29 30 33 34 37 40 41 42 43 44 45, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51262-5 d-SEAMS only, vertices 2 3 5 7 9 14 15 16 17 18 20 23 24 25 26 28 29 30 31 36 39 40 43 45, oxygen bonds missing from the hydrogen-bond graph: 0 of 38
51264 both 0 d-SEAMS-only 0 TRACE-only 0

512-0 has 18 of its 30 oxygen bonds across the periodic boundary. The six 5^12 6^2 cages have 10, 15, 10, 15, 13 and 13 wrapping bonds. The 5^12 that TRACE returns has none, and it fits inside the 12 A cell. Every wrapping bond is still a hydrogen bond, so the drop is the periodic ring, not a missing edge. d-SEAMS is right.

## sI-222

368 molecules. TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A along an out-degree-2 orientation of the four nearest oxygen neighbours

512 both 16 d-SEAMS-only 0 TRACE-only 0
51262 both 48 d-SEAMS-only 0 TRACE-only 0
51264 both 0 d-SEAMS-only 0 TRACE-only 0

## sII

136 molecules. TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A along an out-degree-2 orientation of the four nearest oxygen neighbours

512 both 16 d-SEAMS-only 0 TRACE-only 0
51262 both 0 d-SEAMS-only 0 TRACE-only 0
51264 both 8 d-SEAMS-only 0 TRACE-only 0

## sII-222

1088 molecules. TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A along an out-degree-2 orientation of the four nearest oxygen neighbours

512 both 128 d-SEAMS-only 0 TRACE-only 0
51262 both 0 d-SEAMS-only 0 TRACE-only 0
51264 both 64 d-SEAMS-only 0 TRACE-only 0

## Where 4.29 ms and 12.7 ms come from

Those two numbers are d-SEAMS, not ClavisCage. In `benchmarks/paper-cells/sitting.csv` as committed by `85cf0a24` (`bench: time the quoted kernels in one sitting`, 12 November 2026), the `source_commit` column is `2914162cb8200dc1a5668a54d570130785ccef0c` (`fix(cages): measure tilt-cell edges on the fractional image`). Eight threads, one warm-up, eight calls, `seams` 2.14.0:

- sI 2x2x2, 368 molecules, median 4.285754 ms, which is 4.29 to three significant figures.
- sII 2x2x2, 1088 molecules, median 12.734558 ms, which is 12.7 to three significant figures.

The ClavisCage rows in that same file are 324.664 ms and 739.347 ms. TRACE 0.3.18 has no thread argument (`analyze_waters` takes the water text, the guest text, the box, the kind, a cutoff and an angle; `Cargo.toml` leaves `rayon` commented out), so an eight-thread pin does not speed it up.

The sentence the paper should carry: On the 2x2x2 cells, one thread, d-SEAMS takes 3.72 ms (sI, 368 molecules) and 10.7 ms (sII, 1088 molecules). ClavisCage TRACE 0.3.18 takes 298 ms and 698 ms on those frames. The figures 4.29 ms and 12.7 ms are the eight-thread d-SEAMS medians from commit 2914162c, recorded in sitting.csv by 85cf0a24, and are not TRACE.

## sII-9x9x10

110160 molecules. TRACE cutoff 0.35 nm and angle 35 degrees. Hydrogens sit 0.96 A along an out-degree-2 orientation of the four nearest oxygen neighbours

512 both 12960 d-SEAMS-only 0 TRACE-only 0
51262 both 0 d-SEAMS-only 0 TRACE-only 0
51264 both 6480 d-SEAMS-only 0 TRACE-only 0

TRACE does not expose a thread count. `analyze_waters` has no thread argument, and ClavisCage 0.3.18 leaves `rayon` commented out in `Cargo.toml`. The eight-thread row pins the process to cores 0-7 and times the same serial call: 61957 ms against 61962 ms at one thread. d-SEAMS on that pin is 303 ms.

