# OpenMM

`DseamsReporter` reads positions from an OpenMM context and hands them to
`InSitu`. The engine calls `libseams_lang` for CHILL+ and for primitive
rings. It does not write a trajectory of its own. The driver in
`scripts/openmm_tip4p_ice.py` does, so the same coordinates can be read
back from disk.

The force field is TIP4P/Ice: oxygen Lennard-Jones σ = 0.31668 nm and
ε = 0.882164484 kJ/mol, hydrogen charge 0.5897 e, M-site charge
−1.1794 e, and an M site 0.01577 nm from the oxygen along the H–O–H
bisector. The box in that driver is 4096 molecules. The platform is
CPU. The timestep is 2 fs. The third argument of the driver is the log.

CHILL+ and the primitive rings share one cutoff neighbour list at
3.5 Å, type 1, on the oxygen sites only. Rings run through size 6.
For a ring of at most 7 members the separation around the ring is at
most 3, so a strict shortcut is a walk of one edge or of two. That
question is answered from the sorted adjacency rows, and the per-vertex
balls are not built. The OpenMP team inside the classifier is the same
width as the OpenMM CPU thread count.

A 400-step, 4096-molecule trajectory matched the dump reader and
`seams chill-plus` on every label and every ring count. The same check
on one 32768-molecule frame matched as well. On a saved 4096-oxygen
frame the hop counts equal the ball counts: 1602, 1224, 2059, 2281
rings of size 3, 4, 5, and 6.

## Cost

Release library, 2 fs timestep. The 4096-molecule block is 2000 steps,
four repeats, intervals 100 and 500. The 32768-molecule block is 1000
steps, four repeats, one report every 100 steps. Bare time for 4096
molecules is 50.940 s (6.785 ns/day) on one thread and 12.206 s
(28.313 ns/day) on eight. For 32768 molecules it is 204.945 s
(0.843 ns/day) and 45.846 s (3.769 ns/day).

Clocks inside a report, in milliseconds:

| molecules | threads | interval | copy | neighbours | CHILL+ | rings | sum |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4096 | 1 | 100 | 0.40 | 1.66 | 1.42 | 7.33 | 10.81 |
| 4096 | 1 | 500 | 1.05 | 1.64 | 1.41 | 7.07 | 11.17 |
| 4096 | 8 | 100 | 0.41 | 4.03 | 2.65 | 2.99 | 10.08 |
| 4096 | 8 | 500 | 0.38 | 3.24 | 1.54 | 1.55 | 6.72 |
| 32768 | 1 | 100 | 2.74 | 11.65 | 11.98 | 63.75 | 90.12 |
| 32768 | 8 | 100 | 2.92 | 13.50 | 9.03 | 12.68 | 38.13 |

The ring search is the part that scales. On 4096 oxygens it falls from
7.33 ms on one thread to 2.99 ms on eight (2.4×) at interval 100, and
from 7.07 ms to 1.55 ms (4.6×) at interval 500. On 32768 oxygens it
falls from 63.75 ms to 12.68 ms (5.0×). Adding the 32768 eight-thread
sum onto the bare run is 0.83% of the step time. The measured
wall-clock means at eight threads are 27.665 ns/day for 4096 molecules
at interval 100 (14.3 ms per report, 2.3% slower than bare) and
3.704 ns/day for 32768 molecules (80.8 ms per report, 1.7% slower).
At one thread the reporter's added time is smaller than the spread of
the bare runs.

Set `SEAMS_LANG_LIB` to `libseams_lang` before constructing `InSitu`.
The reporter takes that engine, an unused dump name, the step interval,
and a type array with 1 on the atoms that enter the classifier.

```python
from dseams_insitu import DseamsReporter, InSitu

engine = InSitu(os.environ["SEAMS_LANG_LIB"])
simulation.reporters.append(DseamsReporter(engine, "unused", interval=6, types=oxygen))
```

`engine.frames` then holds the step, the oxygen coordinates, the cell,
the labels, the label counts, and the ring counts.

The 215 K mW frames are not part of the tree. `input/mw215/fetch.py`
writes them under `/var/tmp/dseams/mw215`, or under `MW215_DIR`.
A failed download is replaced by `regenerate.py`, which needs `lmp`
and writes the same directory. The deck is `in.equilibrate` with
velocity seed 215.
