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

The machine is an Intel Xeon, 8 cores, one thread per core. The
release library is built with OpenMP. One process runs the blocks one
after another. Each block is six repeats from one minimized state.
The OpenMP width is the OpenMM thread count, except the cutoff list
itself stays on one thread: waking the team for that list, then again
for the rings, cost more than the list. CHILL+ joins the team from
16384 atoms up. An orthorhombic box with a team already asked for uses
the persistent CSR builder, whose cell tables stay allocated.

```bash
export OMP_NUM_THREADS=1
export OPENMM_TIP4P_STEPS=2000
export OPENMM_TIP4P_REPEATS=6
export OPENMM_TIP4P_INTERVALS=100
export OPENMM_TIP4P_LARGE_STEPS=1000
export OPENMM_TIP4P_LARGE_REPEATS=6
python3 scripts/openmm_tip4p_ice.py "$LIBSEAMS_LANG" "$SEAMS" run.log
```

`$LIBSEAMS_LANG` is the release `libseams_lang` and `$SEAMS` is the
release `seams`. The timestep is 2 fs. 4096 molecules run for 2000
steps. 32768 molecules run for 1000 steps. A report is every 100 steps.

Bare time is 48.616 s (7.109 ns/day) and 12.195 s (28.339 ns/day) for
4096 molecules, and 203.116 s (0.851 ns/day) and 48.178 s
(3.587 ns/day) for 32768 molecules, on one thread and on eight.

Clocks inside a report, in milliseconds:

| molecules | threads | copy | neighbours | CHILL+ | rings | sum |
| --- | --- | --- | --- | --- | --- | --- |
| 4096 | 1 | 0.37 | 1.61 | 1.36 | 7.03 | 10.37 |
| 4096 | 8 | 0.36 | 1.10 | 1.35 | 6.21 | 9.01 |
| 32768 | 1 | 3.03 | 12.19 | 11.98 | 64.00 | 91.19 |
| 32768 | 8 | 2.63 | 8.30 | 9.51 | 13.52 | 33.96 |

On eight threads the neighbour list is no longer the slow case: 1.10 ms
against 1.61 ms at 4096, and 8.30 ms against 12.19 ms at 32768.
CHILL+ at 4096 is 1.35 ms against 1.36 ms. The ring search on 32768
oxygens falls from 64.00 ms to 13.52 ms.

Wall-clock means of the same six repeats, report every 100 steps:

| molecules | threads | reporter ns/day | ms per report | slowdown |
| --- | --- | --- | --- | --- |
| 4096 | 1 | 7.091 | 6.08 | 0.25% |
| 4096 | 8 | 28.148 | 4.14 | 0.67% |
| 32768 | 1 | 0.867 | −379.03 | −1.90% |
| 32768 | 8 | 3.586 | 1.57 | 0.033% |

At eight threads and a report every 100 steps the reporter is 0.67%
slower on 4096 molecules and 0.033% slower on 32768. The one-thread
32768 gap is inside the spread of the bare runs.

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
