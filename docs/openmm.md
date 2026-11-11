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
A positive skin keeps a list built at cutoff plus that skin while every
atom's minimum-image drift stays inside half the skin, and the bonds
are the pairs still inside 3.5 Å. The timed driver passes skin 0: over
100 steps the drift is past a skin that stays cheap to build, so each
report builds one cutoff list and both classifiers use it.

A 400-step trajectory, reported every 100 steps, matched the dump
reader and `seams chill-plus` on every label and every ring count.
A 0.02 Å drift reuses the wide list and still matches a fresh cutoff
list. A 2 Å drift rebuilds.

## Cost

Release library, classifier on one OpenMP thread, 5000 steps, six
repeats from one minimized state. Intervals 100, 500, and 1000.
OpenMM on 1 thread and on 8.

Bare time is 123.476 s (6.997 ns/day) on one thread and 30.368 s
(28.451 ns/day) on eight. The clocks inside a report, in milliseconds:

| threads | interval | copy | neighbours | CHILL+ | rings | sum |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 100 | 0.57 | 1.66 | 1.47 | 9.64 | 13.34 |
| 1 | 500 | 0.42 | 1.60 | 1.44 | 9.34 | 12.80 |
| 1 | 1000 | 0.38 | 1.58 | 1.40 | 9.30 | 12.67 |
| 8 | 100 | 0.30 | 1.38 | 1.31 | 9.29 | 12.29 |
| 8 | 500 | 0.38 | 1.49 | 1.38 | 9.36 | 12.61 |
| 8 | 1000 | 0.38 | 1.56 | 1.39 | 9.38 | 12.71 |

No report reused a neighbour list. Adding that sum to the bare time
slows one thread by 0.54%, 0.10%, and 0.051% at the three intervals,
and eight threads by 2.0%, 0.42%, and 0.21%.

The wall-clock means of the same six repeats:

| threads | interval | ns/day | ms per report | slowdown |
| --- | --- | --- | --- | --- |
| 1 | 100 | 6.789 | 75.78 | 2.98% |
| 1 | 500 | 6.968 | 51.48 | 0.42% |
| 1 | 1000 | 7.068 | −245.77 | −1.01% |
| 8 | 100 | 27.927 | 11.39 | 1.84% |
| 8 | 500 | 28.511 | −6.33 | −0.21% |
| 8 | 1000 | 28.218 | 50.19 | 0.82% |

At eight threads and interval 100 the wall gap is 11.4 ms per report,
next to the 12.3 ms inside the reporter. At one thread and interval 100
three of the six runs took 128–133 s against bare runs of 122–124 s,
so the wall mean sits above the phase sum. At intervals 500 and 1000
the added time is smaller than the spread of the bare runs.

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
