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
CPU with one thread. Twelve steps are taken, and the reporter runs
every sixth step. Six bare runs and six runs with the reporter share
one minimized state. The third argument of the driver is the log.

CHILL+ uses the cutoff neighbour list at 3.5 Å, type 1, on the oxygen
sites only. Rings are the primitive rings of that same list through size
6. A frame matches when every label and every ring count from the
in-situ arrays equals the result of reading the written dump, and the
`seams chill-plus` histogram equals those labels.

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
