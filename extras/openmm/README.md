# OpenMM reporter

In-run neighbours, primitive rings, CHILL, and CHILL+ from `libseams_lang`.

Put `adapters/python` on `PYTHONPATH`. Point `SEAMS_LIB` at
`libseams_lang.so` (`bbdir/bindings/libseams_lang.so` after a Meson build).

```python
from dseams_insitu import DseamsReporter, Engine

with Engine(library, "census.txt", cutoff=3.5) as eng:
    sim.reporters.append(
        DseamsReporter(eng, "frames.lammpstrj", 10, indices=oxygen)
    )
    sim.step(100)
```

The reporter reads positions from the context, on the host, in angstrom,
with the periodic box enforced. `indices` keeps the oxygen atoms of a
TIP4P model. A one-site model omits `indices` and analyses every
particle. The cutoff is angstrom. The bond graph is the cutoff list
`seams chill` uses.

The dump holds those atoms, type 1. On each shared snapshot, these
commands print the lines stored in `census.txt`:

```
seams neighbours frames.lammpstrj --cutoff 3.5 --type 1 --graph cutoff --frame N --jobs 1
seams rings frames.lammpstrj --cutoff 3.5 --type 1 --graph cutoff --frame N --jobs 1
seams chill frames.lammpstrj --cutoff 3.5 --type 1 --graph cutoff --frame N --jobs 1
seams chill-plus frames.lammpstrj --cutoff 3.5 --type 1 --graph cutoff --frame N --jobs 1
```

`N` is the snapshot index in the dump, starting at 1. `NO_COLOR=1` keeps
the ice names as plain text. The census prefixes the CHILL line with
`chill` and the CHILL+ line with `chill-plus`.

The third argument of `DseamsReporter` is the stride in steps. The
OpenMM howto records the call, the match against `seams`, and the
overhead at strides of 10, 100, and 1000 on a few hundred waters and
on about four thousand. The test builds both boxes when it runs and
leaves the dumps out of git. See `data/README`.
