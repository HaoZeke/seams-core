# LAMMPS

Meson builds `run_hook`, the in-run hook. It links `liblammps` and registers a
`fix external` callback. Every step the callback calls the same C API as
`seams chill-plus`. `engines/lammps/compute_dseams.cpp` is a per-atom style.
When the LAMMPS class headers are installed, as in conda-forge and EESSI,
Meson builds it as `dseamsplugin.so` for `plugin load`; otherwise copy it
into a LAMMPS source tree. On several MPI ranks the compute gathers every
atom to rank 0, classifies the whole cell once, and returns each rank its
own labels, so the result does not depend on the domain decomposition.

`in.dseams` is the deck. `ice.data` is 64 oxygens of cubic ice. Forces are
zeroed so the cell does not melt. The hook prints `overhead_s`, the mean
seconds spent inside d-SEAMS per step, and for the last step
`lammps hexagonal N six-rings M q6 Q`. `run_hook.cpp` includes `library.h`,
so the LAMMPS include directory is on the compile line. `meson test lammps_hook`
runs this deck, and `meson test docs_engine_lammps` runs the same hook on
`input/traj/ice1h.lammpstrj`. The copy-in compute style is
`extras/lammps/in.dseams`.

`compute dseams` takes `cutoff` (angstrom), `graph` (`cutoff` or `mutual`),
and `rule` (`published` or `ovito`). The deck's `fix ext` line is
`external pf/callback`.

`meson test lammps_e2e` loads the plugin into `lmp` (or `$LMP`), runs 200
steps of mW ice and liquid with `engines/lammps/in.mw.lmp`, and compares the
`c_ice` column of every dumped frame with `seams chill-plus --per-atom` on
the same dump. `lammps_e2e_np2` repeats it on two ranks.

```{literalinclude} ../../../engines/lammps/in.mw.lmp
```

```{literalinclude} ../../../engines/lammps/in.dseams
:lines: 1-15
```

```{literalinclude} ../../../extras/lammps/in.dseams
:lines: 11-13
```

```{literalinclude} ../../../engines/lammps/run_hook.cpp
:language: cpp
:lines: 42-46
```
