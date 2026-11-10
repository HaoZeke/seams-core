# fix dseams

`in.dseams` reads `water.data` (four oxygen sites) with `pair_style zero`
and analyses timestep 0.

Place `libdseams_lammps.so` beside the input and run it:

```
plugin load ./libdseams_lammps.so
fix sea all dseams record live.txt every 1 cutoff 3.5 depth 8
dump native all custom 1 dump.lammpstrj id type x y z
run 0
```

The plugin is the object from `src/adapters/lammps`. The in-tree package
registers the same fix. A missing `./libdseams_lammps.so` does not abort
the input, so that package binary runs this file as written.

`live.txt` is the census. `seams` on `dump.lammpstrj` must match it.
The keywords are `record`, `every`, `cutoff`, `depth`, and `host`.
`host yes` is the default and keeps the LAMMPS neighbour list.

Build the plugin against the LAMMPS that will load it:

```
cmake -S src/adapters/lammps -B plugin-build \
  -DLAMMPS_SRC=lammps/src \
  -DSEAMS_INCLUDE=src/include/internal \
  -DSEAMS_LIBRARY=bbdir/src/libyodaLib.so
cmake --build plugin-build
cp plugin-build/libdseams_lammps.so extras/lammps/
```

The package build adds `-DPKG_DSEAMS=on` and the same `SEAMS_INCLUDE`
and `SEAMS_LIBRARY`. Its `lmp` runs `in.dseams` without the copy.
