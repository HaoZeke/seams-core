Copyright (c) d-SEAMS developers

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
The keywords are `record`, `every`, `cutoff`, `depth`, `host`, and `neigh`.
`host yes` is the default and keeps the LAMMPS neighbour list.
Without `neigh`, the fix writes no neighbour file.

```
fix sea all dseams record live.txt neigh host.txt every 1 cutoff 3.5
```

`host.txt` holds one block per analysed step:

```
frame
step 0
cutoff 3.5
box xlo ylo zlo lx ly lz xy xz yz
pbc 0 0 0
host 1
atoms 4
atom id type x y z
row id neighbour neighbour
endframe
```

`host 0` means that step did not build a list. Each `row` starts with
the central atom id. The following ids are the neighbours kept for the
fix group. `pbc` is the LAMMPS boundary code on x, y, and z. Zero is
periodic. The nine box numbers are printed with `%.17g`.

The CHILL+ compute accepts `graph`, `cutoff`, and `file`. Without
`file`, it writes no label file.

```
compute ice all dseams file labels.txt cutoff 3.5
```

`labels.txt` holds one block per evaluation:

```
frame
step 0
cutoff 3.5
box xlo ylo zlo lx ly lz xy xz yz
pbc 0 0 0
atoms 4
label id type x y z code
endframe
```

`code` is 0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
5 interClathrate, 6 unclassified, 7 reCubic, 8 reHex. This compute
assigns the first seven from the CHILL+ bond table.

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
