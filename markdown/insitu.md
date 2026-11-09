# In-situ analysis

d-SEAMS 2.13.0 labels one frame while a molecular dynamics calculation
runs.
The live call goes to the library.
It does not call PLUMED.
Chemfiles reads the saved trajectory.
The con format stays on readcon.

## Library

Build seams-core 2.13.0 with Meson.
Chemfiles reads every postprocess format except con.
The Meson build includes Chemfiles.

Run meson setup for bbdir-rel with optimization 2, then meson compile.
The build writes `libyodaLib.so`.
A caller creates an engine, passes one frame, and destroys the engine.

    seams_insitu_options opt;
    seams_insitu_options_default(&opt);
    seams_insitu_engine *eng = seams_insitu_create(&opt);
    seams_insitu_frame frame = {};
    frame.step = 0;
    frame.n = natoms;
    frame.xyz = xyz;
    frame.types = types;
    frame.coord_scale = 1.0;
    seams_insitu_observe(eng, &frame, nullptr);
    seams_insitu_destroy(eng);

Pass positions in Ångström.
Set `coord_scale` to 10 when the engine stores nanometres.
The header is `seams_insitu.h`.

## Python

The package is `adapters/python`.
Put that directory on `PYTHONPATH`.
Point `SEAMS_LIB` at `libyodaLib.so`.

    from dseams_insitu import Engine
    with Engine("bbdir-rel/src/libyodaLib.so", "census.txt") as eng:
        eng.observe(step, positions, box)

`positions` is an array of x, y, z in Ångström.
`box` is the origin, the edge lengths, and the tilt.

## Command

Postprocess the dump with one thread.
`seams rings` takes the dump, cutoff 3.5, type 1, frame 1, and the last
frame.
`seams chill-plus` takes the same arguments.
`seams cages` adds the cutoff graph.
The last frame is the number of frames in the dump.

## LAMMPS plugin

The checked binary is LAMMPS (22 Jul 2025 - Update 6).
Compile `libdseams_lammps.so` against that source tree and against
`libyodaLib`.

    cmake -S plugin -B plugin-build \
      -DLAMMPS_SRC=lammps/src \
      -DSEAMS_LIBRARY=bbdir-rel/src/libyodaLib.so
Compile `plugin-build` after that configure line.

The input loads the object and adds the fix:

    plugin load libdseams_lammps.so
    fix sea all dseams record census.txt every 1 cutoff 3.5

Write a dump of the same atoms on the same stride.
`seams` on that dump must match `census.txt` on every shared step.

## LAMMPS package

The same fix compiles inside LAMMPS.
The package sources are `src/DSEAMS`.

    cmake -S lammps/cmake -B lammps-build \
      -DPKG_DSEAMS=yes \
      -DSEAMS_INCLUDE=seams-core/src/include/internal \
      -DSEAMS_LIBRARY=bbdir-rel/src/libyodaLib.so
Compile `lammps-build` after that configure line.

The input omits `plugin load`.
The fix line is the one above.
The checked package binary reports the same LAMMPS version string.
One rank is the acceptance scope for the package.
The plugin acceptance also runs on two ranks.

## GROMACS

The module lives on the `dseams-insitu` branch, beside the metatomic
module.
GROMACS main is pinned at 2027.0-dev.
The commit is 1ede0115a1c6cc063813d2b13085a86640e4273e.
Configure leaves the module out when either path is empty.

The GROMACS configure line sets `GMX_MPI` to ON.
It passes the same include directory and the same library path.
Compile that build directory next.

`gmx_mpi --version` on the checked build prints
`2027.0-dev-20261108-d067d3f691`.
Add these lines to the mdp file:

    dseams-active = yes
    dseams-input-group = System
    dseams-nstcalc = 1
    dseams-cutoff = 3.5
    dseams-output = census.txt

The cutoff is in Ångström.
The engine stores coordinates in nanometres, and the module multiplies
them by ten before the analysis.
Keep `nstxout` equal to `dseams-nstcalc` so the trajectory has the same
steps as the census.
Run `gmx_mpi grompp` and then `gmx_mpi mdrun` on one rank.

## OpenMM

The checked import is OpenMM 8.6.1.
The reporter is `dseams_insitu.DseamsReporter`.

    from dseams_insitu import DseamsReporter, Engine
    with Engine(library, "census.txt") as eng:
        sim.reporters.append(DseamsReporter(eng, None, 1))
        sim.step(6)

The reporter reads host positions from the context.
Save those positions to a dump.
Postprocess that dump with `seams`.

## Match

On a shared step the ring line matches.
The ice counts match.
The cage counts match.
The comparison checks every step in that list.

## Acceptance driver

The acceptance driver `tests/insitu/engine_e2e.py` runs the three
engines.
Set `SEAMS_LIB`, `SEAMS_CLI`, `SEAMS_INSITU_EXE` and
`SEAMS_GRO_IMD_EXE`.
Also set `DSEAMS_MD_PREFIX`, `DSEAMS_LMP`, `DSEAMS_MPIRUN` and
`DSEAMS_PLUGIN`.
Also set `DSEAMS_LMP_INTREE`, `DSEAMS_GMX` and `DSEAMS_GMX_MPIRUN`.
The driver requires the live step list to equal the trajectory step
list.
The driver fails when a frame is missing, extra, or repeated.
OpenMM writes a disk coordinate dump (DCD).
The driver fails when the census of that dump mismatches.
