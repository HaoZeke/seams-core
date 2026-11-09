# In-situ analysis

d-SEAMS 2.13.0 labels one frame while a molecular dynamics calculation
runs.
The live call goes to the library.
It does not call PLUMED.
Chemfiles reads the saved trajectory.
The con format stays on readcon.

## Library

Build seams-core with the default options.
The shared library exports `seams_insitu_create`,
`seams_insitu_observe`, and `seams_insitu_destroy`.
Pass positions in Ångström.
Set `coord_scale` to 10 when the engine stores nanometres.

## Command

Postprocess the dump with one thread.
`seams rings` takes the dump, cutoff 3.5, type 1, frame 1, and the last
frame.
`seams chill-plus` takes the same arguments.
`seams cages` adds the cutoff graph.
The last frame is the number of frames in the dump.

## LAMMPS plugin

The checked binary is LAMMPS (22 Jul 2025 - Update 6).
Compile `libdseams_lammps.so` against that source and against
`libyodaLib`.
The input loads the object and adds the fix:

    plugin load libdseams_lammps.so
    fix sea all dseams record census.txt every 1 cutoff 3.5

Write a dump of the same atoms on the same stride.
`seams` on that dump must match `census.txt` on every shared step.

## LAMMPS package

The same fix compiles inside LAMMPS.
Configure with `-DPKG_DSEAMS=yes`, and point `SEAMS_INCLUDE` and
`SEAMS_LIBRARY` at seams-core.
The input omits `plugin load`.
The fix line is the one above.
The checked package binary reports the same LAMMPS version string.

## GROMACS

The module lives on the `dseams-insitu` branch, beside the metatomic
module.
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
