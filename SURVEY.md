# In-situ hooks for d-SEAMS 2.0

This note names common neighbour ice labelling (CHILL).
Nanoscale molecular dynamics (NAMD) is one engine in the comparison.
Assisted model building with energy refinement (AMBER) is another.
Visual molecular dynamics (VMD) is the trajectory viewer named beside
them.
The message passing interface (MPI) is the rank collective.
The Massachusetts Institute of Technology (MIT) licence covers i-PI.
The General Public License (GPL) covers LAMMPS, CP2K, and AMBER's
patched build.
The Lesser General Public License (LGPL) covers GROMACS and PLUMED.
The disk coordinate dump (DCD) is the OpenMM trajectory file.
The PLUMED load directive (LOAD) pulls an action into a run.
The interquartile range (IQR) is the third quartile minus the first.

We introduce a census on the running step.
In contrast with a dump that is read after the fact, the hook reads
positions the dynamics code already has.
d-SEAMS reports the primitive ring census and the published CHILL+
labels.
It also reports the hexagonal and double-diamond cage flags for the
chosen frames.
The same numbers must come back from reading the dumped trajectory
afterwards.
CHILL and CHILL+ stay the published classifiers.
The Atomic Simulation Environment is used only as a molecular dynamics
driver.
That use goes through its dynamics callback, and not through an
optimizer.

The architecture is one C core and thin adapters.
The core is `seams_insitu.h` in this repository.
The Meson options are `-Dlammps=enabled|auto|disabled` and the same
pattern for `gromacs`, `openmm`, and `plumed`.
The default is `auto`.
A build that finds none of those codes still compiles and its tests
still pass.
PLUMED is the broad fallback: one PLUMED action calls the core and every
engine that already calls PLUMED picks it up.

Device coordinates are part of the frame the adapters pass in.
This build does not dereference a device pointer.
Copying device memory to the host is left for a later run.

## Recommendation

Call `seams_insitu_step` or `seams_insitu_observe` once per analysed
step.
The call analyses the borrowed positions, box and types.
It does not write a trajectory and it does not read one back.
`debug_dump` appends a LAMMPS dump after that analysis and is off by
default.

The frame carries the step and the positions in ångström.
It also carries optional types and ids.
It also carries the LAMMPS edge box, which is the origin, the edge
lengths, and the tilt.
Pointers are borrowed for the call.
The engine copies them into the point cloud that CHILL+, the ring
search, and the cage search already read.
That copy is required by those published routines.
A single rank makes no second hidden copy when the adapter passes the
engine the buffer it already gathered.

A rank that sees only its local atoms sets `local_only` and leaves
`gathered` clear.
The core returns `SEAMS_INSITU_NEED_GATHER` and does not read the
pointer.
The adapter gathers the full system, because a ring and a CHILL+ label
need the neighbour shell.
The adapter then calls again with `gathered` set.
Message Passing Interface collectives stay in the adapter.
The core does not link that library.

A device pointer sets `gpu_resident`.
The core records the step and the device index and returns
`SEAMS_INSITU_GPU_DEFERRED`.
It does not load the pointer.

Cadence is `every`, `first`, and `last` on the step index.
A step outside that window returns `SEAMS_INSITU_SKIP` and does not
write a frame.

`seams_insitu_observe` appends one LAMMPS dump frame in round-trip
decimal.
The same call reads that frame back with the existing trajectory reader
and analyses those values.
`seams_insitu_replay` reads the same file and writes the same text.
The live record and the postprocessed record are the same bytes.
The text records the ring census by size.
It records the nine CHILL+ counts.
It records the cage counters and, when requested, the per-atom labels.

Adapters, in the order implemented here:

1. LAMMPS fix plugin `dseams`, loaded with `plugin load`.
2. GROMACS interactive molecular dynamics client.
3. OpenMM reporter.
4. PLUMED action `DSEAMS`, loaded with `LOAD`.
5. i-PI socket client and an Atomic Simulation Environment molecular
   dynamics observer.

The PLUMED action covers NAMD.
AMBER `sander` and AMBER `pmemd` use that same action.
CP2K, DL_POLY, and HOOMD-blue use it as well.
Other engines that already call PLUMED use it too.
Those engines do not need a separate compiled adapter for the census.

## Comparison

Each code states the positions the hook is documented to receive.
"Local" means a domain-decomposed rank sees its own atoms plus a halo.
The rank does not see the whole system.
"Copy" is relative to that buffer.

LAMMPS reads local atoms in the plugin callback.
The gather copies those arrays.
A GPU package keeps data on the device unless gathered.
The engine licence is GPLv2.

GROMACS reads local coordinates in the force provider.
Interactive molecular dynamics sends a float copy.
A GPU build keeps nonbonded data on the device.
The engine licence is LGPL-2.1-or-later.

OpenMM reads positions from a reporter state.
The reporter copies a host state.
One context does not split atoms across spatial ranks.
CUDA and OpenCL can keep positions on the device.
The engine licence is MIT and LGPL.

PLUMED reads positions through `plumed_cmd`.
The host copy uses PLUMED units.
The host gathers or passes local atoms.
The engine licence is LGPL.

i-PI reads the cell and the positions at double precision.
The socket copies them.
The server keeps the beads.
The engine licence is MIT on conda-forge.

The Atomic Simulation Environment reads positions on the dynamics
object.
The arrays are a NumPy view or a copy.
One process runs unless the calculator splits it.
The calculator may hold a device array.
The engine licence is LGPL-2.1.

The engine licence for NAMD is separate from the core.
NAMD reads requested atoms in ångström.
Tcl copies that set.
GPU builds exist.
This repository does not link NAMD.

AMBER reads coordinates through the PLUMED patch.
The patch copies coordinates into PLUMED.
`pmemd.cuda` keeps the main loop on the device.
This repository does not link AMBER.

CP2K reads the linked PLUMED library.
The host copies the coordinates.
MPI splits the domain.
GPU builds exist.
The engine licence is GPL.

HOOMD-blue reads a snapshot of positions.
The snapshot is a host copy on rank 0.
GPU arrays stay in the state until a snapshot.
The engine licence is BSD-3-Clause.

DL_POLY reads the linked PLUMED library.
The host copies the coordinates.
The device path depends on the build.
Use PLUMED for the census.

## LAMMPS

LAMMPS loads a shared object with the `plugin` command and calls
`lammpsplugin_init`.
A fix factory is registered on `creator.v2` because a fix takes the
LAMMPS pointer.
The factory also takes the argument count and the argument strings.
The same binary can be extended without rebuilding LAMMPS.
The package variant that provides this command is `plugin`.
That variant has been in the Spack recipe since the 8 April 2021
version.
See [Writing plugins](https://docs.lammps.org/Developer_plugins.html)
and [plugin command](https://docs.lammps.org/plugin.html).

`fix external` is the library-interface hook.
In `pf/callback` mode the fix calls the driver every `Ncall` steps with
the timestep.
The call also passes the local atom count, the local ids, the local
coordinates, and a force array.
`lammps_set_fix_external_callback` installs that function.
The coordinates are the atoms on that processor, so a spatial
decomposition is visible.
See [fix external](https://docs.lammps.org/fix_external.html) and
[Library utility functions](https://docs.lammps.org/Library_utility.html).

`fix python/invoke` calls a Python function at
`post_force` or `end_of_step` every N steps and the `python` command
loads the function.
Both need the PYTHON package.
See [fix python/invoke](https://docs.lammps.org/fix_python_invoke.html)
and [python command](https://docs.lammps.org/python.html).

The adapter in this repository is a fix style named `dseams`.
At the end of the step it packs the atoms in the fix group, gathers them
across ranks, and calls the core.
The gather runs for one rank as well, so the serial test executes that
path.
The plugin source includes LAMMPS headers.
LAMMPS is GPLv2.
An interface that uses LAMMPS source is distributed under GPLv2.
A separate library it calls need not be, unless the distributed binary
requires that library.
The core stays MIT.
The plugin file is GPL-2.0-only, matching the conda-forge package.
See
[LAMMPS open-source license](https://docs.lammps.org/Intro_opensource.html)
and the
[conda-forge lammps feedstock](https://github.com/conda-forge/lammps-feedstock).

GPU packages in LAMMPS can leave coordinates on the device.
The fix reads the host array `atom->x`.
A device-resident path would set `gpu_resident` and is not taken in this
build.

## GROMACS

Anything in `mdrun` beyond bonded forces, nonbonded forces, and
integration is an MDModule.
A module implements `IMDModule` and, typically, a force provider called
every step with the local coordinates.
Colvars and PLUMED are released as such modules.
The module page is
[mdrun modules](https://manual.gromacs.org/current/doxygen/html-full/page_mdmodules.xhtml).
The class page is
[MDModules](https://manual.gromacs.org/documentation/2025-current/doxygen/html-full/classgmx_1_1MDModules.xhtml).
One example is
[plumedMDModule.cpp](https://manual.gromacs.org/documentation/current/doxygen/html-full/plumedMDModule_8cpp.xhtml).

gmxapi exposes `Context.add_mdmodule` so a Python driver can add a
module to a simulation.
See the
[gmxapi Python reference](https://manual.gromacs.org/current/gmxapi/userguide/pythonreference.html).

Colvars is selected in the mdp file with `colvars-active` and
`colvars-configfile` on GROMACS 2024 and later.
See the
[Colvars GROMACS manual](https://colvars.github.io/gromacs-2024/colvars-refman-gromacs.html).

Interactive molecular dynamics speaks the VMD protocol on a TCP port.
`mdrun -imdport` sends positions to a client.
`-imdwait` pauses until a client connects.
The mdp key `IMD-group` selects the atoms, and `System` sends the whole
system.
Coordinates in that stream are a float copy, not a zero-copy view of the
domain.
The client in this repository analyses those values and writes them with
the core's dump.
So the replay sees the same numbers.
Domain decomposition still applies inside `mdrun`.
The main rank assembles the interactive molecular dynamics group before
the send.
See
[Interactive Molecular Dynamics](https://manual.gromacs.org/current/reference-manual/special/vmd-imd.html)
and
[gmx mdrun](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).

GROMACS is LGPL-2.1-or-later on conda-forge.
See the
[gromacs feedstock](https://github.com/conda-forge/gromacs-feedstock).
Spack packages it as well:
[Spack gromacs](https://packages.spack.io/package.html?name=gromacs).

## OpenMM

`CustomForce` classes evaluate a user expression on bonds, angles,
torsions, or nonbonded pairs.
They adjust the force field.
They do not report a census.
Global parameters can change during a simulation through
`Context.setParameter`.
See
[Custom Forces](https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html).

Reporters are the analysis hook.
A reporter implements `describeNextReport` and `report`.
The simulation passes a `State` that can include positions.
`StateDataReporter`, `DCDReporter`, and `PDBReporter` are the included
examples.
See
[Application Layer](https://docs.openmm.org/latest/api-python/app.html).
The adapter here is a reporter.
It converts the box vectors to the LAMMPS edge frame and calls the core.
OpenMM's CUDA and OpenCL platforms keep positions on the device until a
state copy.
The reporter asks for host positions.
A later run can set `gpu_resident` instead of that copy.

A C++ plugin exports `registerPlatforms` and `registerKernelFactories`
and is loaded with `Platform::loadPluginLibrary`.
A force that only forwards positions to an external library can subclass
`CustomCPPForceImpl` and implement `computeForce` without a GPU kernel.
See
[Writing Plugins](https://docs.openmm.org/latest/developerguide/03_writing_plugins.html)
and
[CustomCPPForceImpl](https://docs.openmm.org/latest/developerguide/09_customcppforceimpl.html).

OpenMM is MIT and LGPL.
See the [OpenMM repository](https://github.com/openmm/openmm/).
Spack packages OpenMM.

## PLUMED

`LOAD FILE=extensions.so` registers actions from a shared object.
`plumed mklib` builds that object.
From PLUMED 2.10 the load applies to commands that follow it, and it is
local to one Plumed object.
See
[LOAD](https://www.plumed.org/doc-v2.10/user-doc/html/_l_o_a_d.html).

An action implements `registerKeywords`, a constructor from
`ActionOptions`, and `calculate`.
`PLUMED_REGISTER_ACTION` records the directive.
See
[How to add a new function](https://www.plumed.org/doc-v2.10/developer-doc/html/_adding_a_function.html)
and
[ActionRegister](https://www.plumed.org/doc-v2.10/developer-doc/html/class_p_l_m_d_1_1_action_register.html).

The host calls the C API (`setStep`, `setPositions`, `setBox`, `calc`)
each step.
The action reads positions in nanometres and scales them to ångström
before the core.
The action in this repository is `DSEAMS`.

That DSEAMS action sits in the PLUMED project, which lists the engines
that already call PLUMED.
The list includes AMBER `pmemd` from version 20 and AMBER `sander` from
AmberTools 15.
It also includes CP2K, DL_POLY, GROMACS, LAMMPS.
NAMD, OpenMM through openmm-plumed, i-PI and the Atomic Simulation
Environment.
See the
[PLUMED manual](https://www.plumed.org/doc-v2.10/user-doc/html/index.html).
The project describes the licence as L-GPL on its
[home page](https://www.plumed.org/).
Spack packages PLUMED.

PLUMED therefore needs one action for every patched engine.
That action is the fallback this design uses.

## i-PI

i-PI is a server.
A client connects over an internet or Unix socket.
The server sends the cell and the positions.
The client returns energy, forces, and the virial.
The messages are a small set of 12-byte headers and native doubles.
See [Program Overview](https://docs.ipi-code.org/introduction.html).

The client in this repository is that force loop.
On each `POSDATA` it calls the core, then returns forces from a
short-ranged pair potential.
So the bead actually moves.
It does not call an optimizer.
conda-forge distributes i-PI 3.3.1 as a noarch Python package under the
MIT licence (`i-pi`).

The Atomic Simulation Environment implements the same socket as a
calculator, for codes that already speak it.
See
[socket I/O](https://docs.ase-lib.org/ase/calculators/socketio/socketio.html).
That calculator is not the d-SEAMS hook.
The hook is the client loop above, and, separately, the molecular
dynamics observer below.

## Atomic simulation environment

`MolecularDynamics.attach` runs the callback every n steps.
A dynamics object such as `Langevin` is the usual caller.
A trajectory writer is attached the same way.
See [Molecular dynamics](https://docs.ase-lib.org/ase/md.html).

The observer in this repository attaches only to a molecular dynamics
object.
It refuses classes whose name is an optimizer.
It reads `atoms.get_positions()` and the cell.
Converts the cell with `seams_insitu_box_from_rows`, and calls the core.
ASE is LGPL-2.1 on conda-forge.
See the
[ase feedstock](https://github.com/conda-forge/ase-feedstock/tree/22b3d653b26c31e75b61792570a425442baa1833).

## Other engines that call PLUMED

NAMD runs Tcl on the fly.
`tclForces on` and `tclForcesScript` define `calcforces`.
Inside it, `getstep` is the step.
`loadcoords` fills coordinates in ångström for atoms requested earlier.
The set is whatever the script asked for.
Which is a poor fit for a full-box ring census, and the coordinates are
a copy.
See
[Tcl Forces and Analysis](https://www.ks.uiuc.edu/Research/namd/3.0/ug/node50.html).
NAMD also turns PLUMED on with `plumed on` and `plumedfile`.
The user guide notes that the virial from PLUMED is not applied, so the
safe runs are constant volume.
See
[namd-2.14](https://www.plumed.org/doc-v2.10/user-doc/html/namd-2-14.html).
This repository does not link NAMD.

AMBER `sander` takes `plumed=1` and `plumedfile='plumed.dat'` in the
`&cntrl` namelist.
The same keys are used with `pmemd`, including the CUDA build, when
AMBER was compiled with PLUMED.
See the
[PLUMED manual](https://www.plumed.org/doc-v2.10/user-doc/html/index.html)
and the
[AMBER mailing list note](http://archive.ambermd.org/202006/0036.html).
This repository does not link AMBER.

CP2K calls PLUMED from `&METADYN` with
`USE_PLUMED .TRUE.` and `PLUMED_INPUT_FILE`.
The interface is in CP2K from version 2.7 and does not need a source
patch for PLUMED 2.
See
[How to compile CP2K with PLUMED](https://www.cp2k.org/howto:install_with_plumed).
CP2K is GPL.
The PLUMED action is the adapter.

HOOMD-blue runs a user action from `hoomd.update.CustomUpdater`.
`hoomd.write.CustomWriter`, or `hoomd.tune.CustomTuner`.
`act(timestep)` reads `self._state`.
A snapshot on rank 0 carries positions.
The documentation shows the communicator rank.
So a domain decomposition is visible, and a snapshot is a host copy.
See
[Action](https://hoomd-blue.readthedocs.io/en/stable/hoomd/custom/action.html).
conda-forge packages HOOMD-blue 7.2.0 under BSD-3-Clause.
A writer action would call the core.
It is not built in this run.
The PLUMED list does not name HOOMD-blue, so the Python action is the
direct hook.

DL_POLY 5 enables PLUMED with the `plumed` keyword in the CONTROL file.
The input defaults to a file named `PLUMED`.
The binary must be linked with PLUMED.
See
[DL_POLY PLUMED](https://ccp5.gitlab.io/dl-poly/UserManual/force-fields/plumed.html).

## Build and packaging

conda-forge publishes LAMMPS, GROMACS, OpenMM, PLUMED, ASE, i-PI, and
HOOMD-blue.
The live tests in this run install from that channel.
Spack publishes LAMMPS (with a `plugin` variant and a `plumed` variant),
GROMACS, OpenMM, and PLUMED.
See [Spack lammps](https://packages.spack.io/package.html?name=lammps).
See
[Spack gromacs](https://packages.spack.io/package.html?name=gromacs).

ASE and i-PI are published on conda-forge as pure Python, so a pip
install of those projects is the same code.
OpenMM also publishes wheels.
LAMMPS, GROMACS, PLUMED, CP2K, NAMD, AMBER and DL_POLY are native
builds.
Meson feature options keep this repository optional around them.
`auto` builds the LAMMPS plugin when `lammpsplugin.h` compiles.
The same setting builds the PLUMED action when pkg-config finds
`plumed`.
`enabled` fails configure when that header or that package is missing.
`disabled` skips the search.
A conda-forge LAMMPS package often omits headers.
`-Dlammps_include` then points at a LAMMPS source tree.
`seams-gromacs-imd` is the GROMACS client.
It links no GROMACS library.
A fake socket is the default test, so the client is always built.
A GROMACS build with PLUMED calls the same action through
`mdrun -plumed`.
OpenMM, i-PI and the Atomic Simulation Environment adapters are Python.
Meson does not compile them.
Gathering in the LAMMPS adapter uses MPI_Gatherv onto one rank before
the analysis.
Building always includes the core library.
Building always includes `seams-insitu replay`.

## What this run builds

The results are the overhead quartiles and the census counts.
We implement the core.
The pieces this run builds follow.
The pieces are the replay tool and the LAMMPS plugin.
They include the GROMACS client, the OpenMM reporter, and the PLUMED
action.
They include the i-PI client and the Atomic Simulation Environment
observer.
Analysis reads the borrowed frame.
A debug dump stays off unless requested.
Reading that dump back is not part of the run.
Gathering owned atoms uses `MPI_Gatherv` inside the LAMMPS plugin.
It also gathers the neighbour rows.
Reuse of that list continues while atoms stay inside the skin.
A cutoff shell stays with the core.
Rebuilding starts when an atom moves by more than half of a 0.5 ångström
skin.
Rings and cages update only for the rings that touch a moved atom.
A large set uses one full pass.
The full pass matches the incremental result and finishes sooner.
The PLUMED action passes the position pointer with a coordinate scale.
The action does not copy the array.
The engine builds the cutoff shell.
`NeighborList::getIndexPair` in PLUMED 2.10.2 computes `8*ii` in a
32-bit unsigned integer.
The product wraps from 32769 atoms upward.
The memory guard above that line only refuses lists past about 50000
atoms.
A 40000-atom list walks off the position array.
`tests/insitu/plumed_neighborlist_overflow.cpp` evaluates that
expression.
GROMACS 2026.3 passed the MPI communicator into PLUMED by value.
Open MPI stores a communicator as a pointer.
The call read the wrong address.
The matching build passes the address of a local communicator.
PLUMED and GROMACS are built against the same Open MPI.
From PLUMED 2.10 the action reports that it applies no forces.
Only rank 0 writes the census.
Equality uses each code's own trajectory file through the `seams`
command line.
The comparison does not use a dump the engine wrote.
Overhead uses at least six repeats on one host.
Warmups stay apart from the samples.
The report uses type 7 quartiles.
OpenMP thread counts are the powers of two up to the processors the
process may use.
The LAMMPS route kept for use is the runtime plugin.
An in-tree package of the same fix did not separate from the plugin once
both were built at -O3.
The GROMACS census comes from a d-SEAMS module inside mdrun.
The molecular dynamics parameter keys are dseams-active.
Dseams-input-group, dseams-nstcalc, dseams-cutoff and dseams-output.
The OpenMM reporter reads positions with getState.
It asks for positions and a periodic box.
It does not write a trajectory.

## Other

These manuals are the ones this note opened.

- [Writing plugins](https://docs.lammps.org/Developer_plugins.html)
- [plugin command](https://docs.lammps.org/plugin.html)
- [fix external](https://docs.lammps.org/fix_external.html)
- [Library utility functions](https://docs.lammps.org/Library_utility.html)
- [fix python/invoke](https://docs.lammps.org/fix_python_invoke.html)
- [python command](https://docs.lammps.org/python.html)
- [LAMMPS open-source license](https://docs.lammps.org/Intro_opensource.html)
- [conda-forge lammps feedstock](https://github.com/conda-forge/lammps-feedstock)
- [mdrun modules](https://manual.gromacs.org/current/doxygen/html-full/page_mdmodules.xhtml)
- [MDModules](https://manual.gromacs.org/documentation/2025-current/doxygen/html-full/classgmx_1_1MDModules.xhtml)
- [gmxapi Python reference](https://manual.gromacs.org/current/gmxapi/userguide/pythonreference.html)
- [Colvars GROMACS manual](https://colvars.github.io/gromacs-2024/colvars-refman-gromacs.html)
- [plumedMDModule.cpp](https://manual.gromacs.org/documentation/current/doxygen/html-full/plumedMDModule_8cpp.xhtml)
- [Interactive Molecular Dynamics](https://manual.gromacs.org/current/reference-manual/special/vmd-imd.html)
- [gmx mdrun](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html)
- [conda-forge gromacs feedstock](https://github.com/conda-forge/gromacs-feedstock)
- [Custom Forces](https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html)
- [Application Layer](https://docs.openmm.org/latest/api-python/app.html)
- [Writing Plugins](https://docs.openmm.org/latest/developerguide/03_writing_plugins.html)
- [CustomCPPForceImpl](https://docs.openmm.org/latest/developerguide/09_customcppforceimpl.html)
- [OpenMM repository](https://github.com/openmm/openmm/)
- [LOAD](https://www.plumed.org/doc-v2.10/user-doc/html/_l_o_a_d.html)
- [How to add a new function](https://www.plumed.org/doc-v2.10/developer-doc/html/_adding_a_function.html)
- [ActionRegister](https://www.plumed.org/doc-v2.10/developer-doc/html/class_p_l_m_d_1_1_action_register.html)
- [PLUMED manual](https://www.plumed.org/doc-v2.10/user-doc/html/index.html)
- [PLUMED](https://www.plumed.org/)
- [namd-2.14](https://www.plumed.org/doc-v2.10/user-doc/html/namd-2-14.html)
- [Program Overview](https://docs.ipi-code.org/introduction.html)
- [socket I/O](https://docs.ase-lib.org/ase/calculators/socketio/socketio.html)
- [Molecular dynamics](https://docs.ase-lib.org/ase/md.html)
- [ase feedstock](https://github.com/conda-forge/ase-feedstock/tree/22b3d653b26c31e75b61792570a425442baa1833)
- [Tcl Forces and Analysis](https://www.ks.uiuc.edu/Research/namd/3.0/ug/node50.html)
- [AMBER mailing list note](http://archive.ambermd.org/202006/0036.html)
- [How to compile CP2K with PLUMED](https://www.cp2k.org/howto:install_with_plumed)
- [Action](https://hoomd-blue.readthedocs.io/en/stable/hoomd/custom/action.html)
- [hoomd on conda-forge](https://anaconda.org/conda-forge/hoomd)
- [DL_POLY PLUMED](https://ccp5.gitlab.io/dl-poly/UserManual/force-fields/plumed.html)
- [Spack lammps](https://packages.spack.io/package.html?name=lammps)
- [lammps package.py](https://raw.githubusercontent.com/spack/spack-packages/refs/heads/develop/repos/spack_repo/builtin/packages/lammps/package.py)
- [Spack gromacs](https://packages.spack.io/package.html?name=gromacs)
- [openmm package.py](https://raw.githubusercontent.com/spack/spack-packages/refs/heads/develop/repos/spack_repo/builtin/packages/openmm/package.py)
- [plumed package.py](https://raw.githubusercontent.com/spack/spack-packages/refs/heads/develop/repos/spack_repo/builtin/packages/plumed/package.py)
