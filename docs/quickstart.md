# Quickstart

Build the engine, then classify the cubic cell shipped in this tree.
The Sphinx page is `docs/orgmode/quickstart.org`. The same steps are below.

pixi:

```bash
pixi run setup && pixi run build
./bbdir/seams chill-plus input/traj/mW_cubic.lammpstrj --type 1
```

Compilers and Eigen already on `PATH`:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
./build/seams chill-plus input/traj/mW_cubic.lammpstrj --type 1
```

Both print `nop 4096 cubic 4096`.

`libseams.so` exports `seams.h`. The language fronts and the engine hooks link
that library and do not include the C++ headers.

```bash
export SEAMS_LIB_DIR=$PWD/build/src
export LD_LIBRARY_PATH=$SEAMS_LIB_DIR:$LD_LIBRARY_PATH
export SEAMS_LIB=$SEAMS_LIB_DIR/libseams.so
```

Nix, conda, pip, Spack, and an EESSI toolchain are on the Sphinx install page.
`pip install pydseamslib` installs the Frame package (`import pydseams`).
It does not build `seams` or `import dseams`.

Pages for each front, with the 2.14.0 build command and the test CI runs:

| front | install | tutorial |
| --- | --- | --- |
| C | [api/c.md](api/c.md) | [tutorials/c.md](tutorials/c.md) |
| C++ | [api/cpp.md](api/cpp.md) | [tutorials/cpp.md](tutorials/cpp.md) |
| Fortran | [api/fortran.md](api/fortran.md) | [tutorials/fortran.md](tutorials/fortran.md) |
| Rust | [api/rust.md](api/rust.md) | [tutorials/rust.md](tutorials/rust.md) |
| Python | [api/python.md](api/python.md) | [tutorials/python.md](tutorials/python.md) |
| Julia | [api/julia.md](api/julia.md) | [tutorials/julia.md](tutorials/julia.md) |
| R | [api/r.md](api/r.md) | [tutorials/r.md](tutorials/r.md) |
| Lua / Fennel | [api/lua.md](api/lua.md) | [tutorials/lua.md](tutorials/lua.md) |
| GROMACS | [engines/gromacs.md](engines/gromacs.md) | [tutorials/gromacs.md](tutorials/gromacs.md) |
| LAMMPS | [engines/lammps.md](engines/lammps.md) | [tutorials/lammps.md](tutorials/lammps.md) |
| OpenMM | [engines/openmm.md](engines/openmm.md) | [tutorials/openmm.md](tutorials/openmm.md) |
| PLUMED | [engines/plumed.md](engines/plumed.md) | [tutorials/plumed.md](tutorials/plumed.md) |
| in-situ | [engines/insitu.md](engines/insitu.md) | [tutorials/insitu.md](tutorials/insitu.md) |

Then run one front:

```bash
python -m dseams
cargo run --manifest-path bindings/rust/Cargo.toml --example ice
gfortran -Ibindings/fortran bindings/fortran/test_seams.f90 -L$SEAMS_LIB_DIR -lseams
julia --project=bindings/julia/Dseams bindings/julia/Dseams/test/runtests.jl
```

The eight-atom programs are `docs/examples/quickstart.py` and
`docs/examples/quickstart.c`. `meson test docs_quickstart_py` and
`meson test docs_quickstart_c` both expect 8 cubic molecules.

A frame is a coordinate array and a cell. Orthorhombic boxes pass three
lengths. Triclinic boxes pass the nine LAMMPS numbers. CHILL+ accepts the
published rule and the OVITO rule. q4 and q6 are the local Steinhardt
parameters. Rings are primitive rings with a symmetry label on each ring.
