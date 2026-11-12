#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build libseams and run the C API, the four language fronts, and the
# engine decks. Engine binaries are linked with the system compiler so
# they load the distribution LAMMPS and GROMACS.
set -e
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$root/build
lib=$build/src
eessi_lib=/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/sapphirerapids/software/GCCcore/15.2.0/lib64

if [ ! -x "$lib/test_seams_api" ]; then
  meson setup "$build" --wrap-mode=nofallback -Dwith_tests=false \
    -Dwith_ira=disabled -Dwith_sphericart=disabled -Dwith_nauty=disabled \
    -Dwith_gpulite=disabled -Dwith_openmp_offload=disabled
  meson compile -C "$build" test_seams_api
fi
"$lib/test_seams_api"
python3 "$root/docs/build_docs.py"
python3 "$root/engines/make_ice.py"

# Fortran
mod=$(mktemp -d)
gfortran -J"$mod" -c "$root/bindings/fortran/seams.f90" -o "$mod/seams.o"
gfortran -I"$mod" "$root/bindings/fortran/test_seams.f90" "$mod/seams.o" \
  -o "$mod/test_seams" -L"$lib" -lseams -Wl,-rpath,"$lib"
"$mod/test_seams"

# Rust
SEAMS_LIB_DIR=$lib cargo test --manifest-path "$root/bindings/rust/Cargo.toml" --offline -- --test-threads=1

# Python (nanobind). SEAMS_PY is a virtualenv whose python can import numpy.
if [ -n "$SEAMS_PY" ]; then
  nb=$("$SEAMS_PY" -c 'import nanobind, pathlib; print(pathlib.Path(nanobind.__file__).parent)')
  inc=$("$SEAMS_PY" -c 'import sysconfig; print(sysconfig.get_path("include"))')
  suf=$("$SEAMS_PY" -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
  c++ -shared -fPIC -O2 -std=c++17 -fvisibility=hidden \
    -I"$root/include" -I"$nb/include" -I"$nb/ext/robin_map/include" -I"$inc" \
    "$root/bindings/python/dseams/_core.cpp" "$nb/src/nb_combined.cpp" \
    -L"$lib" -lseams -Wl,-rpath,"$lib" \
    -o "$root/bindings/python/dseams/_core$suf"
  (cd "$root/bindings/python" && PYTHONPATH="$PWD" LD_LIBRARY_PATH="$eessi_lib:$lib" \
    "$SEAMS_PY" -m unittest tests/test_dseams.py)
  (cd "$root/engines/openmm" && PYTHONPATH="$root/bindings/python" LD_LIBRARY_PATH="$eessi_lib:$lib" \
    "$SEAMS_PY" test_reporter.py)
fi

# LAMMPS in-run hook and the compute object.
"$root/engines/link_host.sh" -c -O2 -std=c++17 \
  -I/usr/include/lammps -I/usr/lib/x86_64-linux-gnu/openmpi/include -I"$root/include" \
  "$root/engines/lammps/compute_dseams.cpp" -o "$mod/compute_dseams.o"
"$root/engines/link_host.sh" -O2 -std=c++17 \
  -I/usr/include/lammps -I"$root/include" \
  "$root/engines/lammps/run_hook.cpp" -o "$mod/run_hook" \
  -L/usr/lib/x86_64-linux-gnu -Wl,--allow-shlib-undefined \
  -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,"$lib" \
  -llammps /usr/lib/x86_64-linux-gnu/libcurl.so.4 -L"$lib" -lseams
(cd "$root/engines/lammps" && env -u LIBRARY_PATH LD_LIBRARY_PATH="$lib" "$mod/run_hook" in.dseams)

# GROMACS trajectory tool. analyzeFrame is the per-frame hook.
gmx_cflags=$(PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig pkg-config --cflags libgromacs)
gmx_libs=$(PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig pkg-config --libs libgromacs)
env -i PATH=/usr/bin:/bin /usr/bin/g++ -O2 -std=c++17 "$root/engines/gromacs/seams_ice.cpp" \
  -o "$mod/gmx-seams" $gmx_cflags -I"$root/include" \
  -L"$lib" -Wl,-rpath,"$lib" -lseams $gmx_libs -lfftw3f
env -i PATH=/usr/bin:/bin LD_LIBRARY_PATH="$lib" \
  "$mod/gmx-seams" -f "$root/engines/gromacs/ice.gro" -cutoff 3.5

echo "engines-api ok"
