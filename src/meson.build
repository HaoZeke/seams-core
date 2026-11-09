# d-SEAMS core library (yodaLib)

ydslib_sources = files(
  'absOrientation.cpp',
  'backward.cpp',
  'bond.cpp',
  'bop.cpp',
  'bulkTUM.cpp',
  'cage_affiliation.cpp',
  'cage_enum.cpp',
  'cluster.cpp',
  'domain.cpp',
  'franzblau.cpp',
  'gpu_residency.cpp',
  'gpu_batch.cpp',
  'generic.cpp',
  'ira_sofi.cpp',
  'cage_canon.cpp',
  'chemfiles_io.cpp',
  'chill_offload.cpp',
  'topo_fingerprint.cpp',
  'tum_offload.cpp',
  'mol_sys.cpp',
  'sphericart_ylm.cpp',
  'neighbours.cpp',
  'order_parameter.cpp',
  'ozaki.cpp',
  'phase.cpp',
  'pntCorrespondence.cpp',
  'density.cpp',
  'rdf.cpp',
  'rdf2d.cpp',
  'ring.cpp',
  'seams_config.cpp',
  'seams_input.cpp',
  'seams_output.cpp',
  'selection.cpp',
  'shapeMatch.cpp',
  'seams_c_api.cpp',
  'insitu.cpp',
  'insitu_imd.cpp',
  'site.cpp',
  'steinhardt_parallel.cpp',
  'topo_bulk.cpp',
  'cage_pack.cpp',
  'voronoi_qlm.cpp',
  'structure_desc.cpp',
  'topo_one_dim.cpp',
  'topo_two_dim.cpp',
)

_lib_link_args = _link_args
# nvc++ records the output path (src/libyodaLib.so) as DT_NEEDED.
# The Catch2 tests run with workdir input/, so that relative name
# does not resolve. Name the soname the same way gcc does.
if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
  _lib_link_args += ['-Wl,-soname,libyodaLib.so']
endif

ydslib = library('yodaLib',
  ydslib_sources,
  dependencies: _deps,
  cpp_args: _args,
  link_args: _lib_link_args,
  include_directories: _incdirs,
  build_rpath: _readcondb_rpath != '' ? _readcondb_rpath : meson.current_build_dir(),
  install: true)

_linkto += ydslib

yds_dep = declare_dependency(
  link_with: [ydslib],
  include_directories: _incdirs,
  compile_args: _args,
  link_args: _link_args,
  dependencies: _deps)

# Installed interface.
#
# Building the library is not the same as shipping something another
# project can build against: without headers and a pkg-config file an
# installed copy can only be linked blind, which is why consumers vendor
# this repository as a subproject instead of depending on a packaged one.
#
# The feature defines have to travel with the headers. seams_input.hpp,
# neighbours.hpp and steinhardt_device.hpp change shape under
# SEAMS_HAS_CHEMFILES, SEAMS_HAS_READCON, SEAMS_HAS_LINKCELL and
# SEAMS_HAS_OFFLOAD, so a consumer that compiles without the set this
# library was built with sees a different declaration than the one the
# library exports. yds_dep hands them to in-tree consumers through
# compile_args; extra_cflags is the same promise to an installed one.
# The warning suppressions in _args stay behind: they are this project's
# taste, not part of its interface.
install_subdir(
  'include/internal',
  install_dir: get_option('includedir') / 'seams',
  strip_directory: true)

_public_defines = []
foreach _arg : _args
  if _arg.startswith('-DSEAMS_HAS') or _arg.startswith('-D_LIBCPP_DISABLE_AVAILABILITY')
    _public_defines += _arg
  endif
endforeach

# Naming the library object here would make meson walk the dependency
# graph into Libs.private, and that graph holds absolute paths: the
# build directory of any static subproject, and whichever environment
# supplied chemfiles and highway on the build machine. Those paths do
# not survive the install, let alone reach another machine. yodaLib is
# shared and carries its own dependencies, so the consumer needs the
# link line and nothing more.
pkgconfig = import('pkgconfig')
pkgconfig.generate(
  name: 'seams-core',
  filebase: 'seams-core',
  description: 'd-SEAMS structural elucidation engine',
  version: meson.project_version(),
  subdirs: ['seams'],
  requires: ['eigen3 >= 3.4.0'],
  libraries: ['-L${libdir}', '-lyodaLib'],
  extra_cflags: _public_defines)

# The LAMMPS plugin links this library and nothing from the molecular
# dynamics code. It is built only when that code's headers are visible.
# auto skips a missing dependency. enabled requires it.
# Live analysis does not go through PLUMED, so this tree builds no
# PLUMED action. The neighbour-index fix stays in its own tree.
# The GROMACS client is seams-gromacs-imd and links no GROMACS library,
# so it is built with the core. OpenMM is the Python reporter.
lammps_adapter = false
plumed_adapter = false
if not get_option('lammps').disabled()
  _lammps_inc = []
  if get_option('lammps_include') != ''
    _lammps_inc += include_directories(get_option('lammps_include'))
  endif
  _lammps_ok = cppc.compiles(
    '#include "lammpsplugin.h"\nint main() { return 0; }\n',
    include_directories: _lammps_inc,
    name: 'LAMMPS plugin header')
  if _lammps_ok
    lammps_adapter = true
    _lammps_mpi = dependency('mpi', required: false)
    shared_module(
      'fixdseams',
      'adapters/lammps/fix_dseams.cpp',
      dependencies: [yds_dep, _lammps_mpi],
      include_directories: _lammps_inc,
      name_prefix: 'lib',
      install: true)
  elif get_option('lammps').enabled()
    error('lammps=enabled but lammpsplugin.h was not found. Pass -Dlammps_include')
  endif
endif

if get_option('plumed').enabled()
  error('Live analysis does not use PLUMED. -Dplumed=enabled is an error.')
endif
# The source stays in the tree. -Dplumed=enabled remains an error.
# The action is compiled when the headers are named for a coverage build.
if get_option('plumed_include') != ''
  shared_module(
    'dseams_plumed',
    'adapters/plumed/Dseams.cpp',
    dependencies: [yds_dep],
    include_directories: include_directories(get_option('plumed_include')),
    cpp_args: ['-DPLUMED_VERSION_MINOR=10'],
    name_prefix: 'lib',
    install: false)
endif
