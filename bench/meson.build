# Standalone timing driver. Not part of the Catch2 suite.
bench_link = _link_args
if host_machine.system() == 'linux'
  bench_link += ['-Wl,--enable-new-dtags']
endif

seams_bench = executable(
  'seams_bench',
  files('analyse.cpp'),
  dependencies: [yds_dep],
  include_directories: _incdirs,
  cpp_args: _args,
  link_args: bench_link,
  build_rpath: meson.project_build_root() / 'src',
  install: false,
)

# SQQ's ice rule on d-SEAMS rings, loaded by scripts/h2h_sqq_versus.py.
seams_ring_ice = shared_library(
  'seams_ring_ice',
  files('ring_ice.c'),
  link_with: seams_api_lib,
  dependencies: openmp_found ? [openmp_dep] : [],
  include_directories: include_directories('../include'),
  c_args: ['-I' + meson.project_build_root() / 'src'],
  build_rpath: meson.project_build_root() / 'src',
  install: false,
)
