# Engine end-to-end checks. Each runs a short ice/liquid job in the real
# engine, labels it in situ, and compares every label of every dumped frame
# with `seams chill-plus` on the same frames. A missing engine is a SKIP
# (exit 77). Engine executables come from PATH or LMP, GMX, PLUMED_KERNEL.

e2e_py = find_program('python3', required: false)
if e2e_py.found()
  e2e_src = meson.project_source_root() / 'engines'
  e2e_build = meson.project_build_root()
  e2e_env = environment()
  e2e_env.set('SEAMS_CLI', seams_exe.full_path())
  e2e_deps = [seams_exe, seams_api_lib]
  foreach name : ['lammps_plugin', 'gmx_seams', 'seams_plumed']
    if is_variable(name)
      e2e_deps += get_variable(name)
    endif
  endforeach
  test('lammps_e2e', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_mw.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  # Regression: on two ranks each rank used to label only its own atoms.
  test('lammps_e2e_np2', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_mw.py', e2e_build, '--ranks', '2'],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
endif
