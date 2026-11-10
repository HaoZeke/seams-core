# Copyright (c) d-SEAMS developers
"""Timestep 0 is analysed once, including across two run commands."""

import live


def steps_of(text):
    return [block["step"] for block in live.census_blocks(text)]


def test_a_record_without_step_zero_fails():
    assert 0 not in steps_of("step 1\nrings none\n")


def _script(path, record, plugin, body):
    record = live.Path(record).resolve()
    plugin = str(live.Path(plugin).resolve()) if plugin else plugin
    path.write_text(
        "units real\n"
        "atom_style atomic\n"
        "boundary p p p\n"
        "lattice diamond 6.4\n"
        "region box block 0 1 0 1 0 1\n"
        "create_box 1 box\n"
        "create_atoms 1 box\n"
        "mass 1 18.0154\n"
        "pair_style lj/cut 8.0\n"
        "pair_coeff * * 0.05 2.5\n"
        "neighbor 2.0 bin\n"
        "timestep 2.0\n"
        "velocity all create 30.0 12345\n"
        "fix int all nve\n"
        "plugin load %s\n"
        "fix sea all dseams record %s every 1 cutoff 3.5 depth 8\n"
        "%s"
        % (plugin, record, body)
    )


def _record(name, body):
    folder = live.WORK / name
    folder.mkdir(parents=True, exist_ok=True)
    record = folder / "live.txt"
    script = folder / "in.step0"
    _script(script, record, live.PLUGIN, body)
    live.lammps_bin(script, folder / "log", folder)
    return steps_of(record.read_text())


def test_timestep_zero_is_analysed_once():
    steps = _record("step0", "run 1\n")
    assert steps == [0, 1], steps


def test_a_later_run_does_not_repeat_the_boundary_step():
    steps = _record("step0-twice", "run 1\nrun 1\n")
    assert steps == [0, 1, 2], steps


def run():
    test_a_record_without_step_zero_fails()
    test_timestep_zero_is_analysed_once()
    test_a_later_run_does_not_repeat_the_boundary_step()
