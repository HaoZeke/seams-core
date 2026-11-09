"""Behaviour check for the dseams fix.

LAMMPS writes timestep 0 before any step finishes. The fix used to analyse
a frame only from end_of_step, so that frame was missing. setup now runs
the same census. A record that starts at step 1 fails this test.
"""

import live


def check(ok, message):
    if not ok:
        raise SystemExit(message)


def has_step_zero(text):
    return any(block["step"] == 0 for block in live.census_blocks(text))


def exercise_fail_arm():
    raised = False
    try:
        check(False, "fail arm")
    except SystemExit as exc:
        raised = str(exc) == "fail arm"
    check(raised, "fail arm text")


def run():
    exercise_fail_arm()
    check(not has_step_zero("step 1\nrings none\n"), "skipped record")
    folder = live.WORK / "step0"
    folder.mkdir(parents=True, exist_ok=True)
    dump = folder / "native.lammpstrj"
    record = folder / "live.txt"
    script = folder / "in.step0"
    live.lammps_input(script, 1, dump, record, live.PLUGIN, every=1, ncell=1)
    live.lammps_bin(script, folder / "log", folder)
    check(has_step_zero(record.read_text()), "timestep 0 was not analysed")
