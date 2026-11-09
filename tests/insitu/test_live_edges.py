# Copyright (c) d-SEAMS developers
"""Edges of the live driver that a passing engine run does not take."""

import os
import stat
from pathlib import Path

import pytest

import live


def test_an_empty_census_has_no_blocks():
    assert live.census_blocks("") == []


def test_a_short_time_line_is_skipped():
    assert live.gmx_wall("Time: only\nTime: 1.5 2.5\n") == 2.5


def test_duplicate_live_steps_are_rejected(tmp_path):
    record = tmp_path / "live.txt"
    record.write_text("step 1\nrings a\nstep 1\nrings a\n")
    native = tmp_path / "native.lammpstrj"
    native.write_text("ITEM: TIMESTEP\n0\n")
    with pytest.raises(RuntimeError, match="duplicate live steps"):
        live.match_native_cli(record, native)


def test_a_second_dump_frame_is_kept(tmp_path):
    dump = tmp_path / "two.lammpstrj"
    dump.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
        "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
        "ITEM: ATOMS id type x y z\n1 1 0.1 0 0\n"
    )
    frames = live.parse_dump(dump)
    assert [frame[0] for frame in frames] == [0, 1]


def test_lammps_mpi_reports_absent_tools(tmp_path, monkeypatch):
    monkeypatch.setenv("DSEAMS_LMP_MPI_PREFIX", str(tmp_path))
    monkeypatch.setenv("DSEAMS_LMP_MPI_PLUGIN", str(tmp_path / "missing.so"))
    monkeypatch.setattr(live, "WORK", tmp_path / "work")
    same, note = live.lammps_mpi_identity()
    assert same is False
    assert note == "not run: mpirun or the plugin is absent"


def test_lammps_ranks_that_disagree(tmp_path, monkeypatch):
    root = tmp_path / "mpi"
    (root / "bin").mkdir(parents=True)
    launcher = root / "bin" / "mpirun"
    launcher.write_text(
        "#!/bin/sh\n"
        "np=1\n"
        "prev=\n"
        "for arg in \"$@\"; do\n"
        "  if [ \"$prev\" = \"-np\" ]; then np=$arg; fi\n"
        "  prev=$arg\n"
        "done\n"
        "printf 'step %s\\n' \"$np\" > live.txt\n"
        "exit 0\n"
    )
    launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC)
    (root / "bin" / "lmp_mpi").write_text("#!/bin/sh\nexit 0\n")
    plugin = tmp_path / "plugin.so"
    plugin.write_bytes(b"plugin")
    monkeypatch.setenv("DSEAMS_LMP_MPI_PREFIX", str(root))
    monkeypatch.setenv("DSEAMS_LMP_MPI_PLUGIN", str(plugin))
    monkeypatch.setattr(live, "WORK", tmp_path / "work")
    same, note = live.lammps_mpi_identity()
    assert same is False
    assert note == "The census does not match the one-rank run."


def test_cpuinfo_without_a_model_name(monkeypatch):
    real = Path.read_text

    def read_text(self, *args, **kwargs):
        if str(self) == "/proc/cpuinfo":
            return "processor : 0\n"
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read_text)
    assert read_text(Path(__file__))
    assert live.host_block()["cpu"] == "unknown"


def test_a_scale_lattice_with_the_wrong_count(monkeypatch, tmp_path):
    monkeypatch.setattr(live, "WORK", tmp_path)
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [])
    monkeypatch.setattr(live, "scale_plan", lambda: [(1, 8, 1)])
    with pytest.raises(RuntimeError, match="lattice 8 has 0 atoms"):
        live.scale_gromacs()


def test_an_existing_scale_folder_is_replaced(monkeypatch, tmp_path):
    monkeypatch.setattr(live, "WORK", tmp_path)
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [{}] * 8)
    monkeypatch.setattr(live, "scale_plan", lambda: [(1, 8, 1)])
    seen = {}

    def fake_case(folder, rows, steps, every, repeats):
        seen["path"] = folder
        return {"kept": True}

    monkeypatch.setattr(live, "gromacs_scale_case", fake_case)
    folder = tmp_path / "gscale" / "8" / "1"
    folder.mkdir(parents=True)
    (folder / "old.txt").write_text("stale")
    result = live.scale_gromacs()
    assert result["8"]["1"]["kept"] is True
    assert not (folder / "old.txt").exists()


def test_ipi_that_exits_before_the_client(monkeypatch, tmp_path):
    root = tmp_path / "md"
    (root / "bin").mkdir(parents=True)
    ipi = root / "bin" / "i-pi"
    ipi.write_text("#!/bin/sh\nexit 1\n")
    ipi.chmod(ipi.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "WORK", tmp_path / "work")
    with pytest.raises(RuntimeError):
        live.run_ipi()
    with pytest.raises(RuntimeError):
        live.ipi_elapsed(False)


def test_an_imd_client_that_fails(monkeypatch, tmp_path):
    root = tmp_path / "md"
    (root / "bin").mkdir(parents=True)
    gmx = root / "bin" / "gmx"
    gmx.write_text(
        "#!/bin/sh\n"
        "for arg in \"$@\"; do\n"
        "  if [ \"$arg\" = \"mdrun\" ]; then exit 1; fi\n"
        "done\n"
        "exit 0\n"
    )
    gmx.chmod(gmx.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "IMD", "/bin/false")
    monkeypatch.setenv("SEAMS_GRO_IMD_EXE", "/bin/false")
    folder = tmp_path / "gro"
    folder.mkdir()
    (folder / "solv.gro").write_text("title\n0\n   1.0   1.0   1.0\n")
    (folder / "md.mdp").write_text("nsteps = NSTEPS\n")
    (folder / "topol.top").write_text("[ system ]\n")
    with pytest.raises(RuntimeError):
        live.gromacs_once(folder, 1, 9)


def test_a_bare_live_invocation_exits():
    with pytest.raises(SystemExit, match="pass --ladder-worker"):
        live.module_entry(["live.py"])


def test_main_stops_when_the_ranks_differ(monkeypatch, tmp_path):
    work = tmp_path / "work"
    monkeypatch.setattr(live, "WORK", work)
    monkeypatch.setattr(live, "LIB", str(tmp_path / "lib.so"))
    monkeypatch.setattr(live, "PLUGIN", tmp_path / "plugin.so")
    root = tmp_path / "md"
    (root / "bin").mkdir(parents=True)
    for name in ("lmp", "gmx"):
        path = root / "bin" / name
        path.write_bytes(b"bin")
    (tmp_path / "lib.so").write_bytes(b"lib")
    (tmp_path / "plugin.so").write_bytes(b"plugin")
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "run_lammps", lambda: {"frames": 2, "steps": 4})
    monkeypatch.setattr(live, "time_lammps", lambda *args, **kwargs: [0.1])
    monkeypatch.setattr(live, "run_gromacs", lambda: {"frames": 2, "steps": 4})
    monkeypatch.setattr(live, "time_gromacs", lambda *args, **kwargs: ([0.1], [0.2]))
    monkeypatch.setattr(live, "openmm_identity", lambda: 2)
    monkeypatch.setattr(live, "openmm_time", lambda steps: (0.1, 0.2))
    monkeypatch.setattr(live, "ase_identity", lambda: 2)
    monkeypatch.setattr(live, "ase_time", lambda *args, **kwargs: 0.1)
    monkeypatch.setattr(live, "run_ipi", lambda: 2)
    monkeypatch.setattr(live, "ipi_elapsed", lambda analyse: 0.1)
    monkeypatch.setattr(live, "scale_lammps", lambda: {"8": {"1": {"live": [1.0], "plain": [0.5]}}})
    monkeypatch.setattr(live, "gromacs_trr_identity", lambda: 2)
    monkeypatch.setattr(live, "scale_gromacs", lambda: {})
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [])
    monkeypatch.setattr(live, "frame_clocks", lambda dump, per_frame: {"cli_excluded": 0.1})
    monkeypatch.setattr(live, "lammps_bin", lambda *args, **kwargs: None)
    monkeypatch.setattr(live, "lammps_plain", lambda *args, **kwargs: None)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 1)
    monkeypatch.setattr(live, "thread_ladder", lambda dump, natoms: {"threads": 1})
    monkeypatch.setattr(
        live,
        "lammps_mpi_identity",
        lambda: (False, "The census does not match the one-rank run."),
    )
    monkeypatch.setattr(
        live,
        "gromacs_mpi_identity",
        lambda: (True, "The native module on one, two and four ranks matches the trajectory."),
    )
    with pytest.raises(RuntimeError, match="does not match the one-rank run"):
        live.main()
    saved = work / "samples.json"
    assert saved.is_file()
    assert "does not match" in saved.read_text()


def test_a_clock_that_prints_no_sample(monkeypatch, tmp_path):
    dump = tmp_path / "one.lammpstrj"
    dump.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
        "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
    )
    clock = tmp_path / "clock"
    clock.write_text("#!/bin/sh\nexit 0\n")
    clock.chmod(clock.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("DSEAMS_FRAME_CLOCK", str(clock))
    monkeypatch.setattr(live, "cli_text", lambda *args, **kwargs: "")
    with pytest.raises(RuntimeError):
        live.frame_clocks(dump, 0.1)


def test_the_fallback_clock_is_the_release_binary(monkeypatch, tmp_path):
    dump = tmp_path / "one.lammpstrj"
    dump.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
        "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
    )
    monkeypatch.delenv("DSEAMS_FRAME_CLOCK", raising=False)
    monkeypatch.setattr(live, "cli_text", lambda *args, **kwargs: "")
    seen = {}

    def fake_run(cmd, cwd, env=None):
        seen["clock"] = cmd[0]

        class Proc:
            stdout = "analysis_s 0.25\n"
            stderr = ""
            returncode = 0

        return Proc()

    monkeypatch.setattr(live, "run", fake_run)
    result = live.frame_clocks(dump, 0.2)
    assert seen["clock"].endswith("frame_clock")
    assert result["cli_excluded"] == 0.25
