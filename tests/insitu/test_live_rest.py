# Copyright (c) d-SEAMS developers
"""Branch coverage for the live in-situ driver, with synthetic frames."""

import os
import stat
import struct
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import live  # noqa: E402

SEAMS = REPO / "bbdir-cov" / "seams"


def _exe(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _dump(path, steps=((0, [(1, 0.0, 0.0, 0.0), (2, 2.8, 0.0, 0.0)]),)):
    lines = []
    for step, rows in steps:
        lines.append("ITEM: TIMESTEP")
        lines.append(str(step))
        lines.append("ITEM: NUMBER OF ATOMS")
        lines.append(str(len(rows)))
        lines.append("ITEM: BOX BOUNDS pp pp pp")
        lines.extend(["0 10", "0 10", "0 10"])
        lines.append("ITEM: ATOMS id type x y z")
        for ident, x, y, z in rows:
            lines.append("%d 1 %.17g %.17g %.17g" % (ident, x, y, z))
    path.write_text("\n".join(lines) + "\n")


def _lmp_stub():
    return textwrap.dedent(
        """\
        #!/bin/sh
        in=""
        prev=""
        for arg in "$@"; do
          if [ "$prev" = "-in" ]; then in=$arg; fi
          prev=$arg
        done
        dump=$(awk '/^[[:space:]]*dump /{print $6; exit}' "$in")
        if [ -n "$dump" ]; then
          cat > "$dump" << 'END'
        ITEM: TIMESTEP
        0
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        0 10
        0 10
        0 10
        ITEM: ATOMS id type x y z
        2 1 1 0 0
        1 1 0 0 0
        END
        fi
        echo "Loop time of 2.0 on 1 procs for 4 steps"
        log=""
        prev=""
        for arg in "$@"; do
          if [ "$prev" = "-log" ]; then log=$arg; fi
          prev=$arg
        done
        if [ -n "$log" ]; then
          echo "Loop time of 2.0 on 1 procs for 4 steps" > "$log"
        fi
        exit 0
        """
    )


def _gmx_body(time_line="Time: 0.5 2.0"):
    return textwrap.dedent(
        """\
        #!/bin/sh
        cmd=""
        for arg in "$@"; do
          case "$arg" in
            pdb2gmx|editconf|solvate|grompp|mdrun) cmd=$arg ;;
          esac
        done
        if [ "$cmd" = "solvate" ]; then
          printf 'title\\n1\\n    1SOL     OW    1   0.1   0.1   0.1\\n   2.0   2.0   2.0\\n' > solv.gro
        fi
        if [ "$cmd" = "mdrun" ]; then
          echo "%s"
          deffnm=""
          prev=""
          for arg in "$@"; do
            if [ "$prev" = "-deffnm" ]; then deffnm=$arg; fi
            prev=$arg
          done
          if [ -n "$deffnm" ]; then
            printf 'step 1\\n' > "${deffnm}.txt" 2>/dev/null || true
            : > "${deffnm}.trr"
          fi
        fi
        exit 0
        """
        % time_line
    )


class _Proc:
    def __init__(self, code=0, out="ok\n"):
        self.returncode = code
        self._out = out

    def poll(self):
        return None

    def communicate(self, timeout=None):
        return self._out, ""


@pytest.fixture
def work(tmp_path, monkeypatch):
    monkeypatch.setattr(live, "WORK", tmp_path / "work")
    monkeypatch.setattr(live.time, "sleep", lambda _seconds: None)
    return tmp_path


def test_census_keeps_a_note_and_counts_a_trr_and_a_short_dcd(tmp_path, monkeypatch):
    blocks = live.census_blocks("step 3\nrings 3:0\nnote\nchill water:1\ncage hc:0\n")
    assert blocks == [
        {"step": 3, "rings": "rings 3:0", "chill": {"water": 1}, "cage": {"hc": 0}}
    ]
    trr = tmp_path / "one.trr"
    _pack_trr(trr, 1, 4, [0.1, 0.2, 0.3])
    assert live.native_frame_count(trr) == 1
    short = tmp_path / "short.dcd"
    short.write_bytes(b"1234")
    with pytest.raises(RuntimeError, match="short dcd"):
        live.native_frame_count(short)
    counted = tmp_path / "two.dcd"
    counted.write_bytes(struct.pack("<i4siii", 84, b"CORD", 6, 1, 1) + b"\0" * 16)
    assert live.native_frame_count(counted) == 6
    dump = tmp_path / "one.lammpstrj"
    _dump(dump)
    echo = tmp_path / "echo-seams"
    _exe(echo, "#!/bin/sh\nprintf '%s\\n' \"$*\"\n")
    monkeypatch.setattr(live, "SEAMS", str(echo))
    text = live.cli_text("rings", dump, extra=["--graph", "cutoff"], type_filter=2)
    assert "--graph" in text and "cutoff" in text and "--type 2" in text


def test_dump_steps_ignore_a_timestep_header_with_no_value(tmp_path):
    path = tmp_path / "odd.lammpstrj"
    path.write_text("ITEM: TIMESTEP\n4\nITEM: TIMESTEP\n")
    assert live.dump_steps(path) == [4]


def _pack_trr(path, natoms, step, xyz):
    box = [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0]
    payload = struct.pack(">9f", *box) + struct.pack(">%df" % (natoms * 3), *xyz)
    header = bytearray(84)
    struct.pack_into(">i", header, 0, 1993)
    struct.pack_into(">i", header, 32, 36)
    struct.pack_into(">i", header, 52, natoms * 3 * 4)
    struct.pack_into(">ii", header, 64, natoms, step)
    path.write_bytes(bytes(header) + payload)


def test_a_truncated_trr_is_rejected(tmp_path):
    path = tmp_path / "over.trr"
    _pack_trr(path, 1, 2, [0.0, 0.0, 0.0])
    path.write_bytes(path.read_bytes() + b"\x00")
    with pytest.raises(RuntimeError, match="trr length"):
        live.trr_frames(path)
    clean = tmp_path / "ok.trr"
    _pack_trr(clean, 1, 2, [0.1, 0.2, 0.3])
    frames = live.trr_frames(clean)
    assert frames[0][0] == 2 and frames[0][1] == 1


def test_parse_dump_of_an_empty_file_is_empty(tmp_path):
    path = tmp_path / "empty.lammpstrj"
    path.write_text("comment only\n")
    assert live.parse_dump(path) == []


def test_write_lammps_and_read_multi_xyz_round_trip(tmp_path):
    dest = tmp_path / "out.lammpstrj"
    live.write_lammps(dest, [], [0, 0, 0, 1, 1, 1], [])
    assert dest.read_text() == "\n"
    live.write_lammps(dest, [(3, [[1.0, 2.0, 3.0]])], [0, 0, 0, 5, 5, 5], [1])
    assert "ITEM: TIMESTEP" in dest.read_text() and "\n3\n" in dest.read_text()
    xyz = tmp_path / "mol.xyz"
    xyz.write_text("\n1\ncomment\nO 0.5 0.25 0.125\n")
    frames = live.read_multi_xyz(xyz)
    assert frames == [(0, [[0.5, 0.25, 0.125]])]


def test_a_gro_box_needs_three_lengths(tmp_path):
    path = tmp_path / "bad.gro"
    path.write_text("title\n0\n1.0 2.0\n")
    with pytest.raises(RuntimeError, match="gro box"):
        live.gromacs_box(path)


def test_same_record_counts_steps_and_rejects_a_difference(tmp_path):
    left = tmp_path / "a.txt"
    right = tmp_path / "b.txt"
    left.write_text("step 1\nrings 3:0\n")
    right.write_text(left.read_text())
    assert live.same(left, right) == 1
    right.write_text("step 2\n")
    with pytest.raises(RuntimeError, match="differ"):
        live.same(left, right)
    empty = tmp_path / "e.txt"
    empty.write_text("no census\n")
    with pytest.raises(RuntimeError, match="empty census"):
        live.same(empty, empty)


def test_run_rejects_a_failing_command(tmp_path):
    with pytest.raises(RuntimeError, match="command failed"):
        live.run(["/bin/false"], tmp_path)
    quiet = live.run(["/bin/false"], tmp_path, check=False)
    assert quiet.returncode != 0


def test_a_preamble_before_the_first_step_is_ignored():
    assert live.census_blocks("preamble\nstep 2\nrings 3:1\n")[0]["step"] == 2


def test_dcd_schedule_rejects_a_short_or_illegal_header(tmp_path):
    short = tmp_path / "short.dcd"
    short.write_bytes(b"1234")
    with pytest.raises(RuntimeError, match="short dcd"):
        live.dcd_schedule(short)
    bad = tmp_path / "bad.dcd"
    bad.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 2, 1, 0, 1) + b"\0" * 64)
    with pytest.raises(RuntimeError, match="bad dcd schedule"):
        live.dcd_schedule(bad)


def test_require_equal_and_native_frame_counts():
    live.require_equal(3, 3, "same")
    with pytest.raises(RuntimeError, match="mismatch"):
        live.require_equal(1, 2, "mismatch")
    live.require_native_frames(4, 4, 4)
    with pytest.raises(RuntimeError, match="gromacs native 1"):
        live.require_native_frames(1, 2, 2)


def test_gmx_wall_skips_a_note_and_rejects_silence():
    assert live.gmx_wall("note\nTime: 1.0 3.0\n") == 3.0
    with pytest.raises(RuntimeError, match="no Time line"):
        live.gmx_wall("silent\n")


def test_a_bad_trr_magic_and_a_short_ladder(tmp_path):
    path = tmp_path / "bad.trr"
    path.write_bytes(b"\0" * 84)
    with pytest.raises(RuntimeError, match="unexpected trr magic"):
        live.trr_frames(path)
    with pytest.raises(RuntimeError):
        live.ladder_values("0.1\n0.2\n", "")


def test_loop_time_skips_a_preamble(tmp_path):
    assert live.loop_time("warmup\nLoop time of 2.5 on 1 procs for 4 steps\n") == 2.5
    with pytest.raises(RuntimeError, match="loop time"):
        live.loop_time("no timing\n")


def _record_from_seams(dump):
    saved = live.SEAMS
    live.SEAMS = str(SEAMS)

    def lines(command, extra=None):
        return [
            ln
            for ln in live.cli_text(command, dump, extra=extra).splitlines()
            if ln.startswith("rings ") or ln.startswith("nop ")
        ]

    try:
        rings = lines("rings")[0]
        chill = lines("chill-plus")[0]
        cage = lines("cages", ["--graph", "cutoff"])[0]
    finally:
        live.SEAMS = saved
    chill_words = chill.split()
    chill_bits = [
        "%s:%s" % (chill_words[i], chill_words[i + 1]) for i in range(2, len(chill_words), 2)
    ]
    cage_words = cage.split()
    mapping = {
        cage_words[i]: cage_words[i + 1] for i in range(4, len(cage_words), 2)
    }
    cage_bits = "hc:%s ddc:%s both:%s water:%s" % (
        mapping.get("hexagonal", "0"),
        mapping.get("cubic", "0"),
        mapping.get("both", "0"),
        mapping.get("water", "0"),
    )
    step = live.dump_steps(dump)[0]
    return "step %d\n%s\nchill %s\ncage %s\n" % (step, rings, " ".join(chill_bits), cage_bits)


def test_match_lammps_cli_agrees_with_the_seams_binary(tmp_path, monkeypatch):
    dump = tmp_path / "pair.lammpstrj"
    _dump(dump)
    record = tmp_path / "live.txt"
    monkeypatch.setattr(live, "SEAMS", str(SEAMS))
    record.write_text(_record_from_seams(dump))
    assert live.match_lammps_cli(record, dump) == 1
    dump.write_text(dump.read_text() + dump.read_text().replace("\n0\n", "\n1\n", 1))
    assert live.dump_steps(dump) == [0, 1]
    record.write_text(_record_from_seams_all(dump))
    assert live.match_lammps_cli(record, dump) == 2


def _record_from_seams_all(dump):
    saved = live.SEAMS
    live.SEAMS = str(SEAMS)

    def kept(command, prefix, extra=None):
        return [
            ln
            for ln in live.cli_text(command, dump, extra=extra).splitlines()
            if ln.startswith(prefix)
        ]

    try:
        rings = kept("rings", "rings ")
        chills = kept("chill-plus", "nop ")
        cages = kept("cages", "nop ", ["--graph", "cutoff"])
    finally:
        live.SEAMS = saved
    steps = live.dump_steps(dump)
    blocks = []
    for step, ring, chill, cage in zip(steps, rings, chills, cages):
        chill_words = chill.split()
        chill_bits = [
            "%s:%s" % (chill_words[i], chill_words[i + 1])
            for i in range(2, len(chill_words), 2)
        ]
        cage_words = cage.split()
        mapping = {cage_words[i]: cage_words[i + 1] for i in range(4, len(cage_words), 2)}
        cage_bits = "hc:%s ddc:%s both:%s water:%s" % (
            mapping.get("hexagonal", "0"),
            mapping.get("cubic", "0"),
            mapping.get("both", "0"),
            mapping.get("water", "0"),
        )
        blocks.append("step %d\n%s\nchill %s\ncage %s\n" % (step, ring, " ".join(chill_bits), cage_bits))
    return "".join(blocks)


def test_match_lammps_cli_rejects_each_disagreement(tmp_path, monkeypatch):
    dump = tmp_path / "one.lammpstrj"
    _dump(dump)
    stub = tmp_path / "seams"
    _exe(
        stub,
        "#!/bin/sh\n"
        "case \"$1\" in\n"
        "  rings) echo 'rings 3:0 4:0 5:0 6:0 7:0 8:0' ;;\n"
        "  chill-plus) echo 'nop 2 water 2' ;;\n"
        "  cages) echo 'nop 2 graph cutoff hexagonal 0 cubic 0 both 0 water 2' ;;\n"
        "esac\n",
    )
    monkeypatch.setattr(live, "SEAMS", str(stub))
    good = (
        "step 0\n"
        "rings 3:0 4:0 5:0 6:0 7:0 8:0\n"
        "chill water:2\n"
        "cage hc:0 ddc:0 both:0 water:2\n"
    )
    record = tmp_path / "live.txt"
    record.write_text("step 0\nrings 3:0\nstep 0\nrings 3:0\n")
    with pytest.raises(RuntimeError, match="duplicate live steps"):
        live.match_lammps_cli(record, dump)
    _dump(dump, steps=((1, [(1, 0, 0, 0)]), (1, [(1, 0, 0, 0)])))
    record.write_text("step 1\nrings 3:0\nstep 2\nrings 3:0\n")
    with pytest.raises(RuntimeError, match="duplicate dump steps"):
        live.match_lammps_cli(record, dump)
    _dump(dump)
    record.write_text("step 9\n" + good.split("\n", 1)[1])
    with pytest.raises(RuntimeError, match="step list differs"):
        live.match_lammps_cli(record, dump)
    _exe(stub, "#!/bin/sh\nexit 0\n")
    record.write_text(good)
    with pytest.raises(RuntimeError, match="frame counts differ"):
        live.match_lammps_cli(record, dump)
    _exe(
        stub,
        "#!/bin/sh\n"
        "case \"$1\" in\n"
        "  rings) echo 'rings 3:0 4:0 5:0 6:0 7:0 8:0' ;;\n"
        "  chill-plus) echo 'nop 2 water 2' ;;\n"
        "  cages) echo 'nop 2 graph cutoff hexagonal 0 cubic 0 both 0 water 2' ;;\n"
        "esac\n",
    )
    record.write_text(good.replace("rings 3:0", "rings 3:9", 1))
    with pytest.raises(RuntimeError, match="rings differ"):
        live.match_lammps_cli(record, dump)
    record.write_text(good.replace("water:2", "water:9", 1))
    with pytest.raises(RuntimeError, match="chill water"):
        live.match_lammps_cli(record, dump)
    record.write_text(good.replace("hc:0", "hc:4", 1))
    with pytest.raises(RuntimeError, match="cages differ"):
        live.match_lammps_cli(record, dump)
    record.write_text(good.replace("both:0", "both:3", 1))
    with pytest.raises(RuntimeError, match="cage water differ"):
        live.match_lammps_cli(record, dump)


def test_match_native_cli_writes_the_replay_and_checks_a_dcd(tmp_path, monkeypatch):
    dump = tmp_path / "one.lammpstrj"
    _dump(dump)
    stub = tmp_path / "seams"
    _exe(
        stub,
        "#!/bin/sh\n"
        "case \"$1\" in\n"
        "  rings) echo 'rings 3:0 4:0 5:0 6:0 7:0 8:0' ;;\n"
        "  chill-plus) echo 'nop 2 water 2' ;;\n"
        "  cages) echo 'nop 2 graph cutoff hexagonal 0 cubic 0 both 0 water 2' ;;\n"
        "esac\n",
    )
    monkeypatch.setattr(live, "SEAMS", str(stub))
    record = tmp_path / "live.txt"
    record.write_text(
        "step 0\nrings 3:0 4:0 5:0 6:0 7:0 8:0\nchill water:2\n"
        "cage hc:0 ddc:0 both:0 water:2\n"
    )
    assert live.match_native_cli(record, dump) == 1
    assert "rings 3:0" in (tmp_path / "one.lammpstrj.rings.txt").read_text()
    dcd = tmp_path / "native.dcd"
    dcd.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 1, 0, 1, 0) + b"\0" * 64)
    assert live.match_native_cli(record, dcd) == 1
    shifted = tmp_path / "shift.dcd"
    shifted.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 1, 5, 1, 5) + b"\0" * 64)
    with pytest.raises(RuntimeError, match="DCD step list differs"):
        live.match_native_cli(record, shifted)
    record.write_text(
        "step 0\nrings 3:9 4:0 5:0 6:0 7:0 8:0\nchill water:2\n"
        "cage hc:0 ddc:0 both:0 water:2\n"
    )
    with pytest.raises(RuntimeError, match="native rings differ"):
        live.match_native_cli(record, dump)
    record.write_text(
        "step 0\nrings 3:0 4:0 5:0 6:0 7:0 8:0\nchill water:9\n"
        "cage hc:0 ddc:0 both:0 water:2\n"
    )
    with pytest.raises(RuntimeError, match="native chill"):
        live.match_native_cli(record, dump)
    record.write_text(
        "step 0\nrings 3:0 4:0 5:0 6:0 7:0 8:0\nchill water:2\n"
        "cage hc:4 ddc:0 both:0 water:2\n"
    )
    with pytest.raises(RuntimeError, match="native cages differ"):
        live.match_native_cli(record, dump)
    record.write_text(
        "step 0\nrings 3:0 4:0 5:0 6:0 7:0 8:0\nchill water:2\n"
        "cage hc:0 ddc:0 both:3 water:2\n"
    )
    with pytest.raises(RuntimeError, match="native cage water differ"):
        live.match_native_cli(record, dump)
    _exe(stub, "#!/bin/sh\nexit 0\n")
    record.write_text("step 0\nrings 3:0\n")
    with pytest.raises(RuntimeError, match="native frame counts differ"):
        live.match_native_cli(record, dump)


def test_openmm_dcd_round_trip_and_an_empty_buffer(tmp_path):
    empty = tmp_path / "empty.dcd"
    empty.write_bytes(b"\0" * 276)
    assert live.read_openmm_dcd(empty, 2) == []
    path = tmp_path / "one.dcd"
    natoms = 2
    blob = bytearray(276)
    struct.pack_into("<iii", blob, 8, 1, 1, 1)
    blob += b"\0" * 56
    for axis in range(3):
        blob += struct.pack("<i", 4 * natoms)
        blob += struct.pack("<ff", float(axis), float(axis) + 0.5)
        blob += b"\0\0\0\0"
    path.write_bytes(bytes(blob))
    frames = live.read_openmm_dcd(path, natoms)
    assert len(frames) == 1
    assert frames[0][0] == 1
    assert frames[0][1][0][0] == pytest.approx(0.0)
    assert frames[0][1][1][1] == pytest.approx(1.5)


def _install_openmm(monkeypatch):
    import types

    import numpy as np

    class _Unit:
        def __mul__(self, other):
            return self

        def __rmul__(self, other):
            return self

        def __rtruediv__(self, other):
            return self

    unit = types.ModuleType("openmm.unit")
    unit.kelvin = _Unit()
    unit.nanometer = _Unit()
    unit.picosecond = _Unit()
    unit.angstrom = _Unit()
    app = types.ModuleType("openmm.app")

    class Element:
        def __init__(self, symbol):
            self.symbol = symbol

    class Atom:
        def __init__(self, symbol):
            self.element = Element(symbol)

    class Topology:
        def __init__(self, symbols):
            self._atoms = [Atom(s) for s in symbols]

        def atoms(self):
            return iter(self._atoms)

    class PDB:
        def __init__(self, path):
            symbols = []
            for line in Path(path).read_text().splitlines():
                if not line.startswith("ATOM"):
                    continue
                symbols.append("O" if " O " in line else "H")
            self.topology = Topology(symbols)
            self.positions = np.zeros((len(symbols), 3))

    class Modeller:
        def __init__(self, topology, positions):
            self.topology = topology
            self.positions = positions

    class System:
        pass

    class ForceField:
        def __init__(self, name):
            self.name = name

        def createSystem(self, topology, **kwargs):
            return System()

    class Context:
        def setPositions(self, positions):
            self.positions = positions

        def getState(self, getPositions=False, enforcePeriodicBox=False):
            return _State(np.zeros((3, 3)))

    class _State:
        def __init__(self, positions):
            self._positions = positions

        def getPositions(self, asNumpy=False):
            return self

        def getPeriodicBoxVectors(self, asNumpy=False):
            return self

        def value_in_unit(self, unit_name):
            return self._positions

    class Simulation:
        def __init__(self, topology, system, integrator, platform):
            self.topology = topology
            self.reporters = []
            self.context = Context()
            self.currentStep = 0

        def step(self, count):
            for _ in range(count):
                self.currentStep += 1
                state = _State(np.zeros((len(list(self.topology.atoms())), 3)))
                for reporter in list(self.reporters):
                    if hasattr(reporter, "describeNextReport"):
                        if reporter.describeNextReport(self)["steps"] != 1:
                            continue
                    if hasattr(reporter, "report"):
                        reporter.report(self, state)

    class DCDReporter:
        def __init__(self, path, interval):
            self.path = path
            self.interval = interval
            self.steps = []

        def describeNextReport(self, simulation):
            step = int(simulation.currentStep)
            remain = self.interval - (step % self.interval)
            return {"steps": remain}

        def report(self, simulation, state):
            self.steps.append(int(simulation.currentStep))
            natoms = state.getPositions(asNumpy=True).value_in_unit("A").shape[0]
            blob = bytearray(276)
            n = len(self.steps)
            stride = self.steps[1] - self.steps[0] if n > 1 else 1
            struct.pack_into("<iii", blob, 8, n, self.steps[0], stride)
            for _frame in range(n):
                blob += b"\0" * 56
                for _axis in range(3):
                    blob += struct.pack("<i", 4 * natoms)
                    blob += struct.pack("<%df" % natoms, *([0.0] * natoms))
                    blob += b"\0\0\0\0"
            Path(self.path).write_bytes(bytes(blob))

    app.ForceField = ForceField
    app.Modeller = Modeller
    app.PDBFile = PDB
    app.Simulation = Simulation
    app.CutoffPeriodic = "cutoff"
    app.HBonds = "hbonds"
    app.DCDReporter = DCDReporter
    package = types.ModuleType("openmm")
    package.unit = unit
    package.LangevinIntegrator = lambda *args: object()
    package.Platform = types.SimpleNamespace(getPlatformByName=lambda name: name)
    monkeypatch.setitem(sys.modules, "openmm", package)
    monkeypatch.setitem(sys.modules, "openmm.unit", unit)
    monkeypatch.setitem(sys.modules, "openmm.app", app)


class _Engine:
    def __init__(self, lib, path):
        self.path = Path(path)
        self.path.write_text("")
        self.observed = 0

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def observe(self, step, pos, box, dump, types=None):
        self.observed += 1
        self.path.write_text(self.path.read_text() + "step %s\n" % step)

    def box_from_rows(self, rows):
        return [0, 0, 0, 6.4, 6.4, 6.4, 0, 0, 0]

    def close(self):
        return None


def test_openmm_identity_replays_the_saved_steps(work, monkeypatch):
    _install_openmm(monkeypatch)
    monkeypatch.setattr(live, "load", lambda path=None: object())
    monkeypatch.setattr(live, "Engine", _Engine)

    def matched(record, native, type_filter=1):
        path = Path(native)
        if path.suffix == ".dcd":
            return len(live.dcd_schedule(path))
        return len(live.dump_steps(path))

    monkeypatch.setattr(live, "match_lammps_cli", matched)
    monkeypatch.setattr(live, "match_native_cli", matched)
    frames = live.openmm_identity()
    assert frames == 8
    log = (live.WORK / "openmm" / "dcd-equal.log").read_text()
    assert log == "frames 8 native 8\n"
    folder, sim, types = live.run_openmm(4)
    assert types.count(1) == 27 and types.count(2) == 54
    plain, timed = live.openmm_time(4)
    assert plain >= 0 and timed >= 0
    assert (folder / "time.txt").is_file()


def test_diamond_positions_are_read_once(work, monkeypatch):
    root = work / "md"
    _exe(root / "bin" / "lmp", _lmp_stub())
    monkeypatch.setattr(live, "PREFIX", root)
    rows = live.diamond_positions(1)
    assert [int(row["id"]) for row in rows] == [1, 2]
    _exe(root / "bin" / "lmp", "#!/bin/sh\nexit 1\n")
    again = live.diamond_positions(1)
    assert [row["x"] for row in again] == [row["x"] for row in rows]


def test_lammps_timing_and_the_rejected_fix(work, monkeypatch):
    root = work / "md"
    _exe(root / "bin" / "lmp", _lmp_stub())
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "PLUGIN", work / "plugin.so")
    (work / "plugin.so").write_bytes(b"plugin")
    monkeypatch.setattr(live, "SEAMS", str(SEAMS))
    plain = live.time_lammps("plain", 4, 2, live.PLUGIN)
    live_samples = live.time_lammps("live", 4, 1, live.PLUGIN)
    assert plain == [0.5, 0.5]
    assert live_samples == [0.5]
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 6)
    result = live.run_lammps()
    assert result == {"frames": 6, "steps": 40}
    assert (live.WORK / "lammps" / "native.lammpstrj").is_file()
    assert "fix sea all dseams dump" in (live.WORK / "lammps" / "in.bad").read_text()


def test_gromacs_plain_and_live_samples(work, monkeypatch):
    root = work / "md"
    _exe(root / "bin" / "gmx", _gmx_body())
    monkeypatch.setattr(live, "PREFIX", root)
    imd = work / "imd"
    _exe(imd, "#!/bin/sh\nprintf 'step 1\\nstep 2\\n' > \"$4\"\nexit 0\n")
    monkeypatch.setattr(live, "IMD", str(imd))
    prepared = live.WORK / "gromacs"
    # prepare_gromacs creates the folder. Run the public entry.
    result = live.run_gromacs()
    assert result == {"frames": 2, "steps": 20}
    plain, timed = live.time_gromacs(2, 4)
    assert len(plain) == 2 and len(timed) == 2
    assert min(plain) >= 0 and min(timed) >= 0
    bad = work / "imd-bad"
    _exe(bad, "#!/bin/sh\nexit 1\n")
    monkeypatch.setattr(live, "IMD", str(bad))
    with pytest.raises(RuntimeError):
        live.gromacs_once(live.WORK / "gromacs", 4, 9)


def test_scale_plan_and_the_lammps_ladder(work, monkeypatch):
    monkeypatch.delenv("DSEAMS_SCALE_TINY", raising=False)
    full = live.scale_plan()
    assert [row[1] for row in full] == [1000, 8000, 64000]
    monkeypatch.setenv("DSEAMS_SCALE_TINY", "1")
    assert live.scale_plan() == [(5, 1000, 2)]
    root = work / "md"
    _exe(root / "bin" / "lmp", _lmp_stub())
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "PLUGIN", work / "plugin.so")
    (work / "plugin.so").write_bytes(b"p")
    monkeypatch.setattr(live, "SEAMS", str(SEAMS))
    result = live.scale_lammps()
    block = result["1000"]["1"]
    assert block["steps"] == 2
    assert block["plain"] == [1.0] * 6
    assert block["live"] == [1.0] * 6
    assert block["warmup_plain"] == 1.0
    assert len(block["post_wall"]) == 6
    assert block["dump_bytes"] > 0


def test_gromacs_scale_uses_the_lattice_count(work, monkeypatch):
    root = work / "md"
    _exe(root / "bin" / "mpirun", "#!/bin/sh\necho 'Time: 0.5 2.0'\nexit 0\n")
    _exe(root / "bin" / "gmx", "#!/bin/sh\nexit 0\n")
    monkeypatch.setenv("DSEAMS_GMX", str(root / "bin" / "gmx"))
    monkeypatch.setenv("DSEAMS_GMX_MPIRUN", str(root / "bin" / "mpirun"))
    rows = [{"id": str(i), "x": "1.5", "y": "1.0", "z": "1.0"} for i in range(1, 9)]
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: rows)
    monkeypatch.setattr(live, "scale_plan", lambda: [(1, 8, 2)])
    result = live.scale_gromacs()
    block = result["8"]["1"]
    assert block["steps"] == 2
    assert block["plain"] == [1.0] * 6
    assert block["warmup_live"] == 1.0
    assert len(block["live"]) == 6


def test_gromacs_trr_identity_requires_the_native_match(work, monkeypatch):
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [
        {"id": "1", "x": "1", "y": "1", "z": "1"}
    ])
    monkeypatch.setattr(live, "run", lambda *args, **kwargs: _Proc())
    monkeypatch.setattr(live, "trr_frames", lambda path: [(0, 1, (1,) * 9, (0, 0, 0))])
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 1)
    monkeypatch.setattr(live, "match_native_cli", lambda *args, **kwargs: 1)
    monkeypatch.setenv("DSEAMS_GMX", "gmx")
    monkeypatch.setenv("DSEAMS_GMX_MPIRUN", "mpirun")
    assert live.gromacs_trr_identity() == 1
    log = (live.WORK / "gromacs-trr" / "equal.log").read_text()
    assert log == "frames 1\n"
    monkeypatch.setattr(live, "match_native_cli", lambda *args, **kwargs: 0)
    with pytest.raises(RuntimeError, match="native trr matched 0"):
        live.gromacs_trr_identity()


def test_ladder_worker_prints_seven_samples(work, monkeypatch, capsys):
    dump = work / "frames.lammpstrj"
    _dump(dump, steps=((0, [(1, 0, 0, 0)]), (3, [(1, 0.2, 0, 0)])))

    class Eng(_Engine):
        seen = []

        def __init__(self, lib, path):
            Eng.seen = []

        def observe(self, step, pos, edge, dump, types=None):
            Eng.seen.append(step)

    import dseams_insitu.engine as engine

    monkeypatch.setattr(engine, "load", lambda path=None: object())
    monkeypatch.setattr(engine, "Engine", Eng)
    live.ladder_worker(dump)
    printed = [float(line) for line in capsys.readouterr().out.splitlines() if line.strip()]
    assert len(printed) == 7
    assert Eng.seen == [0, 3] * 7
    assert live.module_entry(["live.py", "--ladder-worker", str(dump)]) == 0


def test_thread_ladder_doubles_until_the_affinity(work, monkeypatch):
    dump = work / "one.lammpstrj"
    _dump(dump)
    empty = work / "empty.lammpstrj"
    empty.write_text("")
    with pytest.raises(RuntimeError, match="empty dump"):
        live.thread_ladder(empty, 1)

    class Proc:
        stdout = "\n".join(["0.2"] * 7) + "\n"
        stderr = ""
        returncode = 0

    monkeypatch.setattr(live, "run", lambda *args, **kwargs: Proc())
    out = live.thread_ladder(dump, 4)
    assert out["waters"] == 4 and out["frames"] == 1
    assert "1" in out["threads"]
    assert out["threads"]["1"]["samples"] == [0.2] * 6


def test_acceptance_reports_a_missing_engine_and_a_failure(monkeypatch, capsys):
    monkeypatch.delenv("SEAMS_INSITU_EXE", raising=False)
    monkeypatch.delenv("SEAMS_GRO_IMD_EXE", raising=False)
    assert live.acceptance() is None
    text = capsys.readouterr().out
    assert "not run: missing" in text
    assert "SEAMS_INSITU_EXE" in text
    monkeypatch.setenv("SEAMS_INSITU_EXE", "present")
    monkeypatch.setenv("SEAMS_GRO_IMD_EXE", "present")
    monkeypatch.setattr(live, "main", lambda: (_ for _ in ()).throw(RuntimeError("boom")))
    assert live.acceptance() is None
    assert "not run: boom" in capsys.readouterr().out
    called = {}
    monkeypatch.setattr(live, "main", lambda: called.setdefault("ok", True))
    assert live.acceptance() is None
    assert called["ok"] is True
    assert live.module_entry(["live.py", "--acceptance"]) == 0


def test_a_bare_module_entry_is_the_script(monkeypatch):
    import runpy

    argv = sys.argv[:]
    sys.argv = ["live.py"]
    try:
        with pytest.raises(SystemExit, match="pass --ladder-worker"):
            runpy.run_path(str(Path(live.__file__).resolve()), run_name="__main__")
    finally:
        sys.argv = argv


def test_frame_clock_ignores_a_preamble(tmp_path, monkeypatch):
    dump = tmp_path / "one.lammpstrj"
    _dump(dump)
    clock = tmp_path / "clock"
    _exe(clock, "#!/bin/sh\necho note\necho 'analysis_s 0.4'\n")
    monkeypatch.setenv("DSEAMS_FRAME_CLOCK", str(clock))
    monkeypatch.setattr(live, "cli_text", lambda *args, **kwargs: "")
    result = live.frame_clocks(dump, 0.3)
    assert result["cli_excluded"] == pytest.approx(0.4)
    assert result["insitu_per_frame"] == 0.3


def test_lammps_ranks_that_agree_replace_an_old_folder(work, monkeypatch):
    root = work / "mpi"
    _exe(
        root / "bin" / "mpirun",
        "#!/bin/sh\nprintf 'step 1\\nrings 3:0\\n' > live.txt\nexit 0\n",
    )
    _exe(root / "bin" / "lmp_mpi", "#!/bin/sh\nexit 0\n")
    plugin = work / "plugin.so"
    plugin.write_bytes(b"p")
    monkeypatch.setenv("DSEAMS_LMP_MPI_PREFIX", str(root))
    monkeypatch.setenv("DSEAMS_LMP_MPI_PLUGIN", str(plugin))
    stale = live.WORK / "mpi" / "1"
    stale.mkdir(parents=True)
    (stale / "old.txt").write_text("stale")
    same, note = live.lammps_mpi_identity()
    assert same is True
    assert note == "The census matches the one-rank run."
    assert not (stale / "old.txt").exists()
    log = (live.WORK / "mpi" / "equal.log").read_text()
    assert "1==2 True" in log and "1==4 True" in log


def test_gromacs_ranks(work, monkeypatch):
    monkeypatch.delenv("DSEAMS_GMX", raising=False)
    monkeypatch.delenv("DSEAMS_GMX_MPIRUN", raising=False)
    monkeypatch.setattr(live, "PREFIX", work / "no-gromacs")
    missing_same, missing_note = live.gromacs_mpi_identity()
    assert missing_same is False
    assert "absent" in missing_note
    root = work / "gro"
    _exe(root / "gmx", "#!/bin/sh\nexit 0\n")
    _exe(root / "mpirun", "#!/bin/sh\necho 'Time: 1 1'\nexit 0\n")

    def fake_run(cmd, cwd, env=None, check=True):
        text = " ".join(cmd)
        if "mdrun" in text:
            for part in cmd:
                if part.startswith("run-"):
                    (Path(cwd) / ("live-%s.txt" % part.split("-")[1])).write_text("step 1\n")
                    (Path(cwd) / (part + ".trr")).write_bytes(b"")
        return _Proc()

    monkeypatch.setattr(live, "run", fake_run)
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [
        {"id": "1", "x": "1", "y": "1", "z": "1"}
    ])
    monkeypatch.setattr(live, "trr_frames", lambda path: [object(), object()])
    monkeypatch.setattr(live, "match_native_cli", lambda *args, **kwargs: 2)
    monkeypatch.setenv("DSEAMS_GMX", str(root / "gmx"))
    monkeypatch.setenv("DSEAMS_GMX_MPIRUN", str(root / "mpirun"))
    stale = live.WORK / "gromacs-mpi"
    stale.mkdir(parents=True)
    (stale / "old.txt").write_text("stale")
    same, note = live.gromacs_mpi_identity()
    assert same is True
    assert "matches the trajectory" in note
    assert not (stale / "old.txt").exists()
    monkeypatch.setattr(live, "match_native_cli", lambda *args, **kwargs: 1)
    with pytest.raises(RuntimeError, match="matched 1 of 2"):
        live.gromacs_mpi_identity()

    def disagree(cmd, cwd, env=None, check=True):
        text = " ".join(cmd)
        if "mdrun" in text:
            for part in cmd:
                if part.startswith("run-"):
                    rank = part.split("-")[1]
                    (Path(cwd) / ("live-%s.txt" % rank)).write_text("step %s\n" % rank)
                    (Path(cwd) / (part + ".trr")).write_bytes(b"")
        return _Proc()

    monkeypatch.setattr(live, "run", disagree)
    monkeypatch.setattr(live, "match_native_cli", lambda *args, **kwargs: 2)
    import shutil

    shutil.rmtree(live.WORK / "gromacs-mpi")
    same, note = live.gromacs_mpi_identity()
    assert same is False
    assert "does not match" in note


def _install_ase(monkeypatch):
    import types

    import numpy as np

    ase = types.ModuleType("ase")
    calculators = types.ModuleType("ase.calculators")
    lj = types.ModuleType("ase.calculators.lj")
    md = types.ModuleType("ase.md")
    langevin = types.ModuleType("ase.md.langevin")
    units = types.ModuleType("ase.units")
    io = types.ModuleType("ase.io")
    units.fs = 1.0

    class Cell:
        def __init__(self, lengths):
            self.array = np.diag(np.asarray(lengths, dtype=float))

    class Atoms:
        def __init__(self, formula, positions, cell, pbc):
            self.positions = np.asarray(positions, dtype=float)
            self.cell = Cell(cell)
            self.calc = None

        def get_positions(self):
            return self.positions.copy()

        def get_atomic_numbers(self):
            return np.ones(len(self.positions), dtype=np.int32)

    class LJ:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class Dyn:
        def __init__(self, atoms, *args, **kwargs):
            self.atoms = atoms
            self.nsteps = 0
            self.callbacks = []

        def attach(self, callback, interval=1):
            self.callbacks.append(callback)

        def run(self, steps):
            for index in range(steps + 1):
                self.nsteps = index
                for callback in list(self.callbacks):
                    callback()

    class Image:
        def __init__(self, positions):
            self.positions = positions

        def get_positions(self):
            return self.positions

    class Trajectory:
        def __init__(self, path, mode):
            self.path = Path(path)
            self.frames = []

        def write(self, atoms):
            self.frames.append(atoms.get_positions().copy())

        def close(self):
            self.path.write_text("frames %d\n" % len(self.frames))

    def read(path, index=":"):
        count = int(Path(path).read_text().split()[1])
        return [Image(np.zeros((8, 3))) for _ in range(count)]

    lj.LennardJones = LJ
    langevin.Langevin = Dyn
    ase.Atoms = Atoms
    ase.units = units
    io.Trajectory = Trajectory
    io.read = read
    for name, module in (
        ("ase", ase),
        ("ase.calculators", calculators),
        ("ase.calculators.lj", lj),
        ("ase.md", md),
        ("ase.md.langevin", langevin),
        ("ase.units", units),
        ("ase.io", io),
    ):
        monkeypatch.setitem(sys.modules, name, module)


def test_ase_identity_and_both_timings(work, monkeypatch):
    _install_ase(monkeypatch)
    monkeypatch.setattr(live, "load", lambda path=None: object())
    monkeypatch.setattr(live, "Engine", _Engine)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 7)
    frames = live.ase_identity()
    assert frames == 7
    assert (live.WORK / "ase" / "from-traj.lammpstrj").is_file()
    plain = live.ase_time(4, False)
    timed = live.ase_time(4, True)
    assert plain >= 0 and timed >= 0
    assert (live.WORK / "ase" / "time.txt").is_file()


def test_ipi_success_and_a_failing_driver(work, monkeypatch):
    import dseams_insitu.ipi as ipi

    root = work / "md"
    _exe(root / "bin" / "i-pi", "#!/bin/sh\nexit 0\n")
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "load", lambda path=None: object())
    monkeypatch.setattr(live, "Engine", _Engine)

    def popen(cmd, cwd=None, **kwargs):
        folder = Path(cwd)
        (folder / "out.pos_0.xyz").write_text("1\ncomment\nO 0 0 0\n")
        return _Proc(0, "ipi done\n")

    monkeypatch.setattr(live.subprocess, "Popen", popen)

    def run_client(self, on_frame, max_frames=None):
        on_frame(0, __import__("numpy").zeros((8, 3)), __import__("numpy").eye(3) * 6.4)
        return 1

    monkeypatch.setattr(ipi.IpiClient, "run", run_client)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 1)
    assert live.run_ipi() == 1
    assert "ipi done" in (live.WORK / "ipi" / "ipi.out").read_text()
    assert (live.WORK / "ipi" / "input.xml").read_text().count("<total_steps>6</total_steps>") == 1

    def popen_bad(cmd, cwd=None, **kwargs):
        return _Proc(3, "failed\n")

    monkeypatch.setattr(live.subprocess, "Popen", popen_bad)
    with pytest.raises(RuntimeError, match="failed"):
        live.run_ipi()

    seen = []

    class Eng(_Engine):
        def observe(self, step, pos, box, dump, types=None):
            seen.append(dump)

    monkeypatch.setattr(live, "Engine", Eng)
    monkeypatch.setattr(live.subprocess, "Popen", popen)
    elapsed = live.ipi_elapsed(False)
    assert elapsed >= 0 and seen == []
    elapsed = live.ipi_elapsed(True)
    assert elapsed >= 0 and seen and str(seen[0]).endswith("dump.lammpstrj")


def test_main_removes_an_old_work_tree_and_seals_the_samples(work, monkeypatch):
    stale = live.WORK / "stale.txt"
    stale.parent.mkdir(parents=True)
    stale.write_text("old")
    monkeypatch.setattr(live, "LIB", str(work / "lib.so"))
    monkeypatch.setattr(live, "PLUGIN", work / "plugin.so")
    (work / "lib.so").write_bytes(b"lib")
    (work / "plugin.so").write_bytes(b"plugin")
    root = work / "md"
    (root / "bin").mkdir(parents=True)
    for name in ("lmp", "gmx"):
        (root / "bin" / name).write_bytes(b"bin")
    monkeypatch.setattr(live, "PREFIX", root)
    monkeypatch.setattr(live, "run_lammps", lambda: {"frames": 2, "steps": 4})
    monkeypatch.setattr(live, "time_lammps", lambda *args, **kwargs: [0.25])
    monkeypatch.setattr(live, "run_gromacs", lambda: {"frames": 3, "steps": 5})
    monkeypatch.setattr(live, "time_gromacs", lambda *args, **kwargs: ([0.1], [0.2]))
    monkeypatch.setattr(live, "openmm_identity", lambda: 4)
    monkeypatch.setattr(live, "openmm_time", lambda steps: (0.1, 0.2))
    monkeypatch.setattr(live, "ase_identity", lambda: 5)
    monkeypatch.setattr(live, "ase_time", lambda *args, **kwargs: 0.3)
    monkeypatch.setattr(live, "run_ipi", lambda: 6)
    monkeypatch.setattr(live, "ipi_elapsed", lambda analyse: 0.4)
    monkeypatch.setattr(
        live,
        "scale_lammps",
        lambda: {"8": {"10": {"live": [1.5, 1.5], "plain": [1.0, 1.0]}}},
    )
    monkeypatch.setattr(live, "gromacs_trr_identity", lambda: 7)
    monkeypatch.setattr(live, "scale_gromacs", lambda: {})
    monkeypatch.setattr(live, "frame_clocks", lambda dump, per_frame: {"cli_excluded": per_frame})
    monkeypatch.setattr(live, "lammps_bin", lambda *args, **kwargs: None)
    monkeypatch.setattr(live, "lammps_plain", lambda *args, **kwargs: None)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *args, **kwargs: 9)
    monkeypatch.setattr(live, "thread_ladder", lambda dump, natoms: {"waters": natoms, "threads": {}})
    monkeypatch.setattr(live, "lammps_mpi_identity", lambda: (True, "matched"))
    monkeypatch.setattr(live, "gromacs_mpi_identity", lambda: (True, "matched"))
    monkeypatch.setattr(live, "diamond_positions", lambda ncell: [{"id": "1"}])
    monkeypatch.setenv("DSEAMS_SCALE_TINY", "1")
    live.main()
    saved = live.WORK / "samples.json"
    text = saved.read_text()
    assert "matched" in text
    assert '"frames": 9' in text or '"frames": 2' in text
    assert not stale.exists()
    assert (live.WORK / "equality-lammps" / "equal.log").read_text() == "frames 9\n"
    lattice = live.WORK / "lattice" / "5" / "one.lammpstrj"
    lattice.parent.mkdir(parents=True, exist_ok=True)
    lattice.write_text("ITEM: TIMESTEP\n0\n")
    sealed = {
        "scale": {"8": {"10": {"live": [1.5], "plain": [1.0]}}},
        "codes": {
            "openmm": {"frames": 4},
            "ase": {"frames": 5},
            "ipi": {"frames": 6},
        },
    }
    live.seal_tail(sealed, 7)
    assert sealed["clock"]["cli_excluded"] == 0.5
    assert sealed["ladder"][0]["waters"] == 8
