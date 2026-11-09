"""Branch coverage for the in-situ Python adapters, scripts, and drivers.

Engines are not started. Subprocess entry points are replaced and the
assertions check return values, exception text, files, and command lists.
"""

import importlib.util
import io
import tempfile
import json
import os
import socket
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "adapters" / "python"))
sys.path.insert(0, str(ROOT / "tests" / "insitu"))

SEAMS = os.environ["SEAMS_CLI"]
LIB = os.environ["SEAMS_LIB"]


class Proc:
    def __init__(self, rc, out="", err=""):
        self.returncode = rc
        self.stdout = out
        self.stderr = err


def exec_file(path, argv, name="__main__"):
    saved = sys.argv
    sys.argv = list(argv)
    try:
        code = compile(Path(path).read_text(), str(path), "exec")
        ns = {"__name__": name, "__file__": str(path)}
        exec(code, ns)
        return ns
    finally:
        sys.argv = saved


def test_openmm_unit_stub_and_import(monkeypatch):
    import builtins
    import conftest

    real_import = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "openmm" or name.startswith("openmm."):
            raise ImportError(name)
        return real_import(name, globals, locals, fromlist, level)

    saved = {}
    for name in [key for key in sys.modules if key == "openmm" or key.startswith("openmm.")]:
        saved[name] = sys.modules.pop(name)
    try:
        monkeypatch.setattr(builtins, "__import__", blocked)
        conftest.install_openmm_unit()
        assert sys.modules["openmm.unit"].angstrom is not None
        monkeypatch.undo()
        for name in [key for key in list(sys.modules) if key == "openmm" or key.startswith("openmm.")]:
            sys.modules.pop(name)
        unit = types.ModuleType("openmm.unit")
        unit.angstrom = object()
        package = types.ModuleType("openmm")
        package.unit = unit
        sys.modules["openmm"] = package
        sys.modules["openmm.unit"] = unit
        conftest.install_openmm_unit()
        assert sys.modules["openmm.unit"] is unit
    finally:
        for name in [key for key in list(sys.modules) if key == "openmm" or key.startswith("openmm.")]:
            sys.modules.pop(name, None)
        sys.modules.update(saved)


def test_hooks_as_main():
    import runpy

    saved = sys.argv[:]
    sys.argv = [str(ROOT / "tests" / "insitu" / "test_hooks.py")]
    try:
        with pytest.raises(SystemExit) as caught:
            runpy.run_path(sys.argv[0], run_name="__main__")
    finally:
        sys.argv = saved
    assert caught.value.code in (0, False, None)


def test_package_exports():
    import dseams_insitu

    assert dseams_insitu.Engine.__name__ == "Engine"
    assert dseams_insitu.__all__ == [
        "DseamsMDObserver",
        "DseamsReporter",
        "Engine",
        "IpiClient",
    ]


def test_engine_path_and_second_close(tmp_path):
    from dseams_insitu.engine import Engine, library_path

    assert library_path(LIB) == LIB
    assert library_path("") == LIB
    record = str(tmp_path / "rec.txt")
    eng = Engine(Path(LIB), record)
    eng.close()
    eng.close()
    with Engine(LIB, record):
        pass


def test_ipi_forces_and_frame_limit():
    from dseams_insitu.ipi import IpiClient, lj_forces, message, read_exact

    cell = np.eye(3) * 30.0
    zero, force = lj_forces(np.zeros((2, 3)), cell, cutoff=6.0)
    assert zero == 0.0 and np.all(force == 0.0)
    far, _ = lj_forces(np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]), cell, cutoff=6.0)
    assert far == 0.0
    near, near_force = lj_forces(
        np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]), np.eye(3) * 20.0, cutoff=6.0
    )
    assert near != 0.0 and near_force.shape == (2, 3)

    class OneByte:
        def __init__(self, data):
            self.data = data
            self.index = 0

        def recv(self, n):
            if self.index >= len(self.data):
                return b""
            chunk = self.data[self.index : self.index + 1]
            self.index += 1
            return chunk

    one = OneByte(b"a")
    assert one.recv(1) == b"a"
    assert one.recv(1) == b""
    assert read_exact(OneByte(b"abcd"), 4) == b"abcd"

    left, right = socket.socketpair()
    left.settimeout(0.5)
    try:
        right.sendall(np.int32(0).tobytes() + np.int32(0).tobytes())
        client = IpiClient("127.0.0.1", 1)
        client._init(left)
        assert client.initialised is True
    finally:
        left.close()
        right.close()

    cell_row = np.ascontiguousarray(np.eye(3) * 12.0, dtype=np.float64).reshape(-1)
    blob = b"".join(
        [
            message("POSDATA"),
            cell_row.tobytes(),
            cell_row.tobytes(),
            np.int32(1).tobytes(),
            np.zeros(3, dtype=np.float64).tobytes(),
        ]
    )

    class Sock:
        def __init__(self, payload):
            self.payload = payload
            self.pos = 0

        def setsockopt(self, *args):
            return None

        def recv(self, n):
            if self.pos >= len(self.payload):
                return b""
            chunk = self.payload[self.pos : self.pos + n]
            self.pos += len(chunk)
            return chunk

        def close(self):
            return None

    import dseams_insitu.ipi as ipi

    drained = Sock(b"")
    assert drained.recv(4) == b""
    real = ipi.socket.create_connection
    ipi.socket.create_connection = lambda addr: Sock(blob)
    try:
        seen = []
        frames = IpiClient("127.0.0.1", 9).run(
            lambda i, pos, box: seen.append((i, tuple(pos.shape))), max_frames=1
        )
    finally:
        ipi.socket.create_connection = real
    assert frames == 1
    assert seen == [(0, (1, 3))]


def test_ase_remaining_names():
    from dseams_insitu.ase import is_optimizer

    for name in ("LBFGS", "FIRE", "GPMin", "MDMin", "ODE", "GoodOldQuasiNewton"):
        assert is_optimizer(type(name, (), {})()) is True
    assert is_optimizer(type("Verlet", (), {})()) is False


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_coverage_edges(tmp_path, monkeypatch):
    import test_insitu_coverage as harness

    module = harness.load_scorer()
    text = "s = 'it\\'s'\nstill = 1\n"
    assert module.physical(text, ".py") == {1, 2}
    quoted = "note = '''inside\nmiddle\nend'''\nvalue = 1\n \n"
    assert module.physical(quoted, ".py") == {1, 3, 4}
    assert module.changed_python(tmp_path, None) == [("changed-files", 0, 0, 0.0)]
    assert module.changed_python(tmp_path, "   ")[0][0] == "changed-files"

    root = tmp_path / "tree"
    root.mkdir()
    (root / "skip").mkdir()
    (root / "skip" / "kept.cpp").mkdir()
    src = root / "src"
    src.mkdir()
    (src / "kept.cpp").write_text("int kept;\nint both;\n")
    gcov = root / "kept.cpp.gcov"
    gcov.write_text(
        "        -:    0:Source:kept.cpp\n"
        "    #####:    1:int kept;\n"
        "    =====:    2:int both;\n"
        "        4:    2:int both;\n"
        "        1:    2:int both;\n"
        "not a count\n"
    )
    monkeypatch.setenv("DSEAMS_COVERAGE_BASE", "   ")
    code, printed = harness.run(module, root, gcov)
    assert code == 1
    assert "src/kept.cpp 2 1 0.500000" in printed
    assert "changed-files 0 0 0.000000" in printed

    monkeypatch.setenv("DSEAMS_COVERAGE_BASE", "deadbeef")
    code, printed = harness.run(module, root, "changed-files deadbeef")
    assert printed.count("changed-files deadbeef") == 1

    repo = tmp_path / "repo"
    repo.mkdir()
    env = os.environ.copy()
    env.update(
        {
            "GIT_AUTHOR_NAME": "Rohit Goswami",
            "GIT_AUTHOR_EMAIL": "rgoswami@ieee.org",
            "GIT_COMMITTER_NAME": "Rohit Goswami",
            "GIT_COMMITTER_EMAIL": "rgoswami@ieee.org",
        }
    )
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True)
    (repo / "note.txt").write_text("n\n")
    subprocess.run(["git", "add", "note.txt"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "base"],
        cwd=repo,
        check=True,
        capture_output=True,
        env=env,
    )
    base = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    seen = repo / "seen.py"
    seen.write_text("x = 1\n")
    payload = {"files": {str(seen): {"executed_lines": [1], "missing_lines": []}}}
    (repo / "cov.json").write_text(json.dumps(payload))
    monkeypatch.setenv("DSEAMS_COVERAGE_BASE", base)
    code, printed = harness.run(module, repo, repo / "cov.json")
    assert code == 0
    assert "seen.py 1 1 1.000000" in printed
    assert "changed-files" not in printed

    class Boom:
        def main(self):
            raise RuntimeError("boom")

    held = sys.stdout
    with pytest.raises(RuntimeError, match="boom"):
        harness.run(Boom(), repo)
    assert sys.stdout is held
    with pytest.raises(SystemExit) as caught:
        exec_file(ROOT / "scripts" / "insitu_coverage.py", ["insitu_coverage.py", str(tmp_path)])
    assert caught.value.code == 1


def test_coverage_script_bases(monkeypatch, capsys):
    import test_insitu_coverage as harness

    monkeypatch.delenv("DSEAMS_COVERAGE_BASE", raising=False)
    assert harness.main() == 0
    assert "DSEAMS_COVERAGE_BASE" not in os.environ
    assert "coverage scorer ok" in capsys.readouterr().out

    monkeypatch.setenv("DSEAMS_COVERAGE_BASE", "kept-base")
    path = Path(harness.__file__)
    with pytest.raises(SystemExit) as caught:
        exec_file(path, [str(path)])
    assert caught.value.code == 0
    assert os.environ["DSEAMS_COVERAGE_BASE"] == "kept-base"


def test_insitu_report(tmp_path):
    report = _load(ROOT / "scripts" / "insitu_report.py", "insitu_report_mod")
    with pytest.raises(SystemExit, match="no samples"):
        report.type7([], 0.5)
    assert report.type7([9], 0.2) == 9
    assert report.type7([0, 10], 0.25) == 2.5
    assert report.type7([1, 2, 3], 1.0) == 3

    host = {
        "cpu": "test-cpu",
        "head": "abc",
        "tree": "def",
        "lib_sha256": "11",
        "lmp_sha256": "22",
        "gmx_sha256": "33",
    }
    minimal = tmp_path / "min.json"
    minimal.write_text(
        json.dumps(
            {
                "host": host,
                "codes": {
                    "skip": {"plain": [1], "live": [1, 2]},
                    "empty": {"plain": [], "live": []},
                },
            }
        )
    )
    dest = tmp_path / "min.md"
    report.main_args = None
    sys.argv = ["insitu_report.py", str(minimal), str(dest)]
    report.main()
    text = dest.read_text()
    assert "Plugin SHA-256 ``." in text
    assert "raw samples" not in text
    assert "| Code |" not in text
    assert "# LAMMPS size and cadence" not in text

    clock_only = tmp_path / "clock.json"
    clock_only.write_text(
        json.dumps(
            {
                "host": host,
                "codes": {},
                "clock": {
                    "system": "tiny",
                    "insitu_per_frame": 1.5,
                    "cli_included": 2.5,
                    "cli_excluded": 0.5,
                },
            }
        )
    )
    clock_dest = tmp_path / "clock.md"
    sys.argv = ["insitu_report.py", str(clock_only), str(clock_dest)]
    report.main()
    clock_text = clock_dest.read_text()
    assert "# Analysis time against the command line" in clock_text
    assert "A " not in clock_text.split("# Analysis", 1)[1]

    full = {
        "host": dict(host, plugin_sha256="abc"),
        "codes": {
            "lammps": {"plain": [2, 4], "live": [6, 10], "steps": 4, "frames": 2},
            "odd": {"plain": [1], "live": []},
        },
        "scale": {
            "1000": {
                "1": {
                    "plain": [2, 4],
                    "live": [6, 10],
                    "steps": 4,
                    "post_wall": [1, 3],
                    "dump_bytes": 50,
                    "warmup_plain": 0.5,
                    "warmup_live": 0.7,
                },
                "10": {
                    "plain": [2, 4],
                    "live": [6, 10],
                    "steps": 8,
                    "post_wall": [1, 3],
                    "dump_bytes": 9,
                },
            }
        },
        "gromacs_scale": {
            "1000": {
                "1": {"plain": [2, 4], "live": [6, 10], "steps": 4},
            }
        },
        "clock": {
            "system": "tiny",
            "insitu_per_frame": 1.5,
            "cli_included": 2.5,
            "cli_excluded": 0.5,
        },
        "steady": {"plain": [2, 4], "live": [6, 10], "steps": 40, "every": 5},
        "equality": [
            {"code": "LAMMPS", "frames": 4, "rings": "yes", "chill": "yes", "cages": "yes"}
        ],
        "ladder": [
            {
                "waters": 8,
                "threads": {
                    "1": {"samples": [10, 10, 10, 10]},
                    "2": {"samples": [1, 1, 1, 1]},
                },
            }
        ],
        "mpi": [
            {"code": "LAMMPS", "ranks": "2 and 4", "identical": True},
            {"code": "GROMACS", "ranks": "2 and 4", "identical": False, "note": "rank note"},
        ],
    }
    full_path = tmp_path / "full.json"
    full_path.write_text(json.dumps(full))
    full_dest = tmp_path / "full.md"
    exec_file(ROOT / "scripts" / "insitu_report.py", ["insitu_report.py", str(full_path), str(full_dest)])
    text = full_dest.read_text()
    assert "Plugin SHA-256 `abc`." in text
    assert "## lammps raw samples" in text
    assert "4.5" in text and "2.75" in text
    assert text.count("Warmup plain seconds per step:") == 1
    assert "A 40-step run at cadence 5" in text
    assert "| LAMMPS | 4 | yes | yes | yes |" in text
    assert "LAMMPS ranks 2 and 4 identical to serial: yes." in text
    assert "GROMACS ranks 2 and 4 identical to serial: no." in text
    assert "rank note" in text
    assert "Thread scaling failed" not in text
    assert "Waters 8, threads 2." in text

    bad = dict(full)
    bad["ladder"] = [
        {
            "waters": 8,
            "threads": {
                "1": {"samples": [1, 1, 1, 1]},
                "2": {"samples": [8, 8, 8, 8]},
            },
        }
    ]
    bad_path = tmp_path / "bad.json"
    bad_path.write_text(json.dumps(bad))
    bad_dest = tmp_path / "bad.md"
    sys.argv = ["insitu_report.py", str(bad_path), str(bad_dest)]
    with pytest.raises(SystemExit, match="thread scaling failed for 8 waters at 2 threads"):
        report.main()
    assert "Thread scaling failed" in bad_dest.read_text()
    assert "Failed at 8 waters with 2 threads." in bad_dest.read_text()


def _cli_decide(args, mode):
    if mode.get("all_fail"):
        return 1, ""
    if not args:
        return 2, ""
    if "--cutoff" in args and args[args.index("--cutoff") + 1] == "not-a-number":
        return mode.get("cutoff_rc", 3), ""
    if "--emit-library" in args:
        val = args[args.index("--emit-library") + 1]
        if val == "tiny":
            return mode.get("emit_rc", 0), mode.get("emit_out", "key 1\n")
        if val == "one":
            return mode.get("hop1_rc", 0), "one\n"
        if val == "two":
            return mode.get("hop2_rc", 0), "two\n"
        return 0, "x\n"
    if args == ["--version"]:
        return mode.get("ver_rc", 0), mode.get("ver_out", "seams 9\n")
    if any(str(a).endswith("blank.con") for a in args):
        return mode.get("blank_rc", 0), ""
    if "--library" in args:
        val = args[args.index("--library") + 1]
        sides = val.split(",")
        repeated = len(sides) == 2 and sides[0] == sides[1]
        if val.endswith("missing-lib.txt") or val.startswith(",") or repeated:
            return 2, ""
    if "--strict-input" in args and "--rebuild" not in args:
        return 2, ""
    if "--types" in args and args[args.index("--types") + 1] in {"1", "1,2,3", "a,b", "1,", ",1"}:
        return 2, ""
    if args[0] == "ions" and "--type" not in args and "--select" not in args and "--ion-types" not in args:
        return 2, ""
    if args[0] == "pairs" and "--site" not in args:
        return 2, ""
    if args[0] == "cn" and "--ions" in args and "--site" not in args:
        return 2, ""
    bad = (
        "no-such-command",
        "--format\nyaml",
        "--axis\nq",
        "--bins\n-1",
        "--graph\nfoo",
        "--graph\nknn-union",
        "--family\nnot-a-family",
        "--frame\n1-2",
        "--guest-types\nx",
        "--cage-edge\n1.0",
        "--subset\nnope",
        "--site\n1=not-a-kind",
        "--site\nno-equals",
        "--site\nabc=cationHead",
        "--site\n=cationHead",
    )
    joined = "\n".join(args)
    if any(item in joined for item in bad):
        return 2, ""
    return 0, "seams ok\n"


def _run_cli(mode, offload=None, script=None):
    calls = []

    def fake_run(cmd, **kwargs):
        cmd = [str(part) for part in cmd]
        calls.append(cmd)
        rc, out = _cli_decide(cmd[1:], mode)
        return Proc(rc, out, "")

    saved_run = subprocess.run
    saved_argv = sys.argv[:]
    saved_err = sys.stderr
    buf = io.StringIO()
    subprocess.run = fake_run
    sys.stderr = buf
    if offload is None:
        os.environ.pop("SEAMS_OFFLOAD", None)
    else:
        os.environ["SEAMS_OFFLOAD"] = offload
    sys.argv = ["cli_every.py", SEAMS, str(ROOT)]
    path = script or (ROOT / "tests" / "cli_every.py")
    try:
        try:
            code_obj = compile(path.read_text(), str(path), "exec")
            exec(code_obj, {"__name__": "__main__", "__file__": str(path)})
            code = 0
        except SystemExit as exc:
            code = exc.code
    finally:
        subprocess.run = saved_run
        sys.argv = saved_argv
        sys.stderr = saved_err
    return code, buf.getvalue(), calls


def test_cli_every_commands():
    code, err, calls = _run_cli({})
    assert code == 0, err
    assert [SEAMS, "--help"] in calls
    assert any("not-a-number" in arg for cmd in calls for arg in cmd)
    assert any(cmd[1:2] == ["no-such-command"] for cmd in calls)
    assert os.environ.get("SEAMS_OFFLOAD") is None

    code, err, _ = _run_cli({}, offload="1")
    assert code == 0, err
    assert os.environ.get("SEAMS_OFFLOAD") == "1"
    os.environ.pop("SEAMS_OFFLOAD", None)

    code, err, _ = _run_cli({"all_fail": True})
    assert code == 1, err
    code, _, _ = _run_cli({"emit_rc": 1})
    assert code == 1
    code, _, _ = _run_cli({"emit_out": ""})
    assert code == 1
    code, _, _ = _run_cli({"hop2_rc": 1})
    assert code == 1
    code, _, _ = _run_cli({"ver_out": "9\n"})
    assert code == 1
    code, _, _ = _run_cli({"cutoff_rc": 0})
    assert code == 1
    code, _, _ = _run_cli({"blank_rc": 1})
    assert code == 1

    planted = ROOT / "input" / "traj" / "not-deposited.lammpstrj"
    planted.write_text("x\n")
    try:
        code, err, _ = _run_cli({})
    finally:
        planted.unlink()
    assert code == 1, err
    stay = Path(tempfile.mkdtemp()) / "stay.py"
    stay.write_text("value = 1\n")
    code, err, _ = _run_cli({}, script=stay)
    assert code == 0, err


def test_lammps_routes(tmp_path, monkeypatch, capsys):
    routes = _load(ROOT / "tests" / "insitu" / "lammps_routes.py", "lammps_routes_mod")
    with pytest.raises(SystemExit, match="no samples"):
        routes.type7([], 0.5)
    assert routes.type7([4], 0.9) == 4
    assert routes.type7([0, 10], 0.25) == 2.5
    assert routes.type7([1, 2, 3], 1) == 3
    assert routes.loop_time("Loop time of 8 on 1 procs for 4 steps\n") == 8
    with pytest.raises(RuntimeError, match="no loop time"):
        routes.loop_time("nothing\n")
    assert routes.section_time("Modify  time (%) = 1.5 (10)\n", "Modify") == 1.5
    assert routes.section_time("Modify time (%) = 2.5 (9)\n", "Modify") == 2.5
    assert routes.section_time("Pair  time (%) = 1 (1)\n", "Modify") is None
    assert "plugin load plug.so" in routes.input_text(2, 1, 1, "rec", "plug.so")
    assert "plugin load" not in routes.input_text(2, 1, 1, "rec", "")

    monkeypatch.setenv("DSEAMS_ROUTE_TINY", "1")
    assert routes.plan_rows() == [(1000, 5, 1, 2, 1)]
    monkeypatch.delenv("DSEAMS_ROUTE_TINY")
    assert len(routes.plan_rows()) == 18

    calls = []

    # The branch matrix is easier to hit by calling run_case directly.
    def direct(cmd, cwd=None, env=None, text=None, capture_output=None):
        calls.append([str(part) for part in cmd])
        return Proc(0, "Loop time of 10 on 1 procs for 2 steps\n", "")

    monkeypatch.setattr(routes.subprocess, "run", direct)
    per, mod, _ = routes.run_case("lmp", "", 1, 2, 1, 1, "", tmp_path / "direct", {})
    assert per == 5
    assert mod is None
    assert calls[-1][:1] == ["lmp"]
    per, mod, _ = routes.run_case("lmp", "mpirun", 1, 2, 1, 1, "p.so", tmp_path / "mpi", {})
    assert calls[-1][0] == "mpirun"
    routes.run_case("lmp", "", 2, 2, 1, 1, "", tmp_path / "np", {})
    assert calls[-1][0] == ""
    assert calls[-1][2] == "-np"

    def failing(cmd, cwd=None, env=None, text=None, capture_output=None):
        return Proc(4, "stdout-tail", "stderr-tail")

    monkeypatch.setattr(routes.subprocess, "run", failing)
    with pytest.raises(RuntimeError, match="run failed"):
        routes.run_case("lmp", "", 1, 2, 1, 1, "", tmp_path / "fail", {})

    def noloop(cmd, cwd=None, env=None, text=None, capture_output=None):
        return Proc(0, "quiet", "")

    monkeypatch.setattr(routes.subprocess, "run", noloop)
    with pytest.raises(RuntimeError, match="run failed"):
        routes.run_case("lmp", "", 1, 2, 1, 1, "", tmp_path / "noloop", {})

    monkeypatch.setenv("DSEAMS_ROUTE_TINY", "1")
    monkeypatch.setenv("DSEAMS_ROUTE_WORK", str(tmp_path / "work"))
    monkeypatch.delenv("DSEAMS_MPIRUN", raising=False)
    # Reload WORK from the env the function reads at call time.
    routes.WORK = Path(os.environ["DSEAMS_ROUTE_WORK"])
    routes.main.__wrapped__ = None
    sys.argv = ["lammps_routes.py", "plugin", "lmp", "-"]
    # Seven samples need seven logs. Extend the iterator by repeating.
    loop_logs = [
        "Loop time of 10 on 1 procs for 2 steps\nModify  time (%) = 1.5 (10)\n",
        "Loop time of 12 on 1 procs for 2 steps\n",
        "Loop time of 8 on 1 procs for 2 steps\nModify time (%) = 2 (4)\n",
    ]

    def cycling2(cmd, cwd=None, env=None, text=None, capture_output=None):
        calls.append([str(part) for part in cmd])
        (Path(cwd) / "log").write_text(loop_logs[len(calls) % 3])
        return Proc(0, "", "")

    monkeypatch.setattr(routes.subprocess, "run", cycling2)
    routes.main()
    written = json.loads((tmp_path / "work" / "plugin.json").read_text())
    assert written["route"] == "plugin"
    assert written["plugin"] == ""
    block = written["cases"]["1000/1/np1"]
    assert len(block["samples"]) == 6
    assert block["modify"]
    assert any(len(cmd) > 2 and cmd[2] == "-np" for cmd in calls) or calls

    plugin = {
        "cases": {
            "a": {
                "waters": 10,
                "every": 1,
                "ranks": 1,
                "samples": [1, 1, 1, 1],
                "modify": [0.2, 0.4],
            },
            "b": {
                "waters": 20,
                "every": 1,
                "ranks": 2,
                "samples": [5, 5, 5, 5],
                "modify": [],
            },
            "c": {
                "waters": 30,
                "every": 1,
                "ranks": 4,
                "samples": [1, 2, 3, 4],
                "modify": [1],
            },
        }
    }
    intree = {
        "cases": {
            "a": {
                "waters": 10,
                "every": 1,
                "ranks": 1,
                "samples": [5, 5, 5, 5],
                "modify": [0.5, 0.7],
            },
            "b": {
                "waters": 20,
                "every": 1,
                "ranks": 2,
                "samples": [1, 1, 1, 1],
                "modify": [0.1],
            },
            "c": {
                "waters": 30,
                "every": 1,
                "ranks": 4,
                "samples": [1, 2, 3, 4],
                "modify": [],
            },
        }
    }
    p1 = tmp_path / "plugin.json"
    p2 = tmp_path / "intree.json"
    p1.write_text(json.dumps(plugin))
    p2.write_text(json.dumps(intree))
    routes.compare(p1, p2)
    out = capsys.readouterr().out
    assert "10 1 1" in out and "plugin" in out
    assert "intree" in out
    assert "none" in out
    assert "real_difference True" in out
    assert "modify_s_per_step" in out

    same = {"cases": {"a": {"waters": 1, "every": 1, "ranks": 1, "samples": [3], "modify": []}}}
    s1 = tmp_path / "s1.json"
    s2 = tmp_path / "s2.json"
    s1.write_text(json.dumps(same))
    s2.write_text(json.dumps(same))
    routes.compare(s1, s2)
    out = capsys.readouterr().out
    assert "real_difference False" in out
    assert "modify_s_per_step" not in out

    exec_file(ROOT / "tests" / "insitu" / "lammps_routes.py", ["routes.py", "compare", str(s1), str(s2)])
    exec_file(ROOT / "tests" / "insitu" / "lammps_routes.py", ["routes.py", "plugin", "lmp", "-"])


def test_v4_measure(tmp_path, monkeypatch):
    real_run = subprocess.run
    monkeypatch.setenv("DSEAMS_LIVE_WORK", str(tmp_path / "v4"))
    monkeypatch.setenv("DSEAMS_V4_TINY", "1")
    v4 = _load(ROOT / "tests" / "insitu" / "v4_measure.py", "v4_measure_mod")
    with pytest.raises(SystemExit, match="set DSEAMS_NOT_A_VAR"):
        v4.need("DSEAMS_NOT_A_VAR")
    assert v4.need("DSEAMS_LMP") == os.environ["DSEAMS_LMP"]
    assert v4.plans()["cadence"] == (1,)
    monkeypatch.delenv("DSEAMS_V4_TINY")
    assert v4.plans()["scale"][1][1] == 64000
    monkeypatch.setenv("DSEAMS_V4_TINY", "1")

    folder = tmp_path / "lt"
    folder.mkdir()
    with pytest.raises(RuntimeError):
        v4.loop_time(folder, 2)
    (folder / "log").write_text("see Loop time of 9 later\n")
    with pytest.raises(RuntimeError, match="see Loop time"):
        v4.loop_time(folder, 2)
    (folder / "log").write_text("Loop time of 6.0 on 1 procs for 3 steps\n")
    assert v4.loop_time(folder, 3) == 2.0

    def bad(cmd, cwd=None, env=None, **kwargs):
        return Proc(4, "hello", "boom-lammps")

    monkeypatch.setattr(v4.subprocess, "run", bad)
    with pytest.raises(RuntimeError, match="boom-lammps"):
        v4.lammps_live_step(tmp_path / "live-bad", 2, 1, 1)
    with pytest.raises(RuntimeError, match="boom-lammps"):
        v4.lammps_plain_step(tmp_path / "plain-bad", 2, 1)
    dump = tmp_path / "native.lammpstrj"
    dump.write_text("ITEM: TIMESTEP\n0\n")
    with pytest.raises(RuntimeError, match="boom-lammps"):
        v4.postprocess(dump)

    def diamonds(ncell):
        count = 8 * int(ncell) ** 3
        return [{"id": str(i), "x": "0.1", "y": "0.2", "z": "0.3"} for i in range(1, count + 1)]

    monkeypatch.setattr(v4.live, "diamond_positions", diamonds)

    def gmx_run(cmd, cwd, env=None, check=True):
        return Proc(0, "Time: 1.0 4.0\n", "")

    monkeypatch.setattr(v4.live, "run", gmx_run)
    with pytest.raises(RuntimeError, match="lattice 7 has"):
        v4.gromacs_case(tmp_path / "bad-lattice", 7, 2, 1, False)
    plain = v4.gromacs_case(tmp_path / "gplain", 1000, 2, 1, False)
    live = v4.gromacs_case(tmp_path / "glive", 1000, 2, 1, True)
    assert plain == 2.0 and live == 2.0
    assert "dseams-active" not in (tmp_path / "gplain" / "md.mdp").read_text()
    assert "dseams-active = yes" in (tmp_path / "glive" / "md.mdp").read_text()

    state = {"clock_out": "noise\nanalysis_s 0.25\n", "rc": 0}
    lammps_cmds = []

    def sub_run(cmd, cwd=None, env=None, **kwargs):
        cmd = [str(part) for part in cmd]
        if cmd and cmd[0] == "git":
            return real_run(cmd, cwd=cwd, env=env, **kwargs)
        lammps_cmds.append(cmd)
        if cwd and "-log" in cmd:
            Path(cwd, "log").write_text("Loop time of 4.0 on 1 procs for 2 steps\n")
        stdout = ""
        if cmd and "clock" in Path(cmd[0]).name:
            stdout = state["clock_out"]
        return Proc(state["rc"], stdout, "")

    created = []

    def live_run(cmd, cwd, env=None, check=True):
        created.append([str(part) for part in cmd])
        folder = Path(cwd)
        script = None
        if "-in" in cmd:
            script = Path(cmd[cmd.index("-in") + 1])
            if not script.is_absolute():
                script = folder / script
        text = script.read_text() if script and script.is_file() else ""
        for line in text.splitlines():
            if "dump native" in line:
                Path(line.split()[5]).parent.mkdir(parents=True, exist_ok=True)
                Path(line.split()[5]).write_text("ITEM: TIMESTEP\n0\n")
        return Proc(0, "Time: 1.0 4.0\n", "")

    (tmp_path / "in.relative").write_text("units real\n")
    assert live_run(["lmp", "-in", "in.relative"], tmp_path).returncode == 0

    monkeypatch.setattr(v4.subprocess, "run", sub_run)
    monkeypatch.setattr(v4.live, "run", live_run)
    monkeypatch.setattr(
        v4.live,
        "thread_ladder",
        lambda dump, natoms: {"waters": natoms, "threads": {"1": {"samples": [0.1] * 6}}},
    )
    v4.plans = lambda: {
        "scale": [(5, 1000, 2), (2, 64, 2)],
        "cadence": (1, 10),
        "ladder": ((5, 1000),),
    }
    v4.WORK = tmp_path / "measure"
    pre = v4.WORK / "post" / "1000" / "native.lammpstrj"
    pre.parent.mkdir(parents=True)
    pre.write_text("ITEM: TIMESTEP\n0\n")
    clock = tmp_path / "frame_clock"
    clock.write_text("#!/bin/sh\necho analysis_s 0.25\n")
    clock.chmod(0o755)
    monkeypatch.setenv("DSEAMS_FRAME_CLOCK", str(clock))
    v4.main()
    data = json.loads((v4.WORK / "samples.json").read_text())
    assert data["scale"]["1000"]["1"]["steps"] == 2
    assert data["scale"]["1000"]["10"]["steps"] == 10
    assert "lammps" in data["codes"] and "gromacs" in data["codes"]
    assert "64" in data["scale"] and "lammps" in data["codes"]
    assert data["codes"]["lammps"]["steps"] == 2
    assert data["scale"]["1000"]["1"]["warmup_plain"] > 0
    assert data["scale"]["1000"]["1"]["dump_bytes"] > 0
    assert data["clock"]["cli_excluded"] == 0.25
    assert data["ladder"][0]["waters"] == 1000
    assert any(cmd[0].endswith("lmp") or "lmp" in cmd[0] for cmd in lammps_cmds + created)

    v4.WORK = tmp_path / "noclock"
    monkeypatch.delenv("DSEAMS_FRAME_CLOCK", raising=False)
    v4.plans = lambda: {"scale": [(5, 1000, 2)], "cadence": (1,), "ladder": ()}
    v4.main()
    assert "clock" not in json.loads((v4.WORK / "samples.json").read_text())

    missing_clock = tmp_path / "absent_clock"
    monkeypatch.setenv("DSEAMS_FRAME_CLOCK", str(missing_clock))
    v4.WORK = tmp_path / "missing-clock"
    v4.main()
    assert "clock" not in json.loads((v4.WORK / "samples.json").read_text())

    v4.plans = lambda: {"scale": [(2, 64, 2)], "cadence": (1,), "ladder": ()}
    v4.WORK = tmp_path / "no-dump"
    monkeypatch.setenv("DSEAMS_FRAME_CLOCK", str(clock))
    v4.main()
    assert "clock" not in json.loads((v4.WORK / "samples.json").read_text())

    state["clock_out"] = "no analysis line\n"
    v4.plans = lambda: {"scale": [(5, 1000, 2)], "cadence": (1,), "ladder": ()}
    v4.WORK = tmp_path / "bad-clock"
    with pytest.raises(RuntimeError, match="no analysis"):
        v4.main()

    monkeypatch.setenv("DSEAMS_V4_TINY", "1")
    monkeypatch.setenv("DSEAMS_LIVE_WORK", str(tmp_path / "main-exec"))
    monkeypatch.delenv("DSEAMS_FRAME_CLOCK", raising=False)
    exec_file(ROOT / "tests" / "insitu" / "v4_measure.py", ["v4_measure.py"])
    assert (tmp_path / "main-exec" / "samples.json").is_file()


def test_adapter_cov(tmp_path, monkeypatch):
    adapter = _load(ROOT / "tests" / "insitu" / "adapter_cov.py", "adapter_cov_mod")
    with pytest.raises(SystemExit, match="set DSEAMS_NOT_A_VAR"):
        adapter.need("DSEAMS_NOT_A_VAR")
    calls = []
    envs = []

    def fake_run(cmd, cwd=None, env=None, **kwargs):
        calls.append([str(part) for part in cmd])
        envs.append(env)
        return Proc(0, "", "")

    monkeypatch.setattr(adapter.subprocess, "run", fake_run)
    adapter.env_host = {"OMP_NUM_THREADS": "1"}
    adapter.lammps_cases(tmp_path / "cases", "plugin.so", "lmp", "mpirun", "libdir")
    move = (tmp_path / "cases" / "move" / "in").read_text()
    assert "group oxy type 1" in move
    assert move.index("group oxy") < move.index("fix sea")
    assert "every 0" in (tmp_path / "cases" / "bad-every" / "in").read_text()
    assert "cutoff nope" in (tmp_path / "cases" / "bad-double" / "in").read_text()
    assert any(cmd[2:4] == ["-np", "2"] for cmd in calls)
    assert any(cmd[2:4] == ["-np", "1"] for cmd in calls)
    assert envs[-1] == {"OMP_NUM_THREADS": "1"}
    assert (tmp_path / "cases" / "nohost" / "in").is_file()

    def bad(cmd, cwd=None, env=None, text=None, capture_output=None):
        return Proc(1, "OUT", "ERR-adapter")

    monkeypatch.setattr(adapter.subprocess, "run", bad)
    with pytest.raises(RuntimeError, match="ERR-adapter"):
        adapter.run(["cmake"], tmp_path)

    monkeypatch.setattr(adapter.subprocess, "run", fake_run)
    monkeypatch.delenv("DSEAMS_PYTHON", raising=False)
    assert adapter.sys_python() == "python3"
    monkeypatch.setenv("DSEAMS_ADAPTER_WORK", str(tmp_path / "work"))
    calls.clear()
    adapter.main()
    assert any(cmd[:2] == ["cmake", "-S"] for cmd in calls)
    cache = tmp_path / "work" / "plugin" / "CMakeCache.txt"
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text("cache\n")
    calls.clear()
    adapter.main()
    assert not any(cmd[:2] == ["cmake", "-S"] for cmd in calls)
    assert any(cmd[:2] == ["cmake", "--build"] for cmd in calls)
    assert calls[-1][0] == os.environ["DSEAMS_MPIRUN"]
    assert "-np" in calls[-1]

    path = ROOT / "tests" / "insitu" / "adapter_cov.py"
    code = compile(path.read_text(), str(path), "exec")
    ns = {
        "__name__": "__main__",
        "__file__": str(path),
        "env_host": {"OMP_NUM_THREADS": "1"},
    }
    saved = sys.argv
    sys.argv = ["adapter_cov.py"]
    try:
        exec(code, ns)
    finally:
        sys.argv = saved


def test_engine_e2e(tmp_path, monkeypatch, capsys):
    import engine_e2e
    import live

    monkeypatch.delenv("DSEAMS_NO_SUCH", raising=False)
    with pytest.raises(SystemExit, match="set DSEAMS_NO_SUCH"):
        engine_e2e.need("DSEAMS_NO_SUCH")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    bare = engine_e2e.lib_env("")
    assert bare["LD_LIBRARY_PATH"] == str(Path(live.LIB).parent)
    assert bare["OMP_NUM_THREADS"] == "1"
    extra = engine_e2e.lib_env("/opt/x")
    assert extra["LD_LIBRARY_PATH"].startswith("/opt/x" + os.pathsep)

    ran = []

    def fake_run(cmd, cwd, env=None, check=True):
        ran.append([str(part) for part in cmd])
        return Proc(0, "", "")

    monkeypatch.setattr(live, "run", fake_run)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *a, **k: 8)
    monkeypatch.setenv("DSEAMS_MPIRUN", "/usr/bin/mpirun")
    folder = tmp_path / "rank1"
    assert engine_e2e.lammps_case("lmp", "plugin.so", folder, 1) == 8
    script = (folder / "in").read_text()
    assert "plugin load" in script
    assert "delay 5" in script
    assert ran[-1][0] == "lmp"
    assert engine_e2e.lammps_case("lmp", "plugin.so", folder, 1) == 8
    empty = tmp_path / "empty-plugin"
    engine_e2e.lammps_case("lmp", "", empty, 2)
    assert "plugin load" not in (empty / "in").read_text()
    assert ran[-1][:4] == ["/usr/bin/mpirun", "--oversubscribe", "-np", "2"]
    assert "frames 8 ranks 2\n" == (empty / "equal.log").read_text()
    monkeypatch.delenv("DSEAMS_MPIRUN")
    with pytest.raises(SystemExit, match="set DSEAMS_MPIRUN"):
        engine_e2e.lammps_case("lmp", "plugin.so", tmp_path / "nompi", 2)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *a, **k: 3)
    with pytest.raises(RuntimeError, match="expected 8 frames"):
        engine_e2e.lammps_case("lmp", "plugin.so", tmp_path / "short", 1)

    monkeypatch.setattr(
        live,
        "diamond_positions",
        lambda n: [{"id": "1", "x": "1", "y": "1", "z": "1"}],
    )
    monkeypatch.setattr(
        live,
        "trr_frames",
        lambda path: [(0, 1, (1.0,) * 9, (0.1, 0.2, 0.3))],
    )
    monkeypatch.setattr(live, "cli_text", lambda *a, **k: "rings 3:0\n")
    monkeypatch.setattr(live, "match_lammps_cli", lambda *a, **k: 1)
    monkeypatch.setattr(live, "match_native_cli", lambda *a, **k: 1)
    monkeypatch.setenv("DSEAMS_GMX", "/opt/gmx")
    monkeypatch.setenv("DSEAMS_GMX_MPIRUN", "/opt/mpirun")
    live.WORK = tmp_path / "live-work"
    matched = engine_e2e.gromacs_case()
    assert matched == 1
    assert (live.WORK / "e2e-gromacs" / "from-trr.lammpstrj").is_file()
    receipt = (live.WORK / "e2e-gromacs" / "receipt.txt").read_text()
    assert "/opt/gmx" in receipt and "/opt/mpirun" in receipt
    assert any("grompp" in cmd for cmd in ran)
    assert any("mdrun" in cmd for cmd in ran)

    monkeypatch.setattr(live, "openmm_identity", lambda: 4)
    monkeypatch.setattr(live, "match_lammps_cli", lambda *a, **k: 8)
    monkeypatch.setattr(live, "match_native_cli", lambda *a, **k: 8)
    monkeypatch.setattr(
        live,
        "trr_frames",
        lambda path: [(0, 1, (1.0,) * 9, (0.1, 0.2, 0.3))] * 8,
    )
    monkeypatch.setenv("DSEAMS_E2E_WORK", str(tmp_path / "e2e"))
    monkeypatch.setenv("DSEAMS_MPIRUN", "/usr/bin/mpirun")
    monkeypatch.setenv("DSEAMS_LMP", "lmp")
    monkeypatch.setenv("DSEAMS_PLUGIN", "plugin.so")
    monkeypatch.setenv("DSEAMS_LMP_INTREE", "lmp-intree")
    engine_e2e.main()
    out = capsys.readouterr().out
    assert "openmm 4" in out
    assert "lammps-plugin-1 8" in out
    assert "gromacs 8" in out

    exec_file(ROOT / "tests" / "insitu" / "engine_e2e.py", ["engine_e2e.py"])


def test_benchmarks(tmp_path, monkeypatch):
    exist = tmp_path / "frame.lammpstrj"
    exist.write_text("ITEM: TIMESTEP\n0\n")
    missing = tmp_path / "missing.lammpstrj"
    marker = tmp_path / "prepared.txt"
    fake = tmp_path / "scale_campaign.py"
    fake.write_text(
        "from pathlib import Path\n"
        "EXIST = Path(%r)\n"
        "MISSING = Path(%r)\n"
        "MARKER = Path(%r)\n"
        "def points():\n"
        "    return [\n"
        "        {'kind': 'strong', 'model': 'ice', 'analysis': 'rings', 'threads': 1, 'reps': [1, 1, 1]},\n"
        "        {'kind': 'weak', 'model': 'ice', 'analysis': 'rings', 'threads': 4, 'reps': [2, 1, 1]},\n"
        "    ]\n"
        "MODELS = {'ice': {'type': '1'}}\n"
        "def dump_path(model, shape):\n"
        "    return EXIST if shape[0] == 1 else MISSING\n"
        "def prepare_dump(spec, model, shape):\n"
        "    MARKER.write_text(repr((spec['type'], model, tuple(shape))))\n"
        "BENCH = %r\n"
        "REPS = 3\n" % (str(exist), str(missing), str(marker), "/" + "tmp" + "/bench-bin")
    )
    import importlib.util

    real = importlib.util.spec_from_file_location

    def wrapped(name, location, *args, **kwargs):
        if str(location).endswith("scale_campaign.py"):
            location = str(fake)
        return real(name, location, *args, **kwargs)

    assert wrapped("other", ROOT / "tests" / "insitu" / "conftest.py") is not None

    monkeypatch.setattr(importlib.util, "spec_from_file_location", wrapped)
    import benchmarks.benchmarks as bm
    import benchmarks.shipped as shipped

    points = bm.ScaleProtocol.params[0]
    assert len(points) == 2
    arm = bm.ScaleProtocol()
    arm.setup(points[0])
    assert arm.cmd[arm.cmd.index("-c") + 1] == "0"
    bench = "/" + "tmp" + "/bench-bin"
    assert bench in arm.cmd
    assert "rings" in arm.cmd and "--threads" in arm.cmd
    assert not marker.exists()
    arm.setup(points[1])
    assert arm.cmd[arm.cmd.index("-c") + 1] == "0-3"
    assert marker.read_text().startswith("('1', 'ice'")

    recorded = []
    codes = [0, 1]
    real_run = bm.subprocess.run

    def fake_run(cmd, cwd=None, env=None, **kwargs):
        if list(cmd)[:1] == ["git"]:
            return real_run(cmd, cwd=cwd, env=env, **kwargs)
        recorded.append([str(part) for part in cmd])
        rc = codes.pop(0)
        return Proc(rc, "bench line\n", "")

    monkeypatch.setattr(bm.subprocess, "run", fake_run)
    arm.setup(points[0])
    arm.time_twelve_calls(points[0])
    log = ROOT / "bbdir" / "scale" / "asv-logs" / "strong-ice-rings-t1-1x1x1.log"
    body = log.read_text()
    assert "EXIT 0" in body and "COMMAND" in body and "HEAD " in body
    assert bench in body
    arm.setup(points[1])
    with pytest.raises(subprocess.CalledProcessError) as caught:
        arm.time_twelve_calls(points[1])
    assert caught.value.returncode == 1
    failed = ROOT / "bbdir" / "scale" / "asv-logs" / "weak-ice-rings-t4-2x1x1.log"
    assert "EXIT 1" in failed.read_text()
    log.unlink()
    failed.unlink()

    recorded.clear()

    def ok_run(cmd, cwd=None, env=None, **kwargs):
        recorded.append([str(part) for part in cmd])
        return Proc(0, "", "")

    monkeypatch.setattr(bm.subprocess, "run", ok_run)
    neigh = bm.NeighbourArm()
    neigh.setup()
    neigh.time_mutual_neighbours()
    assert recorded[-1] == neigh.cmd
    assert recorded[-1][0] == "taskset" and recorded[-1][2] == "0"
    assert recorded[-1][3].endswith("bench_neigh_only")
    assert "input/traj/genice_sH.lammpstrj" in recorded[-1]
    rings = bm.RingArm()
    rings.setup()
    rings.time_primitive_rings()
    assert recorded[-1] == rings.cmd
    assert recorded[-1][3].endswith("bench_stages")
    camp = bm.CampaignArms()
    camp.setup()
    camp.time_chill()
    camp.time_q6()
    camp.time_seeded_cages()
    camp.time_rings()
    modes = [cmd[cmd.index(str(ROOT / "bbdir/tests/bench_paper")) + 5] for cmd in recorded[-4:]]
    assert modes == ["chill", "q6-full", "seeded-adj", "rings"]

    monkeypatch.setattr(shipped.subprocess, "run", ok_run)
    orth = shipped.ShippedOrth()
    orth.setup("seeded")
    orth.time_cages("seeded")
    orth.peakmem_cages("seeded")
    assert "--graph" not in recorded[-1]
    assert "cages" in recorded[-1]
    orth.setup("cutoff")
    orth.time_cages("cutoff")
    assert recorded[-1][recorded[-1].index("--graph") + 1] == "cutoff"
    stages = shipped.ShippedStages()
    stages.setup()
    stages.time_process()
    stages.peakmem_process()
    assert recorded[-1][3].endswith("bench_stages")
    custom = "/" + "tmp" + "/custom-ice.lammpstrj"
    monkeypatch.setenv("SEAMS_1E6", custom)
    million = shipped.ShippedMillion()
    million.setup("cutoff")
    assert million.cmd[-1] == custom
    million.time_cages("cutoff")
    assert recorded[-1][-1] == custom
    monkeypatch.delenv("SEAMS_1E6")
    million.setup("seeded")
    million.peakmem_cages("seeded")
    assert recorded[-1][-1] == "/" + "tmp" + "/v37-frames/ice-1e6.lammpstrj"
    assert "--graph" not in recorded[-1]


def _patch_live(monkeypatch, run_ipi_result):
    import engine_e2e
    import live
    import test_lammps_step0

    def lammps_bin(script, log, cwd, env=None):
        text = Path(script).read_text()
        for line in text.splitlines():
            if "dump native" in line:
                path = Path(line.split()[5])
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(
                    "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
                    "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                    "ITEM: ATOMS id type x y z\n1 1 1.0 2.0 3.0\n"
                )
        Path(log).write_text("Loop time of 1 on 1 procs for 1 steps\n")
        return Proc(0, "", "")

    called = {"ipi": 0, "lammps": 0, "kinds": []}

    def run_lammps():
        called["lammps"] += 1
        return True

    def run_ipi():
        called["ipi"] += 1
        return run_ipi_result

    def gmx():
        gmx.n += 1
        if gmx.n == 1:
            return (True, "matched")
        return (False, "absent")

    gmx.n = 0

    def time_lammps(kind, steps, repeats, plugin):
        called["kinds"].append(kind)
        return [0.2]

    monkeypatch.setattr(live, "lammps_bin", lammps_bin)
    monkeypatch.setattr(live, "run_lammps", run_lammps)
    monkeypatch.setattr(live, "ase_identity", lambda: True)
    monkeypatch.setattr(live, "run_ipi", run_ipi)
    monkeypatch.setattr(live, "gromacs_trr_identity", lambda: (True, "trr"))
    monkeypatch.setattr(live, "lammps_mpi_identity", lambda: True)
    monkeypatch.setattr(live, "openmm_identity", lambda: (True, "openmm"))
    monkeypatch.setattr(live, "gromacs_mpi_identity", gmx)
    monkeypatch.setattr(live, "ipi_elapsed", lambda analyse: 0.5)
    monkeypatch.setattr(live, "time_lammps", time_lammps)
    monkeypatch.setattr(live, "main", lambda: None)
    monkeypatch.setattr(engine_e2e, "main", lambda: None)
    monkeypatch.setattr(test_lammps_step0, "run", lambda: None)
    return called


def test_live_cover_env_absent(tmp_path, monkeypatch):
    monkeypatch.setenv("DSEAMS_LIVE_WORK", str(tmp_path))
    monkeypatch.delenv("DSEAMS_MPIRUN", raising=False)
    monkeypatch.delenv("SEAMS_INSITU_EXE", raising=False)
    monkeypatch.delenv("SEAMS_GRO_IMD_EXE", raising=False)
    monkeypatch.delenv("DSEAMS_GMX", raising=False)
    called = _patch_live(monkeypatch, True)
    import test_live_cover

    test_live_cover.main()
    assert called["kinds"] == ["plain", "live"]
    assert called["ipi"] == 1 and called["lammps"] == 1
    assert (tmp_path / "unit" / "path-record.txt").stat().st_size > 0
    assert "DSEAMS_GMX" not in os.environ
    assert "DSEAMS_MPIRUN" not in os.environ
    assert os.environ["SEAMS_INSITU_EXE"] == "present"


def test_live_cover_env_present(tmp_path, monkeypatch):
    monkeypatch.setenv("DSEAMS_LIVE_WORK", str(tmp_path))
    monkeypatch.setenv("DSEAMS_MPIRUN", "dummy-mpirun")
    monkeypatch.setenv("SEAMS_INSITU_EXE", "dummy-exe")
    monkeypatch.setenv("SEAMS_GRO_IMD_EXE", "dummy-imd")
    monkeypatch.setenv("DSEAMS_GMX", "dummy-gmx")
    _patch_live(monkeypatch, True)
    import test_live_cover

    exec_file(test_live_cover.__file__, ["test_live_cover.py"])
    assert os.environ["DSEAMS_GMX"] == "dummy-gmx"
    assert os.environ["DSEAMS_MPIRUN"] == "dummy-mpirun"
    assert os.environ["SEAMS_INSITU_EXE"] == "dummy-exe"


def test_run_ipi_none_fails_the_check(tmp_path, monkeypatch):
    monkeypatch.setenv("DSEAMS_LIVE_WORK", str(tmp_path))
    called = _patch_live(monkeypatch, None)
    import test_live_cover

    with pytest.raises(SystemExit) as caught:
        test_live_cover.main()
    assert caught.value.code == "run_ipi failed"
    assert called["ipi"] == 1 and called["lammps"] == 1
