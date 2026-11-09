"""Line coverage for the Python adapters that do not need a dynamics package."""

import ctypes
import os
import socket
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "adapters", "python"))

from dseams_insitu.ase import DseamsMDObserver, is_optimizer
from dseams_insitu.engine import Engine, library_path, load, replay
from dseams_insitu.ipi import IpiClient, lj_forces, message, read_exact
from dseams_insitu.openmm import DseamsReporter


class _Opt:
    pass


class _BFGS:
    pass


class _Langevin:
    def __init__(self):
        self.nsteps = 0
        self.calls = []

    def attach(self, fn, interval=1):
        self.calls.append((fn, interval))


class _Atoms:
    def get_positions(self):
        return np.array([[0.0, 0.0, 0.0], [2.8, 0.0, 0.0]])

    class cell:
        array = np.eye(3) * 20.0


class _Sim:
    def __init__(self, step):
        self.currentStep = step


class HookTests(unittest.TestCase):
    def test_optimizer_names(self):
        self.assertTrue(is_optimizer(_BFGS()))
        self.assertTrue(is_optimizer(type("QuasiNewton", (), {})()))
        self.assertTrue(is_optimizer(type("Optimizer", (), {})()))
        self.assertFalse(is_optimizer(_Langevin()))
        dyn = _Langevin()
        dyn.nsteps = 2

        class Eng:
            def box_from_rows(self, cell):
                return [0, 0, 0, 20, 20, 20, 0, 0, 0]

            def observe(self, step, pos, box, dump):
                self.seen = (step, pos.shape, dump)

        eng = Eng()
        obs = DseamsMDObserver(eng, _Atoms(), "dump")
        obs.attach(dyn, interval=2)
        self.assertEqual(dyn.calls[0][1], 2)
        obs()
        self.assertEqual(eng.seen[0], 2)
        with self.assertRaises(TypeError):
            obs.attach(_BFGS())

    def test_message_and_lj(self):
        self.assertEqual(len(message("STATUS")), 12)
        with self.assertRaises(ValueError):
            message("this header is far too long")
        pos = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
        cell = np.eye(3) * 12.0
        energy, force = lj_forces(pos, cell)
        self.assertEqual(force.shape, pos.shape)
        self.assertTrue(np.isfinite(energy))
        a, b = socket.socketpair()
        a.sendall(b"abc")
        self.assertEqual(read_exact(b, 3), b"abc")
        a.close()
        with self.assertRaises(ConnectionError):
            read_exact(b, 1)
        b.close()

    def test_ipi_state(self):
        left, right = socket.socketpair()
        client = IpiClient("127.0.0.1", 1)
        client._status(left)
        self.assertEqual(right.recv(12), message("NEEDINIT"))
        right.sendall(np.int32(4).tobytes())
        right.sendall(np.int32(3).tobytes())
        right.sendall(b"xyz")
        client._init(left)
        self.assertTrue(client.initialised)
        client._status(left)
        self.assertEqual(right.recv(12), message("READY"))
        cell = (np.eye(3) * 20.0).reshape(-1)
        right.sendall(np.ascontiguousarray(cell, dtype=np.float64).tobytes())
        right.sendall(np.ascontiguousarray(cell, dtype=np.float64).tobytes())
        right.sendall(np.int32(2).tobytes())
        right.sendall(np.zeros(6, dtype=np.float64).tobytes())
        seen = []
        client._positions(left, lambda i, pos, box: seen.append((i, pos.shape)), 3)
        self.assertEqual(seen, [(3, (2, 3))])
        client._status(left)
        self.assertEqual(right.recv(12), message("HAVEDATA"))
        client._forces(left)
        self.assertEqual(right.recv(12), message("FORCEREADY"))
        right.settimeout(0.5)
        self.assertTrue(right.recv(256))
        with self.assertRaises(RuntimeError):
            client._forces(left)
        left.close()
        right.close()
        import dseams_insitu.ipi as ipi

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

            def sendall(self, data):
                return None

            def close(self):
                return None

        drained = Sock(b"ab")
        self.assertEqual(drained.recv(1), b"a")
        self.assertEqual(drained.recv(1), b"b")
        self.assertEqual(drained.recv(4), b"")

        real = ipi.socket.create_connection
        ipi.socket.create_connection = lambda addr: Sock(message("NOPE"))
        try:
            with self.assertRaises(RuntimeError):
                IpiClient("127.0.0.1", 9).run(lambda *args: None, max_frames=1)
            ipi.socket.create_connection = lambda addr: Sock(message("EXIT"))
            self.assertEqual(IpiClient("127.0.0.1", 9).run(lambda *args: None), 0)
        finally:
            ipi.socket.create_connection = real
        cell = np.ascontiguousarray(np.eye(3) * 20.0, dtype=np.float64).reshape(-1)
        blob = b"".join(
            [
                message("STATUS"),
                message("INIT"),
                np.int32(0).tobytes(),
                np.int32(3).tobytes(),
                b"abc",
                message("POSDATA"),
                cell.tobytes(),
                cell.tobytes(),
                np.int32(1).tobytes(),
                np.zeros(3, dtype=np.float64).tobytes(),
                message("GETFORCE"),
                message("EXIT"),
            ]
        )
        ipi.socket.create_connection = lambda addr: Sock(blob)
        try:
            seen = []
            frames = IpiClient("127.0.0.1", 9).run(lambda i, pos, box: seen.append(pos.shape))
            self.assertEqual(frames, 1)
            self.assertEqual(seen, [(1, 3)])
        finally:
            ipi.socket.create_connection = real

    def test_reporter_interval(self):
        rep = DseamsReporter(None, "dump", interval=4)
        self.assertEqual(rep.describeNextReport(_Sim(1))["steps"], 3)
        self.assertEqual(rep.describeNextReport(_Sim(4))["steps"], 4)
        self.assertEqual(rep.describeNextReport(_Sim(0))["include"], ["positions"])
        import numpy as np
        from openmm.unit import angstrom

        class Qty:
            def value_in_unit(self, unit):
                return np.eye(3) * 20.0 if unit is angstrom and self.box else np.zeros((2, 3))

            def __init__(self, box=False):
                self.box = box

        self.assertEqual(tuple(Qty().value_in_unit(None).shape), (2, 3))

        class State:
            def getPositions(self, asNumpy=True):
                return Qty()

            def getPeriodicBoxVectors(self, asNumpy=True):
                return Qty(box=True)

        class Eng:
            def box_from_rows(self, rows):
                return [0, 0, 0, 20, 20, 20, 0, 0, 0]

            def observe(self, step, pos, box, dump, types=None):
                self.step = step

        eng = Eng()
        sim = _Sim(3)

        class Ctx:
            def getState(self, getPositions=True, enforcePeriodicBox=True):
                return State()

        sim.context = Ctx()
        DseamsReporter(eng, "dump", 1).report(sim, None)
        self.assertEqual(eng.step, 3)

    def test_engine_round_trip(self):
        lib_path = os.environ["SEAMS_LIB"]
        with self.assertRaises(RuntimeError):
            os.environ.pop("SEAMS_LIB", None)
            library_path()
        os.environ["SEAMS_LIB"] = lib_path
        self.assertEqual(library_path(), lib_path)
        lib = load(lib_path)
        import dseams_insitu.engine as engmod

        class Fn:
            def __init__(self, value):
                self.value = value
                self.restype = None
                self.argtypes = None

            def __call__(self):
                return self.value

        class BadLib:
            def __init__(self, frame, options):
                self.seams_insitu_frame_bytes = Fn(frame)
                self.seams_insitu_options_bytes = Fn(options)

        real_cdll = engmod.ctypes.CDLL
        engmod.ctypes.CDLL = lambda path: BadLib(1, ctypes.sizeof(engmod.Options))
        try:
            with self.assertRaises(RuntimeError):
                load("libyodaLib.so")
            engmod.ctypes.CDLL = lambda path: BadLib(ctypes.sizeof(engmod.Frame), 1)
            with self.assertRaises(RuntimeError):
                load("libyodaLib.so")
        finally:
            engmod.ctypes.CDLL = real_cdll
        folder = tempfile.mkdtemp()
        dump = os.path.join(folder, "round.lammpstrj")
        live = os.path.join(folder, "live.txt")
        post = os.path.join(folder, "post.txt")
        with open(dump, "w", encoding="ascii") as handle:
            handle.write("stale\n")
        for path in (dump, live, post):
            if os.path.exists(path):
                os.remove(path)
        xyz = np.array(
            [[0.0, 0.0, 0.0], [2.7, 0.0, 0.0], [1.3, 2.3, 0.0], [1.3, 0.8, 2.2]],
            dtype=np.float64,
        )
        box = [0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0]
        with Engine(lib, live, debug=True) as eng:
            eng.observe(0, xyz, box, dump, types=np.array([1, 1, 1, 1], dtype=np.int32))
            with self.assertRaises(ValueError):
                eng.box_from_rows(np.array([1, 2, 0, 0, 1, 0, 0, 0, 1], dtype=np.float64))
            self.assertEqual(eng.box_from_rows(np.eye(3) * 10.0)[3], 10.0)
        replay(lib, dump, post)
        with open(live) as handle:
            left = handle.read()
        with open(post) as handle:
            right = handle.read()
        self.assertEqual(left, right)
        with self.assertRaises(RuntimeError):
            Engine(lib, "/no/such/dseams-dir/no.txt")
        with self.assertRaises(RuntimeError):
            Engine(lib, live, cutoff=0)
        with Engine(lib, live) as eng:
            with self.assertRaises(RuntimeError):
                eng.observe(1, xyz, [0, 0, 0, 0, 20, 20, 0, 0, 0], None)
            eng.observe(2, xyz, box, None, ids=np.array([1, 2, 3, 4], dtype=np.int32))
        with self.assertRaises(RuntimeError):
            replay(lib, dump, post, cutoff=0)
        with self.assertRaises(RuntimeError):
            replay(lib, os.path.join(folder, "missing.lammpstrj"), post)


if __name__ == "__main__":
    unittest.main()
