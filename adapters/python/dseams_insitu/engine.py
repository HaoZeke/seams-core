"""ctypes binding of the in-situ C core."""

import ctypes
import os


class Cadence(ctypes.Structure):
    _fields_ = [
        ("every", ctypes.c_int),
        ("first", ctypes.c_int),
        ("last", ctypes.c_int),
    ]


class Options(ctypes.Structure):
    _fields_ = [
        ("cutoff", ctypes.c_double),
        ("type", ctypes.c_int),
        ("ring_depth", ctypes.c_int),
        ("write_labels", ctypes.c_int),
        ("cadence", Cadence),
        ("debug_dump", ctypes.c_int),
    ]


class Frame(ctypes.Structure):
    _fields_ = [
        ("step", ctypes.c_int),
        ("n", ctypes.c_int),
        ("xyz", ctypes.POINTER(ctypes.c_double)),
        ("types", ctypes.POINTER(ctypes.c_int)),
        ("ids", ctypes.POINTER(ctypes.c_int)),
        ("xlo", ctypes.c_double),
        ("ylo", ctypes.c_double),
        ("zlo", ctypes.c_double),
        ("lx", ctypes.c_double),
        ("ly", ctypes.c_double),
        ("lz", ctypes.c_double),
        ("xy", ctypes.c_double),
        ("xz", ctypes.c_double),
        ("yz", ctypes.c_double),
        ("local_only", ctypes.c_int),
        ("gathered", ctypes.c_int),
        ("gpu_resident", ctypes.c_int),
        ("gpu_device", ctypes.c_int),
        ("rank", ctypes.c_int),
        ("nranks", ctypes.c_int),
    ]


def library_path(explicit=None):
    if explicit:
        return explicit
    env = os.environ.get("SEAMS_LIB")
    if env:
        return env
    raise RuntimeError("SEAMS_LIB is not set")


def load(path=None):
    lib = ctypes.CDLL(library_path(path))
    lib.seams_insitu_frame_bytes.restype = ctypes.c_int
    lib.seams_insitu_options_bytes.restype = ctypes.c_int
    if lib.seams_insitu_frame_bytes() != ctypes.sizeof(Frame):
        raise RuntimeError("frame layout does not match libyodaLib")
    if lib.seams_insitu_options_bytes() != ctypes.sizeof(Options):
        raise RuntimeError("options layout does not match libyodaLib")
    lib.seams_insitu_options_default.argtypes = [ctypes.POINTER(Options)]
    lib.seams_insitu_create.argtypes = [ctypes.POINTER(Options)]
    lib.seams_insitu_create.restype = ctypes.c_void_p
    lib.seams_insitu_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_insitu_open_record.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.seams_insitu_open_record.restype = ctypes.c_int
    lib.seams_insitu_close_record.argtypes = [ctypes.c_void_p]
    lib.seams_insitu_close_record.restype = ctypes.c_int
    lib.seams_insitu_observe.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(Frame),
        ctypes.c_char_p,
    ]
    lib.seams_insitu_observe.restype = ctypes.c_int
    lib.seams_insitu_replay.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.seams_insitu_replay.restype = ctypes.c_int
    lib.seams_insitu_box_from_rows.argtypes = [ctypes.POINTER(ctypes.c_double)] + [
        ctypes.POINTER(ctypes.c_double)
    ] * 9
    lib.seams_insitu_box_from_rows.restype = ctypes.c_int
    lib.seams_insitu_scale_nm.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
    ]
    return lib


class Engine:
    def __init__(self, lib, record, cutoff=3.5, depth=8, every=1, water_type=1, debug=False):
        self.lib = lib
        opt = Options()
        lib.seams_insitu_options_default(ctypes.byref(opt))
        opt.cutoff = float(cutoff)
        opt.ring_depth = int(depth)
        opt.type = int(water_type)
        opt.cadence.every = int(every)
        opt.debug_dump = 1 if debug else 0
        self._eng = lib.seams_insitu_create(ctypes.byref(opt))
        if not self._eng:
            raise RuntimeError("rejected in-situ options")
        self._record = record
        rc = lib.seams_insitu_open_record(self._eng, record.encode())
        if rc != 0:
            lib.seams_insitu_destroy(self._eng)
            self._eng = None
            raise RuntimeError("cannot open record " + record)

    def close(self):
        if self._eng:
            self.lib.seams_insitu_close_record(self._eng)
            self.lib.seams_insitu_destroy(self._eng)
            self._eng = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def observe(self, step, xyz, box, dump, types=None, ids=None):
        import numpy as np

        xyz = np.ascontiguousarray(xyz, dtype=np.float64)
        n = int(xyz.shape[0])
        frame = Frame()
        frame.step = int(step)
        frame.n = n
        frame.xyz = xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        type_buf = None
        id_buf = None
        if types is not None:
            type_buf = np.ascontiguousarray(types, dtype=np.int32)
            frame.types = type_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        if ids is not None:
            id_buf = np.ascontiguousarray(ids, dtype=np.int32)
            frame.ids = id_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        frame.xlo, frame.ylo, frame.zlo = (float(v) for v in box[0:3])
        frame.lx, frame.ly, frame.lz = (float(v) for v in box[3:6])
        frame.xy, frame.xz, frame.yz = (float(v) for v in box[6:9])
        frame.gathered = 1
        rc = self.lib.seams_insitu_observe(
            self._eng, ctypes.byref(frame), None if dump is None else dump.encode()
        )
        if rc != 0:
            raise RuntimeError("observe failed with status " + str(rc))
        return rc

    def box_from_rows(self, rows):
        import numpy as np

        raw = np.ascontiguousarray(rows, dtype=np.float64).reshape(9)
        buf = raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        out = [ctypes.c_double() for _ in range(9)]
        ptrs = [ctypes.byref(v) for v in out]
        rc = self.lib.seams_insitu_box_from_rows(buf, *ptrs)
        if rc != 0:
            raise ValueError("cell is not in the reduced LAMMPS frame")
        return [v.value for v in out]


def replay(lib, dump, record, cutoff=3.5, depth=8):
    opt = Options()
    lib.seams_insitu_options_default(ctypes.byref(opt))
    opt.cutoff = float(cutoff)
    opt.ring_depth = int(depth)
    eng = lib.seams_insitu_create(ctypes.byref(opt))
    if not eng:
        raise RuntimeError("rejected replay options")
    rc = lib.seams_insitu_replay(eng, dump.encode(), record.encode())
    lib.seams_insitu_destroy(eng)
    if rc != 0:
        raise RuntimeError("replay failed with status " + str(rc))
