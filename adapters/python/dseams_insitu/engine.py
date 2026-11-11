# Copyright (c) d-SEAMS developers
"""ctypes binding of libseams_lang for one in-situ frame."""

import ctypes
import os

import numpy as np

RING_DEPTH = 8
ICE_NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)

_libc = ctypes.CDLL(None)
_libc.strtod.restype = ctypes.c_double
_libc.strtod.argtypes = [ctypes.c_char_p, ctypes.c_void_p]


def c_float(text):
    """The double the LAMMPS reader stores for this token."""
    return _libc.strtod(text.encode("ascii"), None)


def format_coord(value):
    return format(float(value), ".17g")


def library_path(explicit=None):
    if explicit:
        return os.fspath(explicit)
    env = os.environ.get("SEAMS_LIB")
    if env:
        return env
    raise RuntimeError("SEAMS_LIB is not set")


def load(path=None):
    lib = ctypes.CDLL(library_path(path))
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_chill.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_chill.restype = ctypes.c_int
    lib.seams_frame_chill_plus.argtypes = lib.seams_frame_chill.argtypes
    lib.seams_frame_chill_plus.restype = ctypes.c_int
    lib.seams_frame_rings.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_neighbours.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_neighbours.restype = ctypes.c_int
    return lib


def _check(lib, rc):
    if rc == 0:
        return
    err = lib.seams_last_error()
    detail = err.decode() if err else ""
    raise RuntimeError(detail or "libseams_lang call failed")


def chill_line(counts):
    scored = int(sum(counts))
    present = sorted(
        (ICE_NAMES[i], int(counts[i])) for i in range(len(ICE_NAMES)) if counts[i]
    )
    body = " ".join(f"{name} {n}" for name, n in present)
    if body:
        return f"nop {scored} {body}"
    return f"nop {scored}"


def rings_line(counts, depth=RING_DEPTH):
    parts = [f"{size}:{int(counts[size])}" for size in range(3, depth + 1)]
    return "rings " + " ".join(parts)


def neighbour_lines(ids, degrees, flat):
    lines = ["neighbours"]
    cursor = 0
    for i, degree in enumerate(degrees):
        row = [str(int(ids[i]))]
        row.extend(str(int(flat[cursor + k])) for k in range(int(degree)))
        cursor += int(degree)
        lines.append(" ".join(row))
    lines.append("end_neighbours")
    return lines


class Engine:
    """Live census over coordinates borrowed from the caller.

    `lib` is libseams_lang. `record` receives one block per analysed step.
    Positions are angstrom. The written dump, when a path is passed to
    observe, is the frame `seams` reads back.
    """

    def __init__(self, lib, record, cutoff=3.5, depth=RING_DEPTH, every=1,
                 water_type=1):
        if isinstance(lib, (str, os.PathLike)):
            lib = load(os.fspath(lib))
        self.lib = lib
        self.cutoff = float(cutoff)
        self.depth = int(depth)
        self.every = int(every)
        self.water_type = int(water_type)
        self._record_path = record
        self._record = open(record, "w", encoding="ascii")

    def close(self):
        if self._record is not None:
            self._record.close()
            self._record = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def box_from_rows(self, rows):
        """Reduced LAMMPS edge frame from row-major lattice vectors."""
        raw = np.ascontiguousarray(rows, dtype=np.float64).reshape(9)
        ax, ay, az, bx, by, bz, cx, cy, cz = (float(v) for v in raw)
        scale = max(abs(ax), abs(by), abs(cz))
        tol = 1e-8 * (scale if scale > 1.0 else 1.0)
        if not (ax > 0.0 and by > 0.0 and cz > 0.0):
            raise ValueError("cell is not in the reduced LAMMPS frame")
        if abs(ay) > tol or abs(az) > tol or abs(bz) > tol:
            raise ValueError("cell is not in the reduced LAMMPS frame")
        return [0.0, 0.0, 0.0, ax, by, cz, bx, cx, cy]

    def observe(self, step, xyz, box, dump, types=None, ids=None):
        """Analyse one frame and append it to the record.

        `box` is xlo, ylo, zlo, lx, ly, lz, xy, xz, yz in angstrom.
        `dump`, when set, receives the same atoms in LAMMPS dump text.
        """
        step = int(step)
        if self.every > 1 and (step % self.every) != 0:
            return 0
        xyz = np.ascontiguousarray(xyz, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("positions are n by 3")
        n = int(xyz.shape[0])
        if ids is None:
            id_values = np.arange(1, n + 1, dtype=np.int32)
        else:
            id_values = np.ascontiguousarray(ids, dtype=np.int32)
        if types is None:
            type_values = np.full(n, self.water_type, dtype=np.int32)
        elif np.ndim(types) == 0:
            type_values = np.full(n, int(types), dtype=np.int32)
        else:
            type_values = np.ascontiguousarray(types, dtype=np.int32)
        box = [float(v) for v in box]
        coord_text = [format_coord(v) for v in xyz.ravel()]
        box_text = [format_coord(v) for v in box]
        xyz_rt = np.array([c_float(tok) for tok in coord_text], dtype=np.float64)
        box_rt = [c_float(tok) for tok in box_text]
        if dump:
            _append_dump(dump, step, id_values, type_values, coord_text, box_text)
        text = self._analyse(step, xyz_rt.reshape(n, 3), box_rt, type_values, id_values)
        self._record.write(text)
        self._record.flush()
        return 0

    def _analyse(self, step, xyz, box, types, ids):
        n = int(xyz.shape[0])
        if n == 0:
            rings = rings_line([0] * (self.depth + 1), self.depth)
            lines = [
                f"step {step}",
                "neighbours",
                "end_neighbours",
                rings,
                "chill nop 0",
                "chill-plus nop 0",
                "",
            ]
            return "\n".join(lines)
        xyz_c = np.ascontiguousarray(xyz, dtype=np.float64)
        box_c = np.ascontiguousarray(box, dtype=np.float64)
        frame = self.lib.seams_frame_from_arrays(
            xyz_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            box_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if not frame:
            err = self.lib.seams_last_error()
            detail = err.decode() if err else "rejected frame"
            raise RuntimeError(detail)
        try:
            graph = b"cutoff"
            degrees = (ctypes.c_int * n)()
            cap = max(n * 32, 1)
            flat = (ctypes.c_int * cap)()
            n_flat = ctypes.c_int()
            rc = self.lib.seams_frame_neighbours(
                frame, self.cutoff, self.water_type, graph, degrees, flat, cap,
                ctypes.byref(n_flat),
            )
            if rc != 0 and n_flat.value > cap:
                cap = int(n_flat.value)
                flat = (ctypes.c_int * cap)()
                rc = self.lib.seams_frame_neighbours(
                    frame, self.cutoff, self.water_type, graph, degrees, flat,
                    cap, ctypes.byref(n_flat),
                )
            _check(self.lib, rc)
            ring_counts = (ctypes.c_int * (self.depth + 1))()
            _check(
                self.lib,
                self.lib.seams_frame_rings(
                    frame, self.cutoff, self.water_type, graph, self.depth,
                    ring_counts, self.depth + 1,
                ),
            )
            chill_counts = (ctypes.c_int * 9)()
            _check(
                self.lib,
                self.lib.seams_frame_chill(
                    frame, self.cutoff, self.water_type, graph, None, chill_counts
                ),
            )
            plus_counts = (ctypes.c_int * 9)()
            _check(
                self.lib,
                self.lib.seams_frame_chill_plus(
                    frame, self.cutoff, self.water_type, graph, None, plus_counts
                ),
            )
        finally:
            self.lib.seams_frame_free(frame)
        lines = [f"step {step}"]
        lines.extend(neighbour_lines(ids, degrees, flat[: int(n_flat.value)]))
        lines.append(rings_line(ring_counts, self.depth))
        lines.append("chill " + chill_line(chill_counts))
        lines.append("chill-plus " + chill_line(plus_counts))
        lines.append("")
        return "\n".join(lines)


def _append_dump(path, step, ids, types, coord_text, box_text):
    n = len(ids)
    xlo, ylo, zlo, lx, ly, lz, xy, xz, yz = box_text
    tilted = c_float(xy) != 0.0 or c_float(xz) != 0.0 or c_float(yz) != 0.0
    with open(path, "a", encoding="ascii") as out:
        out.write(f"ITEM: TIMESTEP\n{step}\n")
        out.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
        if tilted:
            out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            out.write(f"{xlo} {format_coord(c_float(xlo) + c_float(lx))} {xy}\n")
            out.write(f"{ylo} {format_coord(c_float(ylo) + c_float(ly))} {xz}\n")
            out.write(f"{zlo} {format_coord(c_float(zlo) + c_float(lz))} {yz}\n")
        else:
            out.write("ITEM: BOX BOUNDS pp pp pp\n")
            out.write(f"{xlo} {format_coord(c_float(xlo) + c_float(lx))}\n")
            out.write(f"{ylo} {format_coord(c_float(ylo) + c_float(ly))}\n")
            out.write(f"{zlo} {format_coord(c_float(zlo) + c_float(lz))}\n")
        out.write("ITEM: ATOMS id type x y z\n")
        for i in range(n):
            x, y, z = coord_text[3 * i : 3 * i + 3]
            out.write(f"{int(ids[i])} {int(types[i])} {x} {y} {z}\n")
