# Copyright (c) d-SEAMS developers
"""In-situ CHILL+ and primitive rings through libseams_lang."""

from __future__ import annotations

import ctypes
import os

import numpy as np

_NAMES = (
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


class InSitu:
    """One library handle. Each observe() is one frame."""

    def __init__(self, library=None, cutoff=3.5, graph="cutoff", max_ring=6):
        path = library or os.environ.get("SEAMS_LANG_LIB", "libseams_lang")
        self.lib = ctypes.CDLL(path)
        self.cutoff = float(cutoff)
        self.graph = graph
        self.max_ring = int(max_ring)
        self.frames = []
        lib = self.lib
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.restype = ctypes.c_int
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.seams_frame_rings.restype = ctypes.c_int
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_int,
        ]
        lib.seams_frame_read.restype = ctypes.c_void_p
        lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]

    def _fail(self, code, what):
        if code == 0:
            return
        err = self.lib.seams_last_error()
        text = err.decode() if err else what
        raise RuntimeError(text)

    @staticmethod
    def box_from_rows(rows):
        """OpenMM box vectors, in angstrom, as nine LAMMPS cell numbers."""
        vecs = np.asarray(rows, dtype=np.float64).reshape(3, 3)
        a, b, c = vecs
        raw = np.zeros(9, dtype=np.float64)
        raw[3] = a[0]
        raw[4] = b[1]
        raw[5] = c[2]
        raw[6] = b[0]
        raw[7] = c[0]
        raw[8] = c[1]
        return raw

    def analyse(self, positions, box):
        """CHILL+ labels and primitive-ring counts for one orthorhombic frame."""
        pos = np.ascontiguousarray(positions, dtype=np.float64)
        n = int(pos.shape[0])
        raw = np.ascontiguousarray(box, dtype=np.float64).reshape(9)
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        labels = np.empty(n, dtype=np.int32)
        counts = np.zeros(9, dtype=np.int32)
        rings = np.zeros(self.max_ring + 1, dtype=np.int32)
        lib = self.lib
        frame = lib.seams_frame_from_arrays(
            pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(n),
            raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if not frame:
            self._fail(1, "null frame")
        try:
            self._fail(
                lib.seams_frame_chill_plus(
                    ctypes.c_void_p(frame),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ),
                "CHILL+",
            )
            self._fail(
                lib.seams_frame_rings(
                    ctypes.c_void_p(frame),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    ctypes.c_int(self.max_ring),
                    rings.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_int(self.max_ring + 1),
                ),
                "rings",
            )
        finally:
            lib.seams_frame_free(ctypes.c_void_p(frame))
        return labels.copy(), counts.copy(), rings.copy()

    def analyse_file(self, path, frame_index):
        """Same classifiers, reading a dump from disk."""
        lib = self.lib
        held = lib.seams_frame_read(str(path).encode(), ctypes.c_int(frame_index), ctypes.c_int(1))
        if not held:
            self._fail(1, "unreadable frame")
        n = lib.seams_frame_nop(ctypes.c_void_p(held))
        labels = np.empty(n, dtype=np.int32)
        counts = np.zeros(9, dtype=np.int32)
        rings = np.zeros(self.max_ring + 1, dtype=np.int32)
        try:
            self._fail(
                lib.seams_frame_chill_plus(
                    ctypes.c_void_p(held),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ),
                "CHILL+",
            )
            self._fail(
                lib.seams_frame_rings(
                    ctypes.c_void_p(held),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    ctypes.c_int(self.max_ring),
                    rings.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_int(self.max_ring + 1),
                ),
                "rings",
            )
        finally:
            lib.seams_frame_free(ctypes.c_void_p(held))
        return labels.copy(), counts.copy(), rings.copy()

    def observe(self, step, positions, box, ids, types=None):
        del ids
        pos = np.asarray(positions, dtype=np.float64)
        if types is None:
            kept = pos
        else:
            kind = np.asarray(types)
            kept = pos[kind == 1]
        labels, counts, rings = self.analyse(kept, box)
        self.frames.append(
            {
                "step": int(step),
                "positions": np.array(kept, copy=True),
                "box": np.array(box, copy=True),
                "labels": labels,
                "counts": counts,
                "rings": rings,
            }
        )
        return labels, counts, rings
