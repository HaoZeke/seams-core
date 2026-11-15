# Copyright (c) d-SEAMS developers
"""Rings, chordless cycles, cages, signature cages, and the shell symmetry label.

Set ``SEAMS_LANG_LIB`` to the built ``libseams_lang``.
"""

from __future__ import annotations

import ctypes
import os

import numpy as np


def _library():
    path = os.environ.get("SEAMS_LANG_LIB", "libseams_lang")
    lib = ctypes.CDLL(path)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_read.restype = ctypes.c_void_p
    lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_chordless.restype = ctypes.c_int
    lib.seams_frame_chordless.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_cages.restype = ctypes.c_int
    lib.seams_frame_signature.restype = ctypes.c_int
    lib.seams_frame_symmetry.restype = ctypes.c_int
    return lib


def _box9(box):
    raw = np.zeros(9, dtype=np.float64)
    lengths = np.asarray(box, dtype=np.float64).reshape(-1)
    if lengths.size == 3:
        raw[3:6] = lengths
    elif lengths.size == 9:
        raw[:] = lengths
    else:
        raise ValueError("box must be 3 lengths or 9 LAMMPS cell numbers")
    return raw


def _frame(lib, positions, box):
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    n = pos.shape[0]
    raw = _box9(box)
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    frame = lib.seams_frame_from_arrays(
        pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if not frame:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "null frame")
    return frame, n


def _read(lib, path, frame, type_filter):
    handle = lib.seams_frame_read(os.fsencode(path), frame, type_filter)
    if not handle:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else f"cannot read {path}")
    return handle, lib.seams_frame_nop(ctypes.c_void_p(handle))


def _check(lib, code):
    if code != 0:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "d-SEAMS call failed")


def rings(positions, box, cutoff=3.5, max_depth=6, graph="cutoff"):
    """Primitive-ring histogram. Index ``i`` counts rings of size ``i``."""
    lib = _library()
    frame, _n = _frame(lib, positions, box)
    counts = np.zeros(max_depth + 1, dtype=np.int32)
    try:
        code = lib.seams_frame_rings(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            ctypes.c_char_p(graph.encode()),
            ctypes.c_int(max_depth),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(counts.size),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    _check(lib, code)
    return counts


def rings_chordless(source, box=None, size=6, cutoff=3.5, graph="cutoff", type=0,
                    frame=1, type_filter=0):
    """Chordless cycles of one size, 3 through 12.

    ``source`` is a trajectory path, read at ``frame`` with ``type_filter``
    (0 guesses the oxygen species, -1 keeps every atom), or positions with
    ``box``. Returns an ``(n, size)`` array of 0-based atom indices and the
    number of cycles through each atom. Each row starts at its lowest index
    and runs towards the smaller of its two neighbours; the rows are sorted,
    so they equal ``seams_rings_chordless`` on the same frame.
    """
    lib = _library()
    if isinstance(source, (str, bytes, os.PathLike)):
        handle, n = _read(lib, source, frame, type_filter)
    elif box is None:
        raise ValueError("positions need a box")
    else:
        handle, n = _frame(lib, source, box)
    per = np.zeros(max(n, 1), dtype=np.int32)
    count = ctypes.c_int()
    try:
        code = lib.seams_frame_chordless(
            ctypes.c_void_p(handle), cutoff, type, graph.encode(), size, None, 0,
            ctypes.byref(count), per.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        rows = np.zeros((count.value if code == 0 else 0, size), dtype=np.int32)
        if code == 0 and rows.size:
            code = lib.seams_frame_chordless(
                ctypes.c_void_p(handle), cutoff, type, graph.encode(), size,
                rows.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), rows.size,
                ctypes.byref(count), None)
    finally:
        lib.seams_frame_free(ctypes.c_void_p(handle))
    _check(lib, code)
    return rows, per[:n]


def cages(positions, box, cutoff=3.5, graph="seeded", k=4):
    """Per-atom cage code and the four population counts."""
    lib = _library()
    frame, n = _frame(lib, positions, box)
    per = np.empty(n, dtype=np.int32)
    n_hex = ctypes.c_int()
    n_cubic = ctypes.c_int()
    n_water = ctypes.c_int()
    n_both = ctypes.c_int()
    try:
        code = lib.seams_frame_cages(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            ctypes.c_int(k),
            ctypes.c_char_p(graph.encode()),
            ctypes.c_int(0),
            per.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(n_hex),
            ctypes.byref(n_cubic),
            ctypes.byref(n_water),
            ctypes.byref(n_both),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    _check(lib, code)
    return per, (n_hex.value, n_cubic.value, n_water.value, n_both.value)


def signature_cages(positions, box, spec="512", cutoff=3.5, graph="cutoff"):
    """Closed cages whose ring census matches ``spec``. Returns (n_cages, n_atoms)."""
    lib = _library()
    frame, _n = _frame(lib, positions, box)
    n_cages = ctypes.c_int()
    n_atoms = ctypes.c_int()
    try:
        code = lib.seams_frame_signature(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            ctypes.c_char_p(graph.encode()),
            ctypes.c_char_p(spec.encode()),
            ctypes.byref(n_cages),
            ctypes.byref(n_atoms),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    _check(lib, code)
    return n_cages.value, n_atoms.value


def symmetry(positions, box, shell_rmsd=0.12):
    """Per-atom shell label. 0 ice Ih, 1 ice Ic, 2 stacking-disordered, -1 otherwise."""
    lib = _library()
    frame, n = _frame(lib, positions, box)
    labels = np.empty(n, dtype=np.int32)
    try:
        code = lib.seams_frame_symmetry(
            ctypes.c_void_p(frame),
            ctypes.c_double(shell_rmsd),
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    _check(lib, code)
    return labels
