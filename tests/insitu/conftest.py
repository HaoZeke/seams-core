# Copyright (c) d-SEAMS developers
"""Give the in-situ tests a library path and a stand-in OpenMM unit."""

import os
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIB = ROOT / "bbdir-rel" / "src" / "libyodaLib.so"
SEAMS = ROOT / "bbdir-rel" / "seams"

os.environ.setdefault("SEAMS_LIB", str(LIB))
os.environ.setdefault("SEAMS_CLI", str(SEAMS))
os.environ.setdefault("DSEAMS_LMP", str(SEAMS))
os.environ.setdefault("DSEAMS_MPIRUN", "/bin/true")
os.environ.setdefault("DSEAMS_PLUGIN", str(LIB))
os.environ.setdefault("DSEAMS_GMX", str(SEAMS))
os.environ.setdefault("DSEAMS_GMX_MPIRUN", "/bin/true")
os.environ.setdefault("DSEAMS_LAMMPS_SRC", str(ROOT))
os.environ.setdefault("DSEAMS_MPICXX", "c++")
os.environ.setdefault("DSEAMS_LMP_INTREE", str(SEAMS))

def install_openmm_unit():
    try:
        import openmm.unit  # noqa: F401
    except ImportError:
        unit = types.ModuleType("openmm.unit")

        class _Angstrom:
            pass

        unit.angstrom = _Angstrom()
        package = types.ModuleType("openmm")
        package.unit = unit
        sys.modules["openmm"] = package
        sys.modules["openmm.unit"] = unit


install_openmm_unit()
