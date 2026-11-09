# Copyright (c) d-SEAMS developers
"""Molecular dynamics observer. Optimizers are refused."""

_OPTIMIZERS = {
    "bfgs",
    "lbfgs",
    "fire",
    "gpmin",
    "mdmin",
    "quasinewton",
    "goodoldquasinewton",
    "ode",
}


def is_optimizer(dyn):
    name = type(dyn).__name__.lower().replace("_", "")
    if "optim" in name:
        return True
    return name in _OPTIMIZERS


class DseamsMDObserver:
    def __init__(self, engine, atoms, dump):
        self.engine = engine
        self.atoms = atoms
        self.dump = dump
        self.dyn = None

    def attach(self, dyn, interval=1):
        if is_optimizer(dyn):
            raise TypeError("d-SEAMS attaches to a molecular dynamics driver")
        self.dyn = dyn
        dyn.attach(self, interval=interval)

    def __call__(self):
        import numpy as np

        pos = np.ascontiguousarray(self.atoms.get_positions(), dtype=np.float64)
        cell = np.ascontiguousarray(self.atoms.cell.array, dtype=np.float64)
        box = self.engine.box_from_rows(cell)
        numbers = np.ascontiguousarray(self.atoms.get_atomic_numbers(), dtype=np.int32)
        self.engine.observe(int(self.dyn.nsteps), pos, box, self.dump, types=numbers)
