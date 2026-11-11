# Copyright (c) d-SEAMS developers
"""OpenMM reporter. Positions come from the context, on the host."""


class DseamsReporter:
    """OpenMM reporter that runs d-SEAMS on the host every `interval` steps.

    The analysed atoms are appended to `dump` in angstrom, and the engine
    writes the neighbour, ring, CHILL, and CHILL+ record for those atoms.
    `indices` keeps a subset (TIP4P oxygen). Omit it for a one-site model.
    """

    def __init__(self, engine, dump, interval=1, types=None, indices=None):
        self.engine = engine
        self.dump = dump
        self.interval = int(interval)
        self.types = types
        self.indices = None if indices is None else list(indices)
        if dump:
            with open(dump, "w", encoding="ascii"):
                pass

    def describeNextReport(self, simulation):
        step = int(simulation.currentStep)
        remain = self.interval - (step % self.interval)
        return {"steps": remain, "periodic": None, "include": ["positions"]}

    def report(self, simulation, state):
        del state
        from openmm.unit import angstrom

        pulled = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        )
        pos = pulled.getPositions(asNumpy=True).value_in_unit(angstrom)
        vecs = pulled.getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom)
        if self.indices is not None:
            pos = pos[self.indices]
        types = self.types
        if types is not None and self.indices is not None and getattr(types, "ndim", 0):
            types = types[self.indices]
        box = self.engine.box_from_rows(vecs)
        self.engine.observe(
            int(simulation.currentStep), pos, box, self.dump, types=types
        )
