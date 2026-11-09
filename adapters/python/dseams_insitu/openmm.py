"""OpenMM reporter. The state positions are copied to the host."""


class DseamsReporter:
    def __init__(self, engine, dump, interval=1, types=None):
        self.engine = engine
        self.dump = dump
        self.interval = int(interval)
        self.types = types

    def describeNextReport(self, simulation):
        step = int(simulation.currentStep)
        remain = self.interval - (step % self.interval)
        return {"steps": remain, "periodic": None, "include": ["positions"]}

    def report(self, simulation, state):
        from openmm.unit import angstrom

        pos = state.getPositions(asNumpy=True).value_in_unit(angstrom)
        vecs = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom)
        box = self.engine.box_from_rows(vecs)
        self.engine.observe(
            int(simulation.currentStep), pos, box, self.dump, types=self.types
        )
