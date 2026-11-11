# Copyright (c) d-SEAMS developers
"""TIP4P/Ice in OpenMM, with the in-situ reporter, checked against the file reader.

Usage:
    python3 scripts/openmm_tip4p_ice.py LIBSEAMS_LANG SEAMS LOG

4096 molecules, CPU, six repeats of a bare run and six of the reporter.
Frames written by the reporter are read back from disk. CHILL+ labels and
primitive-ring counts must match. The seams chill-plus histogram is checked
against the same counts.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "adapters" / "python"))

from dseams_insitu.engine import InSitu, _NAMES
from dseams_insitu.openmm import DseamsReporter

MOLECULES = int(os.environ.get("OPENMM_TIP4P_MOLECULES", "4096"))
STEPS = int(os.environ.get("OPENMM_TIP4P_STEPS", "12"))
INTERVAL = int(os.environ.get("OPENMM_TIP4P_INTERVAL", "6"))
REPEATS = int(os.environ.get("OPENMM_TIP4P_REPEATS", "6"))
CUTOFF = 3.5
OUT = Path(os.environ.get("OPENMM_TIP4P_DIR", "/var/tmp/dseams/openmm-tip4p"))


def water_box(n_mol):
    from openmm.app import ForceField, Modeller, PME, Topology
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picosecond

    ff = ForceField(str(ROOT / "adapters" / "python" / "dseams_insitu" / "tip4pice.xml"))
    modeller = Modeller(Topology(), [])
    modeller.addSolvent(ff, model="tip4pew", numAdded=n_mol)
    span = modeller.topology.getUnitCellDimensions().value_in_unit(nanometer)
    # PME needs the cutoff inside half the shortest edge.
    cutoff = min(0.9, 0.49 * min(span))
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff * nanometer,
        constraints=None,
        rigidWater=True,
    )
    integrator = LangevinMiddleIntegrator(260 * kelvin, 1.0 / picosecond, 0.002 * picosecond)
    platform = Platform.getPlatformByName("CPU")
    from openmm.app import Simulation

    simulation = Simulation(modeller.topology, system, integrator, platform, {"Threads": "1"})
    simulation.context.setPositions(modeller.positions)
    return simulation


def oxygen_index(topology):
    index = []
    for atom in topology.atoms():
        index.append(1 if atom.element is not None and atom.element.symbol == "O" else 0)
    return np.asarray(index, dtype=np.int32)


def num(value):
    return format(float(value), ".17g")


def write_dump(path, frames):
    with path.open("w") as out:
        for number, frame in enumerate(frames, start=1):
            pos = frame["positions"]
            box = frame["box"]
            n = pos.shape[0]
            out.write(f"ITEM: TIMESTEP\n{number}\n")
            out.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
            lx, ly, lz = box[3], box[4], box[5]
            xy, xz, yz = box[6], box[7], box[8]
            tilted = abs(xy) > 0.0 or abs(xz) > 0.0 or abs(yz) > 0.0
            if tilted:
                out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
                out.write(f"0 {num(lx)} {num(xy)}\n")
                out.write(f"0 {num(ly)} {num(xz)}\n")
                out.write(f"0 {num(lz)} {num(yz)}\n")
            else:
                out.write("ITEM: BOX BOUNDS pp pp pp\n")
                out.write(f"0 {num(lx)}\n")
                out.write(f"0 {num(ly)}\n")
                out.write(f"0 {num(lz)}\n")
            out.write("ITEM: ATOMS id type x y z\n")
            for i in range(n):
                x, y, z = pos[i]
                out.write(f"{i + 1} 1 {num(x)} {num(y)} {num(z)}\n")


def snapshot(simulation):
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    return state


def restore(simulation, state):
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    simulation.currentStep = 0


def run_steps(simulation, steps, reporter):
    simulation.reporters.clear()
    if reporter is not None:
        simulation.reporters.append(reporter)
    t0 = time.perf_counter()
    simulation.step(steps)
    return time.perf_counter() - t0


def histogram(labels):
    counts = {name: 0 for name in _NAMES}
    for value in labels:
        counts[_NAMES[int(value)]] += 1
    return counts


def cli_counts(seams, dump, frame):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["NO_COLOR"] = "1"
    proc = subprocess.run(
        [
            seams,
            "chill-plus",
            str(dump),
            "--type",
            "1",
            "--cutoff",
            "3.5",
            "--graph",
            "cutoff",
            "--frame",
            str(frame),
            "--last",
            str(frame),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    text = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout)
    parts = text.split()
    found = {}
    i = 0
    while i < len(parts):
        if parts[i] in _NAMES and i + 1 < len(parts) and parts[i + 1].isdigit():
            found[parts[i]] = int(parts[i + 1])
            i += 2
            continue
        i += 1
    return found, text


def main(argv):
    if len(argv) < 4:
        print(__doc__, file=sys.stderr)
        return 2
    library, seams, log_path = argv[1], argv[2], Path(argv[3])
    os.environ["SEAMS_LANG_LIB"] = library
    os.environ["OMP_NUM_THREADS"] = "1"
    OUT.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(text):
        print(text, flush=True)
        lines.append(text)

    say(f"molecules {MOLECULES}")
    say(f"steps {STEPS} interval {INTERVAL} repeats {REPEATS}")
    say("platform CPU threads 1")
    say(f"library {library}")
    say(f"seams {seams}")
    simulation = water_box(MOLECULES)
    say(f"atoms {simulation.system.getNumParticles()}")
    span = simulation.topology.getUnitCellDimensions()
    from openmm.unit import nanometer

    edges = span.value_in_unit(nanometer)
    say("box-nm " + " ".join(f"{float(v):.6f}" for v in edges))
    from openmm import NonbondedForce

    for force in simulation.system.getForces():
        if isinstance(force, NonbondedForce):
            say(f"pme-cutoff-nm {force.getCutoffDistance().value_in_unit(nanometer):.6f}")
    simulation.minimizeEnergy(maxIterations=40)
    state = snapshot(simulation)
    kinds = oxygen_index(simulation.topology)
    say(f"oxygens {int(kinds.sum())}")

    engine = InSitu(library)
    reporter = DseamsReporter(engine, str(OUT / "reporter.dump"), interval=INTERVAL, types=kinds)
    say("match-run")
    elapsed = run_steps(simulation, STEPS, reporter)
    say(f"match-run seconds {elapsed:.6f} reports {len(engine.frames)}")
    dump = OUT / "oxygens.lammpstrj"
    write_dump(dump, engine.frames)
    for index, frame in enumerate(engine.frames, start=1):
        labels, counts, rings = engine.analyse_file(dump, index)
        if not np.array_equal(labels, frame["labels"]):
            raise SystemExit(f"frame {index} CHILL+ labels differ")
        if not np.array_equal(counts, frame["counts"]):
            raise SystemExit(f"frame {index} CHILL+ counts differ")
        if not np.array_equal(rings, frame["rings"]):
            raise SystemExit(f"frame {index} ring counts differ")
        cli, raw = cli_counts(seams, dump, index)
        live = histogram(frame["labels"])
        for name, number in live.items():
            if number == 0:
                continue
            if cli.get(name) != number:
                raise SystemExit(f"frame {index} CLI {name} {cli.get(name)} != {number}\n{raw}")
        say(
            f"frame {index} step {frame['step']} labels-match rings "
            + " ".join(f"{size}:{int(frame['rings'][size])}" for size in range(3, 7))
        )
    say("match exact")

    bare = []
    reported = []
    for trial in range(1, REPEATS + 1):
        restore(simulation, state)
        seconds = run_steps(simulation, STEPS, None)
        bare.append(seconds)
        say(f"bare {trial} {seconds:.6f}")
        restore(simulation, state)
        engine.frames.clear()
        seconds = run_steps(simulation, STEPS, reporter)
        reported.append(seconds)
        say(f"reporter {trial} {seconds:.6f} reports {len(engine.frames)}")
    bare_mean = sum(bare) / len(bare)
    rep_mean = sum(reported) / len(reported)
    say(f"bare-mean {bare_mean:.6f}")
    say(f"reporter-mean {rep_mean:.6f}")
    say(f"overhead-mean {rep_mean - bare_mean:.6f}")
    log_path.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
