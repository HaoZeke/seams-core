"""Live water and ice runs. The in-situ record must match a replay of the dump."""

import json
import os
import socket
import subprocess
import sys
import textwrap
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "adapters" / "python"))

PREFIX = Path(os.environ.get("DSEAMS_MD_PREFIX", "var/run/dseams-md/.pixi/envs/default"))
LIB = os.environ["SEAMS_LIB"]
REPLAY = os.environ["SEAMS_INSITU_EXE"]
IMD = os.environ["SEAMS_GRO_IMD_EXE"]
PLUGIN = Path(os.environ.get("SEAMS_LAMMPS_PLUGIN", "var/run/libfixdseams.so"))
WORK = Path(os.environ.get("DSEAMS_LIVE_WORK", "var/run/dseams-live"))
os.environ["SEAMS_LIB"] = LIB
os.environ["PATH"] = str(PREFIX / "bin") + os.pathsep + os.environ.get("PATH", "")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["PLUMED_KERNEL"] = str(PREFIX / "lib" / "libplumedKernel.so")

from dseams_insitu.ase import DseamsMDObserver
from dseams_insitu.engine import Engine, load, replay
from dseams_insitu.ipi import IpiClient
from dseams_insitu.openmm import DseamsReporter


def run(cmd, cwd, env=None, check=True):
    merged = os.environ.copy()
    if env:
        merged.update(env)
    proc = subprocess.run(
        cmd, cwd=cwd, env=merged, text=True, capture_output=True, input="\n\n\n"
    )
    if check and proc.returncode != 0:
        raise RuntimeError(
            "command failed: "
            + " ".join(cmd)
            + "\n"
            + proc.stdout
            + "\n"
            + proc.stderr
        )
    return proc


def same(live, post):
    left = Path(live).read_text()
    right = Path(post).read_text()
    if left != right:
        raise RuntimeError("live record and replay differ\n" + left[:400] + "\n---\n" + right[:400])
    if "step " not in left:
        raise RuntimeError("empty census")
    return left.count("\nstep ") + (1 if left.startswith("step ") else 0)


def host_block():
    cpu = "unknown"
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.lower().startswith("model name"):
            cpu = line.split(":", 1)[1].strip()
            break
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    tree = subprocess.check_output(["git", "rev-parse", "HEAD^{tree}"], cwd=ROOT, text=True).strip()
    return {
        "hostname": socket.gethostname(),
        "cpu": cpu,
        "head": head,
        "tree": tree,
    }


def sha256(path):
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lammps_input(path, steps, dump, record, plugin):
    path.write_text(
        textwrap.dedent(
            f"""\
            units real
            atom_style atomic
            boundary p p p
            lattice diamond 6.4
            region box block 0 1 0 1 0 1
            create_box 1 box
            create_atoms 1 box
            mass 1 18.0154
            pair_style lj/cut 8.0
            pair_coeff * * 0.05 2.5
            neighbor 2.0 bin
            timestep 2.0
            velocity all create 30.0 12345
            fix int all nve
            plugin load {plugin}
            fix sea all dseams dump {dump} record {record} every 1 cutoff 3.5 depth 8
            run {steps}
            """
        )
    )


def lammps_plain(path, steps):
    path.write_text(
        textwrap.dedent(
            f"""\
            units real
            atom_style atomic
            boundary p p p
            lattice diamond 6.4
            region box block 0 1 0 1 0 1
            create_box 1 box
            create_atoms 1 box
            mass 1 18.0154
            pair_style lj/cut 8.0
            pair_coeff * * 0.05 2.5
            neighbor 2.0 bin
            timestep 2.0
            velocity all create 30.0 12345
            fix int all nve
            run {steps}
            """
        )
    )


def loop_time(log_text):
    for line in log_text.splitlines():
        if line.startswith("Loop time of "):
            return float(line.split()[3])
    raise RuntimeError("LAMMPS did not report a loop time\n" + log_text[-500:])


def time_lammps(kind, steps, repeats, plugin):
    samples = []
    folder = WORK / "lammps-time" / kind
    folder.mkdir(parents=True, exist_ok=True)
    for i in range(repeats):
        script = folder / f"in.{i}"
        log = folder / f"log.{i}"
        if kind == "live":
            lammps_input(script, steps, folder / f"d.{i}", folder / f"r.{i}", plugin)
        else:
            lammps_plain(script, steps)
        proc = run([str(PREFIX / "bin" / "lmp"), "-in", str(script), "-log", str(log)], folder)
        samples.append(loop_time(proc.stdout + Path(log).read_text()) / steps)
    return samples


def run_lammps():
    folder = WORK / "lammps"
    folder.mkdir(parents=True, exist_ok=True)
    dump = folder / "live.lammpstrj"
    record = folder / "live.txt"
    post = folder / "post.txt"
    dump.unlink(missing_ok=True)
    record.unlink(missing_ok=True)
    script = folder / "in.live"
    lammps_input(script, 6, dump, record, PLUGIN)
    run([str(PREFIX / "bin" / "lmp"), "-in", str(script), "-log", str(folder / "log")], folder)
    run([REPLAY, "replay", str(dump), str(post), "3.5"], folder)
    frames = same(record, post)
    # rejected fix arguments must leave the run alive
    bad = folder / "in.bad"
    bad.write_text(
        "units real\natom_style atomic\nlattice diamond 6.4\nregion box block 0 1 0 1 0 1\n"
        "create_box 1 box\ncreate_atoms 1 box\nmass 1 1\npair_style lj/cut 8\n"
        "pair_coeff * * 0.1 2.5\nplugin load " + str(PLUGIN) + "\n"
        "fix sea all dseams dump\nrun 1\n"
    )
    run([str(PREFIX / "bin" / "lmp"), "-in", str(bad), "-log", str(folder / "bad.log")], folder)
    return {"frames": frames, "steps": 40}


def gromacs_box(gro):
    lines = Path(gro).read_text().splitlines()
    nums = [float(x) for x in lines[-1].split()]
    if len(nums) < 3:
        raise RuntimeError("gro box")
    # nanometre edge lengths, origin zero, no tilt in this setup
    return [0, 0, 0, nums[0] * 10.0, nums[1] * 10.0, nums[2] * 10.0, 0, 0, 0]


def prepare_gromacs(folder):
    folder.mkdir(parents=True, exist_ok=True)
    pdb = folder / "water.pdb"
    pdb.write_text(
        "ATOM      1  O   HOH     1       1.000   1.000   1.000  1.00  0.00\n"
        "ATOM      2  H1  HOH     1       1.957   1.000   1.000  1.00  0.00\n"
        "ATOM      3  H2  HOH     1       1.240   1.927   1.000  1.00  0.00\n"
        "END\n"
    )
    gmx = str(PREFIX / "bin" / "gmx")
    run([gmx, "pdb2gmx", "-f", "water.pdb", "-o", "conf.gro", "-p", "topol.top", "-water", "tip3p", "-ff", "amber99sb-ildn", "-ignh"], folder)
    run([gmx, "editconf", "-f", "conf.gro", "-o", "box.gro", "-box", "2.8", "2.8", "2.8", "-c"], folder)
    run([gmx, "solvate", "-cp", "box.gro", "-cs", str(PREFIX / "share" / "gromacs" / "top" / "spc216.gro"), "-o", "solv.gro", "-p", "topol.top", "-maxsol", "6"], folder)
    (folder / "md.mdp").write_text(
        "integrator = md\n"
        "dt = 0.001\n"
        "nsteps = NSTEPS\n"
        "nstlist = 10\n"
        "cutoff-scheme = Verlet\n"
        "nstcalcenergy = 1\n"
        "nstenergy = 1\n"
        "coulombtype = Cut-off\n"
        "rcoulomb = 0.9\n"
        "rvdw = 0.9\n"
        "pbc = xyz\n"
        "constraints = h-bonds\n"
        "constraint-algorithm = LINCS\n"
        "IMD-group = System\n"
        "gen-vel = yes\n"
        "gen-temp = 300\n"
        "gen-seed = 12345\n"
    )


def gromacs_once(folder, steps, port):
    text = (folder / "md.mdp").read_text().replace("NSTEPS", str(steps))
    mdp = folder / f"run-{steps}-{port}.mdp"
    mdp.write_text(text)
    gmx = str(PREFIX / "bin" / "gmx")
    tpr = folder / f"topol-{port}.tpr"
    run([gmx, "grompp", "-f", str(mdp), "-c", "solv.gro", "-p", "topol.top", "-o", str(tpr), "-maxwarn", "2"], folder)
    dump = folder / f"live-{port}.lammpstrj"
    record = folder / f"live-{port}.txt"
    log = folder / f"md-{port}.log"
    proc = subprocess.Popen(
        [gmx, "mdrun", "-s", str(tpr), "-deffnm", f"run-{port}", "-imdwait", "-imdport", str(port), "-ntomp", "1", "-noconfout"],
        cwd=folder,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    time.sleep(1.5)
    box = gromacs_box(folder / "solv.gro")
    started = time.perf_counter()
    client = subprocess.run(
        [IMD, "127.0.0.1", str(port), str(dump), str(record), *[f"{v:.17g}" for v in box], "0"],
        text=True,
        capture_output=True,
    )
    elapsed = time.perf_counter() - started
    out, _ = proc.communicate(timeout=120)
    (folder / f"mdrun-{port}.out").write_text(out)
    if proc.returncode != 0 or client.returncode != 0:
        raise RuntimeError(out + "\n" + client.stdout + "\n" + client.stderr)
    return dump, record, elapsed


def free_port():
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port


def run_gromacs():
    folder = WORK / "gromacs"
    prepare_gromacs(folder)
    port = free_port()
    dump, record, _ = gromacs_once(folder, 6, port)
    post = folder / "post.txt"
    run([REPLAY, "replay", str(dump), str(post), "3.5"], folder)
    frames = same(record, post)
    return {"frames": frames, "steps": 20}


def time_gromacs(repeats, steps):
    folder = WORK / "gromacs"
    plain = []
    live = []
    gmx = str(PREFIX / "bin" / "gmx")
    for i in range(repeats):
        text = (folder / "md.mdp").read_text().replace("NSTEPS", str(steps))
        mdp = folder / f"plain-{i}.mdp"
        mdp.write_text(text.replace("IMD-group = System\n", ""))
        tpr = folder / f"plain-{i}.tpr"
        run([gmx, "grompp", "-f", str(mdp), "-c", "solv.gro", "-p", "topol.top", "-o", str(tpr), "-maxwarn", "2"], folder)
        started = time.perf_counter()
        run([gmx, "mdrun", "-s", str(tpr), "-deffnm", f"plain-{i}", "-ntomp", "1", "-noconfout"], folder)
        plain.append((time.perf_counter() - started) / steps)
        port = free_port()
        _, _, elapsed = gromacs_once(folder, steps, port)
        live.append(elapsed / steps)
    return plain, live


def ice_positions():
    import numpy as np

    # conventional cubic ice oxygen sites, 8 atoms, lattice 6.4 angstrom
    base = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.25, 0.25, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
            [0.75, 0.75, 0.25],
        ]
    )
    return base * 6.4


def run_openmm(steps):
    from openmm.app import ForceField, Modeller, PDBFile, Simulation
    from openmm import LangevinIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picosecond

    folder = WORK / "openmm"
    folder.mkdir(parents=True, exist_ok=True)
    pdb_path = folder / "water.pdb"
    lines = ["CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1           1"]
    serial = 1
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                x = 4.0 + 8.0 * ix
                y = 4.0 + 8.0 * iy
                z = 4.0 + 8.0 * iz
                res = serial
                lines.append(
                    "ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           O" % (serial, res, x, y, z)
                )
                serial += 1
                lines.append(
                    "ATOM  %5d  H1  HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           H" % (serial, res, x + 0.757, y, z)
                )
                serial += 1
                lines.append(
                    "ATOM  %5d  H2  HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           H" % (serial, res, x - 0.239, y + 0.927, z)
                )
                serial += 1
    lines.append("END")
    pdb_path.write_text("\n".join(lines) + "\n")
    pdb = PDBFile(str(pdb_path))
    ff = ForceField("tip3p.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    from openmm.app import CutoffPeriodic, HBonds

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=0.7 * nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picosecond)
    platform = Platform.getPlatformByName("CPU")
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    types = []
    for atom in modeller.topology.atoms():
        types.append(1 if atom.element.symbol == "O" else 2)
    return folder, sim, types


def openmm_identity():
    folder, sim, types = run_openmm(6)
    lib = load(LIB)
    dump = str(folder / "live.lammpstrj")
    record = str(folder / "live.txt")
    post = str(folder / "post.txt")
    with Engine(lib, record) as eng:
        sim.reporters.append(DseamsReporter(eng, dump, 1, types=types))
        sim.step(6)
    replay(lib, dump, post)
    return same(record, post)


def openmm_time(steps):
    folder, sim, _types = run_openmm(steps)
    started = time.perf_counter()
    sim.step(steps)
    plain = (time.perf_counter() - started) / steps
    folder, sim, types = run_openmm(steps)
    lib = load(LIB)
    with Engine(lib, str(folder / "time.txt")) as eng:
        sim.reporters.append(DseamsReporter(eng, str(folder / "time.lammpstrj"), 1, types=types))
        started = time.perf_counter()
        sim.step(steps)
        live = (time.perf_counter() - started) / steps
    return plain, live


def ase_atoms():
    import numpy as np
    from ase import Atoms
    from ase.calculators.lj import LennardJones

    pos = ice_positions()
    atoms = Atoms("O8", positions=pos, cell=[6.4, 6.4, 6.4], pbc=True)
    atoms.calc = LennardJones(epsilon=0.01, sigma=2.5, rc=6.0)
    return atoms


def ase_identity():
    from ase.md.langevin import Langevin
    from ase import units

    folder = WORK / "ase"
    folder.mkdir(parents=True, exist_ok=True)
    atoms = ase_atoms()
    lib = load(LIB)
    dump = str(folder / "live.lammpstrj")
    record = str(folder / "live.txt")
    with Engine(lib, record) as eng:
        dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=20.0, friction=0.01 / units.fs, fixcm=False)
        DseamsMDObserver(eng, atoms, dump).attach(dyn, 1)
        dyn.run(6)
    post = str(folder / "post.txt")
    replay(lib, dump, post)
    return same(record, post)


def ase_time(steps, live):
    from ase.md.langevin import Langevin
    from ase import units

    atoms = ase_atoms()
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=20.0, friction=0.01 / units.fs, fixcm=False)
    if live:
        lib = load(LIB)
        folder = WORK / "ase"
        with Engine(lib, str(folder / "time.txt")) as eng:
            DseamsMDObserver(eng, atoms, str(folder / "time.lammpstrj")).attach(dyn, 1)
            started = time.perf_counter()
            dyn.run(steps)
            return (time.perf_counter() - started) / steps
    started = time.perf_counter()
    dyn.run(steps)
    return (time.perf_counter() - started) / steps


def write_ipi(folder):
    import numpy as np

    pos = ice_positions()
    xyz = folder / "ice.xyz"
    lines = [str(len(pos)), "# CELL(abcABC): 6.4 6.4 6.4 90.0 90.0 90.0"]
    for row in pos:
        lines.append(f"O {row[0]:.8f} {row[1]:.8f} {row[2]:.8f}")
    xyz.write_text("\n".join(lines) + "\n")
    port = free_port()
    xml = folder / "input.xml"
    xml.write_text(
        f"""<simulation verbosity='low'>
  <output prefix='out'>
    <properties stride='100' filename='md'> [ step ] </properties>
  </output>
  <total_steps>6</total_steps>
  <prng><seed>12345</seed></prng>
  <ffsocket name='dseams' mode='inet'>
    <address>127.0.0.1</address>
    <port>{port}</port>
    <timeout>30</timeout>
    <consolidate_messages>False</consolidate_messages>
  </ffsocket>
  <system>
    <initialize nbeads='1'>
      <file mode='xyz' units='angstrom'>ice.xyz</file>
      <velocities mode='thermal' units='kelvin'>10</velocities>
    </initialize>
    <ensemble><temperature units='kelvin'>10</temperature></ensemble>
    <motion mode='dynamics'>
      <dynamics mode='nve'><timestep units='femtosecond'>0.5</timestep></dynamics>
    </motion>
    <forces><force forcefield='dseams'></force></forces>
  </system>
</simulation>
"""
    )
    return port, xml


def run_ipi():
    import numpy as np

    folder = WORK / "ipi"
    folder.mkdir(parents=True, exist_ok=True)
    port, xml = write_ipi(folder)
    proc = subprocess.Popen(
        [str(PREFIX / "bin" / "i-pi"), str(xml)],
        cwd=folder,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    time.sleep(1.5)
    if proc.poll() is not None:
        raise RuntimeError(proc.stdout.read())
    lib = load(LIB)
    dump = str(folder / "live.lammpstrj")
    record = str(folder / "live.txt")
    seen = {"n": 0}

    def on_frame(index, pos, cell):
        if "eng" not in seen:
            return
        box = seen["eng"].box_from_rows(cell)
        seen["eng"].observe(index, pos, box, dump)
        seen["n"] += 1

    with Engine(lib, record) as eng:
        seen["eng"] = eng
        IpiClient("127.0.0.1", port).run(on_frame)
    out, _ = proc.communicate(timeout=60)
    (folder / "ipi.out").write_text(out)
    if proc.returncode != 0:
        raise RuntimeError(out)
    post = str(folder / "post.txt")
    replay(lib, dump, post)
    return same(record, post)


def compile_plumed():
    so = WORK / "libDseams.so"
    WORK.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PREFIX / "bin" / "x86_64-conda-linux-gnu-g++"),
        "-shared",
        "-fPIC",
        "-std=c++17",
        "-O2",
        "-I" + str(PREFIX / "include" / "plumed"),
        "-I" + str(ROOT / "src" / "include" / "internal"),
        str(ROOT / "src" / "adapters" / "plumed" / "Dseams.cpp"),
        "-o",
        str(so),
        "-L" + str(Path(LIB).parent),
        "-lyodaLib",
        "-L" + str(PREFIX / "lib"),
        "-lplumed",
        "-Wl,-rpath," + str(Path(LIB).parent),
        "-Wl,-rpath," + str(Path(LIB).parent.parent / "subprojects" / "readcon-core"),
        "-Wl,-rpath," + str(Path(LIB).parent.parent / "subprojects" / "readcon-db"),
        "-Wl,-rpath," + str(PREFIX / "lib"),
    ]
    run(cmd, WORK)
    return so


def run_plumed(so):
    import ctypes
    import numpy as np

    libp = ctypes.CDLL(str(PREFIX / "lib" / "libplumed.so"))

    class P(ctypes.Structure):
        _fields_ = [("p", ctypes.c_void_p)]

    libp.plumed_create.restype = P
    libp.plumed_cmd.argtypes = [P, ctypes.c_char_p, ctypes.c_void_p]
    libp.plumed_finalize.argtypes = [P]
    folder = WORK / "plumed"
    folder.mkdir(parents=True, exist_ok=True)
    pos = ice_positions()
    n = len(pos)
    p = libp.plumed_create()
    if not p.p:
        raise RuntimeError("plumed_create failed")
    nat = ctypes.c_int(n)
    dt = ctypes.c_double(0.001)
    libp.plumed_cmd(p, b"setNatoms", ctypes.byref(nat))
    libp.plumed_cmd(p, b"setTimestep", ctypes.byref(dt))
    log = ctypes.create_string_buffer(str(folder / "plumed.log").encode())
    libp.plumed_cmd(p, b"setLogFile", log)
    libp.plumed_cmd(p, b"init", None)
    load_line = ctypes.create_string_buffer(f"LOAD FILE={so}".encode())
    action = ctypes.create_string_buffer(
        (
            f"dseams: DSEAMS ATOMS=1-{n} DUMP={folder / 'live.lammpstrj'} "
            f"RECORD={folder / 'live.txt'} CUTOFF=3.5"
        ).encode()
    )
    libp.plumed_cmd(p, b"readInputLine", load_line)
    libp.plumed_cmd(p, b"readInputLine", action)
    masses = np.ones(n, dtype=np.float64)
    forces = np.zeros(3 * n, dtype=np.float64)
    virial = np.zeros(9, dtype=np.float64)
    from ase.md.langevin import Langevin
    from ase import units

    atoms = ase_atoms()
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=20.0, friction=0.01 / units.fs, fixcm=False)

    def tick():
        nm = np.ascontiguousarray(atoms.get_positions() / 10.0, dtype=np.float64)
        cell = np.ascontiguousarray(atoms.cell.array / 10.0, dtype=np.float64).reshape(9)
        istep = ctypes.c_int(int(dyn.nsteps))
        libp.plumed_cmd(p, b"setStep", ctypes.byref(istep))
        libp.plumed_cmd(p, b"setPositions", nm.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setMasses", masses.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setForces", forces.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setBox", cell.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setVirial", virial.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"calc", None)

    dyn.attach(tick, interval=1)
    dyn.run(6)
    libp.plumed_finalize(p)
    post = folder / "post.txt"
    run([REPLAY, "replay", str(folder / "live.lammpstrj"), str(post), "3.5"], folder)
    return same(folder / "live.txt", post)


def repeats(fn, n, warmup=1):
    for _ in range(warmup):
        fn()
    return [fn() for _ in range(n)]


def ipi_elapsed(analyse):
    folder = WORK / "ipi-bench" / ("live" if analyse else "plain")
    folder.mkdir(parents=True, exist_ok=True)
    port, xml = write_ipi(folder)
    proc = subprocess.Popen(
        [str(PREFIX / "bin" / "i-pi"), str(xml)],
        cwd=folder,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    time.sleep(1.0)
    if proc.poll() is not None:
        raise RuntimeError(proc.stdout.read())
    lib = load(LIB)
    started = time.perf_counter()
    with Engine(lib, str(folder / "rec.txt")) as eng:
        def on_frame(index, pos, cell):
            if analyse:
                eng.observe(index, pos, eng.box_from_rows(cell), str(folder / "dump.lammpstrj"))
        IpiClient("127.0.0.1", port).run(on_frame)
    elapsed = time.perf_counter() - started
    proc.communicate(timeout=60)
    return elapsed / 6


def plumed_elapsed(so, analyse):
    import ctypes
    import numpy as np
    from ase.md.langevin import Langevin
    from ase import units

    if not analyse:
        return ase_time(6, False)
    libp = ctypes.CDLL(str(PREFIX / "lib" / "libplumed.so"))

    class P(ctypes.Structure):
        _fields_ = [("p", ctypes.c_void_p)]

    libp.plumed_create.restype = P
    libp.plumed_cmd.argtypes = [P, ctypes.c_char_p, ctypes.c_void_p]
    libp.plumed_finalize.argtypes = [P]
    atoms = ase_atoms()
    n = len(atoms)
    p = libp.plumed_create()
    nat = ctypes.c_int(n)
    dt = ctypes.c_double(0.001)
    libp.plumed_cmd(p, b"setNatoms", ctypes.byref(nat))
    libp.plumed_cmd(p, b"setTimestep", ctypes.byref(dt))
    libp.plumed_cmd(p, b"init", None)
    action = ctypes.create_string_buffer(
        f"dseams: DSEAMS ATOMS=1-{n} DUMP={WORK / 'plumed-bench.lammpstrj'} RECORD={WORK / 'plumed-bench.txt'} CUTOFF=3.5".encode()
    )
    libp.plumed_cmd(p, b"readInputLine", ctypes.create_string_buffer(f"LOAD FILE={so}".encode()))
    libp.plumed_cmd(p, b"readInputLine", action)
    masses = np.ones(n, dtype=np.float64)
    forces = np.zeros(3 * n, dtype=np.float64)
    virial = np.zeros(9, dtype=np.float64)
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=20.0, friction=0.01 / units.fs, fixcm=False)

    def tick():
        nm = np.ascontiguousarray(atoms.get_positions() / 10.0, dtype=np.float64)
        cell = np.ascontiguousarray(atoms.cell.array / 10.0, dtype=np.float64).reshape(9)
        istep = ctypes.c_int(int(dyn.nsteps))
        libp.plumed_cmd(p, b"setStep", ctypes.byref(istep))
        libp.plumed_cmd(p, b"setPositions", nm.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setMasses", masses.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setForces", forces.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setBox", cell.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"setVirial", virial.ctypes.data_as(ctypes.c_void_p))
        libp.plumed_cmd(p, b"calc", None)

    dyn.attach(tick, interval=1)
    started = time.perf_counter()
    dyn.run(6)
    elapsed = time.perf_counter() - started
    libp.plumed_finalize(p)
    return elapsed / 6


def main():
    import shutil

    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True, exist_ok=True)
    info = host_block()
    info["lib_sha256"] = sha256(LIB)
    info["lmp_sha256"] = sha256(PREFIX / "bin" / "lmp")
    info["gmx_sha256"] = sha256(PREFIX / "bin" / "gmx")
    info["plugin_sha256"] = sha256(PLUGIN)
    out = {"host": info, "codes": {}}
    print("lammps", flush=True)
    ident = run_lammps()
    plain = time_lammps("plain", ident["steps"], 6, PLUGIN)
    live = time_lammps("live", ident["steps"], 6, PLUGIN)
    out["codes"]["lammps"] = {"frames": ident["frames"], "steps": ident["steps"], "plain": plain, "live": live}
    print("gromacs", flush=True)
    g = run_gromacs()
    gp, gl = time_gromacs(6, g["steps"])
    out["codes"]["gromacs"] = {"frames": g["frames"], "steps": g["steps"], "plain": gp, "live": gl}
    print("openmm", flush=True)
    frames = openmm_identity()
    op = repeats(lambda: openmm_time(20)[0], 6)
    ol = repeats(lambda: openmm_time(20)[1], 6)
    out["codes"]["openmm"] = {"frames": frames, "steps": 20, "plain": op, "live": ol}
    print("ase", flush=True)
    frames = ase_identity()
    ap = repeats(lambda: ase_time(20, False), 6)
    al = repeats(lambda: ase_time(20, True), 6)
    out["codes"]["ase"] = {"frames": frames, "steps": 20, "plain": ap, "live": al}
    print("ipi", flush=True)
    frames = run_ipi()
    ip = repeats(lambda: ipi_elapsed(False), 6)
    il = repeats(lambda: ipi_elapsed(True), 6)
    out["codes"]["ipi"] = {"frames": frames, "steps": 6, "plain": ip, "live": il}
    print("plumed", flush=True)
    so = compile_plumed()
    frames = run_plumed(so)
    pp = repeats(lambda: plumed_elapsed(so, False), 6)
    pl = repeats(lambda: plumed_elapsed(so, True), 6)
    out["codes"]["plumed"] = {"frames": frames, "steps": 6, "plain": pp, "live": pl}
    dest = WORK / "samples.json"
    dest.write_text(json.dumps(out))
    print(dest)


if __name__ == "__main__":
    main()
