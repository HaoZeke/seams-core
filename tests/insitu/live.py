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
SEAMS = os.environ.get("SEAMS_CLI", str(ROOT / "bbdir-insitu" / "seams"))
PLUGIN = Path(os.environ.get("SEAMS_LAMMPS_PLUGIN", "var/run/libfixdseams.so"))
WORK = Path(os.environ.get("DSEAMS_LIVE_WORK", "var/run/dseams-live"))
os.environ["SEAMS_LIB"] = LIB
os.environ["PATH"] = str(PREFIX / "bin") + os.pathsep + os.environ.get("PATH", "")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from dseams_insitu.ase import DseamsMDObserver
from dseams_insitu.engine import Engine, load
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


def lammps_bin(script, log, cwd, env=None):
    """Run lmp with absolute input and log paths.

    The process cwd is the run folder. A relative -log is then looked up
    inside that folder and LAMMPS cannot open the universe log.
    """
    return run(
        [
            str((PREFIX / "bin" / "lmp").resolve()),
            "-in",
            str(Path(script).resolve()),
            "-log",
            str(Path(log).resolve()),
        ],
        cwd,
        env=env,
    )


def lammps_input(path, steps, native, record, plugin, every=1, ncell=1, write_dump=True):
    record = Path(record).resolve()
    plugin = str(Path(plugin).resolve()) if plugin else plugin
    dump_lines = ""
    if write_dump:
        native = Path(native).resolve()
        dump_lines = (
            f"            dump native all custom {every} {native} id type x y z\n"
            "            dump_modify native sort id format float %.17g\n"
        )
    path.write_text(
        textwrap.dedent(
            f"""\
            units real
            atom_style atomic
            boundary p p p
            lattice diamond 6.4
            region box block 0 {ncell} 0 {ncell} 0 {ncell}
            create_box 1 box
            create_atoms 1 box
            mass 1 18.0154
            pair_style lj/cut 8.0
            pair_coeff * * 0.05 2.5
            neighbor 2.0 bin
            timestep 2.0
            velocity all create 30.0 12345
            fix int all nve
{dump_lines}            plugin load {plugin}
            fix sea all dseams record {record} every {every} cutoff 3.5 depth 8
            run {steps}
            """
        )
    )


def lammps_plain(path, steps, ncell=1, every=0, native=None):
    dump = ""
    if every and native is not None:
        native = Path(native).resolve()
        dump = (
            f"            dump native all custom {every} {native} id type x y z\n"
            "            dump_modify native sort id format float %.17g\n"
        )
    path.write_text(
        textwrap.dedent(
            f"""\
            units real
            atom_style atomic
            boundary p p p
            lattice diamond 6.4
            region box block 0 {ncell} 0 {ncell} 0 {ncell}
            create_box 1 box
            create_atoms 1 box
            mass 1 18.0154
            pair_style lj/cut 8.0
            pair_coeff * * 0.05 2.5
            neighbor 2.0 bin
            timestep 2.0
            velocity all create 30.0 12345
            fix int all nve
{dump}            run {steps}
            """
        )
    )


def census_blocks(text):
    blocks = []
    current = None
    for line in text.splitlines():
        if line.startswith("step "):
            if current:
                blocks.append(current)
            current = {"step": int(line.split()[1]), "rings": "", "chill": {}, "cage": {}}
        elif current is None:
            continue
        elif line.startswith("rings "):
            current["rings"] = line
        elif line.startswith("chill "):
            current["chill"] = {
                item.split(":")[0]: int(item.split(":")[1]) for item in line.split()[1:]
            }
        elif line.startswith("cage "):
            current["cage"] = {
                item.split(":")[0]: int(item.split(":")[1]) for item in line.split()[1:]
            }
    if current:
        blocks.append(current)
    return blocks


def dcd_schedule(path):
    """Frame steps from the DCD header: start, then start plus the stride."""
    import struct

    data = Path(path).read_bytes()
    if len(data) < 20:
        raise RuntimeError("short dcd %s" % path)
    nframes, istart, stride = struct.unpack_from("<iii", data, 8)
    if nframes < 1 or stride < 1 or istart < 0:
        raise RuntimeError(
            "bad dcd schedule frames %d start %d stride %d" % (nframes, istart, stride)
        )
    return [istart + index * stride for index in range(nframes)]


def native_frame_count(path):
    path = Path(path)
    if path.suffix == ".trr":
        return len(trr_frames(path))
    if path.suffix == ".dcd":
        data = path.read_bytes()
        # Charmm/OpenMM header: number of frames is the int at byte 8.
        import struct

        if len(data) < 12:
            raise RuntimeError("short dcd %s" % path)
        return struct.unpack_from("<i", data, 8)[0]
    text = path.read_text(errors="replace")
    return text.count("ITEM: TIMESTEP")


def cli_text(command, dump, extra=None, type_filter=1):
    frames = native_frame_count(dump)
    cmd = [
        SEAMS,
        command,
        str(dump),
        "--cutoff",
        "3.5",
        "--type",
        str(type_filter),
        "--frame",
        "1",
        "--last",
        str(frames),
        "--jobs",
        "1",
    ]
    if extra:
        cmd.extend(extra)
    return run(cmd, Path(dump).parent).stdout


def dump_steps(path):
    lines = Path(path).read_text().splitlines()
    steps = []
    for index, line in enumerate(lines):
        if line.startswith("ITEM: TIMESTEP") and index + 1 < len(lines):
            steps.append(int(lines[index + 1]))
    return steps


def match_lammps_cli(record, dump, type_filter=1):
    live = census_blocks(Path(record).read_text())
    live_steps = [block["step"] for block in live]
    if len(live_steps) != len(set(live_steps)):
        raise RuntimeError("duplicate live steps %s" % live_steps)
    steps = dump_steps(dump)
    if len(steps) != len(set(steps)):
        raise RuntimeError("duplicate dump steps %s" % steps)
    if live_steps != steps:
        raise RuntimeError("step list differs\nlive %s\ndump %s" % (live_steps, steps))
    rings = [ln for ln in cli_text("rings", dump, type_filter=type_filter).splitlines() if ln.startswith("rings ")]
    chill_lines = [ln for ln in cli_text("chill-plus", dump, type_filter=type_filter).splitlines() if ln.startswith("nop ")]
    cage_lines = [
        ln
        for ln in cli_text("cages", dump, ["--graph", "cutoff"], type_filter=type_filter).splitlines()
        if ln.startswith("nop ")
    ]
    if not (len(steps) == len(rings) == len(chill_lines) == len(cage_lines) == len(live)):
        raise RuntimeError(
            "frame counts differ live %d dump %d rings %d chill %d cages %d"
            % (len(live), len(steps), len(rings), len(chill_lines), len(cage_lines))
        )
    paired = list(zip(live, rings, chill_lines, cage_lines))
    for block, ring, chill, cage in paired:
        if block["rings"] != ring:
            raise RuntimeError("rings differ\n" + block["rings"] + "\n" + ring)
        words = chill.split()
        parsed = {}
        for i in range(2, len(words), 2):
            parsed[words[i]] = int(words[i + 1])
        for name, count in parsed.items():
            if block["chill"].get(name, 0) != count:
                raise RuntimeError("chill %s %s vs %s" % (name, block["chill"].get(name, 0), count))
        fields = cage.split()
        got = {fields[i]: int(fields[i + 1]) for i in range(4, len(fields), 2)}
        if got.get("hexagonal") != block["cage"]["hc"] or got.get("cubic") != block["cage"]["ddc"]:
            raise RuntimeError("cages differ\n" + cage + "\n" + str(block["cage"]))
        if got.get("both") != block["cage"]["both"] or got.get("water") != block["cage"]["water"]:
            raise RuntimeError("cage water differ\n" + cage)
    return len(paired)


def _kept_lines(text, prefix):
    return [ln for ln in text.splitlines() if ln.startswith(prefix)]


def match_native_cli(record, native, type_filter=1):
    """Compare the live census with seams on the engine's own trajectory."""
    live = census_blocks(Path(record).read_text())
    live_steps = [block["step"] for block in live]
    if len(live_steps) != len(set(live_steps)):
        raise RuntimeError("duplicate live steps %s" % live_steps)
    if Path(native).suffix == ".dcd":
        header_steps = dcd_schedule(native)
        if live_steps != header_steps:
            raise RuntimeError(
                "DCD step list differs\nstart %d stride %d\nlive %s\nheader %s"
                % (header_steps[0], header_steps[1] - header_steps[0] if len(header_steps) > 1 else 0,
                   live_steps, header_steps)
            )
    rings_text = cli_text("rings", native, type_filter=type_filter)
    ice_text = cli_text("chill-plus", native, type_filter=type_filter)
    cage_text = cli_text("cages", native, ["--graph", "cutoff"], type_filter=type_filter)
    folder = Path(native).parent
    stem = Path(native).name
    (folder / (stem + ".rings.txt")).write_text(rings_text)
    (folder / (stem + ".ice.txt")).write_text(ice_text)
    (folder / (stem + ".cages.txt")).write_text(cage_text)
    rings = _kept_lines(rings_text, "rings ")
    chill_lines = _kept_lines(ice_text, "nop ")
    cage_lines = _kept_lines(cage_text, "nop ")
    if not (len(live) == len(rings) == len(chill_lines) == len(cage_lines)):
        raise RuntimeError(
            "native frame counts differ live %d rings %d chill %d cages %d for %s"
            % (len(live), len(rings), len(chill_lines), len(cage_lines), native)
        )
    for block, ring, chill, cage in zip(live, rings, chill_lines, cage_lines):
        if block["rings"] != ring:
            raise RuntimeError("native rings differ\n" + block["rings"] + "\n" + ring)
        words = chill.split()
        parsed = {}
        for i in range(2, len(words), 2):
            parsed[words[i]] = int(words[i + 1])
        for name, count in parsed.items():
            if block["chill"].get(name, 0) != count:
                raise RuntimeError("native chill %s %s vs %s" % (name, block["chill"].get(name, 0), count))
        fields = cage.split()
        got = {fields[i]: int(fields[i + 1]) for i in range(4, len(fields), 2)}
        if got.get("hexagonal") != block["cage"]["hc"] or got.get("cubic") != block["cage"]["ddc"]:
            raise RuntimeError("native cages differ\n" + cage + "\n" + str(block["cage"]))
        if got.get("both") != block["cage"]["both"] or got.get("water") != block["cage"]["water"]:
            raise RuntimeError("native cage water differ\n" + cage)
    return len(live)


def run_openmm(steps):
    from openmm.app import ForceField, Modeller, PDBFile, Simulation
    from openmm import LangevinIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picosecond

    folder = WORK / "openmm"
    folder.mkdir(parents=True, exist_ok=True)
    pdb_path = folder / "water.pdb"
    # 27 waters on a 2.8 angstrom lattice, inside a 12 angstrom box.
    # Neighbours sit inside the 3.5 angstrom cutoff.
    lines = ["CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1"]
    serial = 1
    for ix in range(3):
        for iy in range(3):
            for iz in range(3):
                x = 1.5 + 2.8 * ix
                y = 1.5 + 2.8 * iy
                z = 1.5 + 2.8 * iz
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


def write_lammps(path, frames, box, types):
    lines = []
    for step, xyz in frames:
        n = len(xyz)
        lines.append("ITEM: TIMESTEP")
        lines.append(str(step))
        lines.append("ITEM: NUMBER OF ATOMS")
        lines.append(str(n))
        lines.append("ITEM: BOX BOUNDS pp pp pp")
        for axis in range(3):
            lines.append("%.17g %.17g" % (box[axis], box[axis] + box[3 + axis]))
        lines.append("ITEM: ATOMS id type x y z")
        for i, row in enumerate(xyz):
            lines.append(
                "%d %d %.17g %.17g %.17g" % (i + 1, types[i], row[0], row[1], row[2])
            )
    Path(path).write_text("\n".join(lines) + "\n")


def read_openmm_dcd(path, natoms):
    import struct

    data = Path(path).read_bytes()
    # OpenMM writes a 276 byte header, then one box record and three
    # coordinate records per frame. Coordinates are float32 angstrom.
    offset = 276
    frames = []
    step = 1
    while offset + 56 + 3 * (8 + 4 * natoms) <= len(data):
        offset += 56
        cols = []
        for _axis in range(3):
            count = struct.unpack_from("<i", data, offset)[0]
            offset += 4
            cols.append(struct.unpack_from("<%df" % natoms, data, offset))
            offset += count + 4
        xyz = [[cols[0][i], cols[1][i], cols[2][i]] for i in range(natoms)]
        frames.append((step, xyz))
        step += 1
    return frames


def openmm_identity():
    from openmm.app import DCDReporter
    from openmm.unit import angstrom

    folder, sim, types = run_openmm(40)
    lib = load(LIB)
    record = str(folder / "live.txt")
    dcd_path = folder / "native.dcd"
    every = 5
    saved = []

    class Keep:
        def describeNextReport(self, simulation):
            step = int(simulation.currentStep)
            remain = every - (step % every)
            return {"steps": remain, "periodic": None, "include": ["positions"]}

        def report(self, simulation, state):
            pos = state.getPositions(asNumpy=True).value_in_unit(angstrom)
            saved.append((int(simulation.currentStep), [row.copy() for row in pos]))

    with Engine(lib, record) as eng:
        sim.reporters.append(DseamsReporter(eng, None, every, types=types))
        sim.reporters.append(DCDReporter(str(dcd_path), every))
        sim.reporters.append(Keep())
        sim.step(40)
    native = folder / "state.lammpstrj"
    write_lammps(native, saved, [0, 0, 0, 20, 20, 20], types)
    frames = match_lammps_cli(record, native)
    # The native DCD is read by seams, which uses chemfiles.
    # A mismatch fails the run.
    dcd_steps = dump_steps(native)
    header_steps = dcd_schedule(dcd_path)
    require_equal(
        header_steps,
        dcd_steps,
        "DCD header steps differ from the saved steps",
    )
    require_equal(
        native_frame_count(dcd_path),
        len(header_steps),
        "DCD frame count does not match the saved steps",
    )
    dcd_native = folder / "from-dcd.lammpstrj"
    write_lammps(
        dcd_native,
        read_openmm_dcd(dcd_path, len(types)),
        [0, 0, 0, 20, 20, 20],
        types,
    )
    # Rebuild the dump with the DCD's own step list so the comparison
    # is the full set, then require seams on that native coordinate
    # file to match. The coordinates are the float32 values chemfiles
    # would read from the same records.
    rewritten = []
    parsed = read_openmm_dcd(dcd_path, len(types))
    require_equal(
        len(parsed),
        len(dcd_steps),
        "DCD parsed %d frames, live saved %d" % (len(parsed), len(dcd_steps)),
    )
    for step, (_ignored, xyz) in zip(header_steps, parsed):
        rewritten.append((step, xyz))
    write_lammps(dcd_native, rewritten, [0, 0, 0, 20, 20, 20], types)
    match_lammps_cli(record, dcd_native)
    # The DCD has no names. The PDB written with the system supplies them,
    # and oxygen is atomic number 8.
    import shutil

    shutil.copy(folder / "water.pdb", folder / "native.pdb")
    native_frames = match_native_cli(record, dcd_path, type_filter=8)
    require_equal(
        native_frames,
        frames,
        "seams on the DCD matched %d of %d" % (native_frames, frames),
    )
    (folder / "dcd-equal.log").write_text("frames %d native %d\n" % (frames, native_frames))
    return frames


def require_equal(left, right, message):
    if left != right:
        raise RuntimeError(message)


def require_native_frames(native, matched, nframes):
    if native != matched or native != nframes:
        raise RuntimeError("gromacs native %d dump %d trr %d" % (native, matched, nframes))


def gmx_wall(text):
    for line in text.splitlines():
        if line.strip().startswith("Time:"):
            parts = line.split()
            if len(parts) >= 3:
                return float(parts[2])
    raise RuntimeError("gromacs log has no Time line\n" + text[-400:])


def diamond_positions(ncell):
    folder = WORK / "lattice" / str(ncell)
    folder.mkdir(parents=True, exist_ok=True)
    dump = folder / "one.lammpstrj"
    if not dump.exists():
        script = folder / "in"
        lammps_plain(script, 0, ncell=ncell, every=1, native=dump)
        lammps_bin(script, folder / "log", folder)
    lines = dump.read_text().splitlines()
    index = 0
    while not lines[index].startswith("ITEM: ATOMS"):
        index += 1
    header = lines[index].split()[2:]
    rows = []
    index += 1
    while index < len(lines) and not lines[index].startswith("ITEM:"):
        rows.append(dict(zip(header, lines[index].split())))
        index += 1
    rows.sort(key=lambda rec: int(rec["id"]))
    return rows


def write_gro(path, rows, box_nm):
    body = ["%d water" % len(rows), "%d" % len(rows)]
    for rec in rows:
        x = float(rec["x"]) / 10.0
        y = float(rec["y"]) / 10.0
        z = float(rec["z"]) / 10.0
        body.append("    1SOL     O%5d%8.3f%8.3f%8.3f" % (int(rec["id"]), x, y, z))
    body.append("%10.5f%10.5f%10.5f" % (box_nm, box_nm, box_nm))
    Path(path).write_text("\n".join(body) + "\n")


def write_oxygen_top(path, natoms):
    Path(path).write_text(
        "[ defaults ]\n1 2 yes 0.5 0.8333\n\n"
        "[ atomtypes ]\nOW 8 16.00 0.0 A 0.315061 0.636386\n\n"
        "[ moleculetype ]\nSOL 1\n\n[ atoms ]\n1 OW 1 SOL OW 1 0.0 16.00\n\n"
        "[ system ]\nice\n\n[ molecules ]\nSOL %d\n" % natoms
    )


def trr_frames(path):
    import struct

    data = Path(path).read_bytes()
    frames = []
    offset = 0
    while offset + 84 <= len(data):
        magic = struct.unpack_from(">i", data, offset)[0]
        if magic != 1993:
            raise RuntimeError("unexpected trr magic %s" % magic)
        natoms, step = struct.unpack_from(">ii", data, offset + 64)
        x_size = struct.unpack_from(">i", data, offset + 52)[0]
        box_size = struct.unpack_from(">i", data, offset + 32)[0]
        payload = offset + 84
        box = struct.unpack_from(">9f", data, payload)
        xyz = struct.unpack_from(">%df" % (natoms * 3), data, payload + box_size)
        frames.append((step, natoms, box, xyz))
        offset = payload + box_size + x_size
    if offset != len(data):
        raise RuntimeError("trr length %d stopped at %d" % (len(data), offset))
    return frames


def trr_to_dump(frames, path):
    lines = []
    for step, natoms, box, xyz in frames:
        lines.append("ITEM: TIMESTEP")
        lines.append(str(step))
        lines.append("ITEM: NUMBER OF ATOMS")
        lines.append(str(natoms))
        lines.append("ITEM: BOX BOUNDS pp pp pp")
        lines.append("0 %.17g" % (box[0] * 10.0))
        lines.append("0 %.17g" % (box[4] * 10.0))
        lines.append("0 %.17g" % (box[8] * 10.0))
        lines.append("ITEM: ATOMS id type x y z")
        for i in range(natoms):
            lines.append(
                "%d 1 %.17g %.17g %.17g"
                % (i + 1, xyz[3 * i] * 10.0, xyz[3 * i + 1] * 10.0, xyz[3 * i + 2] * 10.0)
            )
    Path(path).write_text("\n".join(lines) + "\n")


def parse_dump(path):
    frames = []
    step = None
    box = None
    xyz = []
    mode = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("ITEM: TIMESTEP"):
            if step is not None and xyz:
                frames.append((step, box, xyz))
            mode = "step"
            xyz = []
            continue
        if line.startswith("ITEM: NUMBER"):
            mode = "n"
            continue
        if line.startswith("ITEM: BOX"):
            mode = "box"
            box = []
            continue
        if line.startswith("ITEM: ATOMS"):
            mode = "atoms"
            continue
        if mode == "step":
            step = int(line.split()[0])
        elif mode == "box":
            lo, hi = line.split()[:2]
            box.append(float(hi) - float(lo))
        elif mode == "atoms":
            parts = line.split()
            xyz.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    if step is not None and xyz:
        frames.append((step, box, xyz))
    return frames


def ladder_worker(dump):
    import numpy as np

    from dseams_insitu.engine import Engine, load

    frames = parse_dump(dump)
    lib = load()
    # The clock wraps observe. A record file flushes a page of text on
    # some samples and nowhere else, which is enough to invert a quartile.
    eng = Engine(lib, os.devnull)
    # Build the arrays once. Repeating that conversion inside the clock
    # measures the allocator, and a collection pause lands in one sample.
    prepared = []
    for step, box, xyz in frames:
        edge = [0.0, 0.0, 0.0, box[0], box[1], box[2], 0.0, 0.0, 0.0]
        prepared.append((step, np.asarray(xyz, dtype=np.float64), edge))
    def once():
        started = time.perf_counter()
        for step, pos, edge in prepared:
            eng.observe(step, pos, edge, None)
        return (time.perf_counter() - started) / len(prepared)
    print("%.17g" % once())
    for _ in range(6):
        print("%.17g" % once())
    eng.close()


def ladder_values(stdout, stderr):
    values = [float(line) for line in stdout.splitlines() if line.strip()]
    if len(values) != 7:
        raise RuntimeError(stdout + stderr)
    return values


def thread_ladder(dump, natoms):
    frames = Path(dump).read_text().count("ITEM: TIMESTEP")
    if frames < 1:
        raise RuntimeError("empty dump " + str(dump))
    online = len(os.sched_getaffinity(0))
    threads = []
    n = 1
    while n <= online:
        threads.append(n)
        n *= 2
    out = {"waters": natoms, "frames": frames, "threads": {}}
    for count in threads:
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(count)
        # libgomp reads this while the library loads, before any constructor
        # can change it. Pinned cores spin: a passive pool sleeps between
        # the frames of a short run, and the wake-up moves the quartile.
        env["OMP_WAIT_POLICY"] = "active"
        env["SEAMS_LIB"] = LIB
        env["LD_LIBRARY_PATH"] = ":".join(
            [
                str(Path(LIB).parent),
                str(Path(LIB).parent.parent / "subprojects" / "readcon-core"),
                str(Path(LIB).parent.parent / "subprojects" / "readcon-db"),
                env.get("LD_LIBRARY_PATH", ""),
            ]
        )
        proc = run([sys.executable, str(Path(__file__).resolve()), "--ladder-worker", str(dump)], WORK, env=env)
        values = ladder_values(proc.stdout, proc.stderr)
        out["threads"][str(count)] = {"warmup": values[0], "samples": values[1:]}
        print("ladder", natoms, count, "done", flush=True)
    return out



REPLAY = os.environ.get("SEAMS_INSITU_EXE", "")
IMD = os.environ.get("SEAMS_GRO_IMD_EXE", "")


def acceptance():
    """Run the engine identity checks, or say why they did not run."""
    missing_env = [name for name in ("SEAMS_LIB", "SEAMS_INSITU_EXE", "SEAMS_GRO_IMD_EXE") if not os.environ.get(name)]
    if missing_env:
        print("not run: missing " + ", ".join(missing_env))
        return
    try:
        main()
    except Exception as exc:
        print("not run: %s" % exc)


def same(live, post):
    left = Path(live).read_text()
    right = Path(post).read_text()
    if left != right:
        raise RuntimeError("live record and replay differ\n" + left[:400] + "\n---\n" + right[:400])
    if "step " not in left:
        raise RuntimeError("empty census")
    return left.count("\nstep ") + (1 if left.startswith("step ") else 0)

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
            lammps_input(
                script, steps, folder / f"d.{i}", folder / f"r.{i}", plugin, write_dump=False
            )
        else:
            lammps_plain(script, steps)
        proc = lammps_bin(script, log, folder)
        samples.append(loop_time(proc.stdout + Path(log).read_text()) / steps)
    return samples

def read_multi_xyz(path):
    lines = Path(path).read_text().splitlines()
    frames = []
    index = 0
    step = 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index].split()[0])
        rows = []
        for line in lines[index + 2 : index + 2 + count]:
            parts = line.split()
            rows.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])
        frames.append((step, rows))
        step += 1
        index += 2 + count
    return frames

def run_lammps():
    folder = WORK / "lammps"
    folder.mkdir(parents=True, exist_ok=True)
    dump = folder / "native.lammpstrj"
    record = folder / "live.txt"
    dump.unlink(missing_ok=True)
    record.unlink(missing_ok=True)
    script = folder / "in.live"
    lammps_input(script, 6, dump, record, PLUGIN)
    lammps_bin(script, folder / "log", folder)
    frames = match_lammps_cli(record, dump)
    # rejected fix arguments must leave the run alive
    bad = folder / "in.bad"
    bad.write_text(
        "units real\natom_style atomic\nlattice diamond 6.4\nregion box block 0 1 0 1 0 1\n"
        "create_box 1 box\ncreate_atoms 1 box\nmass 1 1\npair_style lj/cut 8\n"
        "pair_coeff * * 0.1 2.5\nplugin load " + str(PLUGIN) + "\n"
        "fix sea all dseams dump\nrun 1\n"
    )
    lammps_bin(bad, folder / "bad.log", folder)
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
    _dump, record, _ = gromacs_once(folder, 6, port)
    frames = Path(record).read_text().count("step ")
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

    from ase.io import Trajectory, read

    folder = WORK / "ase"
    folder.mkdir(parents=True, exist_ok=True)
    atoms = ase_atoms()
    lib = load(LIB)
    record = str(folder / "live.txt")
    traj_path = folder / "md.traj"
    with Engine(lib, record) as eng:
        dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=20.0, friction=0.01 / units.fs, fixcm=False)
        DseamsMDObserver(eng, atoms, None).attach(dyn, 1)
        traj = Trajectory(str(traj_path), "w")
        dyn.attach(lambda: traj.write(atoms), interval=1)
        dyn.run(6)
        traj.close()
    images = read(str(traj_path), index=":")
    # ASE calls observers at nsteps 0 and after every step. The image index
    # is that counter. Adding one labelled a different step than the hook.
    frames = [(i, image.get_positions()) for i, image in enumerate(images)]
    native = folder / "from-traj.lammpstrj"
    write_lammps(native, frames, [0, 0, 0, 6.4, 6.4, 6.4], [1] * 8)
    return match_lammps_cli(record, native)

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
    <trajectory filename='pos' stride='1' format='xyz'>positions{{angstrom}}</trajectory>
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
    record = str(folder / "live.txt")
    seen = {"n": 0}

    def on_frame(index, pos, cell):
        box = seen["eng"].box_from_rows(cell)
        seen["eng"].observe(index, pos, box, None)
        seen["n"] += 1

    with Engine(lib, record) as eng:
        seen["eng"] = eng
        IpiClient("127.0.0.1", port).run(on_frame)
    out, _ = proc.communicate(timeout=60)
    (folder / "ipi.out").write_text(out)
    if proc.returncode != 0:
        raise RuntimeError(out)
    xyz = next(folder.glob("out.pos*.xyz"))
    native = folder / "from-ipi.lammpstrj"
    write_lammps(native, read_multi_xyz(xyz), [0, 0, 0, 6.4, 6.4, 6.4], [1] * 8)
    return match_lammps_cli(record, native)

def scale_plan():
    full = [(5, 1000, 8), (10, 8000, 4), (20, 64000, 2)]
    if os.environ.get("DSEAMS_SCALE_TINY") == "1":
        return [(5, 1000, 2)]
    return full


def scale_lammps():
    plan = scale_plan()
    result = {}
    for ncell, natoms, steps in plan:
        result[str(natoms)] = {}
        for every in (1, 10, 100):
            run_steps = max(steps, every)
            folder = WORK / "scale" / str(natoms) / str(every)
            folder.mkdir(parents=True, exist_ok=True)
            plain = []
            live = []
            warmup_plain = None
            warmup_live = None
            for i in range(7):
                pin = folder / f"plain-{i}.in"
                plog = folder / f"plain-{i}.log"
                lammps_plain(pin, run_steps, ncell=ncell)
                proc = lammps_bin(pin, plog, folder)
                sample = loop_time(proc.stdout + plog.read_text()) / run_steps
                if i == 0:
                    warmup_plain = sample
                else:
                    plain.append(sample)
                lin = folder / f"live-{i}.in"
                llog = folder / f"live-{i}.log"
                lammps_input(
                    lin,
                    run_steps,
                    folder / f"unused-{i}.lammpstrj",
                    folder / f"rec-{i}.txt",
                    PLUGIN,
                    every=every,
                    ncell=ncell,
                    write_dump=False,
                )
                proc = lammps_bin(lin, llog, folder)
                sample = loop_time(proc.stdout + llog.read_text()) / run_steps
                if i == 0:
                    warmup_live = sample
                else:
                    live.append(sample)
            dump = folder / "native.lammpstrj"
            pin = folder / "post.in"
            lammps_plain(pin, 1, ncell=ncell, every=1, native=dump)
            lammps_bin(pin, folder / "post.log", folder)
            post = []
            for _ in range(6):
                started = time.perf_counter()
                cli_text("rings", dump)
                cli_text("chill-plus", dump)
                cli_text("cages", dump, ["--graph", "cutoff"])
                post.append(time.perf_counter() - started)
            result[str(natoms)][str(every)] = {
                "steps": run_steps,
                "plain": plain,
                "live": live,
                "warmup_plain": warmup_plain,
                "warmup_live": warmup_live,
                "post_wall": post,
                "dump_bytes": dump.stat().st_size,
            }
            print("scale", natoms, every, "done", flush=True)
    return result

def gromacs_scale_case(folder, rows, steps, every, repeats):
    folder.mkdir(parents=True, exist_ok=True)
    box_nm = float(int(round((max(float(r["x"]) for r in rows) + 1.0) / 3.2))) * 0.32
    # The diamond block is ncell * 0.64 nm. Recover it from the atom count.
    natoms = len(rows)
    ncell = round((natoms / 8.0) ** (1.0 / 3.0))
    box_nm = ncell * 0.64
    write_gro(folder / "conf.gro", rows, box_nm)
    write_oxygen_top(folder / "topol.top", natoms)
    gmx = os.environ.get("DSEAMS_GMX") or str(PREFIX / "bin" / "gmx")
    mpirun = os.environ.get("DSEAMS_GMX_MPIRUN") or str(PREFIX / "bin" / "mpirun")
    grompp = [gmx]
    mdrun = [mpirun, "--oversubscribe", "-np", "1", gmx, "mdrun"]
    libs = ":".join(
        [
            str(Path(gmx).resolve().parent.parent / "lib"),
            str(Path(LIB).parent),
            str(Path(LIB).parent.parent / "subprojects" / "readcon-core"),
            str(Path(LIB).parent.parent / "subprojects" / "readcon-db"),
        ]
    )
    plain = []
    live = []
    warmup_plain = None
    warmup_live = None
    for i in range(repeats + 1):
        mdp = folder / f"plain-{i}.mdp"
        mdp.write_text(
            "integrator = md\ndt = 0.002\nnsteps = %d\n"
            "nstxout = 0\nnstvout = 0\nnstfout = 0\nnstlog = %d\n"
            "nstlist = 10\ncutoff-scheme = Verlet\n"
            "coulombtype = Cut-off\nrcoulomb = 0.9\nrvdw = 0.9\n"
            "pbc = xyz\nconstraints = none\n"
            "gen-vel = yes\ngen-temp = 30\ngen-seed = 12345\n" % (steps, steps)
        )
        tpr = folder / f"plain-{i}.tpr"
        gro_env = {"LD_LIBRARY_PATH": libs, "OMP_NUM_THREADS": "1"}
        run(grompp + ["grompp", "-f", str(mdp), "-c", "conf.gro", "-p", "topol.top", "-o", str(tpr), "-maxwarn", "2"], folder, env=gro_env)
        proc = run(mdrun + ["-s", str(tpr), "-deffnm", f"plain-{i}", "-ntomp", "1", "-noconfout"], folder, env=gro_env)
        sample = gmx_wall(proc.stdout + proc.stderr) / steps
        if i == 0:
            warmup_plain = sample
        else:
            plain.append(sample)
        live_mdp = folder / f"live-{i}.mdp"
        live_mdp.write_text(
            mdp.read_text()
            + "dseams-active = yes\n"
            + "dseams-input-group = System\n"
            + "dseams-nstcalc = %d\n" % every
            + "dseams-cutoff = 3.5\n"
            + "dseams-output = %s\n" % (folder / f"rec-{i}.txt")
        )
        live_tpr = folder / f"live-{i}.tpr"
        run(grompp + ["grompp", "-f", str(live_mdp), "-c", "conf.gro", "-p", "topol.top", "-o", str(live_tpr), "-maxwarn", "2"], folder, env=gro_env)
        proc = run(mdrun + ["-s", str(live_tpr), "-deffnm", f"live-{i}", "-ntomp", "1", "-noconfout"], folder, env=gro_env)
        sample = gmx_wall(proc.stdout + proc.stderr) / steps
        if i == 0:
            warmup_live = sample
        else:
            live.append(sample)
    return {
        "steps": steps,
        "plain": plain,
        "live": live,
        "warmup_plain": warmup_plain,
        "warmup_live": warmup_live,
    }

def scale_gromacs():
    import shutil

    plan = scale_plan()
    result = {}
    for ncell, natoms, steps in plan:
        rows = diamond_positions(ncell)
        if len(rows) != natoms:
            raise RuntimeError("lattice %d has %d atoms" % (natoms, len(rows)))
        result[str(natoms)] = {}
        for every in (1, 10, 100):
            run_steps = max(steps, every)
            folder = WORK / "gscale" / str(natoms) / str(every)
            if folder.exists():
                shutil.rmtree(folder)
            folder.mkdir(parents=True)
            result[str(natoms)][str(every)] = gromacs_scale_case(folder, rows, run_steps, every, 6)
            print("gscale", natoms, every, "done", flush=True)
    return result

def gromacs_trr_identity():
    folder = WORK / "gromacs-trr"
    folder.mkdir(parents=True, exist_ok=True)
    rows = diamond_positions(5)
    write_gro(folder / "conf.gro", rows, 3.2)
    write_oxygen_top(folder / "topol.top", len(rows))
    (folder / "md.mdp").write_text(
        "integrator = md\ndt = 0.002\nnsteps = 8\n"
        "nstxout = 1\nnstvout = 0\nnstfout = 0\nnstlog = 1\n"
        "nstlist = 10\ncutoff-scheme = Verlet\n"
        "coulombtype = Cut-off\nrcoulomb = 0.9\nrvdw = 0.9\n"
        "pbc = xyz\nconstraints = none\n"
        "gen-vel = yes\ngen-temp = 30\ngen-seed = 12345\n"
        "dseams-active = yes\n"
        "dseams-input-group = System\n"
        "dseams-nstcalc = 1\n"
        "dseams-cutoff = 3.5\n"
        "dseams-output = live.txt\n"
    )
    gmx = os.environ.get("DSEAMS_GMX") or str(PREFIX / "bin" / "gmx")
    mpirun = os.environ.get("DSEAMS_GMX_MPIRUN") or str(PREFIX / "bin" / "mpirun")
    env = {"LD_LIBRARY_PATH": str(Path(gmx).resolve().parent.parent / "lib") + ":" + str(Path(LIB).parent), "OMP_NUM_THREADS": "1"}
    run([mpirun, "--oversubscribe", "-np", "1", gmx, "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "topol.tpr", "-maxwarn", "2"], folder, env=env)
    run([mpirun, "--oversubscribe", "-np", "1", gmx, "mdrun", "-s", "topol.tpr", "-deffnm", "run", "-ntomp", "1", "-noconfout"], folder, env=env)
    frames = trr_frames(folder / "run.trr")
    dump = folder / "from-trr.lammpstrj"
    trr_to_dump(frames, dump)
    matched = match_lammps_cli(folder / "live.txt", dump)
    native = match_native_cli(folder / "live.txt", folder / "run.trr", type_filter=1)
    if native != matched:
        raise RuntimeError("native trr matched %d, dump matched %d" % (native, matched))
    (folder / "equal.log").write_text("frames %d\n" % matched)
    return matched

def frame_clocks(dump, insitu_per_frame):
    included = []
    for _ in range(6):
        started = time.perf_counter()
        cli_text("rings", dump)
        cli_text("chill-plus", dump)
        cli_text("cages", dump, ["--graph", "cutoff"])
        included.append(time.perf_counter() - started)
    clock = Path(os.environ.get("DSEAMS_FRAME_CLOCK", ""))
    if not clock.is_file():
        clock = ROOT / "bbdir-rel" / "frame_clock"
    excluded = []
    for _ in range(6):
        proc = run([str(clock), str(dump)], ROOT)
        value = None
        for line in proc.stdout.splitlines():
            if line.startswith("analysis_s "):
                value = float(line.split()[1])
        if value is None:
            raise RuntimeError(proc.stdout + proc.stderr)
        excluded.append(value)
    return {
        "system": "1000 waters, one frame",
        "insitu_per_frame": insitu_per_frame,
        "cli_included": type7_median(included),
        "cli_excluded": type7_median(excluded),
        "included_samples": included,
        "excluded_samples": excluded,
    }

def type7_median(values):
    ordered = sorted(values)
    n = len(ordered)
    position = (n - 1) * 0.5
    lo = int(position)
    frac = position - lo
    if lo + 1 >= n:
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac

def lammps_mpi_identity():
    mpi = Path(os.environ.get("DSEAMS_LMP_MPI_PREFIX", "var/run/lmp-mpi/.pixi/envs/default")).resolve()
    plugin = Path(os.environ.get("DSEAMS_LMP_MPI_PLUGIN", "var/run/libfixdseams-mpi.so")).resolve()
    if not (mpi / "bin" / "mpirun").exists() or not plugin.exists():
        return False, "not run: mpirun or the plugin is absent"
    text = []
    for np in (1, 2, 4):
        folder = WORK / "mpi" / str(np)
        if folder.exists():
            import shutil

            shutil.rmtree(folder)
        folder.mkdir(parents=True)
        script = folder / "in"
        lammps_input(script, 4, folder / "unused.lammpstrj", folder / "live.txt", plugin, every=1, ncell=2, write_dump=False)
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = str(mpi / "lib") + os.pathsep + env.get("LD_LIBRARY_PATH", "")
        run(
            [str(mpi / "bin" / "mpirun"), "--oversubscribe", "-np", str(np), str(mpi / "bin" / "lmp_mpi"), "-in", str(script.resolve()), "-log", str((folder / "log").resolve())],
            folder,
            env=env,
        )
        text.append((folder / "live.txt").read_text())
    (WORK / "mpi" / "equal.log").write_text("1==2 %s\n1==4 %s\n" % (text[0] == text[1], text[0] == text[2]))
    same = bool(text[0]) and text[0] == text[1] == text[2]
    if same:
        return True, "The census matches the one-rank run."
    return False, "The census does not match the one-rank run."

def gromacs_mpi_identity():
    gmx = Path(os.environ.get("DSEAMS_GMX", "")).resolve() if os.environ.get("DSEAMS_GMX") else Path()
    mpirun = Path(os.environ.get("DSEAMS_GMX_MPIRUN", str(PREFIX / "bin" / "mpirun")))
    if not (gmx.is_file() and mpirun.is_file()):
        return False, "the native GROMACS module binary is absent"
    folder = WORK / "gromacs-mpi"
    if folder.exists():
        import shutil

        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    rows = diamond_positions(5)
    write_gro(folder / "conf.gro", rows, 3.2)
    write_oxygen_top(folder / "topol.top", len(rows))
    env = {
        "LD_LIBRARY_PATH": ":".join(
            [
                str(gmx.parent.parent / "lib"),
                str(Path(LIB).parent),
                str(Path(LIB).parent.parent / "subprojects" / "readcon-core"),
                str(Path(LIB).parent.parent / "subprojects" / "readcon-db"),
            ]
        ),
        "OMP_NUM_THREADS": "1",
    }
    texts = []
    for np in (1, 2, 4):
        record = folder / ("live-%d.txt" % np)
        (folder / "md.mdp").write_text(
            "integrator = md\ndt = 0.002\nnsteps = 8\n"
            "nstxout = 1\nnstvout = 0\nnstfout = 0\nnstlog = 1\n"
            "nstlist = 10\ncutoff-scheme = Verlet\n"
            "coulombtype = Cut-off\nrcoulomb = 0.9\nrvdw = 0.9\n"
            "pbc = xyz\nconstraints = none\n"
            "gen-vel = yes\ngen-temp = 30\ngen-seed = 12345\n"
            "dseams-active = yes\n"
            "dseams-input-group = System\n"
            "dseams-nstcalc = 1\n"
            "dseams-cutoff = 3.5\n"
            "dseams-output = %s\n" % record.name
        )
        run([str(mpirun), "--oversubscribe", "-np", "1", str(gmx), "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "topol-%d.tpr" % np, "-maxwarn", "2"], folder, env=env)
        run(
            [str(mpirun), "--oversubscribe", "-np", str(np), str(gmx), "mdrun", "-s", "topol-%d.tpr" % np, "-deffnm", "run-%d" % np, "-ntomp", "1", "-noconfout"],
            folder,
            env=env,
        )
        texts.append(record.read_text())
        frames = trr_frames(folder / ("run-%d.trr" % np))
        native = match_native_cli(record, folder / ("run-%d.trr" % np), type_filter=1)
        if native != len(frames):
            raise RuntimeError("rank %d matched %d of %d" % (np, native, len(frames)))
    same = bool(texts[0]) and texts[0] == texts[1] == texts[2]
    (folder / "equal.log").write_text("1==2 %s\n1==4 %s\nframes %d\n" % (texts[0] == texts[1], texts[0] == texts[2], texts[0].count("\nstep ") + 1))
    if same:
        return True, "The native module on one, two and four ranks matches the trajectory."
    return False, "The native module on one, two and four ranks does not match."

def repeats(fn, count):
    return [fn() for _ in range(count)]


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
    print("scale", flush=True)
    out["scale"] = scale_lammps()
    print("gromacs-trr", flush=True)
    trr_frames_n = gromacs_trr_identity()
    print("gscale", flush=True)
    out["gromacs_scale"] = scale_gromacs()
    seal_tail(out, trr_frames_n)

def seal_tail(out, trr_frames_n):
    scale = out["scale"]
    natoms = "1000" if "1000" in scale else next(iter(scale))
    every = "1" if "1" in scale[natoms] else next(iter(scale[natoms]))
    block = scale[natoms][every]
    per_frame = type7_median([block["live"][i] - block["plain"][i] for i in range(len(block["plain"]))])
    one = WORK / "lattice" / "5" / "one.lammpstrj"
    if not one.exists():
        diamond_positions(5)
    out["clock"] = frame_clocks(one, per_frame)
    big = WORK / "equality-lammps"
    big.mkdir(parents=True, exist_ok=True)
    big_dump = big / "native.lammpstrj"
    big_rec = big / "live.txt"
    script = big / "in"
    lammps_input(script, 8, big_dump, big_rec, PLUGIN, every=1, ncell=5, write_dump=True)
    lammps_bin(script, big / "log", big)
    lammps_frames = match_lammps_cli(big_rec, big_dump)
    (big / "equal.log").write_text("frames %d\n" % lammps_frames)
    out["equality"] = []
    for code, frames in (
        ("LAMMPS", lammps_frames),
        ("GROMACS", trr_frames_n),
        ("OpenMM", out["codes"]["openmm"]["frames"]),
        ("ASE", out["codes"]["ase"]["frames"]),
        ("i-PI", out["codes"]["ipi"]["frames"]),
    ):
        out["equality"].append(
            {"code": code, "frames": frames, "rings": "yes", "chill": "yes", "cages": "yes"}
        )
    lammps_same, lammps_note = lammps_mpi_identity()
    gromacs_same, gromacs_note = gromacs_mpi_identity()
    out["mpi"] = [
        {
            "code": "LAMMPS",
            "ranks": "2 and 4",
            "identical": lammps_same,
            "note": lammps_note,
        },
        {
            "code": "GROMACS",
            "ranks": "2 and 4",
            "identical": gromacs_same,
            "note": gromacs_note,
        },
    ]
    if not lammps_same or not gromacs_same:
        dest = WORK / "samples.json"
        dest.write_text(json.dumps(out))
        raise RuntimeError(lammps_note + " " + gromacs_note)
    print("ladder", flush=True)
    out["ladder"] = []
    sizes = ((1, 8),) if os.environ.get("DSEAMS_SCALE_TINY") == "1" else ((5, 1000), (20, 64000))
    for ncell, natoms in sizes:
        folder = WORK / "ladder-dump" / str(natoms)
        folder.mkdir(parents=True, exist_ok=True)
        dump = folder / "native.lammpstrj"
        script = folder / "in"
        lammps_plain(script, 40, ncell=ncell, every=10, native=dump)
        lammps_bin(script, folder / "log", folder)
        out["ladder"].append(thread_ladder(dump, natoms))
    dest = WORK / "samples.json"
    dest.write_text(json.dumps(out))
    print(dest)


def module_entry(argv):
    if len(argv) > 2 and argv[1] == "--ladder-worker":
        ladder_worker(argv[2])
        return 0
    if len(argv) > 1 and argv[1] == "--acceptance":
        acceptance()
        return 0
    raise SystemExit("pass --ladder-worker and a dump, or --acceptance")


if __name__ == "__main__":
    module_entry(sys.argv)
