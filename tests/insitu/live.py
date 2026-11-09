"""Live water and ice runs. The in-situ record must match a replay of the dump."""

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

from dseams_insitu.engine import Engine, load
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


def lammps_input(path, steps, native, record, plugin, every=1, ncell=1, write_dump=True):
    dump_lines = ""
    if write_dump:
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
    require_equal(
        native_frame_count(dcd_path),
        len(dcd_steps),
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
    for step, (_ignored, xyz) in zip(dcd_steps, parsed):
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
        run([str(PREFIX / "bin" / "lmp"), "-in", str(script), "-log", str(folder / "log")], folder)
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


def module_entry(argv):
    if len(argv) > 2 and argv[1] == "--ladder-worker":
        ladder_worker(argv[2])
        return 0
    raise SystemExit("pass --ladder-worker and a dump")


if __name__ == "__main__": module_entry(sys.argv)
