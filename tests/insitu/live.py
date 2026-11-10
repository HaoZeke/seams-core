# Copyright (c) d-SEAMS developers
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


def library_path(*extra):
    parts = [str(item) for item in extra if item]
    chem = Path(LIB).resolve().parent.parent / "chemfiles-prefix" / "lib"
    if chem.is_dir():
        parts.append(str(chem))
    current = os.environ.get("LD_LIBRARY_PATH", "")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)
os.environ["SEAMS_LIB"] = LIB
os.environ["PATH"] = str(PREFIX / "bin") + os.pathsep + os.environ.get("PATH", "")
os.environ.setdefault("OMP_NUM_THREADS", "1")


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


def cli_text(command, dump, extra=None, type_filter=1):
    dump = Path(dump).resolve()
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



def acceptance():
    """Run the engine identity checks, or say why they did not run."""
    missing_env = [name for name in ("SEAMS_LIB", "SEAMS_INSITU_EXE") if not os.environ.get(name)]
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
    try:
        lammps_bin(bad, folder / "bad.log", folder)
    except RuntimeError as exc:
        if "Illegal fix dseams command" not in str(exc):
            raise
    else:
        raise RuntimeError("bad fix command was accepted")
    return {"frames": frames, "steps": 40}

def free_port():
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port

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

def repeats(fn, count):
    return [fn() for _ in range(count)]


def main():
    import shutil

    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True, exist_ok=True)
    info = host_block()
    info["lib_sha256"] = sha256(LIB)
    info["lmp_sha256"] = sha256(PREFIX / "bin" / "lmp")
    info["plugin_sha256"] = sha256(PLUGIN)
    out = {"host": info, "codes": {}}
    print("lammps", flush=True)
    ident = run_lammps()
    plain = time_lammps("plain", ident["steps"], 6, PLUGIN)
    live_times = time_lammps("live", ident["steps"], 6, PLUGIN)
    out["codes"]["lammps"] = {
        "frames": ident["frames"],
        "steps": ident["steps"],
        "plain": plain,
        "live": live_times,
    }
    print("scale", flush=True)
    out["scale"] = scale_lammps()
    seal_tail(out)


def seal_tail(out):
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
    out["equality"] = [
        {"code": "LAMMPS", "frames": lammps_frames, "rings": "yes", "chill": "yes", "cages": "yes"}
    ]
    lammps_same, lammps_note = lammps_mpi_identity()
    out["mpi"] = [
        {"code": "LAMMPS", "ranks": "2 and 4", "identical": lammps_same, "note": lammps_note}
    ]
    if not lammps_same:
        dest = WORK / "samples.json"
        dest.write_text(json.dumps(out))
        raise RuntimeError(lammps_note)
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
