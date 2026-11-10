# Copyright (c) d-SEAMS developers
"""Build the LAMMPS adapter with coverage and exercise it.

Tool paths come from the environment. Sources are taken from this repo.
"""

import os
import subprocess
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def need(name):
    value = os.environ.get(name, "")
    if not value:
        raise SystemExit("set %s" % name)
    return value


def command_text(proc, cwd):
    log_path = Path(cwd) / "log.lammps"
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    return (proc.stdout or "") + (proc.stderr or "") + log


def run(cmd, cwd, env=None):
    proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        raise SystemExit("exit %d want 0\n%s" % (proc.returncode, command_text(proc, cwd)[-1500:]))
    return proc


def expect(cmd, cwd, env, code, message):
    proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True)
    text = command_text(proc, cwd)
    if proc.returncode != code:
        raise SystemExit("exit %d want %d\n%s" % (proc.returncode, code, text[-1500:]))
    if message not in text:
        raise SystemExit("missing %s\n%s" % (message, text[-1500:]))
    return proc


def lammps_cases(work, plugin, lmp, mpirun, libdir):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = libdir + ":" + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    header = textwrap.dedent(
        """\
        units real
        atom_style atomic
        boundary p p p
        lattice diamond 6.4
        region box block 0 2 0 2 0 2
        create_box 2 box
        create_atoms 1 box
        create_atoms 2 single 6.0 6.0 6.0
        mass * 18.0
        pair_style lj/cut 8.0
        pair_coeff * * 0.05 2.5
        neighbor 2.0 bin
        timestep 2.0
        velocity all set 0 0 0
        plugin load %s
        """
        % plugin
    )
    rejects = {
        "bad-every": (header + "fix sea all dseams record rec.txt every 0\nrun 1\n", 1, "fix dseams: every must be a positive integer"),
        "bad-cutoff": (header + "fix sea all dseams record rec.txt cutoff 0\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "bad-depth": (header + "fix sea all dseams record rec.txt depth 0\nrun 1\n", 1, "fix dseams: depth must be a positive integer"),
        "thin-depth": (header + "fix sea all dseams record rec.txt depth 1\nrun 1\n", 1, "fix dseams: record did not open"),
        "bad-double": (header + "fix sea all dseams record rec.txt cutoff nope\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "bad-word": (header + "fix sea all dseams record rec.txt nope 1\nrun 1\n", 1, "Illegal fix dseams command"),
        "bad-host": (header + "fix sea all dseams record rec.txt host maybe\nrun 1\n", 1, "fix dseams: host is yes or no"),
        "no-record": (header + "fix sea all dseams every 1 cutoff 3.5\nrun 1\n", 1, "fix dseams: record is required"),
        "bad-record": (header + "fix sea all dseams record /no/such/dir/rec.txt every 1\nrun 1\n", 1, "fix dseams: record did not open"),
        "big-every": (header + "fix sea all dseams record rec.txt every 1000001\nrun 1\n", 1, "fix dseams: every must be a positive integer"),
        "junk-every": (header + "fix sea all dseams record rec.txt every 2x\nrun 1\n", 1, "fix dseams: every must be a positive integer"),
        "neg-cutoff": (header + "fix sea all dseams record rec.txt cutoff -3\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "big-depth": (header + "fix sea all dseams record rec.txt depth 1000001\nrun 1\n", 1, "fix dseams: depth must be a positive integer"),
        "bare-dump": (header + "fix sea all dseams dump\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-every": (header + "fix sea all dseams record rec.txt every\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-record": (header + "fix sea all dseams record\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-cutoff": (header + "fix sea all dseams record rec.txt cutoff\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-depth": (header + "fix sea all dseams record rec.txt depth\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-host": (header + "fix sea all dseams record rec.txt host\nrun 1\n", 1, "Illegal fix dseams command"),
        "bare-neigh": (header + "fix sea all dseams record rec.txt neigh\nrun 1\n", 1, "Illegal fix dseams command"),
        "word-every": (header + "fix sea all dseams record rec.txt every abc\nrun 1\n", 1, "fix dseams: every must be a positive integer"),
        "word-cutoff": (header + "fix sea all dseams record rec.txt cutoff abc\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
        "trail-cutoff": (header + "fix sea all dseams record rec.txt cutoff 2x\nrun 1\n", 1, "fix dseams: cutoff must be positive"),
    }
    cases = {
        "shift": (
            "units real\natom_style atomic\nboundary p p p\n"
            "lattice diamond 6.4\nregion box block 0 2 0 2 0 2\n"
            "create_box 1 box\ncreate_atoms 1 box\nmass 1 18.0\n"
            "pair_style zero 4.0\npair_coeff * *\nneighbor 2.0 bin\n"
            "timestep 1.0\nvelocity all set 0 0 0\nplugin load %s\n"
            "fix mov all move linear 5.0 0.0 0.0\n"
            "fix sea all dseams record rec.txt every 1 cutoff 3.5\nrun 2\n"
        ) % plugin,
        "live": header
        + "fix hold all setforce 0 0 0\n"
        + "fix sea all dseams dump dump.lammpstrj record rec.txt every 1 cutoff 3.2 depth 4 host yes\n"
        + "run 2\n",
        "move": header
        + "velocity all create 30.0 7\nfix int all nve\n"
        + "fix sea oxy dseams record rec.txt every 1 cutoff 3.5\n"
        + "group oxy type 1\n"
        + "run 2\n",
    }
    # group must exist before the fix that names it
    cases["move"] = header + textwrap.dedent(
        """\
        group oxy type 1
        velocity all set 0 0 0
        fix hold all setforce 0.0 0.0 0.0
        fix sea oxy dseams record rec.txt every 1 cutoff 3.5
        run 1
        displace_atoms all move 1.0 0.0 0.0
        run 1
        """
    )
    for name, (text, code, message) in rejects.items():
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in").write_text(text)
        expect(
            [mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"],
            folder, env, code, message,
        )
        print("lammps", name, "exit", code, flush=True)
    for name, text in cases.items():
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in").write_text(text)
        run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], folder, env)
        print("lammps", name, "exit 0", flush=True)
    both = work / "ranks"
    both.mkdir(parents=True, exist_ok=True)
    (both / "in").write_text(
        header
        + "fix sea all dseams record rec.txt every 2 cutoff 3.5\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "2", lmp, "-in", "in", "-log", "log.lammps"], both, env)
    print("lammps ranks ok", flush=True)
    still = work / "nohost"
    still.mkdir(parents=True, exist_ok=True)
    (still / "in").write_text(
        header + "fix sea all dseams record rec.txt neigh host.txt every 1 cutoff 3.5 host no\nrun 2\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], still, env)
    print("lammps nohost ok", flush=True)
    empty = work / "empty-group"
    empty.mkdir(parents=True, exist_ok=True)
    (empty / "in").write_text(
        header
        + "group none empty\n"
        + "fix sea none dseams record rec.txt every 1 cutoff 3.5\n"
        + "run 1\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], empty, env)
    print("lammps empty-group exit 0", flush=True)
    twice = work / "twice"
    twice.mkdir(parents=True, exist_ok=True)
    (twice / "in").write_text(
        header
        + "fix hold all setforce 0.0 0.0 0.0\n"
        + "fix sea all dseams record rec.txt every 1 cutoff 3.5\n"
        + "run 1\nrun 1\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], twice, env)
    print("lammps twice exit 0", flush=True)
    moving = work / "moving"
    moving.mkdir(parents=True, exist_ok=True)
    (moving / "in").write_text(
        header
        + "fix sea all dseams record rec.txt every 1 cutoff 3.5\n"
        + "run 0\n"
        + "displace_atoms all move 2.0 0.0 0.0\n"
        + "run 1\n"
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], moving, env)
    print("lammps moving exit 0", flush=True)
    subset_case(work, lmp, mpirun, plugin, env)
    print("lammps subset exit 0", flush=True)
    low, high, _lengths, expected = odd_box()
    cutoff = odd_cutoff()
    boxed = work / "box-digits"
    boxed.mkdir(parents=True, exist_ok=True)
    (boxed / "data").write_text(odd_data(low, high))
    (boxed / "in").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load %s\nread_data data\nmass 1 1.0\n"
        "pair_style zero 20.0\npair_coeff * *\nneighbor 2.0 bin\n"
        "fix sea all dseams record live.txt neigh host.txt every 1 cutoff %.17g\n"
        "run 0\n" % (plugin, cutoff)
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], boxed, env)
    check_dump(boxed / "host.txt", expected, cutoff)
    print("lammps box-digits exit 0", flush=True)
    intree = need("DSEAMS_LMP_INTREE")
    native = work / "box-intree"
    native.mkdir(parents=True, exist_ok=True)
    (native / "data").write_text(odd_data(low, high))
    (native / "in").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "read_data data\nmass 1 1.0\n"
        "pair_style zero 20.0\npair_coeff * *\nneighbor 2.0 bin\n"
        "fix sea all dseams record live.txt neigh host.txt every 1 cutoff %.17g\n"
        "run 0\n" % cutoff
    )
    run([mpirun, "--oversubscribe", "-np", "1", intree, "-in", "in", "-log", "log.lammps"], native, env)
    check_dump(native / "host.txt", expected, cutoff)
    print("lammps intree box-digits exit 0", flush=True)
    must_raise(
        lambda: run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "no-such-input"], work, env),
        "exit ",
    )
    bad = work / "bad-every"
    cmd = [mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"]
    must_raise(
        lambda: expect(cmd, bad, env, 0, "fix dseams: every must be a positive integer"),
        "want 0",
    )
    must_raise(
        lambda: expect(cmd, bad, env, 1, "not-the-message"),
        "missing not-the-message",
    )



def require_same(label, left, right):
    if left != right:
        raise SystemExit("%s\n%s\n%s" % (label, left, right))


def odd_box():
    low = (0.12345678901234568, -0.98765432109876543, 1.1111111111111112)
    high = (12.345678901234567, 9.8765432109876543, 8.1415926535897931)
    lengths = tuple(high[i] - low[i] for i in range(3))
    return low, high, lengths, list(low) + list(lengths) + [0.0, 0.0, 0.0]


def odd_cutoff():
    return 3.141592653589793


def odd_data(low, high):
    return (
        "LAMMPS data file\n\n1 atoms\n1 atom types\n\n"
        "%.17g %.17g xlo xhi\n%.17g %.17g ylo yhi\n%.17g %.17g zlo zhi\n\n"
        "Atoms\n\n1 1 %.17g %.17g %.17g\n"
        % (low[0], high[0], low[1], high[1], low[2], high[2], low[0], low[1], low[2])
    )


def require_numbers(kind, line, expected):
    got = [float(token) for token in line.split()[1:]]
    if got != expected:
        raise SystemExit("%s %s" % (kind, line))


def header_line(path, prefix):
    for line in Path(path).read_text().splitlines():
        if line.startswith(prefix):
            return line
    raise SystemExit("missing %s" % prefix)


def expect_exit(fn, code, message):
    try:
        fn()
    except SystemExit as err:
        status = err.code if isinstance(err.code, int) else 1
        text = err.code if isinstance(err.code, str) else ""
        if status != code or message not in text:
            raise SystemExit("exit %s want %d %s" % (err.code, code, message))
        return
    raise SystemExit("missing exit %d %s" % (code, message))


def require_text(err, needle):
    text = err.code if isinstance(err.code, str) else str(err.code)
    if needle not in text:
        raise SystemExit("missing %s" % needle)


def raise_exit(value):
    raise SystemExit(value)


def require_raised(fn, needle):
    try:
        fn()
    except SystemExit as err:
        require_text(err, needle)
    else:
        raise SystemExit("no exit")


def demand(caught, label):
    if not caught:
        raise SystemExit(label)


def drill_exit_checks():
    require_raised(lambda: require_text(SystemExit(3), "nope"), "missing")
    try:
        demand(False, "broken variant was not caught")
    except SystemExit as err:
        require_text(err, "broken variant")
    caught = False
    try:
        require_raised(lambda: None, "unused")
    except SystemExit as err:
        require_text(err, "no exit")
        caught = True
    demand(caught, "broken variant was not caught")
    caught = False
    try:
        expect_exit(lambda: raise_exit("box ok"), 1, "nope")
    except SystemExit as err:
        require_text(err, "want")
        caught = True
    demand(caught, "wrong message was accepted")
    caught = False
    try:
        expect_exit(lambda: None, 1, "box ")
    except SystemExit as err:
        require_text(err, "missing exit")
        caught = True
    demand(caught, "a command that stays up was accepted")
    caught = False
    try:
        expect_exit(lambda: raise_exit(0), 1, "box ")
    except SystemExit as err:
        require_text(err, "want")
        caught = True
    demand(caught, "wrong status was accepted")
    print("broken variant caught", flush=True)


def check_dump(path, box, cutoff):
    require_numbers("box", header_line(path, "box "), box)
    require_numbers("cutoff", header_line(path, "cutoff "), [cutoff])
    expect_exit(lambda: require_numbers("box", "box 0 0 0 1 1 1 0 0 0", box), 1, "box ")
    expect_exit(lambda: require_numbers("cutoff", "cutoff 1", [cutoff]), 1, "cutoff ")
    blank = Path(path).with_name("no-header.txt")
    blank.write_text("step 0\n")
    expect_exit(lambda: header_line(blank, "box "), 1, "missing box")
    drill_exit_checks()


def census_line(text, key):
    for line in text.splitlines():
        if line.split()[:1] == [key]:
            return line
    return ""


def require_rings_prefix(text):
    line = census_line(text, "rings")
    if not line.startswith("rings "):
        raise SystemExit("subset rings missing")
    return line


def type_dump(src, dest, type_id):
    lines = Path(src).read_text().splitlines()
    out = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            # rewritten after the atom block is known
            out.append(lines[i])
            out.append(lines[i + 1])
            count_at = len(out) - 1
            i += 2
            continue
        if lines[i].startswith("ITEM: ATOMS"):
            header = lines[i].split()
            type_col = header.index("type") - 2
            out.append(lines[i])
            i += 1
            kept = []
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                cols = lines[i].split()
                if int(cols[type_col]) == type_id:
                    kept.append(lines[i])
                i += 1
            out[count_at] = str(len(kept))
            out.extend(kept)
            continue
        out.append(lines[i])
        i += 1
    Path(dest).write_text("\n".join(out) + "\n")


def cli_census(seams, dump, command, extra):
    proc = subprocess.run(
        [seams, command, str(dump), "--cutoff", "3.5", "--jobs", "1", "--frame", "1", "--last", "1", *extra],
        text=True, capture_output=True,
    )
    if proc.returncode != 0:
        raise SystemExit("%s exit %d" % (command, proc.returncode))
    for line in proc.stdout.splitlines():
        if line.startswith("rings ") or line.startswith("nop "):
            return line
    raise SystemExit(command + " wrote no census line")


def subset_case(work, lmp, mpirun, plugin, env):
    seams = need("SEAMS_CLI")
    folder = work / "subset"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "mix.data").write_text(
        "LAMMPS data file for fix dseams subset\n\n"
        "4 atoms\n2 atom types\n\n"
        "0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\n"
        "Atoms\n\n"
        "1 1 0.0 0.0 0.0\n"
        "2 1 2.0 0.0 0.0\n"
        "3 1 1.0 1.8 0.0\n"
        "4 2 0.5 0.2 1.5\n"
    )
    (folder / "in").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load %s\nread_data mix.data\n"
        "mass * 18.0\npair_style zero 4.0\npair_coeff * *\n"
        "neighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\n"
        "group oxy type 1\n"
        "fix sea oxy dseams record live.txt neigh host.txt every 1 cutoff 3.5\n"
        "dump native all custom 1 dump.lammpstrj id type x y z\n"
        "dump_modify native sort id format float %%.17g\n"
        "run 0\n" % plugin
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], folder, env)
    type_dump(folder / "dump.lammpstrj", folder / "type1.lammpstrj", 1)
    live = (folder / "live.txt").read_text()
    live_rings = require_rings_prefix(live)
    host = (folder / "host.txt").read_text()
    require_absent(host, "rings")
    must_raise(lambda: require_rings_prefix(host), "subset rings missing")
    must_raise(lambda: require_absent(live_rings, "rings"), "census line")
    cli_rings = cli_census(seams, folder / "type1.lammpstrj", "rings", [])
    require_same("subset rings", live_rings, cli_rings)
    other = cli_census(seams, folder / "type1.lammpstrj", "rings", ["--cutoff", "1.0"])
    must_raise(lambda: require_same("subset rings", live_rings, other), "subset rings")
    must_raise(
        lambda: cli_census(seams, folder / "type1.lammpstrj", "rings", ["--cutoff", "not-a-number"]),
        "rings exit",
    )
    must_raise(
        lambda: cli_census(seams, folder / "type1.lammpstrj", "--features", []),
        "--features wrote no census line",
    )


def compute_cases(work, lmp, cxx, include, env):
    out = work / "compute"
    out.mkdir(parents=True, exist_ok=True)
    src = ROOT / "extras" / "lammps"
    so = out / "chill_plus_compute.so"
    common = [
        cxx, "-c", "-fPIC", "-std=c++17", "--coverage", "-O0",
        "-fprofile-update=single",
        "-I" + include, "-I" + str(ROOT / "src" / "include" / "internal"),
        "-I" + str(src),
    ]
    run(common + [str(src / "compute_dseams.cpp"), "-o", str(out / "compute_dseams.o")], out, env)
    run(common + [str(ROOT / "tests" / "insitu" / "narg_drive.cpp"), "-o", str(out / "narg_drive.o")], out, env)
    libdir = str(Path(lmp).parent)
    chem = str(Path(env["SEAMS_LIB"]).parents[1] / "chemfiles-prefix" / "lib")
    run(
        [
            cxx, "--coverage",
            str(out / "narg_drive.o"), str(out / "compute_dseams.o"),
            "-L" + libdir, "-llammps",
            "-L" + str(Path(env["SEAMS_LIB"]).parent), "-lyodaLib",
            "-L" + chem, "-lchemfiles",
            "-o", str(out / "narg_drive"),
        ],
        out,
        env,
    )
    drive_env = dict(env)
    drive_env["LD_LIBRARY_PATH"] = libdir + ":" + chem + ":" + drive_env.get("LD_LIBRARY_PATH", "")
    run([str(out / "narg_drive")], out, drive_env)
    print("compute narg guard exit 0", flush=True)
    run(common + [str(src / "plugin_dseams.cpp"), "-o", str(out / "plugin_dseams.o")], out, env)
    run(
        [
            cxx, "-shared", "--coverage",
            str(out / "compute_dseams.o"), str(out / "plugin_dseams.o"),
            "-L" + str(Path(env["SEAMS_LIB"]).parent), "-lyodaLib",
            "-o", str(so),
        ],
        out,
        env,
    )
    water = ROOT / "extras" / "lammps" / "water.data"
    base = (
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load ./chill_plus_compute.so\nread_data water.data\n"
        "mass 1 18.015\npair_style zero 4.0\npair_coeff * *\n"
    )
    dump = "dump d all custom 1 dump.txt id c_ice\n"
    goods = {
        "ok": base + "compute ice all dseams\n" + dump + "run 1\n",
        "graph-cutoff": base + "compute ice all dseams graph cutoff\n" + dump + "run 0\n",
        "mutual": base + "compute ice all dseams graph mutual\n" + dump + "run 0\n",
        "tri": base + "change_box all triclinic\ncompute ice all dseams cutoff 3.2\n" + dump + "run 0\n",
        "labels": base + "compute ice all dseams file labels.txt\n" + dump + "run 1\n",
    }
    rejects = {
        "bad-graph": (base + "compute ice all dseams graph nope\nrun 0\n", "compute dseams: graph is cutoff or mutual"),
        "bad-cut": (base + "compute ice all dseams cutoff 0\nrun 0\n", "compute dseams: cutoff must be positive"),
        "bad-word": (base + "compute ice all dseams nope 1\nrun 0\n", "Illegal compute dseams command"),
        "bare": (base + "compute ice all dseams graph\nrun 0\n", "Illegal compute dseams command"),
        "bare-file": (base + "compute ice all dseams file\nrun 0\n", "Illegal compute dseams command"),
    }
    empty = (
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load ./chill_plus_compute.so\nlattice sc 4.0\n"
        "region box block 0 1 0 1 0 1\ncreate_box 1 box\nmass 1 1\n"
        "pair_style zero 4.0\npair_coeff * *\ncompute ice all dseams\n"
        "dump d all custom 1 dump.txt id c_ice\nrun 0\n"
    )
    for name, text in goods.items():
        folder = out / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "chill_plus_compute.so").write_bytes(so.read_bytes())
        (folder / "water.data").write_bytes(water.read_bytes())
        (folder / "in").write_text(text)
        run([lmp, "-in", "in", "-log", "log.lammps"], folder, env)
        print("compute", name, "exit 0", flush=True)
    for name, (text, message) in rejects.items():
        folder = out / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "chill_plus_compute.so").write_bytes(so.read_bytes())
        (folder / "water.data").write_bytes(water.read_bytes())
        (folder / "in").write_text(text)
        expect([lmp, "-in", "in", "-log", "log.lammps"], folder, env, 1, message)
        print("compute", name, "exit 1", flush=True)
    folder = out / "empty"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "chill_plus_compute.so").write_bytes(so.read_bytes())
    (folder / "in").write_text(empty)
    expect([lmp, "-in", "in", "-log", "log.lammps"], folder, env, 1, "compute dseams: seams_chill_plus_graph failed")
    print("compute empty exit 1", flush=True)
    low, high, _lengths, expected = odd_box()
    cutoff = odd_cutoff()
    boxed = out / "box-digits"
    boxed.mkdir(parents=True, exist_ok=True)
    (boxed / "chill_plus_compute.so").write_bytes(so.read_bytes())
    (boxed / "data").write_text(odd_data(low, high))
    (boxed / "in").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "plugin load ./chill_plus_compute.so\nread_data data\n"
        "mass 1 1.0\npair_style zero 20.0\npair_coeff * *\nneighbor 2.0 bin\n"
        "compute ice all dseams file labels.txt cutoff %.17g\n"
        "dump d all custom 1 dump.txt id c_ice\nrun 0\n" % cutoff
    )
    run([lmp, "-in", "in", "-log", "log.lammps"], boxed, env)
    check_dump(boxed / "labels.txt", expected, cutoff)
    print("compute box-digits exit 0", flush=True)
    mpirun = need("DSEAMS_MPIRUN")
    paired = out / "np2"
    paired.mkdir(parents=True, exist_ok=True)
    (paired / "chill_plus_compute.so").write_bytes(so.read_bytes())
    (paired / "water.data").write_bytes(water.read_bytes())
    (paired / "in").write_text(
        base + "compute ice all dseams file labels.txt\n" + dump + "run 0\n"
    )
    expect(
        [mpirun, "--oversubscribe", "-np", "2", lmp, "-in", "in", "-log", "log.lammps"],
        paired,
        env,
        1,
        "compute dseams: one rank only",
    )
    print("compute np2 exit 1", flush=True)
    rejected = out / "np2-bad-graph"
    rejected.mkdir(parents=True, exist_ok=True)
    (rejected / "chill_plus_compute.so").write_bytes(so.read_bytes())
    (rejected / "water.data").write_bytes(water.read_bytes())
    (rejected / "in").write_text(base + "compute ice all dseams graph nope\nrun 0\n")
    expect(
        [mpirun, "--oversubscribe", "-np", "2", lmp, "-in", "in", "-log", "log.lammps"],
        rejected,
        env,
        1,
        "compute dseams: graph is cutoff or mutual",
    )
    print("compute np2 bad-graph exit 1", flush=True)


def require_absent(text, key):
    if census_line(text, key) != "":
        raise SystemExit("census line")


def must_raise(fn, needle):
    try:
        fn()
    except SystemExit as exc:
        text = exc.code if isinstance(exc.code, str) else str(exc)
        if needle not in text:
            raise SystemExit("raised %r" % text)
        return
    raise SystemExit("did not raise")


def main():
    must_raise(lambda: need("DSEAMS_NOT_SET"), "set DSEAMS_NOT_SET")
    lmp_src = need("DSEAMS_LAMMPS_SRC")
    seams_lib = need("SEAMS_LIB")
    cxx = need("DSEAMS_MPICXX")
    lmp = need("DSEAMS_LMP")
    mpirun = need("DSEAMS_MPIRUN")
    work = Path(os.environ.get("DSEAMS_ADAPTER_WORK", "var/run/dseams-adapter-cov"))
    work.mkdir(parents=True, exist_ok=True)
    build = work / "plugin"
    build.mkdir(parents=True, exist_ok=True)
    (build / "CMakeCache.txt").unlink(missing_ok=True)
    run(
        [
            "cmake", "-S", str(ROOT / "src" / "adapters" / "lammps"), "-B", str(build),
            "-DLAMMPS_SRC=" + lmp_src,
            "-DSEAMS_INCLUDE=" + str(ROOT / "src" / "include" / "internal"),
            "-DSEAMS_LIBRARY=" + seams_lib,
            "-DDSEAMS_LAMMPS_MPI=ON",
            "-DCMAKE_BUILD_TYPE=Debug",
            "-DCMAKE_CXX_COMPILER=" + cxx,
            "-DCMAKE_CXX_FLAGS=--coverage -O0 -fprofile-update=single",
            "-DCMAKE_EXE_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
            "-DCMAKE_SHARED_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
            "-DCMAKE_MODULE_LINKER_FLAGS=-Wl,--allow-shlib-undefined",
        ],
        work,
    )
    run(["cmake", "--build", str(build), "-j", "4"], work)
    plugin = str(build / "libdseams_lammps.so")
    libdir = str(Path(seams_lib).parent)
    lammps_cases(work / "lmp", plugin, lmp, mpirun, libdir + ":" + str(Path(seams_lib).parent.parent / "subprojects" / "readcon-core"))
    compute_cases(work, lmp, cxx, lmp_src, os.environ.copy())
    print(work)


if __name__ == "__main__":
    main()
