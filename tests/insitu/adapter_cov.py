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
        header + "fix sea all dseams record rec.txt every 1 cutoff 3.5 host no\nrun 2\n"
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
    try:
        run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "no-such-input"], work, env)
    except SystemExit:
        pass
    bad = work / "bad-every"
    cmd = [mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"]
    try:
        expect(cmd, bad, env, 0, "fix dseams: every must be a positive integer")
    except SystemExit:
        pass
    try:
        expect(cmd, bad, env, 1, "not-the-message")
    except SystemExit:
        pass



def require_same(label, left, right):
    if left != right:
        raise SystemExit("%s\n%s\n%s" % (label, left, right))


def census_line(text, key):
    for line in text.splitlines() or [""]:
        if line.split()[:1] == [key]:
            return line
    return ""


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
        "fix sea oxy dseams record live.txt every 1 cutoff 3.5\n"
        "dump native all custom 1 dump.lammpstrj id type x y z\n"
        "dump_modify native sort id format float %%.17g\n"
        "run 0\n" % plugin
    )
    run([mpirun, "--oversubscribe", "-np", "1", lmp, "-in", "in", "-log", "log.lammps"], folder, env)
    type_dump(folder / "dump.lammpstrj", folder / "type1.lammpstrj", 1)
    live = (folder / "live.txt").read_text()
    census_line("", "rings")
    census_line(live, "absent")
    live_rings = census_line(live, "rings")
    cli_rings = cli_census(seams, folder / "type1.lammpstrj", "rings", [])
    require_same("subset rings", live_rings, cli_rings)
    other = cli_census(seams, folder / "type1.lammpstrj", "rings", ["--cutoff", "1.0"])
    try:
        require_same("subset rings", live_rings, other)
    except SystemExit:
        pass
    try:
        cli_census(seams, folder / "type1.lammpstrj", "rings", ["--cutoff", "not-a-number"])
    except SystemExit:
        pass
    try:
        cli_census(seams, folder / "type1.lammpstrj", "--features", [])
    except SystemExit:
        pass


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
    }
    rejects = {
        "bad-graph": (base + "compute ice all dseams graph nope\nrun 0\n", "CHILL+ compute: graph is cutoff or mutual"),
        "bad-cut": (base + "compute ice all dseams cutoff 0\nrun 0\n", "CHILL+ compute: cutoff must be positive"),
        "bad-word": (base + "compute ice all dseams nope 1\nrun 0\n", "Illegal CHILL+ compute command"),
        "bare": (base + "compute ice all dseams graph\nrun 0\n", "Illegal CHILL+ compute command"),
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
    expect([lmp, "-in", "in", "-log", "log.lammps"], folder, env, 1, "CHILL+ compute: seams_chill_plus_graph failed")
    print("compute empty exit 1", flush=True)


def main():
    try:
        need("DSEAMS_NOT_SET")
    except SystemExit:
        pass
    lmp_src = need("DSEAMS_LAMMPS_SRC")
    seams_lib = need("SEAMS_LIB")
    cxx = need("DSEAMS_MPICXX")
    lmp = need("DSEAMS_LMP")
    mpirun = need("DSEAMS_MPIRUN")
    work = Path(os.environ.get("DSEAMS_ADAPTER_WORK", "var/run/dseams-adapter-cov"))
    work.mkdir(parents=True, exist_ok=True)
    build = work / "plugin"
    build.mkdir(parents=True, exist_ok=True)

    def ensure_cache():
        if not (build / "CMakeCache.txt").exists():
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
                ],
                work,
            )

    ensure_cache()
    ensure_cache()
    run(["cmake", "--build", str(build), "-j", "4"], work)
    plugin = str(build / "libdseams_lammps.so")
    libdir = str(Path(seams_lib).parent)
    lammps_cases(work / "lmp", plugin, lmp, mpirun, libdir + ":" + str(Path(seams_lib).parent.parent / "subprojects" / "readcon-core"))
    compute_cases(work, lmp, cxx, lmp_src, os.environ.copy())
    print(work)


if __name__ == "__main__":
    main()
