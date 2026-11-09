"""Cadence, postprocess and OpenMP samples for the v4 kit.

LAMMPS live samples are the plugin route already measured at one rank.
This driver adds the plain loop, the native GROMACS module, the command
line postprocess, and the OpenMP ladder. One warmup is stored apart from
the six samples. Quartiles are computed later.
"""

import json
import os
import socket
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("SEAMS_LIB", str(ROOT / "bbdir-rel" / "src" / "libyodaLib.so"))
os.environ.setdefault("SEAMS_INSITU_EXE", str(ROOT / "bbdir-rel" / "seams-insitu"))
os.environ.setdefault("SEAMS_GRO_IMD_EXE", str(ROOT / "bbdir-rel" / "seams-gromacs-imd"))
os.environ.setdefault("SEAMS_CLI", str(ROOT / "bbdir-rel" / "seams"))
os.environ.setdefault("DSEAMS_LIVE_WORK", str(ROOT / "build" / "insitu-measure"))

sys.path.insert(0, str(Path(__file__).resolve().parent))
import live  # noqa: E402


def need(name):
    value = os.environ.get(name, "")
    if not value:
        raise SystemExit("set %s" % name)
    return value


WORK = Path(os.environ["DSEAMS_LIVE_WORK"])
LMP = need("DSEAMS_LMP")
MPIRUN = need("DSEAMS_MPIRUN")
PLUGIN = need("DSEAMS_PLUGIN")
GMX = need("DSEAMS_GMX")
GMX_MPIRUN = need("DSEAMS_GMX_MPIRUN")
SEAMS = os.environ["SEAMS_CLI"]
LIBS = ":".join(
    [
        os.environ.get("DSEAMS_EXTRA_LIBS", ""),
        str(Path(os.environ["SEAMS_LIB"]).parent),
        str(ROOT / "bbdir-rel" / "subprojects" / "readcon-core"),
        str(ROOT / "bbdir-rel" / "subprojects" / "readcon-db"),
    ]
)


def env():
    merged = os.environ.copy()
    merged["OMP_NUM_THREADS"] = "1"
    merged["LD_LIBRARY_PATH"] = LIBS
    return merged


def loop_time(folder, steps):
    log = (folder / "log").read_text() if (folder / "log").exists() else ""
    if "Loop time of " not in log:
        raise RuntimeError(log[-800:])
    for line in log.splitlines():
        if line.startswith("Loop time of "):
            return float(line.split()[3]) / steps
    raise RuntimeError(log[-400:])


def lammps_live_step(folder, steps, ncell, every):
    folder.mkdir(parents=True, exist_ok=True)
    live.lammps_input(
        folder / "in", steps, folder / "unused.lammpstrj", folder / "record.txt",
        PLUGIN, every=every, ncell=ncell, write_dump=False,
    )
    proc = subprocess.run(
        [MPIRUN, "--oversubscribe", "-np", "1", LMP, "-in", "in", "-log", "log"],
        cwd=folder, env=env(), text=True, capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError((proc.stdout + proc.stderr)[-800:])
    return loop_time(folder, steps)


def lammps_plain_step(folder, steps, ncell):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "in").write_text(
        "units real\natom_style atomic\nboundary p p p\n"
        "lattice diamond 6.4\n"
        "region box block 0 %d 0 %d 0 %d\n"
        "create_box 1 box\ncreate_atoms 1 box\n"
        "mass 1 18.0154\npair_style lj/cut 8.0\npair_coeff * * 0.05 2.5\n"
        "neighbor 2.0 bin\ntimestep 2.0\n"
        "velocity all create 30.0 12345\nfix int all nve\n"
        "run %d\n" % (ncell, ncell, ncell, steps)
    )
    proc = subprocess.run(
        [MPIRUN, "--oversubscribe", "-np", "1", LMP, "-in", "in", "-log", "log"],
        cwd=folder,
        env=env(),
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError((proc.stdout + proc.stderr)[-800:])
    return loop_time(folder, steps)


def postprocess(dump):
    samples = []
    for _ in range(6):
        started = __import__("time").perf_counter()
        for command, extra in (
            ("rings", []),
            ("chill-plus", []),
            ("cages", ["--graph", "cutoff"]),
        ):
            cmd = [SEAMS, command, str(dump), "--cutoff", "3.5", "--type", "1", "--jobs", "1", *extra]
            proc = subprocess.run(cmd, cwd=dump.parent, env=env(), text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr[-400:])
        samples.append(__import__("time").perf_counter() - started)
    return samples


def gromacs_case(folder, natoms, steps, every, live_run):
    folder.mkdir(parents=True, exist_ok=True)
    rows = live.diamond_positions(round((natoms / 8.0) ** (1.0 / 3.0)))
    if len(rows) != natoms:
        raise RuntimeError("lattice %d has %d" % (natoms, len(rows)))
    ncell = round((natoms / 8.0) ** (1.0 / 3.0))
    live.write_gro(folder / "conf.gro", rows, ncell * 0.64)
    live.write_oxygen_top(folder / "topol.top", natoms)
    mdp = (
        "integrator = md\ndt = 0.002\nnsteps = %d\n"
        "nstxout = 0\nnstvout = 0\nnstfout = 0\nnstlog = %d\n"
        "nstlist = 10\ncutoff-scheme = Verlet\n"
        "coulombtype = Cut-off\nrcoulomb = 0.9\nrvdw = 0.9\n"
        "pbc = xyz\nconstraints = none\n"
        "gen-vel = yes\ngen-temp = 30\ngen-seed = 12345\n"
    ) % (steps, steps)
    if live_run:
        mdp += (
            "dseams-active = yes\n"
            "dseams-input-group = System\n"
            "dseams-nstcalc = %d\n"
            "dseams-cutoff = 3.5\n"
            "dseams-output = census.txt\n" % every
        )
    (folder / "md.mdp").write_text(mdp)
    live.run(
        [GMX, "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "topol.tpr", "-maxwarn", "2"],
        folder,
        env=env(),
    )
    proc = live.run(
        [GMX_MPIRUN, "--oversubscribe", "-np", "1", GMX, "mdrun", "-s", "topol.tpr", "-deffnm", "run", "-ntomp", "1", "-noconfout"],
        folder,
        env=env(),
    )
    return live.gmx_wall(proc.stdout + proc.stderr) / steps


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    host = live.host_block()
    host.pop("hostname", None)
    out = {
        "host": host,
        "codes": {},
        "scale": {},
        "gromacs_scale": {},
        "ladder": [],
    }
    out["host"]["lib_sha256"] = live.sha256(os.environ["SEAMS_LIB"])
    out["host"]["lmp_sha256"] = live.sha256(LMP)
    out["host"]["gmx_sha256"] = live.sha256(GMX)
    out["host"]["plugin_sha256"] = live.sha256(PLUGIN)
    out["host"]["cpu_online"] = len(os.sched_getaffinity(0))
    plan = [(5, 1000, 40), (20, 64000, 8)]
    for ncell, natoms, base in plan:
        out["scale"][str(natoms)] = {}
        out["gromacs_scale"][str(natoms)] = {}
        for every in (1, 10, 100):
            steps = max(base, every)
            key = "%d/%d/np1" % (natoms, every)
            slug = key.replace("/", "_")
            plain = []
            live_samples = []
            warm = lammps_plain_step(WORK / "plain" / slug / "warm", steps, ncell)
            warm_live = lammps_live_step(WORK / "live" / slug / "warm", steps, ncell, every)
            for i in range(6):
                plain.append(lammps_plain_step(WORK / "plain" / slug / str(i), steps, ncell))
                live_samples.append(lammps_live_step(WORK / "live" / slug / str(i), steps, ncell, every))
            dump = WORK / "post" / str(natoms) / "native.lammpstrj"
            if not dump.exists():
                dump.parent.mkdir(parents=True, exist_ok=True)
                script = dump.parent / "in"
                live.lammps_plain(script, 1, ncell=ncell, every=1, native=dump)
                live.run([str(live.PREFIX / "bin" / "lmp"), "-in", str(script), "-log", str(dump.parent / "log")], dump.parent)
            post = postprocess(dump)
            out["scale"][str(natoms)][str(every)] = {
                "steps": steps,
                "plain": plain,
                "live": live_samples,
                "warmup_plain": warm,
                "warmup_live": warm_live,
                "post_wall": post,
                "dump_bytes": dump.stat().st_size,
            }
            print("lammps", natoms, every, "done", flush=True)
            gplain = []
            glive = []
            gwarm_p = gromacs_case(WORK / "gmx" / key.replace("/", "_") / "plain-warm", natoms, steps, every, False)
            gwarm_l = gromacs_case(WORK / "gmx" / key.replace("/", "_") / "live-warm", natoms, steps, every, True)
            for i in range(6):
                gplain.append(gromacs_case(WORK / "gmx" / key.replace("/", "_") / ("plain-%d" % i), natoms, steps, every, False))
                glive.append(gromacs_case(WORK / "gmx" / key.replace("/", "_") / ("live-%d" % i), natoms, steps, every, True))
            out["gromacs_scale"][str(natoms)][str(every)] = {
                "steps": steps,
                "plain": gplain,
                "live": glive,
                "warmup_plain": gwarm_p,
                "warmup_live": gwarm_l,
            }
            print("gromacs", natoms, every, "done", flush=True)
            if natoms == 1000 and every == 1:
                out["codes"]["lammps"] = {
                    "steps": steps, "frames": steps, "plain": plain, "live": live_samples,
                }
                out["codes"]["gromacs"] = {
                    "steps": steps, "frames": steps, "plain": gplain, "live": glive,
                }
            dest = WORK / "samples.json"
            dest.write_text(json.dumps(out))
    for ncell, natoms in ((5, 1000), (20, 64000)):
        folder = WORK / "ladder-dump" / str(natoms)
        folder.mkdir(parents=True, exist_ok=True)
        dump = folder / "native.lammpstrj"
        script = folder / "in"
        live.lammps_plain(script, 40, ncell=ncell, every=10, native=dump)
        live.run([str(live.PREFIX / "bin" / "lmp"), "-in", str(script), "-log", str(folder / "log")], folder)
        out["ladder"].append(live.thread_ladder(dump, natoms))
        (WORK / "samples.json").write_text(json.dumps(out))
        print("ladder", natoms, "done", flush=True)
    clock = Path(os.environ["DSEAMS_FRAME_CLOCK"]) if os.environ.get("DSEAMS_FRAME_CLOCK") else None
    dump = WORK / "post" / "1000" / "native.lammpstrj"
    if clock is not None and clock.exists() and dump.exists():
        proc = subprocess.run([str(clock), str(dump)], cwd=ROOT, env=env(), text=True, capture_output=True)
        excluded = None
        for line in (proc.stdout + proc.stderr).splitlines():
            if line.startswith("analysis_s "):
                excluded = float(line.split()[1])
        included = postprocess(dump)
        insitu = out["scale"]["1000"]["1"]
        delta = [insitu["live"][i] - insitu["plain"][i] for i in range(6)]
        ordered = sorted(included)
        out["clock"] = {
            "system": "1000 waters, one frame",
            "insitu_per_frame": sorted(delta)[2] * 0.5 + sorted(delta)[3] * 0.5,
            "cli_included": ordered[2] * 0.5 + ordered[3] * 0.5,
            "cli_excluded": excluded,
            "cli_included_samples": included,
        }
        (WORK / "samples.json").write_text(json.dumps(out))
    print(WORK / "samples.json", flush=True)


if __name__ == "__main__":
    main()
