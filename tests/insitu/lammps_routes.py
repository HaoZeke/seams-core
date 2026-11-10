"""Head-to-head of the d-SEAMS LAMMPS plugin and the in-tree package.

Same host, same LAMMPS version string, one thread. Six samples after one
warmup. Quartiles are Hyndman-Fan type 7. A route is faster when its third
quartile is below the other route's first quartile.
"""

import json
import os
import socket
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(os.environ.get("DSEAMS_ROUTE_WORK", "var/run/dseams-routes"))


def type7(values, probability):
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise SystemExit("no samples")
    if n == 1:
        return ordered[0]
    position = (n - 1) * probability
    lo = int(position)
    frac = position - lo
    if lo + 1 >= n:
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac


def loop_time(text):
    for line in text.splitlines():
        if line.startswith("Loop time of "):
            return float(line.split()[3])
    raise RuntimeError("LAMMPS reported no loop time\n" + text[-400:])


def section_time(text, name):
    needle = name + "  time (%)"
    for line in text.splitlines():
        if line.startswith(needle) or line.startswith(name + " time (%)"):
            return float(line.split("=")[1].split("(")[0])
    return None


def input_text(steps, every, ncell, record, plugin):
    load = ""
    if plugin:
        load = "plugin load %s\n" % plugin
    return (
        "units real\natom_style atomic\nboundary p p p\n"
        "lattice diamond 6.4\n"
        "region box block 0 %d 0 %d 0 %d\n"
        "create_box 1 box\ncreate_atoms 1 box\n"
        "mass 1 18.0154\npair_style lj/cut 8.0\npair_coeff * * 0.05 2.5\n"
        "neighbor 2.0 bin\ntimestep 2.0\n"
        "velocity all create 30.0 12345\nfix int all nve\n"
        "%sfix ds all dseams record %s every %d cutoff 3.5\n"
        "run %d\n"
    ) % (ncell, ncell, ncell, load, record, every, steps)


def run_case(binary, mpirun, np, steps, every, ncell, plugin, folder, env):
    folder.mkdir(parents=True, exist_ok=True)
    record = folder / "rec.txt"
    (folder / "in").write_text(input_text(steps, every, ncell, record, plugin))
    cmd = []
    if np > 1 or mpirun:
        cmd = [mpirun, "--oversubscribe", "-np", str(np), binary, "-in", "in", "-log", "log"]
    else:
        cmd = [binary, "-in", "in", "-log", "log"]
    proc = subprocess.run(cmd, cwd=folder, env=env, text=True, capture_output=True)
    log = (folder / "log").read_text() if (folder / "log").exists() else ""
    blob = proc.stdout + proc.stderr + log
    if proc.returncode != 0 or "Loop time of " not in blob:
        raise RuntimeError("run failed %s\n%s" % (" ".join(cmd), blob[-800:]))
    return loop_time(blob) / steps, section_time(blob, "Modify"), blob


def plan_rows():
    if os.environ.get("DSEAMS_ROUTE_TINY") == "1":
        return [(1000, 5, 1, 2, 1)]
    rows = []
    for ncell, natoms, base in ((5, 1000, 40), (20, 64000, 8)):
        for every in (1, 10, 100):
            steps = max(base, every)
            for np in (1, 2, 4):
                rows.append((natoms, ncell, every, steps, np))
    return rows


def main():
    route = sys.argv[1]
    binary = sys.argv[2]
    plugin = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != "-" else ""
    mpirun = os.environ.get("DSEAMS_MPIRUN", "")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    plan = plan_rows()
    out = {
        "route": route,
        "binary": binary,
        "plugin": plugin,
        "hostname": socket.gethostname(),
        "cases": {},
    }
    for natoms, ncell, every, steps, np in plan:
        key = "%d/%d/np%d" % (natoms, every, np)
        folder = WORK / route / key.replace("/", "_")
        warmup, _, _ = run_case(binary, mpirun, np, steps, every, ncell, plugin, folder / "warm", env)
        samples = []
        modify = []
        for i in range(6):
            per, mod, _ = run_case(
                binary, mpirun, np, steps, every, ncell, plugin, folder / str(i), env
            )
            samples.append(per)
            if mod is not None:
                modify.append(mod / steps)
        out["cases"][key] = {
            "waters": natoms,
            "every": every,
            "steps": steps,
            "ranks": np,
            "warmup": warmup,
            "samples": samples,
            "modify": modify,
        }
        print(route, key, "median", "%.17g" % type7(samples, 0.5), flush=True)
    dest = WORK / ("%s.json" % route)
    dest.write_text(json.dumps(out))
    print(dest)


def compare(plugin_path, intree_path):
    plugin = json.loads(Path(plugin_path).read_text())
    intree = json.loads(Path(intree_path).read_text())
    print("waters every ranks plugin_q1 plugin_q3 intree_q1 intree_q3 faster")
    real = False
    for key, block in plugin["cases"].items():
        other = intree["cases"][key]
        pq1, pq3 = type7(block["samples"], 0.25), type7(block["samples"], 0.75)
        iq1, iq3 = type7(other["samples"], 0.25), type7(other["samples"], 0.75)
        faster = "none"
        if pq3 < iq1:
            faster = "plugin"
            real = True
        elif iq3 < pq1:
            faster = "intree"
            real = True
        print(
            "%s %s %s %.17g %.17g %.17g %.17g %s"
            % (block["waters"], block["every"], block["ranks"], pq1, pq3, iq1, iq3, faster)
        )
        pm = block.get("modify") or []
        im = other.get("modify") or []
        if pm and im:
            print(
                "  modify_s_per_step plugin %.17g intree %.17g"
                % (type7(pm, 0.5), type7(im, 0.5))
            )
    print("real_difference", real)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        main()
