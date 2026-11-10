#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Rerun both routes from the clean tip and stream one log line at a time."""

import csv
import hashlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "tools"))
import streamlog

SRC_SEAMS = Path(os.environ["SRC_SEAMS"])
SRC_LAMMPS = Path(os.environ["SRC_LAMMPS"])
BUILD = Path(os.environ["BUILD"])
SEAMS = BUILD / "seams" / "seams"
LMP_STOCK = BUILD / "lammps-stock" / "lmp"
LMP_INTREE = BUILD / "lammps-intree" / "lmp"
PLUGIN = BUILD / "plugin" / "libdseams_lammps.so"
WATER = SRC_SEAMS / "extras" / "lammps" / "water.data"
LITERAL = SRC_SEAMS / "extras" / "lammps" / "in.dseams"
STAGE = BUILD / "kit"
RUNS = STAGE / "runs"
WORK = BUILD / "work"
MPIRUN = os.environ["MPIRUN"]


def stamp():
    return streamlog.date_is()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


HEAD_SEAMS = git(SRC_SEAMS, "rev-parse", "HEAD")
TREE_SEAMS = git(SRC_SEAMS, "rev-parse", "HEAD^{tree}")
HEAD_LAMMPS = git(SRC_LAMMPS, "rev-parse", "HEAD")
TREE_LAMMPS = git(SRC_LAMMPS, "rev-parse", "HEAD^{tree}")
PORCELAIN_SEAMS = git(SRC_SEAMS, "status", "--porcelain")
PORCELAIN_LAMMPS = git(SRC_LAMMPS, "status", "--porcelain")
CXX = subprocess.check_output(["g++", "--version"], text=True).splitlines()[0]
SEAMS_SHA = sha256(SEAMS)
YODA_SHA = sha256(BUILD / "seams" / "src" / "libyodaLib.so")
PLUGIN_SHA = sha256(PLUGIN)
STOCK_SHA = sha256(LMP_STOCK)
INTREE_SHA = sha256(LMP_INTREE)


def clean(text):
    return streamlog.clean(text)


def env_base():
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    env["PRTE_MCA_rmaps_default_mapping_policy"] = ":oversubscribe"
    env.pop("DSEAMS_NO_HOST", None)
    return env


def timed(cmd, cwd):
    started = stamp()
    script = "TIMEFORMAT='wall_s %R'; time \"$@\""
    proc = subprocess.run(
        ["bash", "-c", script, "bash", *cmd],
        cwd=cwd,
        env=env_base(),
        text=True,
        capture_output=True,
    )
    ended = stamp()
    return started, ended, proc


def open_run_log(path, route, launch, system, binary_name, command_text):
    log = streamlog.Log(path)
    log.line("start: " + stamp())
    log.line("HEAD_SEAMS: " + HEAD_SEAMS)
    log.line("TREE_SEAMS: " + TREE_SEAMS)
    log.line("HEAD_LAMMPS: " + HEAD_LAMMPS)
    log.line("TREE_LAMMPS: " + TREE_LAMMPS)
    log.line("git_status_seams: empty" if PORCELAIN_SEAMS == "" else "git_status_seams: dirty")
    log.line("git_status_lammps: empty" if PORCELAIN_LAMMPS == "" else "git_status_lammps: dirty")
    log.line("compiler: " + CXX)
    log.line("cxx_flags: -O3 -DNDEBUG")
    log.line("seams_buildtype: release")
    log.line("lammps_build: Release MPI shared PKG_PLUGIN no GPU no KOKKOS")
    log.line("binary: " + binary_name)
    log.line("command: " + command_text)
    log.line("route: " + route)
    log.line("launch: " + launch)
    log.line("system: " + system)
    log.line("omp_num_threads: 1")
    log.line("host_keyword: absent")
    log.line("dseams_no_host: unset")
    log.line("substitution: live, longest prefix first, see PATHS.txt")
    log.line("---")
    return log


def require_fix_input(text):
    if "fix sea all dseams" not in text or "compute ice all dseams" in text:
        raise SystemExit("shipped input is not the fix")
    if "plugin load ./libdseams_lammps.so" not in text:
        raise SystemExit("shipped input does not load the fix plugin")


def require_ok(code, message):
    if code != 0:
        raise SystemExit(message)


def quantile7(values, probability):
    ordered = sorted(values)
    count = len(ordered)
    if count == 1:
        return ordered[0]
    position = 1.0 + (count - 1) * probability
    lower = int(position)
    if lower >= count:
        return ordered[-1]
    weight = position - lower
    return (1.0 - weight) * ordered[lower - 1] + weight * ordered[lower]


def wall_of(stderr):
    for line in stderr.splitlines():
        if line.startswith("wall_s "):
            return line.split()[1]
    return ""


def run_lammps(binary, binary_token, binary_sha, folder, launch, route, system, log_name, in_name="in"):
    folder.mkdir(parents=True, exist_ok=True)
    tracer = "$BUILD/tools/trace_libs.py"
    if launch == "serial":
        cmd = [sys.executable, str(BUILD / "tools" / "trace_libs.py"), "stdbuf", "-oL", "-eL", str(binary), "-in", in_name, "-log", "log.lammps"]
        shown = "python3 " + tracer + " stdbuf -oL -eL " + binary_token + " -in " + in_name + " -log log.lammps"
    else:
        np = launch[2:]
        cmd = [sys.executable, str(BUILD / "tools" / "trace_libs.py"), "stdbuf", "-oL", "-eL", MPIRUN, "--oversubscribe", "-np", np, str(binary), "-in", in_name, "-log", "log.lammps"]
        shown = "python3 " + tracer + " stdbuf -oL -eL $MPIRUN --oversubscribe -np %s %s -in %s -log log.lammps" % (np, binary_token, in_name)
    log = open_run_log(RUNS / log_name, route, launch, system, binary_token, shown)
    bash = "TIMEFORMAT='wall_s %R'; time \"$@\""
    code = log.run(["bash", "-c", bash, "bash", *cmd], folder, env_base())
    text = (RUNS / log_name).read_text(errors="replace")
    live = folder / "live.txt"
    if live.is_file():
        log.line("--- census ---")
        for line in live.read_text(errors="replace").splitlines():
            head = line.split()[:1]
            if head and head[0] in ("step", "rings", "chill", "cage", "nop"):
                log.line("census " + line)
    log.line("---")
    log.line("end: " + stamp())
    log.line("binary_sha256: %s %s" % (binary, streamlog.sha256(binary)))
    wall = wall_of(text)
    log.line("exit_code: %d" % code)
    log.close()
    loop = None
    for line in text.splitlines():
        if line.startswith("Loop time of "):
            loop = float(line.split()[3])
    if code != 0:
        raise SystemExit("%s %s %s exit %d" % (route, system, launch, code))
    return code, loop, wall


def script_small(folder, route, steps):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "water.data").write_bytes(WATER.read_bytes())
    load = ""
    if route == "plugin":
        (folder / "libdseams_lammps.so").write_bytes(PLUGIN.read_bytes())
        load = "plugin load ./libdseams_lammps.so\n"
    text = (
        "units real\n"
        "atom_style atomic\n"
        "boundary p p p\n"
        "read_data water.data\n"
        "mass 1 18.015\n"
        "pair_style zero 4.0\n"
        "pair_coeff * *\n"
        "neighbor 2.0 bin\n"
        "neigh_modify delay 0 every 1 check yes\n"
        "%s"
        "fix sea all dseams record live.txt every 1 cutoff 3.5 depth 8\n"
        "dump native all custom 1 dump.lammpstrj id type x y z\n"
        "dump_modify native sort id format float %%.17g\n"
        "%s"
    ) % (load, steps)
    (folder / "in").write_text(text)


def script_large(folder, route, steps, dump_on):
    folder.mkdir(parents=True, exist_ok=True)
    load = ""
    if route == "plugin":
        (folder / "libdseams_lammps.so").write_bytes(PLUGIN.read_bytes())
        load = "plugin load ./libdseams_lammps.so\n"
    dump = ""
    if dump_on:
        dump = (
            "dump native all custom 1 dump.lammpstrj id type x y z\n"
            "dump_modify native sort id format float %.17g\n"
        )
    fix = ""
    if route != "plain":
        fix = "fix sea all dseams record live.txt every 1 cutoff 3.5 depth 8\n"
    (folder / "in").write_text(
        "units real\n"
        "atom_style atomic\n"
        "boundary p p p\n"
        "lattice diamond 6.4\n"
        "region box block 0 4 0 4 0 4\n"
        "create_box 1 box\n"
        "create_atoms 1 box\n"
        "mass 1 18.0154\n"
        "pair_style lj/cut 8.0\n"
        "pair_coeff * * 0.05 2.5\n"
        "neighbor 2.0 bin\n"
        "neigh_modify delay 0 every 1 check yes\n"
        "timestep 2.0\n"
        "velocity all create 30.0 12345\n"
        "fix int all nve\n"
        "%s%s%s%s" % (load, fix, dump, steps)
    )


def census(text):
    blocks = []
    current = None
    for line in text.splitlines():
        if line.startswith("step "):
            if current is not None:
                blocks.append(current)
            current = {"step": int(line.split()[1]), "lines": {}}
        elif current is None:
            continue
        else:
            key = line.split()[0] if line else ""
            current["lines"][key] = line
    if current is not None:
        blocks.append(current)
    return blocks


def dump_steps(path):
    lines = Path(path).read_text(errors="replace").splitlines()
    steps = []
    for index, line in enumerate(lines):
        if line.startswith("ITEM: TIMESTEP") and index + 1 < len(lines):
            steps.append(int(lines[index + 1]))
    return steps


def cli_lines(command, dump, extra, prefix):
    frames = dump_steps(dump)
    cmd = [
        str(SEAMS),
        command,
        "dump.lammpstrj",
        "--cutoff",
        "3.5",
        "--type",
        "1",
        "--frame",
        "1",
        "--last",
        str(len(frames)),
        "--jobs",
        "1",
        *extra,
    ]
    proc = subprocess.run(cmd, cwd=dump.parent, text=True, capture_output=True, env=env_base())
    if proc.returncode != 0:
        raise SystemExit("%s exit %d" % (command, proc.returncode))
    kept = [line for line in proc.stdout.splitlines() if line.startswith(prefix)]
    return kept, proc.stdout


def chill_map(line):
    out = {}
    for item in line.split()[1:]:
        if ":" not in item:
            continue
        name, value = item.split(":", 1)
        out[name] = int(value)
    return out


def cli_pairs(line, start):
    words = line.split()
    return {words[i]: int(words[i + 1]) for i in range(start, len(words), 2)}


def compare_live_cli(system, launch, route, folder, rows):
    live = census((folder / "live.txt").read_text())
    steps = dump_steps(folder / "dump.lammpstrj")
    raw = STAGE / "raw" / system / launch / route
    raw.mkdir(parents=True, exist_ok=True)
    (raw / "live.txt").write_bytes((folder / "live.txt").read_bytes())
    rings, rings_out = cli_lines("rings", folder / "dump.lammpstrj", [], "rings ")
    chills, chill_out = cli_lines("chill-plus", folder / "dump.lammpstrj", ["--graph", "cutoff"], "nop ")
    cages, cages_out = cli_lines("cages", folder / "dump.lammpstrj", ["--graph", "cutoff"], "nop ")
    (raw / "cli-rings.txt").write_text(rings_out)
    (raw / "cli-chill.txt").write_text(chill_out)
    (raw / "cli-cages.txt").write_text(cages_out)
    mismatches = 0
    live_steps = [block["step"] for block in live]
    if live_steps != steps or len(rings) != len(live) or len(chills) != len(live) or len(cages) != len(live):
        mismatches += 1
        rows.append(row(system, launch, route, "", "frame_list", " ".join(map(str, live_steps)), " ".join(map(str, steps)), 1))
    for index, block in enumerate(live):
        if index >= len(rings):
            break
        step = block["step"]
        ring_bad = 0 if block["lines"].get("rings", "") == rings[index] else 1
        mismatches += ring_bad
        rows.append(row(system, launch, route, step, "rings", block["lines"].get("rings", ""), rings[index], ring_bad))
        live_chill = chill_map(block["lines"].get("chill", "chill"))
        cli_chill = cli_pairs(chills[index], 2)
        chill_bad = 0
        for name, count in cli_chill.items():
            if live_chill.get(name, 0) != count:
                chill_bad = 1
        for name, count in live_chill.items():
            if count and cli_chill.get(name, 0) != count:
                chill_bad = 1
        mismatches += chill_bad
        rows.append(row(system, launch, route, step, "chill", block["lines"].get("chill", ""), chills[index], chill_bad))
        cage_fields = {}
        for item in block["lines"].get("cage", "").split()[1:]:
            if ":" in item:
                name, value = item.split(":", 1)
                cage_fields[name] = int(value)
        words = cages[index].split()
        cli_cage = {words[i]: int(words[i + 1]) for i in range(4, len(words), 2)}
        mapping = {"hc": "hexagonal", "ddc": "cubic", "both": "both", "water": "water"}
        cage_bad = 0
        for left, right in mapping.items():
            if cage_fields.get(left) != cli_cage.get(right):
                cage_bad = 1
        mismatches += cage_bad
        rows.append(row(system, launch, route, step, "cage", block["lines"].get("cage", ""), cages[index], cage_bad))
    return {
        "frames": len(live),
        "mismatches": mismatches,
        "live_steps": live_steps,
        "dump_steps": steps,
        "duplicate_live": len(live_steps) != len(set(live_steps)),
        "duplicate_dump": len(steps) != len(set(steps)),
        "blocks": live,
    }


def row(system, launch, route, step, field, live, other, mismatch):
    return {
        "comparison": "live_vs_cli",
        "system": system,
        "launch": launch,
        "route": route,
        "step": step,
        "field": field,
        "live": live,
        "other": other,
        "mismatch": mismatch,
    }


def compare_routes(system, launch, left_blocks, right_blocks, rows):
    mismatches = 0
    fields = ("rings", "chill", "cage", "n", "n_rings", "labels", "cages")
    count = max(len(left_blocks), len(right_blocks))
    for index in range(count):
        left = left_blocks[index] if index < len(left_blocks) else None
        right = right_blocks[index] if index < len(right_blocks) else None
        step = "" if left is None else left["step"]
        if right is not None and step == "":
            step = right["step"]
        if left is None or right is None or left["step"] != right["step"]:
            mismatches += 1
            rows.append(
                {
                    "comparison": "plugin_vs_intree",
                    "system": system,
                    "launch": launch,
                    "route": "both",
                    "step": step,
                    "field": "step",
                    "live": "" if left is None else left["step"],
                    "other": "" if right is None else right["step"],
                    "mismatch": 1,
                }
            )
            continue
        for field in fields:
            a = left["lines"].get(field, "")
            b = right["lines"].get(field, "")
            bad = 0 if a == b else 1
            mismatches += bad
            rows.append(
                {
                    "comparison": "plugin_vs_intree",
                    "system": system,
                    "launch": launch,
                    "route": "both",
                    "step": left["step"],
                    "field": field,
                    "live": a,
                    "other": b,
                    "mismatch": bad,
                }
            )
    return mismatches


def binary_for(route):
    if route == "intree":
        return LMP_INTREE, "$BUILD/lammps-intree/lmp", INTREE_SHA
    return LMP_STOCK, "$BUILD/lammps-stock/lmp", STOCK_SHA


def main():
    RUNS.mkdir(parents=True, exist_ok=True)
    if PORCELAIN_SEAMS or PORCELAIN_LAMMPS:
        raise SystemExit("checkout is not empty")
    launches = ("serial", "np1", "np2", "np4")
    systems = {
        "small-run0": ("small", "run 0\n"),
        "small-run4": ("small", "run 4\n"),
        "small-twice": ("small", "run 2\nrun 2\n"),
        "large-run10": ("large", "run 10\n"),
        "large-twice": ("large", "run 5\nrun 5\n"),
    }
    summary = []
    cli_rows = []
    route_rows = []
    saved = {}
    for system, (kind, steps) in systems.items():
        for launch in launches:
            for route in ("plugin", "intree"):
                binary, token, digest = binary_for(route)
                folder = WORK / route / system / launch
                if kind == "small":
                    script_small(folder, route, steps)
                else:
                    script_large(folder, route, steps, True)
                code, loop, wall = run_lammps(
                    binary, token, digest, folder, launch, route, system,
                    "%s-%s-%s.log" % (route, system, launch),
                )
                entry = {
                    "system": system,
                    "launch": launch,
                    "route": route,
                    "exit_code": code,
                    "wall_s": wall,
                    "loop_s": loop,
                }
                if not (folder / "live.txt").exists():
                    raise SystemExit("%s %s %s wrote no live record" % (route, system, launch))
                compared = compare_live_cli(system, launch, route, folder, cli_rows)
                entry.update(
                    frames=compared["frames"],
                    mismatches=compared["mismatches"],
                    duplicate_live=compared["duplicate_live"],
                    duplicate_dump=compared["duplicate_dump"],
                    live_steps=compared["live_steps"],
                )
                summary.append(entry)
                saved[(system, launch, route)] = compared
                print(route, system, launch, compared["frames"], compared["mismatches"], compared["duplicate_live"], flush=True)
    route_summary = []
    for system in systems:
        for launch in launches:
            left = saved.get((system, launch, "plugin"))
            right = saved.get((system, launch, "intree"))
            if left is None or right is None:
                route_summary.append({"system": system, "launch": launch, "frames_plugin": 0, "frames_intree": 0, "mismatches": ""})
                continue
            mismatches = compare_routes(system, launch, left["blocks"], right["blocks"], route_rows)
            route_summary.append(
                {
                    "system": system,
                    "launch": launch,
                    "frames_plugin": left["frames"],
                    "frames_intree": right["frames"],
                    "mismatches": mismatches,
                }
            )
            print("routes", system, launch, mismatches, flush=True)

    shipped = LITERAL.read_text()
    require_fix_input(shipped)
    for route, binary, token, digest, plugin_bytes in (
        ("plugin", LMP_STOCK, "$BUILD/lammps-stock/lmp", STOCK_SHA, PLUGIN.read_bytes()),
        ("intree", LMP_INTREE, "$BUILD/lammps-intree/lmp", INTREE_SHA, None),
    ):
        folder = WORK / "literal" / route
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "in.dseams").write_text(shipped)
        (folder / "water.data").write_bytes(WATER.read_bytes())
        so = folder / "libdseams_lammps.so"
        if plugin_bytes is None:
            if so.exists():
                so.unlink()
        else:
            so.write_bytes(plugin_bytes)
        code, loop, wall = run_lammps(
            binary, token, digest, folder, "serial", route, "literal-in.dseams",
            "%s-literal-in.dseams-serial.log" % route, in_name="in.dseams",
        )
        require_ok(code, "literal %s exit %d" % (route, code))
        compared = compare_live_cli("literal-in.dseams", "serial", route, folder, cli_rows)
        if compared["mismatches"] != 0 or compared["live_steps"] != [0]:
            raise SystemExit("literal %s mismatch %s" % (route, compared))
        summary.append({
            "system": "literal-in.dseams",
            "launch": "serial",
            "route": route,
            "exit_code": code,
            "wall_s": wall,
            "loop_s": loop,
            "frames": compared["frames"],
            "mismatches": compared["mismatches"],
            "duplicate_live": compared["duplicate_live"],
            "duplicate_dump": compared["duplicate_dump"],
            "live_steps": compared["live_steps"],
        })
        print("literal", route, code, compared["mismatches"], flush=True)
        saved[("literal-in.dseams", "serial", route)] = compared
    left = saved[("literal-in.dseams", "serial", "plugin")]
    right = saved[("literal-in.dseams", "serial", "intree")]
    literal_routes = compare_routes("literal-in.dseams", "serial", left["blocks"], right["blocks"], route_rows)
    if literal_routes != 0:
        raise SystemExit("literal routes mismatch %d" % literal_routes)
    route_summary.append({
        "system": "literal-in.dseams",
        "launch": "serial",
        "frames_plugin": left["frames"],
        "frames_intree": right["frames"],
        "mismatches": literal_routes,
    })

    fields = ["comparison", "system", "launch", "route", "step", "field", "live", "other", "mismatch"]
    eq = STAGE / "equality"
    eq.mkdir(parents=True, exist_ok=True)
    with (eq / "live_vs_cli.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cli_rows)
    with (eq / "plugin_vs_intree.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(route_rows)
    with (eq / "counts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["system", "launch", "route", "exit_code", "frames", "mismatches", "duplicate_live", "duplicate_dump", "live_steps", "wall_s", "loop_s"],
        )
        writer.writeheader()
        writer.writerows(summary)
    with (eq / "route_counts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["system", "launch", "frames_plugin", "frames_intree", "mismatches"])
        writer.writeheader()
        writer.writerows(route_summary)

    lscpu = subprocess.check_output(["lscpu"], text=True)
    print("lscpu_begin", flush=True)
    print(lscpu, end="" if lscpu.endswith("\n") else "\n", flush=True)
    print("lscpu_end", flush=True)
    arms = ["plain", "plugin", "intree"]
    for route in arms:
        binary, token, digest = binary_for(route)
        folder = WORK / "timing" / "warmup" / route
        script_large(folder, route, "run 400\n", False)
        code, loop, wall = run_lammps(
            binary, token, digest, folder, "serial", route, "timing-warmup",
            "timing-warmup-%s-serial.log" % route,
        )
        print("warmup", route, format(loop, ".17g"), wall, code, flush=True)
    rng = random.Random(20261109)
    samples = []
    measured = {route: [] for route in arms}
    for index in range(1, 7):
        order = arms[:]
        rng.shuffle(order)
        sample = {"i": index, "order": order}
        for route in order:
            binary, token, digest = binary_for(route)
            folder = WORK / "timing" / ("%02d" % index) / route
            script_large(folder, route, "run 400\n", False)
            code, loop, wall = run_lammps(
                binary, token, digest, folder, "serial", route, "timing-ncell4-run400",
                "timing-%02d-%s-serial.log" % (index, route),
            )
            sample[route] = {"exit_code": code, "loop_s": format(loop, ".17g"), "wall_s": wall}
            measured[route].append(loop)
            print("time", index, route, format(loop, ".17g"), wall, code, " ".join(order), flush=True)
        samples.append(sample)
    summary_times = {}
    for route, values in measured.items():
        median = quantile7(values, 0.5)
        iqr = quantile7(values, 0.75) - quantile7(values, 0.25)
        summary_times[route] = {
            "n": len(values),
            "median_loop_s": format(median, ".17g"),
            "iqr_type7_loop_s": format(iqr, ".17g"),
        }
        print("summary", route, format(median, ".17g"), format(iqr, ".17g"), flush=True)
    (STAGE / "timing").mkdir(parents=True, exist_ok=True)
    (STAGE / "timing" / "lscpu.txt").write_text(lscpu)
    (STAGE / "timing" / "samples.json").write_text(
        json.dumps(
            {
                "HEAD_SEAMS": HEAD_SEAMS,
                "TREE_SEAMS": TREE_SEAMS,
                "HEAD_LAMMPS": HEAD_LAMMPS,
                "TREE_LAMMPS": TREE_LAMMPS,
                "n": 6,
                "launch": "serial",
                "omp_num_threads": 1,
                "warmup": "one discarded run of each arm, same input, before the measured replicates",
                "arm_order": "random.Random(20261109) shuffle of plain, plugin, intree inside each replicate",
                "quantile": "type 7",
                "input": {
                    "ncell": 4,
                    "atoms": 512,
                    "steps": 400,
                    "pair": "lj/cut",
                    "integrate": "nve",
                    "dseams": "every 1 cutoff 3.5 depth 8",
                    "dump": False,
                    "host_keyword": "absent",
                },
                "note": "wall_s is the wall_s line of the matching timed log; loop_s is printed with .17g",
                "summary": summary_times,
                "samples": samples,
            },
            indent=2,
        )
        + "\n"
    )
    print("campaign done", flush=True)


if __name__ == "__main__":
    main()
