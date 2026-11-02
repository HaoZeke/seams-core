#!/usr/bin/env python3
"""One-session timing harness.

Pins the process to one core, attempts to pin the CPU frequency where the
kernel exposes a governor, runs warmup repetitions and then the recorded
repetitions inside one session, and writes raw samples as JSON. Medians and
intervals are not stored; scripts/render_timings.py computes them.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
H2H = Path(os.environ.get("H2H_ROOT", str(ROOT / "h2h")))
CORE = 3
WARMUP = 5
REPS = 30
CUTOFF_A = 3.5
CUT = CUTOFF_A
CUTOFF_NM = CUTOFF_A / 10

RINGS_OPTIONS = """#######################################
       R.I.N.G.S. options file        #
#######################################
 PBC             .true.
 Frac            .false.
 g(r)            .false.
 S(q)            .false.
 S(k)            .false.
 gfft(r)         .false.
 MSD             .false.
 atMSD           .false.
 Bonds           .false.
 Angles          .false.
 Chains          .false.
---- ! chain statistics options ! -----
 Species           0
 AAAA              .false.
 ABAB              .false.
 01                .false.
---------------------------------------
 Rings           .true.
---- ! ring statistics options ! -----
 Species           0
 ABAB              .false.
 Rings0            .false.
 Rings1            .false.
 Rings2            .false.
 Rings3            .false.
 Rings4            .false.
 Prim_Rings        .true.
 Str_Rings         .false.
 BarycRings        .false.
 Prop-1            .false.
 Prop-2            .false.
 Prop-3            .false.
 Prop-4            .false.
 Prop-5            .false.
---------------------------------------
 Vacuum          .false.
#######################################
         ! OUTPUT OPTIONS !           #
#######################################
 Evol            .false.
 Dxout           .false.
! OPTION !       ! DX output options !
 RadOut          .false.
 RingsOut        .false.
 DRngOut         .false.
 VoidsOut        .false.
 TetraOut        .false.
 TrajOut         .false.
! OPTION !       ! Output options !
 Output        h2h.out
#######################################
"""


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="replace")
    except OSError as err:
        return f"<unreadable {path}: {err}>"


def command_version(cmd: list[str]) -> str:
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as err:
        return f"<unavailable: {err}>"
    text = (proc.stdout or proc.stderr).strip()
    return text.splitlines()[0] if text else f"exit {proc.returncode}"


def load_average() -> str:
    return Path("/proc/loadavg").read_text().strip()


def cpuinfo_core(core: int) -> dict:
    blocks = Path("/proc/cpuinfo").read_text().split("\n\n")
    for block in blocks:
        fields = {}
        for line in block.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key.strip()] = value.strip()
        if fields.get("processor") == str(core):
            keep = (
                "processor",
                "vendor_id",
                "cpu family",
                "model",
                "model name",
                "stepping",
                "microcode",
                "cpu MHz",
                "cache size",
                "bogomips",
            )
            return {key: fields.get(key, "") for key in keep}
    return {}


def try_frequency_pin(core: int) -> dict:
    """Record every attempt. This guest has no cpufreq directory; do not write MSRs."""
    attempts = []
    cpu = Path(f"/sys/devices/system/cpu/cpu{core}/cpufreq")
    attempts.append(
        {
            "action": "stat cpufreq directory",
            "path": str(cpu),
            "exists": cpu.is_dir(),
        }
    )
    pstate = Path("/sys/devices/system/cpu/intel_pstate")
    attempts.append(
        {
            "action": "stat intel_pstate",
            "path": str(pstate),
            "exists": pstate.is_dir(),
        }
    )
    msr = Path(f"/dev/cpu/{core}/msr")
    attempts.append(
        {
            "action": "stat msr node; not written",
            "path": str(msr),
            "exists": msr.exists(),
        }
    )
    cpupower = shutil.which("cpupower")
    attempts.append({"action": "which cpupower", "path": cpupower or "", "exists": bool(cpupower)})
    pinned = False
    detail = "no cpufreq sysfs and no cpupower; frequency left at the hypervisor default"
    if cpu.is_dir():
        governor = cpu / "scaling_governor"
        available = cpu / "scaling_available_governors"
        speeds = cpu / "scaling_available_frequencies"
        info = {
            "governor_before": read_text(governor).strip() if governor.exists() else "",
            "available_governors": read_text(available).strip() if available.exists() else "",
            "available_frequencies": read_text(speeds).strip() if speeds.exists() else "",
        }
        try:
            governor.write_text("performance\n")
            info["governor_after"] = read_text(governor).strip()
            pinned = info["governor_after"] == "performance"
            detail = "governor set to performance"
        except OSError as err:
            info["write_error"] = f"{type(err).__name__}: {err}"
            detail = "cpufreq exists but the governor could not be changed"
        attempts.append({"action": "set performance governor", **info})
        if pinned and speeds.exists():
            # Highest listed frequency, if userspace or performance exposes a setter.
            values = read_text(speeds).split()
            setter = cpu / "scaling_setspeed"
            if values and setter.exists():
                try:
                    setter.write_text(values[0] + "\n")
                    attempts.append(
                        {
                            "action": "set scaling_setspeed",
                            "value": values[0],
                            "readback": read_text(setter).strip(),
                        }
                    )
                except OSError as err:
                    attempts.append(
                        {
                            "action": "set scaling_setspeed",
                            "value": values[0],
                            "write_error": f"{type(err).__name__}: {err}",
                        }
                    )
    return {
        "pinned": pinned,
        "detail": detail,
        "attempts": attempts,
        "cpu_mhz_at_record": cpuinfo_core(core).get("cpu MHz", ""),
    }


def blas_lib() -> str:
    """Directory to prepend to LD_LIBRARY_PATH, or empty when unset.

    The baseline driver is linked to FlexiBLAS. That directory depends on
    the machine, so it comes from H2H_BLAS_LIB and is not hard-coded.
    """
    return os.environ.get("H2H_BLAS_LIB", "").strip()


def record_blas_env(path: Path | None = None) -> None:
    """Write the H2H_BLAS_LIB value into raw/env.txt of the record."""
    if path is None:
        raw = Path(os.environ.get("H2H_RAW", str(H2H / "raw")))
        path = raw / "env.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    fb = blas_lib()
    line = "H2H_BLAS_LIB=" + (fb if fb else "unset")
    existing = path.read_text() if path.exists() else ""
    kept = [item for item in existing.splitlines() if not item.startswith("H2H_BLAS_LIB=")]
    kept.append(line)
    path.write_text("\n".join(kept) + "\n")


def child_env() -> dict:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "true",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    fb = blas_lib()
    if fb:
        prev = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = fb if not prev else fb + ":" + prev
    return env


def run_once(cmd: list[str], cwd: Path | None, env: dict) -> tuple[float, subprocess.CompletedProcess]:
    start = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    return (time.perf_counter() - start) * 1e3, proc


def session_shell(session_id: str, role: str, comparison: str, scope: str, frame: str,
                  side: str, code_name: str, command: list[str], cwd: Path | None,
                  repeats: int, warmup: int, prepare, parse, log) -> dict:
    """Warmup then recorded repeats, with no other work between them."""
    env = child_env()
    log.write(f"\n== {session_id} {utc_now()} load {load_average()}\n")
    log.write("command: " + " ".join(command) + "\n")
    started = utc_now()
    load_start = load_average()
    samples = []
    last = None
    failures = 0
    for index in range(warmup + repeats):
        if prepare is not None:
            prepare()
        elapsed, proc = run_once(command, cwd, env)
        last = proc
        if proc.returncode != 0:
            failures += 1
            log.write(f"rep {index} exit {proc.returncode}\n{proc.stderr[-1500:]}\n")
        if index >= warmup:
            samples.append([elapsed])
            log.write(f"rep {index} {elapsed:.6f} ms exit {proc.returncode}\n")
    finished = utc_now()
    parsed = parse(last) if last is not None else {}
    if failures:
        parsed["failed_reps"] = failures
    return {
        "id": session_id,
        "role": role,
        "comparison": comparison,
        "scope": scope,
        "frame": frame,
        "side": side,
        "code_name": code_name,
        "command": command,
        "cwd": str(cwd) if cwd else "",
        "warmup": warmup,
        "columns": ["wall_ms"],
        "samples": samples,
        "counts": parsed,
        "agreement": parsed.get("agreement", ""),
        "load_start": load_start,
        "load_end": load_average(),
        "started": started,
        "finished": finished,
    }


def protocol_driver(binary: str, traj: str, mode: str, session_id: str, log) -> dict:
    command = [
        "taskset",
        "-c",
        str(CORE),
        binary,
        traj,
        "1",
        str(WARMUP),
        str(REPS),
        mode,
    ]
    env = child_env()
    log.write(f"\n== {session_id} {utc_now()} load {load_average()}\n")
    log.write("command: " + " ".join(command) + "\n")
    started = utc_now()
    load_start = load_average()
    elapsed, proc = run_once(command, None, env)
    log.write(proc.stdout)
    if proc.returncode != 0:
        log.write(proc.stderr)
        raise RuntimeError(f"{session_id} failed: {proc.stderr[-500:]}")
    meta = {}
    columns = []
    samples = []
    for line in proc.stdout.splitlines():
        if line.startswith("columns "):
            columns = line.split()[1:]
            continue
        parts = line.split()
        if columns and parts and _is_number(parts[0]) and len(parts) == len(columns):
            samples.append([float(part) for part in parts])
            continue
        if line.startswith("sig "):
            meta.setdefault("signatures", [])
            meta["signatures"].append(line[4:])
            continue
        if " " in line:
            key, value = line.split(" ", 1)
            meta[key] = value
    return {
        "id": session_id,
        "role": "quoted-path" if mode == "path" else "head-to-head",
        "comparison": "mW-cage-path" if mode == "path" else mode,
        "scope": "kernel",
        "frame": traj,
        "side": "dseams",
        "code_name": "d-SEAMS library",
        "command": command,
        "cwd": "",
        "warmup": WARMUP,
        "columns": columns,
        "samples": samples,
        "counts": meta,
        "agreement": "",
        "driver_wall_ms": elapsed,
        "load_start": load_start,
        "load_end": load_average(),
        "started": started,
        "finished": utc_now(),
    }


def read_lammps(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[3])
    box = []
    for line in lines[5:8]:
        left, right = line.split()[:2]
        box.append(float(right) - float(left))
    pos = []
    for line in lines[9 : 9 + n]:
        parts = line.split()
        pos.append((float(parts[2]), float(parts[3]), float(parts[4])))
    return pos, box


def replicate(pos, box, rep: int):
    coords = []
    for ix in range(rep):
        for iy in range(rep):
            for iz in range(rep):
                for x, y, z in pos:
                    coords.append((x + ix * box[0], y + iy * box[1], z + iz * box[2]))
    return coords, [length * rep for length in box]


def write_lammpstrj(path: Path, pos, box) -> None:
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(pos)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0.0 {box[0]:.6f}",
        f"0.0 {box[1]:.6f}",
        f"0.0 {box[2]:.6f}",
        "ITEM: ATOMS id type x y z",
    ]
    for index, (x, y, z) in enumerate(pos, 1):
        lines.append(f"{index} 1 {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")


def write_gro_oxygen(path: Path, pos, box) -> None:
    lines = ["oxygen", str(len(pos))]
    for index, (x, y, z) in enumerate(pos, 1):
        residue = index % 100000
        lines.append(
            f"{residue:5d}{'SOL':<5s}{'O':>5s}{residue:5d}{x / 10:8.3f}{y / 10:8.3f}{z / 10:8.3f}"
        )
    lines.append(f"{box[0] / 10:10.5f}{box[1] / 10:10.5f}{box[2] / 10:10.5f}")
    path.write_text("\n".join(lines) + "\n")


def write_gro_grade(path: Path, pos, box) -> None:
    """GRADE steps through solvent atoms four at a time and treats the first as oxygen."""
    lines = ["grade four-site records", str(4 * len(pos))]
    names = ("OW", "HW1", "HW2", "MW")
    for index, (x, y, z) in enumerate(pos, 1):
        residue = index % 100000
        xyz = f"{x / 10:8.3f}{y / 10:8.3f}{z / 10:8.3f}"
        for name in names:
            lines.append(f"{residue:5d}{'SOL':<5s}{name:>5s}{residue:5d}{xyz}")
    lines.append(f"{box[0] / 10:10.5f}{box[1] / 10:10.5f}{box[2] / 10:10.5f}")
    path.write_text("\n".join(lines) + "\n")


def write_xyz(path: Path, pos, box) -> None:
    lines = [str(len(pos)), f"box {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}"]
    for x, y, z in pos:
        lines.append(f"O {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")


def parse_seams_cages(proc: subprocess.CompletedProcess) -> dict:
    counts = {}
    for line in proc.stdout.splitlines():
        if "signature" in line and "cages" in line:
            parts = line.split()
            # nop N graph cutoff signature 5:12 cages 16 atoms 320
            if "signature" in parts and "cages" in parts:
                sig = parts[parts.index("signature") + 1]
                counts[sig] = int(parts[parts.index("cages") + 1])
    return {"cages": counts, "agreement": cage_agreement(counts)}


def cage_agreement(counts: dict) -> str:
    ordered = []
    for key in ("5:12", "0-12-0", "512", "5¹²"):
        if key in counts:
            ordered.append(f"5^12={counts[key]}")
            break
    for key in ("5:12,6:2", "0-12-2", "51262", "5¹²6²"):
        if key in counts:
            ordered.append(f"5^12 6^2={counts[key]}")
            break
    for key in ("5:12,6:4", "0-12-4", "51264"):
        if key in counts:
            ordered.append(f"5^12 6^4={counts[key]}")
            break
    return ", ".join(ordered)


def parse_trace(proc: subprocess.CompletedProcess, work: Path) -> dict:
    cage = work / "cage.txt"
    counts = {}
    if cage.exists():
        header = []
        data = ""
        for line in cage.read_text().splitlines():
            if line.startswith("#") and "0,12,0" in line:
                header = line.replace("#", " ").split()
            elif line and not line.startswith("#"):
                data = line
        values = data.split()
        # header begins at the type names; the data line begins with the frame index
        # so values[i+1] matches header[i] when the header has no frame column.
        for name, value in zip(header, values[1:]):
            if name in {"0,12,0", "0,12,2", "0,12,4"}:
                counts[{"0,12,0": "5:12", "0,12,2": "5:12,6:2", "0,12,4": "5:12,6:4"}[name]] = int(value)
    internal = ""
    for line in proc.stdout.splitlines():
        if "Total time consumed" in line or line.strip().startswith("Current frame"):
            internal += line.strip() + "; "
    return {"cages": counts, "agreement": cage_agreement(counts), "tool_timer": internal.strip()}


def parse_htr(directory: Path) -> dict:
    counts = {}
    files = list(directory.glob("*/numOfCages/0.txt"))
    text = files[0].read_text() if files else ""
    for line in text.splitlines():
        if ":" not in line:
            continue
        name, value = line.split(":", 1)
        mapped = {"0-12-0": "5:12", "0-12-2": "5:12,6:2", "0-12-4": "5:12,6:4"}.get(name.strip())
        if mapped:
            counts[mapped] = int(value.strip())
    return {"cages": counts, "agreement": cage_agreement(counts)}


def parse_sqq(directory: Path) -> dict:
    counts = {}
    cage = directory / "summary" / "cage.csv"
    if cage.exists():
        rows = cage.read_text().splitlines()
        if len(rows) >= 2:
            header = rows[0].split(",")
            values = rows[1].split(",")
            for name, value in zip(header, values):
                if name in {"5¹²", "5¹²6²", "5¹²6⁴"}:
                    mapped = {"5¹²": "5:12", "5¹²6²": "5:12,6:2", "5¹²6⁴": "5:12,6:4"}[name]
                    counts[mapped] = int(float(value))
    return {"cages": counts, "agreement": cage_agreement(counts)}


def parse_grade(proc: subprocess.CompletedProcess) -> dict:
    counts = {}
    for line in proc.stdout.splitlines():
        if "5¹²" in line and "cage:" in line and "6" not in line.split("cage:")[0]:
            counts["5:12"] = int(line.split("cage:")[1].strip())
        elif "6²5¹²" in line and "cage:" in line:
            counts["5:12,6:2"] = int(line.split("cage:")[1].strip())
        elif "6⁴5¹²" in line and "cage:" in line:
            counts["5:12,6:4"] = int(line.split("cage:")[1].strip())
    return {"cages": counts, "agreement": cage_agreement(counts)}


def parse_rings_census(work: Path) -> dict:
    census = {}
    path = work / "rstat" / "evol-RINGS-5.dat"
    if not path.exists():
        return census
    seen = False
    for line in path.read_text(errors="replace").splitlines():
        if "Total number of rings" in line:
            seen = True
            continue
        parts = line.split()
        if seen and len(parts) >= 2:
            try:
                census[int(float(parts[0]))] = int(float(parts[1]))
            except ValueError:
                continue
    return census


def prepare_rings(work: Path, pos, box) -> None:
    if work.exists():
        shutil.rmtree(work)
    (work / "data").mkdir(parents=True)
    write_xyz(work / "data" / "conf.xyz", pos, box)
    (work / "options").write_text(RINGS_OPTIONS)
    # 200, 500, 25.0, 0.125, 90 and 20 are NUMBER_OF_DELTA_R,
    # NUMBER_OF_Q_POINTS, QMAX, SIGMA_LISS, NUMBER_OF_DELTA_ANG and NDV.
    # They feed the g(r), S(q), angle and void grids. The options file sets
    # g(r), S(q), S(k), angles and voids to .false., so these six numbers
    # are read and have no effect. 3 is TAILLE_REF: the search depth, and
    # the maximum ring size is 2*3 = 6. 10 is TAILLC, the maximum chain
    # length. Chains are .false., so 10 is not used either.
    (work / "input").write_text(
        "#######################################\n"
        "       R.I.N.G.S. input file          #\n"
        "#######################################\n"
        f"h2h-si\n{len(pos)}\n1\nO\n1\n0\n"
        f"{box[0]:.6f} {box[1]:.6f} {box[2]:.6f}\n90.0 90.0 90.0\n"
        "1.0\nANI\nconf.xyz\n200\n500\n25.0\n0.125\n90\n20\n"
        "3\n10\n"
        "#######################################\n"
        f"O O {CUTOFF_A:.3f}\nGrtot {CUTOFF_A:.3f}\n"
        "#######################################\n"
    )


def clean_rings(work: Path) -> None:
    for name in ("rstat", "tmp", "gr", "sq", "cstat", "bonds", "angles"):
        shutil.rmtree(work / name, ignore_errors=True)


def main() -> int:
    out_path = Path(os.environ.get("TIMING_JSON", str(H2H / "out" / "timings.json")))
    log_path = Path(os.environ.get("TIMING_LOG", str(H2H / "out" / "timing-harness.log")))
    work = Path(os.environ.get("H2H_WORK", str(H2H / "work"))) / "timing-frames"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    os.sched_setaffinity(0, {CORE})
    log = log_path.open("w")
    seams = os.environ.get("SEAMS", str(ROOT / "build" / "seams"))
    protocol = os.environ.get("BENCH_PROTOCOL", str(ROOT / "build" / "tests" / "bench_protocol"))
    baseline = os.environ.get("BENCH_BASELINE", str(ROOT / "build" / "tests" / "bench_protocol_base"))
    trace = os.environ.get("TRACE", str(H2H / "trace" / "TRACE.exe"))
    grade = os.environ.get("GRADE", str(H2H / "grade" / "GRADE"))
    htr_jar = os.environ.get("HTR_JAR", str(H2H / "htr" / "HTR-pro"))
    sqq = os.environ.get("SQQ", str(H2H / "sqq" / "bin" / "sqq"))
    rings = os.environ.get("RINGS", str(H2H / "rings" / "bin" / "rings"))

    unit_pos, unit_box = read_lammps(ROOT / "input" / "traj" / "genice_sI.lammpstrj")
    frames = {}
    for rep, label in ((2, "sI 2x2x2, 368"), (4, "sI 4x4x4, 2944"), (8, "sI 8x8x8, 23552")):
        pos, box = replicate(unit_pos, unit_box, rep)
        tag = f"si{rep}"
        write_lammpstrj(work / f"{tag}.lammpstrj", pos, box)
        write_gro_oxygen(work / f"{tag}.gro", pos, box)
        write_gro_grade(work / f"{tag}-grade.gro", pos, box)
        frames[tag] = {"label": label, "pos": pos, "box": box, "n": len(pos)}

    record_blas_env()
    frequency = try_frequency_pin(CORE)
    hardware = {
        "uname": platform.uname()._asdict(),
        "lscpu": subprocess.run(["lscpu"], check=False, capture_output=True, text=True).stdout,
        "core": cpuinfo_core(CORE),
        "meminfo_total": next(
            line for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith("MemTotal")
        ),
        "load_at_start": load_average(),
        "affinity": sorted(os.sched_getaffinity(0)),
    }
    software = {
        "g++": command_version(["g++", "--version"]),
        "gfortran": command_version(["gfortran", "--version"]),
        "meson": command_version(["meson", "--version"]),
        "python": platform.python_version(),
        "java": command_version(["java", "-version"]),
        "sqq": command_version([sqq, "--version"]),
        "index_tree": subprocess.run(
            ["git", "-C", str(ROOT), "write-tree"], check=False, capture_output=True, text=True
        ).stdout.strip(),
        "baseline_commit": "1c8efa382c4ae88ebc68e40a991260c4d57edbd5",
        "baseline_version_string": "2.0.0",
        "seams": seams,
        "protocol_driver": protocol,
        "baseline_driver": baseline,
        "trace": trace,
        "grade": grade,
        "htr_jar": htr_jar,
        "sqq_bin": sqq,
        "rings": rings,
        "H2H_BLAS_LIB": blas_lib() if blas_lib() else "unset",
        "trace_compile": "g++ -O3 -fopenmp -std=c++17 -o TRACE.exe TRACE.cpp",
        "grade_compile": "g++ -O3 -std=c++17 -o GRADE GRADE.cpp MyFunctions.cpp",
        "rings_compile": "the 1.3.5 configure/make default, gfortran -g -O2",
        "htr_note": "public repository publishes HTR-pro, a prebuilt JAR, not Java source",
    }
    sessions = []
    mw = str(ROOT / "input" / "traj" / "mW_cubic.lammpstrj")
    sessions.append(protocol_driver(protocol, mw, "path", "path-current", log))
    sessions.append(protocol_driver_baseline(baseline, mw, log))

    def add_whole(session_id, comparison, frame_label, side, code_name, cmd, cwd, prepare, parse):
        sessions.append(
            session_shell(
                session_id,
                "head-to-head",
                comparison,
                "whole-process",
                frame_label,
                side,
                code_name,
                cmd,
                cwd,
                REPS,
                WARMUP,
                prepare,
                parse,
                log,
            )
        )

    for tag, info in frames.items():
        label = info["label"]
        lmp = str(work / f"{tag}.lammpstrj")
        gro = str(work / f"{tag}.gro")
        grade_gro = str(work / f"{tag}-grade.gro")
        # Library kernel, one process, all repeats.
        kernel = protocol_driver(protocol, lmp, "hydrate", f"hydrate-kernel-{tag}", log)
        kernel["comparison"] = "hydrate cages"
        kernel["frame"] = label
        kernel["scope"] = "kernel"
        kernel["agreement"] = _kernel_agreement(kernel["counts"])
        sessions.append(kernel)

        add_whole(
            f"seams-cli-{tag}",
            "hydrate cages",
            label,
            "dseams",
            "d-SEAMS 2.0 seams cages",
            [
                "taskset",
                "-c",
                str(CORE),
                seams,
                "cages",
                lmp,
                "--graph",
                "cutoff",
                "--signature",
                "512+51262+51264",
                "--cutoff",
                "3.5",
                "--type",
                "1",
            ],
            None,
            None,
            parse_seams_cages,
        )
        trace_dir = work / f"trace-{tag}"
        trace_dir.mkdir(exist_ok=True)
        add_whole(
            f"trace-{tag}",
            "hydrate cages",
            label,
            "trace",
            "TRACE",
            [
                "taskset",
                "-c",
                str(CORE),
                trace,
                "-w",
                gro,
                "-r",
                f"{CUTOFF_NM}",
                "-th",
                "-1",
                "-nt",
                "1",
                "-mr",
                "6",
                "-b",
                "0",
                "-e",
                "0",
            ],
            trace_dir,
            None,
            lambda proc, directory=trace_dir: parse_trace(proc, directory),
        )
        htr_dir = work / f"htr-{tag}"

        def prep_htr(directory=htr_dir, source=lmp, name=f"{tag}.lammpstrj"):
            if directory.exists():
                shutil.rmtree(directory)
            directory.mkdir()
            shutil.copy(source, directory / name)

        prep_htr()
        add_whole(
            f"htr-{tag}",
            "hydrate cages",
            label,
            "htr",
            "HTR beta",
            [
                "taskset",
                "-c",
                str(CORE),
                "java",
                "-jar",
                htr_jar,
                str(htr_dir),
                f"{tag}.lammpstrj",
                "-o",
                "1",
                "-b",
                "12",
            ],
            None,
            prep_htr,
            lambda proc, directory=htr_dir: parse_htr(directory),
        )
        if tag != "si8":
            sqq_dir = work / f"sqq-{tag}"

            def prep_sqq(directory=sqq_dir):
                if directory.exists():
                    shutil.rmtree(directory)

            add_whole(
                f"sqq-{tag}",
                "hydrate cages",
                label,
                "sqq",
                "SQQ 0.5.6 cpp",
                [
                    "taskset",
                    "-c",
                    str(CORE),
                    sqq,
                    "analyze",
                    "-i",
                    gro,
                    "-e",
                    "cpp",
                    "-b",
                    "oo",
                    "-w",
                    "1",
                    "-s",
                    "4,5,6",
                    "--find-half",
                    "off",
                    "--find-quasi",
                    "off",
                    "--find-cluster",
                    "off",
                    "--order-parameter",
                    "none",
                    "--output-type",
                    "summary-csv",
                    "-o",
                    str(sqq_dir),
                ],
                None,
                prep_sqq,
                lambda proc, directory=sqq_dir: parse_sqq(directory),
            )
            grade_dir = work / f"grade-{tag}"
            grade_dir.mkdir(exist_ok=True)
            add_whole(
                f"grade-{tag}",
                "hydrate cages",
                label,
                "grade",
                "GRADE",
                [
                    "taskset",
                    "-c",
                    str(CORE),
                    grade,
                    "-i",
                    grade_gro,
                    "-r",
                    f"{CUTOFF_NM}",
                    "-f4",
                    "no",
                    "-o",
                    str(grade_dir / "out"),
                ],
                grade_dir,
                None,
                parse_grade,
            )

    # Primitive rings against R.I.N.G.S. on the 2x2x2 cell.
    ring_kernel = protocol_driver(protocol, str(work / "si2.lammpstrj"), "rings", "rings-kernel-si2", log)
    ring_kernel["comparison"] = "primitive rings D=6"
    ring_kernel["frame"] = frames["si2"]["label"]
    ring_kernel["scope"] = "kernel"
    ring_kernel["agreement"] = ring_kernel["counts"].get("ring_census", "")
    sessions.append(ring_kernel)
    rings_work = work / "rings-si2"
    prepare_rings(rings_work, frames["si2"]["pos"], frames["si2"]["box"])

    def parse_rings(proc, directory=rings_work):
        census = parse_rings_census(directory)
        text = " ".join(f"{size}:{count}" for size, count in sorted(census.items()))
        return {"ring_census": text, "agreement": text, "stdout_tail": proc.stdout[-500:]}

    add_whole(
        "rings-si2",
        "primitive rings D=6",
        frames["si2"]["label"],
        "rings",
        "R.I.N.G.S. 1.3.5",
        ["taskset", "-c", str(CORE), rings, "input"],
        rings_work,
        lambda: clean_rings(rings_work),
        parse_rings,
    )
    # d-SEAMS whole-process ring census on the same frame.
    add_whole(
        "seams-fingerprint-si2",
        "primitive rings D=6",
        frames["si2"]["label"],
        "dseams",
        "d-SEAMS 2.0 seams fingerprint",
        [
            "taskset",
            "-c",
            str(CORE),
            seams,
            "fingerprint",
            str(work / "si2.lammpstrj"),
            "--cutoff",
            "3.5",
            "--type",
            "1",
        ],
        None,
        None,
        lambda proc: {"stdout": proc.stdout.strip(), "agreement": proc.stdout.strip()},
    )

    payload = {
        "created": utc_now(),
        "protocol": {
            "core": CORE,
            "warmup": WARMUP,
            "repeats": REPS,
            "omp_num_threads": 1,
            "statistic": (
                "median; for an even count, the average of the two central order statistics. "
                "The interval is the narrowest central order-statistic interval whose binomial "
                "coverage of the population median is at least 0.95."
            ),
            "frequency_pin": frequency,
        },
        "hardware": hardware,
        "software": software,
        "sessions": sessions,
        "load_at_end": load_average(),
    }
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    log.write(f"\nwrote {out_path}\n")
    log.close()
    print(out_path)
    return 0


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def _kernel_agreement(counts: dict) -> str:
    signatures = counts.get("signatures") or []
    if signatures:
        return "; ".join(signatures)
    return counts.get("ring_census", "")


def protocol_driver_baseline(binary: str, traj: str, log) -> dict:
    """Baseline binary takes no mode argument; it always runs the path."""
    command = ["taskset", "-c", str(CORE), binary, traj, "1", str(WARMUP), str(REPS)]
    env = child_env()
    log.write(f"\n== path-baseline {utc_now()} load {load_average()}\n")
    log.write("command: " + " ".join(command) + "\n")
    started = utc_now()
    load_start = load_average()
    elapsed, proc = run_once(command, None, env)
    log.write(proc.stdout)
    if proc.returncode != 0:
        log.write(proc.stderr)
        raise RuntimeError(proc.stderr[-500:])
    meta = {}
    columns = []
    samples = []
    for line in proc.stdout.splitlines():
        if line.startswith("columns "):
            columns = line.split()[1:]
            continue
        parts = line.split()
        if columns and parts and _is_number(parts[0]) and len(parts) == len(columns):
            samples.append([float(part) for part in parts])
            continue
        if " " in line:
            key, value = line.split(" ", 1)
            meta[key] = value
    return {
        "id": "path-baseline",
        "role": "quoted-path",
        "comparison": "mW-cage-path",
        "scope": "kernel",
        "frame": traj,
        "side": "baseline",
        "code_name": "baseline 1c8efa38",
        "command": command,
        "cwd": "",
        "warmup": WARMUP,
        "columns": columns,
        "samples": samples,
        "counts": meta,
        "agreement": "",
        "driver_wall_ms": elapsed,
        "load_start": load_start,
        "load_end": load_average(),
        "started": started,
        "finished": utc_now(),
    }


if __name__ == "__main__":
    raise SystemExit(main())
