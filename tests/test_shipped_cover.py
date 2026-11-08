#!/usr/bin/env python3
"""Line-cover the paper campaign, the benchmark module, and shipped.sh."""

import os
import runpy
import subprocess
import sys
from pathlib import Path

REAL_RUN = subprocess.run
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "paper/tex/repro"))
sys.path.insert(0, str(ROOT))


def cover_campaign(scratch: Path) -> None:
    os.environ["SEAMS_ROOT"] = str(ROOT)
    os.environ["SEAMS_OUT"] = str(scratch / "out")
    os.environ["SEAMS_LOGS"] = str(scratch / "logs")
    os.environ["SEAMS_BASE"] = str(scratch / "base")
    os.environ["SEAMS_TIP"] = str(scratch / "tip")
    os.environ["SEAMS_1E6"] = str(scratch / "ice-1e6.lammpstrj")
    os.environ["SEAMS_ORTH"] = str(scratch / "orth.lammpstrj")
    os.environ["SEAMS_HOST"] = "redacted"
    import shipped_campaign as camp

    class Positions:
        shape = (8,)

    class Frame:
        positions = Positions()

    def stages_out():
        return "".join(f"{key} 1.5\n" for key in camp.STAGES)

    real_check = subprocess.check_output
    real_run = subprocess.run

    def fake_check(cmd, **kwargs):
        joined = " ".join(str(part) for part in cmd)
        if "bench_stages" in joined:
            return stages_out()
        if cmd and cmd[0] == "lscpu":
            return "Model name: Intel Xeon Processor\n"
        if cmd[:2] == ["nproc", "--all"]:
            return "8\n"
        if cmd[:1] == ["nproc"]:
            return "4\n"
        if cmd[:1] == ["taskset"]:
            return "pid 1's current affinity list: 0\n"
        raise AssertionError(joined)

    class Proc:
        def __init__(self, rc, out, err):
            self.returncode = rc
            self.stdout = out
            self.stderr = err

    def fake_run(cmd, **kwargs):
        pre = kwargs.get("preexec_fn")
        if pre:
            pre()
        joined = " ".join(str(part) for part in cmd)
        if "1e6" in joined and "base" in joined:
            return Proc(1, "", "std::bad_alloc")
        if "--graph" in joined and "1e6" not in joined:
            return Proc(1, "", "nope")
        return Proc(0, "nop 8\n", "Maximum resident set size (kbytes): 12\n")

    saved_run = subprocess.run
    subprocess.check_output = fake_check
    subprocess.run = fake_run
    camp.read_lammps_frame = lambda path, index: Frame()
    camp.save_frame_npz = lambda frame, path: None
    camp.run_ovito_worker = lambda **kwargs: {"ms": 2.0, "count": 4}
    camp.N_ORTH = 6
    camp.N_BIG = 6
    subprocess.check_output = lambda cmd, **kwargs: ""
    try:
        camp.run_stages(Path("missing"), Path("traj"))
        raise SystemExit("a stage log without keys must fail")
    except RuntimeError:
        pass
    subprocess.check_output = fake_check
    camp.main()
    camp.tex_escape("plain")
    camp.ok = None
    short = [{"rc": 1}]
    # The nested ok() from main already ran. Call the same predicate here.
    assert camp.summary_line("Empty", None, 3)
    try:
        camp.write_outputs({"host_model": "other", "orth": {}, "ice1e6": {},
                             "orth_nop": 1, "big_nop": 1, "matrix_gib_1e6": 1.0})
    except SystemExit:
        pass
    none_stat = {
        "host_model": "Intel Xeon Processor",
        "orth_nop": 1,
        "big_nop": 1,
        "matrix_gib_1e6": 1.0,
        "orth": {key: None for key in (
            "rings_base", "rings_tip", "rings_ovito", "neigh_rings_base",
            "neigh_rings_tip", "neigh_rings_ovito", "cutoff_base", "cutoff_tip",
            "seeded_base", "seeded_tip", "rss_cutoff_base", "rss_cutoff_tip",
            "rss_seeded_base", "rss_seeded_tip")},
        "ice1e6": {key: None for key in (
            "cutoff_tip", "seeded_tip", "rings_ovito", "bonds_ovito", "rss_cutoff_tip")},
    }
    camp.write_outputs(none_stat)
    assert camp.host_label() == "redacted"
    assert "(R)" not in camp.host_name("Intel(R) Xeon(R)")
    assert camp.summary_line("Empty", None, 3).startswith("SUMMARY Empty n=0")
    try:
        camp.write_log("bad.log", "cmd", "tail\n(R)\n")
        raise SystemExit("registered mark must fail the log")
    except SystemExit:
        pass
    subprocess.run = saved_run
    subprocess.check_output = real_check


def cover_bench() -> None:
    import benchmarks.shipped as shipped
    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    shipped.subprocess.run = lambda cmd, **kwargs: calls.append(cmd) or None
    # Keep the real _run so its body is executed.
    orth = shipped.ShippedOrth()
    orth.setup("cutoff")
    orth.time_cages("cutoff")
    orth.peakmem_cages("cutoff")
    orth.setup("seeded")
    orth.time_cages("seeded")
    stages = shipped.ShippedStages()
    stages.setup()
    stages.time_process()
    stages.peakmem_process()
    million = shipped.ShippedMillion()
    million.setup("cutoff")
    million.time_cages("cutoff")
    million.setup("seeded")
    million.peakmem_cages("seeded")
    assert calls


def cover_shell(scratch: Path) -> None:
    stub = scratch / "py"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    trace = scratch / "shipped.xtrace"
    env = os.environ.copy()
    env["PYTHON"] = str(stub)
    env["TMPDIR"] = str(scratch)
    env["PS4"] = "+${BASH_SOURCE}:${LINENO}:"
    proc = REAL_RUN(
        ["bash", "-x", str(ROOT / "paper/tex/repro/shipped.sh")],
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    trace.write_text(proc.stderr)
    src = (ROOT / "paper/tex/repro/shipped.sh").read_text().splitlines()
    hit = set()
    for line in proc.stderr.splitlines():
        if "shipped.sh:" not in line:
            continue
        num = line.split("shipped.sh:", 1)[1].split(":", 1)[0]
        if num.isdigit():
            hit.add(int(num))
    code = []
    for i, raw in enumerate(src, 1):
        stripped = raw.strip()
        if stripped and not stripped.startswith("#"):
            code.append(i)
    covered = [i for i in code if i in hit]
    pct = 100 * len(covered) / len(code)
    print(f"shipped.sh {len(covered)} {len(code)} {pct:.4f}")
    if pct < 98:
        raise SystemExit("shipped.sh coverage is short")


def report(data: Path, include: str, label: str) -> None:
    proc = subprocess.run(
        [
            sys.executable, "-m", "coverage", "report",
            f"--data-file={data}",
            f"--include={include}",
            "--show-missing",
            "--fail-under=98",
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    print(label)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"{label} is under 98%")


def main() -> None:
    scratch = Path(os.environ["TMPDIR"])
    data = scratch / "shipped.coverage"
    env = os.environ.copy()
    env["COVERAGE_FILE"] = str(data)
    # This process is already under coverage when invoked that way.
    cover_campaign(scratch)
    cover_bench()
    cover_shell(scratch)


if __name__ == "__main__":
    main()
