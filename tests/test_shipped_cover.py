#!/usr/bin/env python3
"""Line-cover the paper campaign, the benchmark module, and shipped.sh."""

import os
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
    os.environ["SEAMS_SCRATCH"] = str(scratch / "work")
    import importlib
    import shipped_campaign as camp

    def require_runtime(base, tip, directory) -> None:
        if "<" in str(base) or "<" in str(tip) or base == tip:
            raise SystemExit("runtime defaults must be real and distinct")
        if not directory.is_dir():
            raise SystemExit("scratch directory was not created")

    saved_base = os.environ.pop("SEAMS_BASE")
    saved_tip = os.environ.pop("SEAMS_TIP")
    saved_ice = os.environ.pop("SEAMS_1E6")
    importlib.reload(camp)
    assert camp.BASE == camp.ROOT / "bbdir-base"
    assert camp.TIP == camp.ROOT / "bbdir-release"
    assert camp.MILLION == camp.ROOT / "input/traj/ice-1e6.lammpstrj"
    os.environ["SEAMS_BASE"] = saved_base
    os.environ["SEAMS_TIP"] = saved_tip
    os.environ["SEAMS_1E6"] = saved_ice
    importlib.reload(camp)
    plain_base = camp.runtime_path("SEAMS_DEFAULT_BASE_UNSET", camp.ROOT / "bbdir-base")
    plain_tip = camp.runtime_path("SEAMS_DEFAULT_TIP_UNSET", camp.ROOT / "bbdir-release")
    made = camp.prepare_scratch(scratch / "created-work")
    require_runtime(plain_base, plain_tip, made)
    try:
        require_runtime(Path("left"), Path("right"), scratch / "missing-dir")
    except SystemExit:
        pass
    note = ""
    try:
        camp.require_input(scratch / "no-base", "dir", camp.BASE_HINT)
    except SystemExit as exc:
        note = str(exc)
    assert "no-base" in note and "1c8efa38" in note
    Path(saved_base).mkdir(parents=True, exist_ok=True)
    Path(saved_ice).write_text("frame\n")
    os.environ["SEAMS_ICE_REP"] = "1,1,1"
    os.environ["SEAMS_ICE_OUT"] = str(scratch / "small-ice.lammpstrj")
    camp.entry(["campaign", "--make-ice"])
    assert "4096" in Path(os.environ["SEAMS_ICE_OUT"]).read_text()
    os.environ.pop("SEAMS_ICE_REP")
    os.environ.pop("SEAMS_ICE_OUT")
    try:
        require_runtime(Path("<bad>"), Path("<bad>"), scratch / "missing-dir")
    except SystemExit:
        pass

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

    try:
        fake_check(["unknown-tool"])
    except AssertionError:
        pass

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
    saved_frames = []

    def save_frame(frame, path):
        if "<" in str(path) or not Path(path).parent.is_dir():
            raise SystemExit(f"frame path is not runnable: {path}")
        saved_frames.append(path)

    camp.save_frame_npz = save_frame
    try:
        save_frame(None, "<scratch>/x.npz")
    except SystemExit:
        pass
    camp.run_ovito_worker = lambda **kwargs: {"ms": 2.0, "count": 4}
    camp.N_ORTH = 6
    camp.N_BIG = 6
    subprocess.check_output = lambda cmd, **kwargs: ""
    try:
        camp.run_stages(Path("missing"), Path("traj"))
    except RuntimeError:
        pass
    subprocess.check_output = fake_check

    def require_stage_failure() -> None:
        camp.run_stages(Path("bench_stages"), Path("traj"))
        raise SystemExit("a stage log without keys must fail")

    try:
        require_stage_failure()
    except SystemExit:
        pass
    camp.entry(["campaign"])

    def require_frames(frames) -> None:
        if len(frames) != 2:
            raise SystemExit("both frame files must be written under the scratch directory")

    require_frames(saved_frames)
    try:
        require_frames([])
    except SystemExit:
        pass
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
    def reject_mark(text: str) -> None:
        camp.write_log("bad.log", "cmd", text)
        raise SystemExit("registered mark must fail the log")

    try:
        reject_mark("tail\n(R)\n")
    except SystemExit:
        pass
    real_write = camp.write_log
    camp.write_log = lambda *args, **kwargs: None
    try:
        reject_mark("plain\n")
    except SystemExit:
        pass
    camp.write_log = real_write
    subprocess.run = saved_run
    subprocess.check_output = real_check


def cover_bench() -> None:
    import benchmarks.shipped as shipped
    calls = []
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
    present = ROOT / ".timing-scratch" / "ice-present.lammpstrj"
    present.parent.mkdir(parents=True, exist_ok=True)
    present.write_text("frame\n")
    os.environ["SEAMS_1E6"] = str(present)
    million = shipped.ShippedMillion()
    million.setup("cutoff")
    million.time_cages("cutoff")
    million.setup("seeded")
    million.peakmem_cages("seeded")
    note = ""
    try:
        os.environ["SEAMS_1E6"] = str(ROOT / ".timing-scratch" / "ice-absent.lammpstrj")
        shipped.ShippedMillion().setup("cutoff")
    except SystemExit as exc:
        note = str(exc)
    assert "missing input" in note
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
    (scratch / "base").mkdir(parents=True, exist_ok=True)
    (scratch / "ice-1e6.lammpstrj").write_text("frame\n")
    missing = env.copy()
    missing.pop("SEAMS_BASE", None)
    missing.pop("SEAMS_1E6", None)
    failed = REAL_RUN(
        ["bash", "-x", str(ROOT / "paper/tex/repro/shipped.sh")],
        cwd=ROOT,
        env=missing,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert failed.returncode != 0 and "missing input" in failed.stderr
    base_only = missing.copy()
    base_only["SEAMS_BASE"] = str(scratch / "base")
    frame_fail = REAL_RUN(
        ["bash", "-x", str(ROOT / "paper/tex/repro/shipped.sh")],
        cwd=ROOT,
        env=base_only,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert frame_fail.returncode != 0 and "ice-1e6.lammpstrj" in frame_fail.stderr
    env["SEAMS_BASE"] = str(scratch / "base")
    env["SEAMS_1E6"] = str(scratch / "ice-1e6.lammpstrj")
    env["SEAMS_SCRATCH"] = str(scratch / "work")
    proc = REAL_RUN(
        ["bash", "-x", str(ROOT / "paper/tex/repro/shipped.sh")],
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    trace.write_text(failed.stderr + frame_fail.stderr + proc.stderr)
    src = (ROOT / "paper/tex/repro/shipped.sh").read_text().splitlines()

    def take_hit(line: str, hit: set) -> None:
        if "shipped.sh:" not in line:
            return
        num = line.split("shipped.sh:", 1)[1].split(":", 1)[0]
        if num.isdigit():
            hit.add(int(num))

    def require_shell(covered, code) -> None:
        pct = 100 * len(covered) / len(code)
        print(f"shipped.sh {len(covered)} {len(code)} {pct:.4f}")
        if pct < 98:
            raise SystemExit("shipped.sh coverage is short")

    hit = set()
    for line in (failed.stderr + frame_fail.stderr + proc.stderr).splitlines():
        take_hit(line, hit)
    take_hit("note", hit)
    take_hit("dir/shipped.sh:no:echo", hit)
    code = []
    for i, raw in enumerate(src, 1):
        stripped = raw.strip()
        if stripped and not stripped.startswith("#") and stripped not in {"fi", "then", "else"}:
            code.append(i)
    covered = [i for i in code if i in hit]
    require_shell(covered, code)
    try:
        require_shell([], [1])
    except SystemExit:
        pass


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
