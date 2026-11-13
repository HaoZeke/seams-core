"""Tests of the brine replica workflow: run temperature, the T_K table and the Snakemake rules.

    python3 -m pytest paper/data/brine-frames/test_brine.py

The Snakemake tests dry-run a copy of the workflow in a scratch tree and
are skipped when snakemake is not installed. No frame or restart is read
from the dataset; the table tests write small labelled frames.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import regenerate  # noqa: E402
import replicas  # noqa: E402

REVISION = "0123456789abcdef0123456789abcdef01234567"


def test_tag_and_deck_temperature():
    assert regenerate.tag_of(10, 1) == "T230_m10_r1"
    assert regenerate.tag_of(20, 0, 245.0) == "T245_m20_r0"
    flat = regenerate.brine_vars(0, "1:00:00", 240.0)
    assert flat[flat.index("T") + 1] == "240.00"
    assert flat[flat.index("seed") + 1] == "03"
    default = regenerate.brine_vars(1, "1:00:00")
    assert default[default.index("T") + 1] == "230.00"


@pytest.mark.parametrize("extra, expected", [([], 230.0), (["--temperature", "245"], 245.0)])
def test_run_one_cli_passes_temperature(monkeypatch, tmp_path, extra, expected):
    seen = {}

    def fake_run_one(work, pairs, replica, genice, lmp, threads, timeout, restarts, temperature=None):
        seen.update(pairs=pairs, replica=replica, threads=threads, temperature=temperature)

    monkeypatch.setattr(regenerate, "run_one", fake_run_one)
    argv = ["regenerate.py", "run-one", "--work", str(tmp_path), "--pairs", "10", "--replica", "0",
            "--threads", "2", "--restarts", str(tmp_path / "r"), "--lmp", "lmp", "--genice", "genice2", *extra]
    monkeypatch.setattr(sys, "argv", argv)
    regenerate.main()
    assert seen == {"pairs": 10, "replica": 0, "threads": 2, "temperature": expected}


def write_frame(path: Path, step: int, pairs: int) -> None:
    """A 20-water slab, crystal below z = 25 A, liquid above, and the ions in the liquid."""
    atoms = [("O", 1.0, 1.0, 1.0 + 2.5 * i, "hexagonal" if i < 10 else "water") for i in range(20)]
    atoms += [("Na", 2.0, 2.0, 40.0, "ion")] * pairs + [("Cl", 3.0, 3.0, 42.0, "ion")] * pairs
    header = f'Lattice="10 0 0 0 10 0 0 0 50" origin="0 0 0" step={step} Properties=species:S:1:pos:R:3:chill:S:1'
    lines = [str(len(atoms)), header] + [" ".join(str(v) for v in atom) for atom in atoms]
    path.write_text("\n".join(lines) + "\n")


def manifest_for(frames: Path, out: Path, names: list[str]) -> Path:
    body = "".join(
        f"{hashlib.sha256((frames / n).read_bytes()).hexdigest()}  data/frames/{n}\n" for n in names
    )
    out.write_text(f"# d-SEAMS/brine-front-230K at {REVISION}\n" + body)
    return out


def read_table(path: Path) -> list[dict]:
    return replicas.read_rows(path)


def test_table_writes_temperature_and_source(tmp_path):
    frames = tmp_path / "frames"
    frames.mkdir()
    names = []
    for frame, step in replicas.FRAMES:
        name = replicas.frame_name("r0", 0, frame, 240)
        write_frame(frames / name, step, 0)
        names.append(name)
    local = replicas.frame_name("r0", 10, "t0", 245)
    write_frame(frames / local, 0, 10)
    manifest = manifest_for(frames, tmp_path / "manifest.txt", names)
    out = tmp_path / "warm.csv"
    replicas.table(frames, out, None, manifest, (local,), (240, 245), ("r0",), (0, 10))
    rows = read_table(out)
    assert out.read_text().splitlines()[0].startswith(f"# d-SEAMS/brine-front-230K at {REVISION}, data/frames/T{{T_K}}_")
    assert len(rows) == 2 * 2 * len(replicas.FRAMES)
    by = {(r["T_K"], r["pairs"], r["frame"]): r for r in rows}
    ok = by[("240", "0", "tend")]
    assert ok["status"] == "ok" and ok["source"] == f"hf@{REVISION[:8]}"
    assert ok["crystal"] == "10" and ok["liquid"] == "10" and ok["velocity_seed"] == "3"
    assert by[("245", "10", "t0")]["source"] == "local"
    assert by[("245", "10", "t0")]["na"] == "10"
    missing = by[("240", "10", "tend")]
    assert missing["status"] == "not_in_dataset" and missing["T_K"] == "240"
    assert by[("245", "0", "t0")]["status"] == "not_in_dataset"


def test_table_refuses_a_frame_that_does_not_match_the_revision(tmp_path):
    frames = tmp_path / "frames"
    frames.mkdir()
    name = replicas.frame_name("r0", 0, "t0", 240)
    write_frame(frames / name, 0, 0)
    manifest = manifest_for(frames, tmp_path / "manifest.txt", [name])
    write_frame(frames / name, 0, 1)
    with pytest.raises(SystemExit, match="does not match"):
        replicas.table(frames, tmp_path / "t.csv", None, manifest, (), (240,), ("r0",), (0,))


def snakemake_tree(tmp_path: Path) -> Path:
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not installed")
    root = tmp_path / "tree"
    rules = root / "paper" / "data" / "brine-frames"
    rules.mkdir(parents=True)
    for name in ("Snakefile", "config.yaml", "replicas.py", "regenerate.py"):
        shutil.copy2(HERE / name, rules / name)
    shutil.copytree(HERE / "deck", rules / "deck")
    config = (rules / "config.yaml").read_text()
    revision = config.split("hf_revision:")[1].split()[0]
    dataset = root / "d-SEAMS"
    (dataset / "brine-front-230K" / "data").mkdir(parents=True)
    (dataset / f"brine-front-230K.{revision}.txt").write_text(f"# d-SEAMS/brine-front-230K at {revision}\n")
    return root


def dry_run(root: Path, *args: str) -> str:
    proc = subprocess.run(
        ["snakemake", "-s", "paper/data/brine-frames/Snakefile", "-n", "-p", "--cores", "8", *args],
        cwd=root, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout + proc.stderr


def test_brine_warm_rule_reads_every_warm_loading(tmp_path):
    root = snakemake_tree(tmp_path)
    out = dry_run(root, "paper/data/brine-frames/warm.csv", "--config", "with_runs=true")
    assert "--temperatures 240 245 --replicas r0 --pairs 0 10 20 40" in out
    for temp in (240, 245):
        for pairs in (0, 10, 20, 40):
            assert f"T{temp}_m{pairs}_r0_tend.xyz" in out
            assert f"--pairs {pairs} --replica 0 --threads 2 --temperature {temp}" in out
    table = next(line for line in out.splitlines() if "--temperatures 240 245" in line)
    assert "T230_" not in table and "T245_m10_r0_t0.xyz" in table


def test_finished_run_ignores_a_newer_manifest(tmp_path):
    root = snakemake_tree(tmp_path)
    restarts = root / "d-SEAMS" / "brine-front-230K" / "data" / "restarts" / "T240_m0_r0"
    restarts.mkdir(parents=True)
    for step in (0, 830000, 1670000, 2500000):
        (restarts / f"step-{step}.restart").write_text("restart\n")
    (restarts / "end.restart").write_text("restart\n")
    later = time.time() + 60
    manifest = next((root / "d-SEAMS").glob("brine-front-230K.*.txt"))
    os.utime(manifest, (later, later))
    os.utime(root / "paper" / "data" / "brine-frames" / "replicas.py", (later + 60, later + 60))
    out = dry_run(root, str(restarts.relative_to(root) / "end.restart"), "--rerun-triggers", "mtime")
    assert "rule brine_run:" not in out
