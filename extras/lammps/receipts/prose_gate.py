#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run snapper and both prose profiles on kit markdown and commit messages."""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlog


def messages(repo):
    base = os.environ["SERIES_BASE"]
    text = subprocess.check_output(
        ["git", "-C", str(repo), "log", "--format=%B%x00", "%s..HEAD" % base],
        text=True,
    )
    return [part.strip() + "\n" for part in text.split("\0") if part.strip()]


def guard_cmd(profile, fields, zone=None):
    cmd = [
        os.environ["PROSEGUARD"], "field", "--no-llm",
        "--profile", profile,
        "--styles-path", os.environ["PROSE_STYLES"],
    ]
    if zone:
        cmd.extend(["--zone", zone])
    for spec in fields:
        cmd.extend(["--field", spec])
    return cmd


def other_findings(text, skip_we):
    found = []
    for line in text.splitlines():
        if not line.startswith("field:"):
            continue
        if skip_we and "voice.WeRate" in line:
            continue
        found.append(line)
    return found


def run_guard(cmd):
    return subprocess.run(
        cmd,
        text=True,
        capture_output=True,
        cwd=os.environ.get("PROSE_CWD") or None,
    )


def log_proc(log, label, proc):
    log.line("%s exit %d" % (label, proc.returncode))
    if proc.stdout:
        log.line(proc.stdout.rstrip())
    if proc.stderr:
        log.line(proc.stderr.rstrip())


def note_core(log):
    core = os.environ.get("PAPER_CORE", "")
    if core and Path(core).is_file():
        return core
    if not getattr(note_core, "said", False):
        note_core.said = True
        log.line(
            "paper-core is not on this machine. The profile file is absent, so this VM does not run that check."
        )
    return ""


def main():
    snapper = os.environ["SNAPPER"]
    paper = os.environ["ROHIT_PAPER"]
    commit_profile = os.environ["COMMIT_PROFILE"]
    log = streamlog.Log(os.environ["PROSE_LOG"])
    log.line("start: " + streamlog.date_is())
    root = Path(os.environ["SRC_SEAMS"])
    log.line("HEAD_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True).strip())
    lammps = Path(os.environ["SRC_LAMMPS"])
    log.line("HEAD_LAMMPS: " + subprocess.check_output(["git", "-C", str(lammps), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_LAMMPS: " + subprocess.check_output(["git", "-C", str(lammps), "rev-parse", "HEAD^{tree}"], text=True).strip())
    version = subprocess.check_output([snapper, "--version"], text=True).strip()
    log.line("snapper: " + version)
    kit = Path(os.environ["KIT"])
    failed = 0
    for path in sorted(kit.rglob("*.md")):
        proc = subprocess.run([snapper, "-f", "plaintext", "-w", "72", "--check", str(path)], text=True, capture_output=True)
        log_proc(log, "snapper %s" % path.name, proc)
        if "would reformat" in (proc.stdout or ""):
            failed += 1
        if proc.returncode != 0:
            failed += 1
        for zone in (None, "results"):
            proc = run_guard(guard_cmd(paper, ["body=@%s" % path], zone))
            label = "rohit-paper %s" % path.name if zone is None else "rohit-paper results %s" % path.name
            log_proc(log, label, proc)
            if other_findings((proc.stdout or "") + (proc.stderr or ""), False) or proc.returncode != 0:
                failed += 1
        core = note_core(log)
        if core:
            proc = run_guard(guard_cmd(core, ["body=@%s" % path]))
            log.line("paper-core %s exit %d" % (path.name, proc.returncode))
            if proc.returncode != 0:
                failed += 1
    repos = [root, lammps]
    for repo in repos:
        os.environ["SERIES_BASE"] = os.environ["SEAMS_BASE"] if repo == root else os.environ["LAMMPS_BASE"]
        for index, body in enumerate(messages(repo), 1):
            target = Path(os.environ["BUILD"]) / "work" / ("msg-%s-%d.txt" % (repo.name, index))
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(body)
            lines = body.splitlines()
            subject = lines[0] if lines else ""
            rest = "\n".join(lines[2:]) + "\n"
            body_path = target.with_suffix(".body")
            body_path.write_text(rest)
            proc = subprocess.run([snapper, "-f", "plaintext", "-w", "72", "--check", str(target)], text=True, capture_output=True)
            log_proc(log, "snapper %s %d" % (repo.name, index), proc)
            if proc.returncode != 0 or "would reformat" in (proc.stdout or ""):
                failed += 1
                log.line(subject)
            proc = run_guard(guard_cmd(commit_profile, ["Subject=%s" % subject, "Body=@%s" % body_path]))
            log_proc(log, "commit-msg %s %d" % (repo.name, index), proc)
            if proc.returncode != 0:
                failed += 1
                log.line(subject)
            proc = run_guard(guard_cmd(paper, ["body=@%s" % target]))
            log_proc(log, "rohit-paper %s %d" % (repo.name, index), proc)
            if other_findings((proc.stdout or "") + (proc.stderr or ""), True):
                failed += 1
                log.line(subject)
            core = note_core(log)
            if core:
                proc = run_guard(guard_cmd(core, ["body=@%s" % target]))
                log.line("paper-core %s %d exit %d" % (repo.name, index, proc.returncode))
                if proc.returncode != 0:
                    failed += 1
    log.line("failures: %d" % failed)
    log.line("end: " + streamlog.date_is())
    log.close()
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
