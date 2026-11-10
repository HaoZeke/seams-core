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


def main():
    snapper = os.environ["SNAPPER"]
    guard = os.environ["PROSEGUARD"]
    paper = os.environ["ROHIT_PAPER"]
    core = os.environ.get("PAPER_CORE", "")
    commit_profile = os.environ["COMMIT_PROFILE"]
    log = streamlog.Log(os.environ["PROSE_LOG"])
    log.line("start: " + streamlog.date_is())
    root = Path(os.environ["SRC_SEAMS"])
    log.line("HEAD_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip())
    log.line("TREE_SEAMS: " + subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True).strip())
    version = subprocess.check_output([snapper, "--version"], text=True).strip()
    log.line("snapper: " + version)
    kit = Path(os.environ["KIT"])
    failed = 0
    for path in sorted(kit.rglob("*.md")):
        proc = subprocess.run([snapper, "-f", "plaintext", "-w", "72", "--check", str(path)], text=True, capture_output=True)
        log.line("snapper %s exit %d" % (path.name, proc.returncode))
        if proc.stdout:
            log.line(proc.stdout.rstrip())
        if proc.returncode != 0:
            failed += 1
        for profile, label in ((paper, "rohit-paper"), (core, "paper-core")):
            if not profile:
                log.line("%s %s skipped" % (label, path.name))
                continue
            proc = subprocess.run(
                [guard, "field", "--no-llm", "--profile", profile, "--", str(path)],
                text=True, capture_output=True,
            )
            log.line("%s %s exit %d" % (label, path.name, proc.returncode))
            if proc.returncode != 0:
                failed += 1
                log.line(proc.stdout[-500:])
    repos = [root, Path(os.environ["SRC_LAMMPS"])]
    for repo in repos:
        os.environ["SERIES_BASE"] = os.environ["SEAMS_BASE"] if repo == root else os.environ["LAMMPS_BASE"]
        for index, body in enumerate(messages(repo), 1):
            target = Path(os.environ["BUILD"]) / "work" / ("msg-%s-%d.txt" % (repo.name, index))
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(body)
            for tool, args, label in (
                (snapper, ["-f", "plaintext", "-w", "72", "--check", str(target)], "snapper"),
                (guard, ["field", "--no-llm", "--profile", commit_profile, "--", str(target)], "commit-msg"),
                (guard, ["field", "--no-llm", "--profile", paper, "--", str(target)], "rohit-paper"),
            ):
                proc = subprocess.run([tool, *args], text=True, capture_output=True)
                log.line("%s %s %d exit %d" % (label, repo.name, index, proc.returncode))
                if proc.returncode != 0:
                    failed += 1
                    log.line(body.splitlines()[0] if body else "")
                    log.line((proc.stdout or proc.stderr)[-400:])
            if core:
                proc = subprocess.run(
                    [guard, "field", "--no-llm", "--profile", core, "--", str(target)],
                    text=True, capture_output=True,
                )
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
