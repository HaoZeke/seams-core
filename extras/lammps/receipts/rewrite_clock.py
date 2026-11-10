# Copyright (c) d-SEAMS developers
"""Recreate a linear series, reading the clock once per commit.

Each commit date is one line from date(1) in America/Chicago. The same
reading is the author date and the committer date. When that reading
matches the previous commit, the next line is another reading. The
second is not computed from a start time, and the process does not sleep.
"""

import os
import subprocess
import sys


NAME = "Rohit Goswami"
EMAIL = "rgoswami@ieee.org"


def run(repo, args, env=None, text=True):
    proc = subprocess.run(
        ["git", "-C", repo, *args],
        check=True,
        capture_output=True,
        text=text,
        env=env,
    )
    return proc.stdout


def sample():
    env = os.environ.copy()
    env["TZ"] = "America/Chicago"
    raw = subprocess.check_output(["date", "+%Y-%m-%d %H:%M:%S.%N %z"], env=env, text=True)
    left, zone = raw.strip().rsplit(" ", 1)
    stamp, frac = left.split(".", 1)
    return stamp + " " + zone, int(frac) / 1000000000.0


def write_read(log, line):
    log.write("read " + line + "\n")
    log.flush()


def reading(log):
    line, _frac = sample()
    write_read(log, line)
    return line


def clock(log, previous):
    while True:
        line, frac = sample()
        write_read(log, line)
        if line != previous and frac < 0.2:
            log.write("use " + line + "\n")
            log.flush()
            return line


def load_messages(path):
    found = {}
    if not path:
        return found
    text = open(path, encoding="utf-8").read()
    for block in text.split("\n---\n"):
        block = block.strip("\n")
        if not block:
            continue
        old, message = block.split("\n", 1)
        found[old] = message if message.endswith("\n") else message + "\n"
    return found


def commit_tree(repo, tree, parent, message, stamp):
    env = os.environ.copy()
    env["GIT_AUTHOR_NAME"] = NAME
    env["GIT_AUTHOR_EMAIL"] = EMAIL
    env["GIT_COMMITTER_NAME"] = NAME
    env["GIT_COMMITTER_EMAIL"] = EMAIL
    env["GIT_AUTHOR_DATE"] = stamp
    env["GIT_COMMITTER_DATE"] = stamp
    proc = subprocess.run(
        ["git", "-C", repo, "commit-tree", tree, "-p", parent],
        input=message,
        text=True,
        capture_output=True,
        check=True,
        env=env,
    )
    return proc.stdout.strip()


def recreate(repo, base, tip, ref, messages, log, previous):
    commits = run(repo, ["rev-list", "--reverse", "%s..%s" % (base, tip)]).split()
    current = run(repo, ["rev-parse", ref]).strip()
    if current != base:
        run(repo, ["update-ref", "-m", "clock rewrite base", ref, base, current])
    parent = base
    for old in commits:
        tree = run(repo, ["rev-parse", old + "^{tree}"]).strip()
        message = messages.get(old, run(repo, ["log", "-1", "--format=%B", old]))
        if not message.endswith("\n"):
            message = message + "\n"
        stamp = clock(log, previous)
        previous = stamp
        new = commit_tree(repo, tree, parent, message, stamp)
        subject = message.splitlines()[0]
        log.write("old %s\nnew %s\nsubject %s\n" % (old, new, subject))
        log.flush()
        run(repo, ["update-ref", "-m", "clock %s" % stamp, ref, new, parent])
        parent = new
    return parent, previous


def commit_head(repo, ref, message_path, log, previous):
    message = open(message_path, encoding="utf-8").read()
    if not message.endswith("\n"):
        message += "\n"
    tree = run(repo, ["write-tree"]).strip()
    parent = run(repo, ["rev-parse", "HEAD"]).strip()
    stamp = clock(log, previous)
    new = commit_tree(repo, tree, parent, message, stamp)
    subject = message.splitlines()[0]
    log.write("old %s\nnew %s\nsubject %s\n" % (parent, new, subject))
    log.flush()
    run(repo, ["update-ref", "-m", "clock %s" % stamp, "HEAD", new, parent])
    branch = run(repo, ["rev-parse", ref]).strip()
    run(repo, ["update-ref", "-m", "clock %s" % stamp, ref, new, branch])
    return new, stamp


def main():
    if sys.argv[1] == "commit":
        repo, ref, log_path, message_path, previous = sys.argv[2:7]
        with open(log_path, "a", encoding="utf-8") as log:
            log.write("repo %s\n" % repo)
            new, last = commit_head(repo, ref, message_path, log, previous)
            log.write("rewritten %s\nlast %s\n" % (new, last))
        print(new)
        print(last)
        return
    repo, base, tip, ref, log_path = sys.argv[1:6]
    messages = load_messages(sys.argv[6] if len(sys.argv) > 6 else "")
    previous = sys.argv[7] if len(sys.argv) > 7 else ""
    with open(log_path, "a", encoding="utf-8") as log:
        log.write("repo %s\nbase %s\ntip %s\nref %s\n" % (repo, base, tip, ref))
        tip_new, last = recreate(repo, base, tip, ref, messages, log, previous)
        log.write("rewritten %s\nlast %s\n" % (tip_new, last))
    print(tip_new)
    print(last)


if __name__ == "__main__":
    main()
