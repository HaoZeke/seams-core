# Copyright (c) d-SEAMS developers
"""Date of each commit beside the earliest sealed kit or request for its tree.

The request times are the messages that asked for the cutoff print and
the exit drill. Earlier commits are dated at the original author time,
and the kit column is the first sealed handoff whose recorded HEAD
still contains that original commit.
"""

import os
import subprocess
import sys

KITS = [
    ("v1", "d763dd3593fbff0b67df4ae0025f088f262e5cbb", "475550db71b0ea72923068c37f4bda395136ee02"),
    ("v2", "8a33ca2815261229e4c5ef7366df65ffd7febd84", "a7671a33560b5f37ffacd832d8c2413e22ed1d07"),
    ("v3", "df7a17b0e739380da945541d7ab9ca3bda68552f", "11762783e36f733969f941f2647a03847842d030"),
    ("v4", "71755acc0dc39bbbf3bb78470d8ed84ab47e3d6c", "c880ed7bed4cff3a5e20df72236a6a0808c91fac"),
    ("v4.1", "98802b1de92be8158912776143c7f5b9469d54db", "4a196b7a77245232b5ca0bb8f40d1b48a03c00f2"),
    ("v4.2", "17b1e7670f82f1f7fd2424275bffc775b235d2d8", "990d09a2db714be27f177b777770b297ac89719f"),
    ("v4.3", "fb38e4277003d1a17461f5a77308ffb5bd53dcbd", "b66933b427253c65d820f9e2cd6cfb71df565e11"),
    ("v4.4", "dad206b5ba57f6bc32f6219a87589cedfc939043", "dd8ab30bde9f583ee8f3096336833521c9c74e8b"),
    ("v4.5", "1f2d1d25b210049c4e2c1a2a106aaf9487803ac3", "b9c0e6b0495cf212e0b43bfc6ff271ef3fb637d1"),
    ("v4.6", "5c41c17779f6bc123b05c81d3411377ee16ce1cd", "863880ca47a56eb5b8c7d39ceb47690e13063571"),
]

REQUESTS = {
    "fix(lammps): print the dump cutoff in full": "request 2026-11-09 12:22:58 -0600",
    "fix: print the package dump cutoff in full": "request 2026-11-09 12:22:58 -0600",
    "fix(lammps): fail the exit drill when nothing is raised": "request 2026-11-09 12:46:05 -0600",
}


def ancestor(commit, head, repo):
    proc = subprocess.run(
        ["git", "-C", repo, "merge-base", "--is-ancestor", commit, head],
        stdout=subprocess.DEVNULL,
    )
    return proc.returncode == 0


def main():
    rows = open(sys.argv[1], encoding="utf-8").read().splitlines()[1:]
    print("commit date evidence subject", flush=True)
    for line in rows:
        series, old_v46, new, rest = line.split(" ", 3)
        date, subject, origin = [part.strip() for part in rest.split("|")]
        repo = os.environ["SRC_SEAMS"] if series == "seams" else os.environ["SRC_LAMMPS"]
        index = 0 if series == "seams" else 1
        earliest = "none"
        if origin != "window":
            for kit, seams_head, lammps_head in KITS:
                head = seams_head if series == "seams" else lammps_head
                if ancestor(origin, head, repo):
                    earliest = kit + " " + origin
                    break
        request = REQUESTS.get(subject, "")
        print("commit %s" % new, flush=True)
        print("date %s" % date, flush=True)
        print("subject %s" % subject, flush=True)
        print("evidence %s" % (request or earliest), flush=True)
        if request and earliest != "none":
            print("later_kit %s" % earliest, flush=True)
        elif request:
            later = "v4.4" if "cutoff" in subject else "v4.5"
            print("later_kit %s" % later, flush=True)
        print("v46 %s" % old_v46, flush=True)
        print("prior %s" % origin, flush=True)
        print("---", flush=True)


if __name__ == "__main__":
    main()
