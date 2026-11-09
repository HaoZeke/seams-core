# Copyright (c) d-SEAMS developers
"""Final counts on exampleTraj through the chemfiles read.

The mol column is kept, so these lines match the native reader.
"""

import subprocess
import sys

seams, traj = sys.argv[1], sys.argv[2]
one = ["--frame", "1", "--last", "1"]


def run(args):
    proc = subprocess.run([seams, *args], capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        return None
    return proc.stdout


checks = [
    (["hbonds", traj, *one, "--type", "2", "--htype", "1"],
     "nop 250 hbonds 460\n"),
    (["f4", traj, *one, "--type", "2", "--htype", "1"],
     "nop 750 f4 0.221753 n 250\n"),
    (["ions", traj, *one, "--type", "2", "--ion-types", "1"],
     "nop 750 ice 0 ions 500 in-ice 0 front 0 liquid 500 shell 5.478 shell-ice 0 clusters 0 ions-in-ice 0\n"),
    (["pairs", traj, *one, "--site", "2=cationHead,1=anion"],
     "contact-pair count 0 nCation 250 nAnion 0\n"),
    (["domains", traj, "--site", "2=cationHead,1=anion", "--cutoff", "6", "--subset", "polar"],
     "subset polar n 750 largest 750 P 1\n"),
    (["cn", traj, "--ions", "--site", "2=cationHead,1=anion", "--cutoff", "6"],
     "cage ionCloud types 1 2 cutoff 6 degree 0 nCation 250 nAnion 0\n"),
]
failed = False
for args, expect in checks:
    got = run(args)
    if got != expect:
        sys.stderr.write("got " + repr(got) + "\nwant " + repr(expect) + "\n")
        failed = True
sys.exit(1 if failed else 0)
