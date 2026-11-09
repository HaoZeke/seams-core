# Copyright (c) d-SEAMS developers
"""seams rejects bad arguments with the named error and exit 2.

usage: cli_rejects.py SEAMS TRAJECTORY
"""

import subprocess
import sys
from pathlib import Path

seams = sys.argv[1]
traj = sys.argv[2]


def run(args):
    return subprocess.run([seams, *args], capture_output=True, text=True)


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(name + " " + detail)


def reject(args, message):
    got = run(args)
    first = got.stderr.splitlines()[0] if got.stderr else ""
    check("exit " + message, got.returncode == 2, str(got.returncode))
    check("message " + message, first == message, first)
    check("no frame " + message, "nop " not in got.stdout)


def main():
    reject([], "A command and a trajectory file are required")
    got = run(["no-such-command", traj])
    check("unknown exit", got.returncode == 2, str(got.returncode))
    check("unknown text", got.stdout.strip() == "unknown command: no-such-command",
          got.stdout.strip())
    reject(["--format", "yaml", "read", traj], "bad --format (want text|json)")
    reject(["density-z", traj, "--axis", "q"], "bad --axis (want x|y|z)")
    reject(["rdf", traj, "--bins", "-1"], "bad --bins (want N > 0)")
    reject(["rdf", traj, "--types", "1"], "bad --types (want I,J)")
    reject(["rdf", traj, "--types", "1,2,3"], "bad --types (want I,J)")
    reject(["rdf", traj, "--types", "a,b"], "bad --types (want I,J)")
    reject(["rdf", traj, "--types", "1,"], "bad --types (want I,J)")
    reject(["rdf", traj, "--types", ",1"], "bad --types (want I,J)")
    reject(["ions", traj], "ions needs --ion-types")
    reject(["pairs", traj], "pairs needs --site")
    reject(["cn", traj, "--ions"], "cn --ions needs --site")
    reject(["domains", traj], "domains needs --site")
    reject(["domains", traj, "--site", "1=polar", "--subset", "sideways"],
           "bad --subset (want polar|apolar)")
    reject(["chill", traj, "--family", "ionicLiquid"],
           "CHILL/TUM refused for family ionicLiquid")
    reject(["chill", traj, "--family", "not-a-family"],
           "unknown family 'not-a-family' (want waterIce, ionicLiquid, moltenSalt, des, electrolyte, confinedIL, confinedWater, networkFormer)")
    reject(["pairs", traj, "--site", "1=not-a-kind"],
           "unknown site kind 'not-a-kind'")
    reject(["pairs", traj, "--site", "no-equals"],
           "site spec token 'no-equals' is not key=value")
    reject(["cages", traj, "--cage-edge", "nope"],
           "bad cage edge window (want LO,HI or none)")
    reject(["rings", traj, "--frame", "1-2"],
           "invalid arguments: value \"1-2\" is not a number")
    empty = Path(traj).parent / "cli-empty.xyz"
    empty.write_text("")
    try:
        reject(["read", str(empty)], "XYZ read failed: the file has no frames")
    finally:
        empty.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
