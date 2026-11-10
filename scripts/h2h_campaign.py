"""Head-to-head campaign driver used by paper/tex/repro/campaign.sh.

The stage functions live in scripts/h2h.py. This module is the command
line that script calls, and the names paper/tex/figures/chill_neighbours.py
reads.
"""

import argparse
import sys

import h2h

SOURCE_FRAME = 1
DRIVER_FRAME = 1
CUTOFF_A = h2h.CUT
CHILL_NEIGHBOURS = 4
LABELS = (
    "other",
    "hexagonal",
    "cubic",
    "interfacial",
    "hydrate",
    "interfacial-hydrate",
    "mixed",
)
LABEL_ID = {name: index for index, name in enumerate(LABELS)}
# OVITO ChillPlusModifier structure-type integers.
OVITO_CHILL_LABELS = (
    "other",
    "hexagonal",
    "cubic",
    "interfacial",
    "hydrate",
)


def read_lammps(path, frame=SOURCE_FRAME):
    positions, box, _info = h2h.read_lammps_frame(path, frame=frame, otype=1)
    return positions, box


def write_lammps(path, positions, box):
    h2h.write_lammpstrj(path, positions, box)


def ovito_data(positions, box):
    return h2h.ovito_data(positions, box)


def prepare(args):
    h2h.stage_prep()
    print("prepare", args.data)


def case(args):
    print("case", args.system, args.scale, "repeats", h2h.REPS)


def rings(args):
    print("rings", args.system, args.scale)


def summarize(args):
    print("summarize", args.output)


def environment(args):
    h2h.stage_env()
    print("environment", args.output)


def main(argv):
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    prep = sub.add_parser("prepare")
    prep.add_argument("--data", required=True)
    prep.add_argument("--liquid", required=True)
    prep.add_argument("--confined", required=True)
    prep.set_defaults(func=prepare)

    one = sub.add_parser("case")
    one.add_argument("system")
    one.add_argument("scale")
    one.add_argument("--data", required=True)
    one.add_argument("--output", required=True)
    one.add_argument("--bench", required=True)
    one.add_argument("--palmer", required=True)
    one.add_argument("--timeout", required=True)
    one.add_argument("--ring-timeout", required=True)
    one.add_argument("--reference-timeout", required=True)
    one.set_defaults(func=case)

    ring = sub.add_parser("rings")
    ring.add_argument("system")
    ring.add_argument("scale")
    ring.add_argument("--data", required=True)
    ring.add_argument("--output", required=True)
    ring.add_argument("--bench", required=True)
    ring.add_argument("--rings", required=True)
    ring.add_argument("--timeout", required=True)
    ring.set_defaults(func=rings)

    summary = sub.add_parser("summarize")
    summary.add_argument("--output", required=True)
    summary.set_defaults(func=summarize)

    env = sub.add_parser("environment")
    env.add_argument("--output", required=True)
    env.add_argument("--bench", required=True)
    env.add_argument("--library", required=True)
    env.add_argument("--source", required=True)
    env.add_argument("--palmer-source", required=True)
    env.add_argument("--palmer-library", required=True)
    env.add_argument("--rings-archive", required=True)
    env.add_argument("--rings", required=True)
    env.set_defaults(func=environment)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
