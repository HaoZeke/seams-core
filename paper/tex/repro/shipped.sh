#!/bin/sh
# Release comparison on one core. Twelve calls on the 110592-molecule frame
# and six on the 995328-molecule frame. Writes the JSON record and the TeX
# the paper inputs.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OVITO_THREAD_COUNT=1
ROOT=$(CDPATH= cd -- "$(dirname "$0")/../../.." && pwd)
export SEAMS_BASE="${SEAMS_BASE:-$ROOT/bbdir-base}"
export SEAMS_TIP="${SEAMS_TIP:-$ROOT/bbdir-release}"
export SEAMS_1E6="${SEAMS_1E6:-$ROOT/input/traj/ice-1e6.lammpstrj}"
export SEAMS_SCRATCH="${SEAMS_SCRATCH:-$ROOT/paper/tex/repro/work}"
if [ ! -d "$SEAMS_BASE" ]; then
  echo "missing input $SEAMS_BASE" >&2
  echo "SEAMS_BASE is a build of commit 1c8efa382c4ae88ebc68e40a991260c4d57edbd5." >&2
  echo "docs/orgmode/paper/versions.org names that commit as the d-SEAMS 1 baseline." >&2
  exit 1
fi
if [ ! -f "$SEAMS_1E6" ]; then
  echo "missing input $SEAMS_1E6" >&2
  echo "SEAMS_1E6 is the 995328-molecule frame." >&2
  echo "pixi run paper-ice-1e6 writes it from input/traj/mW_cubic.lammpstrj." >&2
  exit 1
fi
mkdir -p "$SEAMS_SCRATCH"
PY="${PYTHON:-python3}"
"$PY" paper/tex/repro/shipped_campaign.py
"$PY" paper/tex/repro/check_shipped_numbers.py
