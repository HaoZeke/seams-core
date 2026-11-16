#!/bin/sh
# Copyright (c) d-SEAMS developers
# Time findDDC on the ideal-ice ladder.
# Release build, one thread on core 0, one warm-up, six timed calls.
# scripts/pin_run.py sets OMP_NUM_THREADS and RAYON_NUM_THREADS to 1 and
# the mask to cpu 0; its line with the cpus each run used ends
# OUTDIR/ic-N.err and is appended to OUTDIR/cpus.txt.
# Usage: ddc_ladder.sh BENCH OUTDIR
set -eu
BENCH=$1
OUT=$2
WARMUP=1
REPEATS=6
CELLS="3 4 5 6 7 8"
mkdir -p "$OUT/frames"
python3 paper/tex/figures/ladder_frames.py "$OUT/frames" $CELLS
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
PIN_RUN="$(dirname "$0")/pin_run.py"
rm -f "$OUT/cpus.txt"
for m in $CELLS; do
  n=$((8 * m * m * m))
  python3 "$PIN_RUN" --cpus 0 --width 1 --record "$OUT/cpus.txt" -- \
    "$BENCH" "$OUT/frames/ic-$n.lammpstrj" $WARMUP $REPEATS \
    > "$OUT/ic-$n.log" 2> "$OUT/ic-$n.err"
done
