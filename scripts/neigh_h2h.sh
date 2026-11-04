#!/bin/bash
# One-thread neighbour head-to-head. The log records the commit and every sample.
set -o pipefail
status=0
{
  echo "git rev-parse HEAD"
  git rev-parse HEAD
  echo "git rev-parse HEAD^{tree}"
  git rev-parse "HEAD^{tree}"
  echo "command: OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OVITO_THREAD_COUNT=1 .pixi/envs/default/bin/python tests/bench_h2h_neighbours.py --bench-bin bbdir/tests/bench_neigh_only --output-json paper/tex/figures/data/neigh-head2head.json --dseams-threads 1 --freud-threads 1 --ovito-threads 1 --scipy-threads 1"
  OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OVITO_THREAD_COUNT=1 \
    taskset -c 0 .pixi/envs/default/bin/python tests/bench_h2h_neighbours.py \
    --bench-bin bbdir/tests/bench_neigh_only \
    --output-json paper/tex/figures/data/neigh-head2head.json \
    --dseams-threads 1 --freud-threads 1 --ovito-threads 1 --scipy-threads 1 \
    || status=$?
  echo "EXIT:${status}"
} > bbdir/neigh-h2h.log 2>&1
exit "${status}"
