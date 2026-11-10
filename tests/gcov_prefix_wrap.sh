#!/bin/bash
# Copyright (c) d-SEAMS developers
# One profile directory per process. gcov's own merge of a shared
# .gcda stores a negative branch count on a hot untaken arc.
# Fuzz binaries write outside the suite root. Their counters are not
# merged into the coverage report.
root=${GCOV_SUITE_ROOT:?}
base=$(basename "$1")
case "$base" in
  fuzz_*) export GCOV_PREFIX="${root}/../gcov-fuzz/${base}.$$" ;;
  *) export GCOV_PREFIX="${root}/${base}.$$" ;;
esac
mkdir -p "$GCOV_PREFIX"
export GCOV_PREFIX_STRIP="${GCOV_PREFIX_STRIP:?}"
exec "$@"
