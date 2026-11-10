#!/bin/bash
# One profile directory per process. gcov's own merge of a shared
# .gcda stores a negative branch count on a hot untaken arc.
root=${GCOV_SUITE_ROOT:?}
base=$(basename "$1")
export GCOV_PREFIX="${root}/${base}.$$"
mkdir -p "$GCOV_PREFIX"
export GCOV_PREFIX_STRIP="${GCOV_PREFIX_STRIP:?}"
exec "$@"
