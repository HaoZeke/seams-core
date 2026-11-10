#!/bin/bash
# Build the fail-before parent and the fix, and keep both Catch logs.
# Usage: scripts/prism_start_logs.sh BUILD_DIR FAIL_REV FIX_REV OUT_DIR
# The logs name the test and the exit status. Source paths are rewritten
# to paths inside the repository.
set -eu
if [ -z "${PRISM_LOGS_COPIED:-}" ]; then
  copy=$(mktemp)
  cp "$0" "$copy"
  chmod +x "$copy"
  PRISM_LOGS_COPIED=1 exec "$copy" "$@"
fi
BUILD=$1
FAIL_REV=$2
FIX_REV=$3
OUT=$4
ROOT=$(git rev-parse --show-toplevel)
start=$(git rev-parse --abbrev-ref HEAD)
if [ "$start" = HEAD ]; then
  start=$(git rev-parse HEAD)
fi
mkdir -p "$OUT"
run_one() {
  local rev=$1
  local name=$2
  git checkout --detach "$rev"
  ninja -C "$BUILD" -j1 tests/test_cover_orient src/libyodaLib.so
  set +e
  "$BUILD/tests/test_cover_orient" "a two-step shift*" --reporter compact \
    > "$OUT/$name.log" 2>&1
  local status=$?
  set -e
  printf 'exit %s\n' "$status" >> "$OUT/$name.log"
  # Keep the Catch lines. Drop profiler noise and machine prefixes.
  python3 - "$OUT/$name.log" "$ROOT" << 'PY'
import pathlib, sys
path, root = sys.argv[1], sys.argv[2]
lines = []
for line in pathlib.Path(path).read_text().splitlines():
    if "libgcov profiling error" in line:
        continue
    line = line.replace(root + "/", "")
    token = "tests/cover_orient.cpp"
    at = line.find(token)
    if at > 0:
        line = line[at:]
    lines.append(line)
pathlib.Path(path).write_text("\n".join(lines) + "\n")
PY
}
run_one "$FAIL_REV" fail-before
run_one "$FIX_REV" fixed
git checkout "$start"
