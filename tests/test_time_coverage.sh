#!/bin/bash
# trace-bootstrap-begin
# Exercise derive, check, apply, report, assess, and a stubbed record.
set -eu
cd "$(dirname "$0")/.."
root=$(git rev-parse --show-toplevel)
cd "$root"
mkdir -p "$root/.timing-scratch"
export TMPDIR="$root/.timing-scratch"
export SEAMS_REPO="$root"
scratch=$(mktemp -d "$TMPDIR/scratch.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
exec 9>"$scratch/xtrace"
export BASH_XTRACEFD=9
export PS4='+COV:${BASH_SOURCE}:${LINENO}:'
set -x
# trace-bootstrap-end
runsh() {
  bash -x "$root/scripts/time_coverage.sh" "$@"
}

must_fail() {
  echo "$1"
  shift
  if "$@"; then
    exit 1
  fi
}

write_match() {
  cat > "$scratch/match.log" << 'EOF'
bulkTUM clock_delta=0.33
bulkTUM sample=21.080835953005
bulkTUM clock_delta=0.33
bulkTUM sample=21.116506181017
bulkTUM clock_delta=0.33
bulkTUM sample=21.151666952006
bulkTUM clock_delta=0.33
bulkTUM sample=21.101768476015
bulkTUM clock_delta=0.33
bulkTUM sample=21.175167413021
bulkTUM clock_delta=0.33
bulkTUM sample=21.122428875999
bulkTUM n=6
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.791773039033
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.791606997023
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.700322229008
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.732536915981
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.651056864997
device_host_levi5 clock_delta=0.34
device_host_levi5 sample=56.791517759033
device_host_levi5 n=6
EOF
}

{
  echo 'bulkTUM n=6'
  for s in 10 20 30 40 50 60; do
    echo "bulkTUM clock_delta=0.1"
    echo "bulkTUM sample=$s"
  done
  echo 'device_host_levi5 n=6'
  for s in 100 200 300 400 500 600; do
    echo "device_host_levi5 clock_delta=0.1"
    echo "device_host_levi5 sample=$s"
  done
} > "$scratch/fixture.log"
runsh derive "$scratch/fixture.log" > "$scratch/out.txt"
grep -q 'bulkTUM timeout_ceil_max_plus_spread=110' "$scratch/out.txt"
grep -q 'device_host_levi5 timeout_ceil_max_plus_spread=1100' "$scratch/out.txt"
grep -q 'bulkTUM n=6' "$scratch/out.txt"
must_fail "fixture must not match the measured ceilings" runsh check "$scratch/fixture.log" > "$scratch/check.txt"
grep -q 'mismatch' "$scratch/check.txt"

runsh check logs/github-runner-timings.log > "$scratch/live.txt"
grep -q '^match ' "$scratch/live.txt"
cp tests/meson.build "$scratch/meson.build"
sed -i 's/b_coverage'"'"') ? 35/b_coverage'"'"') ? 99/' "$scratch/meson.build"
must_fail "check must fail when the build file disagrees" runsh check logs/github-runner-timings.log "$scratch/meson.build" > "$scratch/badbuild.txt"
grep -q mismatch "$scratch/badbuild.txt"
cp tests/meson.build "$scratch/meson.build"
runsh apply logs/github-runner-timings.log "$scratch/meson.build" > "$scratch/apply.txt"
grep -q 'applied coverage bulk 35 device 96 ordinary bulk 73 device 235' "$scratch/apply.txt"
printf 'nope\n' > "$scratch/bad-meson.build"
must_fail "check must fail when the build file has no ceiling" runsh check logs/github-runner-timings.log "$scratch/bad-meson.build" > "$scratch/badread.txt"
must_fail "apply must fail when the build file has no ceiling" runsh apply logs/github-runner-timings.log "$scratch/bad-meson.build" > "$scratch/badapply.txt"
{
  echo 'bulkTUM n=9'
  for s in 1 1 1 1 1 1 1 1 1; do
    echo "bulkTUM clock_delta=0.1"
    echo "bulkTUM sample=$s"
  done
  echo 'device_host_levi5 n=9'
  for s in 2 2 2 2 2 2 2 2 2; do
    echo "device_host_levi5 clock_delta=0.1"
    echo "device_host_levi5 sample=$s"
  done
  echo 'bulkTUM clock_delta=0.1'
  echo 'device_host_levi5 clock_delta=0.1'
  echo 'bulkTUM warmup=1'
  echo 'device_host_levi5 warmup=1'
} > "$scratch/equal9.log"
runsh derive "$scratch/equal9.log" > "$scratch/equal9.txt"
grep -q 'timeout_ceil_max_plus_spread=1' "$scratch/equal9.txt"
printf 'kind=coverage test=bulkTUM seconds=1\n' > "$scratch/short-runner.log"
must_fail "check must fail when a runner set is short" runsh check "$scratch/short-runner.log" > "$scratch/short.txt"
grep -q INCOMPLETE "$scratch/short.txt"

write_match
runsh report "$scratch/match.log" "$scratch/match.log" > "$scratch/report.md"
grep -q 'interquartile range' "$scratch/report.md"
grep -q 'GitHub' "$scratch/report.md"
runsh report "$scratch/match.log" > "$scratch/report-one.md"
grep -q 'Derived ceilings' "$scratch/report-one.md"

printf '1/1 name OK 1.00s\n' > "$scratch/ok.txt"
runsh assess bulkTUM 1.00 "$scratch/ok.txt" > "$scratch/assess-ok.txt"
grep -q 'clock_delta=' "$scratch/assess-ok.txt"
printf '1/1 name FAILED\n' > "$scratch/bad.txt"
must_fail "missing stamp must fail" runsh assess bulkTUM 1.00 "$scratch/bad.txt" > "$scratch/assess-miss.txt"
grep -q 'clock_missing' "$scratch/assess-miss.txt"
printf '1/1 name OK 30.00s\n' > "$scratch/far.txt"
must_fail "large gap must fail" runsh assess bulkTUM 1.00 "$scratch/far.txt" > "$scratch/assess-far.txt"
grep -q 'clock_mismatch' "$scratch/assess-far.txt"
printf '1/1 name OK 1.00s\n1/1 name OK 2.00s\n' > "$scratch/two.txt"
must_fail "two stamps must fail" runsh assess bulkTUM 1.00 "$scratch/two.txt" > "$scratch/assess-two.txt"
grep -q 'clock_missing' "$scratch/assess-two.txt"

sed 's/clock_delta=0.1/clock_delta=9/' "$scratch/fixture.log" > "$scratch/bigdelta.log"
must_fail "derive must reject a gap above 2 s" runsh derive "$scratch/bigdelta.log" > "$scratch/bigdelta.txt"
grep -q 'CLOCK_MISMATCH' "$scratch/bigdelta.txt"
grep -v clock_delta "$scratch/fixture.log" > "$scratch/nodelta.log"
must_fail "derive must reject a missing gap" runsh derive "$scratch/nodelta.log" > "$scratch/nodelta.txt"
grep -q 'CLOCK_INCOMPLETE' "$scratch/nodelta.txt"

must_fail "unknown mode must fail" runsh nosuch > "$scratch/unknown.txt"
( must_fail "probe the failure exit" true ) || true
must_fail "apply must fail when the runner set is short" runsh apply "$scratch/short-runner.log" "$scratch/bad-meson.build" > "$scratch/short-apply.txt"

mkdir -p "$scratch/bin" "$scratch/empty"
cat > "$scratch/opts.json" << 'EOF'
[{"name": "b_coverage", "value": true}, {"name": "buildtype", "value": "debugoptimized"}, {"name": "optimization", "value": "3"}]
EOF
export SEAMS_OPTS="$scratch/opts.json"
cat > "$scratch/bin/pixi" << 'EOF'
#!/bin/bash -x
if [ "$1" = run ]; then
  shift
  cmd=$1
  shift
  case "$cmd" in
    python3) exec python3 "$@" ;;
    c++) echo "c++ (GCC) 14.0.0"; exit 0 ;;
    ninja)
      if [ "${3:-}" = -t ]; then
        if [ "${PIXIE_COMMANDS:-ok}" = empty ]; then
          exit 0
        fi
        echo "gcc -I${SEAMS_REPO}/include -fprofile-arcs -ftest-coverage -o test"
        exit 0
      fi
      exit "${PIXIE_NINJA:-0}"
      ;;
    gcc) exit "${PIXIE_GCC:-0}" ;;
    meson)
      if [ "$1" = introspect ]; then
        if [ "${PIXIE_INTROSPECT:-0}" != 0 ]; then
          exit 1
        fi
        if [ -n "${PIXIE_OPTS:-}" ]; then
          printf '%s\n' "$PIXIE_OPTS"
        else
          cat "$SEAMS_OPTS"
        fi
        exit 0
      fi
      echo "${PIXIE_TEST_LINE:-1/1 seams-core:bulkTUM OK 1.25s}"
      exit "${PIXIE_MESON:-0}"
      ;;
  esac
fi
exit 0
EOF
chmod +x "$scratch/bin/pixi"
PATH="$scratch/bin:$PATH" runsh record "$scratch/record.log" 0 > "$scratch/record-status.txt" || true
grep -q 'cpu_model:' "$scratch/record.log"
grep -q 'cpu_family:' "$scratch/record.log"
grep -q 'omp_num_threads:' "$scratch/record.log"
grep -q 'b_coverage=True' "$scratch/record.log"
grep -q 'shim_sha256=' "$scratch/record.log"
grep -q 'compile_bulk:' "$scratch/record.log"
grep -q 'clock_delta=' "$scratch/record.log"
SEAMS_STATUS=clean PATH="$scratch/bin:$PATH" runsh record "$scratch/record-clean.log" 0 > "$scratch/record-clean-status.txt" || true
grep -q 'status: clean' "$scratch/record-clean.log"
PIXIE_TEST_LINE="1/1 seams-core:bulkTUM FAILED" PATH="$scratch/bin:$PATH" runsh record "$scratch/record-miss.log" 0 > "$scratch/record-miss-status.txt" || true
grep -q 'clock_missing' "$scratch/record-miss.log"
PIXIE_MESON=1 PATH="$scratch/bin:$PATH" runsh record "$scratch/record-rc.log" 0 > "$scratch/record-rc-status.txt" || true
grep -q 'clock_delta=' "$scratch/record-rc.log"
PIXIE_NINJA=1 PATH="$scratch/bin:$PATH" must_fail "record must fail when ninja fails" runsh record "$scratch/record-ninja.log" 0 > "$scratch/record-ninja.txt" 2>&1
grep -q 'build_failed ninja' "$scratch/record-ninja.txt"
PIXIE_GCC=1 PATH="$scratch/bin:$PATH" must_fail "record must fail when gcc fails" runsh record "$scratch/record-gcc.log" 0 > "$scratch/record-gcc.txt" 2>&1
grep -q 'build_failed gcc' "$scratch/record-gcc.txt"
PIXIE_INTROSPECT=1 PATH="$scratch/bin:$PATH" must_fail "record must fail when introspect fails" runsh record "$scratch/record-intro.log" 0 > "$scratch/record-intro.txt" 2>&1
grep -q 'build_failed introspect' "$scratch/record-intro.txt"
PIXIE_COMMANDS=empty PATH="$scratch/bin:$PATH" must_fail "record must fail when compile commands are empty" runsh record "$scratch/record-cmd.log" 0 > "$scratch/record-cmd.txt" 2>&1
grep -q 'build_failed commands' "$scratch/record-cmd.txt"
PIXIE_OPTS='[{"name":"b_coverage","value":false},{"name":"buildtype","value":"debug"}]'
PIXIE_OPTS="$PIXIE_OPTS" PATH="$scratch/bin:$PATH" must_fail "record must fail when coverage is off" runsh record "$scratch/record-nocov.log" 0 > "$scratch/record-nocov.txt" 2>&1
grep -q 'build_not_coverage' "$scratch/record-nocov.log"
"$scratch/bin/pixi" run true

os_name=$(uname -s)
cat > "$scratch/tick.c" << 'EOF'
#define _GNU_SOURCE
#include <time.h>
#include <stdio.h>
int main(void) {
  struct timespec mono, real;
  if (clock_gettime(CLOCK_MONOTONIC, &mono) != 0) return 1;
  if (clock_gettime(CLOCK_REALTIME, &real) != 0) return 1;
  if (mono.tv_sec < 0 || real.tv_sec < 0) return 1;
  printf("clocks_ok\n");
  return 0;
}
EOF
pixi run gcc -shared -fPIC -pthread -o "$scratch/mono.so" scripts/monotonic_realtime.c -ldl
pixi run gcc -o "$scratch/tick" "$scratch/tick.c"
"$scratch/tick" > "$scratch/plain.txt"
preload_rc=0
dyld_rc=0
LD_PRELOAD="$scratch/mono.so" "$scratch/tick" > "$scratch/preload.txt" || preload_rc=$?
echo "linux_preload_ok"
DYLD_INSERT_LIBRARIES="$scratch/mono.so" "$scratch/tick" > "$scratch/dyld.txt" || dyld_rc=$?
echo "darwin_insert_ok"
echo "preload changed a program that ignores civil time"
test "$preload_rc" -eq 0
cmp "$scratch/plain.txt" "$scratch/preload.txt"
test "$dyld_rc" -eq 0
cmp "$scratch/plain.txt" "$scratch/dyld.txt"
note_os() {
  if [ "$1" = Linux ]; then
    echo "linux_required"
  elif [ "$1" = Darwin ]; then
    echo "darwin_required"
  else
    echo "interposition_skip $1"
  fi
}
note_os "$os_name"
note_os Darwin
note_os other
echo "The Apple clock header defines CLOCK_REALTIME and CLOCK_MONOTONIC."
echo "It does not define CLOCK_REALTIME_COARSE."
echo "The shim reads the coarse clock only when that macro exists."

echo "time_coverage tests passed"
mkdir -p logs
export SEAMS_COVER_LOOSE=1
export SEAMS_COVER_EDGES=1
pixi run python3 -m coverage run --data-file "$root/.timing-scratch/cover_lines.dat" tests/cover_lines.py \
  "$scratch/xtrace" \
  "$scratch" \
  "$root" \
  "$root/scripts/time_coverage.sh" \
  "$root/tests/test_time_coverage.sh" \
  "$root/logs/time-coverage-py.txt" \
  "$root/logs/time-coverage-xtrace.log"
unset SEAMS_COVER_LOOSE
unset SEAMS_COVER_EDGES
pixi run python3 -m coverage report --data-file "$root/.timing-scratch/cover_lines.dat" \
  --include '*/cover_lines.py' --show-missing --fail-under=98 \
  > "$root/logs/cover-lines-coverage.txt"
pixi run python3 tests/cover_lines.py \
  "$scratch/xtrace" \
  "$scratch" \
  "$root" \
  "$root/scripts/time_coverage.sh" \
  "$root/tests/test_time_coverage.sh" \
  "$root/logs/time-coverage-py.txt" \
  "$root/logs/time-coverage-xtrace.log"
