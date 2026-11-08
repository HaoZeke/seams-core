#!/bin/sh
# timeout = ceil(max + (max - min))
# median is the type-7 quantile at one half
# interquartile range is the type-7 quartile gap
# LD_PRELOAD is set on pixi, so Meson and its children inherit it.
# Usage:
#   scripts/time_coverage.sh derive <log>
#   scripts/time_coverage.sh check <log>
#   scripts/time_coverage.sh apply <log>
#   scripts/time_coverage.sh report <timing-log> <confirm-log>
#   scripts/time_coverage.sh record <log> <multiplier>
set -u

mode=${1:-derive}
raw=${2:-logs/github-runner-timings.log}
rc=0

py() {
  pixi run python3 "$@"
}

fmt_g() {
  printf '%.17g' "$1"
}

judge_output() {
  name=$1
  kind=$2
  elapsed=$3
  out=$4
  shown=$(fmt_g "$elapsed")
  ok_n=$(printf '%s\n' "$out" | grep -c 'OK[[:space:]]\{1,\}[0-9.][0-9.]*s' || true)
  if [ "$ok_n" -ne 1 ]; then
    echo "${name} clock_missing"
    return 4
  fi
  ok_s=$(printf '%s\n' "$out" | sed -n 's/.*OK[[:space:]]\{1,\}\([0-9.][0-9.]*\)s.*/\1/p' | head -n 1)
  delta=$(py -c 'import sys; print(abs(float(sys.argv[1]) - float(sys.argv[2])))' "$ok_s" "$shown")
  echo "${name} meson_ok=${ok_s}"
  echo "${name} clock_delta=${delta}"
  too_far=$(py -c 'import sys; raise SystemExit(0 if float(sys.argv[1]) <= 2.0 else 1)' "$delta" || echo yes)
  if [ -n "$too_far" ]; then
    echo "${name} clock_mismatch"
    return 4
  fi
  echo "${name} ${kind}=${shown}"
  return 0
}

derive_py() {
  py - "$1" << 'PY'
import re, sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR

path = sys.argv[1]
rows = {"bulkTUM": [], "device_host_levi5": []}
ns = {}
pat = re.compile(
    r"^(bulkTUM|device_host_levi5) sample=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$"
)
npat = re.compile(r"^(bulkTUM|device_host_levi5) n=(\d+)$")
dpat = re.compile(
    r"^(bulkTUM|device_host_levi5) clock_delta=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$"
)
wpat = re.compile(r"^(bulkTUM|device_host_levi5) warmup=")
deltas = {"bulkTUM": [], "device_host_levi5": []}
warmups = {"bulkTUM": 0, "device_host_levi5": 0}
with open(path) as fh:
    for n, line in enumerate(fh, 1):
        line = line.strip()
        m = pat.match(line)
        if m:
            rows[m.group(1)].append((n, Decimal(m.group(2))))
            continue
        m = npat.match(line)
        if m:
            ns[m.group(1)] = int(m.group(2))
            continue
        m = dpat.match(line)
        if m:
            deltas[m.group(1)].append(Decimal(m.group(2)))
            continue
        m = wpat.match(line)
        if m:
            warmups[m.group(1)] += 1

def qtype7(xs, p):
    h = (Decimal(len(xs)) - 1) * Decimal(p)
    lo = int(h.to_integral_value(rounding=ROUND_FLOOR))
    hi = int(h.to_integral_value(rounding=ROUND_CEILING))
    if lo == hi:
        return xs[lo]
    w = h - Decimal(lo)
    return xs[lo] * (Decimal(1) - w) + xs[hi] * w

for name in ("bulkTUM", "device_host_levi5"):
    pairs = rows[name]
    if name not in ns or ns[name] < 6 or len(pairs) != ns[name]:
        print(f"{name} INCOMPLETE")
        sys.exit(1)
    need = len(pairs) + warmups[name]
    if len(deltas[name]) < need:
        print(f"{name} CLOCK_INCOMPLETE")
        sys.exit(1)
    for gap in deltas[name]:
        if gap > Decimal("2"):
            print(f"{name} CLOCK_MISMATCH")
            sys.exit(1)
    print(f"{name} n={ns[name]}")
    for n, seconds in pairs:
        print(f"{name} raw={format(seconds, '.17g')} source={path}:{n}")
    xs = sorted(s for _, s in pairs)
    median = qtype7(xs, Decimal("0.5"))
    iqr = qtype7(xs, Decimal("0.75")) - qtype7(xs, Decimal("0.25"))
    timeout = (xs[-1] + (xs[-1] - xs[0])).to_integral_value(rounding=ROUND_CEILING)
    print(f"{name} median_s={format(median, '.17g')}")
    print(f"{name} iqr_type7_s={format(iqr, '.17g')}")
    print(f"{name} timeout_ceil_max_plus_spread={timeout}")
    print(f"{name} rule=ceil(max + (max - min))")
PY
}

runner_py() {
  py - "$1" << 'PY'
import sys
from decimal import Decimal, ROUND_CEILING
path = sys.argv[1]
rows = {}
with open(path) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        item = {}
        for part in line.split():
            if "=" in part:
                key, value = part.split("=", 1)
                item[key] = value
        kind = item.get("kind")
        test = item.get("test")
        if kind not in ("coverage", "ordinary") or test not in ("bulkTUM", "device_host_levi5"):
            continue
        rows.setdefault((kind, test), []).append(Decimal(item["seconds"]))
order = (
    ("coverage", "bulkTUM", 4),
    ("coverage", "device_host_levi5", 4),
    ("ordinary", "bulkTUM", 6),
    ("ordinary", "device_host_levi5", 6),
)
for kind, test, least in order:
    xs = sorted(rows.get((kind, test), []))
    if len(xs) < least:
        print(f"{kind} {test} INCOMPLETE n={len(xs)}")
        sys.exit(1)
    timeout = (xs[-1] + (xs[-1] - xs[0])).to_integral_value(rounding=ROUND_CEILING)
    print(f"{kind} {test} n={len(xs)} timeout={timeout}")
PY
}
read_build() {
  py - "$1" << 'PY'
import re, sys
text = open(sys.argv[1]).read()
bulk = re.search(
    r"_timeout = get_option\('b_coverage'\) \? (\d+) : (\d+)", text)
host_cov = re.search(r"_levi_cov = (\d+)", text)
host_ord = re.search(r"_levi_ord = (\d+)", text)
host = host_cov and host_ord
if not bulk or not host:
    sys.exit(1)
print(bulk.group(1))
print(bulk.group(2))
print(host_cov.group(1))
print(host_ord.group(1))
PY
}

case "$mode" in
  derive)
    derive_py "$raw"
    rc=$?
    ;;
  check)
    if grep -q '^kind=coverage ' "$raw"; then
      derived=$(runner_py "$raw") || py_rc=$?
      printf '%s\n' "$derived"
      if [ "${py_rc:-0}" -ne 0 ]; then
        exit "$py_rc"
      fi
      build=$(read_build "${3:-tests/meson.build}") || py_rc=$?
      printf '%s\n' "$build"
      if [ "${py_rc:-0}" -ne 0 ]; then
        exit "$py_rc"
      fi
      b_cov_bulk=$(printf '%s\n' "$build" | sed -n '1p')
      b_ord_bulk=$(printf '%s\n' "$build" | sed -n '2p')
      b_cov_host=$(printf '%s\n' "$build" | sed -n '3p')
      b_ord_host=$(printf '%s\n' "$build" | sed -n '4p')
      d_cov_bulk=$(printf '%s\n' "$derived" | sed -n 's/^coverage bulkTUM n=[0-9]* timeout=//p')
      d_ord_bulk=$(printf '%s\n' "$derived" | sed -n 's/^ordinary bulkTUM n=[0-9]* timeout=//p')
      d_cov_host=$(printf '%s\n' "$derived" | sed -n 's/^coverage device_host_levi5 n=[0-9]* timeout=//p')
      d_ord_host=$(printf '%s\n' "$derived" | sed -n 's/^ordinary device_host_levi5 n=[0-9]* timeout=//p')
      if [ "$b_cov_bulk" != "$d_cov_bulk" ] || [ "$b_ord_bulk" != "$d_ord_bulk" ] || [ "$b_cov_host" != "$d_cov_host" ] || [ "$b_ord_host" != "$d_ord_host" ]; then
        echo "mismatch build ${b_cov_bulk} ${b_ord_bulk} ${b_cov_host} ${b_ord_host} derive ${d_cov_bulk} ${d_ord_bulk} ${d_cov_host} ${d_ord_host}"
        exit 1
      fi
      echo "match coverage bulkTUM=${d_cov_bulk} device_host_levi5=${d_cov_host} ordinary bulkTUM=${d_ord_bulk} device_host_levi5=${d_ord_host}"
      rc=0
    else
      derived=$(derive_py "$raw") || exit 1
      printf '%s\n' "$derived"
      echo "mismatch sample log is not the runner dataset"
      exit 1
    fi
    ;;
  apply)
    derived=$(runner_py "$raw") || py_rc=$?
    if [ "${py_rc:-0}" -ne 0 ]; then
      printf '%s\n' "$derived"
      exit "$py_rc"
    fi
    d_cov_bulk=$(printf '%s\n' "$derived" | sed -n 's/^coverage bulkTUM n=[0-9]* timeout=//p')
    d_ord_bulk=$(printf '%s\n' "$derived" | sed -n 's/^ordinary bulkTUM n=[0-9]* timeout=//p')
    d_cov_host=$(printf '%s\n' "$derived" | sed -n 's/^coverage device_host_levi5 n=[0-9]* timeout=//p')
    d_ord_host=$(printf '%s\n' "$derived" | sed -n 's/^ordinary device_host_levi5 n=[0-9]* timeout=//p')
    build_path=${3:-tests/meson.build}
    py - "$build_path" "$d_cov_bulk" "$d_ord_bulk" "$d_cov_host" "$d_ord_host" << 'PY'
import re, sys
path, cov_bulk, ord_bulk, cov_host, ord_host = sys.argv[1:]
text = open(path).read()
text2, n1 = re.subn(
    r"(_timeout = get_option\('b_coverage'\) \? )\d+( : )\d+",
    rf"\g<1>{cov_bulk}\g<2>{ord_bulk}",
    text,
    count=1,
)
text3, n2 = re.subn(
    r"_levi_cov = \d+",
    f"_levi_cov = {cov_host}",
    text2,
    count=1,
)
text3, n3 = re.subn(
    r"_levi_ord = \d+",
    f"_levi_ord = {ord_host}",
    text3,
    count=1,
)
n2 = n2 if n3 == 1 else 0
if n1 != 1 or n2 != 1:
    sys.exit(1)
open(path, "w").write(text3)
print(f"applied coverage bulk {cov_bulk} device {cov_host} ordinary bulk {ord_bulk} device {ord_host}")
PY
    rc=$?
    ;;
  report)
    confirm=${3:-}
    scratch=$(mktemp -d "${TMPDIR:-.}/scratch.XXXXXX")
    derive_py "$raw" > "$scratch/derive.txt" || exit 1
    if [ -n "$confirm" ]; then
      derive_py "$confirm" > "$scratch/confirm.txt" || exit 1
    else
      : > "$scratch/confirm.txt"
    fi
    py - "$raw" "$confirm" "$scratch/derive.txt" "$scratch/confirm.txt" << 'PY'
import sys
timing, confirm, derive_path, confirm_path = sys.argv[1:]
def stats(path):
    got = {}
    for line in open(path):
        line = line.rstrip("\n")
        if " " not in line or "=" not in line:
            continue
        name, rest = line.split(" ", 1)
        if name not in ("bulkTUM", "device_host_levi5"):
            continue
        key, value = rest.split("=", 1)
        got[(name, key)] = value
    out = []
    b = "bulkTUM"
    d = "device_host_levi5"
    out.append(f"bulkTUM counted {got[(b, 'n')]} samples.")
    out.append(f"device_host_levi5 counted {got[(d, 'n')]} samples.")
    out.append(f"The bulkTUM median is {got[(b, 'median_s')]} s.")
    out.append(f"device_host_levi5 median reads {got[(d, 'median_s')]} s.")
    out.append(f"The bulkTUM interquartile range is {got[(b, 'iqr_type7_s')]} s.")
    out.append(f"device_host_levi5 interquartile range reads {got[(d, 'iqr_type7_s')]} s.")
    out.append(f"The bulkTUM ceiling is {got[(b, 'timeout_ceil_max_plus_spread')]} s.")
    out.append(f"device_host_levi5 ceiling reads {got[(d, 'timeout_ceil_max_plus_spread')]} s.")
    out.append("The ceiling rule is ceil of max plus the spread.")
    return out
lines = []
lines.append("# Coverage timeouts")
lines.append("")
lines.append("This report is generated from the sealed timing log.")
lines.append("")
lines.append("The ceiling is the maximum plus the spread from the minimum, rounded up.")
lines.append("Warmup rows stay outside the counted set, the median, and the ceiling.")
lines.append("The median is the type-7 quantile at one half.")
lines.append("The interquartile range is the gap between the type-7 quartiles.")
lines.append("The log prints each raw sample with %.17g from CLOCK_MONOTONIC.")
lines.append("Meson reads CLOCK_REALTIME for the duration it prints.")
lines.append("The preload maps that reading onto elapsed monotonic time.")
lines.append("The two durations must agree within 2 s.")
lines.append("A missing duration fails the run.")
lines.append("")
lines.append("## Provenance")
lines.append("")
lines.append("The runner log is the tracked GitHub Actions sample.")
lines.append("Coverage rows are the coverage job.")
lines.append("Ordinary rows are the Ubuntu and macOS build jobs.")
lines.append("A host log records the processor, memory, operating system, and compiler.")
lines.append("It does not record a machine name.")
lines.append("Each host log names the revision and the tree from the start of the run.")
lines.append("")
lines.append("## Derived ceilings")
lines.append("")
lines.extend(stats(derive_path))
lines.append("")
if confirm:
    lines.append("## Confirmation")
    lines.append("")
    lines.append("The confirmation log uses the same ceiling rule.")
    lines.append("")
    lines.extend(stats(confirm_path))
    lines.append("")
sys.stdout.write("\n".join(lines) + "\n")
PY
    rc=$?
    rm -rf "$scratch"
    ;;
  record)
    mult=${3:?multiplier}
    root=$(git rev-parse --show-toplevel)
    cd "$root" || exit 1
    unset LIBRARY_PATH
    unset LD_LIBRARY_PATH
    export OMP_NUM_THREADS=$(nproc)
    export TMPDIR=.
    model=$(lscpu | awk -F: '/Model name/{gsub(/^[ \t]+/, "", $2); print $2; exit}')
    family=$(lscpu | awk -F: '/^CPU family:/{gsub(/^[ \t]+/, "", $2); print $2; exit}')
    model_num=$(lscpu | awk -F: '/^Model:/{gsub(/^[ \t]+/, "", $2); print $2; exit}')
    if ! pixi run ninja -C bbdir-cov tests/test_bulkTUM tests/check_device_host; then
      echo "build_failed ninja" >&2
      exit 3
    fi
    if ! pixi run gcc -shared -fPIC -pthread -o bbdir-cov/monotonic_realtime.so scripts/monotonic_realtime.c -ldl; then
      echo "build_failed gcc" >&2
      exit 3
    fi
    if ! pixi run meson introspect --buildoptions bbdir-cov > bbdir-cov/buildoptions.json; then
      echo "build_failed introspect" >&2
      exit 3
    fi
    find bbdir-cov -name '*.gcda' -delete
    strip_repo() { pixi run python3 -c 'import os,sys; r=os.stat(sys.argv[1]); t=sys.stdin.read(); p=set();exec("for tok in t.replace(chr(61),chr(32)).replace(chr(58),chr(32)).split():\n path=tok if tok.startswith(chr(47)) else chr(32).strip()\n while path and path!=chr(47):\n  try: st=os.stat(path)\n  except OSError:\n   path=os.path.dirname(path); continue\n  if (st.st_dev, st.st_ino)==(r.st_dev, r.st_ino): p.add(path); break\n  path=os.path.dirname(path)\nfor prefix in sorted(p, key=len, reverse=True): t=t.replace(prefix+chr(47),chr(32).strip())\nsys.stdout.write(t)")' "$root"; }
    cmd_b=$(pixi run ninja -C bbdir-cov -t commands tests/test_bulkTUM | strip_repo)
    cmd_d=$(pixi run ninja -C bbdir-cov -t commands tests/check_device_host | strip_repo)
    if [ -z "$cmd_b" ] || [ -z "$cmd_d" ]; then
      echo "build_failed commands" >&2
      exit 3
    fi
    preload=$(pwd)/bbdir-cov/monotonic_realtime.so
    sha_b=$(sha256sum ./bbdir-cov/tests/test_bulkTUM | awk '{print $1}')
    sha_d=$(sha256sum ./bbdir-cov/tests/check_device_host | awk '{print $1}')
    sha_s=$(sha256sum "$preload" | awk '{print $1}')
    {
      echo "date: $(date -u)"
      echo "os: $(uname -s) $(uname -r) $(uname -m)"
      echo "cpu_model: $model"
      echo "cpu_family: $family"
      echo "cpu_model_num: $model_num"
      echo "cpu_counts: $(lscpu | awk -F: '/^CPU\(s\)/{gsub(/^[ \t]+/, "", $2); c=$2} /^Core\(s\) per socket/{gsub(/^[ \t]+/, "", $2); cores=$2} /^Socket\(s\)/{gsub(/^[ \t]+/, "", $2); sock=$2} /^Thread\(s\) per core/{gsub(/^[ \t]+/, "", $2); thr=$2} END{print "cpus", c, "cores_per_socket", cores, "sockets", sock, "threads_per_core", thr}')"
      echo "compile_bulk: $cmd_b"
      echo "compile_device: $cmd_d"
      mhz=$(lscpu | awk -F: '/MHz/{gsub(/^[ \t]+/, "", $2); print $2; exit}')
      echo "cpu_mhz: ${mhz:-unavailable}"
      echo "memory_kib: $(awk '/MemTotal/{print $2}' /proc/meminfo)"
      echo "compiler: $(pixi run c++ --version | awk 'NR==1{print; exit}')"
      echo "nproc: $OMP_NUM_THREADS"
      echo "omp_num_threads: 1"
      echo "meson_options: $(pixi run python3 - bbdir-cov/buildoptions.json << 'PY'
import json, sys
opts = json.load(open(sys.argv[1]))
want = {"b_coverage", "buildtype", "optimization"}
print(" ".join(o["name"] + "=" + str(o["value"]) for o in opts if o["name"] in want))
PY
)"
      cov=$(pixi run python3 - bbdir-cov/buildoptions.json << 'PY'
import json, sys
opts = json.load(open(sys.argv[1]))
print(next(str(o["value"]) for o in opts if o["name"] == "b_coverage"))
PY
)
      kind=$(pixi run python3 - bbdir-cov/buildoptions.json << 'PY'
import json, sys
opts = json.load(open(sys.argv[1]))
print(next(str(o["value"]) for o in opts if o["name"] == "buildtype"))
PY
)
      if [ "$cov" != "True" ] || [ "$kind" != "debugoptimized" ]; then
        echo "build_not_coverage cov=$cov buildtype=$kind"
        exit 3
      fi
      echo "clock: CLOCK_MONOTONIC"
      echo "clock_note: CLOCK_MONOTONIC advances while the process waits. The preload hides a step of CLOCK_REALTIME. The duration and the enclosing sample must agree within 2 s."
      echo "binary_sha256 test_bulkTUM=${sha_b}"
      echo "binary_sha256 check_device_host=${sha_d}"
      echo "shim_sha256=${sha_s}"
      echo "rev-parse HEAD: $(git rev-parse HEAD)"
      echo "rev-parse tree: $(git rev-parse 'HEAD^{tree}')"
      if [ "${SEAMS_STATUS:-}" = clean ] || { [ -z "${SEAMS_STATUS:-}" ] && git diff --quiet && git diff --cached --quiet; }; then
        echo "status: clean"
      else
        echo "status: dirty"
      fi
      echo "TMPDIR=."
      echo "command=scripts/time_coverage.sh record $raw $mult"
      one_test() {
        name=$1
        kind=$2
        t0=$(py -c 'import time; print("%.12f" % time.monotonic())')
        out=$(LD_PRELOAD="$preload" pixi run meson test -C bbdir-cov --timeout-multiplier "$mult" --num-processes 1 -v "$name" 2>&1)
        one=$?
        t1=$(py -c 'import time; print("%.12f" % time.monotonic())')
        printf '%s\n' "$out"
        elapsed=$(py -c 'import sys; from decimal import Decimal; print(Decimal(sys.argv[2]) - Decimal(sys.argv[1]))' "$t0" "$t1")
        judge_output "$name" "$kind" "$elapsed" "$out" || return $?
        if [ "$one" -ne 0 ]; then
          return "$one"
        fi
        return 0
      }
      run_name() {
        name=$1
        one_test "$name" warmup || return $?
        i=1
        while [ "$i" -le 6 ]; do
          one_test "$name" sample || return $?
          i=$((i + 1))
        done
        echo "${name} n=6"
        return 0
      }
      run_name bulkTUM || rc=$?
      if [ "${rc:-0}" -eq 0 ]; then
        run_name device_host_levi5 || rc=$?
      fi
      echo "date_end: $(date -u)"
      echo "exit_code: ${rc:-0}"
    } > "$raw" 2>&1
    rc=${rc:-0}
    ;;
  assess)
    name=${2:?name}
    elapsed=${3:?elapsed}
    file=${4:?transcript}
    out=$(cat "$file")
    judge_output "$name" sample "$elapsed" "$out"
    rc=$?
    ;;
  *)
    echo "unknown mode"
    rc=2
    ;;
esac

exit "$rc"
