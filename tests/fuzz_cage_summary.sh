#!/bin/bash
# Copyright (c) d-SEAMS developers
# args: expected-summary binary [binary args...]
# The summary line is the accounting. A degree-74 frame must be
# skipped, not checked and not defined-empty.
set -u
expected=$1
shift
out=$("$@")
rc=$?
printf '%s\n' "$out"
if [[ "$rc" -ne 0 ]]; then
  exit "$rc"
fi
summary=$(printf '%s\n' "$out" | awk '/^frames /{line=$0} END{print line}')
case "$summary" in
  *" $expected "*) ;;
  *)
    printf 'summary mismatch: %s\n' "$summary" >&2
    exit 1
    ;;
esac
if printf '%s\n' "$out" | grep -q 'CHECKED'; then
  if printf '%s\n' "$expected" | grep -q 'checked 0'; then
    printf 'a skipped file printed CHECKED\n' >&2
    exit 1
  fi
fi
exit 0
