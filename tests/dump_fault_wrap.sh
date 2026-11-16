#!/bin/sh
# Copyright (c) d-SEAMS developers
# The fault library has to be preloaded, and libasan has to stay ahead of
# it. Meson's test environment would replace LD_PRELOAD with only the
# fault library, and the driver then aborts before main.
# dyld ignores LD_PRELOAD; the library carries __interpose entries there.
fault=$1
shift
if [ "$(uname -s)" = Darwin ]; then
  export DYLD_INSERT_LIBRARIES="${DYLD_INSERT_LIBRARIES:+$DYLD_INSERT_LIBRARIES:}$fault"
elif [ -n "${ASAN_TEST_PRELOAD:-}" ]; then
  export LD_PRELOAD="${ASAN_TEST_PRELOAD}:$fault"
else
  export LD_PRELOAD="${LD_PRELOAD:+$LD_PRELOAD:}$fault"
fi
exec "$@"
