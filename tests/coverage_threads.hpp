// Copyright (c) d-SEAMS developers
#pragma once

#include <cstdlib>

// A coverage run sets OMP_NUM_THREADS=1. Raising the team inside a test
// races the arc counters and gcov stores a negative branch count.
inline int seams_cap_threads(int requested) {
  if (const char *env = std::getenv("OMP_NUM_THREADS")) {
    char *end = nullptr;
    const long cap = std::strtol(env, &end, 10);
    if (end != env && cap > 0 && requested > static_cast<int>(cap)) {
      return static_cast<int>(cap);
    }
  }
  return requested;
}
