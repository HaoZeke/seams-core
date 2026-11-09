// Copyright (c) d-SEAMS developers
// The process environment sets DSEAMS_PIN and DSEAMS_PHASE before this
// image is loaded. The worker placement runs from a global constructor,
// and the phase clocks are latched on the first analysis.

#include <seams_insitu.h>

int main() {
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  seams_insitu_engine *eng = seams_insitu_create(&opt);
  if (eng == nullptr) {
    return 1;
  }
  const double xyz[6] = {0.0, 0.0, 0.0, 2.0, 0.0, 0.0};
  const int types[2] = {1, 1};
  const int ids[2] = {1, 2};
  seams_insitu_frame frame{};
  frame.step = 0;
  frame.n = 2;
  frame.xyz = xyz;
  frame.types = types;
  frame.ids = ids;
  frame.lx = 10.0;
  frame.ly = 10.0;
  frame.lz = 10.0;
  const int status = seams_insitu_step(eng, &frame);
  seams_insitu_destroy(eng);
  return status == 0 ? 0 : 1;
}
