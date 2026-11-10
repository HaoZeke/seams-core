// Copyright (c) d-SEAMS developers
// A fresh process, so the phase flag is read on its first call.
#include <seams_insitu.h>

#include <cstdio>
#include <cstdlib>
#include <thread>

namespace {

int once() {
  const double xyz[3] = {1.0, 2.0, 3.0};
  seams_insitu_frame frame{};
  frame.step = 1;
  frame.n = 1;
  frame.xyz = xyz;
  frame.lx = 10;
  frame.ly = 10;
  frame.lz = 10;
  auto *eng = seams_insitu_create(nullptr);
  if (eng == nullptr) {
    return 2;
  }
  const int rc = seams_insitu_step(eng, &frame);
  seams_insitu_destroy(eng);
  return rc == SEAMS_INSITU_OK ? 0 : 1;
}

} // namespace

int main() {
  if (std::getenv("DSEAMS_PHASE") == nullptr) {
    std::fprintf(stderr, "phase flag missing\n");
    return 2;
  }
  int left = 1;
  int right = 1;
  std::thread a([&]() { left = once(); });
  std::thread b([&]() { right = once(); });
  a.join();
  b.join();
  return left == 0 && right == 0 ? 0 : 1;
}
