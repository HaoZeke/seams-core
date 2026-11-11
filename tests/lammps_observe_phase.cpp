// Copyright (c) d-SEAMS developers
// A fresh process, so the phase flag is read on its first call.
#include <seams_insitu.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <unistd.h>

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
  const auto note = std::filesystem::temp_directory_path() / "seams-observe-phase.txt";
  FILE *saved = std::fopen(note.string().c_str(), "w");
  if (saved == nullptr) {
    return 2;
  }
  const int fd = fileno(stderr);
  const int copy = dup(fd);
  if (copy < 0) {
    return 2;
  }
  fflush(stderr);
  dup2(fileno(saved), fd);
  int left = 1;
  int right = 1;
  std::thread a([&]() { left = once(); });
  std::thread b([&]() { right = once(); });
  a.join();
  b.join();
  fflush(stderr);
  dup2(copy, fd);
  close(copy);
  std::fclose(saved);
  std::ifstream in(note);
  std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (left != 0 || right != 0) {
    return 1;
  }
  if (text.find("phase filt") == std::string::npos || text.find("phase neigh") == std::string::npos) {
    std::fprintf(stderr, "phase lines missing\n");
    return 3;
  }
  return 0;
}
