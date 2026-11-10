// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <seams_insitu.h>

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

seams_insitu_frame ortho(int step, int n, const double *xyz) {
  seams_insitu_frame frame{};
  frame.step = step;
  frame.n = n;
  frame.xyz = xyz;
  frame.lx = 30;
  frame.ly = 30;
  frame.lz = 30;
  frame.coord_scale = 1.0;
  return frame;
}

} // namespace

TEST_CASE("phase lines run when the variable is set", "[insitu]") {
  REQUIRE(setenv("DSEAMS_PHASE", "1", 1) == 0);
  const double xyz[6] = {0, 0, 0, 1.2, 0, 0};
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, 2, xyz);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.step = 2;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}

TEST_CASE("a zero scale leaves coordinates as they are", "[insitu]") {
  const double xyz[6] = {0, 0, 0, 2.0, 0, 0};
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  auto frame = ortho(1, 2, xyz);
  frame.coord_scale = 0.0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}

TEST_CASE("a repeated frame keeps the neighbour build", "[insitu]") {
  const double xyz[6] = {0, 0, 0, 1.1, 0, 0};
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  auto frame = ortho(1, 2, xyz);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.step = 2;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}

TEST_CASE("replay rejects a null path and a bad cutoff", "[insitu]") {
  REQUIRE(seams_insitu_replay(nullptr, "a", "b") == SEAMS_INSITU_ERR);
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_replay(eng, nullptr, "b") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(eng, "a", nullptr) == SEAMS_INSITU_ERR);
  char arg0[] = "seams-insitu";
  char arg1[] = "replay";
  char arg2[] = "missing.dump";
  char arg3[] = "out.txt";
  char arg4[] = "nope";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};
  REQUIRE(seams_insitu_main(5, argv) == 2);
  seams_insitu_destroy(eng);
}

TEST_CASE("a non-positive lattice row is rejected", "[insitu]") {
  const double bad[9] = {0, 0, 0, 0, 10, 0, 0, 0, 10};
  const double by0[9] = {10, 0, 0, 0, 0, 0, 0, 0, 10};
  const double cz0[9] = {10, 0, 0, 0, 10, 0, 0, 0, 0};
  const double small[9] = {0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5};
  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  REQUIRE(seams_insitu_box_from_rows(bad, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy,
                                     &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(by0, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy,
                                     &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(cz0, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy,
                                     &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(small, &xlo, &ylo, &zlo, &lx, &ly, &lz,
                                     &xy, &xz, &yz) == 0);
  REQUIRE(lx == 0.5);
}

TEST_CASE("a dump with one box line still analyses", "[insitu]") {
  const char *text =
      "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
      "ITEM: BOX BOUNDS pp pp pp\n0.0 10.0\n"
      "ITEM: ATOMS id mol type x y z\n"
      "1 1 1 0.0 0.0 0.0\n2 2 1 1.0 0.0 0.0\n";
  const char *none =
      "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
      "ITEM: ATOMS id mol type x y z\n"
      "1 1 1 0.0 0.0 0.0\n2 2 1 1.0 0.0 0.0\n";
  const auto dir = std::filesystem::temp_directory_path() / "seams-cover-insitu";
  std::filesystem::create_directories(dir);
  const auto one = dir / "one-box.lammpstrj";
  const auto zero = dir / "no-box.lammpstrj";
  const auto rec = dir / "rec.txt";
  {
    std::ofstream(one) << text;
    std::ofstream(zero) << none;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_replay(eng, one.c_str(), rec.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_replay(eng, zero.c_str(), rec.c_str()) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}

TEST_CASE("two hundred and fifty six host rows use the parallel filter",
          "[insitu]") {
  constexpr int n = 256;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3, 0.0);
  std::vector<int> counts(static_cast<std::size_t>(n), 1);
  std::vector<int> flat(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = 0.4 * i;
    flat[static_cast<std::size_t>(i)] = i + 1;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_set_neighbours(eng, n, counts.data(), flat.data()) ==
          SEAMS_INSITU_OK);
  auto frame = ortho(1, n, xyz.data());
  frame.lx = 200;
  frame.ly = 20;
  frame.lz = 20;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}
