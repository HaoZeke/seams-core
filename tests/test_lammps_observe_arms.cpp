// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <seams_insitu.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

seams_insitu_frame ortho(int step, int n, const double *xyz) {
  seams_insitu_frame frame{};
  frame.step = step;
  frame.n = n;
  frame.xyz = xyz;
  frame.lx = 30;
  frame.ly = 30;
  frame.lz = 30;
  return frame;
}

std::string scratch(const std::string &name) {
  return std::string("/tmp/seams-observe-arms-") + name;
}

void require_box_nulls(const double rows[9]) {
  double xlo = 0, ylo = 0, zlo = 0, lx = 0, ly = 0, lz = 0, xy = 0, xz = 0, yz = 0;
  REQUIRE(seams_insitu_box_from_rows(rows, nullptr, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, nullptr, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, nullptr, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, nullptr, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, nullptr, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, nullptr, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, nullptr, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, nullptr, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, nullptr) == 1);
}

} // namespace

TEST_CASE("box rows reject one missing output at a time", "[insitu]") {
  const double rows[9] = {10, 0, 0, 1, 12, 0, 2, 3, 14};
  require_box_nulls(rows);
  double xlo = 0, ylo = 0, zlo = 0, lx = 0, ly = 0, lz = 0, xy = 0, xz = 0, yz = 0;
  const double small[9] = {0.5, 0, 0, 0, 0.4, 0, 0, 0, 0.3};
  REQUIRE(seams_insitu_box_from_rows(small, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 0);
  REQUIRE(lx == 0.5);
  double bad[9] = {0, 0, 0, 0, 10, 0, 0, 0, 10};
  REQUIRE(seams_insitu_box_from_rows(bad, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  bad[0] = 10;
  bad[4] = 0;
  REQUIRE(seams_insitu_box_from_rows(bad, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  bad[4] = 10;
  bad[8] = -1;
  REQUIRE(seams_insitu_box_from_rows(bad, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  const double ay[9] = {10, 1, 0, 0, 10, 0, 0, 0, 10};
  const double az[9] = {10, 0, 1, 0, 10, 0, 0, 0, 10};
  const double bz[9] = {10, 0, 0, 0, 10, 1, 0, 0, 10};
  REQUIRE(seams_insitu_box_from_rows(ay, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(az, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(bz, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
}

TEST_CASE("pack and tally and gather null one argument", "[insitu]") {
  const double xyz[3] = {1, 2, 3};
  double out[3] = {};
  int types_out[1] = {};
  int ids_out[1] = {};
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, nullptr, nullptr, nullptr, out, types_out, ids_out) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, xyz, nullptr, nullptr, nullptr, types_out, ids_out) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, xyz, nullptr, nullptr, out, nullptr, ids_out) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, xyz, nullptr, nullptr, out, types_out, nullptr) == -1);

  int hc = 0, ddc = 0, both = 0, water = 0;
  int types[2] = {1, 2};
  int h[2] = {1, 0};
  int d[2] = {0, 1};
  int codes[2] = {};
  seams_insitu_tally_cages(2, types, 1, h, d, codes, nullptr, &ddc, &both, &water);
  seams_insitu_tally_cages(2, types, 1, h, d, codes, &hc, nullptr, &both, &water);
  seams_insitu_tally_cages(2, types, 1, h, d, codes, &hc, &ddc, nullptr, &water);
  seams_insitu_tally_cages(2, types, 1, h, d, codes, &hc, &ddc, &both, nullptr);
  seams_insitu_tally_cages(2, nullptr, 1, h, d, codes, &hc, &ddc, &both, &water);
  seams_insitu_tally_cages(2, types, 1, nullptr, d, codes, &hc, &ddc, &both, &water);
  seams_insitu_tally_cages(2, types, 1, h, nullptr, codes, &hc, &ddc, &both, &water);
  seams_insitu_tally_cages(2, types, 1, h, d, nullptr, &hc, &ddc, &both, &water);
  seams_insitu_tally_cages(2, types, 1, h, d, codes, &hc, &ddc, &both, &water);
  REQUIRE(hc == 1);
  REQUIRE(ddc == 0);
  REQUIRE(water == 0);

  int counts[1] = {2};
  int displs[1] = {};
  int total = 0;
  REQUIRE(seams_insitu_gather_displs(1, nullptr, displs, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(1, counts, nullptr, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(1, counts, displs, &total) == 0);
  REQUIRE(total == 2);
}

TEST_CASE("a negative id and a zero scale still analyse", "[insitu]") {
  const double xyz[3] = {1, 1, 1};
  int ids[1] = {-5};
  auto *eng = seams_insitu_create(nullptr);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, 1, xyz);
  frame.ids = ids;
  frame.coord_scale = 2;
  const int one = 1;
  const int flat[2] = {-5, 40};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, &one, flat) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const double trio[9] = {0, 0, 0, 8, 0, 0, 0, 8, 0};
  int trio_ids[3] = {1, 2, 3};
  auto *rows = seams_insitu_create(nullptr);
  auto row_frame = ortho(1, 3, trio);
  row_frame.ids = trio_ids;
  const int row_counts[1] = {2};
  const int row_flat[2] = {99, 1};
  REQUIRE(seams_insitu_set_neighbours(rows, 1, row_counts, row_flat) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(rows, &row_frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(rows);
  const int row[1] = {1};
  const int only[1] = {1};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, row, only) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_set_neighbours(eng, 1, nullptr, only) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_set_neighbours(eng, 1, row, nullptr) == SEAMS_INSITU_ERR);
  seams_insitu_destroy(eng);
}

TEST_CASE("debug dump writes and a full device refuses the record", "[insitu]") {
  const double xyz[6] = {0, 0, 0, 1, 0, 0};
  seams_insitu_options opt{};
  seams_insitu_options_default(&opt);
  opt.debug_dump = 1;
  auto *eng = seams_insitu_create(&opt);
  auto frame = ortho(1, 2, xyz);
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_observe(eng, &frame, "") == SEAMS_INSITU_OK);
  const auto path = scratch("debug.lammpstrj");
  std::remove(path.c_str());
  REQUIRE(seams_insitu_observe(eng, &frame, path.c_str()) == SEAMS_INSITU_OK);
  frame.types = nullptr;
  frame.ids = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, path.c_str()) == SEAMS_INSITU_OK);
  std::vector<double> bulk(static_cast<std::size_t>(2500) * 3, 0.0);
  for (int i = 0; i < 2500; ++i) {
    bulk[static_cast<std::size_t>(3 * i)] = static_cast<double>(i) * 10.0;
  }
  auto wide = ortho(2, 2500, bulk.data());
  wide.lx = 30000;
  wide.ly = 30;
  wide.lz = 30;
  REQUIRE(seams_insitu_append_dump(&wide, "/dev/full") == SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_open_record(eng, "/dev/full") == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(eng, &wide) == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(eng);
}

TEST_CASE("replay skips an off-cadence step and rejects a torn timestep", "[insitu]") {
  const double xyz[3] = {4, 5, 6};
  auto frame = ortho(1, 1, xyz);
  const auto dump = scratch("replay.lammpstrj");
  const auto torn = scratch("torn.lammpstrj");
  const auto record = scratch("record.txt");
  std::remove(dump.c_str());
  std::remove(torn.c_str());
  std::remove(record.c_str());
  for (int step = 1; step <= 200; ++step) {
    frame.step = step;
    REQUIRE(seams_insitu_append_dump(&frame, dump.c_str()) == SEAMS_INSITU_OK);
  }
  {
    std::ofstream out(torn);
    out << "ITEM: TIMESTEP\n";
  }
  seams_insitu_options opt{};
  seams_insitu_options_default(&opt);
  opt.cadence.every = 2;
  opt.cadence.first = 0;
  opt.cadence.last = -1;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_replay(nullptr, dump.c_str(), record.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(eng, nullptr, record.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(eng, dump.c_str(), nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(eng, torn.c_str(), record.c_str()) == SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_replay(eng, dump.c_str(), record.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_replay(eng, dump.c_str(), "/dev/full") == SEAMS_INSITU_ERR_IO);
  char arg0[] = "seams-insitu";
  char arg1[] = "replay";
  char arg4[] = "3.5x";
  char *bad[] = {arg0, arg1, const_cast<char *>(dump.c_str()), const_cast<char *>(record.c_str()), arg4};
  REQUIRE(seams_insitu_main(5, bad) == 2);
  char arg5[] = "4.0";
  char *ok[] = {arg0, arg1, const_cast<char *>(dump.c_str()), const_cast<char *>(record.c_str()), arg5};
  REQUIRE(seams_insitu_main(5, ok) == 0);
  seams_insitu_destroy(eng);
}

TEST_CASE("a dilute row of 256 atoms still returns a rings line", "[insitu]") {
  std::vector<double> xyz(static_cast<std::size_t>(256) * 3);
  for (int i = 0; i < 256; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = i * 10.0;
    xyz[static_cast<std::size_t>(3 * i + 1)] = 0;
    xyz[static_cast<std::size_t>(3 * i + 2)] = 0;
  }
  auto *eng = seams_insitu_create(nullptr);
  auto frame = ortho(1, 256, xyz.data());
  frame.lx = 3000;
  frame.ly = 30;
  frame.lz = 30;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(std::strstr(result->text, "rings ") != nullptr);
  frame.lx = 2500;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.lx = 3000;
  frame.ly = 40;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.ly = 30;
  frame.lz = 40;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.lz = 30;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.lz = 50;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
#ifdef _OPENMP
  int nested = SEAMS_INSITU_ERR;
  #pragma omp parallel num_threads(2)
  {
    if (omp_get_thread_num() == 0) {
      nested = seams_insitu_step(eng, &frame);
    }
  }
  REQUIRE(nested == SEAMS_INSITU_OK);
#endif
  seams_insitu_destroy(eng);
}

TEST_CASE("twelve cage oxygens record a hexagonal count", "[insitu]") {
  const double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999},     {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997},  {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996},     {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  double xyz[36];
  for (int i = 0; i < 12; ++i) {
    xyz[3 * i] = coords[i][0];
    xyz[3 * i + 1] = coords[i][1];
    xyz[3 * i + 2] = coords[i][2];
  }
  auto *eng = seams_insitu_create(nullptr);
  auto frame = ortho(1, 12, xyz);
  frame.lx = 53.9690018;
  frame.ly = 54.5289994;
  frame.lz = 51.257;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(result->n_hc + result->n_ddc + result->n_both >= 0);
  REQUIRE(std::strstr(result->text, "hc:") != nullptr);
  seams_insitu_destroy(eng);
}

TEST_CASE("a step issued from an OpenMP worker still returns", "[insitu]") {
  const double xyz[3] = {1, 1, 1};
  auto *eng = seams_insitu_create(nullptr);
  auto frame = ortho(1, 1, xyz);
  int rc = SEAMS_INSITU_ERR;
#ifdef _OPENMP
  #pragma omp parallel num_threads(2)
  {
    if (omp_get_thread_num() == 0) {
      rc = seams_insitu_step(eng, &frame);
    }
  }
#else
  rc = seams_insitu_step(eng, &frame);
#endif
  REQUIRE(rc == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}
