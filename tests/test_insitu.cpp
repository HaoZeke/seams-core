#include <catch2/catch_test_macros.hpp>
#include <seams_c_api.h>
#include <seams_insitu.h>

#include <arpa/inet.h>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <netinet/in.h>
#include <sstream>
#include <string>
#include <sys/socket.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace {

std::string readFile(const std::string &path) {
  std::ifstream in(path);
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}

std::string scratch(const std::string &name) {
  return std::string("var/run/seams-insitu-") + std::to_string(static_cast<long long>(::getpid())) +
         "-" + name;
}

int run(const std::string &cmd) {
  const int rc = std::system(cmd.c_str());
  REQUIRE(WIFEXITED(rc));
  return WEXITSTATUS(rc);
}

seams_insitu_frame ortho(int step, int n, const double *xyz) {
  seams_insitu_frame frame{};
  frame.step = step;
  frame.n = n;
  frame.xyz = xyz;
  frame.lx = 20;
  frame.ly = 20;
  frame.lz = 20;
  return frame;
}

const seams_insitu_result *mustLast(seams_insitu_engine *eng) {
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(result->text != nullptr);
  return result;
}

} // namespace

TEST_CASE("defaults and rejected options", "[insitu]") {
  seams_insitu_options opt{};
  seams_insitu_options_default(nullptr);
  seams_insitu_options_default(&opt);
  REQUIRE(opt.cutoff == 3.5);
  REQUIRE(opt.type == 1);
  REQUIRE(opt.ring_depth == 8);
  REQUIRE(opt.write_labels == 1);
  REQUIRE(opt.cadence.every == 1);
  REQUIRE(opt.cadence.first == 0);
  REQUIRE(opt.cadence.last == -1);
  auto *defaults = seams_insitu_create(nullptr);
  REQUIRE(defaults != nullptr);
  seams_insitu_destroy(defaults);
  seams_insitu_destroy(nullptr);
  opt.cutoff = 0;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.cutoff = std::nan("");
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  seams_insitu_options_default(&opt);
  opt.ring_depth = 2;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.ring_depth = 13;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.ring_depth = 3;
  auto *shallow = seams_insitu_create(&opt);
  REQUIRE(shallow != nullptr);
  seams_insitu_destroy(shallow);
  opt.ring_depth = 12;
  auto *deep = seams_insitu_create(&opt);
  REQUIRE(deep != nullptr);
  seams_insitu_destroy(deep);
  seams_insitu_options_default(&opt);
  opt.cadence.every = 0;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.cadence.every = 1;
  opt.cadence.first = -1;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.cadence.first = 0;
  opt.cadence.last = -2;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.cadence.last = -1;
  opt.write_labels = 2;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.write_labels = 0;
  auto *quiet = seams_insitu_create(&opt);
  REQUIRE(quiet != nullptr);
  seams_insitu_destroy(quiet);
  REQUIRE(seams_insitu_frame_bytes() > 0);
  REQUIRE(seams_insitu_options_bytes() > 0);
  REQUIRE(seams_insitu_result_bytes() > 0);
}

TEST_CASE("cadence boundaries", "[insitu]") {
  REQUIRE(seams_insitu_due(nullptr, 0) == -1);
  seams_insitu_cadence bad{0, 0, -1};
  REQUIRE(seams_insitu_due(&bad, 0) == -1);
  seams_insitu_cadence cadence{10, 5, 20};
  REQUIRE(seams_insitu_due(&cadence, 4) == 0);
  REQUIRE(seams_insitu_due(&cadence, 5) == 1);
  REQUIRE(seams_insitu_due(&cadence, 10) == 0);
  REQUIRE(seams_insitu_due(&cadence, 15) == 1);
  REQUIRE(seams_insitu_due(&cadence, 20) == 0);
  REQUIRE(seams_insitu_due(&cadence, 25) == 0);
  seams_insitu_cadence open{1, 0, -1};
  REQUIRE(seams_insitu_due(&open, 0) == 1);
}

TEST_CASE("step gate skips, defers, and rejects", "[insitu]") {
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cadence.every = 2;
  opt.cadence.first = 2;
  opt.cadence.last = 4;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_status(nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_last(nullptr) == nullptr);
  REQUIRE(seams_insitu_last(eng) == nullptr);
  REQUIRE(seams_insitu_set_sink(nullptr, nullptr, nullptr) == SEAMS_INSITU_ERR);
  const double xyz[3] = {0, 0, 0};
  auto frame = ortho(2, 1, xyz);
  REQUIRE(seams_insitu_step(nullptr, &frame) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_step(eng, nullptr) == SEAMS_INSITU_ERR);
  frame.step = 0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_SKIP);
  REQUIRE(seams_insitu_last(eng) == nullptr);
  frame.step = 2;
  frame.gpu_resident = 1;
  frame.gpu_device = 3;
  frame.xyz = nullptr;
  frame.n = 8;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_GPU_DEFERRED);
  const auto *deferred = mustLast(eng);
  REQUIRE(deferred->status == SEAMS_INSITU_GPU_DEFERRED);
  REQUIRE(std::string(deferred->text) == "step 2\ngpu deferred\ndevice 3\nn 8\n");
  frame.step = 1;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_SKIP);
  REQUIRE(seams_insitu_status(eng) == SEAMS_INSITU_SKIP);
  REQUIRE(mustLast(eng)->status == SEAMS_INSITU_GPU_DEFERRED);
  frame.step = 2;
  frame.gpu_resident = 0;
  frame.xyz = xyz;
  frame.n = 1;
  frame.local_only = 1;
  frame.gathered = 0;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_NEED_GATHER);
  REQUIRE(mustLast(eng)->text[0] == 's');
  frame.xyz = xyz;
  frame.n = -1;
  frame.local_only = 0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR);
  frame.n = 1;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR);
  frame.xyz = xyz;
  frame.lx = 0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR_BOX);
  frame.lx = 20;
  frame.local_only = 1;
  frame.gathered = 1;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const auto kept = std::string(mustLast(eng)->text);
  frame.rank = 4;
  frame.nranks = 9;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(std::string(mustLast(eng)->text) == kept);
  frame.step = 3;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_SKIP);
  REQUIRE(seams_insitu_status(eng) == SEAMS_INSITU_SKIP);
  REQUIRE(mustLast(eng)->step == 2);
  REQUIRE(mustLast(eng)->status == SEAMS_INSITU_OK);
  frame.step = 6;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_SKIP);
  seams_insitu_destroy(eng);
}

TEST_CASE("direct analysis matches the published CHILL+ call", "[insitu]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int api[4] = {-1, -1, -1, -1};
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "cutoff", api) == 0);
  auto *eng = seams_insitu_create(nullptr);
  auto frame = ortho(1, 4, xyz);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const auto *result = mustLast(eng);
  REQUIRE(result->n == 4);
  REQUIRE(result->n_rings >= 0);
  int chillSum = 0;
  for (int i = 0; i < SEAMS_INSITU_CHILL_N; ++i) {
    chillSum += result->chill_count[i];
  }
  REQUIRE(chillSum == 4);
  REQUIRE(result->n_hc + result->n_ddc + result->n_both + result->n_water ==
          result->n_water_atoms);
  REQUIRE(result->n_water_atoms == 4);
  for (int i = 0; i < 4; ++i) {
    REQUIRE(result->chill_labels[i] == api[i]);
    REQUIRE(result->cage_codes[i] >= 0);
    REQUIRE(result->cage_codes[i] <= 3);
  }
  const auto again = std::string(result->text);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(std::string(mustLast(eng)->text) == again);
  seams_insitu_options quiet;
  seams_insitu_options_default(&quiet);
  quiet.write_labels = 0;
  auto *bare = seams_insitu_create(&quiet);
  REQUIRE(seams_insitu_step(bare, &frame) == SEAMS_INSITU_OK);
  const auto bareText = std::string(mustLast(bare)->text);
  REQUIRE(bareText.find("labels") == std::string::npos);
  REQUIRE(bareText.find("chill ") != std::string::npos);
  frame.n = 0;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(std::string(mustLast(eng)->text).find("n 0\n") != std::string::npos);
  const int ids[2] = {7, 7};
  const double pair[6] = {0, 0, 0, 1, 0, 0};
  frame = ortho(1, 2, pair);
  frame.ids = ids;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR);
  REQUIRE(mustLast(eng)->n == 0);
  seams_insitu_destroy(bare);
  seams_insitu_destroy(eng);
}

TEST_CASE("type selection and a failing sink", "[insitu]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const int types[4] = {2, 2, 2, 2};
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.type = 0;
  auto *eng = seams_insitu_create(&opt);
  auto frame = ortho(0, 4, xyz);
  frame.types = types;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(mustLast(eng)->n_water_atoms == 4);
  const int mixed[4] = {1, 1, 1, 8};
  opt.type = 1;
  auto *typed = seams_insitu_create(&opt);
  frame.types = mixed;
  REQUIRE(seams_insitu_step(typed, &frame) == SEAMS_INSITU_OK);
  REQUIRE(mustLast(typed)->n_water_atoms == 3);
  REQUIRE(mustLast(typed)->chill_labels[3] == 6);
  REQUIRE(mustLast(typed)->cage_codes[3] == 0);
  int calls = 0;
  auto boom = [](const seams_insitu_result *, void *user) -> int {
    *static_cast<int *>(user) += 1;
    return 4;
  };
  REQUIRE(seams_insitu_set_sink(typed, boom, &calls) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(typed, &frame) == SEAMS_INSITU_ERR_IO);
  REQUIRE(calls == 1);
  REQUIRE(seams_insitu_status(typed) == SEAMS_INSITU_ERR_IO);
  REQUIRE(mustLast(typed)->status == SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_set_sink(typed, nullptr, nullptr) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(typed, &frame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(typed);
  seams_insitu_destroy(eng);
}

TEST_CASE("dump round trip matches replay", "[insitu]") {
  const double xyz[12] = {1.0 / 3.0, -2.5, 1e-12, 2.8, 0, 0, 1.4, 2.42, 0.1, 1.1, 0.7, 2.3};
  const int types[4] = {1, 1, 1, 1};
  const int ids[4] = {4, 1, 9, 2};
  auto *eng = seams_insitu_create(nullptr);
  const auto dump = scratch("round.lammpstrj");
  const auto live = scratch("live.txt");
  const auto post = scratch("post.txt");
  std::remove(dump.c_str());
  REQUIRE(seams_insitu_open_record(nullptr, live.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, "") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, "/no/such/dseams-dir/live.txt") == SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_close_record(nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_close_record(eng) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_open_record(eng, live.c_str()) == SEAMS_INSITU_OK);
  auto frame = ortho(0, 4, xyz);
  frame.types = types;
  frame.ids = ids;
  frame.xy = 1.25;
  frame.lx = 12;
  frame.ly = 10;
  frame.lz = 11;
  frame.xlo = 0.5;
  frame.ylo = -0.25;
  frame.zlo = 0.125;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const double second[12] = {0.2, 0.2, 0.2, 3.1, 0.4, 0.1, 1.2, 2.7, 0.3, 0.4, 1.1, 2.6};
  frame.step = 2;
  frame.xyz = second;
  REQUIRE(seams_insitu_observe(eng, &frame, dump.c_str()) == SEAMS_INSITU_OK);
  frame.step = 4;
  frame.xyz = xyz;
  REQUIRE(seams_insitu_observe(eng, &frame, dump.c_str()) == SEAMS_INSITU_OK);
  frame.step = 5;
  frame.lx = 0;
  REQUIRE(seams_insitu_observe(eng, &frame, dump.c_str()) == SEAMS_INSITU_ERR_BOX);
  REQUIRE(seams_insitu_close_record(eng) == SEAMS_INSITU_OK);
  auto *replay = seams_insitu_create(nullptr);
  REQUIRE(seams_insitu_replay(nullptr, dump.c_str(), post.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(replay, nullptr, post.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(replay, dump.c_str(), "/no/such/dseams-dir/post.txt") ==
          SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_replay(replay, dump.c_str(), post.c_str()) == SEAMS_INSITU_OK);
  const auto liveText = readFile(live);
  const auto postText = readFile(post);
  const auto marker = liveText.find("step 2\n");
  REQUIRE(marker != std::string::npos);
  REQUIRE(liveText.substr(marker) == postText);
  const auto bad = scratch("bad.txt");
  {
    std::ofstream out(bad);
    out << "hello\n";
  }
  REQUIRE(seams_insitu_replay(replay, bad.c_str(), post.c_str()) == SEAMS_INSITU_ERR);
  {
    std::ofstream out(bad);
    out << "ITEM: TIMESTEP\nnope\n";
  }
  REQUIRE(seams_insitu_replay(replay, bad.c_str(), post.c_str()) == SEAMS_INSITU_ERR_IO);
  const auto missing = scratch("missing.lammpstrj");
  std::remove(missing.c_str());
  REQUIRE(seams_insitu_replay(replay, missing.c_str(), post.c_str()) == SEAMS_INSITU_ERR_IO);
  const char *exe = std::getenv("SEAMS_INSITU_EXE");
  REQUIRE(exe != nullptr);
  REQUIRE(run(std::string(exe)) == 2);
  REQUIRE(run(std::string(exe) + " no") == 2);
  REQUIRE(run(std::string(exe) + " replay " + missing + " " + post) == 1);
  REQUIRE(run(std::string(exe) + " replay " + dump + " " + post + " abc") == 2);
  REQUIRE(run(std::string(exe) + " replay " + dump + " " + post + " 0") == 2);
  const auto cli = scratch("cli.txt");
  REQUIRE(run(std::string(exe) + " replay " + dump + " " + cli + " 3.5") == 0);
  REQUIRE(readFile(cli) == postText);
  seams_insitu_destroy(replay);
  seams_insitu_destroy(eng);
  std::remove(dump.c_str());
  std::remove(live.c_str());
  std::remove(post.c_str());
  std::remove(bad.c_str());
  std::remove(cli.c_str());
}

TEST_CASE("append rejects a bad frame and accepts an empty one", "[insitu]") {
  const auto path = scratch("append.lammpstrj");
  std::remove(path.c_str());
  const double xyz[3] = {1, 2, 3};
  auto frame = ortho(8, 1, xyz);
  REQUIRE(seams_insitu_append_dump(nullptr, path.c_str()) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_append_dump(&frame, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_append_dump(&frame, "") == SEAMS_INSITU_ERR);
  frame.n = -3;
  REQUIRE(seams_insitu_append_dump(&frame, path.c_str()) == SEAMS_INSITU_ERR);
  frame.n = 1;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, path.c_str()) == SEAMS_INSITU_ERR);
  frame.xyz = xyz;
  frame.lx = -1;
  REQUIRE(seams_insitu_append_dump(&frame, path.c_str()) == SEAMS_INSITU_ERR_BOX);
  frame.lx = 20;
  REQUIRE(seams_insitu_append_dump(&frame, "/no/such/dseams-dir/a.lammpstrj") ==
          SEAMS_INSITU_ERR_IO);
  frame.n = 0;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, path.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(readFile(path).find("NUMBER OF ATOMS\n0\n") != std::string::npos);
  auto *eng = seams_insitu_create(nullptr);
  REQUIRE(seams_insitu_observe(eng, &frame, "/no/such/dseams-dir/b.lammpstrj") ==
          SEAMS_INSITU_ERR_IO);
  const auto second = scratch("other.txt");
  REQUIRE(seams_insitu_open_record(eng, path.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_open_record(eng, second.c_str()) == SEAMS_INSITU_OK);
  frame.step = 1;
  REQUIRE(seams_insitu_observe(eng, &frame, "") == SEAMS_INSITU_OK);
  auto *full = seams_insitu_create(nullptr);
  REQUIRE(seams_insitu_open_record(full, "/dev/full") == SEAMS_INSITU_OK);
  frame.step = 2;
  REQUIRE(seams_insitu_observe(full, &frame, "") == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(full);
  REQUIRE(seams_insitu_close_record(eng) == SEAMS_INSITU_OK);
  REQUIRE(readFile(second).find("step 1\n") != std::string::npos);
  seams_insitu_destroy(eng);
  std::remove(path.c_str());
  std::remove(second.c_str());
}

TEST_CASE("pack, gather, scale, box rows, and cage tally", "[insitu]") {
  REQUIRE(seams_insitu_pack_owned(-1, nullptr, 1, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr) == -1);
  REQUIRE(seams_insitu_pack_owned(0, nullptr, 1, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr) == 0);
  const double xyz[9] = {0, 0, 0, 1, 0, 0, 2, 0, 0};
  double out[9];
  int typesOut[3];
  int idsOut[3];
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, nullptr, nullptr, nullptr, out, typesOut,
                                  idsOut) == -1);
  const int mask[3] = {1, 2, 1};
  const int types[3] = {3, 4, 5};
  const int ids[3] = {10, 11, 12};
  REQUIRE(seams_insitu_pack_owned(3, mask, 1, xyz, types, ids, out, typesOut, idsOut) == 2);
  REQUIRE(idsOut[0] == 10);
  REQUIRE(idsOut[1] == 12);
  REQUIRE(typesOut[1] == 5);
  REQUIRE(out[3] == 2);
  REQUIRE(seams_insitu_pack_owned(3, mask, 0, xyz, nullptr, nullptr, out, typesOut, idsOut) ==
          0);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 0, xyz, nullptr, nullptr, out, typesOut, idsOut) ==
          2);
  REQUIRE(typesOut[0] == 1);
  REQUIRE(idsOut[1] == 2);
  int total = -1;
  REQUIRE(seams_insitu_gather_displs(-1, nullptr, nullptr, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(0, nullptr, nullptr, nullptr) == 1);
  REQUIRE(seams_insitu_gather_displs(0, nullptr, nullptr, &total) == 0);
  REQUIRE(total == 0);
  const int counts[3] = {2, 0, 3};
  int displs[3] = {-1, -1, -1};
  REQUIRE(seams_insitu_gather_displs(3, nullptr, displs, &total) == 1);
  const int badCounts[1] = {-2};
  REQUIRE(seams_insitu_gather_displs(1, badCounts, displs, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(3, counts, displs, &total) == 0);
  REQUIRE(total == 5);
  REQUIRE(displs[0] == 0);
  REQUIRE(displs[1] == 2);
  REQUIRE(displs[2] == 2);
  seams_insitu_scale_nm(nullptr, 1, out);
  seams_insitu_scale_nm(xyz, 0, out);
  seams_insitu_scale_nm(xyz, 1, nullptr);
  const double nm[3] = {0.1, -0.2, 0.3};
  seams_insitu_scale_nm(nm, 1, out);
  REQUIRE(out[0] == 1.0);
  REQUIRE(out[1] == -2.0);
  REQUIRE(out[2] == 3.0);
  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  REQUIRE(seams_insitu_box_from_rows(nullptr, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  const double rows[9] = {10, 0, 0, 1.5, 8, 0, 0.25, -0.5, 9};
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 0);
  REQUIRE(xlo == 0);
  REQUIRE(lx == 10);
  REQUIRE(ly == 8);
  REQUIRE(lz == 9);
  REQUIRE(xy == 1.5);
  REQUIRE(xz == 0.25);
  REQUIRE(yz == -0.5);
  double tilted[9] = {10, 0.2, 0, 0, 8, 0, 0, 0, 9};
  REQUIRE(seams_insitu_box_from_rows(tilted, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  tilted[1] = 0;
  tilted[0] = -10;
  REQUIRE(seams_insitu_box_from_rows(tilted, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  int hcCount = 1, ddcCount = 1, both = 1, water = 1;
  seams_insitu_tally_cages(4, nullptr, 1, nullptr, nullptr, nullptr, nullptr, &ddcCount, &both,
                           &water);
  seams_insitu_tally_cages(0, types, 1, types, types, typesOut, &hcCount, &ddcCount, &both,
                           &water);
  REQUIRE(hcCount == 0);
  REQUIRE(water == 0);
  const int atomTypes[4] = {1, 1, 1, 2};
  const int hc[4] = {1, 0, 1, 1};
  const int ddc[4] = {0, 1, 1, 1};
  int codes[4] = {-1, -1, -1, -1};
  seams_insitu_tally_cages(4, atomTypes, 1, hc, ddc, codes, &hcCount, &ddcCount, &both, &water);
  REQUIRE(codes[0] == 1);
  REQUIRE(codes[1] == 2);
  REQUIRE(codes[2] == 3);
  REQUIRE(codes[3] == 0);
  const double hex[18] = {2.5, 0, 0, 1.25, 2.165, 0, -1.25, 2.165, 0,
                          -2.5, 0, 0, -1.25, -2.165, 0, 1.25, -2.165, 0};
  seams_insitu_options hexOpt;
  seams_insitu_options_default(&hexOpt);
  hexOpt.type = 0;
  auto *hexEng = seams_insitu_create(&hexOpt);
  auto hexFrame = ortho(1, 0, nullptr);
  REQUIRE(seams_insitu_step(hexEng, &hexFrame) == SEAMS_INSITU_OK);
  hexFrame = ortho(2, 6, hex);
  REQUIRE(seams_insitu_step(hexEng, &hexFrame) == SEAMS_INSITU_OK);
  REQUIRE(mustLast(hexEng)->n_rings >= 1);
  seams_insitu_destroy(hexEng);
  REQUIRE(hcCount == 1);
  REQUIRE(ddcCount == 1);
  REQUIRE(both == 1);
  REQUIRE(water == 0);
}

namespace {

int listen_port(int *fd) {
  *fd = ::socket(AF_INET, SOCK_STREAM, 0);
  REQUIRE(*fd >= 0);
  int yes = 1;
  ::setsockopt(*fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));
  sockaddr_in addr{};
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
  addr.sin_port = 0;
  REQUIRE(::bind(*fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) == 0);
  REQUIRE(::listen(*fd, 1) == 0);
  socklen_t len = sizeof(addr);
  REQUIRE(::getsockname(*fd, reinterpret_cast<sockaddr *>(&addr), &len) == 0);
  return ntohs(addr.sin_port);
}

void send_all(int fd, const void *data, int n) {
  const auto *p = static_cast<const char *>(data);
  int sent = 0;
  while (sent < n) {
    const ssize_t w = ::send(fd, p + sent, static_cast<std::size_t>(n - sent), 0);
    REQUIRE(w > 0);
    sent += static_cast<int>(w);
  }
}

} // namespace

TEST_CASE("interactive molecular dynamics client round trip", "[insitu]") {
  const double box[9] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  REQUIRE(seams_insitu_imd_session(nullptr, 1, "a", "b", box, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("abc", 1, "a", "b", box, 1) == SEAMS_INSITU_ERR);
  int closed = -1;
  const int closed_port = listen_port(&closed);
  ::close(closed);
  REQUIRE(seams_insitu_imd_session("127.0.0.1", closed_port, scratch("x").c_str(),
                                   scratch("y").c_str(), box, 1) == SEAMS_INSITU_ERR_IO);

  int badfd = -1;
  const int bad_port = listen_port(&badfd);
  std::thread bad([&] {
    const int client = ::accept(badfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 9;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    ::close(client);
    ::close(badfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", bad_port, scratch("bad.lammpstrj").c_str(),
                                   scratch("bad.txt").c_str(), box, 1) == SEAMS_INSITU_ERR);
  bad.join();

  int earlyfd = -1;
  const int early_port = listen_port(&earlyfd);
  std::thread early([&] {
    const int client = ::accept(earlyfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    REQUIRE(::recv(client, go, 8, MSG_WAITALL) == 8);
    ::close(client);
    ::close(earlyfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", early_port, scratch("early.lammpstrj").c_str(),
                                   scratch("early.txt").c_str(), box, 1) == SEAMS_INSITU_ERR_IO);
  early.join();

  int srv = -1;
  const int port = listen_port(&srv);
  const auto dump = scratch("imd.lammpstrj");
  const auto live = scratch("imd-live.txt");
  const auto post = scratch("imd-post.txt");
  std::remove(dump.c_str());
  std::thread server([&] {
    const int client = ::accept(srv, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    REQUIRE(::recv(client, go, 8, MSG_WAITALL) == 8);
    type = htonl(1);
    std::uint32_t length = htonl(1);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    unsigned char energy[40] = {};
    const std::int32_t step = 11;
    std::memcpy(energy, &step, 4);
    send_all(client, energy, 40);
    type = htonl(2);
    length = htonl(2);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    const float xyz[6] = {0.25f, 0.5f, 0.75f, 2.5f, 0.f, 1.f};
    send_all(client, xyz, static_cast<int>(sizeof(xyz)));
    ::close(client);
    ::close(srv);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_OK);
  server.join();
  auto *eng = seams_insitu_create(nullptr);
  REQUIRE(seams_insitu_replay(eng, dump.c_str(), post.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(readFile(live) == readFile(post));
  REQUIRE(readFile(live).find("step 11\n") == 0);
  seams_insitu_options rare;
  seams_insitu_options_default(&rare);
  rare.cadence.every = 100;
  rare.cadence.first = 100;
  auto *skip = seams_insitu_create(&rare);
  const auto skipped = scratch("imd-skip.txt");
  REQUIRE(seams_insitu_replay(skip, dump.c_str(), skipped.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(readFile(skipped).empty());
  seams_insitu_destroy(skip);
  const auto cut = scratch("cut.lammpstrj");
  {
    std::ofstream out(cut);
    out << "ITEM: TIMESTEP\n";
  }
  REQUIRE(seams_insitu_replay(eng, cut.c_str(), post.c_str()) == SEAMS_INSITU_ERR_IO);
  REQUIRE(seams_insitu_replay(eng, dump.c_str(), ".") == SEAMS_INSITU_ERR_IO);
  int shortfd = -1;
  const int short_port = listen_port(&shortfd);
  std::thread shorter([&] {
    const int client = ::accept(shortfd, nullptr, nullptr);
    ::close(client);
    ::close(shortfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", short_port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_ERR_IO);
  shorter.join();
  int oddfd = -1;
  const int odd_port = listen_port(&oddfd);
  std::thread odd([&] {
    const int client = ::accept(oddfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    REQUIRE(::recv(client, go, 8, MSG_WAITALL) == 8);
    type = htonl(9);
    std::uint32_t length = htonl(0);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    ::close(client);
    ::close(oddfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", odd_port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_ERR);
  odd.join();
  int recfd = -1;
  const int rec_port = listen_port(&recfd);
  std::thread recorder([&] {
    const int client = ::accept(recfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    ::close(client);
    ::close(recfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", rec_port, dump.c_str(),
                                   "/no/such/dseams-dir/imd.txt", box, 1) == SEAMS_INSITU_ERR_IO);
  recorder.join();
  int negfd = -1;
  const int neg_port = listen_port(&negfd);
  std::thread neg([&] {
    const int client = ::accept(negfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    type = htonl(2);
    std::uint32_t length = htonl(0xffffffffu);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    ::close(client);
    ::close(negfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", neg_port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_ERR);
  neg.join();
  int thin_fd = -1;
  const int thin_port = listen_port(&thin_fd);
  std::thread thin([&] {
    const int client = ::accept(thin_fd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    type = htonl(1);
    std::uint32_t length = htonl(1);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    ::close(client);
    ::close(thin_fd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", thin_port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_ERR_IO);
  thin.join();
  int trunc_fd = -1;
  const int trunc_port = listen_port(&trunc_fd);
  std::thread short_coord([&] {
    const int client = ::accept(trunc_fd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    type = htonl(2);
    std::uint32_t length = htonl(2);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    ::close(client);
    ::close(trunc_fd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", trunc_port, dump.c_str(), live.c_str(), box, 1) ==
          SEAMS_INSITU_ERR_IO);
  short_coord.join();
  const double bad_box[9] = {0, 0, 0, 0, 20, 20, 0, 0, 0};
  int boxfd = -1;
  const int box_port = listen_port(&boxfd);
  std::thread boxer([&] {
    const int client = ::accept(boxfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    type = htonl(2);
    std::uint32_t length = htonl(1);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    const float xyz[3] = {1.f, 2.f, 3.f};
    send_all(client, xyz, static_cast<int>(sizeof(xyz)));
    ::close(client);
    ::close(boxfd);
  });
  REQUIRE(seams_insitu_imd_session("127.0.0.1", box_port, scratch("badbox.lammpstrj").c_str(),
                                   scratch("badbox.txt").c_str(), bad_box, 1) != SEAMS_INSITU_OK);
  boxer.join();
  int mainfd = -1;
  const int main_port = listen_port(&mainfd);
  std::thread main_srv([&] {
    const int client = ::accept(mainfd, nullptr, nullptr);
    std::uint32_t type = htonl(4);
    std::uint32_t version = 2;
    send_all(client, &type, 4);
    send_all(client, &version, 4);
    char go[8];
    ::recv(client, go, 8, MSG_WAITALL);
    type = htonl(2);
    std::uint32_t length = htonl(1);
    send_all(client, &type, 4);
    send_all(client, &length, 4);
    const float xyz[3] = {0.f, 0.f, 0.f};
    send_all(client, xyz, static_cast<int>(sizeof(xyz)));
    ::close(client);
    ::close(mainfd);
  });
  char port_text[16];
  std::snprintf(port_text, sizeof(port_text), "%d", main_port);
  char arg_host[] = "127.0.0.1";
  char arg_dump[] = "var/run/seams-imd-main.lammpstrj";
  char arg_rec[] = "var/run/seams-imd-main.txt";
  char a0[] = "0";
  char a20[] = "20";
  char a1[] = "1";
  char *imd_argv[] = {arg_host, arg_host, port_text, arg_dump, arg_rec, a0, a0, a0,
                      a20, a20, a20, a0, a0, a0, a1};
  REQUIRE(seams_insitu_imd_main(15, imd_argv) == 0);
  main_srv.join();
  char bad_num[] = "zz";
  char *bad_box_argv[] = {arg_host, arg_host, port_text, arg_dump, arg_rec, bad_num, a0, a0,
                          a20, a20, a20, a0, a0, a0, a1};
  REQUIRE(seams_insitu_imd_main(15, bad_box_argv) == 2);
  char *bad_max_argv[] = {arg_host, arg_host, port_text, arg_dump, arg_rec, a0, a0, a0,
                          a20, a20, a20, a0, a0, a0, bad_num};
  REQUIRE(seams_insitu_imd_main(15, bad_max_argv) == 2);
  seams_insitu_destroy(eng);
  const char *exe = std::getenv("SEAMS_GRO_IMD_EXE");
  REQUIRE(exe != nullptr);
  REQUIRE(run(std::string(exe)) == 2);
  REQUIRE(run(std::string(exe) + " 127.0.0.1 nope a b 0 0 0 1 1 1 0 0 0 1") == 2);
  std::remove(dump.c_str());
  std::remove(live.c_str());
  std::remove(post.c_str());
}
