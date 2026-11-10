// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <seams_c_api.h>
#include <seams_insitu.h>

#include <arpa/inet.h>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <new>
#include <netinet/in.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sstream>
#include <string>
#include <thread>

extern "C" void __gcov_dump();
#include <csignal>
#include <sys/resource.h>
#include <sys/socket.h>
#include <sys/wait.h>
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
  const auto dir = std::filesystem::temp_directory_path() / "seams-insitu";
  std::filesystem::create_directories(dir);
  return (dir / (std::to_string(static_cast<long long>(::getpid())) + "-" + name)).string();
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

long virtual_kilobytes() {
  std::ifstream status("/proc/self/status");
  std::string line;
  long value = 0;
  while (std::getline(status, line)) {
    if (line.rfind("VmSize:", 0) == 0) {
      value = std::strtol(line.c_str() + 7, nullptr, 10);
    }
  }
  return value;
}

void write_lattice(const std::string &path, int n) {
  const int side = 14;
  std::ofstream out(path);
  out << "ITEM: TIMESTEP\n0\n";
  out << "ITEM: NUMBER OF ATOMS\n" << n << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n0 40\n0 40\n0 40\n";
  out << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < n; ++i) {
    const int x = i % side;
    const int y = (i / side) % side;
    const int z = i / (side * side);
    out << (i + 1) << " 1 " << (x * 2.2) << " " << (y * 2.2) << " " << (z * 2.2) << "\n";
  }
}

int replay_at_limit(const std::string &dump, const std::string &record) {
  const long used = virtual_kilobytes();
  struct rlimit saved {};
  getrlimit(RLIMIT_AS, &saved);
  struct rlimit tight {};
  tight.rlim_cur = static_cast<rlim_t>(used) * 1024 + (static_cast<rlim_t>(8) << 20);
  tight.rlim_max = saved.rlim_max;
  setrlimit(RLIMIT_AS, &tight);
  auto *eng = seams_insitu_create(nullptr);
  const int rc = seams_insitu_replay(eng, dump.c_str(), record.c_str());
  const int status = seams_insitu_status(eng);
  seams_insitu_destroy(eng);
  setrlimit(RLIMIT_AS, &saved);
  const bool reported = rc == SEAMS_INSITU_ERR && status == SEAMS_INSITU_ERR;
  return reported ? 0 : 1;
}

TEST_CASE("0 allocation failure is reported by replay", "[insitu]") {
  const auto dump = scratch("alloc.lammpstrj");
  const auto record = scratch("alloc.txt");
  write_lattice(dump, 2200);
  int fds[2];
  REQUIRE(::pipe(fds) == 0);
  const pid_t pid = ::fork();
  REQUIRE(pid >= 0);
  if (pid == 0) {
    ::close(fds[0]);
    const char outcome = static_cast<char>(replay_at_limit(dump, record));
    ::write(fds[1], &outcome, 1);
    ::close(fds[1]);
    __gcov_dump();
    _exit(0);
  }
  ::close(fds[1]);
  int status = 0;
  REQUIRE(::waitpid(pid, &status, 0) == pid);
  REQUIRE(WIFEXITED(status));
  REQUIRE(WEXITSTATUS(status) == 0);
  char outcome = 9;
  const ssize_t got = ::read(fds[0], &outcome, 1);
  ::close(fds[0]);
  REQUIRE(got == 1);
  REQUIRE(outcome == 0);
  std::remove(dump.c_str());
  std::remove(record.c_str());
}

void throw_bad_alloc() { throw std::bad_alloc(); }

TEST_CASE("bad_alloc can be caught", "[insitu]") {
  bool saw = false;
  try {
    throw_bad_alloc();
  } catch (const std::bad_alloc &) {
    saw = true;
  }
  REQUIRE(saw);
}

TEST_CASE("a signalled child fails the allocation check", "[negative]") {
  const pid_t pid = ::fork();
  if (pid == 0) {
    std::abort();
    _exit(1);
  }
  int status = 0;
  REQUIRE(::waitpid(pid, &status, 0) == pid);
  REQUIRE(WIFEXITED(status));
  REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("an uncaught bad_alloc fails the allocation check", "[negative]") {
  const pid_t pid = ::fork();
  if (pid == 0) {
    std::thread([] { throw_bad_alloc(); }).join();
    _exit(1);
  }
  int status = 0;
  REQUIRE(::waitpid(pid, &status, 0) == pid);
  REQUIRE(WIFEXITED(status));
  REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("a thread pool panic fails the allocation check", "[negative]") {
  const pid_t pid = ::fork();
  if (pid == 0) {
    struct rlimit lim {};
    getrlimit(RLIMIT_NPROC, &lim);
    lim.rlim_cur = 1;
    setrlimit(RLIMIT_NPROC, &lim);
    const auto dump = scratch("panic.lammpstrj");
    const auto record = scratch("panic.txt");
    write_lattice(dump, 64);
    auto *eng = seams_insitu_create(nullptr);
    seams_insitu_replay(eng, dump.c_str(), record.c_str());
    seams_insitu_destroy(eng);
    _exit(0);
  }
  int status = 0;
  REQUIRE(::waitpid(pid, &status, 0) == pid);
  REQUIRE(WIFEXITED(status));
  REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("defaults and rejected options", "[insitu]") {
  seams_insitu_options opt{};
  opt.cutoff = 9;
  seams_insitu_options_default(nullptr);
  REQUIRE(opt.cutoff == 9);
  seams_insitu_options_default(&opt);
  REQUIRE(opt.cutoff == 3.5);
  REQUIRE(opt.type == 1);
  REQUIRE(opt.ring_depth == 8);
  REQUIRE(opt.write_labels == 1);
  REQUIRE(opt.debug_dump == 0);
  REQUIRE(opt.cadence.every == 1);
  REQUIRE(opt.cadence.first == 0);
  REQUIRE(opt.cadence.last == -1);
  auto *defaults = seams_insitu_create(nullptr);
  REQUIRE(defaults != nullptr);
  seams_insitu_destroy(defaults);
  seams_insitu_destroy(nullptr);
  auto *afterNull = seams_insitu_create(nullptr);
  REQUIRE(afterNull != nullptr);
  seams_insitu_destroy(afterNull);
  opt.cutoff = 0;
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  opt.cutoff = std::nan("");
  REQUIRE(seams_insitu_create(&opt) == nullptr);
  seams_insitu_options_default(&opt);
  opt.debug_dump = 2;
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

TEST_CASE("step above 2^31 is kept", "[insitu]") {
  auto *eng = seams_insitu_create(nullptr);
  const double xyz[3] = {0, 0, 0};
  auto frame = ortho(2, 1, xyz);
  frame.step = (static_cast<int64_t>(1) << 31) + 5;
  frame.local_only = 1;
  frame.gathered = 1;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(mustLast(eng)->step == frame.step);
  REQUIRE(std::string(mustLast(eng)->text).find("2147483653") != std::string::npos);
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
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.debug_dump = 1;
  auto *eng = seams_insitu_create(&opt);
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
          SEAMS_INSITU_OK);
  seams_insitu_options dbg;
  seams_insitu_options_default(&dbg);
  dbg.debug_dump = 1;
  auto *writing = seams_insitu_create(&dbg);
  REQUIRE(seams_insitu_observe(writing, &frame, "/no/such/dseams-dir/b.lammpstrj") ==
          SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(writing);
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

TEST_CASE("a host neighbour list is the cutoff graph for one frame", "[insitu]") {
  REQUIRE(seams_insitu_set_neighbours(nullptr, 0, nullptr, nullptr) == SEAMS_INSITU_ERR);
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.type = 0;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  REQUIRE(seams_insitu_set_neighbours(eng, -1, nullptr, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_set_neighbours(eng, 1, nullptr, nullptr) == SEAMS_INSITU_ERR);
  const int zeroCount[1] = {0};
  const int oneId[1] = {1};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, zeroCount, oneId) == SEAMS_INSITU_ERR);
  const double hex[18] = {2.5, 0, 0, 1.25, 2.165, 0, -1.25, 2.165, 0,
                          -2.5, 0, 0, -1.25, -2.165, 0, 1.25, -2.165, 0};
  const int alone[6] = {1, 1, 1, 1, 1, 1};
  const int ids[6] = {1, 2, 3, 4, 5, 6};
  REQUIRE(seams_insitu_set_neighbours(eng, 6, alone, ids) == SEAMS_INSITU_OK);
  auto frame = ortho(1, 6, hex);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(mustLast(eng)->n_rings == 0);
  REQUIRE(seams_insitu_set_neighbours(eng, 0, nullptr, nullptr) == SEAMS_INSITU_OK);
  frame.step = 2;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const int built = mustLast(eng)->n_rings;
  REQUIRE(built >= 1);
  double moved[18];
  for (int i = 0; i < 18; ++i) {
    moved[i] = hex[i] + ((i % 3) == 0 ? 0.02 : 0.0);
  }
  auto *fresh = seams_insitu_create(&opt);
  auto movedFrame = ortho(1, 6, moved);
  REQUIRE(seams_insitu_step(fresh, &movedFrame) == SEAMS_INSITU_OK);
  const std::string freshText(mustLast(fresh)->text);
  frame.step = 3;
  frame.xyz = moved;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  const std::string kept(mustLast(eng)->text);
  const auto body = [](const std::string &text) {
    const auto nl = text.find('\n');
    return nl == std::string::npos ? text : text.substr(nl);
  };
  REQUIRE(body(kept) == body(freshText));
  seams_insitu_destroy(fresh);
  seams_insitu_destroy(eng);
  char arg0[] = "seams-insitu";
  char arg1[] = "nope";
  char arg2[] = "a";
  char arg3[] = "b";
  char *argv[] = {arg0, arg1, arg2, arg3};
  REQUIRE(seams_insitu_main(4, argv) == 2);
  // Unknown ids, a repeated image and a non-water neighbour are ignored.
  const int counts[6] = {4, 1, 1, 1, 1, 1};
  const int mixed[9] = {1, 99, 2, 2, 2, 3, 4, 5, 6};
  int types[6] = {1, 2, 1, 1, 1, 1};
  auto *mixedEng = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_set_neighbours(mixedEng, 6, counts, mixed) == SEAMS_INSITU_OK);
  auto mixedFrame = ortho(1, 6, hex);
  mixedFrame.types = types;
  REQUIRE(seams_insitu_step(mixedEng, &mixedFrame) == SEAMS_INSITU_OK);
  seams_insitu_destroy(mixedEng);
  int dupIds[2] = {1, 1};
  const double two[6] = {0, 0, 0, 1, 0, 0};
  auto *dup = seams_insitu_create(&opt);
  auto dupFrame = ortho(1, 2, two);
  dupFrame.ids = dupIds;
  REQUIRE(seams_insitu_step(dup, &dupFrame) == SEAMS_INSITU_ERR);
  seams_insitu_destroy(dup);
}

TEST_CASE("neighbour timing entry points", "[insitu]") {
  REQUIRE(seams_bench_open(0, nullptr, nullptr) == nullptr);
  const double xyz[6] = {0, 0, 0, 1, 0, 0};
  const double box[3] = {10, 10, 10};
  auto *bench = seams_bench_open(2, xyz, box);
  REQUIRE(bench != nullptr);
  REQUIRE(seams_bench_cutoff(nullptr, 3.5) == -1);
  REQUIRE(seams_bench_cutoff(bench, 0.0) == -1);
  REQUIRE(seams_bench_cutoff(bench, 3.5) >= 0);
  REQUIRE(seams_bench_knn(nullptr, 4) == -1);
  REQUIRE(seams_bench_knn(bench, 0) == -1);
  REQUIRE(seams_bench_knn(bench, 1) >= 0);
  REQUIRE(seams_bench_chill(nullptr, 3.5) == -1);
  REQUIRE(seams_bench_chill(bench, 0.0) == -1);
  REQUIRE(seams_bench_chill(bench, 3.5) >= 0);
  seams_bench_close(nullptr);
  REQUIRE(seams_bench_cutoff(bench, 3.5) >= 0);
  seams_bench_close(bench);

  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  const int counts[1] = {1};
  const int flat[1] = {1};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, counts, flat) == SEAMS_INSITU_OK);
  auto frame = ortho(1, 2, xyz);
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const auto dir = scratch("not-a-file");
  std::filesystem::create_directories(dir);
  REQUIRE(seams_insitu_append_dump(&frame, dir.c_str()) == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(eng);
}

TEST_CASE("rejected arguments on the public entry points", "[insitu]") {
  REQUIRE(seams_insitu_append_dump(nullptr, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_append_dump(nullptr, "") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_set_neighbours(nullptr, 1, nullptr, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(nullptr, nullptr, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_pack_owned(-1, nullptr, 1, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr) == -1);
  REQUIRE(seams_insitu_pack_owned(0, nullptr, 1, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr) == 0);
  double xyz[3] = {1, 2, 3};
  double xyzOut[3];
  int typeOut = 0;
  int idOut = 0;
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, nullptr, nullptr, nullptr, &xyzOut[0], &typeOut,
                                  &idOut) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, xyz, nullptr, nullptr, xyzOut, &typeOut, &idOut) ==
          1);
  REQUIRE(typeOut == 1);
  REQUIRE(idOut == 1);
  int total = 0;
  REQUIRE(seams_insitu_gather_displs(-1, nullptr, nullptr, nullptr) == 1);
  REQUIRE(seams_insitu_gather_displs(1, nullptr, nullptr, &total) == 1);
  const int counts[1] = {-1};
  int displs[1] = {0};
  REQUIRE(seams_insitu_gather_displs(1, counts, displs, &total) == 1);
  const int okCounts[2] = {1, 2};
  int okDispls[2] = {0, 0};
  REQUIRE(seams_insitu_gather_displs(2, okCounts, okDispls, &total) == 0);
  REQUIRE(total == 3);
  seams_insitu_scale_nm(nullptr, 1, nullptr);
  double nm[3] = {1, 2, 3};
  double ang[3] = {0, 0, 0};
  seams_insitu_scale_nm(nm, 1, ang);
  REQUIRE(ang[0] == 10);
  REQUIRE(seams_insitu_box_from_rows(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                     nullptr, nullptr, nullptr) == 1);
  double rows[9] = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  rows[0] = 2;
  rows[1] = 1;
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  int nh = 0, nd = 0, nb = 0, nw = 0;
  seams_insitu_tally_cages(1, nullptr, 1, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                           nullptr);
  seams_insitu_tally_cages(1, nullptr, 1, nullptr, nullptr, nullptr, &nh, &nd, &nb, &nw);
  const int types[2] = {1, 2};
  const int hc[2] = {1, 1};
  const int ddc[2] = {1, 0};
  int codes[2] = {0, 0};
  seams_insitu_tally_cages(2, types, 1, hc, ddc, codes, &nh, &nd, &nb, &nw);
  REQUIRE(codes[0] == 3);
  REQUIRE(codes[1] == 0);
}

TEST_CASE("the other side of each public check", "[insitu]") {
  const double pairXyz[6] = {0, 0, 0, 3, 0, 0};
  const int pairIds[2] = {-5, 2};
  const int pairTypes[2] = {1, 1};
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, 2, pairXyz);
  frame.ids = pairIds;
  frame.types = pairTypes;
  frame.coord_scale = 10;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  frame.lx = 21;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  frame.ly = 21;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  frame.lz = 21;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const int wideCounts[1] = {3};
  const int wideFlat[3] = {2, -1, 100000};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, wideCounts, wideFlat) == SEAMS_INSITU_OK);
  frame.step = 2;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const int missingCounts[1] = {1};
  const int missingFlat[1] = {50};
  const double moved[6] = {0.4, 0.2, 0, 3.2, 0.1, 0};
  frame.xyz = moved;
  REQUIRE(seams_insitu_set_neighbours(eng, 1, missingCounts, missingFlat) == SEAMS_INSITU_OK);
  frame.step = 3;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_set_neighbours(eng, 1, missingCounts, nullptr) == SEAMS_INSITU_ERR);
  seams_insitu_destroy(eng);

  opt.debug_dump = 1;
  auto *dumping = seams_insitu_create(&opt);
  frame = ortho(4, 2, pairXyz);
  frame.ids = pairIds;
  frame.types = pairTypes;
  REQUIRE(seams_insitu_observe(dumping, &frame, "") == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_append_dump(&frame, "/dev/full") == SEAMS_INSITU_ERR_IO);
  const auto dump = scratch("one.lammpstrj");
  std::remove(dump.c_str());
  frame.ids = nullptr;
  frame.types = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, dump.c_str()) == SEAMS_INSITU_OK);
  frame.ids = pairIds;
  frame.types = pairTypes;
  REQUIRE(seams_insitu_append_dump(&frame, dump.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_replay(dumping, dump.c_str(), nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(dumping, dump.c_str(), "/dev/full") == SEAMS_INSITU_ERR_IO);
  const auto truncated = scratch("truncated.lammpstrj");
  {
    std::ofstream out(truncated);
    out << "ITEM: TIMESTEP";
  }
  REQUIRE(seams_insitu_replay(dumping, truncated.c_str(), dump.c_str()) == SEAMS_INSITU_ERR_IO);
  seams_insitu_options slow;
  seams_insitu_options_default(&slow);
  slow.cadence.every = 2;
  auto *skipping = seams_insitu_create(&slow);
  const auto odd = scratch("odd.lammpstrj");
  std::remove(odd.c_str());
  frame.step = 1;
  frame.ids = nullptr;
  frame.types = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, odd.c_str()) == SEAMS_INSITU_OK);
  const auto skipped = scratch("skipped.txt");
  REQUIRE(seams_insitu_replay(skipping, odd.c_str(), skipped.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(readFile(skipped).empty());
  seams_insitu_destroy(skipping);
  seams_insitu_destroy(dumping);

  double lattice[9] = {2, 0, 0, 0, 2, 0, 0, 0, 2};
  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  REQUIRE(seams_insitu_box_from_rows(lattice, nullptr, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, nullptr, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, nullptr, &lx, &ly, &lz, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, nullptr, &ly, &lz, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, nullptr, &lz, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, nullptr, &xy, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, nullptr, &xz, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, nullptr, &yz) ==
          1);
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, nullptr) ==
          1);
  lattice[4] = 0;
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  lattice[4] = 2;
  lattice[8] = 0;
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  lattice[8] = 2;
  lattice[2] = 1;
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  lattice[2] = 0;
  lattice[5] = 1;
  REQUIRE(seams_insitu_box_from_rows(lattice, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);

  int hcCount = 0, ddcCount = 0, both = 0, water = 0;
  const int oneType[1] = {1};
  const int oneHc[1] = {1};
  const int oneDdc[1] = {0};
  int oneCode[1] = {0};
  seams_insitu_tally_cages(1, oneType, 1, oneHc, oneDdc, oneCode, &hcCount, nullptr, &both, &water);
  seams_insitu_tally_cages(1, oneType, 1, oneHc, oneDdc, oneCode, &hcCount, &ddcCount, nullptr,
                           &water);
  seams_insitu_tally_cages(1, oneType, 1, oneHc, oneDdc, oneCode, &hcCount, &ddcCount, &both,
                           nullptr);
  seams_insitu_tally_cages(1, nullptr, 1, oneHc, oneDdc, oneCode, &hcCount, &ddcCount, &both,
                           &water);
  seams_insitu_tally_cages(1, oneType, 1, nullptr, oneDdc, oneCode, &hcCount, &ddcCount, &both,
                           &water);
  seams_insitu_tally_cages(1, oneType, 1, oneHc, nullptr, oneCode, &hcCount, &ddcCount, &both,
                           &water);
  seams_insitu_tally_cages(1, oneType, 1, oneHc, oneDdc, nullptr, &hcCount, &ddcCount, &both,
                           &water);

  double owned[3] = {1, 2, 3};
  double ownedOut[3];
  int typeOut = 0;
  int idOut = 0;
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, owned, nullptr, nullptr, nullptr, &typeOut,
                                  &idOut) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, owned, nullptr, nullptr, ownedOut, nullptr,
                                  &idOut) == -1);
  REQUIRE(seams_insitu_pack_owned(1, nullptr, 1, owned, nullptr, nullptr, ownedOut, &typeOut,
                                  nullptr) == -1);
  int total = 0;
  const int rankCounts[1] = {1};
  REQUIRE(seams_insitu_gather_displs(1, rankCounts, nullptr, &total) == 1);

  const double benchXyz[6] = {0, 0, 0, 1, 0, 0};
  const double benchBox[3] = {10, 10, 10};
  REQUIRE(seams_bench_open(2, nullptr, benchBox) == nullptr);
  REQUIRE(seams_bench_open(2, benchXyz, nullptr) == nullptr);
  const double zeroBox[3] = {0, 0, 0};
  REQUIRE(seams_bench_open(2, benchXyz, zeroBox) == nullptr);

  const char *exe = std::getenv("SEAMS_INSITU_EXE");
  REQUIRE(exe != nullptr);
  const auto cli = scratch("trail.txt");
  REQUIRE(run(std::string(exe) + " replay " + dump + " " + cli + " 3.5x") == 2);

  std::remove(dump.c_str());
  std::remove(truncated.c_str());
  std::remove(odd.c_str());
  std::remove(skipped.c_str());
}

TEST_CASE("a wide frame takes the parallel neighbour guard", "[insitu]") {
  constexpr int n = 256;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  for (int i = 0; i < n; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(i % 16);
    xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>((i / 16) % 16);
    xyz[static_cast<std::size_t>(3 * i + 2)] = 0;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 1.2;
  opt.ring_depth = 3;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, n, xyz.data());
  frame.lx = 40;
  frame.ly = 40;
  frame.lz = 40;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  seams_insitu_destroy(eng);
}

TEST_CASE("a nested parallel region skips the inner neighbour team", "[insitu]") {
#ifdef _OPENMP
  constexpr int n = 256;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3, 0.0);
  for (int i = 0; i < n; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(i % 16);
    xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>(i / 16);
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 1.2;
  opt.ring_depth = 3;
  int rc = SEAMS_INSITU_ERR;
#pragma omp parallel num_threads(2)
  {
#pragma omp single
    {
      auto *eng = seams_insitu_create(&opt);
      auto frame = ortho(1, n, xyz.data());
      frame.lx = 40;
      frame.ly = 40;
      frame.lz = 40;
      rc = seams_insitu_observe(eng, &frame, nullptr);
      seams_insitu_destroy(eng);
    }
  }
  REQUIRE(rc == SEAMS_INSITU_OK);
#else
  REQUIRE(true);
#endif
}

TEST_CASE("a long record fails inside fwrite on a full device", "[insitu]") {
  constexpr int n = 1500;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3, 0.0);
  for (int i = 0; i < n; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(i % 20) * 3.0;
    xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>((i / 20) % 20) * 3.0;
    xyz[static_cast<std::size_t>(3 * i + 2)] = static_cast<double>(i / 400) * 3.0;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 1.1;
  opt.ring_depth = 3;
  opt.write_labels = 1;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  REQUIRE(seams_insitu_open_record(eng, "/dev/full") == SEAMS_INSITU_OK);
  auto frame = ortho(1, n, xyz.data());
  frame.lx = 80;
  frame.ly = 80;
  frame.lz = 80;
  const int rc = seams_insitu_observe(eng, &frame, nullptr);
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(result->text_len > 8192);
  REQUIRE(rc == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(eng);
}

TEST_CASE("bench chill counts labels for two atoms", "[insitu]") {
  const double xyz[6] = {0, 0, 0, 1, 0, 0};
  const double box[3] = {10, 10, 10};
  auto *bench = seams_bench_open(2, xyz, box);
  REQUIRE(bench != nullptr);
  REQUIRE(seams_bench_chill(bench, 0.0) == -1);
  const int labels = seams_bench_chill(bench, 3.5);
  REQUIRE(labels == 2);
  seams_bench_close(bench);
}

TEST_CASE("a diamond cell sets a cubic cage flag", "[insitu]") {
  // Conventional diamond cell, 6.4 A, two cells on a side.
  const double basis[8][3] = {{0, 0, 0},
                              {0.5, 0.5, 0},
                              {0.5, 0, 0.5},
                              {0, 0.5, 0.5},
                              {0.25, 0.25, 0.25},
                              {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75},
                              {0.25, 0.75, 0.75}};
  constexpr int cells = 2;
  constexpr double scale = 6.4;
  constexpr int n = cells * cells * cells * 8;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  int cursor = 0;
  for (int iz = 0; iz < cells; ++iz) {
    for (int iy = 0; iy < cells; ++iy) {
      for (int ix = 0; ix < cells; ++ix) {
        for (const auto &site : basis) {
          xyz[static_cast<std::size_t>(cursor++)] = (ix + site[0]) * scale;
          xyz[static_cast<std::size_t>(cursor++)] = (iy + site[1]) * scale;
          xyz[static_cast<std::size_t>(cursor++)] = (iz + site[2]) * scale;
        }
      }
    }
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 3.5;
  opt.ring_depth = 6;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, n, xyz.data());
  frame.lx = cells * scale;
  frame.ly = cells * scale;
  frame.lz = cells * scale;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(result->n == n);
  REQUIRE(result->n_ddc > 0);
  seams_insitu_destroy(eng);
}

TEST_CASE("hexagonal ice sets a hexagonal cage flag", "[insitu]") {
  // Lonsdaleite oxygen network, two cells on each side.
  const double a = 2.76 * std::sqrt(8.0 / 3.0);
  const double c = a * std::sqrt(8.0 / 3.0);
  const double b = a * std::sqrt(3.0);
  const double basis[4][3] = {{1.0 / 3, 2.0 / 3, 0.0},
                              {2.0 / 3, 1.0 / 3, 0.5},
                              {1.0 / 3, 2.0 / 3, 3.0 / 8},
                              {2.0 / 3, 1.0 / 3, 7.0 / 8}};
  constexpr int ncell = 2;
  std::vector<double> xyz;
  for (int ix = 0; ix < ncell; ++ix) {
    for (int iy = 0; iy < ncell; ++iy) {
      for (int iz = 0; iz < ncell; ++iz) {
        for (const double shift : {0.0, 0.5}) {
          for (const auto &site : basis) {
            const double x = std::fmod(shift * a + site[0] * a - site[1] * a / 2 + ix * a, ncell * a);
            const double y = std::fmod(shift * b + site[1] * b / 2 + iy * b, ncell * b);
            const double z = (site[2] + iz) * c;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
          }
        }
      }
    }
  }
  const int n = static_cast<int>(xyz.size() / 3);
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 3.2;
  opt.ring_depth = 6;
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  auto frame = ortho(1, n, xyz.data());
  frame.lx = ncell * a;
  frame.ly = ncell * b;
  frame.lz = ncell * c;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  const auto *result = seams_insitu_last(eng);
  REQUIRE(result != nullptr);
  REQUIRE(result->n_hc > 0);
  seams_insitu_destroy(eng);
}



