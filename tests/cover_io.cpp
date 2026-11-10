// Copyright (c) d-SEAMS developers

#include <catch2/catch_test_macros.hpp>

#include <bond.hpp>
#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_insitu.h>
#include <seams_output.hpp>
#include <selection.hpp>

#ifdef SEAMS_HAS_CHEMFILES
#include <chemfiles.hpp>
#endif

#include <arpa/inet.h>
#include <netinet/in.h>
#include <poll.h>
#include <spawn.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <new>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

// Five coverage profiles run at once. A fixed temp directory is shared
// across those processes, so each process keeps its own tree.
fs::path coverTemp(const char *name) {
  auto dir = fs::temp_directory_path() /
             (std::string(name) + "-" +
              std::to_string(static_cast<long long>(::getpid())));
  fs::create_directories(dir);
  return dir;
}

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

namespace {

std::atomic<int> g_alloc_budget{-1};

struct PhaseEnv {
  PhaseEnv() { ::setenv("DSEAMS_PHASE", "1", 0); }
};

PhaseEnv g_phase_env;

struct AllocLimit {
  explicit AllocLimit(int budget) { g_alloc_budget.store(budget); }
  ~AllocLimit() { g_alloc_budget.store(-1); }
};

struct Cwd {
  fs::path saved;
  explicit Cwd(const fs::path &next) : saved(fs::current_path()) {
    fs::current_path(next);
  }
  ~Cwd() { fs::current_path(saved); }
};

fs::path exePath() { return fs::read_symlink("/proc/self/exe"); }

fs::path repoRoot() {
  fs::path cursor = exePath().parent_path();
  for (int hop = 0; hop < 5; ++hop) {
    if (fs::exists(cursor / "input" / "traj")) {
      return cursor;
    }
    cursor = cursor.parent_path();
  }
  return fs::current_path();
}

fs::path seamsBinary() {
  if (const char *env = std::getenv("SEAMS_CLI")) {
    return env;
  }
  const fs::path beside = exePath().parent_path().parent_path() / "seams";
  if (fs::exists(beside)) {
    return beside;
  }
  return "seams";
}

fs::path langLibrary() {
  if (const char *env = std::getenv("SEAMS_LANG")) {
    return env;
  }
  if (const char *cli = std::getenv("SEAMS_CLI")) {
    return fs::path(cli).parent_path() / "bindings" / "libseams_lang.so";
  }
  return exePath().parent_path().parent_path() / "bindings" / "libseams_lang.so";
}

Cloud makeCloud(int n, int mol = 0) {
  Cloud cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  cloud.nop = n;
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.atomID = i + 1;
    point.molID = mol == 0 ? i + 1 : mol;
    point.type = 1;
    point.x = 0.4 * i;
    point.y = 0.2 * i;
    point.z = 0.1 * i;
    point.inSlice = true;
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

std::string readFile(const fs::path &path) {
  std::ifstream in(path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

void writeText(const fs::path &path, const std::string &text) {
  std::ofstream out(path);
  out << text;
}

struct Proc {
  int status = -1;
  std::string out;
  std::string err;
  bool exited = false;
};

Proc runSeams(const std::vector<std::string> &args) {
  int outPipe[2];
  int errPipe[2];
  REQUIRE(::pipe(outPipe) == 0);
  REQUIRE(::pipe(errPipe) == 0);
  const std::string bin = seamsBinary().string();
  posix_spawn_file_actions_t actions;
  REQUIRE(::posix_spawn_file_actions_init(&actions) == 0);
  REQUIRE(::posix_spawn_file_actions_adddup2(&actions, outPipe[1], STDOUT_FILENO) == 0);
  REQUIRE(::posix_spawn_file_actions_adddup2(&actions, errPipe[1], STDERR_FILENO) == 0);
  REQUIRE(::posix_spawn_file_actions_addclose(&actions, outPipe[0]) == 0);
  REQUIRE(::posix_spawn_file_actions_addclose(&actions, outPipe[1]) == 0);
  REQUIRE(::posix_spawn_file_actions_addclose(&actions, errPipe[0]) == 0);
  REQUIRE(::posix_spawn_file_actions_addclose(&actions, errPipe[1]) == 0);
  std::vector<char *> argv;
  argv.push_back(const_cast<char *>(bin.c_str()));
  for (const auto &arg : args) {
    argv.push_back(const_cast<char *>(arg.c_str()));
  }
  argv.push_back(nullptr);
  std::vector<std::string> envStore;
  for (char **row = environ; row != nullptr && *row != nullptr; ++row) {
    const std::string item(*row);
    if (item.rfind("DSEAMS_PHASE=", 0) == 0 || item.rfind("OMP_NUM_THREADS=", 0) == 0) {
      continue;
    }
    envStore.emplace_back(*row);
  }
  envStore.emplace_back("OMP_NUM_THREADS=1");
  std::vector<char *> envp;
  envp.reserve(envStore.size() + 1);
  for (auto &row : envStore) {
    envp.push_back(row.data());
  }
  envp.push_back(nullptr);
  pid_t pid = 0;
  const int spawned = ::posix_spawn(&pid, bin.c_str(), &actions, nullptr, argv.data(), envp.data());
  ::posix_spawn_file_actions_destroy(&actions);
  ::close(outPipe[1]);
  ::close(errPipe[1]);
  REQUIRE(spawned == 0);
  Proc proc;
  pollfd fds[2] = {{outPipe[0], POLLIN, 0}, {errPipe[0], POLLIN, 0}};
  bool openOut = true;
  bool openErr = true;
  char buffer[512];
  while (openOut || openErr) {
    const int ready = ::poll(fds, 2, -1);
    REQUIRE(ready >= 0);
    for (int slot = 0; slot < 2; ++slot) {
      if ((fds[slot].revents & (POLLIN | POLLHUP | POLLERR)) == 0) {
        continue;
      }
      const ssize_t got = ::read(fds[slot].fd, buffer, sizeof(buffer));
      if (got > 0) {
        (slot == 0 ? proc.out : proc.err).append(buffer, static_cast<std::size_t>(got));
      } else {
        if (slot == 0) {
          openOut = false;
        } else {
          openErr = false;
        }
        fds[slot].fd = -1;
      }
    }
  }
  ::close(outPipe[0]);
  ::close(errPipe[0]);
  int waitStatus = 0;
  REQUIRE(::waitpid(pid, &waitStatus, 0) == pid);
  proc.exited = WIFEXITED(waitStatus) != 0;
  proc.status = proc.exited ? WEXITSTATUS(waitStatus) : 128 + WTERMSIG(waitStatus);
  return proc;
}

std::string dumpText(int nAtoms, const std::string &rows, int step = 0) {
  std::ostringstream out;
  out << "ITEM: TIMESTEP\n"
      << step << "\nITEM: NUMBER OF ATOMS\n"
      << nAtoms << "\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS id type x y z\n"
      << rows;
  return out.str();
}

seams_insitu_frame frameOf(int n, const double *xyz, const int *types, const int *ids,
                           int step = 0) {
  seams_insitu_frame frame{};
  frame.step = step;
  frame.n = n;
  frame.xyz = xyz;
  frame.types = types;
  frame.ids = ids;
  frame.xlo = 0;
  frame.ylo = 0;
  frame.zlo = 0;
  frame.lx = 10;
  frame.ly = 10;
  frame.lz = 10;
  frame.gathered = 1;
  return frame;
}

int listenLocal(uint16_t *port) {
  const int fd = ::socket(AF_INET, SOCK_STREAM, 0);
  REQUIRE(fd >= 0);
  sockaddr_in addr{};
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
  REQUIRE(::bind(fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) == 0);
  REQUIRE(::listen(fd, 1) == 0);
  socklen_t len = sizeof(addr);
  REQUIRE(::getsockname(fd, reinterpret_cast<sockaddr *>(&addr), &len) == 0);
  *port = ntohs(addr.sin_port);
  return fd;
}

void sendAll(int fd, const void *data, std::size_t bytes) {
  const auto *in = static_cast<const unsigned char *>(data);
  std::size_t sent = 0;
  while (sent < bytes) {
    const ssize_t n = ::send(fd, in + sent, bytes - sent, 0);
    REQUIRE(n > 0);
    sent += static_cast<std::size_t>(n);
  }
}

void sendHeader(int fd, int kind, int length) {
  const std::uint32_t words[2] = {htonl(static_cast<std::uint32_t>(kind)),
                                  htonl(static_cast<std::uint32_t>(length))};
  sendAll(fd, words, sizeof(words));
}

// The handshake compares the version word as a host integer. Later
// headers swap both words.
void sendHandshake(int fd, int version) {
  const std::uint32_t words[2] = {htonl(4), static_cast<std::uint32_t>(version)};
  sendAll(fd, words, sizeof(words));
}

enum class ImdScript { handshakeBad, close, energyCoord, zeroCoord, shortCoord, negative, badKind };

void imdServer(int listenFd, ImdScript script) {
  const int fd = ::accept(listenFd, nullptr, nullptr);
  if (fd < 0) {
    return;
  }
  if (script == ImdScript::handshakeBad) {
    sendHandshake(fd, 9);
    ::close(fd);
    return;
  }
  sendHandshake(fd, 2);
  unsigned char go[8];
  std::size_t got = 0;
  while (got < sizeof(go)) {
    const ssize_t n = ::recv(fd, go + got, sizeof(go) - got, 0);
    if (n <= 0) {
      ::close(fd);
      return;
    }
    got += static_cast<std::size_t>(n);
  }
  if (script == ImdScript::close) {
    ::close(fd);
    return;
  }
  if (script == ImdScript::energyCoord) {
    sendHeader(fd, 1, 40);
    unsigned char body[40] = {};
    sendAll(fd, body, sizeof(body));
    sendHeader(fd, 2, 1);
    const float xyz[3] = {1.f, 2.f, 3.f};
    sendAll(fd, xyz, sizeof(xyz));
  } else if (script == ImdScript::zeroCoord) {
    sendHeader(fd, 2, 0);
  } else if (script == ImdScript::shortCoord) {
    sendHeader(fd, 2, 1);
  } else if (script == ImdScript::negative) {
    sendHeader(fd, 2, -1);
  } else if (script == ImdScript::badKind) {
    sendHeader(fd, 9, 0);
  }
  ::close(fd);
}

int imdSession(ImdScript script, int maxFrames, const fs::path &dir) {
  uint16_t port = 0;
  const int listenFd = listenLocal(&port);
  std::thread server(imdServer, listenFd, script);
  const double box[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
  const auto dump = (dir / "imd.dump").string();
  const auto record = (dir / "imd.rec").string();
  const int rc = seams_insitu_imd_session("127.0.0.1", port, dump.c_str(), record.c_str(), box,
                                          maxFrames);
  server.join();
  ::close(listenFd);
  return rc;
}

struct Lang {
  void *handle = nullptr;
  const char *(*version)() = nullptr;
  const char *(*last_error)() = nullptr;
  int (*nframes)(const char *) = nullptr;
  void *(*read)(const char *, int, int) = nullptr;
  void *(*from_arrays)(const double *, int, const double *, const int *, const int *) = nullptr;
  void (*free_frame)(void *) = nullptr;
  int (*nop)(const void *) = nullptr;
  int (*index)(const void *) = nullptr;
  void (*box)(const void *, double *) = nullptr;
  const double *(*positions)(const void *) = nullptr;
  const int *(*types)(const void *) = nullptr;
  const int *(*ids)(const void *) = nullptr;
  int (*chill)(void *, double, int, const char *, int *, int *) = nullptr;
  int (*chill_plus)(void *, double, int, const char *, int *, int *) = nullptr;
  int (*rings)(void *, double, int, const char *, int, int *, int) = nullptr;
  int (*steinhardt)(void *, double, int, int, double *, double *) = nullptr;
  int (*cages)(void *, double, int, int, const char *, int, int *, int *, int *, int *,
               int *) = nullptr;
  int (*signature)(void *, double, int, const char *, const char *, int *, int *) = nullptr;
  int (*signature_edge)(void *, double, int, const char *, const char *, const char *, int *,
                        int *) = nullptr;
  int (*rdf_nbins)(double, int) = nullptr;
  int (*rdf)(const void *, double, int, int, int, double *, double *, int *, double *) = nullptr;
  int (*rdf2d_nbin)(double, double) = nullptr;
  int (*rdf2d)(void *, double, double, const char *, double *, double *) = nullptr;
  int (*prisms)(void *, double, int, int, const char *, void *, int *) = nullptr;
  int (*monolayer)(void *, double, int, int, double, const char *, int *, int) = nullptr;
};

template <typename T> void loadSym(Lang &lang, T &fn, const char *name) {
  fn = reinterpret_cast<T>(::dlsym(lang.handle, name));
  REQUIRE(fn != nullptr);
}

Lang &langApi() {
  static Lang lang = [] {
    Lang out;
    out.handle = ::dlopen(langLibrary().c_str(), RTLD_NOW);
    REQUIRE(out.handle != nullptr);
    loadSym(out, out.version, "seams_lang_version");
    loadSym(out, out.last_error, "seams_last_error");
    loadSym(out, out.nframes, "seams_path_nframes");
    loadSym(out, out.read, "seams_frame_read");
    loadSym(out, out.from_arrays, "seams_frame_from_arrays");
    loadSym(out, out.free_frame, "seams_frame_free");
    loadSym(out, out.nop, "seams_frame_nop");
    loadSym(out, out.index, "seams_frame_index");
    loadSym(out, out.box, "seams_frame_box");
    loadSym(out, out.positions, "seams_frame_positions");
    loadSym(out, out.types, "seams_frame_types");
    loadSym(out, out.ids, "seams_frame_ids");
    loadSym(out, out.chill, "seams_frame_chill");
    loadSym(out, out.chill_plus, "seams_frame_chill_plus");
    loadSym(out, out.rings, "seams_frame_rings");
    loadSym(out, out.steinhardt, "seams_frame_steinhardt");
    loadSym(out, out.cages, "seams_frame_cages");
    loadSym(out, out.signature, "seams_frame_signature");
    loadSym(out, out.signature_edge, "seams_frame_signature_edge");
    loadSym(out, out.rdf_nbins, "seams_frame_rdf_nbins");
    loadSym(out, out.rdf, "seams_frame_rdf");
    loadSym(out, out.rdf2d_nbin, "seams_frame_rdf2d_nbin");
    loadSym(out, out.rdf2d, "seams_frame_rdf2d");
    loadSym(out, out.prisms, "seams_frame_prisms");
    loadSym(out, out.monolayer, "seams_frame_monolayer");
    return out;
  }();
  return lang;
}

bool threwBadAlloc(int budget, const std::function<void()> &call) {
  try {
    AllocLimit limit(budget);
    call();
    return false;
  } catch (const std::bad_alloc &) {
    return true;
  }
}

} // namespace

void *operator new(std::size_t bytes) {
  const int budget = g_alloc_budget.load();
  if (budget == 0) {
    g_alloc_budget.store(-1);
    throw std::bad_alloc();
  }
  if (budget > 0) {
    g_alloc_budget.store(budget - 1);
  }
  if (void *memory = std::malloc(bytes == 0 ? 1 : bytes)) {
    return memory;
  }
  throw std::bad_alloc();
}

void *operator new[](std::size_t bytes) { return ::operator new(bytes); }

void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

TEST_CASE("map builders throw through the return frame", "[cover_io]") {
  auto cloud = makeCloud(2);
  bool sawMap = false;
  bool sawMulti = false;
  bool sawHydrogens = false;
  for (int budget = 0; budget < 8; ++budget) {
    sawMap = sawMap || threwBadAlloc(budget, [&] { auto map = molSys::createIDMolIDmap(cloud); REQUIRE(map.size() == 2); });
    sawMulti = sawMulti || threwBadAlloc(budget, [&] {
      auto map = molSys::createMolIDAtomIDMultiMap(cloud);
      REQUIRE(map.size() == 2);
    });
    auto hydrogens = makeCloud(2, 7);
    auto oxygens = makeCloud(1, 7);
    sawHydrogens = sawHydrogens || threwBadAlloc(budget, [&] {
      auto rows = molSys::hAtomMolList(hydrogens, oxygens);
      REQUIRE(rows.size() == 1);
    });
  }
  REQUIRE(sawMap);
  REQUIRE(sawMulti);
  REQUIRE(sawHydrogens);
  auto map = molSys::createIDMolIDmap(cloud);
  REQUIRE(map[1] == 1);
  REQUIRE(map[2] == 2);
}

TEST_CASE("slice selection keeps mates and a second cloud", "[cover_io]") {
  auto cloud = makeCloud(3);
  cloud.pts[0].molID = 4;
  cloud.pts[1].molID = 4;
  cloud.pts[2].molID = 9;
  cloud.pts[0].x = cloud.pts[0].y = cloud.pts[0].z = 1.0;
  cloud.pts[1].x = cloud.pts[1].y = cloud.pts[1].z = 8.0;
  cloud.pts[2].x = cloud.pts[2].y = cloud.pts[2].z = 8.0;
  cloud.pts[1].inSlice = true;
  const std::array<double, 3> low{0, 0, 0};
  const std::array<double, 3> high{2, 2, 2};
  gen::moleculesInSingleSlice(cloud, false, low, high);
  REQUIRE(cloud.pts[0].inSlice);
  REQUIRE(cloud.pts[1].inSlice);
  REQUIRE(cloud.pts[2].inSlice);
  gen::moleculesInSingleSlice(cloud, true, low, high);
  REQUIRE_FALSE(cloud.pts[2].inSlice);

  auto oxygen = cloud;
  auto full = cloud;
  full.pts.push_back(cloud.pts[0]);
  full.pts.back().atomID = 40;
  full.pts.back().molID = 4;
  full.nop = 4;
  full.idIndexMap[40] = 3;
  const std::vector<std::vector<int>> rings{{0, 1}, {2}};
  ring::getEdgeMoleculesInRings(rings, oxygen, full, low, high, false);
  REQUIRE(full.pts[0].inSlice);
  REQUIRE(full.pts[3].inSlice);
  REQUIRE(oxygen.pts[1].inSlice);
}

TEST_CASE("output paths and padded data files", "[cover_io]") {
  const auto root = coverTemp("dseams-cover-io");
  fs::remove_all(root);
  const auto work = root / "work";
  fs::create_directories(work);
  fs::create_directories(root / "output" / "prisms");
  fs::create_directories(root / "output" / "cages" / "hexCages");
  fs::create_directories(root / "output" / "cages" / "doubleDiaCages");
  fs::create_directories(root / "output" / "cages" / "hexBasalRings");
  fs::create_directories(root / "output" / "deformed" / "basalRings");
  fs::create_directories(root / "output" / "perfect" / "basalRings");
  Cwd guard(work);

  const auto kept = root / "already";
  fs::create_directories(kept);
  REQUIRE(sout::makePath(kept.string()) == 0);
  REQUIRE(sout::makePath((root / "fresh" / "nested").string()) == 0);
  REQUIRE(fs::is_directory(root / "fresh" / "nested"));
  const auto blocking = root / "blocking";
  writeText(blocking, "not a directory\n");
  REQUIRE(sout::makePath(blocking.string()) == 1);
  const auto raced = root / "race-dir";
  std::atomic<bool> stopRace{false};
  std::thread racer([&] {
    while (!stopRace.load()) {
      std::error_code ignore;
      fs::create_directory(raced, ignore);
    }
  });
  int sawDirectory = 0;
  int lastCode = -1;
  for (int attempt = 0; attempt < 20000; ++attempt) {
    std::error_code ignore;
    fs::remove(raced, ignore);
    lastCode = sout::makePath(raced.string());
    sawDirectory += lastCode == 0 ? 1 : 0;
  }
  stopRace.store(true);
  racer.join();
  REQUIRE((lastCode == 0 || lastCode == 1));
  REQUIRE(sawDirectory > 0);

  auto cloud = makeCloud(4);
  cloud.box = {10, 10, 10, 0.2, 0.1, 0.0};
  std::vector<std::vector<int>> rings{{1, 3}};
  std::vector<std::vector<int>> bonds{{1, 3}};
  REQUIRE(sout::writeLAMMPSdata(cloud, {}, bonds, "empty.data") == 1);
  REQUIRE(sout::writeLAMMPSdata(cloud, rings, bonds, "gapped.data") == 0);
  const auto gapped = readFile(root / "output" / "gapped.data");
  REQUIRE(gapped.find("2 1.0 # dummy") != std::string::npos);
  REQUIRE(gapped.find("\n2 ") != std::string::npos);

  std::vector<int> basal1{0, 1, 2};
  std::vector<int> basal2{1, 2, 0};
  REQUIRE(sout::writePrisms(basal1, basal2, 1, cloud) == 0);
  REQUIRE(fs::exists(root / "output" / "prisms" / "singleRing.dat"));
  REQUIRE(sout::writePrisms(basal1, basal2, 2, cloud) == 0);
  REQUIRE(fs::file_size(root / "output" / "prisms" / "prism2.dat") > 0);

  std::vector<std::vector<int>> cageRings{{1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};
  auto wide = makeCloud(12);
  std::vector<std::vector<int>> nList(12);
  for (int i = 0; i < 12; ++i) {
    nList[static_cast<std::size_t>(i)] = {i + 1};
  }
  for (int i = 0; i < 6; ++i) {
    nList[static_cast<std::size_t>(i)].push_back(i + 7);
  }
  std::vector<cage::Cage> cages{
      {cage::cageType::HexC, {0, 1}},
      {cage::cageType::DoubleDiaC, {0}},
      {static_cast<cage::cageType>(9), {0}},
  };
  REQUIRE(sout::writeAllCages("ignored", cages, cageRings, nList, wide, 1) == 0);
  REQUIRE(sout::writeAllCages("ignored", cages, cageRings, nList, wide, 1) == 0);
  std::vector<cage::Cage> none;
  REQUIRE(sout::writeAllCages("ignored", none, cageRings, nList, wide, 1) == 1);
  REQUIRE(sout::writeEachCage({0, 1}, 3, static_cast<cage::cageType>(4), cageRings, wide) == 1);
  REQUIRE(sout::writeEachCage({0}, 1, cage::cageType::DoubleDiaC, cageRings, wide) == 0);

  REQUIRE(sout::writeBasalRingsHex({0, 1}, 1, nList, cageRings) == 0);
  std::vector<std::vector<int>> rotated{{1, 2, 3, 4, 5, 6}, {9, 10, 11, 12, 7, 8}};
  nList[0] = {1};
  nList[1] = {2, 8};
  nList[3] = {4, 10};
  REQUIRE(sout::writeBasalRingsHex({0, 1}, 2, nList, rotated) == 0);
  std::vector<std::vector<int>> reversed{{1, 2, 3, 4, 5, 6}, {7, 12, 8, 11, 9, 10}};
  for (int i = 0; i < 6; ++i) {
    nList[static_cast<std::size_t>(i)] = {i + 1, i + 7};
  }
  REQUIRE(sout::writeBasalRingsHex({0, 1}, 3, nList, reversed) == 0);
  nList[2] = {3, 8};
  REQUIRE(sout::writeBasalRingsHex({0, 1}, 4, nList, cageRings) == 1);

  auto prismCloud = makeCloud(8);
  for (int i = 0; i < 8; ++i) {
    prismCloud.pts[static_cast<std::size_t>(i)].atomID = i;
    prismCloud.idIndexMap.clear();
  }
  prismCloud.pts[1].x = prismCloud.pts[1].y = prismCloud.pts[1].z = 0;
  prismCloud.pts[4].x = 0.2;
  prismCloud.pts[4].y = 0;
  prismCloud.pts[4].z = 0;
  prismCloud.pts[6].x = 9;
  prismCloud.pts[6].y = 9;
  prismCloud.pts[6].z = 9;
  prismCloud.pts[7].x = 0.1;
  std::vector<std::vector<int>> prismList(8);
  for (int i = 0; i < 8; ++i) {
    prismList[static_cast<std::size_t>(i)] = {i};
  }
  prismList[1].push_back(7);
  std::vector<int> ringA{0, 1, 2, 3};
  std::vector<int> ringB{4, 5, 6, 7};
  REQUIRE(sout::writeBasalRingsPrism(ringA, ringB, 1, prismList, prismCloud, true) == 0);
  REQUIRE(fs::file_size(root / "output" / "deformed" / "basalRings" / "basalRings1.dat") > 0);
  prismCloud.pts[4].x = 9;
  prismCloud.pts[6].x = 0.2;
  prismCloud.pts[6].y = 0;
  prismCloud.pts[6].z = 0;
  REQUIRE(sout::writeBasalRingsPrism(ringA, ringB, 2, prismList, prismCloud, false) == 0);
  prismCloud.pts[4] = prismCloud.pts[6];
  REQUIRE(sout::writeBasalRingsPrism(ringA, ringB, 3, prismList, prismCloud, true) == 1);
  prismList[1] = {1};
  REQUIRE(sout::writeBasalRingsPrism(ringA, ringB, 4, prismList, prismCloud, false) == 1);

  auto sliced = makeCloud(4);
  sliced.pts[0].molID = 1;
  sliced.pts[1].molID = 2;
  sliced.pts[2].molID = 5;
  sliced.pts[3].molID = 1;
  sliced.pts[3].inSlice = false;
  const auto sliceDir = (root / "slice").string() + "/";
  REQUIRE(sout::writeMoleculeIDsInSlice(sliceDir, sliced) == 0);
  const auto sliceText =
      readFile(root / "slice" / "selection" / "IDtextFiles" / "molID-1.dat");
  REQUIRE(sliceText.find(":2 5") != std::string::npos);
  sliced.pts[1].molID = 3;
  sliced.pts[2].molID = 4;
  sliced.currentFrame = 2;
  REQUIRE(sout::writeMoleculeIDsInSlice(sliceDir, sliced) == 0);
  const auto runText =
      readFile(root / "slice" / "selection" / "IDtextFiles" / "molID-2.dat");
  REQUIRE(runText.find("1 3 4") != std::string::npos);
  REQUIRE(sout::writeMoleculeIDsExpressionSelectOVITO(sliceDir, sliced) == 0);

  const std::vector<int> counts{1, 0};
  const std::vector<double> heights{0.5, 0.25};
  const auto statDir = (root / "stats").string() + "/";
  REQUIRE(sout::writePrismNum(statDir, counts, counts, heights, 4, 3, 1) == 0);
  REQUIRE(sout::writeClusterStats(statDir, 3, 4, 2, 1, 2.5, 1) == 0);
  REQUIRE(sout::writeRingNum(statDir, 3, counts, heights, heights, heights, 4, 1) == 0);
  REQUIRE(sout::writeRingNumBulk(statDir, 3, counts, 4, 1) == 0);
  REQUIRE(sout::writeTopoBulkData(statDir, 3, 1, 1, 0, 0, 0, 1) == 0);
  const std::string later = readFile(root / "stats" / "topoINT" / "nPrisms.dat");
  REQUIRE(later.find("Frame") == std::string::npos);

  std::vector<std::vector<int>> oneRing{{1, 2}};
  std::vector<int> which{0};
  std::vector<std::vector<int>> bondsByIndex{{0, 1, 3}, {1, 0}, {2}, {3}};
  const auto bondFile = root / "bonds.txt";
  writeText(bondFile, "1 Bonds\n1 2\n");
  REQUIRE(sout::writeLAMMPSdataPrisms(cloud, oneRing, true, bondFile.string(), which, bondsByIndex,
                                      "from-file.data") == 0);
  REQUIRE(sout::writeLAMMPSdataPrisms(cloud, oneRing, false, "", which, bondsByIndex,
                                      "from-nlist.data") == 0);
  const auto prismData = readFile(root / "output" / "from-nlist.data");
  REQUIRE(prismData.find("2 1.0 # dummy") != std::string::npos);
  REQUIRE(prismData.find(" 2 ") != std::string::npos);
  REQUIRE(sout::writeLAMMPSdataPrisms(cloud, {}, false, "", which, bondsByIndex) == 1);

  std::vector<std::vector<int>> cageRing{{1, 2, 3, 4, 5, 6}};
  auto cageCloud = makeCloud(7);
  std::vector<cage::Cage> only{{cage::cageType::HexC, {0}}};
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, none, cage::cageType::HexC, 1) == 1);
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, only, cage::cageType::HexC, 0) == 1);
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, only, cage::cageType::HexC, 1,
                                     "hc.data") == 0);
  only[0].type = cage::cageType::DoubleDiaC;
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, only, cage::cageType::DoubleDiaC, 1,
                                     "ddc.data") == 0);
  REQUIRE(readFile(root / "output" / "ddc.data").find("# DDC") != std::string::npos);
  only[0].type = static_cast<cage::cageType>(3);
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, only, only[0].type, 1, "bad.data") == 0);
  REQUIRE(readFile(root / "output" / "bad.data").find("# error") != std::string::npos);
  std::vector<cage::Cage> mixedTypes{{cage::cageType::DoubleDiaC, {0}},
                                     {cage::cageType::HexC, {0}}};
  REQUIRE(sout::writeLAMMPSdataCages(cageCloud, cageRing, mixedTypes, cage::cageType::HexC, 1,
                                     "skip.data") == 0);

  std::vector<cage::iceType> kinds{cage::iceType::dummy, cage::iceType::hc, cage::iceType::ddc,
                                   cage::iceType::mixed, cage::iceType::pnc, cage::iceType::mixed2};
  auto bulk = makeCloud(6);
  std::vector<std::vector<int>> bulkList{{0, 1}, {1, 0}, {2, 3}, {3, 2}, {4}, {5}};
  const auto bulkDir = (root / "bulk").string() + "/";
  REQUIRE(sout::writeLAMMPSdataTopoBulk(bulk, bulkList, kinds, bulkDir, true) == 0);
  REQUIRE(sout::writeLAMMPSdataTopoBulk(bulk, bulkList, kinds, bulkDir, false) == 0);
  REQUIRE(fs::file_size(root / "bulk" / "bulkTopo" / "dataFiles" / "system-1.data") > 0);

  auto histoCloud = makeCloud(2);
  histoCloud.pts[0].type = 2;
  histoCloud.pts[1].type = 1;
  histoCloud.pts[1].c_ij.push_back({molSys::bond_type::staggered, 0.4});
  std::vector<std::vector<int>> histo{{0}, {1, 0}};
  std::vector<double> q6{0.1, 0.2};
  REQUIRE(sout::writeHisto(histoCloud, histo, q6) == 0);
  REQUIRE(readFile(work / "q6.txt").find("0.2") != std::string::npos);

  fs::remove_all(root);
}

#ifdef SEAMS_HAS_CHEMFILES
TEST_CASE("chemfiles names, residues, and frame slices", "[cover_io]") {
  const auto dir = coverTemp("dseams-cover-chem");
  fs::remove_all(dir);
  fs::create_directories(dir);
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D{12, 12, 12}));
  const std::pair<const char *, const char *> named[] = {
      {"Xx", "O"}, {"Xx", "OH"}, {"Xx", "H"}, {"Xx", "HO"}, {"OW", "Q"}, {"HW", "Q"}, {"HO", "Q"},
  };
  for (int i = 0; i < 7; ++i) {
    chemfiles::Atom atom(named[i].second);
    atom.set_type(named[i].first);
    frame.add_atom(atom, {0.1 * i, 0.2 * i, 0.3 * i});
  }
  chemfiles::Residue numbered("WAT", 7);
  numbered.add_atom(0);
  frame.add_residue(numbered);
  chemfiles::Residue bare("SOL");
  bare.add_atom(1);
  frame.add_residue(bare);
  const auto path = dir / "names.pdb";
  chemfiles::Trajectory(path.string(), 'w').write(frame);
  chemfiles::Frame symbols;
  symbols.set_cell(chemfiles::UnitCell(chemfiles::Vector3D{12, 12, 12}));
  for (const char *token : {"OW", "OH", "HW", "HO", "Xx"}) {
    chemfiles::Atom atom(token);
    atom.set_type(token);
    symbols.add_atom(atom, {1.0, 1.0, 1.0});
  }
  const auto symbolPath = dir / "symbols.xyz";
  chemfiles::Trajectory(symbolPath.string(), 'w', "XYZ").write(symbols);

  Cloud cloud;
  sinp::ChemfilesRead request;
  request.frame = 1;
  request.guessOxygen = false;
  request.typeFilter = -1;
  Cloud symbolsCloud;
  symbolsCloud = sinp::readChemfiles(symbolPath.string(), symbolsCloud, request);
  REQUIRE(symbolsCloud.nop == 5);
  bool sawOxygen = false;
  bool sawHydrogen = false;
  bool sawOther = false;
  for (const auto &point : symbolsCloud.pts) {
    sawOxygen = sawOxygen || point.type == 8;
    sawHydrogen = sawHydrogen || point.type == 1;
    sawOther = sawOther || point.type == 0;
  }
  REQUIRE(sawOxygen);
  REQUIRE(sawHydrogen);
  REQUIRE(sawOther);

  request.guessOxygen = true;
  cloud = sinp::readChemfiles(path.string(), cloud, request);
  REQUIRE(cloud.nop >= 1);
  bool oxygen = false;
  for (const auto &point : cloud.pts) {
    oxygen = oxygen || point.type == 8;
  }
  REQUIRE(oxygen);

  request.guessOxygen = false;
  request.typeFilter = -1;
  cloud = sinp::readChemfiles(path.string(), cloud, request);
  REQUIRE(cloud.nop == 7);
  bool hydrogen = false;
  for (const auto &point : cloud.pts) {
    hydrogen = hydrogen || point.type == 1;
  }
  REQUIRE(hydrogen);

  const auto xyz = dir / "range.xyz";
  writeText(xyz, "1\n\nO 0 0 0\n");
  int seen = 0;
  REQUIRE(sinp::forEachChemfilesFrame(xyz.string(), 1, 1, -1, "", -1, 2,
                                      [&](int, Cloud &) { seen += 1; }));
  REQUIRE(seen == 0);
  REQUIRE(sinp::forEachChemfilesFrame(xyz.string(), 1, 1, -1, "", 3, 2,
                                      [&](int, Cloud &) { seen += 1; }));
  REQUIRE(seen == 0);

  const auto dumped = dir / "native.dump";
  writeText(dumped, dumpText(1, "1 2 0.0 0.0 0.0\n"));
  Cloud native = sinp::openFrame(dumped.string(), 1, 0);
  REQUIRE(native.nop >= 0);
  const auto junk = dir / "junk.dump";
  writeText(junk, "not a trajectory\n");
  Cloud rejected = sinp::openFrame(junk.string(), 1, 0);
  REQUIRE(rejected.nop == 0);

  const auto dcd = dir / "movie.dcd";
  writeText(dcd, "not a dcd");
  writeText(dir / "movie.pdb", "END\n");
  Cloud dcdCloud;
  request.guessOxygen = false;
  request.typeFilter = -1;
  dcdCloud = sinp::readChemfiles(dcd.string(), dcdCloud, request);
  REQUIRE(dcdCloud.nop == 0);

  fs::remove_all(dir);
}
#endif

TEST_CASE("in-situ nulls, phase trace, and imd client", "[cover_io]") {
  const auto dir = coverTemp("dseams-cover-insitu");
  fs::remove_all(dir);
  fs::create_directories(dir);

  REQUIRE(seams_insitu_create(nullptr) != nullptr);
  seams_insitu_options bad{};
  seams_insitu_options_default(&bad);
  bad.cutoff = 0;
  REQUIRE(seams_insitu_create(&bad) == nullptr);
  bad.cutoff = 3.5;
  bad.ring_depth = 2;
  REQUIRE(seams_insitu_create(&bad) == nullptr);
  bad.ring_depth = 8;
  bad.write_labels = 2;
  REQUIRE(seams_insitu_create(&bad) == nullptr);
  bad.write_labels = 1;
  bad.debug_dump = 3;
  REQUIRE(seams_insitu_create(&bad) == nullptr);

  seams_insitu_options opt{};
  seams_insitu_options_default(&opt);
  auto *eng = seams_insitu_create(&opt);
  REQUIRE(eng != nullptr);
  REQUIRE(seams_insitu_set_sink(nullptr, nullptr, nullptr) == SEAMS_INSITU_ERR);
  int sinkCalls = 0;
  REQUIRE(seams_insitu_set_sink(
              eng, [](const seams_insitu_result *, void *user) {
                *static_cast<int *>(user) += 1;
                return 0;
              },
              &sinkCalls) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_open_record(nullptr, "x") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, "") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_open_record(eng, (dir / "nope" / "rec.txt").c_str()) == SEAMS_INSITU_ERR_IO);
  const auto record = dir / "rec.txt";
  REQUIRE(seams_insitu_open_record(eng, record.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_open_record(eng, (dir / "rec2.txt").c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_close_record(nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_close_record(eng) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_close_record(eng) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_set_sink(
              eng, [](const seams_insitu_result *, void *user) {
                *static_cast<int *>(user) += 1;
                return 0;
              },
              &sinkCalls) == SEAMS_INSITU_OK);

  REQUIRE(seams_insitu_due(nullptr, 0) == -1);
  seams_insitu_cadence cadence{0, 0, -1};
  REQUIRE(seams_insitu_due(&cadence, 0) == -1);
  cadence = {2, 1, 4};
  REQUIRE(seams_insitu_due(&cadence, 0) == 0);
  REQUIRE(seams_insitu_due(&cadence, 5) == 0);
  REQUIRE(seams_insitu_due(&cadence, 2) == 0);
  REQUIRE(seams_insitu_due(&cadence, 1) == 1);

  const double xyz[6] = {0, 0, 0, 1.2, 0, 0};
  const int types[2] = {1, 1};
  int ids[2] = {1, -5};
  auto frame = frameOf(2, xyz, types, ids, 0);
  REQUIRE(seams_insitu_step(nullptr, &frame) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_step(eng, nullptr) == SEAMS_INSITU_ERR);
  frame.gpu_resident = 1;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_GPU_DEFERRED);
  frame.gpu_resident = 0;
  frame.local_only = 1;
  frame.gathered = 0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_NEED_GATHER);
  frame.local_only = 0;
  frame.gathered = 1;
  frame.n = -1;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR);
  frame.n = 2;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR);
  frame.xyz = xyz;
  frame.lx = 0;
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_ERR_BOX);
  frame.lx = 10;
  const int countsBad[1] = {1};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, countsBad, nullptr) == SEAMS_INSITU_ERR);
  const int flatUnknown[2] = {99, 1};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, countsBad, flatUnknown) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  REQUIRE(sinkCalls >= 1);
  const int flatKnown[3] = {1, 2, 1};
  const int countsRow[1] = {3};
  REQUIRE(seams_insitu_set_neighbours(eng, 1, countsRow, flatKnown) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_step(eng, &frame) == SEAMS_INSITU_OK);
  frame.step = 1;
  frame.lx = 11;
  REQUIRE(seams_insitu_observe(eng, &frame, "") == SEAMS_INSITU_OK);
  frame.lx = 10;
  frame.ly = 11;
  REQUIRE(seams_insitu_observe(eng, &frame, nullptr) == SEAMS_INSITU_OK);
  opt.debug_dump = 1;
  auto *dumping = seams_insitu_create(&opt);
  const auto dumpPath = dir / "live.dump";
  REQUIRE(seams_insitu_observe(dumping, &frame, dumpPath.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_observe(dumping, &frame, (dir).c_str()) == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(dumping);

  frame.n = 0;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, (dir / "empty.dump").c_str()) == SEAMS_INSITU_OK);
  frame.n = 2;
  frame.xyz = xyz;
  REQUIRE(seams_insitu_append_dump(nullptr, "x") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_append_dump(&frame, "") == SEAMS_INSITU_ERR);
  frame.n = -1;
  REQUIRE(seams_insitu_append_dump(&frame, (dir / "x.dump").c_str()) == SEAMS_INSITU_ERR);
  frame.n = 2;
  frame.xyz = nullptr;
  REQUIRE(seams_insitu_append_dump(&frame, (dir / "x.dump").c_str()) == SEAMS_INSITU_ERR);
  frame.xyz = xyz;
  frame.lz = 0;
  REQUIRE(seams_insitu_append_dump(&frame, (dir / "badbox.dump").c_str()) == SEAMS_INSITU_ERR_BOX);
  frame.lz = 10;
  REQUIRE(seams_insitu_append_dump(&frame, dir.c_str()) == SEAMS_INSITU_ERR_IO);

  const auto replayDump = dir / "replay.dump";
  frame.ids = ids;
  frame.step = 0;
  REQUIRE(seams_insitu_append_dump(&frame, replayDump.c_str()) == SEAMS_INSITU_OK);
  frame.step = 3;
  REQUIRE(seams_insitu_append_dump(&frame, replayDump.c_str()) == SEAMS_INSITU_OK);
  REQUIRE(seams_insitu_replay(nullptr, "a", "b") == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_replay(eng, (dir / "missing.dump").c_str(),
                              (dir / "out.txt").c_str()) == SEAMS_INSITU_ERR_IO);
  writeText(dir / "badstep.dump", "ITEM: TIMESTEP\nnope\n");
  REQUIRE(seams_insitu_replay(eng, (dir / "badstep.dump").c_str(),
                              (dir / "out.txt").c_str()) == SEAMS_INSITU_ERR_IO);
  writeText(dir / "nostep.dump", "no timesteps here\n");
  REQUIRE(seams_insitu_replay(eng, (dir / "nostep.dump").c_str(),
                              (dir / "out.txt").c_str()) == SEAMS_INSITU_ERR);
  opt.cadence.every = 2;
  opt.cadence.first = 0;
  opt.cadence.last = 10;
  auto *paced = seams_insitu_create(&opt);
  REQUIRE(seams_insitu_replay(paced, replayDump.c_str(), (dir / "paced.txt").c_str()) ==
          SEAMS_INSITU_OK);
  REQUIRE(fs::file_size(dir / "paced.txt") > 0);
  REQUIRE(seams_insitu_replay(paced, replayDump.c_str(), dir.c_str()) == SEAMS_INSITU_ERR_IO);
  seams_insitu_destroy(paced);
  seams_insitu_destroy(eng);
  seams_insitu_destroy(nullptr);
  REQUIRE(seams_insitu_status(nullptr) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_last(nullptr) == nullptr);

  const double one[3] = {0, 0, 0};
  const int oneType[1] = {1};
  const int oneId[1] = {7};
  auto single = frameOf(1, one, oneType, oneId, 4);
  const auto kept = dir / "one.dump";
  REQUIRE(seams_insitu_append_dump(&single, kept.c_str()) == SEAMS_INSITU_OK);
  const auto keptText = kept.string();
  const auto mainRecord = (dir / "main.txt").string();
  char *mainArgv[] = {const_cast<char *>("seams-insitu"), const_cast<char *>("replay"),
                      const_cast<char *>(keptText.c_str()),
                      const_cast<char *>(mainRecord.c_str()), const_cast<char *>("3.2")};
  REQUIRE(seams_insitu_main(5, mainArgv) == 0);
  REQUIRE(fs::file_size(dir / "main.txt") > 0);
  REQUIRE(seams_insitu_main(1, mainArgv) == 2);
  mainArgv[1] = const_cast<char *>("nope");
  REQUIRE(seams_insitu_main(4, mainArgv) == 2);
  mainArgv[1] = const_cast<char *>("replay");
  mainArgv[4] = const_cast<char *>("abc");
  REQUIRE(seams_insitu_main(5, mainArgv) == 2);
  mainArgv[4] = const_cast<char *>("0");
  REQUIRE(seams_insitu_main(5, mainArgv) == 2);

  double outXyz[6];
  int outTypes[2];
  int outIds[2];
  const int mask[2] = {1, 0};
  REQUIRE(seams_insitu_pack_owned(-1, nullptr, 1, xyz, types, ids, outXyz, outTypes, outIds) ==
          -1);
  REQUIRE(seams_insitu_pack_owned(0, nullptr, 1, xyz, types, ids, outXyz, outTypes, outIds) == 0);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 1, nullptr, types, ids, outXyz, outTypes, outIds) ==
          -1);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 1, xyz, types, ids, nullptr, outTypes, outIds) == -1);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 1, xyz, types, ids, outXyz, nullptr, outIds) == -1);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 1, xyz, types, ids, outXyz, outTypes, nullptr) == -1);
  REQUIRE(seams_insitu_pack_owned(2, mask, 1, xyz, nullptr, nullptr, outXyz, outTypes, outIds) ==
          1);
  REQUIRE(outTypes[0] == 1);
  REQUIRE(outIds[0] == 1);
  REQUIRE(seams_insitu_pack_owned(2, nullptr, 1, xyz, types, ids, outXyz, outTypes, outIds) == 2);

  int displs[2];
  int total = 0;
  const int rankCounts[2] = {1, 2};
  REQUIRE(seams_insitu_gather_displs(-1, rankCounts, displs, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(1, rankCounts, displs, nullptr) == 1);
  REQUIRE(seams_insitu_gather_displs(1, nullptr, displs, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(1, rankCounts, nullptr, &total) == 1);
  const int negative[1] = {-1};
  REQUIRE(seams_insitu_gather_displs(1, negative, displs, &total) == 1);
  REQUIRE(seams_insitu_gather_displs(0, nullptr, nullptr, &total) == 0);
  REQUIRE(total == 0);
  REQUIRE(seams_insitu_gather_displs(2, rankCounts, displs, &total) == 0);
  REQUIRE(total == 3);
  REQUIRE(displs[1] == 1);

  double scaled[3];
  seams_insitu_scale_nm(nullptr, 1, scaled);
  seams_insitu_scale_nm(xyz, 1, nullptr);
  seams_insitu_scale_nm(xyz, 0, scaled);
  seams_insitu_scale_nm(xyz, 1, scaled);
  REQUIRE(scaled[0] == 0);
  REQUIRE(scaled[1] == 0);

  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  const double rows[9] = {2, 0.2, 0, 0.3, 4, 0.2, 0.1, 0.2, 5};
  REQUIRE(seams_insitu_box_from_rows(nullptr, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, nullptr, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, nullptr, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, nullptr, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, nullptr, &ly, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, nullptr, &lz, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, nullptr, &xy, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, nullptr, &xz, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, nullptr, &yz) == 1);
  REQUIRE(seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, nullptr) == 1);
  const double negativeAxis[9] = {-1, 0, 0, 0, 2, 0, 0, 0, 2};
  REQUIRE(seams_insitu_box_from_rows(negativeAxis, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz,
                                     &yz) == 1);
  const double byLarger[9] = {1, 0, 0, 0, 3, 0, 0, 0, 2};
  REQUIRE(seams_insitu_box_from_rows(byLarger, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) ==
          0);
  REQUIRE(lx == 1);
  REQUIRE(ly == 3);
  const double ayBig[9] = {2, 1, 0, 0, 2, 0, 0, 0, 2};
  REQUIRE(seams_insitu_box_from_rows(ayBig, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  const double azBig[9] = {2, 0, 1, 0, 2, 0, 0, 0, 2};
  REQUIRE(seams_insitu_box_from_rows(azBig, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);
  const double bzBig[9] = {2, 0, 0, 0, 2, 1, 0, 0, 2};
  REQUIRE(seams_insitu_box_from_rows(bzBig, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) == 1);

  int nHc = 0, nDdc = 0, nBoth = 0, nWater = 0;
  int codes[3] = {};
  const int tallyTypes[3] = {1, 2, 1};
  const int hc[3] = {1, 1, 1};
  const int ddc[3] = {1, 0, 0};
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, codes, nullptr, &nDdc, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, codes, &nHc, nullptr, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, codes, &nHc, &nDdc, nullptr, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, codes, &nHc, &nDdc, &nBoth, nullptr);
  seams_insitu_tally_cages(0, tallyTypes, 1, hc, ddc, codes, &nHc, &nDdc, &nBoth, &nWater);
  REQUIRE(nHc == 0);
  seams_insitu_tally_cages(3, nullptr, 1, hc, ddc, codes, &nHc, &nDdc, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, nullptr, ddc, codes, &nHc, &nDdc, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, nullptr, codes, &nHc, &nDdc, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, nullptr, &nHc, &nDdc, &nBoth, &nWater);
  seams_insitu_tally_cages(3, tallyTypes, 1, hc, ddc, codes, &nHc, &nDdc, &nBoth, &nWater);
  REQUIRE(nBoth == 1);
  REQUIRE(nHc == 1);
  REQUIRE(nWater == 0);
  REQUIRE(codes[0] == 3);
  REQUIRE(codes[1] == 0);
  REQUIRE(codes[2] == 1);

  REQUIRE(seams_insitu_imd_session(nullptr, 1, "d", "r", byLarger, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("127.0.0.1", 1, nullptr, "r", byLarger, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("127.0.0.1", 1, "d", nullptr, byLarger, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("127.0.0.1", 1, "d", "r", nullptr, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("127.0.0.1", 0, "d", "r", byLarger, 1) == SEAMS_INSITU_ERR);
  REQUIRE(seams_insitu_imd_session("not-an-ip", 9, "d", "r", byLarger, 1) == SEAMS_INSITU_ERR);
  REQUIRE(imdSession(ImdScript::handshakeBad, 1, dir) == SEAMS_INSITU_ERR);
  REQUIRE(imdSession(ImdScript::close, 0, dir) == SEAMS_INSITU_ERR_IO);
  REQUIRE(imdSession(ImdScript::energyCoord, 1, dir) == SEAMS_INSITU_OK);
  REQUIRE(imdSession(ImdScript::zeroCoord, 1, dir) == SEAMS_INSITU_OK);
  REQUIRE(imdSession(ImdScript::shortCoord, 1, dir) == SEAMS_INSITU_ERR_IO);
  REQUIRE(imdSession(ImdScript::negative, 1, dir) == SEAMS_INSITU_ERR);
  REQUIRE(imdSession(ImdScript::badKind, 1, dir) == SEAMS_INSITU_ERR);

  char *imdArgv[15];
  std::string host = "127.0.0.1";
  std::string portText = "1";
  std::string dumpArg = (dir / "imd.dump").string();
  std::string recArg = (dir / "imd.rec").string();
  std::string nums[10] = {"0", "0", "0", "10", "10", "10", "0", "0", "0", "1"};
  imdArgv[0] = const_cast<char *>("seams-gromacs-imd");
  imdArgv[1] = host.data();
  imdArgv[2] = portText.data();
  imdArgv[3] = dumpArg.data();
  imdArgv[4] = recArg.data();
  for (int i = 0; i < 10; ++i) {
    imdArgv[5 + i] = nums[static_cast<std::size_t>(i)].data();
  }
  REQUIRE(seams_insitu_imd_main(2, imdArgv) == 2);
  portText = "nope";
  REQUIRE(seams_insitu_imd_main(15, imdArgv) == 2);
  portText = "1";
  nums[0] = "x";
  REQUIRE(seams_insitu_imd_main(15, imdArgv) == 2);
  nums[0] = "0";
  nums[9] = "y";
  REQUIRE(seams_insitu_imd_main(15, imdArgv) == 2);

  fs::remove_all(dir);
}

TEST_CASE("language bindings nulls and graph names", "[cover_io]") {
  Lang &lang = langApi();
  REQUIRE(lang.version() != nullptr);
  REQUIRE(std::string(lang.version()).size() > 0);
  REQUIRE(lang.nframes(nullptr) == -1);
  const auto dir = coverTemp("dseams-cover-lang");
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "tiny.dump";
  writeText(dump, dumpText(2, "1 1 0 0 0\n2 1 1.5 0 0\n"));
  const auto twin = dir / "tiny.lammpstrj";
  writeText(twin, dumpText(2, "1 1 0 0 0\n2 1 1.5 0 0\n"));
  REQUIRE(lang.nframes(dump.c_str()) >= 1);
  REQUIRE(lang.nframes(twin.c_str()) >= 1);
  const auto missingXyz = dir / "missing.xyz";
  const int missingFrames = lang.nframes(missingXyz.c_str());
  REQUIRE(missingFrames < 0);
  REQUIRE(lang.read(nullptr, 1, -1) == nullptr);
  void *read = lang.read(dump.c_str(), 1, -1);
  REQUIRE(read != nullptr);
  REQUIRE(lang.nop(read) == 2);
  REQUIRE(lang.positions(read) != nullptr);
  REQUIRE(lang.types(read) != nullptr);
  REQUIRE(lang.ids(read) != nullptr);
  double box[9];
  lang.box(read, box);
  REQUIRE(box[3] > 0);
  lang.free_frame(read);
  REQUIRE(lang.nop(nullptr) == 0);
  REQUIRE(lang.positions(nullptr) == nullptr);
  REQUIRE(lang.types(nullptr) == nullptr);
  REQUIRE(lang.ids(nullptr) == nullptr);
  lang.box(nullptr, box);

  const double xyz[6] = {0, 0, 0, 1.4, 0.2, 0};
  const double ortho[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
  const double shear[9] = {0, 0, 0, 10, 10, 10, 0, 0, 1.5};
  const int types[2] = {1, 2};
  REQUIRE(lang.from_arrays(nullptr, 2, ortho, types, nullptr) == nullptr);
  REQUIRE(lang.from_arrays(xyz, 2, nullptr, types, nullptr) == nullptr);
  REQUIRE(lang.from_arrays(xyz, 0, ortho, types, nullptr) == nullptr);
  void *bare = lang.from_arrays(xyz, 2, ortho, nullptr, nullptr);
  REQUIRE(bare != nullptr);
  REQUIRE(lang.types(bare)[0] == 1);
  REQUIRE(lang.ids(bare)[1] == 2);
  lang.free_frame(bare);
  void *sheared = lang.from_arrays(xyz, 2, shear, types, types);
  REQUIRE(sheared != nullptr);
  lang.box(sheared, box);
  REQUIRE(box[8] == 1.5);

  int labels[2];
  int counts9[9];
  REQUIRE(lang.chill(nullptr, 3.5, 1, "cutoff", labels, counts9) == 1);
  REQUIRE(lang.chill(sheared, 3.5, 0, nullptr, nullptr, nullptr) == 0);
  REQUIRE(lang.chill_plus(sheared, 3.5, -1, "", labels, counts9) == 0);
  REQUIRE(lang.chill(sheared, 3.5, 1, "mutual", labels, counts9) == 0);
  REQUIRE(lang.chill(sheared, 3.5, 1, "nope", labels, counts9) == 1);
  REQUIRE(std::string(lang.last_error()).find("unknown") != std::string::npos);

  int ringCounts[13] = {};
  REQUIRE(lang.rings(nullptr, 3.5, 1, "cutoff", 6, ringCounts, 13) == 1);
  REQUIRE(lang.rings(sheared, 3.5, 1, "cutoff", 2, ringCounts, 13) == 1);
  REQUIRE(lang.rings(sheared, 3.5, 1, "seeded", 6, ringCounts, 7) == 1);
  REQUIRE(lang.rings(sheared, 3.5, 1, "union", 6, ringCounts, 13) == 0);
  REQUIRE(lang.rings(sheared, 3.5, 1, "knn-union", 6, ringCounts, 13) == 0);
  REQUIRE(lang.rings(sheared, 3.5, 1, "knn", 6, ringCounts, 13) == 0);
  REQUIRE(lang.rings(sheared, 3.5, 1, "knn-mutual", 6, ringCounts, 13) == 0);
  REQUIRE(lang.rings(sheared, 3.5, 1, "mutual", 6, ringCounts, 13) == 0);

  double ql[2];
  double qlbar[2];
  REQUIRE(lang.steinhardt(nullptr, 3.5, 1, 6, ql, qlbar) == 1);
  REQUIRE(lang.steinhardt(sheared, 3.5, 1, 6, nullptr, nullptr) == 0);
  REQUIRE(lang.steinhardt(sheared, 3.5, 1, 4, ql, qlbar) == 0);

  int perAtom[2];
  int nHex = -1, nCubic = -1, nWater = -1, nBoth = -1;
  REQUIRE(lang.cages(nullptr, 3.5, 1, 4, "cutoff", 0, perAtom, &nHex, &nCubic, &nWater, &nBoth) ==
          1);
  REQUIRE(lang.cages(sheared, 3.5, 1, 0, "cutoff", 1, perAtom, &nHex, &nCubic, &nWater, &nBoth) ==
          0);
  REQUIRE(nHex >= 0);
  REQUIRE(lang.cages(sheared, 3.5, 1, 4, "seeded", 0, nullptr, nullptr, nullptr, nullptr,
                     nullptr) == 0);
  int nCages = 0;
  int nAtoms = 0;
  REQUIRE(lang.signature(nullptr, 3.5, 1, "cutoff", "512", &nCages, &nAtoms) == 1);
  REQUIRE(lang.signature(sheared, 3.5, 1, "cutoff", "512", &nCages, &nAtoms) == 0);
  REQUIRE(lang.signature_edge(sheared, 3.5, 0, "mutual", "hc", "none", &nCages, &nAtoms) == 0);
  REQUIRE(lang.signature(sheared, 3.5, 1, "nope", "512", &nCages, &nAtoms) == 1);

  REQUIRE(lang.rdf_nbins(5.0, 10) > 0);
  double radius[4];
  double gee[4];
  int rdfCount[4];
  double volume = 0;
  REQUIRE(lang.rdf(nullptr, 5, 4, 1, 1, radius, gee, rdfCount, &volume) == 1);
  REQUIRE(lang.rdf(sheared, 5, 4, 1, 2, radius, gee, rdfCount, &volume) == 0);
  REQUIRE(volume > 0);
  const int bins = lang.rdf2d_nbin(6.0, 0.5);
  REQUIRE(bins > 0);
  std::vector<double> r2(static_cast<std::size_t>(bins));
  std::vector<double> g2(static_cast<std::size_t>(bins));
  REQUIRE(lang.rdf2d(nullptr, 6, 0.5, nullptr, r2.data(), g2.data()) == 1);
  const auto rdfDir = (dir / "rdf").string();
  REQUIRE(lang.rdf2d(sheared, 6, 0.5, rdfDir.c_str(), r2.data(), g2.data()) == 0);
  REQUIRE(lang.rdf2d(sheared, 6, 0.5, (rdfDir + "/").c_str(), r2.data(), g2.data()) == 0);

  int nRows = 0;
  REQUIRE(lang.prisms(sheared, 3.5, 1, 6, nullptr, nullptr, &nRows) == 0);
  REQUIRE(nRows == 4);
  std::vector<int> prismBuf(static_cast<std::size_t>(nRows) * 3);
  const auto prismDir = dir / "prisms";
  fs::create_directories(prismDir);
  REQUIRE(lang.prisms(sheared, 3.5, 1, 6, (prismDir.string() + "/").c_str(), prismBuf.data(),
                      &nRows) == 0);
  REQUIRE(lang.prisms(nullptr, 3.5, 1, 6, nullptr, nullptr, &nRows) == 1);
  int mono[13] = {};
  const auto sheet = (dir / "sheet").string();
  REQUIRE(lang.monolayer(sheared, 3.5, 1, 6, 100.0, sheet.c_str(), mono, 13) == 0);
  REQUIRE(lang.monolayer(sheared, 3.5, 1, 2, 100.0, nullptr, mono, 13) == 1);
  REQUIRE(lang.monolayer(sheared, 3.5, 1, 6, 100.0, (sheet + "/").c_str(), mono, 13) == 0);

  const auto xyzPath = dir / "plain.xyz";
  writeText(xyzPath, "1\n\nXx 0 0 0\n");
  void *loose = lang.read(xyzPath.c_str(), 1, -1);
  REQUIRE(loose != nullptr);
  lang.box(loose, box);
  lang.free_frame(loose);
  lang.free_frame(sheared);
  fs::remove_all(dir);
}

TEST_CASE("cli edges for ranges, axes, and unclassified ice", "[cover_io]") {
  const auto dir = coverTemp("dseams-cover-cli");
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto one = dir / "one.lammpstrj";
  writeText(one, dumpText(1, "1 1 1 1 1\n"));
  auto chill = runSeams({"chill", one.string(), "--frame", "1", "--type", "1", "--cutoff", "3.5"});
  REQUIRE(chill.exited);
  REQUIRE(chill.status == 0);
  REQUIRE(chill.out.find("nop") != std::string::npos);

  const auto empty = dir / "empty.lammpstrj";
  writeText(empty, dumpText(0, ""));
  auto rdf = runSeams({"rdf", empty.string(), "--type", "0", "--cutoff", "3", "--bins", "4"});
  REQUIRE(rdf.exited);
  REQUIRE(rdf.status == 0);
  REQUIRE(rdf.out.find("# types") != std::string::npos);

  const auto oxy = dir / "oxy.lammpstrj";
  writeText(oxy, dumpText(2, "1 1 0 0 0\n2 1 1 0 0\n"));
  auto f4 = runSeams({"f4", oxy.string(), "--type", "0", "--htype", "1", "--cutoff", "3.5"});
  REQUIRE(f4.exited);
  REQUIRE(f4.status == 0);
  REQUIRE(f4.out.find("f4") != std::string::npos);

  const auto xyz = dir / "late.xyz";
  writeText(xyz, "1\n\nO 0 0 0\n");
  auto late = runSeams({"read", xyz.string(), "--frame", "2", "--last", "3"});
  REQUIRE(late.exited);
  REQUIRE(late.status == 2);
  REQUIRE(late.err.find("past the last frame") != std::string::npos);
  // A start after the requested end is an empty window, not a missing frame.
  auto inverted = runSeams({"read", xyz.string(), "--frame", "3", "--last", "1"});
  REQUIRE(inverted.exited);
  REQUIRE(inverted.status == 0);
  REQUIRE(inverted.out.empty());
  REQUIRE(inverted.err.empty());

  const auto nodot = dir / "nodot";
  writeText(nodot, dumpText(1, "1 1 0 0 0\n"));
  auto plain = runSeams({"read", nodot.string(), "--frame", "1", "--last", "2"});
  REQUIRE(plain.exited);

  const auto pdb = dir / "one.pdb";
  writeText(pdb, "END\n");
  auto pdbRun = runSeams({"read", pdb.string(), "--frame", "1", "--last", "2"});
  REQUIRE(pdbRun.exited);
  const auto gro = dir / "one.gro";
  writeText(gro, "title\n0\n   1.0 1.0 1.0\n");
  auto groRun = runSeams({"read", gro.string(), "--frame", "1", "--last", "2"});
  REQUIRE(groRun.exited);
  const auto dcd = dir / "one.dcd";
  writeText(dcd, "nope");
  auto dcdRun = runSeams({"read", dcd.string(), "--frame", "1", "--last", "2"});
  REQUIRE(dcdRun.exited);

  const auto junk = dir / "junk.lammpstrj";
  writeText(junk, "this is not a dump\n");
  auto junkRun = runSeams({"rings", junk.string(), "--frame", "1", "--last", "2"});
  REQUIRE(junkRun.exited);

  const auto three = dir / "three.lammpstrj";
  writeText(three, dumpText(3, "1 1 0 0 0\n2 2 1 0 0\n3 3 2 0 0\n"));
  auto hint = runSeams({"density-z", three.string(), "--type", "0", "--bins", "2"});
  REQUIRE(hint.exited);
  REQUIRE(hint.status == 0);
  REQUIRE(hint.err.find("more than two") != std::string::npos);

  for (const char *axis : {"X", "Y", "Z", "1", "2"}) {
    auto density = runSeams({"density-z", three.string(), "--axis", axis, "--bins", "2"});
    REQUIRE(density.exited);
    REQUIRE(density.status == 0);
  }
  auto waters = runSeams(
      {"rings", three.string(), "--water-types", "1,,2", "--frame", "1", "--last", "2"});
  REQUIRE(waters.exited);

  const auto ions = dir / "ions.lammpstrj";
  writeText(ions, dumpText(2, "1 1 0 0 0\n2 2 1 0 0\n"));
  auto domains = runSeams({"domains", ions.string(), "--site", "1=cationHead,2=anion", "--cutoff",
                           "4"});
  REQUIRE(domains.exited);
  REQUIRE(domains.status == 0);

  auto sig = runSeams({"cages", one.string(), "--signature", "", "--cutoff", "3.5"});
  REQUIRE(sig.exited);

  const auto pair = dir / "pair.lammpstrj";
  writeText(pair, "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
                  "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                  "ITEM: ATOMS id mol type x y z\n"
                  "1 1 1 0 0 0\n2 1 2 0.1 0 0\n3 2 1 2.8 0 0\n");
  auto bonds = runSeams({"hbonds", pair.string(), "--type", "1", "--htype", "2", "--cutoff", "4",
                         "--hdist", "3", "--hangle", "50", "--donors"});
  REQUIRE(bonds.exited);
  REQUIRE(bonds.status == 0);

  const auto blanks = dir / "blanks.lammpstrj";
  writeText(blanks, dumpText(0, "", 0) + dumpText(0, "", 1));
  auto strict = runSeams(
      {"read", blanks.string(), "--frame", "1", "--last", "2", "--strict-input", "--jobs", "2"});
  REQUIRE(strict.exited);
  REQUIRE(strict.status == 2);

  const auto root = repoRoot();
  const auto ice = root / "input" / "traj" / "genice_sI.lammpstrj";
  if (fs::exists(ice)) {
    auto census = runSeams({"ions", ice.string(), "--frame", "1", "--type", "1", "--ion-types",
                            "1", "--ion-cutoff", "4", "--cutoff", "3.5"});
    REQUIRE(census.exited);
    REQUIRE(census.out.find("ions-in-ice") != std::string::npos);
  }

  fs::remove_all(dir);
}
