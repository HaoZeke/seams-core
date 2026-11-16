#include <catch2/catch_test_macros.hpp>

#include "fresh_temp.hpp"

#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_fingerprint.hpp>

#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace topo {
FrameFingerprint fingerprintPaths(const Rows &rows, int hops, int maxRingSize,
                                  const std::vector<int> &colors, int which);
}

namespace {

topo::Rows cycle(int n) {
  topo::Rows rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)] = {i, (i + 1) % n, (i + n - 1) % n};
  }
  return rows;
}

topo::Rows path(int n) {
  topo::Rows rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    if (i > 0) {
      rows[static_cast<std::size_t>(i)].push_back(i - 1);
    }
    if (i + 1 < n) {
      rows[static_cast<std::size_t>(i)].push_back(i + 1);
    }
  }
  return rows;
}

void expectEqual(const topo::FrameFingerprint &a, const topo::FrameFingerprint &b) {
  REQUIRE(a.method == b.method);
  REQUIRE(a.key == b.key);
  REQUIRE(a.atomKeys == b.atomKeys);
  REQUIRE(a.classes == b.classes);
  REQUIRE(a.ringCensus == b.ringCensus);
  REQUIRE(a.wlAtom == b.wlAtom);
  REQUIRE(a.hops == b.hops);
  REQUIRE(a.colored == b.colored);
}

// The ball the local key hashes: breadth-first neighbourhood, self dropped,
// neighbours outside the ball dropped, a short row sorted.
topo::Rows inducedBall(const topo::Rows &rows, int atom, int hops, std::vector<int> &localColors,
                       const std::vector<int> &colors) {
  const auto ball = topo::hopNeighbourhood(rows, atom, hops);
  std::map<int, int> loc;
  for (int i = 0; i < static_cast<int>(ball.size()); i++) {
    loc[ball[static_cast<std::size_t>(i)]] = i;
  }
  topo::Rows adj(ball.size());
  localColors.clear();
  if (colors.size() == rows.size()) {
    localColors.resize(ball.size());
    for (int i = 0; i < static_cast<int>(ball.size()); i++) {
      localColors[static_cast<std::size_t>(i)] = colors[static_cast<std::size_t>(ball[static_cast<std::size_t>(i)])];
    }
  }
  for (int v : ball) {
    const int i = loc[v];
    const auto &row = rows[static_cast<std::size_t>(v)];
    auto &out = adj[static_cast<std::size_t>(i)];
    for (std::size_t k = 1; k < row.size(); k++) {
      const auto it = loc.find(row[k]);
      if (it == loc.end() || it->second == i) {
        continue;
      }
      out.push_back(it->second);
    }
    if (out.size() <= 16) {
      std::sort(out.begin(), out.end());
    }
  }
  return adj;
}

void expectKeyMatchesHash(const topo::Rows &rows, int atom, int hops,
                          const std::vector<int> &colors) {
  std::vector<int> local;
  const auto adj = inducedBall(rows, atom, hops, local, colors);
  const auto expect = topo::hex(topo::wlHash(adj, 0, hops + 2, local));
  const auto first = topo::localKey(rows, atom, hops, colors);
  const auto second = topo::localKey(rows, atom, hops, colors);
  REQUIRE(first.method == "wl");
  REQUIRE(first.key == expect);
  REQUIRE(second.key == expect);
  REQUIRE(second.edges == first.edges);
  REQUIRE(second.vertices == static_cast<int>(adj.size()));
}

std::string frameText(const topo::FrameFingerprint &fp) {
  std::ostringstream out;
  out << "nop " << fp.atomKeys.size() << " method " << fp.method << " key " << fp.key
      << " classes " << fp.classes.size() << " rings";
  for (std::size_t s = 3; s < fp.ringCensus.size(); s++) {
    out << ' ' << s << ':' << fp.ringCensus[s];
  }
  out << '\n';
  for (const auto &key : fp.atomKeys) {
    out << key << '\n';
  }
  return out.str();
}

void writeText(const std::filesystem::path &path, const std::string &text) {
  std::ofstream out(path);
  REQUIRE(out);
  out << text;
}

std::string readText(const std::filesystem::path &path) {
  std::ifstream in(path);
  REQUIRE(in);
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}

#ifdef SEAMS_HAS_OPENMP
struct ThreadWidth {
  int saved;
  explicit ThreadWidth(int n) : saved(omp_get_max_threads()) { omp_set_num_threads(n); }
  ~ThreadWidth() { omp_set_num_threads(saved); }
};
#endif

} // namespace

TEST_CASE("an empty or odd colour vector is the uncolored key", "[fingerprint]") {
  const auto rows = cycle(6);
  const auto plain = topo::fingerprint(rows, 2, 7);
  REQUIRE_FALSE(plain.colored);
  expectEqual(plain, topo::fingerprint(rows, 2, 7, {}));
  // One colour, and three colours, neither matching the six rows.
  expectEqual(plain, topo::fingerprint(rows, 2, 7, {1}));
  expectEqual(plain, topo::fingerprint(rows, 2, 7, {1, 2, 3}));
  REQUIRE_FALSE(topo::fingerprint(rows, 2, 7, {1, 2, 3}).colored);
}

TEST_CASE("odd colour types agree on the shared ball and the per-atom ball", "[fingerprint]") {
  constexpr int n = 6;
  topo::Rows rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    for (int j = 0; j < n; j++) {
      if (j != i) {
        rows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
  }
  // Three types, an odd set, one per atom.
  std::vector<int> colors(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    colors[static_cast<std::size_t>(i)] = (i % 3) + 1;
  }
  const auto shared = topo::fingerprintPaths(rows, 1, 7, colors, 1);
  const auto perAtom = topo::fingerprintPaths(rows, 1, 7, colors, 2);
  expectEqual(shared, perAtom);
  REQUIRE(shared.colored);
  REQUIRE(shared.key != topo::fingerprint(rows, 1, 7).key);

  const auto pathRows = path(5);
  const auto automatic = topo::fingerprint(pathRows, 2, 7, {1, 2, 3, 1, 2});
  const auto reference = topo::fingerprintPaths(pathRows, 2, 7, {1, 2, 3, 1, 2}, 2);
  expectEqual(automatic, reference);
  REQUIRE(automatic.colored);
}

TEST_CASE("a repeated ball hits the memo and still matches the hash", "[fingerprint]") {
  const auto rows = cycle(6);
  const std::vector<int> colors = {1, 2, 3, 1, 2, 3};
  expectKeyMatchesHash(rows, 0, 2, {});
  expectKeyMatchesHash(rows, 3, 2, {});
  expectKeyMatchesHash(rows, 0, 2, colors);
  const auto left = topo::localKey(rows, 0, 2, colors);
  const auto right = topo::localKey(rows, 3, 2, colors);
  REQUIRE(left.key == right.key);
  REQUIRE(left.edges == right.edges);

  const auto other = path(4);
  const auto different = topo::localKey(other, 0, 2, {});
  REQUIRE(different.key != left.key);
  expectKeyMatchesHash(rows, 1, 2, colors);
  expectKeyMatchesHash(other, 1, 2, {});
}

TEST_CASE("one, three and eight threads write the same fingerprint", "[fingerprint]") {
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "two-frames.lammpstrj";
  writeText(dump,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n4\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 0 0 0\n2 2 1 0 0\n3 3 0 1 0\n4 1 1 1 0\n"
            "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n4\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 0.1 0 0\n2 2 1.1 0 0\n3 3 0 1.1 0\n4 1 1 1.1 0\n");

  auto report = [&](int threads) {
    std::vector<std::string> lines(2);
#ifdef SEAMS_HAS_OPENMP
    ThreadWidth hold(threads);
#endif
    sinp::forEachLammpsFrame(
        dump.string(), 1, 2, 0,
        [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
          std::vector<int> colors;
          colors.reserve(static_cast<std::size_t>(cloud.nop));
          topo::Rows rows(static_cast<std::size_t>(cloud.nop));
          for (int i = 0; i < cloud.nop; i++) {
            colors.push_back(cloud.pts[static_cast<std::size_t>(i)].type);
            rows[static_cast<std::size_t>(i)].push_back(i);
            for (int j = 0; j < cloud.nop; j++) {
              if (j != i && gen::periodicDist(cloud, i, j) < 1.5) {
                rows[static_cast<std::size_t>(i)].push_back(j);
              }
            }
          }
          const auto fp = topo::fingerprint(rows, 2, 7, colors);
          lines[static_cast<std::size_t>(frame - 1)] = frameText(fp);
        },
        threads);
    std::ostringstream all;
    for (const auto &line : lines) {
      all << line;
    }
    const auto out = dir / ("threads-" + std::to_string(threads) + ".txt");
    writeText(out, all.str());
    return readText(out);
  };

  const auto one = report(1);
  const auto three = report(3);
  const auto eight = report(8);
  REQUIRE_FALSE(one.empty());
  REQUIRE(three == one);
  REQUIRE(eight == one);
  sinp::dropLammpsDumpIndex(dump.string());
  std::filesystem::remove_all(dir);
}

TEST_CASE("a bad atom line is dropped and a truncated frame errors", "[fingerprint]") {
  const auto dir = dseams_test::freshDir();
  const auto bad = dir / "bad.lammpstrj";
  writeText(bad,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
            "not-an-atom\n");
  // A token that is not a coordinate is not an atom. The following
  // well-formed line is kept, and the read does not throw.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  REQUIRE_NOTHROW(cloud = sinp::readLammpsTrj(bad.string(), 1, cloud));
  REQUIRE(cloud.nop == 0);

  const auto mixed = dir / "mixed.lammpstrj";
  writeText(mixed,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 abc 0 0\n2 1 1 0 0\n");
  REQUIRE_NOTHROW(cloud = sinp::readLammpsTrj(mixed.string(), 1, cloud));
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.pts.size() == 1);
  REQUIRE(cloud.pts[0].atomID == 2);

  const auto empty = dir / "empty.lammpstrj";
  writeText(empty, "");
  molSys::PointCloud<molSys::Point<double>, double> none;
  REQUIRE_NOTHROW(none = sinp::readLammpsTrj(empty.string(), 1, none));
  REQUIRE(none.nop == 0);

  const auto shortDump = dir / "short.lammpstrj";
  writeText(shortDump,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 0 0 0\n");
  try {
    sinp::readLammpsTrj(shortDump.string(), 1, cloud);
    FAIL("expected a truncated dump");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) ==
            "LAMMPS dump frame is truncated: expected 3 atoms, read 1.");
  }
  sinp::dropLammpsDumpIndex(bad.string());
  sinp::dropLammpsDumpIndex(mixed.string());
  sinp::dropLammpsDumpIndex(empty.string());
  sinp::dropLammpsDumpIndex(shortDump.string());
  std::filesystem::remove_all(dir);
}

TEST_CASE("parking the linkcell team does not change the neighbour rows", "[fingerprint]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = 6;
  cloud.currentFrame = 1;
  cloud.box = {12.0, 8.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double pos[6][3] = {{0.2, 0.2, 0.2}, {1.0, 0.2, 0.2}, {0.2, 1.1, 0.2},
                            {11.7, 0.3, 0.2}, {0.4, 7.6, 0.3}, {0.3, 0.3, 9.6}};
  for (int i = 0; i < 6; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.x = pos[i][0];
    p.y = pos[i][1];
    p.z = pos[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  std::vector<int> offPark, colPark, offRun, colRun;
  nneigh::kNearestCsr(cloud, 4, 3.0, 1, true, offPark, colPark, true);
  nneigh::kNearestCsr(cloud, 4, 3.0, 1, true, offRun, colRun, false);
  REQUIRE(offPark == offRun);
  REQUIRE(colPark == colRun);
  const auto listed = nneigh::neighListO(3.0, cloud, 1);
  REQUIRE(listed.size() == static_cast<std::size_t>(cloud.nop));
}
