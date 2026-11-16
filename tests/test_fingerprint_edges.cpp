#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <topo_fingerprint.hpp>

namespace primitive {
std::vector<int> ringCensusCounts(const std::vector<std::vector<int>> &nList, int maxDepth);
}

#include <algorithm>
#include <array>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

// (max degree + 63) / 64. This is the shell width the ring search picks
// for the whole frame.
int shellWords(const topo::Rows &rows) {
  int maxDeg = 0;
  const int n = static_cast<int>(rows.size());
  for (int i = 0; i < n; i++) {
    int deg = 0;
    for (std::size_t k = 1; k < rows[static_cast<std::size_t>(i)].size(); k++) {
      const int j = rows[static_cast<std::size_t>(i)][k];
      if (j >= 0 && j < n && j != i) {
        deg++;
      }
    }
    maxDeg = std::max(maxDeg, deg);
  }
  return (maxDeg + 63) >> 6;
}

topo::Rows withSelf(topo::Rows rows) {
  for (int i = 0; i < static_cast<int>(rows.size()); i++) {
    auto &row = rows[static_cast<std::size_t>(i)];
    row.insert(row.begin(), i);
  }
  return rows;
}

// A hub of `leaves` neighbours, plus one chord so a vertex has two shortest
// paths. The chord is a triangle; a disjoint 4-cycle adds rings that do not
// touch the hub.
topo::Rows hubFrame(int leaves) {
  const int hubN = leaves + 1;
  const int n = hubN + 4;
  topo::Rows rows(static_cast<std::size_t>(n));
  for (int j = 1; j <= leaves; j++) {
    rows[0].push_back(j);
    rows[static_cast<std::size_t>(j)].push_back(0);
  }
  rows[1].push_back(2);
  rows[2].push_back(1);
  const int c = hubN;
  const int cycle[4] = {c, c + 1, c + 2, c + 3};
  for (int k = 0; k < 4; k++) {
    const int a = cycle[k];
    const int b = cycle[(k + 1) % 4];
    rows[static_cast<std::size_t>(a)].push_back(b);
    rows[static_cast<std::size_t>(b)].push_back(a);
  }
  return withSelf(std::move(rows));
}

std::vector<int> tally(const std::vector<std::vector<int>> &rings, int maxDepth) {
  std::vector<int> hist(static_cast<std::size_t>(maxDepth) + 1, 0);
  for (const auto &ring : rings) {
    if (static_cast<int>(ring.size()) <= maxDepth) {
      hist[ring.size()]++;
    }
  }
  return hist;
}

void expectCensus(const topo::Rows &rows, int words) {
  REQUIRE(shellWords(rows) == words);
  const int maxDepth = 7;
  const auto fast = primitive::ringCensusCounts(rows, maxDepth);
  const auto listed = primitive::ringNetwork(rows, maxDepth);
  const auto legacy = primitive::ringNetworkBalls(rows, maxDepth);
  REQUIRE(fast == tally(listed, maxDepth));
  REQUIRE(fast == tally(legacy, maxDepth));
  const auto fp = topo::fingerprint(rows, 2, maxDepth);
  REQUIRE(fp.ringCensus == fast);
  REQUIRE(fp.atomKeys.size() == rows.size());
  REQUIRE(fp.wlAtom.size() == rows.size());
}

void expectSameFingerprint(const topo::FrameFingerprint &a, const topo::FrameFingerprint &b) {
  REQUIRE(a.method == b.method);
  REQUIRE(a.key == b.key);
  REQUIRE(a.atomKeys == b.atomKeys);
  REQUIRE(a.classes == b.classes);
  REQUIRE(a.ringCensus == b.ringCensus);
  REQUIRE(a.wlAtom == b.wlAtom);
  REQUIRE(a.hops == b.hops);
  REQUIRE(a.colored == b.colored);
}

#ifdef SEAMS_HAS_OPENMP
struct ThreadWidth {
  int saved;
  explicit ThreadWidth(int n) : saved(omp_get_max_threads()) { omp_set_num_threads(n); }
  ~ThreadWidth() { omp_set_num_threads(saved); }
};
#endif

void expectThreadsAgree(const topo::Rows &rows) {
#ifdef SEAMS_HAS_OPENMP
  topo::FrameFingerprint one;
  topo::FrameFingerprint three;
  topo::FrameFingerprint many;
  {
    ThreadWidth hold(1);
    one = topo::fingerprint(rows, 2, 7);
  }
  {
    ThreadWidth hold(3);
    three = topo::fingerprint(rows, 2, 7);
  }
  {
    ThreadWidth hold(std::max(8, static_cast<int>(rows.size()) + 3));
    many = topo::fingerprint(rows, 2, 7);
  }
  expectSameFingerprint(one, three);
  expectSameFingerprint(one, many);
#else
  const auto fp = topo::fingerprint(rows, 2, 7);
  REQUIRE(fp.atomKeys.size() == rows.size());
#endif
}

molSys::PointCloud<molSys::Point<double>, double>
cloudOf(const std::vector<std::array<double, 3>> &pos, std::vector<double> box,
        std::vector<double> low) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = static_cast<int>(pos.size());
  cloud.currentFrame = 1;
  cloud.box = std::move(box);
  cloud.boxLow = std::move(low);
  for (std::size_t i = 0; i < pos.size(); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = static_cast<int>(i) + 1;
    p.molID = p.atomID;
    p.x = pos[i][0];
    p.y = pos[i][1];
    p.z = pos[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = static_cast<int>(i);
  }
  return cloud;
}

// One neighbour when the minimum image lies inside the cutoff. A second
// periodic image is not a second neighbour.
topo::Rows minimumImageRows(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                            double cutoff) {
  topo::Rows rows(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    auto &row = rows[static_cast<std::size_t>(i)];
    row.push_back(i);
    for (int j = 0; j < cloud.nop; j++) {
      if (j != i && gen::periodicDist(cloud, i, j) < cutoff) {
        row.push_back(j);
      }
    }
    std::sort(row.begin() + 1, row.end());
  }
  return rows;
}

topo::Rows cutoffRows(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                      double cutoff) {
  const auto ids = nneigh::neighListO(cutoff, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, ids);
  for (auto &row : rows) {
    if (row.size() > 1) {
      std::sort(row.begin() + 1, row.end());
    }
  }
  return rows;
}

} // namespace

TEST_CASE("an empty frame and a single atom have no rings", "[fingerprint]") {
  const topo::Rows empty;
  const auto none = topo::fingerprint(empty, 2, 7);
  REQUIRE(none.atomKeys.empty());
  REQUIRE(none.wlAtom.empty());
  REQUIRE(none.classes.empty());
  REQUIRE(none.ringCensus == primitive::ringCensusCounts(empty, 7));
  for (std::size_t s = 3; s < none.ringCensus.size(); s++) {
    REQUIRE(none.ringCensus[s] == 0);
  }
  expectSameFingerprint(none, topo::fingerprint(empty, 2, 7));

  // Degree 0 is below a one-word shell. The census is still the legacy list.
  const topo::Rows lone = {{0}};
  REQUIRE(shellWords(lone) == 0);
  const auto fast = primitive::ringCensusCounts(lone, 7);
  REQUIRE(fast == tally(primitive::ringNetwork(lone, 7), 7));
  REQUIRE(fast == tally(primitive::ringNetworkBalls(lone, 7), 7));
  const auto one = topo::fingerprint(lone, 2, 7);
  REQUIRE(one.ringCensus == fast);
  for (std::size_t s = 3; s < one.ringCensus.size(); s++) {
    REQUIRE(one.ringCensus[s] == 0);
  }
}

TEST_CASE("a path has no rings on the one-word shell", "[fingerprint]") {
  topo::Rows path(4);
  path[0] = {1};
  path[1] = {0, 2};
  path[2] = {1, 3};
  path[3] = {2};
  path = withSelf(std::move(path));
  REQUIRE(shellWords(path) == 1);
  expectCensus(path, 1);
  for (std::size_t s = 3; s < 8; s++) {
    REQUIRE(primitive::ringCensusCounts(path, 7)[s] == 0);
  }
}

TEST_CASE("one-word, two-word and three-word shells match the legacy search",
          "[fingerprint]") {
  // Degree 20 still fits in one word and is past the path cap, so the
  // chordless one-word loop runs. Degree 70 is two words. Degree 130 is
  // three words.
  expectCensus(hubFrame(20), 1);
  expectCensus(hubFrame(70), 2);
  expectCensus(hubFrame(130), 3);
}

TEST_CASE("a one-way triangle matches the legacy ring list", "[fingerprint]") {
  topo::Rows rows(3);
  rows[0] = {1};
  rows[1] = {2};
  rows[2] = {0};
  rows = withSelf(std::move(rows));
  REQUIRE(shellWords(rows) == 1);
  expectCensus(rows, 1);
}

TEST_CASE("odd thread counts and more threads than atoms agree", "[fingerprint]") {
  expectThreadsAgree(hubFrame(20));
  expectThreadsAgree(hubFrame(70));
  const topo::Rows triangle = withSelf(topo::Rows{{1, 2}, {0, 2}, {0, 1}});
  expectThreadsAgree(triangle);
}

TEST_CASE("atoms on the cell faces use the minimum image", "[fingerprint]") {
  // 0 and 9.8 are 0.2 apart through the face. The corner sees three faces.
  const double box = 10.0;
  const double cutoff = 0.5;
  const auto cloud = cloudOf({{0.0, 0.0, 0.0},
                              {9.8, 0.0, 0.0},
                              {0.0, 9.8, 0.0},
                              {0.0, 0.0, 9.8},
                              {9.8, 9.8, 9.8}},
                             {box, box, box}, {0.0, 0.0, 0.0});
  const auto fast = cutoffRows(cloud, cutoff);
  const auto ref = minimumImageRows(cloud, cutoff);
  REQUIRE(fast == ref);
  expectSameFingerprint(topo::fingerprint(fast, 2, 7), topo::fingerprint(ref, 2, 7));
  expectCensus(fast, shellWords(fast));
}

TEST_CASE("a non-cubic box and a sheared box match the minimum image", "[fingerprint]") {
  const auto ortho = cloudOf({{0.0, 0.0, 0.0}, {11.7, 0.2, 0.1}, {0.2, 4.6, 0.1}},
                             {12.0, 5.0, 8.0}, {0.0, 0.0, 0.0});
  const auto orthoFast = cutoffRows(ortho, 0.6);
  REQUIRE(orthoFast == minimumImageRows(ortho, 0.6));
  expectSameFingerprint(topo::fingerprint(orthoFast, 2, 7),
                        topo::fingerprint(minimumImageRows(ortho, 0.6), 2, 7));

  const double xy = 1.5;
  const double xz = 0.4;
  const double yz = 0.25;
  const auto shear = cloudOf({{0.2, 0.2, 0.2}, {9.6, 0.3, 0.2}, {1.0, 4.7, 0.4}},
                             {10.0, 5.0, 6.0, xy, xz, yz}, {0.0, 0.0, 0.0});
  const auto shearFast = cutoffRows(shear, 0.8);
  REQUIRE(shearFast == minimumImageRows(shear, 0.8));
  expectSameFingerprint(topo::fingerprint(shearFast, 2, 7),
                        topo::fingerprint(minimumImageRows(shear, 0.8), 2, 7));
}

TEST_CASE("a cutoff larger than half the box lists each neighbour once", "[fingerprint]") {
  // Both the short gap and the long way around the box are inside the
  // cutoff. The list keeps the minimum image, once.
  const auto cloud =
      cloudOf({{0.0, 0.0, 0.0}, {1.5, 0.0, 0.0}, {0.2, 0.2, 0.2}}, {4.0, 4.0, 4.0},
              {0.0, 0.0, 0.0});
  const double cutoff = 3.0;
  REQUIRE(cutoff > 2.0);
  const auto fast = cutoffRows(cloud, cutoff);
  const auto ref = minimumImageRows(cloud, cutoff);
  REQUIRE(fast == ref);
  for (const auto &row : fast) {
    auto copy = row;
    std::sort(copy.begin(), copy.end());
    REQUIRE(std::unique(copy.begin(), copy.end()) == copy.end());
  }
  expectSameFingerprint(topo::fingerprint(fast, 2, 7), topo::fingerprint(ref, 2, 7));
#ifdef SEAMS_HAS_OPENMP
  topo::Rows odd;
  topo::Rows crowd;
  {
    ThreadWidth hold(3);
    odd = cutoffRows(cloud, cutoff);
  }
  {
    ThreadWidth hold(8);
    crowd = cutoffRows(cloud, cutoff);
  }
  REQUIRE(odd == ref);
  REQUIRE(crowd == ref);
#endif
}
