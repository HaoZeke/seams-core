#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <domain.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <thread_team.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

namespace {

molSys::PointCloud<molSys::Point<double>, double> jitteredIce(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double basis[8][3] = {{0, 0, 0},       {0.5, 0.5, 0},   {0.5, 0, 0.5},
                              {0, 0.5, 0.5},   {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  unsigned s = 1;
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          const double xyz[3] = {i + b[0], j + b[1], k + b[2]};
          double *out[3] = {&p.x, &p.y, &p.z};
          for (int c = 0; c < 3; c++) {
            s = s * 1664525u + 1013904223u;
            const double jitter = (static_cast<int>(s % 1000) / 1000.0 - 0.5) * 0.3;
            *out[c] = std::fmod((xyz[c] * lattice + jitter), side);
            if (*out[c] < 0.0) {
              *out[c] += side;
            }
          }
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

struct Snapshot {
  std::vector<std::vector<int>> neighbours;
  std::vector<int> published;
  std::vector<int> ovito;
  std::vector<double> q6;
  std::vector<double> q6bar;
  std::vector<int> rings;
  std::vector<std::vector<int>> ringList;
};

std::vector<std::vector<int>> rowsByIndex(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const std::vector<std::vector<int>> &byId) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
    if (static_cast<int>(byId.size()) <= i) {
      continue;
    }
    for (std::size_t k = 1; k < byId[static_cast<std::size_t>(i)].size(); k++) {
      const auto it = cloud.idIndexMap.find(byId[static_cast<std::size_t>(i)][k]);
      if (it != cloud.idIndexMap.end()) {
        rows[static_cast<std::size_t>(i)].push_back(it->second);
      }
    }
  }
  return rows;
}

std::vector<int> census(const std::vector<std::vector<int>> &rings) {
  std::vector<int> counts(7, 0);
  for (const auto &ring : rings) {
    if (ring.size() < counts.size()) {
      counts[ring.size()]++;
    }
  }
  return counts;
}

Snapshot runOnce(const molSys::PointCloud<molSys::Point<double>, double> &base, int threads) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#else
  (void)threads;
#endif
  auto cloud = base;
  const auto nList = nneigh::neighListO(3.5, cloud, 1);
  Snapshot snap;
  snap.neighbours = nList;
  snap.published.resize(static_cast<std::size_t>(cloud.nop));
  snap.ovito.resize(static_cast<std::size_t>(cloud.nop));
  {
    auto labeled = cloud;
    chill::applyChillPlus(labeled, nList, chill::ChillPlusRule::Published);
    for (int i = 0; i < cloud.nop; i++) {
      snap.published[static_cast<std::size_t>(i)] =
          static_cast<int>(labeled.pts[static_cast<std::size_t>(i)].iceType);
    }
  }
  {
    auto labeled = cloud;
    chill::applyChillPlus(labeled, nList, chill::ChillPlusRule::Ovito);
    for (int i = 0; i < cloud.nop; i++) {
      snap.ovito[static_cast<std::size_t>(i)] =
          static_cast<int>(labeled.pts[static_cast<std::size_t>(i)].iceType);
    }
  }
  const auto q = chill::steinhardtQl(cloud, nList, 6);
  snap.q6 = q.ql;
  snap.q6bar = q.qlBar;
  snap.ringList = primitive::ringNetwork(rowsByIndex(cloud, nList), 6);
  snap.rings = census(snap.ringList);
  return snap;
}

} // namespace

TEST_CASE("one, two and four threads match the serial CHILL+, q6 and rings",
          "[thread]") {
  // Above the threaded cell-list threshold, so two and four threads take
  // that walk and must still match the one-thread linkcell rows.
  const auto base = jitteredIce(8);
  const Snapshot serial = runOnce(base, 1);
  REQUIRE(lc_num_threads() == 1);
  for (const int threads : {2, 4}) {
    const Snapshot got = runOnce(base, threads);
    REQUIRE(lc_num_threads() == threads);
    REQUIRE(got.neighbours == serial.neighbours);
    REQUIRE(got.published == serial.published);
    REQUIRE(got.ovito == serial.ovito);
    REQUIRE(got.rings == serial.rings);
    REQUIRE(got.ringList == serial.ringList);
    REQUIRE(got.q6.size() == serial.q6.size());
    for (std::size_t i = 0; i < serial.q6.size(); i++) {
      REQUIRE_THAT(got.q6[i], Catch::Matchers::WithinAbs(serial.q6[i], 1e-15));
      REQUIRE_THAT(got.q6bar[i], Catch::Matchers::WithinAbs(serial.q6bar[i], 1e-15));
    }
  }
}

TEST_CASE("one, two and four ranks list the same rings as one rank", "[thread]") {
  const auto base = jitteredIce(8);
  const auto whole = nneigh::getNewNeighbourListByIndex(base, 3.5);
  const auto serial = primitive::ringNetwork(whole, 6);
  for (const int nRanks : {1, 2, 4, 8}) {
    std::vector<std::vector<int>> joined;
    for (int rank = 0; rank < nRanks; rank++) {
      const auto share = seams::domain::decompose(base, rank, nRanks, 3.5 * 4);
      const auto nbrs = seams::domain::haloRanks(share);
      for (const int other : nbrs) {
        REQUIRE(other != rank);
        REQUIRE(other >= 0);
        REQUIRE(other < nRanks);
      }
      const auto mine = seams::domain::rings(base, 3.5, 6, rank, nRanks);
      joined.insert(joined.end(), mine.begin(), mine.end());
    }
    std::stable_sort(joined.begin(), joined.end(),
                     [](const auto &a, const auto &b) { return a[0] < b[0]; });
    auto expect = serial;
    std::stable_sort(expect.begin(), expect.end(),
                     [](const auto &a, const auto &b) { return a[0] < b[0]; });
    REQUIRE(joined == expect);
  }
}

TEST_CASE("spatial ranks match the serial neighbour list, CHILL+ labels and q6",
          "[thread]") {
  const auto base = jitteredIce(8);
  const double cutoff = 3.5;
  const auto serialRows = nneigh::neighListO(cutoff, base, 1);
  auto serialCloud = base;
  chill::applyChillPlus(serialCloud, serialRows, chill::ChillPlusRule::Published);
  auto serialOvito = base;
  chill::applyChillPlus(serialOvito, serialRows, chill::ChillPlusRule::Ovito);
  const auto serialQ = chill::steinhardtQl(base, serialRows, 6);
  for (const int nRanks : {1, 2, 4, 8}) {
    for (const double halo : {cutoff, 2.0 * cutoff}) {
      for (int rank = 0; rank < nRanks; rank++) {
        const auto share = seams::domain::decompose(base, rank, nRanks, halo);
        const auto local = seams::domain::subCloud(base, share.local);
        const auto rows = nneigh::neighListO(cutoff, local, 1);
        std::vector<std::vector<int>> window;
        const bool blocked =
            nneigh::blockWindowRows(local, cutoff, &share.owned, window);
        if (nRanks > 1) {
          REQUIRE(blocked);
        }
        for (std::size_t k = 0; k < share.local.size(); k++) {
          if (!share.owned[k]) {
            continue;
          }
          const int frame = share.local[k];
          if (halo == cutoff) {
            REQUIRE(rows[k] == serialRows[static_cast<std::size_t>(frame)]);
            if (blocked) {
              REQUIRE(window[k] == serialRows[static_cast<std::size_t>(frame)]);
            }
          }
        }
        // q and c_ij at a neighbour need that neighbour's shell, one
        // cutoff further out. One cutoff is not enough; two is.
        if (halo != 2.0 * cutoff) {
          continue;
        }
        auto published = local;
        chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published);
        auto ovito = local;
        chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito);
        const auto q = chill::steinhardtQl(local, rows, 6);
        for (std::size_t k = 0; k < share.local.size(); k++) {
          if (!share.owned[k]) {
            continue;
          }
          const int frame = share.local[k];
          REQUIRE(published.pts[k].iceType ==
                  serialCloud.pts[static_cast<std::size_t>(frame)].iceType);
          REQUIRE(ovito.pts[k].iceType ==
                  serialOvito.pts[static_cast<std::size_t>(frame)].iceType);
          REQUIRE_THAT(q.ql[k], Catch::Matchers::WithinAbs(
                                    serialQ.ql[static_cast<std::size_t>(frame)], 1e-15));
          REQUIRE_THAT(q.qlBar[k],
                       Catch::Matchers::WithinAbs(
                           serialQ.qlBar[static_cast<std::size_t>(frame)], 1e-15));
        }
      }
    }
  }
}

#ifdef __linux__
namespace {

bool sameMask(const cpu_set_t &a, const cpu_set_t &b) {
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, &a) != CPU_ISSET(cpu, &b)) {
      return false;
    }
  }
  return true;
}

struct AffinityGuard {
  cpu_set_t saved{};
#ifdef SEAMS_HAS_OPENMP
  int threads = 1;
#endif
  AffinityGuard() {
    CPU_ZERO(&saved);
    sched_getaffinity(0, sizeof(saved), &saved);
#ifdef SEAMS_HAS_OPENMP
    threads = omp_get_max_threads();
#endif
  }
  ~AffinityGuard() {
    sched_setaffinity(0, sizeof(saved), &saved);
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(threads);
#endif
  }
};

} // namespace

int maskCount(const cpu_set_t &set) {
  int n = 0;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, &set)) {
      ++n;
    }
  }
  return n;
}

// The pool follows the mask the process was given. Proc-bind pins the
// calling thread to one place, so the place count is that inherited
// width when it is wider than this thread.
int expectedPool(int requested, int maskCpus) {
  int limit = maskCpus;
#ifdef SEAMS_HAS_OPENMP
  const int places = omp_get_num_procs();
  const char *bind = std::getenv("OMP_PROC_BIND");
  const bool procBind = bind != nullptr && bind[0] != '\0' && bind[0] != 'f' &&
                        bind[0] != 'F' && bind[0] != '0';
  if (places > 0 && places < limit) {
    limit = places;
  } else if (procBind && places > limit) {
    limit = places;
  }
#else
  (void)requested;
  return 1;
#endif
  if (limit > 0 && requested > limit) {
    requested = limit;
  }
  return requested < 1 ? 1 : requested;
}

TEST_CASE("a restricted affinity mask is kept and sizes the pool", "[thread]") {
  AffinityGuard guard;
  cpu_set_t saved = guard.saved;
  int allowed[8];
  int nAllowed = 0;
  for (int cpu = 0; cpu < CPU_SETSIZE && nAllowed < 8; ++cpu) {
    if (CPU_ISSET(cpu, &saved)) {
      allowed[nAllowed++] = cpu;
    }
  }
  REQUIRE(nAllowed >= 1);

  auto check = [](const cpu_set_t &want, int requested) {
    REQUIRE(sched_setaffinity(0, sizeof(want), &want) == 0);
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(requested);
#endif
    seams::adoptThreads();
    cpu_set_t after;
    CPU_ZERO(&after);
    REQUIRE(sched_getaffinity(0, sizeof(after), &after) == 0);
    REQUIRE(sameMask(want, after));
    const int expect = expectedPool(requested, maskCount(after));
    REQUIRE(lc_num_threads() == expect);
#ifdef SEAMS_HAS_OPENMP
    REQUIRE(omp_get_max_threads() == expect);
#endif
    // A mask the runtime still reports as the place count sizes the pool.
#ifdef SEAMS_HAS_OPENMP
    if (omp_get_num_procs() == maskCount(after)) {
      const int capped = requested < maskCount(after) ? requested : maskCount(after);
      REQUIRE(lc_num_threads() == (capped < 1 ? 1 : capped));
    }
#endif
  };

  cpu_set_t one;
  CPU_ZERO(&one);
  CPU_SET(allowed[0], &one);
  check(one, 4);

  if (nAllowed >= 2) {
    cpu_set_t two;
    CPU_ZERO(&two);
    CPU_SET(allowed[0], &two);
    CPU_SET(allowed[1], &two);
    check(two, 8);
  }
}
#endif
