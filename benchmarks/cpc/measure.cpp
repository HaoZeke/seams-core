/*
** Serial measurements for the d-SEAMS 2.0 results section.
**
**   measure_cpc labels chill|chillplus FILE FRAME TYPE
**   measure_cpc time MODE FILE FRAME TYPE WARMUP REPEATS
**
** MODE is chill, chillplus, q3, q6, q6avg, cages, cluster-q6 or cluster-chill.
** The cutoff graph is built once outside the cage timers. CHILL and the
** bond-order modes report a kernel time (neighbours already built) and a
** path time (neighbours inside the timer). Cluster sizes use the same ice
** mask as clump::clusterAnalysis: q6 > 0.5, or CHILL labels other than
** water and unclassified. One thread.
*/

#include <bop.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

volatile std::size_t gSink = 0;

double millis(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

void pinOneThread() {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
}

Cloud loadFrame(const std::string &path, int frame, int type) {
  Cloud cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, type);
  return cloud;
}

void chillLabels(Cloud &cloud, const std::vector<std::vector<int>> &nList,
                 bool plus) {
  if (plus) {
    chill::getCorrelPlus(cloud, nList, false);
    chill::getIceTypePlusNoPrint(cloud, nList, false);
  } else {
    chill::getCorrel(cloud, nList, false);
    chill::getIceTypeNoPrint(cloud, nList, false);
  }
  gSink += static_cast<std::size_t>(cloud.nop);
}

std::vector<bool> iceMaskQ6(Cloud &cloud,
                            const std::vector<std::vector<int>> &nList) {
  const auto q6 = chill::getq6(cloud, nList, false);
  std::vector<bool> mask(static_cast<std::size_t>(cloud.nop), false);
  for (int i = 0; i < cloud.nop; i++) {
    if (q6[static_cast<std::size_t>(i)] > 0.5) {
      mask[static_cast<std::size_t>(i)] = true;
    }
  }
  gSink += q6.size();
  return mask;
}

std::vector<bool> iceMaskChill(Cloud &cloud,
                               const std::vector<std::vector<int>> &nList) {
  chillLabels(cloud, nList, false);
  std::vector<bool> mask(static_cast<std::size_t>(cloud.nop), false);
  for (int i = 0; i < cloud.nop; i++) {
    const auto kind = cloud.pts[static_cast<std::size_t>(i)].iceType;
    if (kind == molSys::atom_state_type::water ||
        kind == molSys::atom_state_type::unclassified) {
      continue;
    }
    mask[static_cast<std::size_t>(i)] = true;
  }
  return mask;
}

void printDomain(const Cloud &cloud, const std::vector<std::vector<int>> &nList,
                 const std::vector<bool> &mask) {
  int ice = 0;
  for (bool bit : mask) {
    ice += bit ? 1 : 0;
  }
  const auto domain = clump::largestDomain(cloud, nList, mask);
  std::cout << "count ice_like " << ice << "\n"
            << "count largest " << domain.largest << "\n"
            << "count subset " << domain.subset << "\n"
            << "count percolation " << domain.percolation << "\n";
  gSink += static_cast<std::size_t>(domain.largest);
}

struct CageCount {
  int atoms = 0;
  int rings = 0;
  int six = 0;
  int hcRings = 0;
  int ddcRings = 0;
  int hcCages = 0;
  int ddcCages = 0;
  int equators = 0;
};

std::vector<std::vector<int>>
sixRings(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  return six;
}

int countEquators(const std::vector<std::vector<int>> &six,
                  const std::vector<int> &listHC,
                  const ring::RingSearchIndex &index) {
  std::vector<char> blocked(six.size(), 0);
  for (int id : listHC) {
    if (id >= 0 && id < static_cast<int>(blocked.size())) {
      blocked[static_cast<std::size_t>(id)] = 1;
    }
  }
  int equators = 0;
  for (int iring = 0; iring < static_cast<int>(six.size()); iring++) {
    if (blocked[static_cast<std::size_t>(iring)]) {
      continue;
    }
    if (ring::equatorialConditions(six, iring, index, nullptr)) {
      equators++;
    }
  }
  return equators;
}

void countCages(const std::vector<cage::Cage> &cageList, int &hc, int &ddc) {
  hc = 0;
  ddc = 0;
  for (const auto &cage : cageList) {
    if (cage.type == cage::cageType::HexC) {
      hc++;
    } else if (cage.type == cage::cageType::DoubleDiaC) {
      ddc++;
    }
  }
}

CageCount referenceCages(Cloud &cloud, int type) {
  auto nList = nneigh::neighListO(3.5, cloud, type);
  auto bonds = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(bonds, 7);
  auto six = sixRings(rings);
  const auto index = ring::buildRingSearchIndex(six, cloud);
  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cageList;
  const auto listHC = ring::findHC(six, ringType, bonds, cageList, index);
  const int equators = countEquators(six, listHC, index);
  const auto listDDC = ring::findDDC(six, ringType, listHC, cageList, index);
  CageCount out;
  out.atoms = cloud.nop;
  out.rings = static_cast<int>(rings.size());
  out.six = static_cast<int>(six.size());
  out.hcRings = static_cast<int>(listHC.size());
  out.ddcRings = static_cast<int>(listDDC.size());
  out.equators = equators;
  countCages(cageList, out.hcCages, out.ddcCages);
  gSink += static_cast<std::size_t>(out.hcCages + out.ddcCages);
  return out;
}

void timeCages(Cloud &cloud, int type, int warmup, int repeats) {
  auto nList = nneigh::neighListO(3.5, cloud, type);
  auto bonds = nneigh::neighbourListByIndex(cloud, nList);
  const CageCount counts = referenceCages(cloud, type);
  std::cout << "count atoms " << counts.atoms << "\n"
            << "count rings " << counts.rings << "\n"
            << "count six_rings " << counts.six << "\n"
            << "count hc_rings " << counts.hcRings << "\n"
            << "count ddc_rings " << counts.ddcRings << "\n"
            << "count ddc_equators " << counts.equators << "\n"
            << "count hc_cages " << counts.hcCages << "\n"
            << "count ddc_cages " << counts.ddcCages << "\n";
  auto once = [&]() {
    auto rings = primitive::ringNetwork(bonds, 7);
    auto six = sixRings(rings);
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> ringType(six.size(),
                                          ring::strucType::unclassified);
    std::vector<cage::Cage> cageList;
    auto listHC = ring::findHC(six, ringType, bonds, cageList, index);
    auto listDDC = ring::findDDC(six, ringType, listHC, cageList, index);
    gSink += rings.size() + listHC.size() + listDDC.size();
  };
  for (int i = 0; i < warmup; i++) {
    once();
  }
  for (int rep = 1; rep <= repeats; rep++) {
    const auto t0 = Clock::now();
    auto rings = primitive::ringNetwork(bonds, 7);
    auto six = sixRings(rings);
    const auto t1 = Clock::now();
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> ringType(six.size(),
                                          ring::strucType::unclassified);
    std::vector<cage::Cage> cageList;
    auto listHC = ring::findHC(six, ringType, bonds, cageList, index);
    const auto t2 = Clock::now();
    auto listDDC = ring::findDDC(six, ringType, listHC, cageList, index);
    const auto t3 = Clock::now();
    gSink += rings.size() + listHC.size() + listDDC.size();
    std::cout << "repeat " << rep << " rings_ms " << millis(t0, t1) << "\n"
              << "repeat " << rep << " hc_ms " << millis(t1, t2) << "\n"
              << "repeat " << rep << " ddc_ms " << millis(t2, t3) << "\n"
              << "repeat " << rep << " path_ms " << millis(t0, t3) << "\n";
  }
}

void timeOrder(Cloud &cloud, int type, int order, bool averaged, int warmup,
               int repeats) {
  auto nList = nneigh::neighListO(3.5, cloud, type);
  auto kernel = [&]() {
    if (averaged) {
      auto q = chill::getq6(cloud, nList, false);
      gSink += q.size();
    } else {
      auto q = chill::steinhardtQl(cloud, nList, order);
      gSink += q.ql.size();
    }
  };
  auto path = [&]() {
    auto local = nneigh::neighListO(3.5, cloud, type);
    if (averaged) {
      auto q = chill::getq6(cloud, local, false);
      gSink += q.size();
    } else {
      auto q = chill::steinhardtQl(cloud, local, order);
      gSink += q.ql.size();
    }
  };
  for (int i = 0; i < warmup; i++) {
    kernel();
    path();
  }
  std::cout << "count atoms " << cloud.nop << "\n";
  for (int rep = 1; rep <= repeats; rep++) {
    const auto k0 = Clock::now();
    kernel();
    const auto k1 = Clock::now();
    const auto p0 = Clock::now();
    path();
    const auto p1 = Clock::now();
    std::cout << "repeat " << rep << " kernel_ms " << millis(k0, k1) << "\n"
              << "repeat " << rep << " path_ms " << millis(p0, p1) << "\n";
  }
}

void timeChill(Cloud &cloud, int type, bool plus, int warmup, int repeats) {
  auto nList = nneigh::neighListO(3.5, cloud, type);
  auto kernel = [&]() { chillLabels(cloud, nList, plus); };
  auto path = [&]() {
    auto local = nneigh::neighListO(3.5, cloud, type);
    chillLabels(cloud, local, plus);
  };
  for (int i = 0; i < warmup; i++) {
    kernel();
    path();
  }
  std::cout << "count atoms " << cloud.nop << "\n";
  for (int rep = 1; rep <= repeats; rep++) {
    const auto k0 = Clock::now();
    kernel();
    const auto k1 = Clock::now();
    const auto p0 = Clock::now();
    path();
    const auto p1 = Clock::now();
    std::cout << "repeat " << rep << " kernel_ms " << millis(k0, k1) << "\n"
              << "repeat " << rep << " path_ms " << millis(p0, p1) << "\n";
  }
}

void timeCluster(Cloud &cloud, int type, bool byChill, int warmup, int repeats) {
  auto run = [&]() {
    auto nList = nneigh::neighListO(3.5, cloud, type);
    const auto mask =
        byChill ? iceMaskChill(cloud, nList) : iceMaskQ6(cloud, nList);
    const auto domain = clump::largestDomain(cloud, nList, mask);
    gSink += static_cast<std::size_t>(domain.largest);
    return std::pair<clump::Domain, int>(domain, [&]() {
      int ice = 0;
      for (bool bit : mask) {
        ice += bit ? 1 : 0;
      }
      return ice;
    }());
  };
  for (int i = 0; i < warmup; i++) {
    (void)run();
  }
  for (int rep = 1; rep <= repeats; rep++) {
    const auto t0 = Clock::now();
    const auto found = run();
    const auto t1 = Clock::now();
    std::cout << "repeat " << rep << " path_ms " << millis(t0, t1) << "\n"
              << "repeat " << rep << " largest " << found.first.largest << "\n"
              << "repeat " << rep << " ice_like " << found.second << "\n"
              << "repeat " << rep << " subset " << found.first.subset << "\n";
  }
}

void printLabels(Cloud &cloud, int type, bool plus) {
  auto nList = nneigh::neighListO(3.5, cloud, type);
  chillLabels(cloud, nList, plus);
  for (int i = 0; i < cloud.nop; i++) {
    std::cout << static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType)
              << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  pinOneThread();
  if (argc < 2) {
    std::cerr << "usage: measure_cpc labels|time ...\n";
    return 2;
  }
  const std::string action = argv[1];
  if (action == "labels") {
    if (argc != 6) {
      std::cerr << "usage: measure_cpc labels chill|chillplus FILE FRAME TYPE\n";
      return 2;
    }
    const std::string mode = argv[2];
    Cloud cloud = loadFrame(argv[3], std::atoi(argv[4]), std::atoi(argv[5]));
    if (cloud.nop <= 0) {
      std::cerr << "nop 0\n";
      return 2;
    }
    if (mode == "chill") {
      printLabels(cloud, std::atoi(argv[5]), false);
    } else if (mode == "chillplus") {
      printLabels(cloud, std::atoi(argv[5]), true);
    } else {
      std::cerr << "unknown label mode\n";
      return 2;
    }
    return 0;
  }
  if (action != "time" || argc != 8) {
    std::cerr << "usage: measure_cpc time MODE FILE FRAME TYPE WARMUP REPEATS\n";
    return 2;
  }
  const std::string mode = argv[2];
  const int frame = std::atoi(argv[4]);
  const int type = std::atoi(argv[5]);
  const int warmup = std::atoi(argv[6]);
  const int repeats = std::atoi(argv[7]);
  Cloud cloud = loadFrame(argv[3], frame, type);
  if (cloud.nop <= 0) {
    std::cerr << "nop 0\n";
    return 2;
  }
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  std::cout << "meta mode " << mode << "\n"
            << "meta frame " << frame << "\n"
            << "meta type " << type << "\n"
            << "meta warmup " << warmup << "\n"
            << "meta repeats " << repeats << "\n"
            << "meta threads 1\n"
            << "meta cutoff 3.5\n";
  if (mode == "cages") {
    std::cout << "meta maxring 7\n";
    timeCages(cloud, type, warmup, repeats);
  } else if (mode == "q3") {
    timeOrder(cloud, type, 3, false, warmup, repeats);
  } else if (mode == "q6") {
    timeOrder(cloud, type, 6, false, warmup, repeats);
  } else if (mode == "q6avg") {
    timeOrder(cloud, type, 6, true, warmup, repeats);
  } else if (mode == "chill") {
    timeChill(cloud, type, false, warmup, repeats);
  } else if (mode == "chillplus") {
    timeChill(cloud, type, true, warmup, repeats);
  } else if (mode == "cluster-q6") {
    timeCluster(cloud, type, false, warmup, repeats);
  } else if (mode == "cluster-chill") {
    timeCluster(cloud, type, true, warmup, repeats);
  } else {
    std::cerr << "unknown mode\n";
    return 2;
  }
  std::cout << "meta sink " << gSink << "\n";
  return 0;
}
