// Driver for the published d-SEAMS 1 library (tag v1.0.1).
// Compiled against that tree's headers, not against the current sources.
// Same input as tests/bench_cages: one LAMMPS frame, cutoff 3.5, rings
// through 7, then findHC and findDDC. Prints one JSON object. With a sixth
// argument (warm-up calls) it prints one line per retained call in the
// columns of `bench_paper ... cages` instead.

#include <cage.hpp>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <franzblau.hpp>
#include <iostream>
#include <limits>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <string>
#include <topo_bulk.hpp>
#include <vector>

namespace {

template <typename Fn>
double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "input/traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int reps = argc > 3 ? std::atoi(argv[3]) : 1;
  const int atomType = argc > 4 ? std::atoi(argv[4]) : 1;
  const std::string commit = argc > 5 ? argv[5] : "";

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(traj, frame, &yCloud, atomType);
  if (yCloud.nop == 0) {
    std::cerr << "could not read " << traj << "\n";
    return 1;
  }

  auto nList = nneigh::neighListO(3.5, &yCloud, atomType);
  auto hbond = nneigh::neighbourListByIndex(&yCloud, nList);

  if (argc > 6) {
    const int warmup = std::atoi(argv[6]);
    const auto since = [](std::chrono::steady_clock::time_point start) {
      return std::chrono::duration<double, std::milli>(
                 std::chrono::steady_clock::now() - start)
          .count();
    };
    std::cout << "mode cages\natoms " << yCloud.nop << "\nwarmup " << warmup
              << "\nreps " << reps
              << "\ncolumns ring_ms hc_ms ddc_ms path_ms\n";
    std::vector<std::vector<int>> rings;
    std::vector<std::vector<int>> six;
    std::vector<int> listHC;
    std::vector<int> listDDC;
    std::vector<cage::Cage> cages;
    for (int rep = 0; rep < warmup + reps; rep++) {
      const auto t0 = std::chrono::steady_clock::now();
      rings = primitive::ringNetwork(hbond, 7);
      six.clear();
      for (const auto &ring : rings) {
        if (ring.size() == 6) {
          six.push_back(ring);
        }
      }
      const double ringMs = since(t0);
      const auto t1 = std::chrono::steady_clock::now();
      std::vector<ring::strucType> ringType(six.size(), ring::unclassified);
      cages.clear();
      listHC = ring::findHC(six, &ringType, hbond, &cages);
      const double hcMs = since(t1);
      const auto t2 = std::chrono::steady_clock::now();
      listDDC = ring::findDDC(six, &ringType, &cages);
      const double ddcMs = since(t2);
      if (rep >= warmup) {
        std::cout << ringMs << " " << hcMs << " " << ddcMs << " "
                  << (ringMs + hcMs + ddcMs) << "\n";
      }
    }
    int nHC = 0;
    int nDDC = 0;
    for (const auto &cage : cages) {
      if (cage.type == cage::HexC) {
        nHC++;
      } else if (cage.type == cage::DoubleDiaC) {
        nDDC++;
      }
    }
    std::cerr << "rings " << rings.size() << " six_rings " << six.size()
              << " hc_rings " << listHC.size() << " ddc_rings "
              << listDDC.size() << " hc_cages " << nHC << " ddc_cages "
              << nDDC << "\n";
    return 0;
  }

  std::vector<std::vector<int>> rings;
  const double tRings = bestMillis(
      [&]() { rings = primitive::ringNetwork(hbond, 7); }, reps);

  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }

  std::vector<ring::strucType> ringType;
  std::vector<cage::Cage> cages;
  std::vector<int> listHC;
  const double tHC = bestMillis(
      [&]() {
        ringType.assign(six.size(), ring::unclassified);
        cages.clear();
        listHC = ring::findHC(six, &ringType, hbond, &cages);
      },
      reps);

  // bestMillis leaves the last repetition's HC state in place.
  const auto ringTypeAfter = ringType;
  const auto cagesAfter = cages;
  std::vector<int> listDDC;
  const double tDDC = bestMillis(
      [&]() {
        ringType = ringTypeAfter;
        cages = cagesAfter;
        listDDC = ring::findDDC(six, &ringType, &cages);
      },
      reps);

  int nHC = 0;
  int nDDC = 0;
  for (const auto &cage : cages) {
    if (cage.type == cage::HexC) {
      nHC++;
    } else if (cage.type == cage::DoubleDiaC) {
      nDDC++;
    }
  }

  std::cout << "{\n"
            << "  \"baseline_commit\": \"" << commit << "\",\n"
            << "  \"baseline_tag\": \"v1.0.1\",\n"
            << "  \"input\": \"" << traj << "\",\n"
            << "  \"frame\": " << frame << ",\n"
            << "  \"reps\": " << reps << ",\n"
            << "  \"atoms\": " << yCloud.nop << ",\n"
            << "  \"rings\": " << rings.size() << ",\n"
            << "  \"six_rings\": " << six.size() << ",\n"
            << "  \"hc_rings\": " << listHC.size() << ",\n"
            << "  \"ddc_rings\": " << listDDC.size() << ",\n"
            << "  \"cages\": " << cages.size() << ",\n"
            << "  \"hc_cages\": " << nHC << ",\n"
            << "  \"ddc_cages\": " << nDDC << ",\n"
            << "  \"ring_ms\": " << tRings << ",\n"
            << "  \"findHC_ms\": " << tHC << ",\n"
            << "  \"findDDC_ms\": " << tDDC << ",\n"
            << "  \"path_ms\": " << (tRings + tHC + tDDC) << "\n"
            << "}\n";
  return 0;
}
