/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** In-memory kernels behind the paper's comparison rows. Warmup repetitions
** are not recorded. Mode q6 builds the neighbour list outside its timer;
** q6-full includes that search. Both write the optional sidecar.
*/

#include <bop.hpp>
#include <cage.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

double millisSince(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
}

const char *stateName(molSys::atom_state_type state) {
  switch (state) {
  case molSys::atom_state_type::cubic:
    return "cubic";
  case molSys::atom_state_type::reCubic:
    return "reCubic";
  case molSys::atom_state_type::hexagonal:
    return "hexagonal";
  case molSys::atom_state_type::reHex:
    return "reHex";
  case molSys::atom_state_type::water:
    return "water";
  case molSys::atom_state_type::interfacial:
    return "interfacial";
  case molSys::atom_state_type::clathrate:
    return "clathrate";
  case molSys::atom_state_type::interClathrate:
    return "interClathrate";
  case molSys::atom_state_type::unclassified:
    return "unclassified";
  }
  return "unknown";
}

molSys::PointCloud<molSys::Point<double>, double>
load(const std::string &traj, int frame) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, frame, cloud, 1);
  if (cloud.nop == 0) {
    throw std::runtime_error("could not read " + traj);
  }
  return cloud;
}

void header(const std::string &mode, int atoms, int warmup, int reps) {
  std::cout << "mode " << mode << "\natoms " << atoms << "\nwarmup " << warmup
            << "\nreps " << reps << "\ncolumns kernel_ms\n";
}

void runChill(const std::string &traj, int frame, int warmup, int reps,
              const std::string &sidecar) {
  auto cloud = load(traj, frame);
  header("chill", cloud.nop, warmup, reps);
  for (int rep = 0; rep < warmup + reps; ++rep) {
    auto work = cloud;
    auto t0 = std::chrono::steady_clock::now();
    // 3.5 A is the first-neighbour cutoff. The list is inside the timer
    // so this call is comparable to seededLabels, which builds its lists
    // inside the same kind of timer.
    auto list = nneigh::neighListO(3.5, work, 1);
    chill::getCorrelPlus(work, list, false);
    chill::getIceTypePlusNoPrint(work, list, false);
    const double ms = millisSince(t0);
    if (rep + 1 == warmup + reps) {
      cloud = std::move(work);
    }
    if (rep >= warmup) {
      std::cout << ms << "\n";
    }
  }
  if (!sidecar.empty()) {
    std::ofstream out(sidecar);
    for (int i = 0; i < cloud.nop; ++i) {
      out << stateName(cloud.pts[static_cast<std::size_t>(i)].iceType) << "\n";
    }
  }
}

// The cell edge is the first-shell cutoff, as the command line uses.
// It only sizes the cell list of kNearestNeighbourList. The search
// falls back to brute force when the fourth neighbour lies beyond it,
// so the labels do not depend on this value. tests/test_neighbours.cpp,
// "k-nearest graph is exact under an undersized candidate cutoff".
// The timing does depend on it.
constexpr double kCandidateA = 3.5;

void runQ6(const std::string &traj, int frame, int warmup, int reps,
           bool includeNeighbours, const std::string &sidecar) {
  auto cloud = load(traj, frame);
  std::vector<std::vector<int>> rows;
  if (!includeNeighbours) {
    rows = nneigh::kNearestNeighbourList(cloud, 4, kCandidateA, 1, true);
  }
  header(includeNeighbours ? "q6-full" : "q6", cloud.nop, warmup, reps);
  chill::SteinhardtQl last;
  for (int rep = 0; rep < warmup + reps; ++rep) {
    auto t0 = std::chrono::steady_clock::now();
    if (includeNeighbours) {
      rows = nneigh::kNearestNeighbourList(cloud, 4, kCandidateA, 1, true);
    }
    last = chill::steinhardtQl(cloud, rows, 6);
    const double ms = millisSince(t0);
    if (rep >= warmup) {
      std::cout << ms << "\n";
    }
  }
  if (!sidecar.empty()) {
    std::ofstream out(sidecar);
    out << std::setprecision(17);
    out << "bonds\n";
    for (const auto &row : rows) {
      if (row.empty()) {
        continue;
      }
      const int query = row[0] - 1;
      for (std::size_t k = 1; k < row.size(); ++k) {
        out << query << " " << (row[k] - 1) << "\n";
      }
    }
    out << "values\n";
    for (int i = 0; i < cloud.nop; ++i) {
      out << last.ql[static_cast<std::size_t>(i)] << " "
          << last.qlBar[static_cast<std::size_t>(i)] << "\n";
    }
  }
}

void runRings(const std::string &traj, int frame, int warmup, int reps,
              const std::string &sidecar) {
  auto cloud = load(traj, frame);
  header("rings", cloud.nop, warmup, reps);
  std::vector<std::vector<int>> last;
  for (int rep = 0; rep < warmup + reps; ++rep) {
    auto t0 = std::chrono::steady_clock::now();
    auto neighbours = nneigh::neighListO(3.5, cloud, 1);
    auto indexed = nneigh::neighbourListByIndex(cloud, neighbours);
    last = primitive::ringNetwork(indexed, 6);
    const double ms = millisSince(t0);
    if (rep >= warmup) {
      std::cout << ms << "\n";
    }
  }
  if (!sidecar.empty()) {
    std::array<std::size_t, 7> counts{};
    for (const auto &ring : last) {
      if (ring.size() < counts.size()) {
        ++counts[ring.size()];
      }
    }
    std::ofstream out(sidecar);
    for (std::size_t size = 3; size < counts.size(); ++size) {
      out << size << " " << counts[size] << "\n";
    }
  }
}

std::vector<std::string> seededLabels(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                                      bool ringAdjacent, int depth) {
  // One nomination pass produces both the mutual seed graph and the union
  // completion graph. Building each graph through a separate public call
  // repeats the same nearest-neighbour search.
  auto pair = nneigh::kNearestNeighbourPair(cloud, 4, kCandidateA, 1);
  auto strictRows = std::move(pair.first);
  auto permRows = std::move(pair.second);
  auto strict = nneigh::neighbourListByIndex(cloud, strictRows);
  auto perm = nneigh::neighbourListByIndex(cloud, permRows);
  std::vector<std::vector<int>> sixS;
  std::vector<std::vector<int>> sixP;
  for (const auto &ring : primitive::ringNetwork(strict, depth)) {
    if (ring.size() == 6) {
      sixS.push_back(ring);
    }
  }
  for (const auto &ring : primitive::ringNetwork(perm, depth)) {
    if (ring.size() == 6) {
      sixP.push_back(ring);
    }
  }
  auto labels =
      ring::seededCageAffiliation(sixS, strict, sixP, perm, ringAdjacent);
  std::vector<std::string> out(static_cast<std::size_t>(cloud.nop), "water");
  for (int i = 0; i < cloud.nop; ++i) {
    const bool hc = i < static_cast<int>(labels.hc.size()) && labels.hc[static_cast<std::size_t>(i)];
    const bool ddc = i < static_cast<int>(labels.ddc.size()) && labels.ddc[static_cast<std::size_t>(i)];
    out[static_cast<std::size_t>(i)] =
        hc && ddc ? "mixed" : hc ? "hexagonal" : ddc ? "cubic" : "water";
  }
  return out;
}

void runSeeded(const std::string &traj, int frame, int warmup, int reps,
               bool ringAdjacent, int depth, const std::string &sidecar) {
  auto cloud = load(traj, frame);
  header(ringAdjacent ? "seeded-adj" : "seeded", cloud.nop, warmup, reps);
  std::vector<std::string> last;
  for (int rep = 0; rep < warmup + reps; ++rep) {
    auto t0 = std::chrono::steady_clock::now();
    last = seededLabels(cloud, ringAdjacent, depth);
    const double ms = millisSince(t0);
    if (rep >= warmup) {
      std::cout << ms << "\n";
    }
  }
  if (!sidecar.empty()) {
    std::ofstream out(sidecar);
    for (const auto &label : last) {
      out << label << "\n";
    }
  }
}

// The stages and graph of scripts/baseline_dseams1.cpp: cutoff graph built
// once outside the timers, rings through seven, then findHC and findDDC.
// The ring search index belongs to the cage stage, so it is inside hc_ms.
void runCages(const std::string &traj, int frame, int warmup, int reps) {
  auto cloud = load(traj, frame);
  auto neighbours = nneigh::neighListO(3.5, cloud, 1);
  auto indexed = nneigh::neighbourListByIndex(cloud, neighbours);
  std::cout << "mode cages\natoms " << cloud.nop << "\nwarmup " << warmup
            << "\nreps " << reps << "\ncolumns ring_ms hc_ms ddc_ms path_ms\n";
  std::vector<std::vector<int>> rings;
  std::vector<std::vector<int>> six;
  std::vector<int> listHC;
  std::vector<int> listDDC;
  std::vector<cage::Cage> cages;
  for (int rep = 0; rep < warmup + reps; ++rep) {
    const auto t0 = std::chrono::steady_clock::now();
    rings = primitive::ringNetwork(indexed, 7);
    six.clear();
    for (const auto &ring : rings) {
      if (ring.size() == 6) {
        six.push_back(ring);
      }
    }
    const double ringMs = millisSince(t0);
    const auto t1 = std::chrono::steady_clock::now();
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> ringType(six.size(),
                                          ring::strucType::unclassified);
    cages.clear();
    listHC = ring::findHC(six, ringType, indexed, cages, index);
    const double hcMs = millisSince(t1);
    const auto t2 = std::chrono::steady_clock::now();
    listDDC = ring::findDDC(six, ringType, listHC, cages, index);
    const double ddcMs = millisSince(t2);
    if (rep >= warmup) {
      std::cout << ringMs << " " << hcMs << " " << ddcMs << " "
                << (ringMs + hcMs + ddcMs) << "\n";
    }
  }
  int nHC = 0;
  int nDDC = 0;
  for (const auto &cage : cages) {
    if (cage.type == cage::cageType::HexC) {
      nHC++;
    } else if (cage.type == cage::cageType::DoubleDiaC) {
      nDDC++;
    }
  }
  std::cerr << "rings " << rings.size() << " six_rings " << six.size()
            << " hc_rings " << listHC.size() << " ddc_rings " << listDDC.size()
            << " hc_cages " << nHC << " ddc_cages " << nDDC << "\n";
}

void runCageIds(const std::string &traj, int frame) {
  auto cloud = load(traj, frame);
  auto neighbours = nneigh::neighListO(3.5, cloud, 1);
  auto indexed = nneigh::neighbourListByIndex(cloud, neighbours);
  auto rings = ring::contractibleRings(primitive::ringNetwork(indexed, 6), cloud);
  const auto rank = ring::positionRank(cloud);
  const std::vector<cage::Signature> sigs = {
      cage::Signature::parse("512"),
      cage::Signature::parse("51262"),
      cage::Signature::parse("51264"),
  };
  auto found = cage::findBySignatures(rings, indexed, sigs, rank);
  std::cout << "mode cageids\natoms " << cloud.nop << "\n";
  int id = 0;
  for (std::size_t s = 0; s < sigs.size(); ++s) {
    for (const auto &cage : found[s]) {
      std::cout << "cage " << id << " " << sigs[s].str();
      for (int vertex : cage.vertices) {
        std::cout << " " << vertex;
      }
      std::cout << "\n";
      ++id;
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 6) {
    std::cerr << "usage: bench_paper traj frame warmup reps "
                 "chill|q6|q6-full|rings|seeded|seeded-adj|cages|cageids "
                 "[sidecar]\n";
    return 2;
  }
  const std::string traj = argv[1];
  const int frame = std::atoi(argv[2]);
  const int warmup = std::atoi(argv[3]);
  const int reps = std::atoi(argv[4]);
  const std::string mode = argv[5];
  const std::string sidecar = argc > 6 ? argv[6] : "";
  try {
    if (mode == "chill") {
      runChill(traj, frame, warmup, reps, sidecar);
    } else if (mode == "q6") {
      runQ6(traj, frame, warmup, reps, false, sidecar);
    } else if (mode == "q6-full") {
      runQ6(traj, frame, warmup, reps, true, sidecar);
    } else if (mode == "rings") {
      runRings(traj, frame, warmup, reps, sidecar);
    } else if (mode == "seeded") {
      runSeeded(traj, frame, warmup, reps, false, 6, sidecar);
    } else if (mode == "seeded-adj") {
      runSeeded(traj, frame, warmup, reps, true, 7, sidecar);
    } else if (mode == "cages") {
      runCages(traj, frame, warmup, reps);
    } else if (mode == "cageids") {
      runCageIds(traj, frame);
    } else {
      std::cerr << "unknown mode " << mode << "\n";
      return 2;
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << "\n";
    return 1;
  }
  return 0;
}
