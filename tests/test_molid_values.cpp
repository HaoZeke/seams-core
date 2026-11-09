// Copyright (c) d-SEAMS developers

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bond.hpp>
#include <cluster.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <seams_input.hpp>
#include <site.hpp>

#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

void addPoint(Cloud &cloud, int id, int mol, int type, double x, double y, double z) {
  molSys::Point<double> p;
  p.atomID = id;
  p.molID = mol;
  p.type = type;
  p.x = x;
  p.y = y;
  p.z = z;
  cloud.pts.push_back(p);
  cloud.idIndexMap[id] = static_cast<int>(cloud.pts.size()) - 1;
  cloud.nop = static_cast<int>(cloud.pts.size());
}

Cloud blank() {
  Cloud cloud;
  cloud.box = {30.0, 30.0, 30.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  return cloud;
}

Cloud waterDimer(bool sharedMol) {
  Cloud cloud = blank();
  const int m1 = sharedMol ? 1 : 1;
  const int m2 = sharedMol ? 1 : 2;
  const int m3 = sharedMol ? 1 : 3;
  const int m4 = sharedMol ? 2 : 4;
  const int m5 = sharedMol ? 2 : 5;
  const int m6 = sharedMol ? 2 : 6;
  addPoint(cloud, 1, m1, 2, 5.0, 5.0, 5.0);
  addPoint(cloud, 2, m2, 1, 5.96, 5.0, 5.0);
  addPoint(cloud, 3, m3, 1, 4.5, 5.7, 5.0);
  addPoint(cloud, 4, m4, 2, 7.6, 5.0, 5.0);
  addPoint(cloud, 5, m5, 1, 8.4, 5.4, 5.0);
  addPoint(cloud, 6, m6, 1, 8.4, 4.6, 5.0);
  return cloud;
}

Cloud byType(const Cloud &src, int type) {
  Cloud out = blank();
  out.box = src.box;
  out.boxLow = src.boxLow;
  for (const auto &p : src.pts) {
    if (p.type == type) {
      addPoint(out, p.atomID, p.molID, p.type, p.x, p.y, p.z);
    }
  }
  return out;
}

int hbondCount(Cloud &full) {
  Cloud oxy = byType(full, 2);
  Cloud hyd = byType(full, 1);
  auto nList = nneigh::neighListO(3.5, oxy, 2);
  auto net = bond::populateHbondsWithInputClouds(oxy, hyd, nList, 2.42, 30.0);
  int edges = 0;
  for (const auto &row : net) {
    if (row.size() > 1) {
      edges += static_cast<int>(row.size()) - 1;
    }
  }
  return edges / 2;
}

int finiteF4(const Cloud &full, double &mean) {
  auto nList = nneigh::kNearestNeighbourList(full, 4, 5.0, 2, true);
  const auto f4 = topoparam::rodgerF4(full, nList, 2, 1);
  int n = 0;
  for (double v : f4) {
    n += std::isfinite(v) ? 1 : 0;
  }
  mean = topoparam::meanFinite(f4);
  return n;
}

site::Table ionTable() {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;
  return table;
}

Cloud ionPair(bool sharedMol) {
  Cloud cloud = blank();
  addPoint(cloud, 1, sharedMol ? 10 : 1, 1, 0.0, 0.0, 0.0);
  addPoint(cloud, 2, sharedMol ? 10 : 2, 1, 10.0, 0.0, 0.0);
  addPoint(cloud, 3, sharedMol ? 11 : 3, 2, 1.0, 0.0, 0.0);
  addPoint(cloud, 4, sharedMol ? 11 : 4, 2, 11.0, 0.0, 0.0);
  return cloud;
}

} // namespace

TEST_CASE("shared molID keeps one hydrogen bond on the dimer", "[molid]") {
  Cloud shared = waterDimer(true);
  Cloud split = waterDimer(false);
  REQUIRE(hbondCount(shared) == 1);
  REQUIRE(hbondCount(split) == 0);
}

TEST_CASE("shared molID keeps a finite Rodger F4 on the dimer", "[molid]") {
  Cloud shared = waterDimer(true);
  Cloud split = waterDimer(false);
  double sharedMean = 0.0;
  double splitMean = 0.0;
  REQUIRE(finiteF4(shared, sharedMean) == 2);
  REQUIRE(finiteF4(split, splitMean) == 0);
  REQUIRE_THAT(sharedMean, Catch::Matchers::WithinAbs(1.0, 1e-9));
}

TEST_CASE("shared molID collapses each ion to one site", "[molid]") {
  const auto table = ionTable();
  const auto shared = site::ionCloud(ionPair(true), table);
  const auto split = site::ionCloud(ionPair(false), table);
  REQUIRE(shared.nop == 2);
  REQUIRE(split.nop == 4);
  REQUIRE(nneigh::mutualNearestUnlike(shared, 1, 2).size() == 1);
  REQUIRE(nneigh::mutualNearestUnlike(split, 1, 2).size() == 2);
}

TEST_CASE("shared molID changes the ion-site domain", "[molid]") {
  const auto table = ionTable();
  const auto shared = site::ionCloud(ionPair(true), table);
  const auto split = site::ionCloud(ionPair(false), table);
  auto idNeighbours = [](const Cloud &cloud, double cutoff) {
    const auto idxList = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
    std::vector<std::vector<int>> idList(idxList.size());
    for (std::size_t i = 0; i < idxList.size(); ++i) {
      for (const int idx : idxList[i]) {
        if (idx >= 0 && idx < cloud.nop) {
          idList[i].push_back(cloud.pts[static_cast<std::size_t>(idx)].atomID);
        }
      }
    }
    return idList;
  };
  const std::vector<bool> sharedMask(static_cast<std::size_t>(shared.nop), true);
  const std::vector<bool> splitMask(static_cast<std::size_t>(split.nop), true);
  const auto sharedDomain = clump::largestDomain(shared, idNeighbours(shared, 6.0), sharedMask);
  const auto splitDomain = clump::largestDomain(split, idNeighbours(split, 6.0), splitMask);
  REQUIRE(sharedDomain.subset == 2);
  REQUIRE(sharedDomain.largest == 2);
  REQUIRE(splitDomain.subset == 4);
  REQUIRE(splitDomain.largest == 2);
}

TEST_CASE("the native LAMMPS reader keeps exampleTraj molecule ids", "[molid]") {
  Cloud cloud;
  cloud = sinp::readLammpsTrj("traj/exampleTraj.lammpstrj", 1, cloud);
  REQUIRE(cloud.nop == 750);
  std::set<int> mols;
  int mol1 = 0;
  for (const auto &p : cloud.pts) {
    mols.insert(p.molID);
    mol1 += p.molID == 1 ? 1 : 0;
  }
  REQUIRE(mols.size() == 250);
  REQUIRE(mol1 == 3);
  REQUIRE(cloud.pts[0].atomID == 1);
  REQUIRE(cloud.pts[1].molID == cloud.pts[0].molID);
}
