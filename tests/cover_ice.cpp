//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

// Internal
#include <bop.hpp>
#include <bulkTUM.hpp>
#include <cage.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <tum_offload.hpp>
#include <voronoi_qlm.hpp>

// Standard
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace fs = std::filesystem;

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud makeCloud(const std::vector<std::array<double, 3>> &xyz,
                const std::vector<double> &box, int frame = 1) {
  Cloud cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = frame;
  for (int i = 0; i < static_cast<int>(xyz.size()); i++) {
    molSys::Point<double> point;
    point.type = 1;
    point.atomID = i;
    point.molID = i;
    point.x = xyz[static_cast<std::size_t>(i)][0];
    point.y = xyz[static_cast<std::size_t>(i)][1];
    point.z = xyz[static_cast<std::size_t>(i)][2];
    point.inSlice = true;
    cloud.idIndexMap[i] = i;
    cloud.pts.push_back(point);
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

molSys::Result bondOf(molSys::bond_type kind, double value) {
  molSys::Result bond;
  bond.classifier = kind;
  bond.c_value = value;
  return bond;
}

void pushBonds(molSys::Point<double> &point, int staggered, int eclipsed,
               int other, double value) {
  point.c_ij.clear();
  for (int i = 0; i < staggered; i++) {
    point.c_ij.push_back(bondOf(molSys::bond_type::staggered, value));
  }
  for (int i = 0; i < eclipsed; i++) {
    point.c_ij.push_back(bondOf(molSys::bond_type::eclipsed, value));
  }
  for (int i = 0; i < other; i++) {
    point.c_ij.push_back(bondOf(molSys::bond_type::out_of_range, value));
  }
}

std::vector<std::vector<int>> selfRows(int n) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)] = {i};
  }
  return rows;
}

void link(std::vector<std::vector<int>> &rows, int a, int b) {
  rows[static_cast<std::size_t>(a)].push_back(b);
  rows[static_cast<std::size_t>(b)].push_back(a);
}

void unlink(std::vector<std::vector<int>> &rows, int a, int b) {
  auto &ra = rows[static_cast<std::size_t>(a)];
  auto &rb = rows[static_cast<std::size_t>(b)];
  ra.erase(std::remove(ra.begin(), ra.end(), b), ra.end());
  rb.erase(std::remove(rb.begin(), rb.end(), a), rb.end());
}

Cloud squarePrism(double height, const std::vector<double> &box) {
  return makeCloud({{0, 0, 0},
                    {2.75, 0, 0},
                    {2.75, 2.75, 0},
                    {0, 2.75, 0},
                    {0, 0, height},
                    {2.75, 0, height},
                    {2.75, 2.75, height},
                    {0, 2.75, height}},
                   box);
}

Cloud hexCage() {
  return makeCloud(
      {{8.995, 10.3859997, 15.0939999},
       {6.7459998, 14.2810001, 15.0939999},
       {4.4970002, 10.3859997, 15.0939999},
       {8.995, 12.9829998, 14.1949997},
       {6.7459998, 9.0880003, 14.1949997},
       {4.4970002, 12.9829998, 14.1949997},
       {8.995, 12.9829998, 11.4329996},
       {6.7459998, 9.0880003, 11.4329996},
       {4.4970002, 12.9829998, 11.4329996},
       {8.995, 10.3859997, 10.5340004},
       {6.7459998, 14.2810001, 10.5340004},
       {4.4970002, 10.3859997, 10.5340004}},
      {53.9690018, 54.5289994, 51.257});
}

fs::path scratch(const char *name) {
  auto dir = fs::temp_directory_path() / name;
  fs::remove_all(dir);
  fs::create_directories(dir);
  return dir;
}

} // namespace

TEST_CASE("spherical harmonics accept the degrees outside the hard-coded tables",
          "[cover_ice][bop]") {
  const std::array<double, 2> angles = {0.4, 0.8};
  int degree = 5;
  auto y5 = sph::spheriHarmo(degree, angles);
  REQUIRE(y5.size() == 11);
  REQUIRE(std::abs(y5[5]) > 0.0);

  degree = 8;
  auto y8 = sph::spheriHarmo(degree, angles);
  REQUIRE(y8.size() == 17);
  auto y8Lookup = sph::lookupTableQ8Vec(angles);
  REQUIRE(y8Lookup.size() == 17);
  REQUIRE_THAT(y8[8].real(),
               Catch::Matchers::WithinAbs(y8Lookup[8].real(), 1e-8));

  auto y0 = sph::spheriHarmo(0, angles);
  REQUIRE(y0.size() == 1);
  REQUIRE_THAT(y0[0].real(),
               Catch::Matchers::WithinAbs(0.5 / std::sqrt(M_PI), 1e-8));

  auto empty = sph::spheriHarmo(-1, angles);
  REQUIRE(empty.empty());

  REQUIRE(sph::lookupTableQ3(-1, angles) == std::complex<double>{0.0, 0.0});
  REQUIRE(sph::lookupTableQ6(13, angles) == std::complex<double>{0.0, 0.0});
  REQUIRE(sph::lookupTableQ4(9, angles) == std::complex<double>{0.0, 0.0});
  REQUIRE(sph::lookupTableQ8(17, angles) == std::complex<double>{0.0, 0.0});
}

TEST_CASE("CHILL classifies partial eclipsed counts and ignores a bad neighbour index",
          "[cover_ice][bop]") {
  Cloud cloud = makeCloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                          {20, 20, 20});
  std::vector<std::vector<int>> nList = selfRows(4);

  pushBonds(cloud.pts[0], 4, 0, 0, -0.9);
  pushBonds(cloud.pts[1], 3, 1, 0, -0.5);
  pushBonds(cloud.pts[2], 2, 1, 1, 0.0);
  pushBonds(cloud.pts[3], 0, 0, 2, 0.2);
  chill::getIceTypeNoPrint(cloud, nList, false);
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::cubic);
  REQUIRE(cloud.pts[1].iceType == molSys::atom_state_type::hexagonal);
  REQUIRE(cloud.pts[2].iceType == molSys::atom_state_type::water);
  REQUIRE(cloud.pts[3].iceType == molSys::atom_state_type::water);

  pushBonds(cloud.pts[0], 0, 4, 0, 0.1);
  pushBonds(cloud.pts[1], 0, 3, 1, 0.1);
  pushBonds(cloud.pts[2], 1, 1, 2, 0.0);
  pushBonds(cloud.pts[3], 4, 0, 0, -0.9);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::clathrate);
  REQUIRE(cloud.pts[1].iceType == molSys::atom_state_type::interClathrate);
  REQUIRE(cloud.pts[3].iceType == molSys::atom_state_type::cubic);

  // nop stops short of the stored points, so a mapped index is past nop
  // without being past the point vector the distance uses.
  cloud.nop = 2;
  cloud.idIndexMap[8] = 3;
  std::vector<std::vector<int>> stale = {{0, 8}, {1}, {2}, {3}};
  REQUIRE_FALSE(chill::isInterfacial(cloud, stale, 0, 2, 0, false));
  REQUIRE_FALSE(chill::isInterfacial(cloud, stale, 0, 3, 1, false));
  REQUIRE_FALSE(chill::isInterfacial(cloud, stale, 0, 3, 0, true));

  cloud.nop = 4;
  cloud.pts[2].inSlice = false;
  cloud.pts[2].iceType = molSys::atom_state_type::clathrate;
  chill::getIceTypeNoPrint(cloud, nList, true);
  REQUIRE(cloud.pts[2].iceType == molSys::atom_state_type::clathrate);

  auto printed = scratch("dseams_cover_ice_chill");
  Cloud report = makeCloud({{0, 0, 0}, {1, 0, 0}}, {10, 10, 10}, 1);
  pushBonds(report.pts[0], 2, 1, 1, 0.0);
  pushBonds(report.pts[1], 0, 4, 0, 0.2);
  chill::getIceType(report, selfRows(2), printed.string() + "/", 1, false, "chill.txt");
  REQUIRE(report.pts[0].iceType == molSys::atom_state_type::water);
  REQUIRE(report.pts[1].iceType == molSys::atom_state_type::water);
  fs::remove_all(printed);
}

TEST_CASE("q6 of an isolated atom is zero and a low q3 reclassifies water",
          "[cover_ice][bop]") {
  Cloud cloud = makeCloud({{5, 5, 5}, {5, 6.2, 5}, {12, 5, 5}}, {20, 20, 20});
  std::vector<std::vector<int>> nList = {{0, 1}, {1, 0}, {2}};
  auto q6 = chill::getq6(cloud, nList, false);
  REQUIRE(q6.size() == 3);
  REQUIRE_THAT(q6[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE(std::isfinite(q6[0]));

  cloud.pts[0].iceType = molSys::atom_state_type::water;
  cloud.pts[1].iceType = molSys::atom_state_type::water;
  cloud.pts[2].iceType = molSys::atom_state_type::water;
  pushBonds(cloud.pts[0], 0, 0, 4, -0.9);
  pushBonds(cloud.pts[1], 0, 0, 4, -0.8);
  pushBonds(cloud.pts[2], 0, 0, 4, -0.2);
  q6 = {0.7, 0.6, 0.9};
  chill::reclassifyWater(cloud, q6);
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::reCubic);
  REQUIRE(cloud.pts[1].iceType == molSys::atom_state_type::reHex);
  REQUIRE(cloud.pts[2].iceType == molSys::atom_state_type::water);

  chill::clearIceTypesExcept(cloud, {});
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::reCubic);
  chill::clearIceTypesExcept(cloud, {2});
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::unclassified);
}

TEST_CASE("bond classification accepts non-positive and out-of-range atom ids",
          "[cover_ice][bop]") {
  Cloud sparse = makeCloud({{1, 1, 1}, {2.2, 1, 1}, {1, 2.2, 1}}, {15, 15, 15});
  sparse.pts[0].atomID = 0;
  sparse.pts[1].atomID = -4;
  sparse.pts[2].atomID = 0;
  sparse.idIndexMap.clear();
  sparse.idIndexMap[0] = 0;
  sparse.idIndexMap[-4] = 1;
  chill::getCorrel(sparse, {{0, -4}, {1, 0}, {2}}, false, 4);
  REQUIRE(sparse.pts[0].c_ij.size() == 1);
  REQUIRE(std::isfinite(sparse.pts[0].c_ij[0].c_value));

  Cloud dense = makeCloud({{1, 1, 1}, {2.4, 1, 1}, {1, 2.4, 1}, {1, 1, 2.4}},
                          {15, 15, 15});
  dense.pts[0].atomID = 3;
  dense.pts[1].atomID = -2;
  dense.pts[2].atomID = 1;
  dense.pts[3].atomID = 4;
  dense.idIndexMap.clear();
  dense.idIndexMap[3] = 0;
  dense.idIndexMap[-2] = 1;
  dense.idIndexMap[1] = 2;
  dense.idIndexMap[4] = 3;
  // 99 is above every stored atom id, so the dense map declines it.
  chill::getCorrel(dense, {{0, 3, -2, 99}, {1, 3}, {2, 1}, {3, 4, 1}}, false, 0);
  REQUIRE(dense.pts[2].c_ij.size() == 1);
  REQUIRE(dense.pts[0].c_ij.size() >= 1);

  unsetenv("DSEAMS_PHASE");
  Cloud phase = makeCloud({{0, 0, 0}, {1.5, 0, 0}}, {10, 10, 10});
  chill::getCorrel(phase, {{0, 1}, {1, 0}}, false, 1);
  REQUIRE(phase.pts[0].c_ij.size() == 1);
}

TEST_CASE("a non-hex ring and a negative atom reach the basal candidate tails",
          "[cover_ice][cage]") {
  Cloud cloud = hexCage();
  auto byId = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, byId);
  auto rings = primitive::ringNetwork(rows, 7);
  auto six = rings;
  six.erase(std::remove_if(six.begin(), six.end(),
                           [](const std::vector<int> &ring) {
                             return ring.size() != 6;
                           }),
            six.end());
  REQUIRE(six.size() >= 5);

  auto affiliation = ring::cageAffiliation(six, rows);
  REQUIRE(affiliation.hc.size() == six.size());
  REQUIRE(affiliation.ddc.size() == six.size());
  int hcRings = 0;
  for (std::size_t i = 0; i < six.size(); i++) {
    hcRings += affiliation.hc[i] ? 1 : 0;
  }
  REQUIRE(hcRings >= 5);

  auto planes = ring::stackingPlanes(six, rows);
  REQUIRE(planes.basal.size() == six.size());
  REQUIRE(planes.equatorial.size() == six.size());

  auto none = ring::cageAffiliation({}, rows);
  REQUIRE(none.hc.empty());
  auto noPlanes = ring::stackingPlanes({}, rows);
  REQUIRE(noPlanes.basal.empty());

  // One hexagon plus a shorter ring that shares its atoms. The short ring
  // makes the all-hex search unavailable, and the shared atoms stop the pair
  // before findPrismatic indexes the short ring.
  std::vector<std::vector<int>> mixed = {{0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4}};
  std::vector<std::vector<int>> mixedRows = selfRows(6);
  for (int i = 0; i < 6; i++) {
    link(mixedRows, i, (i + 1) % 6);
  }
  mixedRows[0].push_back(-1);
  mixedRows[2].push_back(20);
  auto mixedAff = ring::cageAffiliation(mixed, mixedRows);
  REQUIRE(mixedAff.hc.size() == 2);
  REQUIRE_FALSE(mixedAff.hc[0]);
  REQUIRE_FALSE(mixedAff.hc[1]);

  std::vector<std::vector<int>> ragged = {
      {0, 1, -4, 8}, {1, 0, 2}, {2, 1, 3}, {3, 2, 4}, {4, 3, 5}, {5, 4, 0}};
  auto taggedAff = ring::cageAffiliation({{-1, 0, 1, 2, 3, 4}}, ragged);
  REQUIRE(taggedAff.hc.size() == 1);
  REQUIRE_FALSE(taggedAff.hc[0]);

  // A neighbour row that names the atom itself, and two copies of one hexagon.
  std::vector<std::vector<int>> looped = selfRows(6);
  for (int i = 0; i < 6; i++) {
    link(looped, i, (i + 1) % 6);
    looped[static_cast<std::size_t>(i)].push_back(i);
  }
  auto copies = ring::cageAffiliation({{0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}}, looped);
  REQUIRE(copies.hc.size() == 2);
  REQUIRE_FALSE(copies.hc[0]);

  std::vector<std::vector<int>> many(1024, {0, 1, 2, 3, 4, 5});
  auto wide = ring::cageAffiliation(many, looped);
  REQUIRE(wide.hc.size() == 1024);
  ring::AffiliationUpdater bulk;
  const auto &bulkFlags = bulk.update(many, looped);
  REQUIRE(bulkFlags.hc.size() == 1024);
  REQUIRE(bulk.lastReclassified() == 1024);
}

TEST_CASE("affiliation update carries a reversed cycle and falls back to a batch",
          "[cover_ice][cage]") {
  std::vector<std::vector<int>> rows = {
      {0, 1, 5}, {1, 0, 2}, {2, 1, 3}, {3, 2, 4}, {4, 3, 5}, {5, 4, 0, -3, 9}};
  std::vector<std::vector<int>> backward = {{0, 5, 4, 3, 2, 1},
                                            {0, 1, 2, 3, 4, 5}};
  ring::AffiliationUpdater updater;
  const auto &first = updater.update(backward, rows);
  REQUIRE(first.hc.size() == 2);
  REQUIRE(updater.lastReclassified() == 2);

  std::vector<std::vector<int>> forward = {{0, 1, 2, 3, 4, 5},
                                           {0, 2, 1, 3, 4, 6}};
  const auto &second = updater.update(forward, rows);
  REQUIRE(second.hc.size() == 2);

  updater.setBatchFraction(0.0);
  rows[0].push_back(2);
  const auto &hot = updater.update(forward, rows);
  REQUIRE(hot.hc.size() == 2);
  rows[1].push_back(3);
  const auto &again = updater.update(forward, rows);
  REQUIRE(again.hc.size() == 2);
  REQUIRE(updater.lastReclassified() == 2);

  const auto &cleared = updater.update({}, rows);
  REQUIRE(cleared.hc.empty());
  REQUIRE(updater.lastReclassified() == 0);
}

TEST_CASE("seeded labels reject an out-of-range neighbour and fill one vacancy",
          "[cover_ice][cage]") {
  std::vector<std::vector<int>> strict = {{0, 1, 2, 3, 4, 5}};
  std::vector<std::vector<int>> permissive = {{0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, -2}};
  auto strictRows = selfRows(6);
  link(strictRows, 0, 1);
  auto permissiveRows = strictRows;
  permissiveRows[0].push_back(-1);
  permissiveRows[0].push_back(8);
  permissiveRows[1].push_back(0);

  ring::CageAffiliation strictAff;
  strictAff.hc = {true};
  strictAff.ddc = {true};
  ring::CageAffiliation permissiveAff;
  permissiveAff.hc = {true, false};
  permissiveAff.ddc = {true, false};

  auto seeded = ring::seededCageAffiliation(strict, strictAff, permissive,
                                            permissiveAff, permissiveRows, false);
  REQUIRE(seeded.hc.size() == 6);
  REQUIRE(seeded.hc[0]);
  REQUIRE(seeded.ddc[0]);

  auto completed = ring::seededCageAffiliation(strict, strictAff, permissive,
                                               permissiveAff, permissiveRows, true);
  REQUIRE(completed.hc.size() == 6);

  ring::SeededAtomLabels labels;
  labels.hc = {true, true, true, true, true, false};
  labels.ddc.assign(6, false);
  auto filled = ring::ringAdjacentCompletion(labels, {{0, 1, 2, 3, 4, 5}});
  REQUIRE(filled.hc[5]);
  REQUIRE_FALSE(filled.ddc[5]);

  auto emptyLabels = ring::ringAdjacentCompletion({}, {});
  REQUIRE(emptyLabels.hc.empty());
  auto noRings = ring::ringAdjacentCompletion(labels, {});
  REQUIRE(noRings.hc == labels.hc);
}

TEST_CASE("NeighBond and basal conditions reject a negative or partial bond",
          "[cover_ice][topo]") {
  std::vector<std::vector<int>> rows = {{0, 1, -1, 7}, {1, 0}, {2}};
  ring::NeighBond bonds;
  bonds.build(rows);
  REQUIRE(bonds(0, 1));
  REQUIRE_FALSE(bonds(0, -1));
  REQUIRE_FALSE(bonds(-1, 0));
  REQUIRE_FALSE(bonds(0, 9));
  REQUIRE_FALSE(bonds(4, 0));

  std::vector<int> basal1 = {0, 1, 2, 3, 4, 5};
  std::vector<int> basal2 = {-1, 1, 2, 3, 4, 5};
  REQUIRE_FALSE(ring::basalConditions(rows, basal1, basal2));
  REQUIRE_FALSE(ring::basalConditions(bonds, basal1, basal2));

  std::vector<ring::strucType> ringType = {ring::strucType::bothBasal,
                                           ring::strucType::unclassified};
  std::vector<int> listDDC = {0, 2};
  std::vector<int> listHC = {0};
  std::vector<std::vector<int>> rings = {basal1, basal2};
  auto mixed = ring::findMixedRings(rings, ringType, listDDC, listHC);
  REQUIRE(mixed.size() == 1);
  REQUIRE(listDDC[0] == -10);
  REQUIRE(listHC[0] == -10);

  Cloud twins = makeCloud({{1, 1, 1}, {11, 1, 1}}, {10, 10, 10});
  REQUIRE_THROWS_AS(ring::positionRank(twins), std::invalid_argument);
  Cloud shortBox = makeCloud({{1, 2, 3}, {1.2, 2.2, 3.2}}, {0, 4});
  auto ranks = ring::positionRank(shortBox);
  REQUIRE(ranks.size() == 2);
  REQUIRE(ranks[0] != ranks[1]);

  Cloud thin = makeCloud({{0, 0, 0}, {1, 0, 0}}, {1.0, 8, 8});
  REQUIRE_THROWS_AS(ring::cageRings({{0, 1}}, thin, 3.0), std::invalid_argument);
  auto kept = ring::contractibleRings({{-1, 0}, {0, 1}}, thin);
  REQUIRE(kept.size() <= 2);

  auto emptyIndex = ring::buildRingSearchIndex(rings, 0);
  REQUIRE(emptyIndex.ringOff.empty());
  std::vector<int> wide(17);
  std::iota(wide.begin(), wide.end(), 0);
  std::vector<std::vector<int>> oddRings = {wide, {}};
  std::vector<ring::strucType> oddType(2, ring::strucType::unclassified);
  std::vector<int> listHCEmpty;
  std::vector<cage::Cage> cages;
  auto ddc = ring::findDDC(oddRings, oddType, listHCEmpty, cages);
  REQUIRE(ddc.empty());
}

TEST_CASE("a wide equator with unpaired triplets stops before a cage is kept",
          "[cover_ice][topo]") {
  std::vector<std::vector<int>> rings;
  rings.push_back({0, 1, 2, 3, 4, 5});
  for (int copy = 0; copy < 7; copy++) {
    const int extra = 100 + copy * 30;
    rings.push_back({0, 1, 2, extra, extra + 1, extra + 2});
    rings.push_back({2, 3, 4, extra + 3, extra + 4, extra + 5});
    rings.push_back({4, 5, 0, extra + 6, extra + 7, extra + 8});
    rings.push_back({1, 2, 3, extra + 9, extra + 10, extra + 11});
    rings.push_back({3, 4, 5, extra + 12, extra + 13, extra + 14});
    rings.push_back({5, 0, 1, extra + 15, extra + 16, extra + 17});
  }
  rings.push_back({0, 1, 2, 3, 4});
  rings.push_back({0, 1, 2, -1, 8, 9});
  const auto index = ring::buildRingSearchIndex(rings, 4000);
  std::vector<int> peripheral;
  REQUIRE_FALSE(ring::equatorialConditions(rings, 0, index, &peripheral));
  REQUIRE(peripheral.empty());

  std::vector<std::vector<int>> paired;
  paired.push_back({0, 1, 2, 3, 4, 5});
  for (int copy = 0; copy < 7; copy++) {
    paired.push_back({0, 1, 2, 10, 11, 12});
    paired.push_back({2, 3, 4, 10, 11, 12});
    paired.push_back({4, 5, 0, 10, 11, 12});
    const int extra = 200 + copy * 3;
    paired.push_back({1, 2, 3, extra, extra + 1, extra + 2});
    paired.push_back({3, 4, 5, extra + 30, extra + 31, extra + 32});
    paired.push_back({5, 0, 1, extra + 60, extra + 61, extra + 62});
  }
  const auto pairedIndex = ring::buildRingSearchIndex(paired, 4000);
  std::vector<int> chosen;
  std::vector<int> rank(paired.size());
  std::iota(rank.begin(), rank.end(), 0);
  ring::EquatorChoice raw;
  REQUIRE_FALSE(ring::equatorialConditions(paired, 0, pairedIndex, &peripheral,
                                           &chosen, &rank, &raw));
  REQUIRE(chosen.empty());
  REQUIRE(raw.n == 0);
}

TEST_CASE("prism blocks cover a missing bond, the other axes, and a lost anchor",
          "[cover_ice][prism]") {
  Cloud cloud = squarePrism(2.75, {10, 10, 50});
  auto rows = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
  unlink(rows, 0, 4);
  auto rings = primitive::ringNetwork(rows, 5);
  std::vector<std::vector<int>> squares;
  for (const auto &ring : rings) {
    if (ring.size() == 4) {
      squares.push_back(ring);
    }
  }
  REQUIRE(squares.size() >= 2);
  std::vector<ring::strucType> ringType(squares.size(), ring::strucType::unclassified);
  std::vector<double> rmsd(cloud.nop, -1.0);
  int perfect = 0;
  int deformed = 0;
  auto found = ring::findPrisms(squares, ringType, perfect, deformed, rows, cloud,
                                rmsd, true);
  REQUIRE(perfect + deformed >= 1);
  REQUIRE_FALSE(found.empty());
  bool sawDeformed = false;
  for (auto kind : ringType) {
    sawDeformed = sawDeformed || kind == ring::strucType::deformedPrism ||
                  kind == ring::strucType::mixedPrismRing;
  }
  REQUIRE(sawDeformed);

  std::vector<ring::strucType> strictType(squares.size(), ring::strucType::unclassified);
  int strictPerfect = 0;
  int strictDeformed = 0;
  auto strict = ring::findPrisms(squares, strictType, strictPerfect, strictDeformed,
                                 rows, cloud, rmsd, false);
  REQUIRE(strictDeformed == 0);
  REQUIRE(static_cast<int>(strict.size()) == strictPerfect);

  std::vector<int> one = {0, 1, 2, 3};
  std::vector<int> two = {4, 5, 6, 7};
  REQUIRE_FALSE(ring::basalPrismConditions(rows, one, two));
  REQUIRE(ring::relaxedPrismConditions(rows, one, two));
  std::vector<int> rotated = {1, 2, 3, 0};
  REQUIRE_FALSE(ring::basalPrismConditions(rows, rotated, two));
  std::vector<std::vector<int>> lonely = selfRows(8);
  REQUIRE_FALSE(ring::relaxedPrismConditions(lonely, one, two));
  REQUIRE_FALSE(ring::basalPrismConditions(lonely, one, two));

  // A strict pair first marks both rings Prism. Dropping one bond on the
  // next call takes the deformed branch while those rings are already Prism.
  Cloud stacked = squarePrism(2.75, {10, 10, 50});
  auto full = nneigh::neighbourListByIndex(stacked, nneigh::neighListO(3.5, stacked, 1));
  std::vector<std::vector<int>> ordered = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<ring::strucType> staged(2, ring::strucType::unclassified);
  std::vector<double> stagedRmsd(stacked.nop, -1.0);
  int nPerfect = 0;
  int nDeformed = 0;
  auto strictIds = ring::findPrisms(ordered, staged, nPerfect, nDeformed, full,
                                    stacked, stagedRmsd, false);
  REQUIRE(nPerfect == 1);
  REQUIRE(staged[0] == ring::strucType::Prism);
  REQUIRE(staged[1] == ring::strucType::Prism);
  auto dropped = full;
  unlink(dropped, 0, 4);
  nPerfect = 0;
  nDeformed = 0;
  auto deformedIds = ring::findPrisms(ordered, staged, nPerfect, nDeformed, dropped,
                                      stacked, stagedRmsd, true);
  REQUIRE(nDeformed == 1);
  REQUIRE(staged[0] == ring::strucType::mixedPrismRing);
  REQUIRE(staged[1] == ring::strucType::mixedPrismRing);
  REQUIRE_FALSE(strictIds.empty());
  // Already-Prism rings are retyped in place and are not appended again.
  REQUIRE(deformedIds.empty());

  // The other way round: a deformed pair, then the restored bond, so the
  // strict branch sees rings that are already deformed.
  std::vector<ring::strucType> reverse(2, ring::strucType::unclassified);
  nPerfect = 0;
  nDeformed = 0;
  REQUIRE_FALSE(ring::findPrisms(ordered, reverse, nPerfect, nDeformed, dropped,
                                 stacked, stagedRmsd, true)
                    .empty());
  REQUIRE(reverse[0] == ring::strucType::deformedPrism);
  nPerfect = 0;
  nDeformed = 0;
  auto restored = ring::findPrisms(ordered, reverse, nPerfect, nDeformed, full,
                                   stacked, stagedRmsd, false);
  REQUIRE(restored.empty());
  REQUIRE(nPerfect == 1);
  REQUIRE(reverse[0] == ring::strucType::mixedPrismRing);
  REQUIRE(reverse[1] == ring::strucType::mixedPrismRing);

  std::vector<std::vector<int>> apart = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<ring::strucType> apartType(2, ring::strucType::unclassified);
  nPerfect = 0;
  nDeformed = 0;
  auto none = ring::findPrisms(apart, apartType, nPerfect, nDeformed, lonely,
                               stacked, stagedRmsd, true);
  REQUIRE(none.empty());
  REQUIRE(nPerfect == 0);
  REQUIRE(nDeformed == 0);

  // One cross bond keeps the relaxed test, and a coplanar second ring fails
  // the prism shape match.
  Cloud flat = makeCloud({{0, 0, 0},
                          {2.75, 0, 0},
                          {2.75, 2.75, 0},
                          {0, 2.75, 0},
                          {6, 0, 0},
                          {8.75, 0, 0},
                          {8.75, 2.75, 0},
                          {6, 2.75, 0}},
                         {10, 10, 50});
  auto flatRows = selfRows(8);
  for (int i = 0; i < 4; i++) {
    link(flatRows, i, (i + 1) % 4);
    link(flatRows, 4 + i, 4 + (i + 1) % 4);
  }
  link(flatRows, 1, 5);
  std::vector<std::vector<int>> flatRings = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<ring::strucType> flatType(2, ring::strucType::unclassified);
  std::vector<double> flatRmsd(8, -1.0);
  nPerfect = 0;
  nDeformed = 0;
  auto unmatched = ring::findPrisms(flatRings, flatType, nPerfect, nDeformed, flatRows,
                                    flat, flatRmsd, true);
  REQUIRE(nPerfect == 0);
  REQUIRE(nDeformed == 1);
  REQUIRE(unmatched.size() == 2);

  std::vector<int> yzA = {0, 1, 2, 3};
  std::vector<int> yzB = {4, 5, 6, 7};
  Cloud alongX = makeCloud({{0, 0, 0},
                            {0, 2, 0},
                            {0, 2, 2},
                            {0, 0, 2},
                            {2, 0, 0},
                            {2, 2, 0},
                            {2, 2, 2},
                            {2, 0, 2}},
                           {50, 10, 10});
  REQUIRE(ring::discardExtraTetragonBlocks(yzA, yzB, alongX));

  Cloud alongY = makeCloud({{0, 0, 0},
                            {2, 0, 0},
                            {2, 0, 2},
                            {0, 0, 2},
                            {0, 2, 0},
                            {2, 2, 0},
                            {2, 2, 2},
                            {0, 2, 2}},
                           {10, 50, 10});
  REQUIRE(ring::discardExtraTetragonBlocks(yzA, yzB, alongY));
  std::vector<int> lateralA = {0, 1, 5, 4};
  std::vector<int> lateralB = {3, 2, 6, 7};
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lateralA, lateralB, cloud));

  std::vector<std::vector<int>> typed = {{0, 1, 2, 3}, {0, 1, 2, 3}};
  std::vector<ring::strucType> kinds = {ring::strucType::deformedPrism,
                                        ring::strucType::Prism};
  std::vector<int> atomTypes = {1, 1, 1, 1};
  std::vector<ring::strucType> atomState(4, ring::strucType::unclassified);
  REQUIRE(ring::assignPrismType(typed, {0, 1}, 4, kinds, atomTypes, atomState) == 0);
  REQUIRE(atomState[0] == ring::strucType::mixedPrismRing);
  REQUIRE(atomTypes[0] == 4);

  kinds[0] = ring::strucType::mixedPrismRing;
  atomState.assign(4, ring::strucType::Prism);
  REQUIRE(ring::assignPrismType(typed, {0}, 4, kinds, atomTypes, atomState) == 0);
  REQUIRE(atomState[0] == ring::strucType::mixedPrismRing);
  kinds[0] = ring::strucType::deformedPrism;
  atomState.assign(4, ring::strucType::Prism);
  REQUIRE(ring::assignPrismType(typed, {0}, 4, kinds, atomTypes, atomState) == 0);
  REQUIRE(atomState[0] == ring::strucType::mixedPrismRing);

  std::vector<ring::strucType> states = {ring::strucType::deformedPrism,
                                         ring::strucType::mixedPrismRing,
                                         ring::strucType::Prism};
  std::vector<int> deformedTypes = {4, 4, 3};
  REQUIRE(ring::deformedPrismTypes(states, deformedTypes, 6) == 0);
  REQUIRE(deformedTypes[0] == 8);
  REQUIRE(deformedTypes[1] == 2);
  REQUIRE(deformedTypes[2] == 3);

  int anchor = 0;
  cloud.pts[0].z = 3.0;
  REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 1) == 0);
  REQUIRE(anchor == 0);
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(0.0, 1e-8));
  cloud.pts[0].z = 1.5;
  REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 2) == 0);
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(0.0, 1e-8));
  anchor = 99;
  cloud.idIndexMap.erase(0);
  REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 2) == 1);
}

TEST_CASE("bulk matching records every ice type and a cage that cannot be ordered",
          "[cover_ice][tum]") {
  Cloud hc = hexCage();
  auto hcRows = nneigh::neighbourListByIndex(hc, nneigh::neighListO(3.5, hc, 1));
  auto hcRings = primitive::ringNetwork(hcRows, 7);

  const double radius = 1.2;
  std::vector<std::array<double, 3>> pentagon;
  for (int i = 0; i < 5; i++) {
    const double angle = 2.0 * M_PI * static_cast<double>(i) / 5.0;
    pentagon.push_back({20.0 + radius * std::cos(angle),
                        20.0 + radius * std::sin(angle), 0.0});
  }
  for (int i = 0; i < 5; i++) {
    auto top = pentagon[static_cast<std::size_t>(i)];
    top[2] = 1.2;
    pentagon.push_back(top);
  }
  pentagon.push_back({40.0, 40.0, 40.0});
  for (const auto &point : pentagon) {
    molSys::Point<double> extra;
    extra.type = 1;
    extra.atomID = hc.nop;
    extra.x = point[0];
    extra.y = point[1];
    extra.z = point[2];
    hc.idIndexMap[hc.nop] = hc.nop;
    hc.pts.push_back(extra);
    hc.nop++;
  }
  hcRows.resize(static_cast<std::size_t>(hc.nop));
  for (int i = 12; i < hc.nop; i++) {
    hcRows[static_cast<std::size_t>(i)] = {i};
  }
  for (int i = 0; i < 5; i++) {
    link(hcRows, 12 + i, 12 + (i + 1) % 5);
    link(hcRows, 17 + i, 17 + (i + 1) % 5);
    link(hcRows, 12 + i, 17 + i);
  }
  hcRings.push_back({12, 13, 14, 15, 16});
  hcRings.push_back({17, 18, 19, 20, 21});

  auto dir = scratch("dseams_cover_ice_tum");
  int rc = tum3::topoUnitMatchingBulk(dir.string() + "/", hcRings, hcRows, hc, 1,
                                      true, false, "../templates");
  REQUIRE(rc == 0);
  REQUIRE(fs::exists(dir / "bulkTopo"));
  fs::remove_all(dir);

  std::vector<ring::strucType> ringType;
  int numHC = 0;
  int numDDC = 0;
  auto emptyCages = tum3::topoBulkCriteria(dir.string(), {}, hcRows, hc, 1, numHC,
                                           numDDC, ringType);
  REQUIRE(emptyCages.empty());
  REQUIRE(numHC == 0);
  REQUIRE(numDDC == 0);

  std::vector<ring::strucType> matchedType(hcRings.size(),
                                           ring::strucType::unclassified);
  std::vector<cage::Cage> matched;
  auto hcIds = ring::findHC(hcRings, matchedType, hcRows, matched);
  REQUIRE_FALSE(hcIds.empty());
  REQUIRE_FALSE(matched.empty());
  Eigen::MatrixXd badRef(11, 3);
  badRef.setZero();
  std::vector<double> quat;
  double rmsd = 5.0;
  REQUIRE(tum3::shapeMatchHC(hc, badRef, matched[0], hcRings, hcRows, quat, rmsd) ==
          0);
  REQUIRE(rmsd < 0.0);

  cage::Cage loose;
  loose.type = cage::cageType::DoubleDiaC;
  loose.rings = {0};
  std::vector<std::vector<int>> two = {{0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}};
  rmsd = 5.0;
  Cloud twelve = makeCloud({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
                            {2, 0, 0}, {2, 1, 0}, {2, 0, 1}, {2, 1, 1}},
                           {10, 10, 10});
  REQUIRE(tum3::shapeMatchDDC(twelve, badRef, {loose}, 0, two, quat, rmsd) == 0);
  REQUIRE(rmsd < 0.0);

  std::vector<double> perAtom = {-1, -1, -1, -1, -1, -1};
  std::vector<int> common = {0, 0, 0, 0, 0, 0};
  std::vector<cage::iceType> kinds = {
      cage::iceType::pnc,     cage::iceType::mixed2, cage::iceType::hc,
      cage::iceType::dummy,   cage::iceType::ddc,    cage::iceType::mixed};
  REQUIRE(tum3::updateRMSDatom(two, loose, 0.25, perAtom, common, kinds) == 0);
  REQUIRE(perAtom[0] == -1.0);
  REQUIRE(perAtom[1] == -1.0);
  REQUIRE_THAT(perAtom[2], Catch::Matchers::WithinAbs(0.25, 1e-12));
  REQUIRE(common[2] == 1);
  REQUIRE(tum3::updateRMSDatom(two, loose, -1.0, perAtom, common, kinds) == 0);
  REQUIRE(common[2] == 1);
  REQUIRE(tum3::averageRMSDatom(perAtom, common) == 0);
  REQUIRE(perAtom[0] == -1.0);
  REQUIRE(perAtom[1] == -1.0);
  REQUIRE_THAT(perAtom[2], Catch::Matchers::WithinAbs(0.25, 1e-12));
  REQUIRE_THAT(perAtom[3], Catch::Matchers::WithinAbs(0.25, 1e-12));

  auto clusterDir = scratch("dseams_cover_ice_cluster");
  Cloud clusterCloud = twelve;
  std::vector<cage::Cage> clusterCages(2);
  clusterCages[0].type = cage::cageType::HexC;
  clusterCages[0].rings = {0};
  clusterCages[1].type = cage::cageType::DoubleDiaC;
  clusterCages[1].rings = {1};
  REQUIRE(tum3::clusterCages(clusterCloud, clusterDir.string() + "/", two,
                             clusterCages, 1, 1) == 0);
  const auto info = clusterDir / "bulkTopo" / "clusterXYZ" /
                    ("frame-" + std::to_string(clusterCloud.currentFrame)) / "info.dat";
  REQUIRE(fs::exists(info));
  fs::remove_all(clusterDir);
}

TEST_CASE("cage counts cover an empty graph, a forced offload, and a long ring row",
          "[cover_ice][tum]") {
  unsetenv("SEAMS_OFFLOAD");
  REQUIRE_FALSE(tum::preferOffload());
  auto absent = tum::hostCageCounts({});
  REQUIRE(absent.nSix == 0);
  REQUIRE_FALSE(absent.usedDevice);
  auto viaHost = tum::cageCounts({});
  REQUIRE(viaHost.nSix == 0);

  setenv("SEAMS_OFFLOAD", "", 1);
  REQUIRE_FALSE(tum::preferOffload());
  setenv("SEAMS_OFFLOAD", "0", 1);
  REQUIRE_FALSE(tum::preferOffload());
  setenv("SEAMS_OFFLOAD", "1", 1);
  REQUIRE(tum::preferOffload());

  std::vector<std::vector<int>> hex = selfRows(6);
  for (int i = 0; i < 6; i++) {
    link(hex, i, (i + 1) % 6);
  }
  auto counted = tum::cageCounts(hex);
  REQUIRE(counted.nSix >= 1);
  REQUIRE(counted.atomSix.size() == 6);
  REQUIRE(counted.atomSix[0] >= 1);

  auto emptySpecial = tum::specializedCageCounts({});
  REQUIRE(emptySpecial.nSix == 0);
  auto blankRow = tum::hostCageCounts({{0, 1}, {}});
  REQUIRE(blankRow.atomSix.size() == 2);
  auto blankSpecial = tum::specializedCageCounts({{}, {1, 0}});
  REQUIRE(blankSpecial.nSix == 0);
  tum::SpecializedCaps tiny;
  tiny.ringCap = 1;
  tiny.pairCap = 1;
  tiny.rowCap = 1;
  auto grown = tum::specializedCageCounts(hex, tiny);
  REQUIRE(grown.nSix >= 1);

  std::vector<std::vector<int>> many = selfRows(1 + 33 * 5);
  for (int copy = 0; copy < 33; copy++) {
    const int base = 1 + copy * 5;
    link(many, 0, base);
    for (int step = 0; step < 4; step++) {
      link(many, base + step, base + step + 1);
    }
    link(many, base + 4, 0);
  }
  auto overflow = tum::specializedCageCounts(many);
  REQUIRE(overflow.nSix >= 33);
  REQUIRE(overflow.atomSix[0] >= 33);
  unsetenv("SEAMS_OFFLOAD");
}

TEST_CASE("Voronoi weights certify a tetrahedron and reject a bad degree",
          "[cover_ice][voronoi]") {
  Cloud emptyBox = makeCloud({{0, 0, 0}, {1, 0, 0}}, {4, 4});
  auto rejected = chill::voronoiFacetWeights(emptyBox, 2.0);
  REQUIRE(rejected.size() == 2);
  REQUIRE(rejected[0].neighbours.empty());
  REQUIRE_FALSE(rejected[0].certified);

  Cloud pair = makeCloud({{0, 0, 0}, {1, 0, 0}}, {20, 20, 20});
  auto noCutoff = chill::voronoiFacetWeights(pair, 0.0);
  REQUIRE(noCutoff[0].neighbours.empty());
  auto coincident = makeCloud({{0, 0, 0}, {0, 0, 0}, {1.2, 0, 0}}, {20, 20, 20});
  auto skipped = chill::voronoiFacetWeights(coincident, 3.0);
  REQUIRE(skipped.size() == 3);

  Cloud openCell = makeCloud({{0, 0, 0},
                              {1, 0, 0},
                              {0.5, std::sqrt(3.0) / 2.0, 0},
                              {0.5, std::sqrt(3.0) / 6.0, std::sqrt(6.0) / 3.0}},
                             {30, 30, 30});
  auto openWeights = chill::voronoiFacetWeights(openCell, 2.0);
  REQUIRE_FALSE(openWeights[0].neighbours.empty());
  REQUIRE_FALSE(openWeights[0].certified);

  Cloud cube = makeCloud({{5, 5, 5},
                          {7, 5, 5},
                          {3, 5, 5},
                          {5, 7, 5},
                          {5, 3, 5},
                          {5, 5, 7},
                          {5, 5, 3}},
                         {20, 20, 20});
  auto weights = chill::voronoiFacetWeights(cube, 6.0);
  REQUIRE(weights[0].neighbours.size() == 6);
  double sum = 0.0;
  for (double weight : weights[0].weights) {
    sum += weight;
  }
  REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-8));
  REQUIRE(weights[0].certified);

  auto ql = chill::steinhardtQlVoronoi(cube, weights, 4);
  REQUIRE(ql.ql.size() == 7);
  REQUIRE(ql.ql[0] > 0.0);
  REQUIRE(ql.qlBar[0] > 0.0);
  auto q3 = chill::steinhardtQlVoronoi(cube, weights, 3);
  REQUIRE(q3.ql[0] > 0.0);
  Cloud alone = makeCloud({{4, 4, 4}}, {20, 20, 20});
  auto lonely = chill::voronoiFacetWeights(alone, 3.0);
  REQUIRE(lonely[0].neighbours.empty());
  REQUIRE_FALSE(lonely[0].certified);
  auto wrong = chill::steinhardtQlVoronoi(cube, weights, 5);
  REQUIRE_THAT(wrong.ql[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THAT(wrong.qlBar[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  auto fromCloud = chill::steinhardtQlVoronoi(cube, 6.0, 6);
  REQUIRE(fromCloud.ql.size() == 7);
  REQUIRE(fromCloud.ql[0] > 0.0);
}
