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

#include <franzblau.hpp>
#include <neighbours.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <unistd.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace fs = std::filesystem;

fs::path pidTemp(const char *name) {
  return fs::temp_directory_path() /
         (std::string(name) + "-" +
          std::to_string(static_cast<long long>(::getpid())));
}


namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud makeCloud(const std::vector<std::array<double, 3>> &xyz,
                const std::vector<double> &box) {
  Cloud cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int i = 0; i < static_cast<int>(xyz.size()); ++i) {
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

template <std::size_t N>
Cloud makeCloud(const double (&xyz)[N][3], const std::vector<double> &box) {
  std::vector<std::array<double, 3>> points;
  points.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    points.push_back({xyz[i][0], xyz[i][1], xyz[i][2]});
  }
  return makeCloud(points, box);
}

std::vector<std::vector<int>> selfRows(int n) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
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

Cloud prismCloud(double height, const std::vector<double> &box) {
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

std::vector<std::vector<int>> prismRows() {
  auto rows = selfRows(8);
  for (int i = 0; i < 4; ++i) {
    link(rows, i, (i + 1) % 4);
    link(rows, 4 + i, 4 + ((i + 1) % 4));
    link(rows, i, 4 + i);
  }
  return rows;
}

std::string capture(const std::function<void()> &body) {
  std::stringstream buffer;
  auto *previous = std::cerr.rdbuf(buffer.rdbuf());
  body();
  std::cerr.rdbuf(previous);
  return buffer.str();
}

const std::array<std::array<double, 3>, 4> kYz = {
    {{0, 0, 0}, {0, 2, 0}, {0, 2, 2}, {0, 0, 2}}};
const std::array<std::array<double, 3>, 4> kXz = {
    {{0, 0, 0}, {2, 0, 0}, {2, 0, 2}, {0, 0, 2}}};
const std::array<std::array<double, 3>, 4> kXy = {
    {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}};
// areaYZ > areaXY, and areaYZ is not greater than areaXZ.
const std::array<std::array<double, 3>, 4> kXPartial = {
    {{0, 0, 0}, {2, 0.4, 0}, {2, 0.2, 2}, {0, 0, 2}}};
// areaXZ > areaXY, and areaXZ is not greater than areaYZ.
const std::array<std::array<double, 3>, 4> kYPartial = {
    {{0, 0, 0}, {1, 2, 0}, {1, 2, 3}, {0, 0, 3}}};
// areaXY > areaXZ, and areaXY is not greater than areaYZ.
const std::array<std::array<double, 3>, 4> kZPartial = {
    {{0, 0, 0}, {1, 0, -2}, {0, 5, -2.5}, {-1, 5, -0.5}}};

Cloud ringCloud(const std::array<std::array<double, 3>, 4> &first,
                const std::array<std::array<double, 3>, 4> &second,
                const std::vector<double> &box) {
  std::vector<std::array<double, 3>> xyz;
  for (const auto &point : first) {
    xyz.push_back(point);
  }
  for (const auto &point : second) {
    xyz.push_back({point[0] + 6.0, point[1] + 6.0, point[2]});
  }
  return makeCloud(xyz, box);
}

} // namespace

TEST_CASE("prism blocks record a perfect cage, a rejected match, and a repeat",
          "[cover_topo][prism]") {
  Cloud cloud = prismCloud(2.75, {10, 10, 50});
  auto rows = prismRows();
  std::vector<std::vector<int>> ordered = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<ring::strucType> strict(2, ring::strucType::unclassified);
  std::vector<double> rmsd(cloud.nop, -1.0);
  int perfect = 0;
  int deformed = 0;
  auto first = ring::findPrisms(ordered, strict, perfect, deformed, rows, cloud,
                                rmsd, false);
  REQUIRE(perfect == 1);
  REQUIRE(deformed == 0);
  REQUIRE(first.size() == 2);
  REQUIRE(strict[0] == ring::strucType::Prism);
  REQUIRE(strict[1] == ring::strucType::Prism);
  perfect = 0;
  deformed = 0;
  auto again = ring::findPrisms(ordered, strict, perfect, deformed, rows, cloud,
                                rmsd, false);
  REQUIRE(again.empty());
  REQUIRE(perfect == 1);
  REQUIRE(strict[0] == ring::strucType::Prism);
  REQUIRE(strict[1] == ring::strucType::Prism);

  auto dropped = rows;
  unlink(dropped, 0, 4);
  std::vector<ring::strucType> bent(2, ring::strucType::unclassified);
  perfect = 0;
  deformed = 0;
  auto matched = ring::findPrisms(ordered, bent, perfect, deformed, dropped, cloud,
                                  rmsd, true);
  REQUIRE(perfect == 0);
  REQUIRE(deformed == 1);
  REQUIRE(matched.size() == 2);
  REQUIRE(bent[0] == ring::strucType::deformedPrism);
  REQUIRE(bent[1] == ring::strucType::deformedPrism);
  perfect = 0;
  deformed = 0;
  auto matchedAgain = ring::findPrisms(ordered, bent, perfect, deformed, dropped,
                                       cloud, rmsd, true);
  REQUIRE(matchedAgain.empty());
  REQUIRE(deformed == 1);
  REQUIRE(bent[0] == ring::strucType::deformedPrism);
  REQUIRE(bent[1] == ring::strucType::deformedPrism);

  std::vector<std::array<double, 3>> star;
  for (int i = 0; i < 5; ++i) {
    const double angle = 2.0 * 3.141592653589793 * static_cast<double>(i) / 5.0;
    star.push_back({3.0 * std::cos(angle), 3.0 * std::sin(angle), 0.0});
  }
  star.push_back({8, 0, 0.4});
  star.push_back({8, 0, 0.4});
  star.push_back({8, 0, 0.4});
  star.push_back({8, 0, 0.4});
  star.push_back({14, 6, 0.4});
  Cloud flat = makeCloud(star, {40, 20, 80});
  auto flatRows = selfRows(10);
  for (int i = 0; i < 5; ++i) {
    link(flatRows, i, (i + 1) % 5);
    link(flatRows, 5 + i, 5 + ((i + 1) % 5));
  }
  link(flatRows, 0, 5);
  std::vector<std::vector<int>> pentagons = {{0, 1, 2, 3, 4}, {5, 6, 7, 8, 9}};
  std::vector<ring::strucType> flatType(2, ring::strucType::unclassified);
  std::vector<double> flatRmsd(10, -1.0);
  perfect = 0;
  deformed = 0;
  auto rejected = ring::findPrisms(pentagons, flatType, perfect, deformed, flatRows,
                                   flat, flatRmsd, true);
  CAPTURE(perfect, deformed, rejected.size());
  REQUIRE(rejected.empty());
  REQUIRE(perfect == 0);
  REQUIRE(deformed == 0);
  REQUIRE(flatType[0] == ring::strucType::unclassified);
  REQUIRE(flatType[1] == ring::strucType::unclassified);

  perfect = 0;
  deformed = 0;
  auto skipped = ring::findPrisms(pentagons, flatType, perfect, deformed, flatRows,
                                  flat, flatRmsd, false);
  REQUIRE(skipped.empty());
  REQUIRE(perfect == 0);
  REQUIRE(deformed == 0);
}

TEST_CASE("projected basal area accepts and rejects each axial direction",
          "[cover_topo][prism]") {
  std::vector<int> lower = {0, 1, 2, 3};
  std::vector<int> upper = {4, 5, 6, 7};
  Cloud alongX = ringCloud(kYz, kYz, {50, 10, 10});
  REQUIRE(ring::discardExtraTetragonBlocks(lower, upper, alongX));
  Cloud xFalse = ringCloud(kXy, kYz, {50, 10, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, xFalse));
  Cloud xPartial = ringCloud(kXPartial, kYz, {50, 10, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, xPartial));
  Cloud xSecond = ringCloud(kYz, kXy, {50, 10, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, xSecond));
  Cloud xSecondPartial = ringCloud(kYz, kXPartial, {50, 10, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, xSecondPartial));

  Cloud alongY = ringCloud(kXz, kXz, {10, 50, 10});
  REQUIRE(ring::discardExtraTetragonBlocks(lower, upper, alongY));
  Cloud yFalse = ringCloud(kYz, kXz, {10, 50, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, yFalse));
  Cloud yPartial = ringCloud(kYPartial, kXz, {10, 50, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, yPartial));
  Cloud ySecond = ringCloud(kXz, kYz, {10, 50, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, ySecond));
  Cloud ySecondPartial = ringCloud(kXz, kYPartial, {10, 50, 10});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, ySecondPartial));

  Cloud alongZ = ringCloud(kXy, kXy, {10, 10, 50});
  REQUIRE(ring::discardExtraTetragonBlocks(lower, upper, alongZ));
  Cloud zPartial = ringCloud(kZPartial, kXy, {10, 10, 50});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, zPartial));
  Cloud zSecondPartial = ringCloud(kXy, kZPartial, {10, 10, 50});
  REQUIRE_FALSE(ring::discardExtraTetragonBlocks(lower, upper, zSecondPartial));
}

TEST_CASE("prism atom states keep an equal type and mix only the opposed pair",
          "[cover_topo][prism]") {
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3}};
  std::vector<int> atomTypes = {1, 1, 1, 1};
  std::vector<ring::strucType> same(4, ring::strucType::Prism);
  std::vector<ring::strucType> kinds = {ring::strucType::Prism};
  REQUIRE(ring::assignPrismType(rings, {0}, 4, kinds, atomTypes, same) == 0);
  REQUIRE(same[0] == ring::strucType::Prism);
  REQUIRE(atomTypes[0] == 4);

  std::vector<ring::strucType> deformed(4, ring::strucType::deformedPrism);
  kinds[0] = ring::strucType::HCbasal;
  REQUIRE(ring::assignPrismType(rings, {0}, 4, kinds, atomTypes, deformed) == 0);
  REQUIRE(deformed[0] == ring::strucType::deformedPrism);

  std::vector<ring::strucType> other(4, ring::strucType::HCbasal);
  kinds[0] = ring::strucType::deformedPrism;
  REQUIRE(ring::assignPrismType(rings, {0}, 4, kinds, atomTypes, other) == 0);
  REQUIRE(other[0] == ring::strucType::HCbasal);

  std::vector<ring::strucType> prismAtoms(4, ring::strucType::Prism);
  kinds[0] = ring::strucType::HCbasal;
  REQUIRE(ring::assignPrismType(rings, {0}, 4, kinds, atomTypes, prismAtoms) == 0);
  REQUIRE(prismAtoms[0] == ring::strucType::Prism);

  kinds[0] = ring::strucType::deformedPrism;
  REQUIRE(ring::assignPrismType(rings, {0}, 4, kinds, atomTypes, prismAtoms) == 0);
  REQUIRE(prismAtoms[0] == ring::strucType::mixedPrismRing);
}

TEST_CASE("a later frame recentres on the saved atom and reports a lost one",
          "[cover_topo][prism]") {
  Cloud cloud = prismCloud(2.75, {10, 10, 50});
  cloud.pts[0].z = 4.0;
  int anchor = -1;
  REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 1) == 0);
  REQUIRE(anchor == 0);
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(0.0, 1e-8));

  cloud.pts[0].z = 5.0;
  REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 4) == 0);
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(0.0, 1e-6));

  anchor = 42;
  const auto text = capture([&] {
    REQUIRE(ring::rmAxialTranslations(cloud, anchor, 1, 4) == 1);
  });
  REQUIRE(text.find("Lost atoms.") != std::string::npos);
}

TEST_CASE("prism analysis writes one tetragonal prism and skips empty sizes",
          "[cover_topo][prism]") {
  Cloud cloud = prismCloud(2.75, {10, 10, 50});
  auto byId = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, byId);
  auto rings = primitive::ringNetwork(rows, 6);
  int squares = 0;
  for (const auto &ring : rings) {
    squares += static_cast<int>(ring.size() == 4);
  }
  REQUIRE(squares >= 2);
  const auto directory = pidTemp("dseams-cover-topo");
  fs::remove_all(directory);
  fs::create_directories(directory);
  const auto path = directory.string() + "/";
  int anchor = 0;
  REQUIRE(ring::prismAnalysis(path, rings, rows, cloud, 6, anchor, 1, 1, false) ==
          0);
  std::ifstream input(path + "topoINT/nPrisms.dat");
  REQUIRE(input.good());
  std::string header;
  std::getline(input, header);
  REQUIRE(header.find("Num_of_prisms") != std::string::npos);
  int frame = 0;
  input >> frame;
  REQUIRE(frame == 1);
  bool sawFour = false;
  for (int size = 3; size <= 6; ++size) {
    int ringSize = 0;
    int prisms = 0;
    int deformed = 0;
    double height = 0;
    input >> ringSize >> prisms >> deformed >> height;
    REQUIRE(ringSize == size);
    if (size == 4) {
      sawFour = true;
      REQUIRE(prisms == 1);
      REQUIRE(deformed == 0);
      REQUIRE(height > 0.0);
    } else {
      REQUIRE(prisms == 0);
      REQUIRE(deformed == 0);
    }
  }
  REQUIRE(sawFour);
  fs::remove_all(directory);
}

TEST_CASE("bulk search rejects a thin or missing cell and a duplicated point",
          "[cover_topo][bulk]") {
  Cloud cloud = makeCloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, {10, 10, 10});
  cloud.boxLow.clear();
  auto rank = ring::positionRank(cloud);
  REQUIRE(rank.size() == 3);
  REQUIRE(rank[0] + rank[1] + rank[2] == 3);

  cloud.pts[1].x = 0;
  cloud.pts[1].y = 0;
  cloud.pts[1].z = 0;
  try {
    ring::positionRank(cloud);
    FAIL("duplicate wrapped positions must be refused");
  } catch (const std::invalid_argument &error) {
    REQUIRE(std::string(error.what()).find("share a wrapped position") !=
            std::string::npos);
  }

  std::vector<std::vector<int>> rings = {{0, 1, 2}};
  cloud.box = {10, -1, 10};
  cloud.boxLow = {0, 0, 0};
  auto missing = ring::buildRingSearchIndex(rings, cloud);
  REQUIRE(missing.ringShift.size() == 1);
  REQUIRE(missing.ringShift[0] == std::array<int, 3>({0, 0, 0}));
  auto kept = ring::cageRings(rings, cloud, 3.0);
  REQUIRE(kept.size() == 1);
  REQUIRE(kept[0] == rings[0]);

  cloud.box = {1.0, 8.0, 8.0};
  try {
    ring::cageRings(rings, cloud, 3.0);
    FAIL("a cell thinner than twice the cutoff must be refused");
  } catch (const std::invalid_argument &error) {
    REQUIRE(std::string(error.what()).find("cage search refused") !=
            std::string::npos);
  }

  Cloud tilted = makeCloud({{0, 0, 0}, {3, 0, 0}, {0.4, 3, 0.2}, {0, 0.2, 3}},
                           {12, 12, 12, 1.5, 0.4, 0.8});
  std::vector<std::vector<int>> tiltedRing = {{0, 1, 2, 3}};
  auto tiltedKept = ring::contractibleRings(tiltedRing, tilted);
  REQUIRE(tiltedKept.size() == 1);
  REQUIRE(tiltedKept[0].size() == 4);

  Cloud wide = makeCloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {20, 20, 20});
  std::vector<std::vector<int>> many;
  for (int i = 0; i < 300; ++i) {
    many.push_back({0, 1, 2});
  }
  auto outside = ring::contractibleRings(many, wide);
  REQUIRE(outside.size() == 300);
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(2);
  std::vector<std::vector<int>> insideKept;
  ring::RingSearchIndex nested;
#pragma omp parallel
  {
#pragma omp single
    {
      insideKept = ring::contractibleRings(many, wide);
      nested = ring::buildRingSearchIndex(rings, 3);
    }
  }
  REQUIRE(insideKept.size() == 300);
  REQUIRE(nested.ringsContainingAtom[0].size() == 1);
  omp_set_num_threads(1);
#endif
}

TEST_CASE("mixed rings survive when the cage lists do not name them",
          "[cover_topo][bulk]") {
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5}, {1, 2, 3, 4, 5, 6}};
  std::vector<ring::strucType> types = {ring::strucType::bothBasal,
                                        ring::strucType::bothPrismatic};
  std::vector<int> listDdc = {9};
  std::vector<int> listHc = {9};
  auto mixed = ring::findMixedRings(rings, types, listDdc, listHc);
  REQUIRE(mixed.size() == 2);
  REQUIRE(mixed[0] == 0);
  REQUIRE(mixed[1] == 1);
  REQUIRE(listDdc[0] == 9);
  REQUIRE(listHc[0] == 9);
}

TEST_CASE("bulk prisms stay classified when the same pair is tested again",
          "[cover_topo][bulk]") {
  Cloud cloud = prismCloud(2.75, {20, 20, 20});
  auto rows = prismRows();
  std::vector<std::vector<int>> squares = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<ring::strucType> types(2, ring::strucType::unclassified);
  std::vector<double> rmsd(8, -1.0);
  REQUIRE(prism3::findBulkPrisms(squares, types, rows, cloud, rmsd, 8.0) == 0);
  REQUIRE(types[0] == ring::strucType::Prism);
  REQUIRE(types[1] == ring::strucType::Prism);
  REQUIRE(prism3::findBulkPrisms(squares, types, rows, cloud, rmsd, 8.0) == 0);
  REQUIRE(types[0] == ring::strucType::Prism);
  REQUIRE(types[1] == ring::strucType::Prism);

  cloud.pts[5].x = 9.0;
  cloud.pts[6].y = 9.0;
  std::vector<ring::strucType> strict(2, ring::strucType::unclassified);
  REQUIRE(prism3::findBulkPrisms(squares, strict, rows, cloud, rmsd, 30.0) == 0);
  REQUIRE(strict[0] == ring::strucType::Prism);
  REQUIRE(strict[1] == ring::strucType::Prism);
  REQUIRE(prism3::findBulkPrisms(squares, strict, rows, cloud, rmsd, 30.0) == 0);
  REQUIRE(strict[0] == ring::strucType::Prism);
}

TEST_CASE("basal hexagon tests reject a short ring, a gap, and a missing bond",
          "[cover_topo][bulk]") {
  auto rows = selfRows(6);
  std::vector<int> full = {0, 1, 2, 3, 4, 5};
  std::vector<int> shortRing = {0, 1, 2, 3, 4};
  REQUIRE_FALSE(ring::basalConditions(rows, full, shortRing));
  REQUIRE_FALSE(ring::basalConditions(rows, full, full));
  std::vector<int> high = {6, 1, 2, 3, 4, 5};
  REQUIRE_FALSE(ring::basalConditions(rows, high, full));
  std::vector<int> negative = {-1, 1, 2, 3, 4, 5};
  REQUIRE_FALSE(ring::basalConditions(rows, negative, full));
  rows[0].clear();
  REQUIRE_FALSE(ring::basalConditions(rows, full, full));
  link(rows, 0, 1);
  REQUIRE_FALSE(ring::basalConditions(rows, full, full));

  ring::NeighBond bonds;
  std::vector<std::vector<int>> wide(70);
  wide[0] = {0, -1, 3, 80};
  wide[3] = {3, 0};
  bonds.build(wide);
  REQUIRE(bonds(0, 3));
  REQUIRE_FALSE(bonds(0, -1));
  REQUIRE_FALSE(bonds(0, 80));
  REQUIRE_FALSE(bonds(-2, 3));
  REQUIRE_FALSE(bonds(3, 90));
  REQUIRE_FALSE(ring::basalConditions(bonds, full, shortRing));
  REQUIRE_FALSE(ring::basalConditions(bonds, high, full));
}

TEST_CASE("an equatorial atom past the index is not a cage", "[cover_topo][bulk]") {
  std::vector<std::vector<int>> rings = {{100, 0, 1, 2, 3, 4}};
  auto index = ring::buildRingSearchIndex(rings, 5);
  std::vector<int> peripheral;
  REQUIRE_FALSE(ring::conditionOneDDC(rings, peripheral, 0, index));
  REQUIRE(peripheral.empty());

  rings[0][0] = -3;
  REQUIRE_FALSE(ring::conditionOneDDC(rings, peripheral, 0, index));
}

TEST_CASE("a hexagonal cage still counts when a neighbour index is out of range",
          "[cover_topo][bulk]") {
  const double hc[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  auto cloud = makeCloud(hc, {40, 40, 40});
  auto listed = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, listed);
  rows[0].push_back(-1);
  rows[0].push_back(11);
  rows[1].push_back(100);
  auto rings = ring::getSingleRingSize(primitive::ringNetwork(rows, 6), 6);
  REQUIRE_FALSE(rings.empty());
  std::vector<ring::strucType> types(rings.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  auto found = ring::findHC(rings, types, rows, cages);
  REQUIRE_FALSE(found.empty());
  REQUIRE_FALSE(cages.empty());
  REQUIRE(cages[0].type == cage::cageType::HexC);
  auto shared = rings;
  shared[0].push_back(-1);
  if (shared.size() > 1) {
    shared[1].push_back(-1);
  }
  std::vector<ring::strucType> sharedTypes(shared.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> sharedCages;
  auto sharedFound = ring::findHC(shared, sharedTypes, rows, sharedCages);
  REQUIRE(sharedFound.size() <= found.size());

  auto small = ring::buildRingSearchIndex(rings, static_cast<int>(rows.size()) - 1);
  rows[2].push_back(static_cast<int>(rows.size()) - 1);
  std::vector<ring::strucType> again(rings.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> more;
  auto second = ring::findHC(rings, again, rows, more, small);
  REQUIRE(second.size() == 4);
  REQUIRE(more.size() == 1);
  REQUIRE(more[0].type == cage::cageType::HexC);

#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(2);
  std::vector<ring::strucType> parallelTypes(rings.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> parallelCages;
  auto parallelFound = ring::findHC(rings, parallelTypes, rows, parallelCages);
  std::vector<ring::strucType> nestedTypes(rings.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> nestedCages;
  std::vector<int> nestedFound;
#pragma omp parallel
  {
#pragma omp single
    nestedFound = ring::findHC(rings, nestedTypes, rows, nestedCages);
  }
  REQUIRE(parallelFound.size() == found.size());
  REQUIRE(parallelCages.size() == cages.size());
  omp_set_num_threads(1);
#pragma omp parallel
  {
#pragma omp single
    {
      std::vector<ring::strucType> nestedTypes(rings.size(),
                                              ring::strucType::unclassified);
      std::vector<cage::Cage> nestedCages;
      auto nested = ring::findHC(rings, nestedTypes, rows, nestedCages);
      REQUIRE(nested.size() == found.size());
    }
  }
#endif
}

TEST_CASE("filler rings keep one double-diamond cage and drop an empty chain",
          "[cover_topo][bulk]") {
  const double ddc[14][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073}};
  auto cloud = makeCloud(ddc, {54, 54, 54});
  auto listed = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, listed);
  auto rings = ring::getSingleRingSize(primitive::ringNetwork(rows, 6), 6);
  REQUIRE(rings.size() >= 7);
  auto index = ring::buildRingSearchIndex(rings, 14);
  int equator = -1;
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  for (int i = 0; i < static_cast<int>(rings.size()); ++i) {
    chosen.clear();
    if (ring::equatorialConditions(rings, i, index, nullptr, &chosen, &rank,
                                   nullptr)) {
      equator = i;
      break;
    }
  }
  REQUIRE(equator >= 0);
  REQUIRE(chosen.size() == 6);
  ring::EquatorChoice plainRaw;
  REQUIRE(ring::equatorialConditions(rings, equator, index, nullptr, nullptr, &rank,
                                     &plainRaw));
  REQUIRE(plainRaw.n == 6);

  auto grow = [&](int copies, bool oddsOnly) {
    auto grown = rings;
    int nextAtom = 14;
    const auto eq = grown[static_cast<std::size_t>(equator)];
    for (int k = 0; k < 6; ++k) {
      if (oddsOnly && k % 2 == 0) {
        continue;
      }
      for (int copy = 0; copy < copies; ++copy) {
        grown.push_back({eq[static_cast<std::size_t>(k)],
                         eq[static_cast<std::size_t>((k + 1) % 6)],
                         eq[static_cast<std::size_t>((k + 2) % 6)], nextAtom,
                         nextAtom + 1, nextAtom + 2});
        nextAtom += 3;
      }
    }
    return grown;
  };

  auto both = grow(6, false);
  auto bothIndex = ring::buildRingSearchIndex(both, 14 + 6 * 6 * 3);
  std::vector<int> bothRank(both.size());
  std::iota(bothRank.begin(), bothRank.end(), 0);
  std::vector<int> bothChosen;
  std::vector<int> peripheral;
  ring::EquatorChoice raw;
  REQUIRE(ring::equatorialConditions(both, equator, bothIndex, &peripheral,
                                     &bothChosen, &bothRank, &raw));
  REQUIRE(bothChosen.size() == 6);
  REQUIRE(raw.n == 6);
  REQUIRE_FALSE(peripheral.empty());

  ring::EquatorChoice rawOnly;
  REQUIRE(ring::equatorialConditions(both, equator, bothIndex, nullptr, nullptr,
                                     &bothRank, &rawOnly));
  REQUIRE(rawOnly.n == 6);
  std::vector<int> peripheryOnly;
  REQUIRE(ring::equatorialConditions(both, equator, bothIndex, &peripheryOnly,
                                     nullptr, nullptr, nullptr));
  REQUIRE_FALSE(peripheryOnly.empty());
  REQUIRE(ring::equatorialConditions(both, equator, bothIndex, nullptr, nullptr,
                                     nullptr, nullptr));

  auto narrow = grow(1, false);
  auto narrowIndex =
      ring::buildRingSearchIndex(narrow, 14 + 6 * 1 * 3);
  std::vector<int> narrowRank(narrow.size());
  std::iota(narrowRank.begin(), narrowRank.end(), 0);
  std::vector<int> narrowChosen;
  ring::EquatorChoice narrowRaw;
  REQUIRE(ring::equatorialConditions(narrow, equator, narrowIndex, nullptr,
                                     &narrowChosen, &narrowRank, &narrowRaw));
  REQUIRE(narrowChosen.size() == 6);
  REQUIRE(narrowRaw.n == 6);
  ring::EquatorChoice narrowRawOnly;
  REQUIRE(ring::equatorialConditions(narrow, equator, narrowIndex, nullptr, nullptr,
                                     &narrowRank, &narrowRawOnly));
  REQUIRE(narrowRawOnly.n == 6);
  REQUIRE(ring::equatorialConditions(narrow, equator, narrowIndex, nullptr, nullptr,
                                     nullptr, nullptr));

  auto extra = both;
  auto addExtra = [&](int a, int b, int c, int x, int y, int z) {
    extra.push_back({a, b, c, x, y, z});
  };
  const auto &eq = rings[static_cast<std::size_t>(equator)];
  auto at = [&](int k) { return eq[static_cast<std::size_t>(k % 6)]; };
  addExtra(at(0), at(1), at(2), 200, 201, 202);
  addExtra(at(2), at(3), at(4), 200, 201, 203);
  addExtra(at(4), at(5), at(0), 200, 202, 204);
  addExtra(at(0), at(1), at(2), 210, 211, 212);
  addExtra(at(2), at(3), at(4), 210, 211, 213);
  addExtra(at(4), at(5), at(0), 210, 212, 214);
  std::unordered_set<int> eqAtoms(eq.begin(), eq.end());
  std::unordered_map<int, int> alias;
  int aliasNext = 400;
  for (int i = 0; i < static_cast<int>(rings.size()); ++i) {
    if (i == equator) {
      continue;
    }
    std::vector<int> copy;
    for (int atom : rings[static_cast<std::size_t>(i)]) {
      if (eqAtoms.count(atom) != 0) {
        copy.push_back(atom);
      } else {
        auto [it, inserted] = alias.emplace(atom, aliasNext);
        if (inserted) {
          ++aliasNext;
        }
        copy.push_back(it->second);
      }
    }
    extra.push_back(std::move(copy));
  }
  auto extraIndex = ring::buildRingSearchIndex(extra, aliasNext);
  std::vector<int> extraRank(extra.size());
  std::iota(extraRank.begin(), extraRank.end(), 0);
  std::vector<int> extraChosen;
  REQUIRE(ring::equatorialConditions(extra, equator, extraIndex, nullptr,
                                     &extraChosen, &extraRank, nullptr));
  REQUIRE(extraChosen.size() == 6);

  auto odds = grow(6, true);
  auto oddsIndex = ring::buildRingSearchIndex(odds, 14 + 3 * 6 * 3);
  std::vector<int> oddsRank(odds.size());
  std::iota(oddsRank.begin(), oddsRank.end(), 0);
  std::vector<int> oddsChosen;
  REQUIRE(ring::equatorialConditions(odds, equator, oddsIndex, nullptr, &oddsChosen,
                                     &oddsRank, nullptr));
  REQUIRE(oddsChosen.size() == 6);

  // Chains that share atoms, then one atom outside the cloud, are removed
  // before a cage is kept.
  std::vector<std::vector<int>> loose = {{0, 1, 2, 3, 4, 5}};
  int fresh = 6;
  auto add = [&](int a, int b, int c, std::initializer_list<int> extra) {
    std::vector<int> ring = {a, b, c};
    ring.insert(ring.end(), extra.begin(), extra.end());
    while (ring.size() < 6) {
      ring.push_back(fresh++);
    }
    loose.push_back(ring);
  };
  for (int k = 0; k < 6; ++k) {
    for (int copy = 0; copy < 7; ++copy) {
      add(k, (k + 1) % 6, (k + 2) % 6, {});
    }
  }
  // Two even chains of middles 1, so the later chain's key is worse.
  add(0, 1, 2, {150, 151, 152});
  add(2, 3, 4, {150, 151, 153});
  add(4, 5, 0, {150, 152, 154});
  add(0, 1, 2, {190, 191, 192});
  add(2, 3, 4, {190, 191, 193});
  add(4, 5, 0, {190, 192, 194});
  // Odd middles 4 is dead against even bucket 1. Odd middles 2 stays live.
  add(1, 2, 3, {180, 181, 182});
  add(5, 0, 1, {180, 181, 183});
  add(3, 4, 5, {183, 181, 185});
  add(1, 2, 3, {170, 173, 172});
  add(3, 4, 5, {170, 173, 174});
  add(5, 0, 1, {174, 173, 177});
  auto looseIndex = ring::buildRingSearchIndex(loose, 200);
  std::vector<int> looseRank(loose.size());
  std::iota(looseRank.begin(), looseRank.end(), 0);
  ring::EquatorChoice looseRaw;
  REQUIRE_FALSE(ring::equatorialConditions(loose, 0, looseIndex, nullptr, nullptr,
                                           &looseRank, &looseRaw));
  REQUIRE(looseRaw.n == 0);
  REQUIRE_FALSE(ring::equatorialConditions(loose, 0, looseIndex, nullptr, nullptr,
                                           nullptr, nullptr));

  std::vector<std::vector<int>> dropped = {{0, 1, 2, 3, 4, 5}};
  int freshDrop = 6;
  auto addDrop = [&](int a, int b, int c, int x, int y, int z) {
    dropped.push_back({a, b, c, x, y, z});
  };
  for (int k = 0; k < 6; ++k) {
    for (int copy = 0; copy < 7; ++copy) {
      dropped.push_back({k, (k + 1) % 6, (k + 2) % 6, freshDrop, freshDrop + 1,
                         freshDrop + 2});
      freshDrop += 3;
    }
  }
  addDrop(0, 1, 2, 150, 151, 4000);
  addDrop(2, 3, 4, 150, 151, 4000);
  addDrop(4, 5, 0, 150, 151, 4000);
  auto droppedIndex = ring::buildRingSearchIndex(dropped, 220);
  std::vector<int> droppedRank(dropped.size());
  std::iota(droppedRank.begin(), droppedRank.end(), 0);
  REQUIRE_FALSE(ring::equatorialConditions(dropped, 0, droppedIndex, nullptr, nullptr,
                                           &droppedRank, nullptr));

  std::vector<std::vector<int>> dead = {{0, 1, 2, 3, 4, 5}};
  int freshDead = 6;
  auto addDead = [&](int a, int b, int c, std::initializer_list<int> extra) {
    std::vector<int> ring = {a, b, c};
    ring.insert(ring.end(), extra.begin(), extra.end());
    while (ring.size() < 6) {
      ring.push_back(freshDead++);
    }
    dead.push_back(std::move(ring));
  };
  for (int k = 0; k < 6; ++k) {
    for (int copy = 0; copy < 7; ++copy) {
      addDead(k, (k + 1) % 6, (k + 2) % 6, {});
    }
  }
  addDead(0, 1, 2, {150, 151, 152});
  addDead(2, 3, 4, {150, 151, 153});
  addDead(4, 5, 0, {150, 152, 154});
  addDead(1, 2, 3, {180, 181, 182});
  addDead(5, 0, 1, {180, 181, 183});
  addDead(3, 4, 5, {183, 181, 185});
  auto deadIndex = ring::buildRingSearchIndex(dead, 220);
  std::vector<int> deadRank(dead.size());
  std::iota(deadRank.begin(), deadRank.end(), 0);
  REQUIRE_FALSE(ring::equatorialConditions(dead, 0, deadIndex, nullptr, nullptr,
                                           &deadRank, nullptr));

  std::vector<std::vector<int>> evenDead = {{0, 1, 2, 3, 4, 5}};
  int freshEven = 6;
  auto addEven = [&](int a, int b, int c, std::initializer_list<int> extra) {
    std::vector<int> ring = {a, b, c};
    ring.insert(ring.end(), extra.begin(), extra.end());
    while (ring.size() < 6) {
      ring.push_back(freshEven++);
    }
    evenDead.push_back(std::move(ring));
  };
  for (int k = 0; k < 6; ++k) {
    for (int copy = 0; copy < 7; ++copy) {
      addEven(k, (k + 1) % 6, (k + 2) % 6, {});
    }
  }
  addEven(0, 1, 2, {150, 151, 152});
  addEven(2, 3, 4, {150, 151, 153});
  addEven(4, 5, 0, {150, 152, 154});
  addEven(1, 2, 3, {170, 173, 172});
  addEven(3, 4, 5, {170, 173, 174});
  addEven(5, 0, 1, {174, 173, 177});
  auto evenDeadIndex = ring::buildRingSearchIndex(evenDead, 220);
  REQUIRE_FALSE(ring::equatorialConditions(evenDead, 0, evenDeadIndex, nullptr,
                                           nullptr, nullptr, nullptr));
}

TEST_CASE("a duplicated peripheral uses the triple loop and a disjoint set does not",
          "[cover_topo][bulk]") {
  const double ddc[14][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073}};
  auto cloud = makeCloud(ddc, {54, 54, 54});
  auto rows = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
  auto rings = ring::getSingleRingSize(primitive::ringNetwork(rows, 6), 6);
  auto index = ring::buildRingSearchIndex(rings, 14);
  int equator = -1;
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  for (int i = 0; i < static_cast<int>(rings.size()); ++i) {
    chosen.clear();
    if (ring::equatorialConditions(rings, i, index, nullptr, &chosen, &rank,
                                   nullptr)) {
      equator = i;
      break;
    }
  }
  REQUIRE(equator >= 0);
  auto twin = rings;
  twin.push_back(rings[static_cast<std::size_t>(chosen[0])]);
  auto twinIndex = ring::buildRingSearchIndex(twin, 14);
  REQUIRE(ring::equatorialConditions(twin, equator, twinIndex, nullptr, nullptr,
                                     nullptr, nullptr));
  std::vector<int> peripheral;
  REQUIRE(ring::equatorialConditions(twin, equator, twinIndex, &peripheral,
                                     nullptr, nullptr, nullptr));
  REQUIRE(peripheral.size() >= 6);

  std::vector<std::vector<int>> apart = {{0, 1, 2, 3, 4, 5}};
  int atom = 6;
  for (int k = 0; k < 6; ++k) {
    apart.push_back({k, (k + 1) % 6, (k + 2) % 6, atom, atom + 1, atom + 2});
    atom += 3;
  }
  auto apartIndex = ring::buildRingSearchIndex(apart, atom);
  REQUIRE_FALSE(ring::equatorialConditions(apart, 0, apartIndex, nullptr, nullptr,
                                           nullptr, nullptr));

  std::vector<std::vector<int>> past = {{80, 1, 2, 3, 4, 5}};
  auto pastIndex = ring::buildRingSearchIndex(past, 6);
  REQUIRE_FALSE(ring::equatorialConditions(past, 0, pastIndex, nullptr, nullptr,
                                           nullptr, nullptr));

  auto hcRows = rows;
  auto hcRings = rings;
  hcRings[static_cast<std::size_t>(equator)].push_back(400);
  std::vector<ring::strucType> hcTypes(hcRings.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> hcCages;
  auto hcFound = ring::findHC(hcRings, hcTypes, hcRows, hcCages);
  REQUIRE(hcFound.size() <= hcRings.size());
}
