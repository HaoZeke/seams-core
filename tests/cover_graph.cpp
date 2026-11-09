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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

namespace {

Cloud makeCloud(const std::vector<std::array<double, 3>> &xyz,
                const std::vector<int> &types, std::vector<double> box,
                std::vector<double> boxLow = {0.0, 0.0, 0.0}) {
  Cloud cloud;
  cloud.box = std::move(box);
  cloud.boxLow = std::move(boxLow);
  cloud.currentFrame = 1;
  for (std::size_t i = 0; i < xyz.size(); i++) {
    molSys::Point<double> pt;
    pt.type = types.empty() ? 1 : types[i];
    pt.atomID = static_cast<int>(i);
    pt.molID = static_cast<int>(i);
    pt.x = xyz[i][0];
    pt.y = xyz[i][1];
    pt.z = xyz[i][2];
    pt.inSlice = true;
    cloud.idIndexMap[pt.atomID] = static_cast<int>(i);
    cloud.pts.push_back(pt);
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

bool listsId(const std::vector<int> &row, int id) {
  return std::find(row.begin(), row.end(), id) != row.end();
}

class CerrCapture {
public:
  CerrCapture() : prev_(std::cerr.rdbuf(sink_.rdbuf())) {}
  ~CerrCapture() { std::cerr.rdbuf(prev_); }
  std::string str() const { return sink_.str(); }

private:
  std::ostringstream sink_;
  std::streambuf *prev_;
};

std::filesystem::path scratchDir(const char *name) {
  auto dir = std::filesystem::temp_directory_path() / name;
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  return dir;
}

} // namespace

// Unreachable without rewriting production code, recorded so a later
// coverage pass can leave them:
//
// cluster.cpp:231. largestIceCluster returns at the empty-nClusters test,
// so totalClusters == 0 is not evaluated.
// cluster.cpp:410. chill::getIceTypeNoPrint overwrites the unclassified
// default for every atom (water when the bond count is not 4, otherwise
// cubic, hexagonal, interfacial, or water).
//
// neighbours.cpp:724. orthoFullList is reached only after hasPeriodicBox,
// which already requires a positive recovered diagonal.
// neighbours.cpp:1065-1066. An atom has one type, so an unlike (typeI !=
// typeJ) pair never has the same cloud index on both sides.
// neighbours.cpp:1368. Every public k-nearest entry returns before
// nominatePacked when nop or k is not positive.
// neighbours.cpp:1410-1411. A non-positive ortho edge makes
// nominateByLinkcell return success before the in-tree fallback, and a
// six-element box takes the triclinic fallback instead.
// neighbours.cpp:1481-1483. stamp is unsigned, starts at 1, and increments
// once per owner. nop is a signed int, so the counter cannot wrap to 0.
// neighbours.cpp:2013 and 2029 and 2043 and 2056. update() rebuilds
// candidates whenever the atom count changes, and the rebuild only stores
// indices inside the current cloud.
//
// franzblau.cpp:55-65 run only when DSEAMS_PIN is set in the process
// environment before libyodaLib's static WorkerPlacement constructor. A
// Catch2 case starts too late to set that. The coverage test environment
// sets the variable.
// franzblau.cpp:301-302. refreshNear runs only after fastEnumerate has
// filled balls for the same vertex count, so off.size() is n + 1.
// franzblau.cpp:373. dist() is called with ring members, which the
// breadth-first walk only admits inside [0, n).
// franzblau.cpp:566. With a null through-mask the walk's minimum vertex
// is the source, so every touched index is >= src.
// franzblau.cpp:1434. nThreads is at least 1 and scratches is resized to
// that count before the empty test.
// franzblau.cpp:1475. fastEnumerate sizes bySource to the vertex count,
// and the masked update keeps that count.
// franzblau.cpp:1745-1747. ballsStale starts true and fastEnumerate clears
// it on the first frame. Nothing sets it true again, and the masked update
// is only reached on a later frame.
// ring.cpp:99 is the exception landing pad of findsCommonElements. The
// body is already executed. std::find is noexcept; the pad runs only if
// vector reallocation throws std::bad_alloc.
//
// WorkerPlacement, BoundedBalls::dist, levelsFrom, and the quartile
// helpers above are not noexcept as a group: vector growth can throw.
// The pads stay dark because the throwing condition is allocation failure.

TEST_CASE("odd quartile average drops a high outlier", "[cover_graph]") {
  const double odd =
      gen::getAverageWithoutOutliers({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1000.0});
  REQUIRE_THAT(odd, Catch::Matchers::WithinAbs(3.5, 1e-12));
  const double low =
      gen::getAverageWithoutOutliers({-100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  REQUIRE_THAT(low, Catch::Matchers::WithinAbs(3.5, 1e-12));
  const double kept =
      gen::getAverageWithoutOutliers({1.0, 2.0, 3.0, 4.0, 100.0});
  REQUIRE_THAT(kept, Catch::Matchers::WithinAbs(22.0, 1e-12));
}

TEST_CASE("prettyPrint reports success when the path is not a file",
          "[cover_graph]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0}}, {}, {4.0, 4.0, 4.0});
  REQUIRE(gen::prettyPrintYoda(cloud, ".") == 0);
}

TEST_CASE("generic box helpers cover ortho, tilt, and rejected cells",
          "[cover_graph]") {
  Cloud bare;
  bare.pts.push_back({});
  bare.nop = 1;
  bare.box = {10.0, 10.0};
  auto bad = gen::makeFracBox(bare);
  REQUIRE_FALSE(bad.ok);

  Cloud zero = bare;
  zero.box = {0.0, 10.0, 10.0};
  zero.boxLow = {1.0};
  auto flat = gen::makeFracBox(zero);
  REQUIRE_FALSE(flat.ok);
  REQUIRE_THAT(flat.ox, Catch::Matchers::WithinAbs(1.0, 1e-12));

  Cloud ortho = makeCloud({{0.0, 0.0, 0.0}, {9.5, 0.0, 0.0}}, {},
                          {10.0, 10.0, 10.0}, {});
  auto frame = gen::makeFracBox(ortho);
  REQUIRE(frame.ok);
  REQUIRE_FALSE(frame.triclinic);
  REQUIRE(gen::smithInside(frame, 1.0));
  const double wrapped = gen::periodicDistSq(frame, ortho, 0, 1);
  REQUIRE_THAT(wrapped, Catch::Matchers::WithinAbs(0.25, 1e-9));

  Cloud neg = ortho;
  neg.box = {-10.0, 10.0, 10.0};
  neg.pts[0] = {};
  neg.pts[1].x = 0.0;
  neg.pts[1].y = -6.0;
  neg.pts[1].z = 0.0;
  auto negFrame = gen::makeFracBox(neg);
  REQUIRE_FALSE(negFrame.ok);
  const double raw = gen::periodicDistSq(negFrame, neg, 0, 1);
  REQUIRE(raw > 0.0);
  const auto rel = gen::relDistFromPoint(negFrame, neg, 0, 9.5, 0.0, 0.0);
  REQUIRE(rel[0] < 0.0);
  const auto down = gen::relDistFromPoint(negFrame, neg, 0, 0.0, -6.0, 0.0);
  REQUIRE(down[1] < 5.0);

  Cloud tilt = makeCloud({{0.2, 0.2, 0.2}, {9.6, 0.3, 0.2}}, {},
                         {12.0, 10.0, 10.0, 2.0, 0.0, 0.0}, {0.5});
  auto tri = gen::makeFracBox(tilt);
  REQUIRE(tri.ok);
  REQUIRE(tri.triclinic);
  REQUIRE(gen::smithInside(tri, 0.01));
  REQUIRE_FALSE(gen::smithInside(tri, 1.0e6));
  const auto near =
      gen::triclinicMinImage(tilt, 0.2, 0.2, 0.2, 0.4, 0.2, 0.2);
  REQUIRE_THAT(near[0], Catch::Matchers::WithinAbs(-0.2, 1e-8));
  const auto far =
      gen::triclinicMinImage(tilt, 0.2, 0.2, 0.2, 9.6, 0.3, 0.2);
  REQUIRE(far[0] * far[0] + far[1] * far[1] + far[2] * far[2] < 4.0);
  const double farSq = gen::periodicDistSq(tri, tilt, 0, 1);
  REQUIRE(farSq < 4.0);
  const auto fromPt =
      gen::relDistFromPoint(tri, tilt, 0, 9.6, 0.3, 0.2);
  REQUIRE(fromPt[0] * fromPt[0] < 4.0);

  Cloud broken = tilt;
  broken.box = {0.0, 10.0, 10.0, 1.0, 0.0, 0.0};
  broken.boxLow = {1.0, 2.0};
  auto brokenFrame = gen::makeFracBox(broken);
  REQUIRE_FALSE(brokenFrame.ok);
  const auto manual = gen::triclinicMinImage(broken, 0.2, 0.2, 0.2, 1.2, 0.2,
                                             0.2);
  REQUIRE(std::isfinite(manual[0]));
  const double manualSq = gen::periodicDistSq(brokenFrame, broken, 0, 1);
  REQUIRE(std::isfinite(manualSq));
  const auto manualRel =
      gen::relDistFromPoint(brokenFrame, broken, 0, 9.0, 0.0, 0.0);
  REQUIRE(std::isfinite(manualRel[0]));

  gen::FracBox dead;
  REQUIRE_FALSE(gen::smithInside(dead, 1.0));
  dead.ok = true;
  REQUIRE(gen::smithInside(dead, 1.0e9));

#ifdef SEAMS_HAS_MINIMAGE
  mi_cell cell{};
  Cloud cellCloud = ortho;
  cellCloud.boxLow.clear();
  REQUIRE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.boxLow = {1.0};
  REQUIRE(gen::pointCloudCell(cellCloud, &cell));
  REQUIRE_THAT(cell.ox, Catch::Matchers::WithinAbs(1.0, 1e-12));
  cellCloud.boxLow = {1.0, 2.0};
  REQUIRE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.boxLow = {1.0, 2.0, 3.0};
  REQUIRE(gen::pointCloudCell(cellCloud, &cell));
  REQUIRE_THAT(cell.oz, Catch::Matchers::WithinAbs(3.0, 1e-12));
  cellCloud.box = {10.0, 10.0, 10.0, 1.0, 0.0, 0.0};
  REQUIRE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.box = {0.0, 10.0, 10.0};
  REQUIRE_FALSE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.box = {-1.0, 10.0, 10.0};
  REQUIRE_FALSE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.box = {10.0, 0.0, 10.0};
  REQUIRE_FALSE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.box = {10.0, 10.0, 0.0};
  REQUIRE_FALSE(gen::pointCloudCell(cellCloud, &cell));
  cellCloud.box = {10.0, 10.0};
  REQUIRE_FALSE(gen::pointCloudCell(cellCloud, &cell));
  double ed[3] = {1.0, 1.0, 1.0};
  REQUIRE_FALSE(gen::euclideanDelta(cellCloud, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ed));
#endif

  int one = 1;
  double slot = -1.0;
  gen::batchPeriodicDistSq(ortho, 0, &one, 0, &slot);
  REQUIRE_THAT(slot, Catch::Matchers::WithinAbs(-1.0, 1e-15));
  gen::batchPeriodicDistSq(ortho, 0, &one, 1, &slot);
  REQUIRE_THAT(slot, Catch::Matchers::WithinAbs(0.25, 1e-8));
  gen::batchPeriodicDistSq(tilt, 0, &one, 1, &slot);
  REQUIRE(slot < 4.0);

  Cloud skew = makeCloud({{0.0, 0.0, 0.0}, {0.0, 6.0, 0.0}}, {},
                         {20.0, 20.0, 20.0, 18.0, 0.0, 0.0});
  auto skewFrame = gen::makeFracBox(skew);
  REQUIRE(skewFrame.ok);
  REQUIRE(skewFrame.triclinic);
  const auto skewDr = gen::fracDelta(skewFrame, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0);
  const double skewR2 =
      skewDr[0] * skewDr[0] + skewDr[1] * skewDr[1] + skewDr[2] * skewDr[2];
  REQUIRE_FALSE(gen::smithInside(skewFrame, skewR2));
  const auto skewImg =
      gen::triclinicMinImage(skew, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0);
  REQUIRE(std::isfinite(skewImg[1]));
  REQUIRE(gen::periodicDistSq(skewFrame, skew, 0, 1) > 0.0);
  const auto skewRel = gen::relDistFromPoint(skewFrame, skew, 0, 0.0, 6.0, 0.0);
  REQUIRE(std::isfinite(skewRel[1]));
}

TEST_CASE("dump box helpers reject bad edges and keep a right angle ortho",
          "[cover_graph]") {
  double H[3][3] = {};
  double origin[3] = {};
  nneigh::dumpBoundsToH({}, {}, H, origin);
  REQUIRE_THAT(H[0][0], Catch::Matchers::WithinAbs(0.0, 1e-15));
  nneigh::dumpBoundsToH({4.0}, {1.0}, H, origin);
  REQUIRE_THAT(origin[0], Catch::Matchers::WithinAbs(1.0, 1e-15));
  REQUIRE_THAT(nneigh::dumpVolume({-2.0, 3.0, 4.0}, {}),
               Catch::Matchers::WithinAbs(24.0, 1e-12));

  std::vector<double> box;
  std::vector<double> boxLow;
  REQUIRE_FALSE(nneigh::setBoxFromLammpsEdges(box, boxLow, 0, 0, 0, 0, 1, 1,
                                             0, 0, 0));
  REQUIRE(box.empty());
  REQUIRE_FALSE(nneigh::setBoxFromLammpsEdges(box, boxLow, 0, 0, 0, 1, -1, 1,
                                             0, 0, 0));
  REQUIRE_FALSE(nneigh::setBoxFromLammpsEdges(box, boxLow, 0, 0, 0, 1, 1, 0,
                                             0, 0, 0));
  REQUIRE(nneigh::setBoxFromLammpsEdges(box, boxLow, 1, 2, 3, 4, 5, 6, 0, 0,
                                       0.5));
  REQUIRE(box.size() == 6);
  REQUIRE_THAT(box[5], Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE(nneigh::setBoxFromLammpsEdges(box, boxLow, 0, 0, 0, 4, 5, 6, 0, 0.5,
                                       0));
  REQUIRE_THAT(box[4], Catch::Matchers::WithinAbs(0.5, 1e-12));

  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 0, 1, 1, 90, 90, 90));
  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 1, 0, 1, 90, 90, 90));
  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 1, 1, 0, 90, 90, 90));
  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 1, 1, 1, 90, 90, 0));
  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 1, 1, 1, 90, 0, 90));
  REQUIRE(nneigh::setBoxFromLengthsAngles(box, boxLow, 10, 10, 10, 90, 90, 90));
  REQUIRE(box.size() == 3);
  REQUIRE(nneigh::setBoxFromLengthsAngles(box, boxLow, 10, 10, 10, 80, 90, 90));
  REQUIRE(box.size() == 6);
  REQUIRE_FALSE(
      nneigh::setBoxFromLengthsAngles(box, boxLow, 10, 10, 10, 90, 90, 190));
  REQUIRE(box.empty());
}

TEST_CASE("brute-force neighbour lists run when the cutoff list rejects",
          "[cover_graph]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0},
                          {0.4, 0.0, 0.0},
                          {3.0, 0.0, 0.0},
                          {0.2, 0.2, 0.0}},
                         {1, 1, 2, 2}, {10.0, 10.0, 10.0});
  CerrCapture err;
  {
    auto full = nneigh::neighListO(-1.0, cloud, 1);
    REQUIRE(full.size() == 4);
    REQUIRE(listsId(full[0], 1));
    REQUIRE_FALSE(listsId(full[0], 2));
    auto half = nneigh::halfNeighList(-1.0, cloud, 1);
    REQUIRE(listsId(half[0], 1));
    REQUIRE_FALSE(listsId(half[1], 0));
    auto mixed = nneigh::neighListPair(-1.0, cloud, 1, 2);
    REQUIRE(listsId(mixed[0], 3));
    auto none = nneigh::neighListPair(-1.0, cloud, 1, 9);
    REQUIRE(none[0].size() == 1);
    auto byIndex = nneigh::getNewNeighbourListByIndex(cloud, -1.0);
    REQUIRE(listsId(byIndex[0], 1));
    REQUIRE_FALSE(listsId(byIndex[0], 2));
  }
  REQUIRE(err.str().find("falling back to brute force") != std::string::npos);

  auto tilt = makeCloud({{0.0, 0.0, 0.0}, {9.6, 0.0, 0.0}}, {},
                        {12.0, 10.0, 10.0, 2.0, 0.0, 0.0});
  auto wrapped = nneigh::neighListO(-1.0, tilt, 1);
  REQUIRE(listsId(wrapped[0], 1));

  cloud.idIndexMap.erase(1);
  CerrCapture missing;
  auto rows = nneigh::neighListO(-1.0, cloud, 1);
  REQUIRE(rows[1].empty());
  REQUIRE(missing.str().find("idIndexMap") != std::string::npos);

  auto far = makeCloud({{0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}, {6.0, 0.0, 0.0}}, {},
                       {10.0, 10.0, 10.0});
  auto halfFar = nneigh::halfNeighList(-1.0, far, 1);
  REQUIRE(listsId(halfFar[0], 1));
  REQUIRE_FALSE(listsId(halfFar[0], 2));

}

TEST_CASE("cell-list rows report an empty subset and a useless cutoff",
          "[cover_graph]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {},
                         {30.0, 30.0, 30.0});
  std::vector<std::vector<int>> rows = {{7}};
  REQUIRE(nneigh::cellListRowsThreaded(cloud, {}, 3.0, rows));
  REQUIRE(rows.empty());
  rows = {{7}};
  REQUIRE_FALSE(nneigh::cellListRowsThreaded(cloud, {0, 1}, 0.0, rows));
  REQUIRE(rows.size() == 1);
  cloud.box = {0.0, 30.0, 30.0};
  REQUIRE_FALSE(nneigh::cellListRowsThreaded(cloud, {0, 1}, 3.0, rows));
  cloud.box = {30.0, 30.0, 30.0};
  REQUIRE(nneigh::cellListRowsThreaded(cloud, {0, 1}, 3.0, rows));
  REQUIRE(listsId(rows[0], 1));
}

TEST_CASE("large ortho list hits the borderline shell and a short cell",
          "[cover_graph]") {
  std::vector<std::array<double, 3>> xyz;
  xyz.reserve(2200);
  for (int ix = 0; ix < 13; ix++) {
    for (int iy = 0; iy < 13; iy++) {
      for (int iz = 0; iz < 13; iz++) {
        xyz.push_back({2.0 * ix, 2.0 * iy, 2.0 * iz});
      }
    }
  }
  xyz[1] = {40.0, 0.0, 0.0};
  xyz[2] = {0.4, 0.0, 0.0};
  auto cloud = makeCloud(xyz, {}, {40.0, 40.0, 40.0});
  auto rows = nneigh::neighListO(4.0, cloud, 1);
  REQUIRE(rows.size() == xyz.size());
  REQUIRE(rows[0].size() > 1);

  auto wide = nneigh::neighListO(20.0, cloud, 1);
  REQUIRE(wide[0].size() > rows[0].size());

  CerrCapture err;
  auto brute = nneigh::getNewNeighbourListByIndex(cloud, -1.0);
  REQUIRE(brute.size() == xyz.size());
  REQUIRE(listsId(brute[0], 2));
  REQUIRE(err.str().find("falling back to brute force") != std::string::npos);
}

TEST_CASE("in-tree k-nearest matches the periodic heap when linkcell refuses",
          "[cover_graph]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.0, 1.2, 0.0},
                          {0.0, 0.0, 3.0},
                          {8.0, 8.0, 8.0},
                          {0.2, 0.2, 0.0}},
                         {1, 1, 1, 2, 1, 1}, {10.0, 10.0, 10.0});
  std::vector<std::vector<int>> rows;
  CerrCapture err;
  // Hint ~0.04 makes linkcell's mesh exceed its cell cap and throws.
  // The in-tree grid is the same edge and still fits in memory.
  rows = nneigh::kNearestNeighbourList(cloud, 2, 0.05, 1, true);
  REQUIRE(err.str().find("falling back") != std::string::npos);
  REQUIRE(listsId(rows[0], 1));
  REQUIRE(listsId(rows[0], 5));
  REQUIRE_FALSE(listsId(rows[0], 3));

  auto both = nneigh::kNearestNeighbourPair(cloud, 2, 0.05, 1);
  REQUIRE(listsId(both.first[0], 1));
  REQUIRE(both.second[0].size() >= both.first[0].size());

  auto typed = nneigh::kNearestNeighbourList(cloud, 2, 0.05,
                                             std::vector<int>{1, 2}, false);
  REQUIRE(typed[0].size() > 1);
  auto typedPair = nneigh::kNearestNeighbourPair(cloud, 1, 0.05,
                                                 std::vector<int>{});
  REQUIRE(typedPair.first.empty());
  REQUIRE(typedPair.second.empty());

  auto tilt = cloud;
  tilt.box = {12.0, 10.0, 10.0, 2.0, 0.5, 0.0};
  CerrCapture tiltErr;
  rows = nneigh::kNearestNeighbourList(tilt, 1, 0.05, 1, false);
  REQUIRE(tiltErr.str().find("falling back") != std::string::npos);
  REQUIRE(rows[0].size() >= 2);

  Cloud emptyBox = cloud;
  emptyBox.box.clear();
  auto noBox = nneigh::kNearestNeighbourList(emptyBox, 1, 1.0, 1, true);
  REQUIRE(noBox[0].size() == 1);
  Cloud flat = cloud;
  flat.box = {0.0, 10.0, 10.0};
  auto flatRows = nneigh::kNearestNeighbourList(flat, 1, 1.0, 1, true);
  REQUIRE(flatRows[0].size() == 1);

  const auto shells = nneigh::shellSeparation(cloud, 1, 1);
  REQUIRE(shells.first > 0.0);
  REQUIRE(shells.second > 0.0);
  const auto lonely = nneigh::shellSeparation(cloud, 4, 2);
  REQUIRE_THAT(lonely.first, Catch::Matchers::WithinAbs(0.0, 1e-15));
  const auto tilted = nneigh::shellSeparation(tilt, 1, 1);
  REQUIRE(tilted.first > 0.0);
  Cloud bad = cloud;
  bad.box = {-10.0, 10.0, 10.0};
  const auto badShell = nneigh::shellSeparation(bad, 1, 1);
  REQUIRE(std::isfinite(badShell.first));

  auto unlike = nneigh::nearestUnlike(cloud, 1, 2);
  REQUIRE(unlike.size() == 5);
  auto missing = nneigh::nearestUnlike(cloud, 1, 9);
  REQUIRE(missing.empty());
  auto mutual = nneigh::mutualNearestUnlike(cloud, 1, 2);
  REQUIRE_FALSE(mutual.empty());
}

TEST_CASE("bond graph names and a skin that misses the cutoff list",
          "[cover_graph]") {
  REQUIRE(std::string(nneigh::bondGraphName(nneigh::BondGraph::Cutoff)) ==
          "cutoff");
  REQUIRE(std::string(nneigh::bondGraphName(nneigh::BondGraph::KnnMutual)) ==
          "knn");
  REQUIRE(std::string(nneigh::bondGraphName(nneigh::BondGraph::KnnUnion)) ==
          "knn-union");
  const auto bogus =
      static_cast<nneigh::BondGraph>(static_cast<int>(nneigh::BondGraph::KnnUnion) +
                                     3);
  REQUIRE(std::string(nneigh::bondGraphName(bogus)) == "cutoff");

  auto cloud = makeCloud({{0.0, 0.0, 0.0},
                          {0.3, 0.0, 0.0},
                          {0.3, 0.4, 0.0},
                          {5.0, 5.0, 5.0}},
                         {1, 1, 1, 2}, {10.0, 10.0, 10.0});
  nneigh::SkinNeighborList skin(-1.0, 0.4, 1, nneigh::BondGraph::Cutoff);
  const auto &bonds = skin.update(cloud);
  REQUIRE(listsId(bonds[0], 1));
  REQUIRE(skin.lastRebuilt());
  cloud.pts[1].x = 0.35;
  const auto &held = skin.update(cloud);
  REQUIRE_FALSE(skin.lastRebuilt());
  REQUIRE(listsId(held[0], 1));
  cloud.box[0] = 14.0;
  REQUIRE(skin.update(cloud).size() == 4);
  REQUIRE(skin.lastRebuilt());
  cloud.box.push_back(0.0);
  REQUIRE(skin.update(cloud).size() == 4);
  REQUIRE(skin.lastRebuilt());

  nneigh::SkinNeighborList knn(1.5, 0.2, 1, nneigh::BondGraph::KnnMutual, 1);
  const auto &mutual = knn.update(cloud);
  REQUIRE(listsId(mutual[0], 1));
  REQUIRE(mutual[2].size() == 1);
  nneigh::SkinNeighborList uni(1.5, 0.2, 1, nneigh::BondGraph::KnnUnion, 1);
  const auto &unionBonds = uni.update(cloud);
  REQUIRE(unionBonds[0].size() >= mutual[0].size());

  cloud.idIndexMap.erase(1);
  nneigh::SkinNeighborList missing(-1.0, 0.4, 1, nneigh::BondGraph::Cutoff);
  CerrCapture skinErr;
  const auto &dropped = missing.update(cloud);
  REQUIRE(dropped.size() == 4);
  REQUIRE(skinErr.str().find("falling back to brute force") != std::string::npos);
}

TEST_CASE("three-ring overlap and wrapped triplets", "[cover_graph]") {
  const std::vector<int> a{1, 2, 3};
  const std::vector<int> b{3, 4, 5};
  const std::vector<int> c{3, 8, 9};
  const std::vector<int> d{7, 8};
  REQUIRE(ring::commonElementsInThreeRings(a, b, c));
  REQUIRE_FALSE(ring::commonElementsInThreeRings(a, b, d));
  REQUIRE_FALSE(ring::commonElementsInThreeRings(a, d, c));
  REQUIRE_FALSE(ring::commonElementsInThreeRings({}, b, c));

  const std::vector<int> hex{0, 1, 2, 3, 4, 5};
  REQUIRE(ring::findTripletInRing(hex, {0, 1, 2}));
  REQUIRE(ring::findTripletInRing(hex, {2, 1, 0}));
  REQUIRE(ring::findTripletInRing(hex, {5, 0, 1}));
  REQUIRE(ring::findTripletInRing(hex, {1, 0, 5}));
  REQUIRE_FALSE(ring::findTripletInRing(hex, {0, 2, 4}));
  REQUIRE_FALSE(ring::findTripletInRing({}, {0, 1, 2}));

  REQUIRE(ring::hasCommonElements(b, a));
  REQUIRE(ring::hasCommonElements(a, b));
  REQUIRE_FALSE(ring::hasCommonElements(a, d));
}

TEST_CASE("empty ice cluster, a short mask, and a chill mix", "[cover_graph]") {
  auto dir = scratchDir("dseams-cover-graph-cluster");
  auto cloud = makeCloud({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}},
                         {}, {10.0, 10.0, 10.0});
  cloud.pts[0].inSlice = false;
  cloud.pts[1].inSlice = false;
  cloud.pts[2].inSlice = true;
  std::vector<bool> isIce{true, true, false};
  std::vector<std::vector<int>> nList{{0}, {1}};
  Cloud ice;
  std::vector<int> linked;
  std::vector<int> counts;
  std::unordered_map<int, int> indexNumber;
  {
    CerrCapture quiet;
    REQUIRE(clump::largestIceCluster(dir.string() + "/", cloud, ice, nList, isIce,
                                     linked, counts, indexNumber, 1) == 0);
    (void)quiet;
  }
  REQUIRE(ice.nop == 0);
  REQUIRE(counts.empty());

  std::vector<bool> shortMask{true};
  cloud.pts[0].inSlice = false;
  cloud.pts[1].inSlice = true;
  isIce = {true, true, true};
  nList = {{0}, {1}, {2}};
  REQUIRE(clump::largestIceCluster(dir.string() + "/", cloud, ice, nList, isIce,
                                   linked, counts, indexNumber, 2) == 0);
  REQUIRE(ice.nop == 0);
  const auto domain = clump::largestDomain(cloud, nList, shortMask);
  REQUIRE(domain.subset == 0);
  REQUIRE_THAT(domain.percolation, Catch::Matchers::WithinAbs(0.0, 1e-15));
  cloud.pts[0].inSlice = true;
  cloud.pts[1].inSlice = true;
  cloud.pts[2].inSlice = true;
  std::vector<bool> hole{true, false, true};
  std::vector<std::vector<int>> linkedRows{{0, 1, 2}, {1, 0}, {2, 0}};
  const auto holed = clump::largestDomain(cloud, linkedRows, hole);
  REQUIRE(holed.subset == 2);
  REQUIRE(holed.largest == 2);
  const auto clipped = clump::largestDomain(cloud, linkedRows, std::vector<bool>{true});
  REQUIRE(clipped.subset == 1);
  REQUIRE(clipped.largest == 1);

  Cloud one = makeCloud({{1.0, 1.0, 1.0}}, {}, {8.0, 8.0, 8.0});
  std::vector<std::vector<int>> selfRow{{0}};
  REQUIRE(clump::recenterClusterCloud(one, selfRow) == 0);
  REQUIRE_THAT(one.pts[0].x, Catch::Matchers::WithinAbs(4.0, 1e-9));

  Cloud bare = makeCloud({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {}, {10.0, 10.0, 10.0});
  std::vector<bool> bothIce{true, true};
  std::vector<std::vector<int>> noRowsIce;
  Cloud iceBare;
  REQUIRE(clump::largestIceCluster(dir.string() + "/", bare, iceBare, noRowsIce, bothIce,
                                   linked, counts, indexNumber, 3) == 0);
  REQUIRE(iceBare.nop == 0);

  Cloud empty;
  empty.box = {8.0, 8.0, 8.0};
  empty.boxLow = {0.0, 0.0, 0.0};
  std::vector<std::vector<int>> noRows;
  CerrCapture none;
  REQUIRE(clump::recenterClusterCloud(empty, noRows) == 1);
  REQUIRE(none.str().find("no particles") != std::string::npos);

  Cloud crystal;
  crystal = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, crystal, 1);
  REQUIRE(crystal.nop > 4);
  crystal.box[0] += 30.0;
  molSys::Point<double> stray = crystal.pts[0];
  stray.atomID = 900000;
  stray.x = crystal.box[0] - 5.0;
  stray.y = 0.1;
  stray.z = 0.1;
  stray.inSlice = true;
  crystal.idIndexMap[stray.atomID] = crystal.nop;
  crystal.pts.push_back(stray);
  crystal.nop += 1;
  auto idList = nneigh::neighListO(3.5, crystal, 1);
  Cloud iceOut;
  std::vector<std::vector<int>> iceRows;
  REQUIRE(clump::clusterAnalysis(dir.string() + "/", iceOut, crystal, idList,
                                 iceRows, 3.5, 1, "q6") == 0);
  REQUIRE(iceOut.nop > 0);
  REQUIRE(iceOut.nop < crystal.nop);
  Cloud iceChill;
  std::vector<std::vector<int>> chillRows;
  REQUIRE(clump::clusterAnalysis(dir.string() + "/", iceChill, crystal, idList,
                                 chillRows, 3.5, 1, "chill") == 0);
  int water = 0;
  int iceLike = 0;
  for (const auto &pt : crystal.pts) {
    if (pt.iceType == molSys::atom_state_type::water) {
      water++;
    } else if (pt.iceType != molSys::atom_state_type::unclassified) {
      iceLike++;
    }
  }
  REQUIRE(water > 0);
  REQUIRE(iceLike > 0);
  std::filesystem::remove_all(dir);
}

TEST_CASE("masked ring update sees a one-sided new neighbour", "[cover_graph]") {
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(1);
#endif
  const int n = 24;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  auto link = [&](int i, int j) {
    rows[static_cast<std::size_t>(i)].push_back(j);
    rows[static_cast<std::size_t>(j)].push_back(i);
  };
  for (int i = 0; i < 6; i++) {
    link(i, (i + 1) % 6);
  }
  for (int i = 5; i < n - 1; i++) {
    link(i, i + 1);
  }
  for (int i = 0; i < n; i++) {
    auto &row = rows[static_cast<std::size_t>(i)];
    row.insert(row.begin(), i);
    std::sort(row.begin() + 1, row.end());
  }

  primitive::RingUpdater updater(8);
  const auto &first = updater.update(rows);
  const auto fresh = primitive::ringNetwork(rows, 8);
  REQUIRE(first == fresh);
  REQUIRE(updater.lastRecomputedSources() == n);
  const auto &same = updater.update(rows);
  REQUIRE(same == fresh);
  REQUIRE(updater.lastChangedBonds() == 0);
  REQUIRE(updater.lastRecomputedSources() == 0);

#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(2);
#endif
  auto edited = rows;
  edited[12].push_back(0);
  edited[0].push_back(20);
  edited[3].push_back(21);
  for (int v : {0, 3, 12}) {
    std::sort(edited[static_cast<std::size_t>(v)].begin() + 1,
              edited[static_cast<std::size_t>(v)].end());
  }
  const auto &masked = updater.update(edited);
  const auto editedFresh = primitive::ringNetwork(edited, 8);
  REQUIRE(masked == editedFresh);
  // One changed row is one directed endpoint. changedBonds is that
  // count divided by two, so it stays 0. The masked walk still runs.
  REQUIRE(updater.lastRecomputedSources() > 0);
  REQUIRE(updater.lastRecomputedSources() < n);
  REQUIRE(updater.lastBallsRefreshed() > 0);

  std::vector<std::vector<int>> dirty = rows;
  dirty[2].push_back(-1);
  dirty[2].push_back(n + 5);
  dirty[2].push_back(3);
  const auto dirtyRings = primitive::ringNetwork(dirty, 8);
  REQUIRE_FALSE(dirtyRings.empty());
  primitive::RingUpdater dirtyUp(8);
  REQUIRE(dirtyUp.update(dirty) == dirtyRings);
  auto dirtyNext = dirty;
  dirtyNext[4].push_back(0);
  std::sort(dirtyNext[4].begin() + 1, dirtyNext[4].end());
  const auto dirtySecond = dirtyUp.update(dirtyNext);
  REQUIRE(dirtySecond == primitive::ringNetwork(dirtyNext, 8));

  std::vector<std::vector<int>> isolated(8192);
  for (int i = 0; i < 8192; i++) {
    isolated[static_cast<std::size_t>(i)] = {i};
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(2);
#endif
  const auto isolatedRings = primitive::ringNetwork(isolated, 6);
  REQUIRE(isolatedRings.empty());
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(1);
#endif

  primitive::RingUpdater shallow(2);
  REQUIRE(shallow.update(rows).empty());
  REQUIRE(shallow.lastRecomputedSources() == 0);
  auto bothSides = edited;
  bothSides[0].push_back(12);
  std::sort(bothSides[0].begin() + 1, bothSides[0].end());
  primitive::RingUpdater full(8);
  full.setBondCrossover(0.0);
  full.update(rows);
  const auto &rebuilt = full.update(bothSides);
  REQUIRE(rebuilt == primitive::ringNetwork(bothSides, 8));
  REQUIRE(full.lastRecomputedSources() == n);
  REQUIRE(full.lastChangedBonds() > 0);
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(1);
#endif
}
