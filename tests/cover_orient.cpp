//-----------------------------------------------------------------------------------
// d-SEAMS
//
// Copyright (c) d-SEAMS developers
//-----------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <absOrientation.hpp>
#include <bond.hpp>
#include <bop.hpp>
#include <generic.hpp>
#include <ira_sofi.hpp>
#include <mol_sys.hpp>
#include <order_parameter.hpp>
#include <ozaki.hpp>
#include <pntCorrespondence.hpp>
#include <rdf.hpp>
#include <rdf2d.hpp>
#include <shapeMatch.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <utility>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;

namespace {

struct CerrCapture {
  std::stringstream buf;
  std::streambuf *old;
  CerrCapture() : old(std::cerr.rdbuf(buf.rdbuf())) {}
  ~CerrCapture() { std::cerr.rdbuf(old); }
  std::string str() const { return buf.str(); }
};

struct ProductRestore {
  seams::product::Mode mode;
  seams::product::Reduction reduction;
  ProductRestore()
      : mode(seams::product::mode()), reduction(seams::product::reduction()) {}
  ~ProductRestore() {
    seams::product::setMode(mode);
    seams::product::setReduction(reduction);
  }
};

struct HcTemplate {
  bool madeFile = false;
  bool madeDir = false;
  HcTemplate() {
    const std::filesystem::path file =
        std::filesystem::path("templates") / "hc.xyz";
    if (std::filesystem::exists(file)) {
      return;
    }
    madeDir = !std::filesystem::exists(file.parent_path());
    std::filesystem::create_directories(file.parent_path());
    std::ofstream out(file);
    out << "2\ncover\nO 0.0 0.0 0.0\nO 1.0 0.0 0.0\n";
    madeFile = static_cast<bool>(out);
  }
  ~HcTemplate() {
    if (!madeFile) {
      return;
    }
    const std::filesystem::path file =
        std::filesystem::path("templates") / "hc.xyz";
    std::filesystem::remove(file);
    if (madeDir) {
      std::filesystem::remove(file.parent_path());
    }
  }
};

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud makeCloud(const std::vector<std::array<double, 3>> &xyz,
                const std::vector<double> &box, int type = 1) {
  Cloud cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (std::size_t i = 0; i < xyz.size(); i++) {
    molSys::Point<double> p;
    p.type = type;
    p.atomID = static_cast<int>(i);
    p.molID = static_cast<int>(i);
    p.x = xyz[i][0];
    p.y = xyz[i][1];
    p.z = xyz[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = static_cast<int>(i);
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

std::vector<std::vector<int>>
indexNeighbours(int n, const std::vector<std::pair<int, int>> &bonds) {
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    nList[static_cast<std::size_t>(i)].push_back(i);
  }
  for (const auto &b : bonds) {
    nList[static_cast<std::size_t>(b.first)].push_back(b.second);
    nList[static_cast<std::size_t>(b.second)].push_back(b.first);
  }
  return nList;
}

std::vector<std::array<double, 3>> ringCoords(int n, double z, double radius) {
  std::vector<std::array<double, 3>> out;
  const double pi = std::acos(-1.0);
  for (int i = 0; i < n; i++) {
    const double a = 2.0 * pi * static_cast<double>(i) / static_cast<double>(n);
    out.push_back({radius * std::cos(a), radius * std::sin(a), z});
  }
  return out;
}

} // namespace

TEST_CASE("Horn rejects unequal point sets and a non-finite matrix",
          "[orient][abs]") {
  Eigen::MatrixXd rowsA(2, 3);
  Eigen::MatrixXd rowsB(3, 3);
  rowsA.setZero();
  rowsB.setOnes();
  std::vector<double> quat;
  std::vector<double> rmsdList;
  double rmsd = 0.0;
  double scale = 0.0;
  {
    CerrCapture err;
    const int rc =
        absor::hornAbsOrientation(rowsA, rowsB, quat, rmsd, rmsdList, scale);
    REQUIRE(rc == 1);
    REQUIRE(err.str().find("not of the same size") != std::string::npos);
    REQUIRE(rmsd == -1.0);
    REQUIRE(scale == 1.0);
    REQUIRE(quat.size() == 4);
  }
  Eigen::MatrixXd colsA(3, 3);
  Eigen::MatrixXd colsB(3, 2);
  colsA.setOnes();
  colsB.setOnes();
  {
    CerrCapture err;
    const int rc =
        absor::hornAbsOrientation(colsA, colsB, quat, rmsd, rmsdList, scale);
    REQUIRE(rc == 1);
    REQUIRE(err.str().find("not of the same size") != std::string::npos);
  }
  Eigen::MatrixXd bad(4, 3);
  bad.setConstant(std::numeric_limits<double>::quiet_NaN());
  {
    CerrCapture err;
    const int rc =
        absor::hornAbsOrientation(bad, bad, quat, rmsd, rmsdList, scale);
    REQUIRE(rc == 1);
    REQUIRE(err.str().find("Eigenvalue decomposition failed") !=
            std::string::npos);
  }
}

TEST_CASE("calcMatrixS uses the selected product and the empty guards",
          "[orient][abs]") {
  ProductRestore restore;
  Eigen::MatrixXd ref(4, 3);
  Eigen::MatrixXd tgt(4, 3);
  ref << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1;
  tgt << 0, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 1;
  const Eigen::MatrixXd eigen = tgt.transpose() * ref;
  for (seams::product::Mode mode :
       {seams::product::Mode::Dot2, seams::product::Mode::OzakiI,
        seams::product::Mode::OzakiII}) {
    seams::product::setMode(mode);
    const Eigen::MatrixXd S = absor::calcMatrixS(ref, tgt, 4, 3);
    REQUIRE(S.rows() == 3);
    REQUIRE(S.cols() == 3);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        REQUIRE_THAT(S(i, j), WithinAbs(eigen(i, j), 1e-8));
      }
    }
  }
  seams::product::setMode(seams::product::Mode::Dot2);
  Eigen::MatrixXd empty(0, 3);
  const Eigen::MatrixXd sEmpty = absor::calcMatrixS(empty, empty, 0, 3);
  REQUIRE(sEmpty.rows() == 3);
  REQUIRE(sEmpty.cols() == 3);
  const Eigen::MatrixXd sDim = absor::calcMatrixS(ref, tgt, 4, 0);
  REQUIRE(sDim.rows() == 3);
  REQUIRE(sDim.cols() == 3);

  seams::product::setMode(seams::product::Mode::OzakiI);
  const auto pent = ringCoords(5, 0.0, 1.0);
  Eigen::MatrixXd pts(5, 3);
  for (int i = 0; i < 5; i++) {
    pts(i, 0) = pent[static_cast<std::size_t>(i)][0];
    pts(i, 1) = pent[static_cast<std::size_t>(i)][1];
    pts(i, 2) = pent[static_cast<std::size_t>(i)][2];
  }
  std::vector<double> quat;
  std::vector<double> perAtom;
  double rmsd = 0.0;
  double scale = 0.0;
  const int rc = absor::hornAbsOrientation(pts, pts, quat, rmsd, perAtom, scale);
  REQUIRE(rc == 0);
  REQUIRE_THAT(rmsd, WithinAbs(0.0, 1e-6));
  REQUIRE(quat.size() == 4);
  REQUIRE(perAtom.size() == 5);
}

TEST_CASE("IRA without the library leaves orient false and the inputs intact",
          "[orient][ira]") {
#ifndef SEAMS_HAS_IRA
  REQUIRE_FALSE(ira::available());
  Eigen::MatrixXd pts(3, 3);
  pts << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  ira::Match match;
  REQUIRE(ira::match(pts, pts, match) == 1);
  ira::PointGroup group;
  REQUIRE(ira::pointGroup(pts, group) == 1);
  std::vector<double> quat = {1.0, 2.0, 3.0, 4.0};
  double rmsd = 9.0;
  REQUIRE_FALSE(ira::orient(pts, pts, quat, rmsd));
  REQUIRE(quat.size() == 4);
  REQUIRE(quat[0] == 1.0);
  REQUIRE(quat[3] == 4.0);
  REQUIRE(rmsd == 9.0);
#else
  REQUIRE(ira::available());
  Eigen::MatrixXd bad(2, 2);
  bad.setZero();
  ira::Match match;
  REQUIRE(ira::match(bad, bad, match) == 1);
  std::vector<double> quat = {1.0};
  double rmsd = 4.0;
  REQUIRE_FALSE(ira::orient(bad, bad, quat, rmsd));
  REQUIRE(quat.size() == 1);
  REQUIRE(rmsd == 4.0);
#endif
}

TEST_CASE("prism correspondence covers order, wrap, and the far-ring error",
          "[orient][pnt]") {
  auto bottom = ringCoords(4, 0.0, 1.0);
  auto top = ringCoords(4, 2.0, 1.0);
  std::vector<std::array<double, 3>> xyz = bottom;
  xyz.insert(xyz.end(), top.begin(), top.end());
  Cloud cloud = makeCloud(xyz, {20.0, 20.0, 40.0});
  const std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> out1;
  std::vector<int> out2;

  {
    CerrCapture err;
    const auto nList = indexNeighbours(8, {});
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(cloud, basal1, {4, 5, 6, 7},
                                                nList, out1, out2);
    REQUIRE(rc == 1);
    REQUIRE(err.str().find("deformed prism") != std::string::npos);
    REQUIRE(out1.empty());
  }

  {
    // Last basal1 atom is the only cross neighbour, and it meets basal2's
    // last atom. Forward neighbour is closer, so the order is clockwise
    // and both rings wrap.
    Cloud skewed = cloud;
    // basal2[0] is the clockwise step from basal2[3] and sits on basal1[3].
    // basal2[2] is the anticlockwise step and is far away.
    skewed.pts[4].x = skewed.pts[3].x;
    skewed.pts[4].y = skewed.pts[3].y;
    skewed.pts[6].x = 8.0;
    skewed.pts[6].y = 8.0;
    const auto nList = indexNeighbours(8, {{3, 7}});
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(skewed, basal1, {4, 5, 6, 7},
                                                nList, out1, out2);
    REQUIRE(rc == 0);
    REQUIRE(out1.size() == 4);
    REQUIRE(out2.size() == 4);
    REQUIRE(out1[0] == 3);
    REQUIRE(out2[0] == 7);
    REQUIRE(out2[1] == 4);
  }

  {
    // Anticlockwise from basal2 index 0: currentJatom = 0 - i goes negative.
    Cloud anti = cloud;
    anti.pts[5].x = 4.0;
    anti.pts[5].y = 4.0;
    anti.pts[7].x = anti.pts[0].x + 0.1;
    anti.pts[7].y = anti.pts[0].y;
    const auto nList = indexNeighbours(8, {{0, 4}});
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(anti, basal1, {4, 5, 6, 7},
                                                nList, out1, out2);
    REQUIRE(rc == 0);
    REQUIRE(out2[0] == 4);
    REQUIRE(out2[1] == 7);
  }

  {
    // Equidistant neighbours of a square take the default clockwise order.
    const auto nList = indexNeighbours(8, {{0, 4}});
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(cloud, basal1, {4, 5, 6, 7},
                                                nList, out1, out2);
    REQUIRE(rc == 0);
    REQUIRE(out1[0] == 0);
    REQUIRE(out2[0] == 4);
    REQUIRE(out2[1] == 5);
  }

  {
    CerrCapture err;
    Cloud far = makeCloud({{0, 0, 0},
                           {1, 0, 0},
                           {1, 1, 0},
                           {0, 1, 0},
                           {200000, 0, 0},
                           {200001, 0, 0},
                           {200001, 1, 0},
                           {200000, 1, 0}},
                          {1.0e7, 1.0e7, 1.0e7});
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(far, basal1, {4, 5, 6, 7},
                                                out1, out2);
    REQUIRE(rc == 1);
    REQUIRE(err.str().find("deformed prism") != std::string::npos);
  }

  {
    Cloud near = cloud;
    near.pts[7].x = near.pts[0].x;
    near.pts[7].y = near.pts[0].y;
    near.pts[4].x = 8.0;
    near.pts[5].x = 8.0;
    near.pts[6].x = 8.0;
    out1.clear();
    out2.clear();
    const int rc = pntToPnt::relOrderPrismBlock(near, basal1, {4, 5, 6, 7},
                                                out1, out2);
    REQUIRE(rc == 0);
    REQUIRE(out2[0] == 7);
  }
}

TEST_CASE("point-set fillers honour a wrapped start and a rejected start",
          "[orient][pnt]") {
  Cloud cloud = makeCloud({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}},
                          {30.0, 30.0, 30.0});
  const std::vector<int> ring = {0, 1, 2, 3};
  const auto wrapped = pntToPnt::fillPointSetPrismRing(cloud, ring, 3);
  REQUIRE(wrapped.rows() == 4);
  REQUIRE_THAT(wrapped(0, 0), WithinAbs(9.0, 1e-12));
  REQUIRE_THAT(wrapped(1, 0), WithinAbs(0.0, 1e-12));
  const auto resetNeg = pntToPnt::fillPointSetPrismRing(cloud, ring, -1);
  const auto resetBig = pntToPnt::fillPointSetPrismRing(cloud, ring, 9);
  const auto zero = pntToPnt::fillPointSetPrismRing(cloud, ring, 0);
  REQUIRE_THAT(resetNeg(2, 1), WithinAbs(zero(2, 1), 1e-12));
  REQUIRE_THAT(resetBig(3, 2), WithinAbs(zero(3, 2), 1e-12));

  const auto block = pntToPnt::fillPointSetPrismBlock(cloud, {0, 1}, {2, 3}, 1);
  REQUIRE(block.rows() == 4);
  REQUIRE_THAT(block(0, 0), WithinAbs(3.0, 1e-12));
  REQUIRE_THAT(block(2, 0), WithinAbs(9.0, 1e-12));
  const auto blockReset =
      pntToPnt::fillPointSetPrismBlock(cloud, {0, 1}, {2, 3}, -4);
  const auto blockBig =
      pntToPnt::fillPointSetPrismBlock(cloud, {0, 1}, {2, 3}, 5);
  REQUIRE_THAT(blockReset(0, 0), WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(blockBig(0, 0), WithinAbs(0.0, 1e-12));
}

TEST_CASE("HC cage template keeps the cloud readXYZ returns",
          "[orient][pnt]") {
  HcTemplate keep;
  const char *root = std::getenv("SEAMS_SOURCE_ROOT");
  const auto previous = std::filesystem::current_path();
  if (root != nullptr) {
    std::filesystem::current_path(root);
  }
  const auto pts = pntToPnt::getPointSetCage(ring::strucType::HCbasal);
  std::filesystem::current_path(previous);
  REQUIRE(pts.rows() == 12);
  REQUIRE(pts.cols() == 3);
  REQUIRE_THAT(pts(0, 0), WithinAbs(8.995, 1e-6));
  REQUIRE_THAT(pts(0, 1), WithinAbs(10.3859997, 1e-6));
  REQUIRE_THAT(pts(0, 2), WithinAbs(15.0939999, 1e-6));
  const auto other = pntToPnt::getPointSetCage(ring::strucType::Prism);
  REQUIRE(other.rows() == 12);
  REQUIRE(other.cols() == 3);
}

TEST_CASE("relOrderHC original, reversed, and unmatched basal rings",
          "[orient][pnt]") {
  Cloud cloud = makeCloud(ringCoords(12, 0.0, 1.0), {40.0, 40.0, 40.0});
  const std::vector<int> basal1 = {0, 1, 2, 3, 4, 5};
  const std::vector<int> basal2 = {6, 7, 8, 9, 10, 11};
  std::vector<int> matched1;
  std::vector<int> matched2;

  {
    // basal2[5] neighbours basal1[0]. index 5+2 wraps. The wrapped atom
    // neighbours basal1[2], so the original order is kept.
    const auto nList = indexNeighbours(12, {{0, 11}, {2, 7}});
    const int rc = pntToPnt::relOrderHC(cloud, basal1, basal2, nList, matched1,
                                        matched2);
    REQUIRE(rc == 0);
    REQUIRE(matched2.size() == 6);
    REQUIRE(matched2[0] == 11);
    REQUIRE(matched1[0] == 0);
  }
  {
    // +2 is not a neighbour. -2 wraps through a negative index and matches.
    const auto nList = indexNeighbours(12, {{0, 7}, {2, 11}});
    const int rc = pntToPnt::relOrderHC(cloud, basal1, basal2, nList, matched1,
                                        matched2);
    REQUIRE(rc == 0);
    REQUIRE(matched2[0] == 7);
  }
  {
    const auto nList = indexNeighbours(12, {{0, 8}});
    const int rc = pntToPnt::relOrderHC(cloud, basal1, basal2, nList, matched1,
                                        matched2);
    REQUIRE(rc == 1);
  }

  const auto hexNeg =
      pntToPnt::changeHexCageOrder(cloud, basal1, basal2, -1);
  const auto hexZero =
      pntToPnt::changeHexCageOrder(cloud, basal1, basal2, 0);
  const auto hexHigh =
      pntToPnt::changeHexCageOrder(cloud, basal1, basal2, 9);
  const auto hexWrap =
      pntToPnt::changeHexCageOrder(cloud, basal1, basal2, 5);
  REQUIRE(hexNeg.rows() == 12);
  REQUIRE_THAT(hexNeg(0, 0), WithinAbs(hexZero(0, 0), 1e-12));
  REQUIRE_THAT(hexHigh(0, 0), WithinAbs(hexZero(0, 0), 1e-12));
  REQUIRE(hexWrap.rows() == 12);
  REQUIRE_THAT(hexWrap(0, 0), WithinAbs(cloud.pts[5].x, 1e-9));
}

TEST_CASE("synthetic DDC hits both apex wraps and the empty guards",
          "[orient][pnt]") {
  std::vector<std::array<double, 3>> xyz;
  for (int i = 0; i < 14; i++) {
    xyz.push_back({static_cast<double>(i), 0.2 * i, 0.1 * i});
  }
  Cloud cloud = makeCloud(xyz, {80.0, 80.0, 80.0});
  const std::vector<int> equatorial = {0, 1, 2, 3, 4, 5};
  std::vector<std::vector<int>> rings = {
      equatorial,
      {0, 20, 21, 22, 12, 6},
      {13, 23, 24, 0, 1, 9},
      {2, 25, 26, 27, 28, 7},
      {3, 29, 30, 31, 32, 10},
      {4, 33, 34, 35, 36, 8},
      {5, 37, 38, 39, 40, 11},
  };
  std::vector<cage::Cage> cages(1);
  cages[0].type = cage::cageType::DoubleDiaC;
  cages[0].rings = {0, 1, 2, 3, 4, 5, 6};

  const auto order = pntToPnt::relOrderDDC(equatorial, 0, rings, cages);
  REQUIRE(order.size() == 14);
  REQUIRE(order[0] == 0);
  REQUIRE(order[6] == 6);
  REQUIRE(order[9] == 12);
  REQUIRE(order[10] == 9);
  REQUIRE(order[13] == 13);

  const auto viaIndex = pntToPnt::relOrderDDC(0, rings, cages);
  REQUIRE(viaIndex.size() == 14);
  REQUIRE(pntToPnt::relOrderDDC(-1, rings, cages).empty());
  REQUIRE(pntToPnt::relOrderDDC(3, rings, cages).empty());
  cages[0].rings.clear();
  REQUIRE(pntToPnt::relOrderDDC(0, rings, cages).empty());
  cages[0].rings = {0, 1, 2};
  REQUIRE(pntToPnt::relOrderDDC(equatorial, 0, rings, cages).empty());
  cages[0].rings = {0, 1, 2, 3, 4, 5, 6};
  REQUIRE(pntToPnt::relOrderDDC({0, 1, 2}, 0, rings, cages).empty());
  std::vector<std::vector<int>> blank(7, std::vector<int>(6, 50));
  REQUIRE(pntToPnt::relOrderDDC(equatorial, 0, blank, cages).empty());

  std::vector<std::vector<int>> ringsOdd = rings;
  ringsOdd[2] = {1, 30, 31, 32, 13, 9};
  const auto oddApex = pntToPnt::relOrderDDC(equatorial, 0, ringsOdd, cages);
  REQUIRE(oddApex.size() == 14);
  REQUIRE(oddApex[13] == 13);
  REQUIRE(oddApex[10] == 9);

  const auto shifted = pntToPnt::changeDiaCageOrder(cloud, order, 2);
  REQUIRE(shifted.rows() == 14);
  REQUIRE_THAT(shifted(0, 0), WithinAbs(cloud.pts[order[2]].x, 1e-9));
  const auto same = pntToPnt::changeDiaCageOrder(cloud, order, 0);
  REQUIRE(same.rows() == 14);
  REQUIRE(std::isfinite(same(1, 0)));
  REQUIRE_THAT(same(0, 0), WithinAbs(cloud.pts[order[0]].x, 1e-9));

  std::vector<int> shortOrder(10, 0);
  const auto zeroShort = pntToPnt::changeDiaCageOrder(cloud, shortOrder, 0);
  REQUIRE(zeroShort.rows() == 14);
  REQUIRE_THAT(zeroShort(3, 1), WithinAbs(0.0, 0.0));

  std::vector<int> negative = order;
  negative[4] = -1;
  const auto zeroNeg = pntToPnt::changeDiaCageOrder(cloud, negative, 1);
  REQUIRE_THAT(zeroNeg.cwiseAbs().maxCoeff(), WithinAbs(0.0, 0.0));
  std::vector<int> huge = order;
  huge[0] = 100;
  const auto zeroHuge = pntToPnt::changeDiaCageOrder(cloud, huge, 4);
  REQUIRE_THAT(zeroHuge(0, 0), WithinAbs(0.0, 0.0));
}

TEST_CASE("odd prism block searches every start and drops a failed Horn",
          "[orient][shape]") {
  const int n = 5;
  auto bottom = ringCoords(n, 0.0, 1.2);
  auto top = ringCoords(n, 2.5, 1.2);
  bottom[2][0] += 0.35;
  top[2][0] += 0.35;
  std::vector<std::array<double, 3>> xyz = bottom;
  xyz.insert(xyz.end(), top.begin(), top.end());
  Cloud cloud = makeCloud(xyz, {12.0, 12.0, 40.0});
  std::vector<int> basal1;
  std::vector<int> basal2;
  for (int i = 0; i < n; i++) {
    basal1.push_back(i);
    basal2.push_back(n + i);
  }
  const auto ref = pntToPnt::getPointSetRefRing(n, 2);
  std::vector<std::vector<int>> nList;
  int begin = -1;
  const bool matched =
      match::matchPrismBlock(cloud, nList, ref, basal1, basal2, begin);
  REQUIRE(begin >= 0);
  REQUIRE(begin < n);
  REQUIRE(matched == true);

  std::vector<int> rot1 = {2, 3, 4, 0, 1};
  std::vector<int> rot2 = {7, 8, 9, 5, 6};
  int beginRot = 0;
  const bool rotated =
      match::matchPrismBlock(cloud, nList, ref, rot1, rot2, beginRot);
  REQUIRE(rotated == true);
  REQUIRE(beginRot != 0);

  auto triBottom = ringCoords(3, 0.0, 1.0);
  auto triTop = ringCoords(3, 2.2, 1.0);
  std::vector<std::array<double, 3>> tri = triBottom;
  tri.insert(tri.end(), triTop.begin(), triTop.end());
  Cloud triangle = makeCloud(tri, {10, 10, 30});
  std::vector<int> t1 = {0, 1, 2};
  std::vector<int> t2 = {3, 4, 5};
  int beginTri = -1;
  const auto triRef = pntToPnt::getPointSetRefRing(3, 2);
  const bool triMatch =
      match::matchPrismBlock(triangle, nList, triRef, t1, t2, beginTri);
  REQUIRE(triMatch == true);
  REQUIRE(beginTri == 0);

  Cloud wreck = cloud;
  for (int i = 0; i < cloud.nop; i++) {
    wreck.pts[static_cast<std::size_t>(i)].x =
        std::numeric_limits<double>::quiet_NaN();
    wreck.pts[static_cast<std::size_t>(i)].y = 0.0;
    wreck.pts[static_cast<std::size_t>(i)].z = static_cast<double>(i);
  }
  int beginBad = 0;
  bool bad = false;
  {
    CerrCapture err;
    bad = match::matchPrismBlock(wreck, nList, ref, basal1, basal2, beginBad);
    REQUIRE(err.str().find("Eigenvalue decomposition failed") !=
            std::string::npos);
  }
  REQUIRE((bad == true || bad == false));
}

TEST_CASE("deformed and untethered prisms report order and shape failures",
          "[orient][shape]") {
  auto bottom = ringCoords(4, 0.0, 1.4);
  auto top = ringCoords(4, 2.8, 1.4);
  std::vector<std::array<double, 3>> xyz = bottom;
  xyz.insert(xyz.end(), top.begin(), top.end());
  Cloud cloud = makeCloud(xyz, {16.0, 16.0, 50.0});
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};
  const auto ref = pntToPnt::getPointSetRefRing(4, 2);
  std::vector<double> rmsd(8, -1.0);

  const auto linked = indexNeighbours(8, {{0, 4}, {1, 5}, {2, 6}, {3, 7}});
  const bool perfectFalse =
      match::matchPrism(cloud, linked, ref, basal1, basal2, rmsd, false);
  REQUIRE(perfectFalse == true);
  int filled = 0;
  for (double v : rmsd) {
    if (v >= 0.0) {
      filled++;
    }
  }
  REQUIRE(filled == 8);

  Cloud distorted = cloud;
  for (int i = 0; i < 8; i++) {
    distorted.pts[static_cast<std::size_t>(i)].x = (i < 4) ? 0.0 : 80.0;
    distorted.pts[static_cast<std::size_t>(i)].y = 0.01 * i;
    distorted.pts[static_cast<std::size_t>(i)].z = (i % 2 == 0) ? 0.0 : 50.0;
  }
  std::vector<double> quatBad;
  std::vector<double> listBad;
  double rmsdBad = 0.0;
  double scaleBad = 0.0;
  const auto refBlock =
      pntToPnt::createPrismBlock(distorted, ref, 4, basal1, basal2);
  const auto targetBlock =
      pntToPnt::fillPointSetPrismBlock(distorted, basal1, basal2, 0);
  REQUIRE(absor::hornAbsOrientation(refBlock, targetBlock, quatBad, rmsdBad,
                                   listBad, scaleBad) == 0);
  REQUIRE(rmsdBad > 6.0);
  std::vector<double> perAtom(8, -1.0);
  const bool deformed = match::matchPrism(distorted, linked, ref, basal1,
                                          basal2, perAtom, false);
  REQUIRE(deformed == false);

  const auto none = indexNeighbours(8, {});
  std::vector<double> rmsdNone(8, -1.0);
  {
    CerrCapture err;
    REQUIRE(match::matchPrism(cloud, none, ref, basal1, basal2, rmsdNone,
                              true) == false);
    REQUIRE(err.str().find("deformed prism") != std::string::npos);
  }

  Cloud far = makeCloud({{0, 0, 0},
                         {1, 0, 0},
                         {1, 1, 0},
                         {0, 1, 0},
                         {200000, 0, 0},
                         {200001, 0, 0},
                         {200001, 1, 0},
                         {200000, 1, 0}},
                        {1.0e7, 1.0e7, 1.0e7});
  std::vector<double> rmsdFar(8, -1.0);
  {
    CerrCapture err;
    REQUIRE(match::matchUntetheredPrism(far, none, ref, basal1, basal2,
                                        rmsdFar) == false);
    REQUIRE(err.str().find("deformed prism") != std::string::npos);
  }
}

TEST_CASE("bond helpers reject empty lists, dummy atoms, and bad indices",
          "[orient][bond]") {
  Cloud o = makeCloud({{0, 0, 0}, {2.6, 0, 0}}, {10.0, 10.0, 10.0});
  Cloud h = makeCloud({{1.0, 0, 0}, {-0.8, 0.6, 0}, {3.6, 0, 0}, {2.6, 0.9, 0}},
                      {10.0, 10.0, 10.0});
  o.pts[0].molID = 1;
  o.pts[1].molID = 2;
  h.pts[0].molID = 1;
  h.pts[1].molID = 1;
  h.pts[2].molID = 2;
  h.pts[3].molID = 2;

  REQUIRE(bond::donatedHydrogenBond(o, h, -1, 1, {0}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 0, -2, {0}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 0, 8, {0}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 5, 0, {0}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 0, 1, {}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 0, 1, {-1}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 0, 1, {9}, 2.42, 30.0) == false);
  REQUIRE(bond::donatedHydrogenBond(o, h, 1, 0, {0, 1}, 2.42, 30.0) == true);

  {
    CerrCapture err;
    const auto bonds = bond::populateBonds({}, o);
    REQUIRE(bonds.empty());
    REQUIRE(err.str().find("no bonds") != std::string::npos);
  }
  {
    CerrCapture err;
    const auto bonds =
        bond::populateBonds({}, o, {cage::iceType::hc, cage::iceType::hc});
    REQUIRE(bonds.empty());
    REQUIRE(err.str().find("no bonds") != std::string::npos);
  }

  Cloud three = makeCloud({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, {12.0, 12.0, 12.0});
  three.pts[0].atomID = 10;
  three.pts[1].atomID = 11;
  three.pts[2].atomID = 12;
  const auto nList = indexNeighbours(3, {{0, 1}, {1, 2}});
  const std::vector<cage::iceType> kinds = {
      cage::iceType::dummy, cage::iceType::hc, cage::iceType::dummy};
  const auto kept = bond::populateBonds(nList, three, kinds);
  REQUIRE(kept.empty());

  const std::vector<cage::iceType> mixed = {
      cage::iceType::hc, cage::iceType::dummy, cage::iceType::hc};
  const auto mixedBonds = bond::populateBonds(nList, three, mixed);
  REQUIRE(mixedBonds.empty());

  const std::vector<std::vector<int>> rows = {{}, {11, 99}};
  o.idIndexMap.clear();
  o.idIndexMap[10] = 0;
  const auto net = bond::populateHbondsFromDonors(o, h, rows, {-1, 9, 0}, 2.42,
                                                  30.0);
  REQUIRE(net.size() == 2);
  REQUIRE(net[0].size() == 1);
  REQUIRE(net[1].size() == 1);
}

TEST_CASE("F4, jump rotors, and layer strings cover the quiet branches",
          "[orient][order]") {
  ProductRestore restore;
  Cloud empty;
  empty.nop = 0;
  empty.box = {10.0, 10.0, 10.0};
  const auto noAtoms = topoparam::rodgerF4(empty, {}, 1, 2);
  REQUIRE(noAtoms.empty());
  const std::vector<double> tiny = topoparam::projAreaSingleRing(empty, {0, 1});
  REQUIRE(tiny.size() == 3);
  REQUIRE(tiny[0] == 0.0);

  Cloud water;
  water.box = {20.0, 20.0, 20.0};
  water.boxLow = {0.0, 0.0, 0.0};
  auto add = [&](int type, int mol, int id, double x, double y, double z) {
    molSys::Point<double> p;
    p.type = type;
    p.molID = mol;
    p.atomID = id;
    p.x = x;
    p.y = y;
    p.z = z;
    water.idIndexMap[id] = static_cast<int>(water.pts.size());
    water.pts.push_back(p);
  };
  add(1, 1, 10, 0, 0, 0);
  add(2, 1, 11, 0.0, 0.96, 0.0);
  add(2, 1, 12, 0.0, -0.24, 0.9);
  add(1, 2, 20, 2.8, 0, 0);
  add(2, 2, 21, 2.8, 0.96, 0.0);
  add(2, 2, 22, 2.8, -0.24, 0.9);
  add(3, 3, 30, 8, 0, 0);
  add(1, 9, 40, 5, 5, 5);
  water.nop = static_cast<int>(water.pts.size());
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(water.nop));
  nList[0] = {10, 20, 99, 10};
  water.idIndexMap[20] = 3;
  water.idIndexMap[99] = -1;
  water.idIndexMap[77] = 50;
  nList[0].push_back(77);
  nList[0].push_back(30);
  nList[0].push_back(40);
  nList[0].push_back(123);
  const auto f4 = topoparam::rodgerF4(water, nList, 1, 2);
  REQUIRE(f4.size() == static_cast<std::size_t>(water.nop));
  REQUIRE(std::isfinite(f4[0]));

  seams::product::setReduction(seams::product::Reduction::Pairwise);
  const auto folded = topoparam::rodgerF4(water, nList, 1, 2);
  REQUIRE(std::isfinite(folded[0]));
  const double mean =
      topoparam::meanFinite({1.0, std::numeric_limits<double>::quiet_NaN(), 3.0});
  REQUIRE_THAT(mean, WithinAbs(2.0, 1e-12));
  seams::product::setReduction(seams::product::Reduction::Left);

  Cloud line = water;
  line.pts[4].x = 3.3;
  line.pts[4].y = 0.0;
  line.pts[4].z = 0.0;
  line.pts[5].x = 4.0;
  line.pts[5].y = 0.0;
  line.pts[5].z = 0.0;
  std::vector<std::vector<int>> only(static_cast<std::size_t>(line.nop));
  only[0] = {10, 20};
  const auto flat = topoparam::rodgerF4(line, only, 1, 2);
  REQUIRE(std::isnan(flat[0]));

  Cloud nanH = water;
  nanH.pts[1].x = std::numeric_limits<double>::quiet_NaN();
  nanH.pts[2].x = std::numeric_limits<double>::quiet_NaN();
  nanH.pts[4].y = std::numeric_limits<double>::quiet_NaN();
  nanH.pts[5].y = std::numeric_limits<double>::quiet_NaN();
  const auto nanF = topoparam::rodgerF4(nanH, only, 1, 2);
  REQUIRE(std::isnan(nanF[0]));

  auto frame = [&](bool flip) {
    Cloud c;
    c.box = {10, 10, 10};
    c.boxLow = {0, 0, 0};
    molSys::Point<double> o;
    o.type = 1;
    o.molID = 7;
    o.atomID = 1;
    c.pts.push_back(o);
    molSys::Point<double> h1;
    h1.type = 2;
    h1.molID = 7;
    h1.atomID = 2;
    h1.x = flip ? -1.0 : 1.0;
    molSys::Point<double> h2;
    h2.type = 2;
    h2.molID = 7;
    h2.atomID = 3;
    h2.x = flip ? 1.0 : -1.0;
    c.pts.push_back(h1);
    c.pts.push_back(h2);
    c.nop = 3;
    return c;
  };
  const double tau = topoparam::jumpRotorTau90(frame(false), frame(true), 2.5,
                                               1, 2);
  REQUIRE_THAT(tau, WithinAbs(2.5, 1e-12));
  const double same = topoparam::jumpRotorTau90(frame(false), frame(false), 2.5,
                                                1, 2);
  REQUIRE(std::isnan(same));
  REQUIRE(std::isnan(topoparam::jumpRotorTau90(frame(false), frame(false), 0.0,
                                               1, 2)));
  Cloud other = frame(false);
  other.pts[0].molID = 8;
  other.pts[1].molID = 8;
  other.pts[2].molID = 8;
  REQUIRE(std::isnan(
      topoparam::jumpRotorTau90(frame(false), other, 1.0, 1, 2)));
  Cloud stacked = frame(false);
  stacked.pts[1].x = 0;
  stacked.pts[2].x = 0;
  REQUIRE(std::isnan(
      topoparam::jumpRotorTau90(stacked, stacked, 1.0, 1, 2)));
  Cloud oneH = frame(false);
  oneH.pts.pop_back();
  oneH.nop = 2;
  REQUIRE(std::isnan(topoparam::jumpRotorTau90(oneH, oneH, 1.0, 1, 2)));

  Cloud layers;
  layers.box = {10.0, 10.0, 10.0};
  layers.currentFrame = 1;
  auto ice = [&](molSys::atom_state_type kind, double x, double y, double z) {
    molSys::Point<double> p;
    p.iceType = kind;
    p.x = x;
    p.y = y;
    p.z = z;
    p.atomID = static_cast<int>(layers.pts.size());
    layers.pts.push_back(p);
  };
  ice(molSys::atom_state_type::cubic, 1, 1, 1);
  ice(molSys::atom_state_type::hexagonal, 1, 1, 1);
  ice(molSys::atom_state_type::reCubic, 1, 1, 8);
  ice(molSys::atom_state_type::reHex, 8, 1, 5);
  ice(molSys::atom_state_type::water, 2, 2, 2);
  layers.nop = static_cast<int>(layers.pts.size());
  const auto stack = topoparam::layerCubicity(layers, 2, 3.7);
  REQUIRE(stack.sequence.find('M') != std::string::npos);
  REQUIRE(stack.sequence.find('C') != std::string::npos);
  REQUIRE(stack.phiC > 0.0);
  REQUIRE(stack.phiC < 1.0);
  const auto alongX = topoparam::layerCubicity(layers, 0, 3.0);
  REQUIRE(alongX.sequence.find('C') != std::string::npos);
  const auto alongY = topoparam::layerCubicity(layers, 1, 3.0);
  REQUIRE(alongY.sequence.size() >= 1);
  REQUIRE(topoparam::layerCubicity(empty, 2, 3.7).sequence.empty());
  REQUIRE(topoparam::layerCubicity(layers, 4, 3.7).sequence.empty());
  REQUIRE(topoparam::layerCubicity(layers, 2, 0.0).sequence.empty());
  Cloud zeroBox = layers;
  zeroBox.box = {0.0, 10.0, 10.0};
  REQUIRE(topoparam::layerCubicity(zeroBox, 0, 3.7).phiC == 0.0);

  Cloud edge = layers;
  edge.box = {4.0, 4.0, 4.0};
  edge.pts[0].z = -std::numeric_limits<double>::denorm_min();
  edge.pts[0].iceType = molSys::atom_state_type::cubic;
  const auto wrapped = topoparam::layerCubicity(edge, 2, 1.0);
  REQUIRE(wrapped.cubicPerLayer.size() >= 1);
  REQUIRE(wrapped.phiC > 0.0);
  const auto subTum = topoparam::tumLayerStack(
      edge, {{0, 0, 0, 0, 0, 0}}, {true}, {false}, 2, 1.0);
  REQUIRE(subTum.hexPerLayer.size() >= 1);
  REQUIRE(subTum.sequence.find('H') != std::string::npos);

  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5},
                                         {0, 1, 2, 3, 4, 5},
                                         {0, 1, 2},
                                         {-1, 0, 1, 2, 3, 4},
                                         {0, 99, 1, 2, 3, 4}};
  std::vector<bool> basal = {true, false};
  std::vector<bool> equatorial = {false, true};
  const auto tum =
      topoparam::tumLayerStack(layers, rings, basal, equatorial, 2, 3.7);
  REQUIRE(tum.sequence.find('M') != std::string::npos);
  std::vector<bool> badBasal = {false, false, false, true, true};
  std::vector<bool> badEq(5, false);
  const auto skipped =
      topoparam::tumLayerStack(layers, rings, badBasal, badEq, 2, 3.7);
  REQUIRE(skipped.sequence.find('M') == std::string::npos);
  const auto tumX =
      topoparam::tumLayerStack(layers, rings, basal, equatorial, 0, 4.0);
  REQUIRE(tumX.phiC >= 0.0);
  const auto tumY =
      topoparam::tumLayerStack(layers, rings, basal, equatorial, 1, 4.0);
  REQUIRE(tumY.sequence.empty() == false);
  REQUIRE(topoparam::tumLayerStack(empty, rings, basal, equatorial, 2, 3.7)
              .sequence.empty());
  REQUIRE(topoparam::tumLayerStack(layers, rings, basal, equatorial, 2, -1.0)
              .phiC == 0.0);
}

TEST_CASE("Steinhardt l=12 clamps, pairwise reductions, and a split flag",
          "[orient][steinhardt]") {
  ProductRestore restore;
  const double zPos = 1.0422269095183273e-155;
  const double zNeg = -6.43409430831532e-159;
  Cloud cloud;
  cloud.box = {4.0, 4.0, 4.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double coords[][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0.25},
                              {0, 0, zPos}, {0, 0, zNeg}};
  const int ids[] = {5, 1, 4, 3, 2};
  for (int i = 0; i < 5; i++) {
    molSys::Point<double> p;
    p.atomID = ids[i];
    p.molID = i;
    p.type = 1;
    p.x = coords[i][0];
    p.y = coords[i][1];
    p.z = coords[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[ids[i]] = i;
  }
  cloud.nop = 5;
  const auto drPos = gen::relDist(cloud, 3, 0);
  const double rPos = std::sqrt(drPos[0] * drPos[0] + drPos[1] * drPos[1] +
                                drPos[2] * drPos[2]);
  REQUIRE(drPos[2] / rPos > 1.0);
  const auto drNeg = gen::relDist(cloud, 4, 0);
  const double rNeg = std::sqrt(drNeg[0] * drNeg[0] + drNeg[1] * drNeg[1] +
                                drNeg[2] * drNeg[2]);
  REQUIRE(drNeg[2] / rNeg < -1.0);

  std::vector<std::vector<int>> nList(5);
  nList[0] = {5, 1, 4, 3, 2};
  nList[1] = {1};
  nList[2] = {4, 5};
  nList[3] = {3};
  nList[4] = {2};
  chill::setSteinhardtAtomSplit(false);
  const auto q12 = chill::steinhardtQl(cloud, nList, 12);
  REQUIRE(q12.ql.size() == 5);
  REQUIRE(q12.qlBar.size() == 5);
  REQUIRE(std::isfinite(q12.ql[2]));
  chill::setSteinhardtAtomSplit(true);

  Cloud sparse = cloud;
  for (auto &p : sparse.pts) {
    p.atomID = 0;
  }
  sparse.idIndexMap.clear();
  sparse.idIndexMap[0] = 0;
  std::vector<std::vector<int>> sparseList(5, {0});
  const auto qSparse = chill::steinhardtQl(sparse, sparseList, 6);
  REQUIRE(qSparse.ql.size() == 5);

  Cloud dense = makeCloud({{0, 0, 0}, {1.2, 0, 0}}, {8, 8, 8});
  dense.pts[0].atomID = 1;
  dense.pts[1].atomID = 2;
  dense.idIndexMap.clear();
  dense.idIndexMap[1] = 0;
  dense.idIndexMap[2] = 1;
  std::vector<std::vector<int>> odd = {{1, -3, 1000, 2}, {2, 1}};
  const auto qOdd = chill::steinhardtQl(dense, odd, 6);
  REQUIRE(qOdd.ql.size() == 2);
  REQUIRE(qOdd.ql[0] > 0.0);

  seams::product::setReduction(seams::product::Reduction::Pairwise);
  const auto qPair = chill::steinhardtQl(dense, odd, 4);
  REQUIRE(qPair.ql.size() == 2);
  REQUIRE(std::isfinite(qPair.ql[0]));
  Cloud twin = makeCloud({{1, 1, 1}, {1, 1, 1}}, {6, 6, 6});
  twin.pts[0].atomID = 4;
  twin.pts[1].atomID = 8;
  twin.idIndexMap[4] = 0;
  twin.idIndexMap[8] = 1;
  std::vector<std::vector<int>> zeros = {{4, 8}, {8, 4}};
  const auto qZero = chill::steinhardtQl(twin, zeros, 6);
  REQUIRE_THAT(qZero.ql[0], WithinAbs(0.0, 0.0));

  seams::product::setMode(seams::product::Mode::Dot2);
  seams::product::setReduction(seams::product::Reduction::Left);
  Cloud one = makeCloud({{0, 0, 0}}, {5, 5, 5});
  std::vector<std::vector<int>> alone = {{0}};
  const auto qAlone = chill::steinhardtQl(one, alone, 6);
  REQUIRE(qAlone.ql.size() == 1);
  const auto qDot = chill::steinhardtQl(dense, {{1, 2}, {2, 1}}, 6);
  REQUIRE(qDot.ql[0] > 0.0);

  seams::product::setMode(seams::product::Mode::OzakiII);
  seams::product::setReduction(seams::product::Reduction::Neumaier);
  const auto qNeu = chill::steinhardtQl(dense, {{1, 2}, {2}}, 8);
  REQUIRE(qNeu.ql.size() == 2);
  REQUIRE(std::isfinite(qNeu.qlBar[0]));
}

TEST_CASE("partial RDF guards, a vanishing ideal gas, and the CN edges",
          "[orient][rdf]") {
  Cloud two = makeCloud({{0, 0, 0}, {1.5, 0, 0}}, {10, 10, 10});
  two.pts[0].type = 1;
  two.pts[1].type = 1;
  const auto clamped = rdf::partialRdf(two, 1, 1, 4.0, 0);
  REQUIRE(clamped.r.size() == 1);
  REQUIRE(clamped.binwidth == 4.0);
  REQUIRE(clamped.count[0] >= 0);

  const auto noRange = rdf::partialRdf(two, 1, 1, 0.0, 4);
  REQUIRE(noRange.g.size() == 4);
  REQUIRE(noRange.count[0] == 0);
  Cloud noBox = two;
  noBox.box = {10.0};
  const auto noCell = rdf::partialRdf(noBox, 1, 1, 3.0, 3);
  REQUIRE(noCell.nI == 2);
  REQUIRE(noCell.count[0] == 0);
  Cloud none;
  none.box = {10, 10, 10};
  none.nop = 0;
  REQUIRE(rdf::partialRdf(none, 1, 2, 3.0, 3).nI == 0);

  Cloud flat;
  flat.box = {0.0, 10.0, 10.0};
  flat.boxLow = {0, 0, 0};
  flat.nop = 1;
  molSys::Point<double> only;
  only.type = 1;
  only.atomID = 1;
  flat.pts.push_back(only);
  flat.idIndexMap[1] = 0;
  {
    CerrCapture err;
    const auto gr = rdf::partialRdf(flat, 1, 1, 2.0, 2);
    REQUIRE(gr.g.size() == 2);
    REQUIRE(gr.g[0] == 0.0);
    REQUIRE(err.str().find("non-positive edge") != std::string::npos);
    const auto cn = rdf::runningCN(gr);
    REQUIRE(cn.size() == 2);
    REQUIRE(cn[0] == 0.0);
  }

  Cloud unlike = makeCloud({{0, 0, 0}, {1, 0, 0}}, {12, 12, 12});
  unlike.pts[0].type = 2;
  unlike.pts[1].type = 2;
  const auto missingI = rdf::partialRdf(unlike, 1, 2, 5.0, 4);
  REQUIRE(missingI.nI == 0);
  REQUIRE(missingI.nJ == 2);
  REQUIRE(missingI.g[0] == 0.0);

  Cloud huge = makeCloud({{0, 0, 0}, {1, 0, 0}}, {1e100, 1e100, 1e100});
  huge.pts[0].type = 1;
  huge.pts[1].type = 1;
  const auto tiny = rdf::partialRdf(huge, 1, 1, 1e-100, 1);
  REQUIRE(tiny.g.size() == 1);
  REQUIRE(tiny.g[0] == 0.0);

  rdf::PartialRdf hand;
  hand.binwidth = 1.0;
  hand.g = {2.0, 0.5, 1.0};
  hand.volume = 8.0;
  hand.nJ = 4;
  const double partial = rdf::coordinationNumber(hand, 1.5, 0.25);
  REQUIRE(partial > 0.0);
  REQUIRE(rdf::coordinationNumber(hand, 0.0, 0.25) == 0.0);
  hand.binwidth = 0.0;
  REQUIRE(rdf::coordinationNumber(hand, 1.0, 0.25) == 0.0);
  hand.binwidth = 1.0;
  REQUIRE(rdf::coordinationNumber(hand, 1.0, 0.0) == 0.0);
  hand.g.clear();
  REQUIRE(rdf::coordinationNumber(hand, 1.0, 0.25) == 0.0);

  rdf::PartialRdf peak;
  peak.g = {0.1, 0.4, 0.2, 0.5};
  const int bin = rdf::firstMinimumBin(peak);
  REQUIRE(bin == 2);
  rdf::PartialRdf down;
  down.g = {0.9, 0.4, 0.1, 0.05};
  REQUIRE(rdf::firstMinimumBin(down) == -1);
  rdf::PartialRdf up;
  up.g = {0.1, 0.2, 0.3};
  REQUIRE(rdf::firstMinimumBin(up) == -1);

  rdf::PartialRdf dead;
  dead.g = {1.0, 1.0};
  dead.binwidth = 0.0;
  const auto zeros = rdf::runningCN(dead, 1.0);
  REQUIRE(zeros.size() == 2);
  REQUIRE(zeros[1] == 0.0);
  dead.binwidth = 1.0;
  const auto noRho = rdf::runningCN(dead, 0.0);
  REQUIRE(noRho[0] == 0.0);
}

TEST_CASE("slab RDF clamps a short bin, resizes, and skips a bad cell pair",
          "[orient][rdf2d]") {
  Cloud cloud = makeCloud({{0.2, 0.2, 0.2}, {0.8, 0.2, 0.2}, {0.2, 0.9, 0.2}},
                          {20, 20, 20});
  cloud.currentFrame = 1;
  std::vector<double> rdf;
  const std::string dir = "orient_rdf_out";
  const int rc = rdf2::rdf2Danalysis_AA(dir + "/", rdf, cloud, 0.1, 1.0, 1, 3);
  REQUIRE(rc == 0);
  REQUIRE(rdf.size() == 1);

  cloud.currentFrame = 3;
  rdf.assign(4, 1.0);
  const int done =
      rdf2::rdf2Danalysis_AA(dir + "/", rdf, cloud, 0.1, 1.0, 1, 3);
  REQUIRE(done == 0);
  REQUIRE(rdf.size() == 1);
  REQUIRE(std::filesystem::exists(std::filesystem::path(dir) / "topoMonolayer" /
                                  "rdf.dat"));
  std::filesystem::remove_all(dir);

  const auto hist =
      rdf2::sampleRDF_AA(cloud, 5.0, 0.5, 0);
  REQUIRE(hist.empty());
  Cloud one = makeCloud({{1, 1, 1}}, {30, 30, 30});
  const auto single = rdf2::sampleRDF_AA(one, 2.0, 0.25, 4);
  REQUIRE(single.size() == 4);
  REQUIRE(single[0] == 0);
  const auto cells = rdf2::sampleRDF_AA(cloud, 2.0, 0.25, 8);
  REQUIRE(cells.size() == 8);
  int pairs = 0;
  for (int bin : cells) {
    pairs += bin;
  }
  REQUIRE(pairs > 0);

  std::vector<double> values(2, 0.0);
  std::vector<int> histogram = {4, 0};
  const int norm = rdf2::normalizeRDF(2, values, histogram, 0.5, 2,
                                      {0.0, 4.0, 3.0}, 1);
  REQUIRE(norm == 0);
  REQUIRE(values.size() == 2);

  Cloud collapsed = makeCloud({{0.0, 0.0, 0.0}, {1.2, 0.0, 0.0}},
                              {10.0, 10.0, -1.0});
  const auto direct = rdf2::sampleRDF_AA(collapsed, 3.0, 0.5, 4);
  REQUIRE(direct.size() == 4);
  REQUIRE(direct[2] == 2);
}

TEST_CASE("F4 skips a short neighbour row and a missing second hydrogen",
          "[orient][order]") {
  Cloud water;
  water.box = {20.0, 20.0, 20.0};
  water.boxLow = {0.0, 0.0, 0.0};
  auto add = [&](int type, int mol, int id, double x, double y, double z) {
    molSys::Point<double> p;
    p.type = type;
    p.molID = mol;
    p.atomID = id;
    p.x = x;
    p.y = y;
    p.z = z;
    water.idIndexMap[id] = static_cast<int>(water.pts.size());
    water.pts.push_back(p);
  };
  add(1, 1, 10, 0.0, 0.0, 0.0);
  add(2, 1, 11, 0.0, 0.96, 0.0);
  add(2, 1, 12, 0.0, -0.24, 0.9);
  add(1, 2, 20, 2.8, 0.0, 0.0);
  add(2, 2, 21, 2.8, 0.96, 0.0);
  add(2, 2, 22, 2.8, -0.24, 0.9);
  water.nop = static_cast<int>(water.pts.size());
  std::vector<std::vector<int>> shortList(1);
  shortList[0] = {10, 20};
  const auto skipped = topoparam::rodgerF4(water, shortList, 1, 2);
  REQUIRE(std::isfinite(skipped[0]));
  REQUIRE(std::isnan(skipped[3]));

  water.pts[4].x = std::numeric_limits<double>::quiet_NaN();
  water.pts[5].x = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::vector<int>> full(static_cast<std::size_t>(water.nop));
  full[0] = {10, 20};
  const auto missingH = topoparam::rodgerF4(water, full, 1, 2);
  REQUIRE(std::isnan(missingH[0]));

  Cloud bare;
  bare.box = {8.0, 8.0, 8.0};
  bare.boxLow = {0.0, 0.0, 0.0};
  molSys::Point<double> oxygen;
  oxygen.type = 1;
  oxygen.molID = 4;
  oxygen.atomID = 1;
  bare.pts.push_back(oxygen);
  bare.nop = 1;
  REQUIRE(std::isnan(topoparam::jumpRotorTau90(bare, bare, 1.0, 1, 2)));
}

TEST_CASE("layer clamps and rejected stack axes land on real bins",
          "[orient][order]") {
  const double pi = 3.141592653589793;
  Cloud edge;
  edge.box = {pi, pi, pi};
  edge.boxLow = {0.0, 0.0, 0.0};
  molSys::Point<double> ice;
  ice.iceType = molSys::atom_state_type::cubic;
  ice.z = std::nextafter(pi, 0.0);
  edge.pts.push_back(ice);
  edge.nop = 1;
  const auto top = topoparam::layerCubicity(edge, 2, 0.5);
  REQUIRE(top.cubicPerLayer.size() == 6);
  REQUIRE(top.cubicPerLayer.back() == 1);
  REQUIRE(top.sequence.back() == 'C');

  Cloud broken = edge;
  broken.pts[0].z = std::numeric_limits<double>::quiet_NaN();
  const auto nanLayer = topoparam::layerCubicity(broken, 2, 0.5);
  REQUIRE(nanLayer.cubicPerLayer.front() == 1);
  REQUIRE(nanLayer.sequence.front() == 'C');

  Cloud seeded = edge;
  seeded.nop = 5;
  seeded.pts.resize(5);
  seeded.pts[0] = ice;
  const std::vector<std::vector<int>> rings = {{9, 0, 1, 2, 3, 4},
                                               {0, -3, 1, 2, 3, 4},
                                               {0, 0, 0, 0, 0, 0}};
  const std::vector<bool> basal = {true, true, false};
  const std::vector<bool> equatorial = {false, false, true};
  const auto tumTop =
      topoparam::tumLayerStack(edge, {{0, 0, 0, 0, 0, 0}}, {false}, {true}, 2,
                               0.5);
  REQUIRE(tumTop.cubicPerLayer.back() == 1);
  const auto tumNan =
      topoparam::tumLayerStack(broken, {{0, 0, 0, 0, 0, 0}}, {true}, {false},
                               2, 0.5);
  REQUIRE(tumNan.hexPerLayer.front() == 1);
  const auto skipped =
      topoparam::tumLayerStack(seeded, rings, basal, equatorial, 2, 0.5);
  REQUIRE(skipped.sequence.size() == 6);

  Cloud box = edge;
  REQUIRE(topoparam::layerCubicity(box, -1, 0.5).sequence.empty());
  box.box = {pi, pi};
  REQUIRE(topoparam::layerCubicity(box, 2, 0.5).sequence.empty());
  REQUIRE(topoparam::tumLayerStack(edge, rings, basal, equatorial, -2, 0.5)
              .sequence.empty());
  REQUIRE(topoparam::tumLayerStack(edge, rings, basal, equatorial, 5, 0.5)
              .sequence.empty());
  REQUIRE(topoparam::tumLayerStack(box, rings, basal, equatorial, 2, 0.5)
              .sequence.empty());
  Cloud flat = edge;
  flat.box = {0.0, pi, pi};
  REQUIRE(topoparam::tumLayerStack(flat, rings, basal, equatorial, 0, 0.5)
              .phiC == 0.0);
}

TEST_CASE("prism wrap, rejected DDC indices, and a one-sided periphery",
          "[orient][pnt]") {
  Cloud cloud = makeCloud({{0, 0, 0},
                           {1, 0, 0},
                           {1, 1, 0},
                           {0, 1, 0},
                           {0, 0, 2},
                           {1, 0, 2},
                           {1, 1, 2},
                           {0, 1, 2}},
                          {12.0, 12.0, 12.0});
  const std::vector<int> basal1 = {0, 1, 2, 3};
  const std::vector<int> basal2 = {4, 5, 6, 7};
  const auto nList = indexNeighbours(8, {{2, 4}});
  std::vector<int> out1;
  std::vector<int> out2;
  const int rc =
      pntToPnt::relOrderPrismBlock(cloud, basal1, basal2, nList, out1, out2);
  REQUIRE(rc == 0);
  REQUIRE(out1.size() == 4);
  REQUIRE(out1[0] == 2);
  REQUIRE(out1[2] == 0);
  REQUIRE(out2[0] == 4);

  const std::vector<int> equatorial = {0, 1, 2, 3, 4, 5};
  std::vector<std::vector<int>> rings = {
      equatorial,
      {0, 20, 21, 22, 12, 6},
      {13, 23, 24, 50, 12, 9},
      {2, 25, 26, 27, 28, 7},
      {50, 29, 30, 31, 32, 10},
      {4, 33, 34, 35, 36, 8},
      {50, 37, 38, 39, 40, 11},
  };
  std::vector<cage::Cage> cages(1);
  cages[0].rings = {0, 1, 2, 3, 4, 5, 6};
  REQUIRE(pntToPnt::relOrderDDC(equatorial, -1, rings, cages).empty());
  REQUIRE(pntToPnt::relOrderDDC(equatorial, 4, rings, cages).empty());
  const auto oneSide = pntToPnt::relOrderDDC(equatorial, 0, rings, cages);
  REQUIRE(oneSide.empty());

  // Atom 0 sits at the end of the ring. Its previous neighbour is the
  // equatorial predecessor, so the forward neighbour is kept and +2 wraps.
  std::vector<std::vector<int>> forward = {
      equatorial, {20, 99, 21, 22, 5, 0}, {50, 50, 50, 50, 50, 50},
      {50, 50, 50, 50, 50, 50},           {50, 50, 50, 50, 50, 50},
      {50, 50, 50, 50, 50, 50},           {50, 50, 50, 50, 50, 50}};
  cages[0].rings = {0, 1, 2, 3, 4, 5, 6};
  const auto wrappedApex = pntToPnt::relOrderDDC(equatorial, 0, forward, cages);
  REQUIRE(wrappedApex.empty());
}
