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

#include <ira_sofi.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <limits>

TEST_CASE("ira::available matches the compile-time flag", "[ira]") {
#ifdef SEAMS_HAS_IRA
  REQUIRE(ira::available());
#else
  REQUIRE_FALSE(ira::available());
  ira::Match m;
  Eigen::MatrixXd pts(3, 3);
  pts.setZero();
  REQUIRE(ira::match(pts, pts, m) == 1);
  ira::PointGroup pg;
  REQUIRE(ira::pointGroup(pts, pg) == 1);
#endif
}

#ifdef SEAMS_HAS_IRA

static Eigen::MatrixXd squareXY() {
  Eigen::MatrixXd p(4, 3);
  p << 1.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, -1.0, 0.0, 1.0, -1.0, 0.0;
  return p;
}

TEST_CASE("IRA residual is ref versus R times assigned target plus t",
          "[ira]") {
  // libira maps structure 2 onto structure 1. 90 deg is not an involution,
  // so R*ref+t and ref-(R*target+t) cannot both be ~0.
  const Eigen::MatrixXd ref = squareXY();
  Eigen::MatrixXd tgt(4, 3);
  tgt.row(0) = Eigen::RowVector3d(-1.0, 1.0, 0.0);
  tgt.row(1) = Eigen::RowVector3d(-1.0, -1.0, 0.0);
  tgt.row(2) = Eigen::RowVector3d(1.0, -1.0, 0.0);
  tgt.row(3) = Eigen::RowVector3d(1.0, 1.0, 0.0);
  ira::Match m;
  REQUIRE(ira::match(ref, tgt, m) == 0);
  REQUIRE_THAT(m.rmsd, Catch::Matchers::WithinAbs(0.0, 1e-6));
  double alt = 0.0;
  int used = 0;
  const int n = 4;
  bool hasZero = false;
  bool hasN = false;
  for (int p : m.assignment) {
    if (p == 0) {
      hasZero = true;
    }
    if (p == n) {
      hasN = true;
    }
  }
  const int shift = (hasN && !hasZero) ? 1 : 0;
  for (int i = 0; i < n; i++) {
    int j = (i < static_cast<int>(m.assignment.size()))
                ? static_cast<int>(m.assignment[static_cast<size_t>(i)]) - shift
                : -1;
    if (j < 0 || j >= n) {
      continue;
    }
    const Eigen::Vector3d a(ref(i, 0), ref(i, 1), ref(i, 2));
    const Eigen::Vector3d b(tgt(j, 0), tgt(j, 1), tgt(j, 2));
    alt += (m.rotation * a + m.translation - b).squaredNorm();
    used++;
  }
  if (used > 0) {
    alt = std::sqrt(alt / static_cast<double>(used));
  }
  REQUIRE(alt > 0.5);
}

TEST_CASE("IRA overlays a rotated and permuted square", "[ira]") {
  const Eigen::MatrixXd ref = squareXY();
  Eigen::MatrixXd tgt(4, 3);
  // 90 deg about z, then swap first two atoms.
  tgt.row(0) = Eigen::RowVector3d(-1.0, -1.0, 0.0);
  tgt.row(1) = Eigen::RowVector3d(-1.0, 1.0, 0.0);
  tgt.row(2) = Eigen::RowVector3d(1.0, 1.0, 0.0);
  tgt.row(3) = Eigen::RowVector3d(1.0, -1.0, 0.0);
  ira::Match m;
  REQUIRE(ira::match(ref, tgt, m) == 0);
  REQUIRE(m.assignment.size() == 4);
  REQUIRE_THAT(m.rmsd, Catch::Matchers::WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(m.hausdorff, Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("IRA assignments are 0-based indices inside the set", "[ira]") {
  ira::Match m;
  const Eigen::MatrixXd ref = squareXY();
  REQUIRE(ira::match(ref, ref, m) == 0);
  REQUIRE(m.assignment.size() == 4);
  bool sawZero = false;
  for (int p : m.assignment) {
    REQUIRE(p >= 0);
    REQUIRE(p < 4);
    sawZero = sawZero || (p == 0);
  }
  REQUIRE(sawZero);
}

TEST_CASE("IRA returns the library error when no basis exists", "[ira]") {
  Eigen::MatrixXd ref(4, 3);
  ref << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  Eigen::MatrixXd target(4, 3);
  target << 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.2, 0.0, 8.0, 8.0, 8.0;
  ira::Match m;
  REQUIRE(ira::match(ref, target, m, 1.01) != 0);
}

TEST_CASE("SOFI returns the library error when the basis cannot be set", "[ira]") {
  ira::PointGroup pg;
  const int rc = ira::pointGroup(squareXY(), pg, 10.0);
  REQUIRE(rc == -3);
  REQUIRE(pg.nOperations == 0);
}

TEST_CASE("SOFI reports a non-C1 group for a square", "[ira]") {
  ira::PointGroup pg;
  REQUIRE(ira::pointGroup(squareXY(), pg) == 0);
  REQUIRE(pg.nOperations >= 2);
  REQUIRE(pg.symbol != "C1");
  REQUIRE(pg.symbol != "C1 ");
}

TEST_CASE("IRA rejects a matrix that is not a point set", "[ira]") {
  ira::Match m;
  Eigen::MatrixXd wide(2, 2);
  wide.setZero();
  REQUIRE(ira::match(wide, wide, m) == 1);
  Eigen::MatrixXd empty(0, 3);
  REQUIRE(ira::match(empty, squareXY(), m) == 1);
  ira::PointGroup pg;
  REQUIRE(ira::pointGroup(wide, pg) == 1);
  REQUIRE(ira::pointGroup(empty, pg) == 1);
}

TEST_CASE("IRA matches a triangle onto a larger set", "[ira]") {
  Eigen::MatrixXd ref(3, 3);
  ref << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
  Eigen::MatrixXd target(4, 3);
  target << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 5.0, 5.0, 0.0;
  ira::Match m;
  const int rc = ira::match(ref, target, m);
  REQUIRE(rc == 0);
  REQUIRE(m.assignment.size() == 4);
  REQUIRE(m.rmsd >= 0.0);
}

TEST_CASE("IRA rejects a matrix that is not n by 3", "[ira]") {
  Eigen::MatrixXd wide(2, 2);
  wide.setZero();
  Eigen::MatrixXd empty(0, 3);
  ira::Match m;
  ira::PointGroup pg;
  REQUIRE(ira::match(wide, squareXY(), m) == 1);
  REQUIRE(ira::match(squareXY(), wide, m) == 1);
  REQUIRE(ira::match(empty, squareXY(), m) == 1);
  REQUIRE(ira::pointGroup(empty, pg) == 1);
}

TEST_CASE("IRA rejects a non-finite point set", "[ira]") {
  ira::Match m;
  for (const double badValue : {std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()}) {
    Eigen::MatrixXd bad = squareXY();
    bad(0, 0) = badValue;
    REQUIRE(ira::match(bad, squareXY(), m) == 1);
    REQUIRE(ira::match(squareXY(), bad, m) == 1);
  }
}

TEST_CASE("SOFI rejects a non-finite cloud before libira", "[ira]") {
  // A NaN coordinate used to enter libira_compute_all and segfault.
  Eigen::MatrixXd bad = squareXY();
  ira::PointGroup pg;
  for (const double badValue : {std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()}) {
    bad(1, 2) = badValue;
    REQUIRE(ira::pointGroup(bad, pg) == 1);
    REQUIRE(pg.nOperations == 0);
  }
}

TEST_CASE("orient copies the quaternion of a matched square", "[ira]") {
  const Eigen::MatrixXd ref = squareXY();
  std::vector<double> quat;
  double rmsd = -1.0;
  REQUIRE(ira::orient(ref, ref, quat, rmsd));
  REQUIRE(quat.size() == 4);
  REQUIRE_THAT(rmsd, Catch::Matchers::WithinAbs(0.0, 1e-6));
  Eigen::MatrixXd bad(1, 2);
  REQUIRE_FALSE(ira::orient(bad, bad, quat, rmsd));
}

#endif
