#include <catch2/catch_test_macros.hpp>

#include <mol_sys.hpp>
#include <site.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// Eight "water" vertices on a cube of edge 4 whose centre sits at (2, 2, 2)
// and a second cube shifted along x by the box length minus its edge, so
// that its vertices straddle the periodic boundary; guests are appended.
Cloud twoCubes(const std::vector<std::array<double, 3>> &guests) {
  Cloud cloud;
  const double L = 20.0;
  int id = 0;
  auto add = [&](double x, double y, double z, int type) {
    molSys::Point<double> p;
    p.x = std::fmod(x + L, L);
    p.y = std::fmod(y + L, L);
    p.z = std::fmod(z + L, L);
    p.type = type;
    p.atomID = id + 1;
    p.molID = id + 1;
    cloud.pts.push_back(p);
    cloud.idIndexMap[id + 1] = id;
    ++id;
  };
  for (int c = 0; c < 2; c++) {
    const double ox = c == 0 ? 0.0 : L - 2.0;  // second cube spans x in [18, 22) -> wraps
    for (int i = 0; i < 8; i++) {
      add(ox + 4.0 * (i & 1), 4.0 * ((i >> 1) & 1), 4.0 * ((i >> 2) & 1), 1);
    }
  }
  for (const auto &g : guests) {
    add(g[0], g[1], g[2], 2);
  }
  cloud.nop = id;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

std::vector<int> range(int lo, int hi) {
  std::vector<int> v;
  for (int i = lo; i < hi; i++) {
    v.push_back(i);
  }
  return v;
}

} // namespace

TEST_CASE("the periodic centroid of a cage across the boundary is inside it", "[guests]") {
  const auto cloud = twoCubes({});
  const auto c0 = site::periodicCentroid(cloud, range(0, 8));
  REQUIRE(std::fabs(c0[0] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[1] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[2] - 2.0) < 1e-9);
  const auto c1 = site::periodicCentroid(cloud, range(8, 16));
  // the wrapped cube spans x in [18, 22): centre at 20 == 0 modulo the box
  const double wrapped = std::fmod(c1[0] + 20.0, 20.0);
  REQUIRE((std::fabs(wrapped) < 1e-9 || std::fabs(wrapped - 20.0) < 1e-9));
  REQUIRE(std::fabs(c1[1] - 2.0) < 1e-9);
}

TEST_CASE("guests are assigned to the cage whose centre they are nearest", "[guests]") {
  // one guest at the first centre, one near the wrapped centre seen from
  // the far side of the box, one in open space
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {19.5, 2.0, 2.2}, {10.0, 10.0, 10.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17, 18}, 3.0);
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, 1, -1});
  REQUIRE(occ.guestsPerCage == std::vector<int>{1, 1});
  REQUIRE(occ.occupied == 2);
  REQUIRE(occ.multiply == 0);
  REQUIRE(occ.free == 1);
  REQUIRE(std::fabs(occ.centreDistance[0]) < 1e-9);
  REQUIRE(std::fabs(occ.centreDistance[1] - std::sqrt(0.25 + 0.04)) < 1e-9);
  REQUIRE(occ.centreDistance[2] == -1.0);
}

std::vector<std::vector<int>> cubeFaces(int base) {
  return {
      {base + 0, base + 1, base + 3, base + 2},
      {base + 4, base + 5, base + 7, base + 6},
      {base + 0, base + 1, base + 5, base + 4},
      {base + 2, base + 3, base + 7, base + 6},
      {base + 0, base + 2, base + 6, base + 4},
      {base + 1, base + 3, base + 7, base + 5},
  };
}

TEST_CASE("ray parity puts a guest inside, on a face, or outside", "[guests]") {
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {2.0, 2.0, 4.0}, {2.0, 2.0, 4.5},
                               {19.5, 2.0, 2.2}, {10.0, 10.0, 10.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const std::vector<std::vector<std::vector<int>>> faces = {cubeFaces(0),
                                                            cubeFaces(8)};
  const auto inside =
      site::guestOccupancyInside(cloud, cages, faces, {16, 17, 18, 19, 20});
  REQUIRE(inside.cageOfGuest == std::vector<int>{0, 0, -1, 1, -1});
  REQUIRE(inside.guestsPerCage == std::vector<int>{2, 1});
  REQUIRE(inside.occupied == 2);
  REQUIRE(inside.multiply == 1);
  REQUIRE(inside.free == 2);
}

TEST_CASE("both predicates agree on centred guests and differ on a near miss",
          "[guests]") {
  // (2,2,4.5) is outside the cube of half-width 2 but within r=3 of the centre
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {2.0, 2.0, 4.5}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const std::vector<std::vector<std::vector<int>>> faces = {cubeFaces(0),
                                                            cubeFaces(8)};
  const auto both =
      site::guestOccupancyBoth(cloud, cages, faces, {16, 17}, 3.0);
  REQUIRE(both.radius.cageOfGuest == std::vector<int>{0, 0});
  REQUIRE(both.inside.cageOfGuest == std::vector<int>{0, -1});
  REQUIRE(both.radius.occupied == 1);
  REQUIRE(both.inside.occupied == 1);
  REQUIRE(both.radius.free == 0);
  REQUIRE(both.inside.free == 1);
}

TEST_CASE("two guests in one cage count as multiple occupancy", "[guests]") {
  const auto cloud = twoCubes({{1.5, 2.0, 2.0}, {2.5, 2.0, 2.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17}, 3.0);
  REQUIRE(occ.guestsPerCage == std::vector<int>{2, 0});
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.multiply == 1);
  REQUIRE(occ.free == 0);
  // a radius too small leaves every guest free
  const auto none = site::guestOccupancy(cloud, cages, {16, 17}, 0.1);
  REQUIRE(none.free == 2);
  REQUIRE(none.occupied == 0);
}
