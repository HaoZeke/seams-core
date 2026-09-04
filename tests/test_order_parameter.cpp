#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>

// Helper to build a PointCloud from a list of (x,y,z) coordinates
static molSys::PointCloud<molSys::Point<double>, double>
makeCloud(const std::vector<std::array<double, 3>> &coords,
          double boxLen = 100.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < static_cast<int>(coords.size()); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i;
    p.molID = 0;
    p.x = coords[i][0];
    p.y = coords[i][1];
    p.z = coords[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = static_cast<int>(coords.size());
  return cloud;
}

// Regression test for projAreaSingleRing return order.
// The function should return {areaXY, areaXZ, areaYZ} (index 0=XY, 1=XZ, 2=YZ).
// The caller calcCoverageArea reads [0]=XY, [1]=XZ, [2]=YZ.
TEST_CASE("projAreaSingleRing returns areas in XY, XZ, YZ order",
          "[order_parameter]") {
  // A 2x3 rectangle in the XZ plane (y=5 for all points).
  // Vertices: (0,5,0), (2,5,0), (2,5,3), (0,5,3)
  // XY projected area: all y are the same, so area = 0
  // XZ projected area: 2 * 3 = 6
  // YZ projected area: all y are the same, so area = 0
  auto cloud = makeCloud({{0, 5, 0}, {2, 5, 0}, {2, 5, 3}, {0, 5, 3}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE(areas.size() == 3);
  // [0] = XY area = 0 (flat in XY since all y identical)
  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  // [1] = XZ area = 6.0
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(6.0, 1e-10));
  // [2] = YZ area = 0 (flat in YZ since all y identical)
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

// Additional test: rectangle in XY plane
TEST_CASE("projAreaSingleRing XY plane rectangle", "[order_parameter]") {
  // A 4x3 rectangle in the XY plane (z=0 for all).
  // Vertices: (0,0,0), (4,0,0), (4,3,0), (0,3,0)
  // XY area: 4*3 = 12
  // XZ area: 0 (all z identical)
  // YZ area: 0 (all z identical)
  auto cloud = makeCloud({{0, 0, 0}, {4, 0, 0}, {4, 3, 0}, {0, 3, 0}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(12.0, 1e-10));
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("rodgerF4 is cos 3 phi on a known H-O-O-H dihedral",
          "[order_parameter]") {
  // Two waters. Outer hydrogens give a 90 degree H-O-O-H dihedral
  // (phi = pi/2, cos 3phi = 0) or a planar 0 degree one (cos 3phi = 1).
  auto twoWaters = [](double h2y, double h2z) {
    // O1 (0,0,0) mol 1; H_outer (0,1,0); H_inner (1,0,0)
    // O2 (3,0,0) mol 2; H_outer (3,h2y,h2z); H_inner (2,0,0)
    auto cloud = makeCloud({{0, 0, 0}, {0, 1, 0}, {1, 0, 0},
                            {3, 0, 0}, {3, h2y, h2z}, {2, 0, 0}});
    cloud.pts[0].type = 1;
    cloud.pts[0].molID = 1;
    cloud.pts[0].atomID = 1;
    cloud.pts[1].type = 2;
    cloud.pts[1].molID = 1;
    cloud.pts[1].atomID = 2;
    cloud.pts[2].type = 2;
    cloud.pts[2].molID = 1;
    cloud.pts[2].atomID = 3;
    cloud.pts[3].type = 1;
    cloud.pts[3].molID = 2;
    cloud.pts[3].atomID = 4;
    cloud.pts[4].type = 2;
    cloud.pts[4].molID = 2;
    cloud.pts[4].atomID = 5;
    cloud.pts[5].type = 2;
    cloud.pts[5].molID = 2;
    cloud.pts[5].atomID = 6;
    cloud.idIndexMap.clear();
    for (int i = 0; i < 6; i++) {
      cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    return cloud;
  };
  // nList by atom ID, leading self
  const std::vector<std::vector<int>> nList = {
      {1, 4}, {2}, {3}, {4, 1}, {5}, {6}};
  auto planar = twoWaters(1.0, 0.0);
  const auto f0 = topoparam::rodgerF4(planar, nList, 1, 2);
  REQUIRE(std::isfinite(f0[0]));
  REQUIRE_THAT(f0[0], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(f0[3], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(topoparam::meanFinite(f0), Catch::Matchers::WithinAbs(1.0, 1e-9));

  auto right = twoWaters(0.0, 1.0);
  const auto f90 = topoparam::rodgerF4(right, nList, 1, 2);
  REQUIRE_THAT(f90[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(f90[3], Catch::Matchers::WithinAbs(0.0, 1e-9));
}

TEST_CASE("rodgerF4 is NaN on mW with no hydrogens", "[order_parameter]") {
  auto cloud = makeCloud({{0, 0, 0}, {3, 0, 0}});
  cloud.pts[0].atomID = 1;
  cloud.pts[1].atomID = 2;
  cloud.idIndexMap = {{1, 0}, {2, 1}};
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  const auto f4 = topoparam::rodgerF4(cloud, nList, 1, 2);
  REQUIRE_FALSE(std::isfinite(f4[0]));
  REQUIRE_FALSE(std::isfinite(f4[1]));
  REQUIRE_FALSE(std::isfinite(topoparam::meanFinite(f4)));
}

TEST_CASE("CHILL+ layerCubicity reproduces the literature I_sd string",
          "[order_parameter]") {
  // Four layers along z at 0, 3.7, 7.4, 11.1 in a 14.8 box: C H C H
  std::vector<std::array<double, 3>> coords;
  for (int k = 0; k < 4; k++) {
    coords.push_back({1.0, 1.0, k * 3.7});
    coords.push_back({2.0, 2.0, k * 3.7});
  }
  auto cloud = makeCloud(coords, 14.8);
  for (int i = 0; i < 8; i++) {
    cloud.pts[static_cast<std::size_t>(i)].iceType =
        (i / 2) % 2 == 0 ? molSys::atom_state_type::cubic
                         : molSys::atom_state_type::hexagonal;
  }
  const auto st = topoparam::layerCubicity(cloud, 2, 3.7);
  REQUIRE(st.sequence == "CHCH");
  REQUIRE_THAT(st.phiC, Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE(st.phiC > 0.0);
  REQUIRE(st.phiC < 1.0);
}

TEST_CASE("TUM stacking uses ring planes and stays empty without basal rings",
          "[order_parameter]") {
  // Two disjoint six-rings stacked along z, marked basal. Their centroids
  // fall in two H layers. A five-ring (clathrate face) does not vote.
  auto cloud = makeCloud({{0, 0, 0}, {2, 0, 0}, {3, 1.5, 0},
                          {2, 3, 0}, {0, 3, 0}, {-1, 1.5, 0},
                          {0, 0, 7.4}, {2, 0, 7.4}, {3, 1.5, 7.4},
                          {2, 3, 7.4}, {0, 3, 7.4}, {-1, 1.5, 7.4},
                          {1, 1, 3.7}, {2, 1, 3.7}, {2.5, 2, 3.7},
                          {1.5, 2.8, 3.7}, {0.5, 2, 3.7}},
                         14.8);
  const std::vector<std::vector<int>> rings = {
      {0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}, {12, 13, 14, 15, 16}};
  std::vector<bool> basal = {true, true, false};
  std::vector<bool> equatorial = {false, false, false};
  const auto tum = topoparam::tumLayerStack(cloud, rings, basal, equatorial, 2, 3.7);
  REQUIRE(tum.sequence.size() == 4);
  REQUIRE(tum.sequence[0] == 'H');
  REQUIRE(tum.sequence[2] == 'H');
  REQUIRE(tum.phiC == 0.0);
  // the five-ring is not a plane: no C letter
  REQUIRE(tum.sequence.find('C') == std::string::npos);

  std::vector<bool> eq = {false, true, false};
  const auto mixed = topoparam::tumLayerStack(cloud, rings, basal, eq, 2, 3.7);
  REQUIRE(mixed.phiC > 0.0);
  REQUIRE(mixed.phiC < 1.0);
  REQUIRE(mixed.sequence.find('C') != std::string::npos);
  REQUIRE(mixed.sequence.find('H') != std::string::npos);
}

static void wrapCoord(double &x, double L) {
  if (!(L > 0.0)) {
    return;
  }
  x -= L * std::floor(x / L);
  if (x < 0.0) {
    x += L;
  }
  if (x >= L) {
    x = 0.0;
  }
}

static bool nearExisting(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double x,
    double y, double z) {
  for (const auto &p : cloud.pts) {
    double dx = p.x - x;
    double dy = p.y - y;
    double dz = p.z - z;
    if (cloud.box.size() >= 3) {
      const double bx = cloud.box[0];
      const double by = cloud.box[1];
      const double bz = cloud.box[2];
      dx -= bx * std::round(dx / bx);
      dy -= by * std::round(dy / by);
      dz -= bz * std::round(dz / bz);
    }
    if (dx * dx + dy * dy + dz * dz < 0.04) {
      return true;
    }
  }
  return false;
}

// Ice Ih oxygens: wurtzite 4f sites tiled in an orthorhombic cell
// (two hexagonal cells: a x a sqrt(3) x c).
static molSys::PointCloud<molSys::Point<double>, double> iceIhOxygens(int nx,
                                                                     int ny,
                                                                     int nz) {
  const double a = 4.5112;
  const double c = 7.351;
  const double u = 0.0622;
  const double fx[4] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0};
  const double fy[4] = {2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0};
  const double fz[4] = {u, 0.5 + u, 0.5 - u, 1.0 - u};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {a * static_cast<double>(nx),
               a * std::sqrt(3.0) * static_cast<double>(ny),
               c * static_cast<double>(nz)};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  int atomID = 1;
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        for (int ja = 0; ja < 2; ja++) {
          for (int k = 0; k < 4; k++) {
            const double hx = static_cast<double>(ix) + fx[k];
            const double hy = static_cast<double>(2 * iy + ja) + fy[k];
            const double hz = static_cast<double>(iz) + fz[k];
            double x = hx * a + hy * (0.5 * a);
            double y = hy * (a * std::sqrt(3.0) * 0.5);
            double z = hz * c;
            wrapCoord(x, cloud.box[0]);
            wrapCoord(y, cloud.box[1]);
            wrapCoord(z, cloud.box[2]);
            if (nearExisting(cloud, x, y, z)) {
              continue;
            }
            molSys::Point<double> p;
            p.type = 1;
            p.atomID = atomID;
            p.molID = atomID;
            p.x = x;
            p.y = y;
            p.z = z;
            cloud.pts.push_back(p);
            cloud.idIndexMap[atomID] = static_cast<int>(cloud.pts.size()) - 1;
            ++atomID;
          }
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

// Ice Ic oxygens on the diamond lattice (same tetrahedral ice-I F4 as Ih).
static molSys::PointCloud<molSys::Point<double>, double> iceIcOxygens(int n) {
  const double a = 6.358;
  const double frac[8][3] = {
      {0.00, 0.00, 0.00}, {0.50, 0.50, 0.00}, {0.50, 0.00, 0.50},
      {0.00, 0.50, 0.50}, {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
      {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {a * static_cast<double>(n), a * static_cast<double>(n),
               a * static_cast<double>(n)};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  int atomID = 1;
  for (int iz = 0; iz < n; iz++) {
    for (int iy = 0; iy < n; iy++) {
      for (int ix = 0; ix < n; ix++) {
        for (int k = 0; k < 8; k++) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = atomID;
          p.molID = atomID;
          p.x = (static_cast<double>(ix) + frac[k][0]) * a;
          p.y = (static_cast<double>(iy) + frac[k][1]) * a;
          p.z = (static_cast<double>(iz) + frac[k][2]) * a;
          wrapCoord(p.x, cloud.box[0]);
          wrapCoord(p.y, cloud.box[1]);
          wrapCoord(p.z, cloud.box[2]);
          cloud.pts.push_back(p);
          cloud.idIndexMap[atomID] = static_cast<int>(cloud.pts.size()) - 1;
          ++atomID;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

// Place two hydrogens per oxygen along an Eulerian orientation of the
// mutual four-neighbour graph (ice rules).
static void addIceRuleHydrogens(
    molSys::PointCloud<molSys::Point<double>, double> &cloud, int oType,
    int hType) {
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.0, oType, true);
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < cloud.nop; i++) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != oType) {
      continue;
    }
    if (static_cast<std::size_t>(i) >= nList.size()) {
      continue;
    }
    for (std::size_t m = 1; m < nList[static_cast<std::size_t>(i)].size();
         m++) {
      const int jid = nList[static_cast<std::size_t>(i)][m];
      const auto it = cloud.idIndexMap.find(jid);
      if (it == cloud.idIndexMap.end()) {
        continue;
      }
      const int j = it->second;
      if (j <= i) {
        continue;
      }
      if (cloud.pts[static_cast<std::size_t>(j)].type != oType) {
        continue;
      }
      edges.emplace_back(i, j);
    }
  }
  const int nE = static_cast<int>(edges.size());
  std::vector<std::vector<std::pair<int, int>>> adj(
      static_cast<std::size_t>(cloud.nop));
  for (int e = 0; e < nE; e++) {
    adj[static_cast<std::size_t>(edges[static_cast<std::size_t>(e)].first)]
        .push_back({edges[static_cast<std::size_t>(e)].second, e});
    adj[static_cast<std::size_t>(edges[static_cast<std::size_t>(e)].second)]
        .push_back({edges[static_cast<std::size_t>(e)].first, e});
  }
  std::vector<char> used(static_cast<std::size_t>(nE), 0);
  std::vector<int> donor(static_cast<std::size_t>(nE), -1);
  for (int start = 0; start < cloud.nop; start++) {
    std::vector<int> st = {start};
    while (!st.empty()) {
      const int v = st.back();
      bool found = false;
      for (const auto &pr : adj[static_cast<std::size_t>(v)]) {
        if (used[static_cast<std::size_t>(pr.second)] != 0) {
          continue;
        }
        used[static_cast<std::size_t>(pr.second)] = 1;
        donor[static_cast<std::size_t>(pr.second)] = v;
        st.push_back(pr.first);
        found = true;
        break;
      }
      if (!found) {
        st.pop_back();
      }
    }
  }
  int nextID = 0;
  for (const auto &p : cloud.pts) {
    nextID = std::max(nextID, p.atomID);
  }
  ++nextID;
  const int nO = cloud.nop;
  for (int e = 0; e < nE; e++) {
    const int d = donor[static_cast<std::size_t>(e)];
    if (d < 0) {
      continue;
    }
    const int a = edges[static_cast<std::size_t>(e)].first == d
                      ? edges[static_cast<std::size_t>(e)].second
                      : edges[static_cast<std::size_t>(e)].first;
    const auto toA = gen::relDist(cloud, a, d);
    const double len =
        std::sqrt(toA[0] * toA[0] + toA[1] * toA[1] + toA[2] * toA[2]);
    if (!(len > 0.5)) {
      continue;
    }
    const double f = 0.957 / len;
    molSys::Point<double> h;
    h.type = hType;
    h.atomID = nextID++;
    h.molID = cloud.pts[static_cast<std::size_t>(d)].molID;
    h.x = cloud.pts[static_cast<std::size_t>(d)].x + f * toA[0];
    h.y = cloud.pts[static_cast<std::size_t>(d)].y + f * toA[1];
    h.z = cloud.pts[static_cast<std::size_t>(d)].z + f * toA[2];
    wrapCoord(h.x, cloud.box[0]);
    wrapCoord(h.y, cloud.box[1]);
    wrapCoord(h.z, cloud.box[2]);
    cloud.pts.push_back(h);
    cloud.idIndexMap[h.atomID] = static_cast<int>(cloud.pts.size()) - 1;
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  REQUIRE(cloud.nop == nO + nE);
}

static bool tetrahedralFourNN(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, int oType) {
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.0, oType, true);
  int nO = 0;
  for (int i = 0; i < cloud.nop; i++) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != oType) {
      continue;
    }
    ++nO;
    if (static_cast<std::size_t>(i) >= nList.size() ||
        nList[static_cast<std::size_t>(i)].size() != 5) {
      return false;
    }
    for (std::size_t m = 1; m < nList[static_cast<std::size_t>(i)].size();
         m++) {
      const int j = cloud.idIndexMap.at(nList[static_cast<std::size_t>(i)][m]);
      const auto dr = gen::relDist(cloud, i, j);
      const double d =
          std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
      if (d < 2.55 || d > 2.95) {
        return false;
      }
    }
  }
  return nO >= 8;
}

static double meanF4OnOxygens(
    molSys::PointCloud<molSys::Point<double>, double> &cloud, int oType,
    int hType) {
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.0, oType, true);
  const auto f4 = topoparam::rodgerF4(cloud, nList, oType, hType);
  return topoparam::meanFinite(f4);
}

// Small random water rotations stand in for thermal proton disorder.
static void rotateWaters(
    molSys::PointCloud<molSys::Point<double>, double> &cloud, int oType,
    int hType, double maxRad, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::uniform_real_distribution<double> uang(-maxRad, maxRad);
  for (int i = 0; i < cloud.nop; i++) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != oType) {
      continue;
    }
    std::vector<int> hs;
    const int mol = cloud.pts[static_cast<std::size_t>(i)].molID;
    for (int j = 0; j < cloud.nop; j++) {
      if (cloud.pts[static_cast<std::size_t>(j)].type == hType &&
          cloud.pts[static_cast<std::size_t>(j)].molID == mol) {
        hs.push_back(j);
      }
    }
    if (hs.size() < 2) {
      continue;
    }
    const double th = uang(rng);
    const double z = 2.0 * u01(rng) - 1.0;
    const double phi = 6.283185307179586 * u01(rng);
    const double rxy = std::sqrt(std::max(0.0, 1.0 - z * z));
    const double ax = rxy * std::cos(phi);
    const double ay = rxy * std::sin(phi);
    const double az = z;
    const double c = std::cos(th);
    const double s = std::sin(th);
    const double ox = cloud.pts[static_cast<std::size_t>(i)].x;
    const double oy = cloud.pts[static_cast<std::size_t>(i)].y;
    const double oz = cloud.pts[static_cast<std::size_t>(i)].z;
    for (int h : hs) {
      auto &p = cloud.pts[static_cast<std::size_t>(h)];
      double vx = p.x - ox;
      double vy = p.y - oy;
      double vz = p.z - oz;
      if (cloud.box.size() >= 3) {
        vx -= cloud.box[0] * std::round(vx / cloud.box[0]);
        vy -= cloud.box[1] * std::round(vy / cloud.box[1]);
        vz -= cloud.box[2] * std::round(vz / cloud.box[2]);
      }
      const double dot = ax * vx + ay * vy + az * vz;
      const double rx = vx * c + (ay * vz - az * vy) * s + ax * dot * (1.0 - c);
      const double ry = vy * c + (az * vx - ax * vz) * s + ay * dot * (1.0 - c);
      const double rz = vz * c + (ax * vy - ay * vx) * s + az * dot * (1.0 - c);
      p.x = ox + rx;
      p.y = oy + ry;
      p.z = oz + rz;
      wrapCoord(p.x, cloud.box[0]);
      wrapCoord(p.y, cloud.box[1]);
      wrapCoord(p.z, cloud.box[2]);
    }
  }
}

TEST_CASE("rodgerF4 is near -0.4 on tetrahedral ice I when hydrogens exist",
          "[order_parameter]") {
  auto ice = iceIhOxygens(3, 2, 2);
  if (!tetrahedralFourNN(ice, 1)) {
    ice = iceIcOxygens(2);
  }
  REQUIRE(tetrahedralFourNN(ice, 1));
  addIceRuleHydrogens(ice, 1, 2);
  rotateWaters(ice, 1, 2, 0.82, 11u);
  const double fIce = meanF4OnOxygens(ice, 1, 2);
  INFO("ice I F4=" << fIce << " nop=" << ice.nop);
  REQUIRE(std::isfinite(fIce));
  REQUIRE(fIce >= -0.55);
  REQUIRE(fIce <= -0.25);
}

TEST_CASE("rodgerF4 is near 0.7 on GenIce sI when hydrogens exist",
          "[order_parameter]") {
  molSys::PointCloud<molSys::Point<double>, double> sI;
  sI = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, sI, 1);
  REQUIRE(sI.nop == 46);
  REQUIRE(tetrahedralFourNN(sI, 1));
  for (int i = 0; i < sI.nop; i++) {
    sI.pts[static_cast<std::size_t>(i)].molID = i + 1;
    sI.pts[static_cast<std::size_t>(i)].type = 1;
  }
  addIceRuleHydrogens(sI, 1, 2);
  rotateWaters(sI, 1, 2, 0.35, 3u);
  const double fSI = meanF4OnOxygens(sI, 1, 2);
  INFO("sI F4=" << fSI);
  REQUIRE(std::isfinite(fSI));
  REQUIRE(fSI >= 0.55);
  REQUIRE(fSI <= 0.85);
}

TEST_CASE("normHeightPercent uses recovered lz not tilt",
          "[order_parameter]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0}});
  cloud.box = {61.0, 12.0, 50.0, 60.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double h = topoparam::normHeightPercent(cloud, 10, 2.5);
  REQUIRE_THAT(h, Catch::Matchers::WithinAbs(50.0, 1e-10));
}
