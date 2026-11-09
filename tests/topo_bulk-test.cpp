// Internal
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_bulk.hpp>

// Standard
#include <algorithm>
#include <cmath>
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

#include <catch2/catch_test_macros.hpp>
#include <rang.hpp>

namespace fs = std::filesystem;

// Helper: build the 12-atom HC geometry (reusable across tests)
static void buildHCCloud(molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud.box = {53.9690018, 54.5289994, 51.257};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  for (int i = 0; i < 12; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0]; iPoint.y = coords[i][1]; iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 12;
}

SCENARIO("Test the HC algorithm for a single hexagonal cage.", "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud;  // pointCloud
    molSys::Point<double> iPoint;                              // A single point
    std::vector<double> box = {53.9690018, 54.5289994, 51.257};  // Box lengths
    yCloud.box = box;  // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listHC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    // --------------------
    // Building the hexagonal cage
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 14.2810001;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 15.0939999;  // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 9.0880003;   // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 14.1949997;  // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;      // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 9.0880003;   // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;      // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 12.9829998;  // y
    iPoint.z = 11.4329996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;      // iatom
    iPoint.x = 8.995;       // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;     // iatom
    iPoint.x = 6.7459998;   // x
    iPoint.y = 14.2810001;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;     // iatom
    iPoint.x = 4.4970002;   // x
    iPoint.y = 10.3859997;  // y
    iPoint.z = 10.5340004;  // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = yCloud.pts.size();
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    // --------------------
    WHEN("Given a pointCloud, and a neighbour list") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Neighbour list by index
      nList = nneigh::neighbourListByIndex(yCloud, nList);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 7);
      THEN("There should be exactly one hexagonal cage.") {
        ringType.resize(
            rings.size());  // Has a value for each ring. init to zero.
        // Find the number of hexagonal cages
        listHC = ring::findHC(rings, ringType, nList, cageList);

        // Assert the number of cages
        REQUIRE(cageList.size() ==
                1);  // Evaluate condition for a single tetragonal prism
      }
    }  // End of getting the neighbour list
  }    // End of given
}  // End of scenario

// DDC determination
SCENARIO("Test the DDC algorithm for a single double-diamond cage.", "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud;  // pointCloud
    molSys::Point<double> iPoint;                              // A single point
    std::vector<double> box = {53.9690018, 54.5289994, 51.257};  // Box lengths
    yCloud.box = box;  // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listDDC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    std::vector<int> listHC;  // Empty
    // --------------------
    // Building the double-diamond cage
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 19.1100006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 25.4799996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 22.2950001;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;      // iatom
    iPoint.x = 25.4799996;  // x
    iPoint.y = 28.6650009;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;      // iatom
    iPoint.x = 22.2950001;  // x
    iPoint.y = 22.2950001;  // y
    iPoint.z = 25.4799996;  // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;      // iatom
    iPoint.x = 22.2950001;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;      // iatom
    iPoint.x = 28.6650009;  // x
    iPoint.y = 25.4799996;  // y
    iPoint.z = 22.2950001;  // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;      // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 20.7029991;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;      // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 27.073;      // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;      // iatom
    iPoint.x = 27.073;      // x
    iPoint.y = 27.073;      // y
    iPoint.z = 20.7029991;  // z
    yCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;     // iatom
    iPoint.x = 20.7029991;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;     // iatom
    iPoint.x = 27.073;      // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 23.8880006;  // z
    yCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 12;     // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 20.7029991;  // z
    yCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 13;     // iatom
    iPoint.x = 23.8880006;  // x
    iPoint.y = 23.8880006;  // y
    iPoint.z = 27.073;      // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = yCloud.pts.size();
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    // --------------------
    WHEN("Given a pointCloud") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Neighbour list by index
      nList = nneigh::neighbourListByIndex(yCloud, nList);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 7);
      THEN("There should be exactly one double-diamond cage.") {
        ringType.resize(
            rings.size());  // Has a value for each ring. init to zero.
        // Find the number of hexagonal cages
        listDDC = ring::findDDC(rings, ringType, listHC, cageList);
        // Assert the number of cages
        REQUIRE(cageList.size() ==
                1);  // Evaluate condition for a single tetragonal prism
      }
    }  // End of getting the neighbour list
  }    // End of given
} // End of scenario

// -- Additional tests for wrapper functions --

TEST_CASE("topoBulkAnalysis integration with HC geometry", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulk_analysis/").string();

  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}


TEST_CASE("getStrucNumbers counts HC/DDC structures correctly", "[topo_bulk]") {
  // Set up known ring types: 5 rings, 2 basal HC, 3 prismatic
  std::vector<ring::strucType> ringType = {
      ring::strucType::HCbasal, ring::strucType::HCbasal,
      ring::strucType::HCprismatic, ring::strucType::HCprismatic,
      ring::strucType::HCprismatic};

  // One HC cage referencing rings 0-4
  cage::Cage hcCage;
  hcCage.type = cage::cageType::HexC;
  hcCage.rings = {0, 1, 2, 3, 4};
  std::vector<cage::Cage> cageList = {hcCage};

  int numHC = 0, numDDC = 0, mixedRings = 0, prismaticRings = 0,
      basalRings = 0;
  int ret = ring::getStrucNumbers(ringType, cageList, numHC, numDDC,
                                   mixedRings, prismaticRings, basalRings);

  REQUIRE(ret == 0);
  REQUIRE(numHC == 1);
  REQUIRE(numDDC == 0);
  REQUIRE(basalRings == 2);
  REQUIRE(prismaticRings == 3);
}

TEST_CASE("getAtomTypesTopoBulk assigns atom types from ring types",
          "[topo_bulk]") {
  // 6-membered rings with known atoms
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5},
                                          {6, 7, 8, 9, 10, 11}};
  std::vector<ring::strucType> ringType = {ring::strucType::HCbasal,
                                            ring::strucType::DDC};
  std::vector<cage::iceType> atomTypes(12, cage::iceType::dummy);

  int ret = ring::getAtomTypesTopoBulk(rings, ringType, atomTypes);
  REQUIRE(ret == 0);

  // Atoms in ring 0 (HC basal) should be hc type
  REQUIRE(atomTypes[0] == cage::iceType::hc);
  // Atoms in ring 1 (DDC) should be ddc type
  REQUIRE(atomTypes[6] == cage::iceType::ddc);
}

TEST_CASE("findMixedRings identifies shared DDC/HC rings", "[topo_bulk]") {
  // Build the HC geometry, find rings, then test
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);
  auto listDDC = ring::findDDC(rings, ringType, listHC, cageList);

  // Find mixed rings (with a single HC, there should be none)
  auto mixedList = ring::findMixedRings(rings, ringType, listDDC, listHC);
  REQUIRE(mixedList.empty());
}

TEST_CASE("buildRingSearchIndex lists every ring that contains an atom",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {
      {0, 1, 2}, {2, 3, 4}, {0, 4, 5}, {9}};
  const auto index = ring::buildRingSearchIndex(rings, 6);

  REQUIRE(index.ringsContainingAtom.size() == 6);
  REQUIRE(index.ringsContainingAtom[0] == std::vector<int>{0, 2});
  REQUIRE(index.ringsContainingAtom[1] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[2] == std::vector<int>{0, 1});
  REQUIRE(index.ringsContainingAtom[3] == std::vector<int>{1});
  REQUIRE(index.ringsContainingAtom[4] == std::vector<int>{1, 2});
  REQUIRE(index.ringsContainingAtom[5] == std::vector<int>{2});
}

TEST_CASE("buildRingSearchIndex ignores atoms outside numAtoms",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {{0, 1, 99}, {-1, 2}};
  const auto index = ring::buildRingSearchIndex(rings, 3);
  REQUIRE(index.ringsContainingAtom.size() == 3);
  REQUIRE(index.ringsContainingAtom[0] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[1] == std::vector<int>{0});
  REQUIRE(index.ringsContainingAtom[2] == std::vector<int>{1});
}

TEST_CASE("conditionOneDDC rejects an atom outside the ring index",
          "[topo_bulk]") {
  const std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5}};
  ring::RingSearchIndex index;
  index.ringsContainingAtom.resize(3);
  std::vector<int> peripheral;
  REQUIRE_FALSE(ring::conditionOneDDC(rings, peripheral, 0, index));
}

TEST_CASE("findHC on the twelve-atom HC reports one HexC of five rings",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);

  REQUIRE(cageList.size() == 1);
  REQUIRE(cageList[0].type == cage::cageType::HexC);
  REQUIRE(cageList[0].rings.size() == 5);
  REQUIRE_FALSE(listHC.empty());
}

TEST_CASE("indexed findHC and findDDC match the convenience overloads",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<std::vector<int>> sixRings;
  int maxAtom = 0;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      sixRings.push_back(r);
      for (int atom : r) {
        maxAtom = std::max(maxAtom, atom);
      }
    }
  }

  const auto index =
      ring::buildRingSearchIndex(sixRings, maxAtom + 1);

  std::vector<ring::strucType> typeWrap(sixRings.size());
  std::vector<ring::strucType> typeIdx(sixRings.size());
  std::vector<cage::Cage> cagesWrap;
  std::vector<cage::Cage> cagesIdx;

  auto hcWrap = ring::findHC(sixRings, typeWrap, nList, cagesWrap);
  auto hcIdx = ring::findHC(sixRings, typeIdx, nList, cagesIdx, index);
  REQUIRE(hcWrap == hcIdx);
  REQUIRE(typeWrap == typeIdx);

  auto ddcWrap = ring::findDDC(sixRings, typeWrap, hcWrap, cagesWrap);
  auto ddcIdx = ring::findDDC(sixRings, typeIdx, hcIdx, cagesIdx, index);
  REQUIRE(ddcWrap == ddcIdx);
  REQUIRE(typeWrap == typeIdx);
}

TEST_CASE("topoBulkAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topobulk_mw/").string();
  int ret = ring::topoBulkAnalysis(tmpPath, rings, nList, yCloud, 1, true);
  REQUIRE(ret == 0);
  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("bulkPolygonRingAnalysis on mW cubic trajectory", "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_bulkpoly_mw/").string();
  int ret = ring::bulkPolygonRingAnalysis(tmpPath, rings, nList, yCloud, 7, 1);
  REQUIRE(ret == 0);
  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

// -- prism3 namespace tests --

TEST_CASE("basalPrismConditions returns true for connected basal rings",
          "[topo_bulk][prism3]") {
  // Build a prism: 8 atoms, two 4-membered basal rings connected
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalPrismConditions(nList, basal1, basal2);
  // Each atom in basal1 should have a neighbour in basal2
  REQUIRE(result == true);
}

TEST_CASE("basalPrismConditions returns false for disconnected rings",
          "[topo_bulk][prism3]") {
  // Two rings far apart with no bonds between them
  // Build nList manually: atoms 0-3 are in one plane, 4-7 in another far away
  // nList where atoms only know their own basal ring neighbours
  std::vector<std::vector<int>> nList = {
      {0, 1, 3}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalPrismConditions(nList, basal1, basal2);
  REQUIRE(result == false);
}

TEST_CASE("relaxedPrismConditions returns true when at least one bond exists",
          "[topo_bulk][prism3]") {
  // nList where only atom 0 has atom 4 as neighbour (one bond between rings)
  std::vector<std::vector<int>> nList = {
      {0, 1, 3, 4}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7, 0}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::relaxedPrismConditions(nList, basal1, basal2);
  REQUIRE(result == true);
}

TEST_CASE("relaxedPrismConditions returns false when no bonds exist",
          "[topo_bulk][prism3]") {
  std::vector<std::vector<int>> nList = {
      {0, 1, 3}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2},
      {4, 5, 7}, {5, 4, 6}, {6, 5, 7}, {7, 4, 6}};

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::relaxedPrismConditions(nList, basal1, basal2);
  REQUIRE(result == false);
}

TEST_CASE("basalRingsSeparation returns true for close rings",
          "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  // Two parallel 4-rings separated by 2.75 in z
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  // Height cutoff 8.0 (default), separation is 2.75
  bool result = prism3::basalRingsSeparation(cloud, basal1, basal2, 8.0);
  REQUIRE(result == true);
}

TEST_CASE("basalRingsSeparation returns false for far apart rings",
          "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {100, 100, 100};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  // Two rings separated by 20 in z
  double coords[8][3] = {{0, 0, 0},    {2, 0, 0},    {2, 2, 0},    {0, 2, 0},
                          {0, 0, 20},   {2, 0, 20},   {2, 2, 20},   {0, 2, 20}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};

  bool result = prism3::basalRingsSeparation(cloud, basal1, basal2, 8.0);
  REQUIRE(result == false);
}

TEST_CASE("findBulkPrisms on ideal prism geometry", "[topo_bulk][prism3]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10, 10, 50};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> pt;
  pt.type = 1;
  double coords[8][3] = {{0, 0, 0},       {2.75, 0, 0},    {2.75, 2.75, 0},
                          {0, 2.75, 0},    {0, 0, 2.75},    {2.75, 0, 2.75},
                          {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    pt.atomID = i;
    pt.x = coords[i][0]; pt.y = coords[i][1]; pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;

  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(nList, 5);

  // Get only 4-membered rings
  auto fourRings = ring::getSingleRingSize(rings, 4);
  if (fourRings.size() >= 2) {
    std::vector<ring::strucType> ringType(fourRings.size(),
                                           ring::strucType::unclassified);
    std::vector<double> rmsdPerAtom(cloud.nop, -1.0);

    int ret = prism3::findBulkPrisms(fourRings, ringType, nList, cloud,
                                      rmsdPerAtom, 8.0);
    REQUIRE(ret == 0);
  }
}

TEST_CASE("findPrismatic identifies prismatic rings between HC basals",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);

  std::vector<ring::strucType> ringType(rings.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(rings, ringType, nList, cageList);

  // The HC should have basal rings; find prismatic between them
  REQUIRE(cageList.size() == 1);
  if (cageList[0].rings.size() >= 2) {
    std::vector<int> prismaticRings;
    int ret = ring::findPrismatic(rings, listHC, ringType,
                                   cageList[0].rings[0], cageList[0].rings[1],
                                   prismaticRings);
    REQUIRE(ret == 0);
    // HC has 3 prismatic rings between its 2 basal rings
    REQUIRE(prismaticRings.size() == 3);
  }
}

TEST_CASE("DDC cage assembly is independent of ring enumeration order",
          "[topo_bulk]") {
  // The greedy assembly claims peripheral rings as it accepts cages, so
  // without a canonical visiting order the cage-object count depends on the
  // order the ring network was enumerated in. The per-ring classification
  // and the cage partition must both survive any permutation of the input.
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(idx, 7);
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }

  size_t referenceCages = 0;
  size_t referenceDDCRings = 0;
  for (int variant = 0; variant < 3; variant++) {
    auto shuffled = six;
    if (variant == 1) {
      std::sort(shuffled.begin(), shuffled.end());
    } else if (variant == 2) {
      std::sort(shuffled.rbegin(), shuffled.rend());
    }
    std::vector<ring::strucType> rt(shuffled.size(),
                                    ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto index = ring::buildRingSearchIndex(shuffled, yCloud);
    auto hc = ring::findHC(shuffled, rt, idx, cl, index);
    auto ddc = ring::findDDC(shuffled, rt, hc, cl, index);
    if (variant == 0) {
      referenceCages = cl.size();
      referenceDDCRings = ddc.size();
      REQUIRE(referenceDDCRings == 8192);
    } else {
      REQUIRE(cl.size() == referenceCages);
      REQUIRE(ddc.size() == referenceDDCRings);
    }
  }
}

TEST_CASE("Cage counts are invariant under atom relabelling",
          "[topo_bulk]") {
  // The frame is relabelled at random and its neighbour list and rings are
  // rebuilt. The cage count is the size of a maximum non-conflicting set,
  // so it does not depend on the labels. The particular packing is the one
  // selected by the atom-index key.
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  Cloud yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  const int n = yCloud.nop;

  using Cages = std::vector<std::pair<int, std::vector<std::vector<int>>>>;
  auto cagesOf = [](Cloud &cloud, const std::vector<int> &label) {
    auto nList = nneigh::neighListO(3.5, cloud, 1);
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    std::vector<std::vector<int>> six;
    for (const auto &r : primitive::ringNetwork(idx, 6)) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto hc = ring::findHC(six, rt, idx, cl, index);
    (void)ring::findDDC(six, rt, hc, cl, index);
    Cages out;
    for (const auto &c : cl) {
      std::vector<std::vector<int>> faces;
      for (const int ri : c.rings) {
        std::vector<int> face;
        for (const int a : six[ri]) {
          face.push_back(label[a]);
        }
        std::sort(face.begin(), face.end());
        faces.push_back(face);
      }
      std::sort(faces.begin(), faces.end());
      out.emplace_back(c.type == cage::cageType::HexC ? 0 : 1, faces);
    }
    std::sort(out.begin(), out.end());
    return out;
  };

  std::vector<int> identity(n);
  std::iota(identity.begin(), identity.end(), 0);
  const Cages reference = cagesOf(yCloud, identity);
  REQUIRE(!reference.empty());
  int refHc = 0;
  int refDdc = 0;
  for (const auto &c : reference) {
    (c.first == 0 ? refHc : refDdc)++;
  }

  for (const unsigned seed : {1U, 2U, 3U}) {
    std::vector<int> perm(identity);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);
    Cloud relabelled = yCloud;
    std::vector<int> back(n);
    for (int a = 0; a < n; a++) {
      relabelled.pts[perm[a]] = yCloud.pts[a];
      back[perm[a]] = a;
    }
    relabelled.idIndexMap.clear();
    for (int i = 0; i < n; i++) {
      relabelled.pts[i].atomID = i + 1;
      relabelled.idIndexMap[i + 1] = i;
    }
    const Cages got = cagesOf(relabelled, back);
    int hc = 0;
    int ddc = 0;
    for (const auto &c : got) {
      (c.first == 0 ? hc : ddc)++;
    }
    // Cardinality is the maximum packing. The particular cages follow the
    // atom-index key, so they are not the same set after a relabeling.
    REQUIRE(hc == refHc);
    REQUIRE(ddc == refDdc);
    REQUIRE(got.size() == reference.size());
  }
}

TEST_CASE("Cage objects are invariant under a periodic translation",
          "[topo_bulk]") {
  // Shifting the same frame by a fraction of the box and rewrapping does
  // not change the bonds. The cage set is a maximum packing keyed by atom
  // indices and the unwrapped shift, so it does not move with the origin.
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  auto cagesOf = [](Cloud &cloud) {
    auto nList = nneigh::neighListO(3.5, cloud, 1);
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    std::vector<std::vector<int>> six;
    for (const auto &r : primitive::ringNetwork(idx, 6)) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    const auto index = ring::buildRingSearchIndex(six, cloud);
    std::vector<ring::strucType> rt(six.size(), ring::strucType::unclassified);
    std::vector<cage::Cage> cl;
    const auto hc = ring::findHC(six, rt, idx, cl, index);
    (void)ring::findDDC(six, rt, hc, cl, index);
    std::vector<std::pair<int, std::vector<std::vector<int>>>> out;
    for (const auto &c : cl) {
      std::vector<std::vector<int>> faces;
      for (const int ri : c.rings) {
        auto face = six[static_cast<size_t>(ri)];
        std::sort(face.begin(), face.end());
        faces.push_back(std::move(face));
      }
      std::sort(faces.begin(), faces.end());
      out.emplace_back(c.type == cage::cageType::HexC ? 0 : 1, faces);
    }
    std::sort(out.begin(), out.end());
    return out;
  };
  Cloud cloud;
  cloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, cloud, 1);
  REQUIRE(cloud.nop == 4096);
  const auto reference = cagesOf(cloud);
  int ddc = 0;
  for (const auto &c : reference) {
    ddc += c.first == 1 ? 1 : 0;
  }
  // One equator for every two molecules: the conflict graph is the line
  // graph of a 4-regular bipartite graph on the molecules' rings, so the
  // maximum packing has n/2 cages.
  REQUIRE(ddc == 2048);
  const double frac[3] = {0.37, 0.21, 0.53};
  for (int k = 0; k < 3; k++) {
    const double lo = k < static_cast<int>(cloud.boxLow.size()) ? cloud.boxLow[static_cast<size_t>(k)] : 0.0;
    const double len = cloud.box[static_cast<size_t>(k)];
    for (int i = 0; i < cloud.nop; i++) {
      double *coord = k == 0 ? &cloud.pts[static_cast<size_t>(i)].x
                              : k == 1 ? &cloud.pts[static_cast<size_t>(i)].y
                                       : &cloud.pts[static_cast<size_t>(i)].z;
      double v = *coord - lo + frac[k] * len;
      v -= len * std::floor(v / len);
      *coord = lo + v;
    }
  }
  REQUIRE(cagesOf(cloud) == reference);
}

TEST_CASE("positionRank rejects two atoms at one wrapped position",
          "[topo_bulk]") {
  // An atom and the periodic image of another wrap onto one point, which no
  // position order separates, so the frame is rejected
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.resize(3);
  const double xs[3] = {2.0, 1.0, 3.0};
  for (int i = 0; i < 3; i++) {
    cloud.pts[i].x = xs[i];
    cloud.pts[i].y = 1.0;
    cloud.pts[i].z = 1.0;
  }
  cloud.nop = 3;
  REQUIRE(ring::positionRank(cloud) == std::vector<int>{1, 0, 2});
  cloud.pts[2].x = 12.0;
  REQUIRE_THROWS_AS(ring::positionRank(cloud), std::invalid_argument);
}

namespace {

std::vector<ring::EquatorChoice>
conflictChoices(int n, const std::vector<std::pair<int, int>> &edges) {
  std::vector<std::vector<int>> adj(static_cast<size_t>(n));
  for (const auto &e : edges) {
    adj[static_cast<size_t>(e.first)].push_back(e.second);
    adj[static_cast<size_t>(e.second)].push_back(e.first);
  }
  std::vector<ring::EquatorChoice> periph(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &row = periph[static_cast<size_t>(i)];
    if (adj[static_cast<size_t>(i)].empty()) {
      row.n = 1;
      row.ring[0] = -1;
      continue;
    }
    row.n = static_cast<int>(adj[static_cast<size_t>(i)].size());
    for (int k = 0; k < row.n; k++) {
      row.ring[k] = adj[static_cast<size_t>(i)][static_cast<size_t>(k)];
    }
  }
  return periph;
}

int bruteIndependent(int n, const std::vector<std::pair<int, int>> &edges) {
  int best = 0;
  const int limit = 1 << n;
  for (int mask = 0; mask < limit; mask++) {
    bool ok = true;
    for (const auto &e : edges) {
      if ((mask & (1 << e.first)) && (mask & (1 << e.second))) {
        ok = false;
        break;
      }
    }
    if (ok) {
      best = std::max(best, __builtin_popcount(mask));
    }
  }
  return best;
}

} // namespace

TEST_CASE("equatorial search is exact on K4 minus an edge", "[topo_bulk]") {
  for (int a = 0; a < 4; a++) {
    for (int b = a + 1; b < 4; b++) {
      std::vector<std::pair<int, int>> edges;
      for (int u = 0; u < 4; u++) {
        for (int v = u + 1; v < 4; v++) {
          if (u == a && v == b) {
            continue;
          }
          edges.push_back({u, v});
        }
      }
      const std::vector<int> order{0, 1, 2, 3};
      ring::PackingCertificate cert;
      const auto kept = ring::maximumEquators(conflictChoices(4, edges), order, &cert);
      REQUIRE(static_cast<int>(kept.size()) == 2);
      REQUIRE(cert.method == 2);
    }
  }
}

TEST_CASE("equatorial search matches brute force on random graphs", "[topo_bulk]") {
  std::mt19937 rng(20261102);
  std::uniform_real_distribution<double> coin(0.0, 1.0);
  int checked = 0;
  for (int n = 4; n <= 8; n++) {
    for (int sample = 0; sample < 20; sample++) {
      std::vector<std::pair<int, int>> edges;
      std::vector<int> deg(static_cast<size_t>(n), 0);
      for (int u = 0; u < n; u++) {
        for (int v = u + 1; v < n; v++) {
          if (coin(rng) < 0.45 && deg[static_cast<size_t>(u)] < 6 &&
              deg[static_cast<size_t>(v)] < 6) {
            edges.push_back({u, v});
            deg[static_cast<size_t>(u)]++;
            deg[static_cast<size_t>(v)]++;
          }
        }
      }
      std::vector<int> order(static_cast<size_t>(n));
      std::iota(order.begin(), order.end(), 0);
      const auto kept = ring::maximumEquators(conflictChoices(n, edges), order, nullptr);
      REQUIRE(static_cast<int>(kept.size()) == bruteIndependent(n, edges));
      checked++;
    }
  }
  REQUIRE(checked == 100);
}

TEST_CASE("a wide parity product searches the large equator", "[topo]") {
  const double xyz[][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073},
  };
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {40.0, 40.0, 40.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  molSys::Point<double> point;
  point.type = 1;
  for (int i = 0; i < 14; ++i) {
    point.atomID = i;
    point.x = xyz[i][0];
    point.y = xyz[i][1];
    point.z = xyz[i][2];
    cloud.pts.push_back(point);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 14;
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  rows = nneigh::neighbourListByIndex(cloud, rows);
  auto rings = primitive::ringNetwork(rows, 6);
  const int base = static_cast<int>(rings.size());
  REQUIRE(base > 6);
  // Eight copies of every ring put more than sixteen candidates on a
  // triplet, so the product leaves the triple loop.
  rings.reserve(static_cast<std::size_t>(base * 9));
  for (int copy = 0; copy < 8; ++copy) {
    for (int i = 0; i < base; ++i) {
      rings.push_back(rings[static_cast<std::size_t>(i)]);
    }
  }
  auto index = ring::buildRingSearchIndex(rings, 14);
  std::vector<int> peripheral;
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  ring::EquatorChoice raw;
  bool any = false;
  for (int iring = 0; iring < base; ++iring) {
    any = ring::equatorialConditions(rings, iring, index, &peripheral, &chosen,
                                     &rank, &raw) ||
          any;
    any = ring::equatorialConditions(rings, iring, index, nullptr) || any;
  }
  // Equal ranks make the first live even chain the pivot, so the even
  // side of the pair search runs as well as the odd side.
  std::fill(rank.begin(), rank.end(), 0);
  for (int iring = 0; iring < base; ++iring) {
    any = ring::equatorialConditions(rings, iring, index, &peripheral, &chosen,
                                     &rank, &raw) ||
          any;
  }
  std::vector<ring::strucType> ringType(rings.size(), ring::strucType::unclassified);
  std::vector<int> listHC;
  std::vector<cage::Cage> cages;
  setenv("DSEAMS_CAGE_TIME", "1", 1);
  const auto ddc = ring::findDDC(rings, ringType, listHC, cages, index);
  unsetenv("DSEAMS_CAGE_TIME");
  REQUIRE_FALSE(ddc.empty());
  REQUIRE(any);
  const auto split = ring::lastDdcSplit();
  REQUIRE(split.equator_ms >= 0.0);
  std::vector<int> peripheralAgain;
  for (int iring = 0; iring < base && iring < static_cast<int>(rings.size()); ++iring) {
    (void)ring::conditionOneDDC(rings, peripheralAgain, iring, index);
  }
  if (base >= 2) {
    const auto &b1 = rings[0];
    const auto &b2 = rings[1];
    (void)ring::basalConditions(rows, b1, b2);
    ring::NeighBond bonds;
    bonds.build(rows);
    (void)ring::basalConditions(bonds, b1, b2);
    std::vector<int> triplet = {b1[0], b1[1], b1[2]};
    (void)ring::basalNeighbours(rows, triplet, b2[0], b2[1]);
    (void)ring::notNeighboursOfRing(rows, triplet, b2);
  }
  std::vector<cage::iceType> atomTypes(14, cage::iceType::dummy);
  (void)ring::getAtomTypesTopoBulk(rings, ringType, atomTypes);
  int numHC = 0, numDDC = 0, mixed = 0, prism = 0, basal = 0;
  (void)ring::getStrucNumbers(ringType, cages, numHC, numDDC, mixed, prism, basal);
  std::vector<int> peripheralBare;
  (void)ring::conditionOneDDC(rings, peripheralBare, 0);
  std::vector<int> listDDC = ddc;
  std::vector<int> listHC2;
  (void)ring::findMixedRings(rings, ringType, listDDC, listHC2);
}

TEST_CASE("constructed rings hit the bulk type and neighbour branches", "[topo]") {
  const double xyz[][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073},
  };
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {40.0, 40.0, 40.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  molSys::Point<double> point;
  point.type = 1;
  for (int i = 0; i < 14; ++i) {
    point.atomID = i;
    point.x = xyz[i][0];
    point.y = xyz[i][1];
    point.z = xyz[i][2];
    cloud.pts.push_back(point);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 14;
  auto idRows = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, idRows);
  auto rings = primitive::ringNetwork(rows, 8);
  REQUIRE(rings.size() > 4);
  auto index = ring::buildRingSearchIndex(rings, 14);
  std::vector<ring::strucType> ringType(rings.size(), ring::strucType::unclassified);
  std::vector<int> listHC;
  std::vector<cage::Cage> cages;
  auto ddc = ring::findDDC(rings, ringType, listHC, cages, index);
  auto hc = ring::findHC(rings, ringType, rows, cages, index);
  std::vector<int> prismatic;
  if (rings.size() > 1) {
    (void)ring::findPrismatic(rings, listHC, ringType, 0, 1, prismatic, index);
  }
  (void)ring::findMixedRings(rings, ringType, ddc, hc);
  std::vector<double> rmsd(14, 0.0);
  (void)prism3::findBulkPrisms(rings, ringType, rows, cloud, rmsd, 2.0);

  std::vector<std::vector<int>> typed = {
      {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}, {1, 2, 3, 4, 5, 0},
      {0, 1, 2, 3, 4, 5}, {2, 3, 4, 5, 0, 1}, {0, 1, 2, 3, 4, 5},
      {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5},
  };
  std::vector<ring::strucType> kinds = {
      ring::strucType::unclassified, ring::strucType::DDC,
      ring::strucType::HCbasal,      ring::strucType::HCprismatic,
      ring::strucType::bothBasal,    ring::strucType::bothPrismatic,
      ring::strucType::Prism,        ring::strucType::deformedPrism,
      ring::strucType::mixedPrismRing,
  };
  std::vector<cage::iceType> atoms(6, cage::iceType::dummy);
  atoms[0] = cage::iceType::pnc;
  atoms[1] = cage::iceType::mixed;
  atoms[2] = cage::iceType::mixed2;
  REQUIRE(ring::getAtomTypesTopoBulk(typed, kinds, atoms) == 0);
  REQUIRE(atoms[0] == cage::iceType::mixed2);

  std::vector<cage::Cage> cageList(2);
  cageList[0].type = cage::cageType::HexC;
  cageList[1].type = cage::cageType::DoubleDiaC;
  int numHC = 0, numDDC = 0, mixed = 0, prism = 0, basal = 0;
  REQUIRE(ring::getStrucNumbers(kinds, cageList, numHC, numDDC, mixed, prism, basal) == 0);
  REQUIRE(numHC == 1);
  REQUIRE(numDDC == 1);
  REQUIRE(mixed == 2);

  std::vector<std::vector<int>> nList = {
      {0, 1, 2}, {1, 0, 3}, {2, 0, 3}, {3, 1, 2}, {4},
  };
  std::vector<int> triplet = {0, 1, 2};
  REQUIRE(ring::basalNeighbours(nList, triplet, 0, 3));
  REQUIRE_FALSE(ring::basalNeighbours(nList, triplet, 1, 2));
  triplet = {0, 1, 8};
  REQUIRE_FALSE(ring::basalNeighbours(nList, triplet, 0, 3));
  triplet = {0, 3, 8};
  REQUIRE_FALSE(ring::basalNeighbours(nList, triplet, 1, 0));
  triplet = {0, 1, 2};
  REQUIRE(ring::notNeighboursOfRing(nList, triplet, std::vector<int>{4}));
  REQUIRE_FALSE(ring::notNeighboursOfRing(nList, triplet, std::vector<int>{0}));

  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};
  std::vector<std::vector<int>> prismList(8);
  for (int i = 0; i < 4; ++i) {
    prismList[i] = {i, basal2[i]};
    prismList[basal2[i]] = {basal2[i], i};
  }
  (void)prism3::basalPrismConditions(prismList, basal1, basal2);
  (void)prism3::basalPrismConditions(nList, basal1, basal1);

  std::vector<int> bare;
  (void)ring::conditionOneDDC(rings, bare, 0);
  rings.push_back({-1, 0, 1, 2, 3, 4});
  (void)ring::conditionOneDDC(rings, bare, static_cast<int>(rings.size()) - 1);
  auto wideIndex = ring::buildRingSearchIndex(rings, 14);
  (void)ring::conditionOneDDC(rings, bare, static_cast<int>(rings.size()) - 1, wideIndex);
  std::vector<int> peripheral;
  (void)ring::equatorialConditions(rings, 0, wideIndex, &peripheral);
  (void)ring::findMixedRings(typed, kinds, ddc, hc);

  // Many small edits of one real ring network. Each edit is a lattice the
  // search can refuse, and together they walk the remaining branches.
  std::vector<std::vector<int>> baseRings = rings;
  if (!baseRings.empty() && baseRings[0].size() == 6) {
    std::vector<std::vector<int>> variants;
    variants.push_back(baseRings[0]);
    auto shorter = baseRings[0];
    shorter.pop_back();
    variants.push_back(shorter);
    auto longer = baseRings[0];
    longer.push_back(longer.back());
    variants.push_back(longer);
    auto negative = baseRings[0];
    negative[0] = -1;
    variants.push_back(negative);
    auto far = baseRings[0];
    far[0] = 40;
    variants.push_back(far);
    auto reversed = baseRings[0];
    std::reverse(reversed.begin(), reversed.end());
    variants.push_back(reversed);
    for (const auto &extra : variants) {
      auto trial = baseRings;
      trial.push_back(extra);
      trial.push_back(baseRings[0]);
      auto trialType = std::vector<ring::strucType>(trial.size(), ring::strucType::unclassified);
      std::vector<int> trialHC;
      std::vector<cage::Cage> trialCages;
      auto trialIndex = ring::buildRingSearchIndex(trial, 48);
      auto trialDDC = ring::findDDC(trial, trialType, trialHC, trialCages, trialIndex);
      auto trialHCFound = ring::findHC(trial, trialType, rows, trialCages, trialIndex);
      std::vector<int> peripheralTrial;
      for (int i = 0; i < std::min(4, static_cast<int>(trial.size())); ++i) {
        (void)ring::equatorialConditions(trial, i, trialIndex, &peripheralTrial);
        (void)ring::conditionOneDDC(trial, peripheralTrial, i, trialIndex);
      }
      bool hexOnly = true;
      for (const auto &ring : trial) {
        if (ring.size() != 6) {
          hexOnly = false;
        }
      }
      if (hexOnly && trial.size() > 1) {
        std::vector<int> prisms;
        (void)ring::findPrismatic(trial, trialHC, trialType, 0, 1, prisms, trialIndex);
      }
      (void)ring::findMixedRings(trial, trialType, trialDDC, trialHCFound);
    }
  }
  std::vector<std::vector<int>> unequal = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5, 6}, {0}};
  (void)ring::buildRingSearchIndex(unequal, 8);
  std::vector<ring::strucType> unequalType(unequal.size(), ring::strucType::unclassified);
  std::vector<int> unequalHC;
  std::vector<cage::Cage> unequalCages;
  (void)ring::findDDC(unequal, unequalType, unequalHC, unequalCages);
  (void)ring::findHC(unequal, unequalType, rows, unequalCages);
}

TEST_CASE("basal bonds, a non-hex atom type and a failed prism match",
          "[topo_bulk]") {
  std::vector<std::vector<int>> nList(12);
  for (int i = 0; i < 12; ++i) {
    nList[static_cast<std::size_t>(i)] = {i};
  }
  nList[0].push_back(6);
  nList[2].push_back(8);
  nList[4].push_back(10);
  const std::vector<int> basal1 = {0, 1, 2, 3, 4, 5};
  const std::vector<int> basal2 = {6, 7, 8, 9, 10, 11};
  REQUIRE(ring::basalConditions(nList, basal1, basal2));
  ring::NeighBond bonds;
  bonds.build(nList);
  REQUIRE(ring::basalConditions(bonds, basal1, basal2));
  REQUIRE_FALSE(ring::basalConditions(nList, basal1, basal1));

  std::vector<std::vector<int>> neigh = {{0, 2}, {1}, {2}, {3, 1}};
  std::vector<int> triplet = {9, 1, 2};
  REQUIRE(ring::basalNeighbours(neigh, triplet, 0, 3));
  triplet = {9, 1, 8};
  REQUIRE_FALSE(ring::basalNeighbours(neigh, triplet, 0, 3));

  std::vector<std::vector<int>> quads = {{0, 1, 2, 3}};
  std::vector<ring::strucType> quadType = {ring::strucType::Prism};
  std::vector<cage::iceType> atoms(4, cage::iceType::dummy);
  REQUIRE(ring::getAtomTypesTopoBulk(quads, quadType, atoms) == 0);
  REQUIRE(atoms[0] == cage::iceType::pnc);

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {40, 10, 20};
  cloud.boxLow = {0, 0, 0};
  cloud.currentFrame = 1;
  double coords[8][3] = {{0, 0, 0},       {20, 0, 0},      {20, 1, 0},
                         {0, 1, 0},       {0, 0, 2.5},     {2.5, 0, 2.5},
                         {2.5, 2.5, 2.5}, {0, 2.5, 2.5}};
  for (int i = 0; i < 8; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;
  std::vector<std::vector<int>> faces = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<std::vector<int>> faceRows(8);
  for (int i = 0; i < 8; ++i) {
    faceRows[static_cast<std::size_t>(i)] = {i};
  }
  std::vector<ring::strucType> faceType(2, ring::strucType::unclassified);
  std::vector<double> rmsd(8, -1.0);
  REQUIRE(prism3::findBulkPrisms(faces, faceType, faceRows, cloud, rmsd, 8.0) ==
          0);
}

TEST_CASE("repeated triplets exercise the large equatorial search",
          "[topo_bulk]") {
  constexpr int copies = 7;
  std::vector<std::vector<int>> rings;
  rings.push_back({0, 1, 2, 0, 1, 2});
  for (int i = 0; i < copies; ++i) {
    rings.push_back({0, 1, 2, 10 + i, 30 + i, 50 + i});
  }
  const int nAtoms = 60;
  const auto index = ring::buildRingSearchIndex(rings, nAtoms);
  std::vector<int> peripheral;
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  ring::EquatorChoice raw;
  (void)ring::equatorialConditions(rings, 0, index, &peripheral, &chosen, &rank,
                                   &raw);

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);
  auto listed = nneigh::neighListO(3.5, yCloud, 1);
  auto byIndex = nneigh::neighbourListByIndex(yCloud, listed);
  auto six = ring::getSingleRingSize(primitive::ringNetwork(byIndex, 6), 6);
  REQUIRE_FALSE(six.empty());
  std::vector<ring::strucType> typed(six.size(), ring::strucType::DDC);
  std::vector<cage::Cage> cages;
  const auto hc = ring::findHC(six, typed, byIndex, cages);
  REQUIRE_FALSE(hc.empty());
  REQUIRE(std::find(typed.begin(), typed.end(), ring::strucType::bothBasal) !=
          typed.end());

  const auto outDir = fs::temp_directory_path() / "dseams_topo_penta";
  fs::create_directories(outDir);
  auto pentagons = six;
  if (!pentagons.empty()) {
    pentagons.push_back({0, 1, 2, 3, 4});
  }
  REQUIRE(ring::topoBulkAnalysis(outDir.string(), pentagons, byIndex, yCloud, 1,
                                 false) == 0);
  fs::remove_all(outDir);

  molSys::PointCloud<molSys::Point<double>, double> prism;
  prism.box = {10, 10, 10};
  prism.boxLow = {0, 0, 0};
  prism.currentFrame = 1;
  double xyz[8][3] = {{0, 0, 0}, {2.75, 0, 0}, {2.75, 2.75, 0}, {0, 2.75, 0},
                      {0, 0, 2.75}, {2.75, 0, 2.75}, {2.75, 2.75, 2.75},
                      {0, 2.75, 2.75}};
  for (int i = 0; i < 8; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = xyz[i][0];
    pt.y = xyz[i][1];
    pt.z = xyz[i][2];
    prism.pts.push_back(pt);
    prism.idIndexMap[i] = i;
  }
  prism.nop = 8;
  std::vector<std::vector<int>> squares = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<std::vector<int>> squareRows(8);
  for (int i = 0; i < 8; ++i) {
    squareRows[static_cast<std::size_t>(i)] = {i, (i + 1) % 4 + (i < 4 ? 0 : 4)};
  }
  std::vector<ring::strucType> squareType(2, ring::strucType::unclassified);
  std::vector<double> squareRmsd(8, -1.0);
  REQUIRE(prism3::findBulkPrisms(squares, squareType, squareRows, prism,
                                 squareRmsd, 8.0) == 0);
}

TEST_CASE("small cages reach the live bulk search returns", "[topo_bulk]") {
  auto cloudOf = [](const double (*xyz)[3], int n, double box) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud.box = {box, box, box};
    cloud.boxLow = {0, 0, 0};
    cloud.currentFrame = 1;
    for (int i = 0; i < n; ++i) {
      molSys::Point<double> pt;
      pt.type = 1;
      pt.atomID = i;
      pt.x = xyz[i][0];
      pt.y = xyz[i][1];
      pt.z = xyz[i][2];
      cloud.pts.push_back(pt);
      cloud.idIndexMap[i] = i;
    }
    cloud.nop = n;
    return cloud;
  };
  auto sixOf = [](molSys::PointCloud<molSys::Point<double>, double> &cloud) {
    auto listed = nneigh::neighListO(3.5, cloud, 1);
    auto rows = nneigh::neighbourListByIndex(cloud, listed);
    auto rings = ring::getSingleRingSize(primitive::ringNetwork(rows, 6), 6);
    return std::pair{rows, rings};
  };

  const double hc[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  auto hcCloud = cloudOf(hc, 12, 40.0);
  auto [hcRows, hcRings] = sixOf(hcCloud);
  REQUIRE_FALSE(hcRings.empty());

  auto shifted = hcRings;
  for (auto &ring : shifted) {
    for (int &atom : ring) {
      atom += 12;
    }
  }
  auto both = hcRings;
  both.insert(both.end(), shifted.begin(), shifted.end());
  std::vector<std::vector<int>> bothRows = hcRows;
  bothRows.resize(24);
  for (int i = 0; i < static_cast<int>(hcRows.size()); ++i) {
    auto row = hcRows[static_cast<std::size_t>(i)];
    for (int &atom : row) {
      atom += 12;
    }
    bothRows[static_cast<std::size_t>(i + 12)] = std::move(row);
  }
  both.push_back(both.front());
  std::vector<ring::strucType> bothType(both.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> bothCages;
#ifdef SEAMS_HAS_OPENMP
  const int threadsBefore = omp_get_max_threads();
  omp_set_num_threads(8);
#endif
  const auto twoHc = ring::findHC(both, bothType, bothRows, bothCages);
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(1);
#endif
  REQUIRE_FALSE(twoHc.empty());

  std::vector<std::vector<int>> shareRows(9);
  for (int i = 0; i < 9; ++i) {
    shareRows[static_cast<std::size_t>(i)] = {i};
    for (int j = 0; j < 9; ++j) {
      if (j != i) {
        shareRows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
  }
  std::vector<std::vector<int>> shareRings = {
      {0, 1, 2, 3, 4, 5},
      {0, 1, 2, 6, 7, 8},
      {},
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7},
      {0, 1, 2, 3, 4, 99},
  };
  std::vector<ring::strucType> shareType(shareRings.size(),
                                         ring::strucType::unclassified);
  std::vector<cage::Cage> shareCages;
  (void)ring::findHC(shareRings, shareType, shareRows, shareCages);

  std::vector<std::vector<int>> orderRings = {{0, 1, 2, 3, 4, 5},
                                              {0, 2, 1, 3, 4, 5},
                                              {0, 1, 3, 2, 4, 5}};
  std::vector<std::vector<int>> orderRows(6);
  for (int i = 0; i < 6; ++i) {
    orderRows[static_cast<std::size_t>(i)] = {i, (i + 1) % 6};
  }
  std::vector<ring::strucType> orderType(orderRings.size(),
                                         ring::strucType::unclassified);
  std::vector<cage::Cage> orderCages;
  (void)ring::findHC(orderRings, orderType, orderRows, orderCages);

  std::vector<std::vector<int>> shortList(6);
  for (int i = 0; i < 6; ++i) {
    shortList[static_cast<std::size_t>(i)] = {i};
  }
  REQUIRE_FALSE(ring::basalConditions(shortList, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}));
  REQUIRE_FALSE(ring::basalConditions(shortList, {-1, 0, 1, 2, 3, 4},
                                      {0, 1, 2, 3, 4, 5}));
  ring::NeighBond bare;
  REQUIRE_FALSE(bare(0, 1));
  ring::NeighBond built;
  built.build(shortList);
  REQUIRE_FALSE(built(-1, 0));
  (void)ring::buildRingSearchIndex(shareRings, 0);

  std::vector<std::vector<int>> one = {{0, 1, 2, 3, 4, 5}};
  std::vector<ring::strucType> mixedType = {ring::strucType::bothBasal};
  std::vector<int> listD = {0};
  std::vector<int> listH = {0};
  const auto mixed = ring::findMixedRings(one, mixedType, listD, listH);
  REQUIRE(mixed.size() == 1);
  REQUIRE(listH[0] == -10);

  const double ddc[14][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073}};
  auto ddcCloud = cloudOf(ddc, 14, 54.0);
  auto [ddcRows, ddcRings] = sixOf(ddcCloud);
  REQUIRE_FALSE(ddcRings.empty());
  auto ddcIndex = ring::buildRingSearchIndex(ddcRings, 14);
  std::vector<int> rank(ddcRings.size());
  std::iota(rank.begin(), rank.end(), 0);
  bool chose = false;
  int equator = -1;
  for (int i = 0; i < static_cast<int>(ddcRings.size()) && !chose; ++i) {
    std::vector<int> peripheral;
    std::vector<int> chosen;
    chose = ring::equatorialConditions(ddcRings, i, ddcIndex, &peripheral, &chosen,
                                       &rank, nullptr);
    if (chose) {
      equator = i;
    }
  }
  REQUIRE(chose);
  auto mutated = ddcRings;
  mutated.push_back(ddcRings[static_cast<std::size_t>((equator + 1) % ddcRings.size())]);
  auto mutIndex = ring::buildRingSearchIndex(mutated, 14);
  std::vector<int> mutRank(mutated.size());
  std::iota(mutRank.begin(), mutRank.end(), 0);
  std::vector<int> mutChosen;
  std::vector<int> mutPeri;
  (void)ring::equatorialConditions(mutated, equator, mutIndex, &mutPeri, &mutChosen,
                                   &mutRank, nullptr);
  (void)ring::equatorialConditions(mutated, equator, mutIndex, nullptr);
  if (equator >= 0) {
    auto shortened = ddcRings;
    const int victim = (equator + 1) % static_cast<int>(shortened.size());
    if (shortened[static_cast<std::size_t>(victim)].size() > 3) {
      shortened[static_cast<std::size_t>(victim)].pop_back();
    }
    auto shortIndex = ring::buildRingSearchIndex(shortened, 14);
    std::vector<int> sink;
    (void)ring::equatorialConditions(shortened, equator, shortIndex, &sink);
  }
  std::vector<ring::strucType> ddcType(ddcRings.size(), ring::strucType::unclassified);
  if (ddcRings.size() > 1) {
    ddcType[1] = ring::strucType::HCprismatic;
  }
  std::vector<int> ddcHc;
  std::vector<cage::Cage> ddcCages;
  (void)ring::findDDC(ddcRings, ddcType, ddcHc, ddcCages, ddcIndex);

  std::vector<std::vector<int>> broken = {{0, 1, 2, 3, 4, 5}};
  for (int k = 0; k < 6; ++k) {
    if (k == 0) {
      broken.push_back({0, 1, 2, 6, 7});
    } else {
      broken.push_back({k % 6, (k + 1) % 6, (k + 2) % 6, 10 + 3 * k, 11 + 3 * k,
                        12 + 3 * k});
    }
  }
  auto brokenIndex = ring::buildRingSearchIndex(broken, 40);
  std::vector<int> ignored;
  (void)ring::equatorialConditions(broken, 0, brokenIndex, &ignored);
  broken[1] = {0, 1, 2, 6, 7, 8};
  brokenIndex = ring::buildRingSearchIndex(broken, 40);
  (void)ring::equatorialConditions(broken, 0, brokenIndex, &ignored);

  molSys::PointCloud<molSys::Point<double>, double> stretched;
  stretched.box = {200, 20, 20};
  stretched.boxLow = {0, 0, 0};
  stretched.currentFrame = 1;
  const double face[8][3] = {{0, 0, 0},       {2.5, 0, 0},       {2.5, 2.5, 0},     {0, 2.5, 0},
                             {80, 0, 2.0},   {82.5, 0, 2.0},    {82.5, 2.5, 2.0},  {80, 2.5, 2.0}};
  for (int i = 0; i < 8; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = face[i][0];
    pt.y = face[i][1];
    pt.z = face[i][2];
    stretched.pts.push_back(pt);
    stretched.idIndexMap[i] = i;
  }
  stretched.nop = 8;
  std::vector<std::vector<int>> faces = {{0, 1, 2, 3}, {4, 5, 6, 7}};
  std::vector<std::vector<int>> faceRows(8);
  for (int i = 0; i < 4; ++i) {
    faceRows[static_cast<std::size_t>(i)] = {i, 4, 5, 6, 7};
    faceRows[static_cast<std::size_t>(i + 4)] = {i + 4, 0, 1, 2, 3};
  }
  std::vector<int> lower = {0, 1, 2, 3};
  std::vector<int> upper = {4, 5, 6, 7};
  REQUIRE(prism3::basalPrismConditions(faceRows, lower, upper));
  std::vector<ring::strucType> faceType(2, ring::strucType::unclassified);
  std::vector<double> faceRmsd(8, -1.0);
  REQUIRE(prism3::findBulkPrisms(faces, faceType, faceRows, stretched, faceRmsd,
                                 100.0) == 0);
  REQUIRE(faceType[0] == ring::strucType::Prism);
  REQUIRE(faceType[1] == ring::strucType::Prism);

  const auto outDir = fs::temp_directory_path() / "dseams_topo_hex_only";
  fs::create_directories(outDir);
  (void)ring::topoBulkAnalysis(outDir.string(), hcRings, hcRows, hcCloud, 1, false);
  fs::remove_all(outDir);
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(threadsBefore);
#endif
}

TEST_CASE("wide equator skips a short ring and an HC basal also closes a DDC",
          "[topo_bulk]") {
  std::vector<std::vector<int>> rings;
  rings.push_back({0, 1, 2, 3, 4, 5});
  for (int i = 0; i < 8; ++i) {
    rings.push_back({0, 1, 2, 10 + i, 30 + i, 50 + i});
    rings.push_back({2, 3, 4, 70 + i, 90 + i, 110 + i});
    rings.push_back({4, 5, 0, 130 + i, 150 + i, 170 + i});
  }
  rings.push_back({0, 1, 2, 3, 4});
  rings.push_back({0, 1, 2, -1, 8, 9});
  rings.push_back({0, 1, 2, 400, 401, 402});
  const int nAtoms = 420;
  const auto index = ring::buildRingSearchIndex(rings, nAtoms);
  std::vector<int> peripheral;
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  ring::EquatorChoice raw;
  (void)ring::equatorialConditions(rings, 0, index, &peripheral, &chosen, &rank,
                                   &raw);

  const double xyz[][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073},
  };
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {40.0, 40.0, 40.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  molSys::Point<double> point;
  point.type = 1;
  for (int i = 0; i < 14; ++i) {
    point.atomID = i;
    point.x = xyz[i][0];
    point.y = xyz[i][1];
    point.z = xyz[i][2];
    cloud.pts.push_back(point);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 14;
  auto idRows = nneigh::neighListO(3.5, cloud, 1);
  auto rows = nneigh::neighbourListByIndex(cloud, idRows);
  auto ddcRings = primitive::ringNetwork(rows, 8);
  REQUIRE(ddcRings.size() > 4);
  auto ddcIndex = ring::buildRingSearchIndex(ddcRings, 14);
  std::vector<ring::strucType> marked(ddcRings.size(), ring::strucType::HCbasal);
  std::vector<int> listHC;
  std::vector<cage::Cage> cages;
  auto ddc = ring::findDDC(ddcRings, marked, listHC, cages, ddcIndex);
  REQUIRE_FALSE(ddc.empty());
  REQUIRE(std::find(marked.begin(), marked.end(), ring::strucType::bothBasal) !=
          marked.end());

  // Seven copies on the even triplets push the parity product into the
  // wide search. A short ring, a negative atom and peripherals that do
  // not share an apex are inputs that search refuses.
  std::vector<std::vector<int>> wide;
  wide.push_back({0, 1, 2, 3, 4, 5});
  auto face = [](int a, int b, int c, int x, int y, int z) {
    return std::vector<int>{a, b, c, x, y, z};
  };
  for (int i = 0; i < 7; ++i) {
    wide.push_back(face(0, 1, 2, 10, 11, 12));
    wide.push_back(face(2, 3, 4, 10, 11, 12));
    wide.push_back(face(4, 5, 0, 10, 11, 12));
    wide.push_back(face(1, 2, 3, 20, 21, 22));
    wide.push_back(face(3, 4, 5, 20, 21, 22));
    wide.push_back(face(5, 0, 1, 20, 21, 22));
  }
  wide.push_back(face(0, 1, 2, -1, 11, 12));
  wide.push_back(face(2, 3, 4, -1, 11, 12));
  wide.push_back(face(4, 5, 0, -1, 11, 12));
  wide.push_back({0, 1, 2, 7, 8});
  const auto wideIndex = ring::buildRingSearchIndex(wide, 800);
  std::vector<int> widePeripheral;
  std::vector<int> wideChosen;
  std::vector<int> wideRank(wide.size());
  std::iota(wideRank.begin(), wideRank.end(), 0);
  ring::EquatorChoice wideRaw;
  (void)ring::equatorialConditions(wide, 0, wideIndex, &widePeripheral,
                                   &wideChosen, &wideRank, &wideRaw);

  std::vector<std::vector<int>> oneEach = {
      {0, 1, 2, 3, 4, 5}, {0, 1, 2, 50, 51, 52}, {1, 2, 3, 60, 61, 62},
      {2, 3, 4, 50, 51, 52}, {3, 4, 5, 60, 61, 62}, {4, 5, 0, 50, 51, 52},
      {5, 0, 1, 60, 61, 62},
  };
  const auto oneIndex = ring::buildRingSearchIndex(oneEach, 40);
  std::vector<int> onePeripheral;
  (void)ring::equatorialConditions(oneEach, 0, oneIndex, &onePeripheral);

  std::vector<std::vector<int>> shifted = {{0, 1, 2, 3, 4, 5},
                                           {0, 1, 2, 3, 4, 5}};
  auto shiftIndex = ring::buildRingSearchIndex(shifted, 6);
  shiftIndex.ringShift = {{0, 0, 0}, {1, 0, 0}};
  std::vector<ring::strucType> shiftType(2, ring::strucType::unclassified);
  std::vector<int> shiftHC;
  std::vector<cage::Cage> shiftCages;
  (void)ring::findDDC(shifted, shiftType, shiftHC, shiftCages, shiftIndex);

  std::vector<std::vector<int>> oddRings = {{-1, 0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9, 10}};
  std::vector<std::vector<int>> oddRows(11);
  for (int i = 0; i < 11; ++i) {
    oddRows[static_cast<std::size_t>(i)] = {i};
  }
  oddRows[0] = {0, 5, 6, 7};
  oddRows[5] = {5, 0};
  std::vector<ring::strucType> oddType(2, ring::strucType::unclassified);
  std::vector<cage::Cage> oddCages;
  (void)ring::findHC(oddRings, oddType, oddRows, oddCages);

  std::vector<std::vector<int>> typed = {{0, 1, 2, 3, 4, 5}};
  std::vector<ring::strucType> kinds = {static_cast<ring::strucType>(99)};
  std::vector<cage::iceType> atoms(6, cage::iceType::dummy);
  (void)ring::getAtomTypesTopoBulk(typed, kinds, atoms);

  std::vector<std::vector<int>> prismRows(8);
  for (int i = 0; i < 8; ++i) {
    prismRows[static_cast<std::size_t>(i)] = {i};
  }
  prismRows[0] = {0, 4};
  std::vector<int> basalA = {0, 1, 2, 3, 4, 5};
  std::vector<int> basalB = {4, 5, 6, 7, 8, 9};
  prismRows.resize(10);
  for (int i = 8; i < 10; ++i) {
    prismRows[static_cast<std::size_t>(i)] = {i};
  }
  REQUIRE_FALSE(prism3::basalPrismConditions(prismRows, basalA, basalB));

  molSys::PointCloud<molSys::Point<double>, double> flat;
  flat.box = {40, 40, 10};
  flat.boxLow = {0, 0, 0};
  flat.currentFrame = 1;
  for (int i = 0; i < 12; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    const double ang = (i % 6) * 1.0471975512;
    pt.x = 8.0 * (i / 6) + 2.0 * std::cos(ang);
    pt.y = 2.0 * std::sin(ang);
    pt.z = 0.0;
    flat.pts.push_back(pt);
    flat.idIndexMap[i] = i;
  }
  flat.nop = 12;
  std::vector<std::vector<int>> flatRows(12);
  for (int i = 0; i < 12; ++i) {
    flatRows[static_cast<std::size_t>(i)] = {i};
    for (int j = 0; j < 6; ++j) {
      flatRows[static_cast<std::size_t>(i)].push_back((i < 6 ? 6 : 0) + j);
    }
  }
  std::vector<std::vector<int>> flatRings = {{0, 1, 2, 3, 4, 5},
                                             {6, 7, 8, 9, 10, 11}};
  std::vector<ring::strucType> flatType(2, ring::strucType::unclassified);
  std::vector<double> flatRmsd(12, -1.0);
  (void)prism3::findBulkPrisms(flatRings, flatType, flatRows, flat, flatRmsd,
                               50.0);
}

TEST_CASE("mixed ring lengths and a negative atom stay inside the cage search",
          "[topo_bulk]") {
  std::mt19937 rng(20261109);
  std::uniform_int_distribution<int> atom(-1, 17);
  std::uniform_int_distribution<int> length(3, 8);
  for (int trial = 0; trial < 30; ++trial) {
    const int nRings = 6 + trial % 8;
    std::vector<std::vector<int>> rings;
    rings.push_back({0, 1, 2, 3, 4, 5});
    for (int r = 1; r < nRings; ++r) {
      const int n = length(rng);
      std::vector<int> ring(static_cast<std::size_t>(n));
      for (int k = 0; k < n; ++k) {
        ring[static_cast<std::size_t>(k)] = atom(rng);
      }
      rings.push_back(std::move(ring));
    }
    const auto index = ring::buildRingSearchIndex(rings, 18);
    std::vector<int> peripheral;
    std::vector<int> chosen;
    std::vector<int> rank(rings.size());
    std::iota(rank.begin(), rank.end(), 0);
    ring::EquatorChoice raw;
    (void)ring::equatorialConditions(rings, 0, index, &peripheral, &chosen,
                                     &rank, &raw);
    std::vector<ring::strucType> types(rings.size(), ring::strucType::unclassified);
    std::vector<int> listHC;
    std::vector<cage::Cage> cages;
    (void)ring::findDDC(rings, types, listHC, cages, index);
    std::vector<std::vector<int>> rows(18);
    for (int a = 0; a < 18; ++a) {
      rows[static_cast<std::size_t>(a)] = {a, (a + 1) % 18, (a + 3) % 18};
    }
    (void)ring::findHC(rings, types, rows, cages, index);
  }
  std::vector<std::vector<int>> neg = {{-1, -1, -1, -1, -1, -1},
                                       {0, 1, 2, 3, 4, 5}};
  std::vector<std::vector<int>> negRows(6);
  for (int a = 0; a < 6; ++a) {
    negRows[static_cast<std::size_t>(a)] = {a, (a + 1) % 6};
  }
  std::vector<ring::strucType> negType(2, ring::strucType::unclassified);
  std::vector<cage::Cage> negCages;
  (void)ring::findHC(neg, negType, negRows, negCages);

  molSys::PointCloud<molSys::Point<double>, double> ddcCloud;
  const double ddcXyz[][3] = {
      {25.4799996, 25.4799996, 19.1100006}, {25.4799996, 25.4799996, 25.4799996},
      {25.4799996, 22.2950001, 22.2950001}, {25.4799996, 28.6650009, 22.2950001},
      {22.2950001, 22.2950001, 25.4799996}, {22.2950001, 25.4799996, 22.2950001},
      {28.6650009, 25.4799996, 22.2950001}, {23.8880006, 20.7029991, 23.8880006},
      {23.8880006, 27.073, 23.8880006},     {27.073, 27.073, 20.7029991},
      {20.7029991, 23.8880006, 23.8880006}, {27.073, 23.8880006, 23.8880006},
      {23.8880006, 23.8880006, 20.7029991}, {23.8880006, 23.8880006, 27.073},
  };
  ddcCloud.box = {40, 40, 40};
  ddcCloud.boxLow = {0, 0, 0};
  for (int i = 0; i < 14; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = ddcXyz[i][0];
    pt.y = ddcXyz[i][1];
    pt.z = ddcXyz[i][2];
    ddcCloud.pts.push_back(pt);
    ddcCloud.idIndexMap[i] = i;
  }
  ddcCloud.nop = 14;
  auto ddcListed = nneigh::neighListO(3.5, ddcCloud, 1);
  auto ddcBy = nneigh::neighbourListByIndex(ddcCloud, ddcListed);
  auto ddcRings = primitive::ringNetwork(ddcBy, 8);
  if (!ddcRings.empty()) {
    ddcRings.push_back(ddcRings[0]);
    ddcRings.push_back(std::vector<int>{0, 1, 2, 3, 4});
    ddcRings.push_back(std::vector<int>{-1, 0, 1, 2, 3, 4});
    auto ddcIndex = ring::buildRingSearchIndex(ddcRings, 14);
    std::vector<int> peripheral;
    std::vector<int> chosen;
    std::vector<int> rank(ddcRings.size());
    std::iota(rank.begin(), rank.end(), 0);
    ring::EquatorChoice raw;
    for (int i = 0; i < static_cast<int>(ddcRings.size()); ++i) {
      if (ddcRings[static_cast<std::size_t>(i)].size() == 6) {
        (void)ring::equatorialConditions(ddcRings, i, ddcIndex, &peripheral,
                                         &chosen, &rank, &raw);
      }
    }
    std::vector<ring::strucType> ddcTypes(ddcRings.size(),
                                          ring::strucType::unclassified);
    std::vector<int> ddcHC;
    std::vector<cage::Cage> ddcCages;
    (void)ring::findDDC(ddcRings, ddcTypes, ddcHC, ddcCages, ddcIndex);
    (void)ring::findHC(ddcRings, ddcTypes, ddcBy, ddcCages, ddcIndex);
  }
}

TEST_CASE("bonded coplanar rings fail the shape match and use the strict basal test",
          "[topo_bulk]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {200, 200, 200};
  cloud.boxLow = {0, 0, 0};
  cloud.currentFrame = 1;
  // Three atoms of the second ring sit over the square and one sits
  // far away, so the scaled block stays beyond the prism RMSD. Every
  // atom of the first ring is still bonded to every atom of the second.
  const double coords[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {80, 0, 1}};
  for (int i = 0; i < 8; ++i) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 8;
  std::vector<std::vector<int>> nList(8);
  for (int i = 0; i < 8; ++i) {
    nList[static_cast<std::size_t>(i)] = {i};
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = 4; j < 8; ++j) {
      nList[static_cast<std::size_t>(i)].push_back(j);
      nList[static_cast<std::size_t>(j)].push_back(i);
    }
  }
  std::vector<int> basal1 = {0, 1, 2, 3};
  std::vector<int> basal2 = {4, 5, 6, 7};
  REQUIRE(prism3::basalPrismConditions(nList, basal1, basal2));
  const std::vector<std::vector<int>> rings = {basal1, basal2};
  std::vector<ring::strucType> ringType(2, ring::strucType::unclassified);
  std::vector<double> rmsd(8, -1.0);
  REQUIRE(prism3::findBulkPrisms(rings, ringType, nList, cloud, rmsd, 50.0) == 0);
  REQUIRE(ringType[0] == ring::strucType::Prism);
  REQUIRE(ringType[1] == ring::strucType::Prism);
}

// A kept hex has three equator atoms and three others. Pairing needs two of
// those others, so the seven rings of a chain cannot pass 14 molecules.
// The search still has to see a B-C partner that is not an A-C partner.
TEST_CASE("a side partner reaches the two-middle chain the triple misses",
          "[topo_bulk]") {
  std::vector<std::vector<int>> rings;
  rings.push_back({0, 1, 2, 3, 4, 5});
  auto add = [&](std::vector<int> ring) { rings.push_back(std::move(ring)); };
  for (int k : {0, 2, 4}) {
    for (int i = 0; i < 7; ++i) {
      const int base = 100 + k * 20 + i * 3;
      add({k % 6, (k + 1) % 6, (k + 2) % 6, base, base + 1, base + 2});
    }
  }
  for (int k : {1, 3, 5}) {
    const int base = 300 + k * 3;
    add({k % 6, (k + 1) % 6, (k + 2) % 6, base, base + 1, base + 2});
  }
  add({0, 1, 2, 10, 11, 12});
  add({2, 3, 4, 10, 11, 13});
  add({4, 5, 0, 10, 11, 14});
  add({4, 5, 0, 10, 13, 15});
  add({2, 3, 4, 10, 14, 19});
  add({1, 2, 3, 20, 21, 22});
  add({3, 4, 5, 20, 21, 23});
  add({5, 0, 1, 20, 21, 24});
  add({5, 0, 1, 20, 23, 25});
  const auto index = ring::buildRingSearchIndex(rings, 400);
  (void)ring::equatorialConditions(rings, 0, index, nullptr);
  std::vector<int> peripheral;
  (void)ring::equatorialConditions(rings, 0, index, &peripheral);
  std::vector<int> chosen;
  std::vector<int> rank(rings.size());
  std::iota(rank.begin(), rank.end(), 0);
  ring::EquatorChoice raw;
  (void)ring::equatorialConditions(rings, 0, index, &peripheral, &chosen, &rank,
                                   &raw);
  rank[0] = 0;
  (void)ring::equatorialConditions(rings, 0, index, nullptr, &chosen, &rank, &raw);
}

TEST_CASE("an equator id past the inline six reads the overflow", "[topo_bulk]") {
  ring::EquatorChoice choice;
  choice.put(6, 9);
  REQUIRE(choice.idAt(0) == 0);
  REQUIRE(choice.idAt(6) == 9);
}
