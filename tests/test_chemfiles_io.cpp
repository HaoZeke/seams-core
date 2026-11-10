// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <chemfiles.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

std::map<int, int> ringSizes(Cloud cloud, int type) {
  std::map<int, int> hist;
  if (cloud.nop <= 0) {
    return hist;
  }
  auto rows = nneigh::neighListO(3.5, cloud, type);
  auto index = nneigh::neighbourListByIndex(cloud, rows);
  for (const auto &ring : primitive::ringNetwork(index, 6)) {
    hist[static_cast<int>(ring.size())] += 1;
  }
  return hist;
}

std::map<int, int> iceCounts(Cloud cloud, int type) {
  std::map<int, int> hist;
  auto rows = nneigh::neighListO(3.5, cloud, type);
  chill::getCorrelPlus(cloud, rows, false);
  chill::getIceTypePlusNoPrint(cloud, rows, false);
  for (const auto &point : cloud.pts) {
    if (point.type == type) {
      hist[static_cast<int>(point.iceType)] += 1;
    }
  }
  return hist;
}

void writeAtoms(const fs::path &path, const std::string &format,
                const std::vector<std::array<double, 3>> &xyz,
                const std::array<double, 3> &lengths, const std::string &name) {
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(
      chemfiles::Vector3D(lengths[0], lengths[1], lengths[2])));
  for (const auto &p : xyz) {
    frame.add_atom(chemfiles::Atom(name),
                   chemfiles::Vector3D(p[0], p[1], p[2]));
  }
  if (format == std::string("XTC")) {
    frame.set("xtc_precision", 1000000.0);
  }
  chemfiles::Trajectory out(path.string(), 'w', format);
  out.write(frame);
}

std::vector<std::array<double, 3>> oxygenCoords(const std::string &dump,
                                                int type) {
  Cloud native;
  native = sinp::readLammpsTrjO(dump, 1, native, type);
  std::vector<std::array<double, 3>> xyz;
  xyz.reserve(static_cast<std::size_t>(native.nop));
  for (const auto &point : native.pts) {
    xyz.push_back({point.x, point.y, point.z});
  }
  return xyz;
}

} // namespace

TEST_CASE("chemfiles matches the native LAMMPS and XYZ readers",
          "[chemfiles]") {
  Cloud native;
  native = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, native, 2);
  Cloud via = sinp::openFrame("traj/exampleTraj.lammpstrj", 1, 2);
  REQUIRE(native.nop > 10);
  REQUIRE(via.nop == native.nop);
  REQUIRE(via.box.size() >= 3);
  for (int k = 0; k < 3; ++k) {
    REQUIRE_THAT(via.box[static_cast<std::size_t>(k)],
                 Catch::Matchers::WithinAbs(native.box[static_cast<std::size_t>(k)],
                                            1e-6));
  }
  for (int i = 0; i < native.nop; ++i) {
    REQUIRE(via.pts[static_cast<std::size_t>(i)].type ==
            native.pts[static_cast<std::size_t>(i)].type);
    REQUIRE(via.pts[static_cast<std::size_t>(i)].atomID ==
            native.pts[static_cast<std::size_t>(i)].atomID);
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].x,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].x,
                                            1e-4));
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].y,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].y,
                                            1e-4));
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].z,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].z,
                                            1e-4));
  }
  REQUIRE(ringSizes(via, 2) == ringSizes(native, 2));
  REQUIRE(iceCounts(via, 2) == iceCounts(native, 2));

  const auto path = fs::temp_directory_path() / "dseams_chem_native.xyz";
  {
    std::ofstream out(path);
    out << native.nop << "\nbox\n";
    for (const auto &point : native.pts) {
      out << "O " << point.x << " " << point.y << " " << point.z << "\n";
    }
  }
  auto xyzNative = sinp::readXYZ(path.string());
  Cloud xyzChem;
  xyzChem = sinp::readChemfiles(path.string(), 1, xyzChem, -1);
  REQUIRE(xyzChem.nop == xyzNative.nop);
  REQUIRE(xyzChem.nop == native.nop);
  for (int i = 0; i < xyzNative.nop; ++i) {
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].x,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].x, 1e-4));
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].y,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].y, 1e-4));
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].z,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].z, 1e-4));
  }
  fs::remove(path);
}

TEST_CASE("XTC DCD and NetCDF carry the same census as the dump",
          "[chemfiles]") {
  const auto xyz = oxygenCoords("traj/exampleTraj.lammpstrj", 2);
  REQUIRE(xyz.size() > 10);
  Cloud native;
  native = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, native, 2);
  const std::array<double, 3> lengths = {
      native.box[0], native.box[1], native.box[2]};
  const auto nativeRings = ringSizes(native, 2);
  const auto dir = fs::temp_directory_path();
  const struct {
    const char *file;
    const char *format;
    double tol;
    bool census;
  } formats[] = {
      {"dseams_chem.dcd", "DCD", 1e-4, true},
      {"dseams_chem.xtc", "XTC", 1e-3, true},
      {"dseams_chem.nc", "Amber NetCDF", 1e-4, true},
      {"dseams_chem.gro", "GRO", 1e-2, false},
      {"dseams_chem.pdb", "PDB", 1e-2, false},
  };
  for (const auto &fmt : formats) {
    const auto path = dir / fmt.file;
    writeAtoms(path, fmt.format, xyz, lengths, "O");
    Cloud got = sinp::openFrame(path.string(), 1, -1);
    INFO(fmt.file);
    REQUIRE(got.nop == native.nop);
    for (int i = 0; i < got.nop; ++i) {
      REQUIRE_THAT(got.pts[static_cast<std::size_t>(i)].x,
                   Catch::Matchers::WithinAbs(
                       native.pts[static_cast<std::size_t>(i)].x, fmt.tol));
    }
    if (fmt.census) {
      REQUIRE(ringSizes(got, got.pts[0].type) == nativeRings);
      REQUIRE(iceCounts(got, got.pts[0].type) == iceCounts(native, 2));
    }
    fs::remove(path);
  }
}

TEST_CASE("water names and a selection split oxygen from hydrogen",
          "[chemfiles]") {
  const auto path = fs::temp_directory_path() / "dseams_chem_water.xyz";
  {
    std::ofstream out(path);
    out << "3\n\nO 0 0 0\nH 1 0 0\nHW 0 1 0\n";
  }
  Cloud all;
  sinp::ChemfilesRead request;
  request.frame = 1;
  request.typeFilter = -1;
  all = sinp::readChemfiles(path.string(), all, request);
  REQUIRE(all.nop == 3);
  REQUIRE(all.pts[0].type == 8);
  REQUIRE(all.pts[1].type == 1);
  REQUIRE(all.pts[2].type == 1);

  Cloud guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 1);
  REQUIRE(guessed.pts[0].type == 8);
  Cloud selected = sinp::openFrame(path.string(), 1, 8, "name O");
  REQUIRE(selected.nop == 1);
  REQUIRE(selected.pts[0].type == 8);

  request.selection = "name H";
  request.guessOxygen = false;
  Cloud hydrogens;
  hydrogens = sinp::readChemfiles(path.string(), hydrogens, request);
  REQUIRE(hydrogens.nop == 1);
  REQUIRE(hydrogens.pts[0].type == 1);

  request.selection = "name O";
  request.typeFilter = 8;
  Cloud both;
  both = sinp::readChemfiles(path.string(), both, request);
  REQUIRE(both.nop == 1);

  request.selection = "name Zz";
  Cloud none;
  none = sinp::readChemfiles(path.string(), none, request);
  REQUIRE(none.nop == 0);

  request.selection = "this is not a selection";
  request.typeFilter = -1;
  Cloud bad;
  bad = sinp::readChemfiles(path.string(), bad, request);
  REQUIRE(bad.nop == 0);

  request.selection.clear();
  request.frame = 0;
  Cloud early;
  early = sinp::readChemfiles(path.string(), early, request);
  REQUIRE(early.nop == 0);

  request.frame = 4;
  Cloud late;
  late = sinp::readChemfiles(path.string(), late, request);
  REQUIRE(late.nop == 0);
  fs::remove(path);
}

TEST_CASE("a numeric dump guesses type 2 and then type 1", "[chemfiles]") {
  const auto path = fs::temp_directory_path() / "dseams_chem_types.xyz";
  {
    std::ofstream out(path);
    out << "3\n\n1 0 0 0\n1 1 0 0\n3 2 0 0\n";
  }
  Cloud guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 2);
  REQUIRE(guessed.pts[0].type == 1);

  {
    std::ofstream out(path);
    out << "2\n\n2 0 0 0\n1 3 0 0\n";
  }
  guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 1);
  REQUIRE(guessed.pts[0].type == 2);

  Cloud missing = sinp::openFrame("no/such/trajectory.lammpstrj", 1, -1);
  REQUIRE(missing.nop == 0);
  const auto junkDump = fs::temp_directory_path() / "dseams_chem_junk.lammpstrj";
  const auto junkXyz = fs::temp_directory_path() / "dseams_chem_junk.xyz";
  {
    std::ofstream out(junkDump);
    out << "this is not a dump\n";
    std::ofstream xyz(junkXyz);
    xyz << "not an xyz either\n";
  }
  REQUIRE(sinp::openFrame(junkDump.string(), 1, -1).nop == 0);
  REQUIRE(sinp::openFrame(junkDump.string(), 1, 2).nop == 0);
  REQUIRE(sinp::openFrame(junkDump.string(), 1, 0).nop == 0);
  REQUIRE(sinp::openFrame(junkXyz.string(), 1, -1).nop == 0);
  fs::remove(junkDump);
  fs::remove(junkXyz);
  REQUIRE(sinp::chemfilesSteps("no/such/trajectory.lammpstrj") == -1);
  REQUIRE_FALSE(sinp::forEachChemfilesFrame(
      "no/such/trajectory.lammpstrj", 1, 1, -1, "", 0, 1,
      [&](int, Cloud &) {}));
  REQUIRE(sinp::openFrame("noextension", 1, -1).nop == 0);
  REQUIRE(sinp::chemfilesSteps(path.string()) == 1);
  fs::remove(path);
}

TEST_CASE("a range walks every chemfiles step and skips CON", "[chemfiles]") {
  const auto path = fs::temp_directory_path() / "dseams_chem_range.xyz";
  {
    std::ofstream out(path);
    out << "1\n\nO 0 0 0\n1\n\nO 1 0 0\n1\n\nO 2 0 0\n";
  }
  std::vector<int> seen;
  const bool walked = sinp::forEachChemfilesFrame(
      path.string(), 1, 3, -1, "", 0, 1,
      [&](int frame, Cloud &cloud) {
        seen.push_back(frame);
        REQUIRE(cloud.nop == 1);
      });
  REQUIRE(walked);
  REQUIRE(seen.size() == 3);
  seen.clear();
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 0, -1, "", 0, 1,
      [&](int frame, Cloud &) { seen.push_back(frame); }));
  REQUIRE(seen.size() == 3);

  seen.clear();
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 3, -1, "", 1, 2,
      [&](int frame, Cloud &) { seen.push_back(frame); }));
  REQUIRE(seen.size() == 1);
  REQUIRE(seen[0] == 2);

  int calls = 0;
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 1, -1, "", 3, 2,
      [&](int, Cloud &) { calls += 1; }));
  REQUIRE(calls == 0);
  calls = 0;
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 1, -1, "", 0, 0,
      [&](int, Cloud &) { calls += 1; }));
  REQUIRE(calls == 1);

  REQUIRE_FALSE(sinp::forEachChemfilesFrame(
      "con/tiny_multi_cuh2.con", 1, 1, -1, "", 0, 1,
      [&](int, Cloud &) {}));
#ifdef SEAMS_HAS_READCON
  Cloud via = sinp::openFrame("con/tiny_multi_cuh2.con", 1, -1);
  Cloud direct;
  direct = sinp::readCon("con/tiny_multi_cuh2.con", 1, direct);
  REQUIRE(via.nop == direct.nop);
  REQUIRE(via.nop > 0);
#endif
  fs::remove(path);
}

TEST_CASE("a triclinic cell is stored with three tilt factors", "[chemfiles]") {
  const auto path = fs::temp_directory_path() / "dseams_chem_tilt.pdb";
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(10, 10, 10),
                                     chemfiles::Vector3D(90, 90, 60)));
  frame.add_atom(chemfiles::Atom("O"), chemfiles::Vector3D(1, 1, 1));
  chemfiles::Trajectory out(path.string(), 'w', "PDB");
  out.write(frame);
  out.close();
  Cloud cloud = sinp::openFrame(path.string(), 1, -1);
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.box.size() == 6);
  fs::remove(path);

  const auto pdb = fs::temp_directory_path() / "dseams_chem_res.pdb";
  {
    std::ofstream out(pdb);
    out << "ATOM      1  OW  SOL     7       1.000   2.000   3.000  1.00  0.00\n"
        << "END\n";
  }
  Cloud residue;
  residue = sinp::openFrame(pdb.string(), 1, -1);
  REQUIRE(residue.nop == 1);
  REQUIRE(residue.pts[0].molID == 7);
  fs::remove(pdb);
}

TEST_CASE("a dcd uses the sibling pdb names", "[chemfiles]") {
  const auto root = fs::temp_directory_path();
  const auto dcd = root / "dseams_chem_names.dcd";
  const auto pdb = root / "dseams_chem_names.pdb";
  {
    chemfiles::Frame frame;
    frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(20, 20, 20)));
    frame.add_atom(chemfiles::Atom("X"), chemfiles::Vector3D(1, 1, 1));
    frame.add_atom(chemfiles::Atom("X"), chemfiles::Vector3D(2, 1, 1));
    chemfiles::Trajectory out(dcd.string(), 'w', "DCD");
    out.write(frame);
  }
  {
    std::ofstream out(pdb);
    out << "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        << "ATOM      1  O   HOH A   1       1.000   1.000   1.000  1.00  0.00           O\n"
        << "ATOM      2  H1  HOH A   1       2.000   1.000   1.000  1.00  0.00           H\n"
        << "END\n";
  }
  Cloud oxygen = sinp::openFrame(dcd.string(), 1, 8);
  REQUIRE(oxygen.nop == 1);
  REQUIRE(oxygen.pts[0].type == 8);
  fs::remove(dcd);
  fs::remove(pdb);
}
