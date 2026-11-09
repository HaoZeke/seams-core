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

#include <seams_input.hpp>

#include <neighbours.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_CHEMFILES
#include <chemfiles.hpp>
#endif

namespace {

std::string extensionOf(const std::string &path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return {};
  }
  return path.substr(dot + 1);
}

#ifdef SEAMS_HAS_CHEMFILES

bool digitsOnly(const std::string &text) {
  return !text.empty() &&
         text.find_first_not_of("0123456789") == std::string::npos;
}

std::string upperToken(std::string text) {
  for (char &ch : text) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return text;
}

bool oxygenToken(const std::string &token) {
  const std::string name = upperToken(token);
  return name == "O" || name == "OW" || name == "OH";
}

bool hydrogenToken(const std::string &token) {
  const std::string name = upperToken(token);
  return name == "H" || name == "HW" || name == "HO";
}

// Numeric LAMMPS types stay numbers. Element names become atomic numbers,
// and the water names OW/HW do too when the table has no entry for them.
int speciesOf(const chemfiles::Atom &atom) {
  if (digitsOnly(atom.type())) {
    return std::stoi(atom.type());
  }
  if (const auto number = atom.atomic_number()) {
    return static_cast<int>(*number);
  }
  if (oxygenToken(atom.type()) || oxygenToken(atom.name())) {
    return 8;
  }
  if (hydrogenToken(atom.type()) || hydrogenToken(atom.name())) {
    return 1;
  }
  return 0;
}

int atomicNumberOf(const chemfiles::Atom &atom) {
  if (const auto number = atom.atomic_number()) {
    return static_cast<int>(*number);
  }
  return 0;
}

void rememberId(molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  const int idx = static_cast<int>(cloud.pts.size()) - 1;
  cloud.idIndexMap[cloud.pts.back().atomID] = idx;
}

void applyCell(molSys::PointCloud<molSys::Point<double>, double> &cloud,
               const chemfiles::UnitCell &cell) {
  if (cell.shape() == chemfiles::UnitCell::INFINITE) {
    cloud.box.clear();
    cloud.boxLow.clear();
    return;
  }
  const auto lengths = cell.lengths();
  const auto angles = cell.angles();
  nneigh::setBoxFromLengthsAngles(cloud.box, cloud.boxLow, lengths[0],
                                 lengths[1], lengths[2], angles[0], angles[1],
                                 angles[2]);
}

#endif

} // namespace

#ifdef SEAMS_HAS_CHEMFILES

// Set when chemfiles rejects the file, so openFrame can try the native reader.
// A selection that keeps no atoms does not set it.
thread_local bool chemfilesUnreadable = false;

int sinp::chemfilesSteps(const std::string &filename) {
  try {
    chemfiles::Trajectory trajectory(filename);
    return static_cast<int>(trajectory.nsteps());
  } catch (const chemfiles::Error &) {
    return -1;
  }
}

molSys::PointCloud<molSys::Point<double>, double>
sinp::readChemfiles(std::string filename, int targetFrame,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    int typeFilter) {
  ChemfilesRead request;
  request.frame = targetFrame;
  request.typeFilter = typeFilter;
  return readChemfiles(std::move(filename), yCloud, request);
}

molSys::PointCloud<molSys::Point<double>, double>
sinp::readChemfiles(std::string filename,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const ChemfilesRead &request) {
  chemfilesUnreadable = false;
  yCloud = molSys::clearPointCloud(yCloud);
  if (request.frame < 1) {
    return yCloud;
  }
  try {
    chemfiles::Trajectory trajectory(filename);
    // A DCD has coordinates and a cell, not atom names. A sibling PDB
    // carries the names the host used to choose the analysed species.
    if (extensionOf(filename) == "dcd") {
      const auto dot = filename.find_last_of('.');
      if (dot != std::string::npos) {
        const auto pdb = filename.substr(0, dot) + ".pdb";
        std::ifstream probe(pdb);
        if (probe) {
          trajectory.set_topology(pdb);
        }
      }
    }
    if (static_cast<std::size_t>(request.frame) > trajectory.nsteps()) {
      std::cerr << "Frame " << request.frame << " does not exist in "
                << filename << " (has " << trajectory.nsteps() << " frames).\n";
      return yCloud;
    }
    auto frame =
        trajectory.read_step(static_cast<std::size_t>(request.frame - 1));
    applyCell(yCloud, frame.cell());
    const auto &topology = frame.topology();
    const auto positions = frame.positions();
    const std::size_t n = frame.size();

    std::vector<char> wanted(n, 1);
    if (!request.selection.empty()) {
      std::fill(wanted.begin(), wanted.end(), 0);
      const chemfiles::Selection selection(request.selection);
      for (const std::size_t index : selection.list(frame)) {
        if (index < n) {
          wanted[index] = 1;
        }
      }
    }

    bool anyOxygen = false;
    bool anyElement = false;
    int numericTwos = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const int number = atomicNumberOf(topology[i]);
      const int species = speciesOf(topology[i]);
      const bool named =
          number > 0 || oxygenToken(topology[i].name()) ||
          oxygenToken(topology[i].type()) || hydrogenToken(topology[i].name()) ||
          hydrogenToken(topology[i].type());
      if (named) {
        anyElement = true;
      }
      if (number == 8 || species == 8) {
        anyOxygen = true;
      }
      if (!named && species == 2) {
        numericTwos += 1;
      }
    }
    const bool keepOxygen = request.guessOxygen && anyOxygen;
    const bool keepTwo =
        request.guessOxygen && !anyElement && numericTwos > 0;
    const bool keepOne =
        request.guessOxygen && !anyElement && numericTwos == 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (!wanted[i]) {
        continue;
      }
      int species = speciesOf(topology[i]);
      const int number = atomicNumberOf(topology[i]);
      // A file with no names and no elements is one species.
      // A GROMACS TRR is that case, and the live module labels it type 1.
      if (!anyElement && species == 0) {
        species = 1;
      }
      if (keepOxygen && number != 8 && species != 8) {
        continue;
      }
      if (keepTwo && species != 2) {
        continue;
      }
      if (keepOne && species != 1) {
        continue;
      }
      if (request.typeFilter > 0 && species != request.typeFilter) {
        continue;
      }
      molSys::Point<double> point;
      point.type = species;
      // Chemfiles places a LAMMPS atom at the index of its id, so the
      // row index plus one is that id. The mol column is a custom
      // property on the atom. A residue id covers formats that have one.
      point.atomID = static_cast<int>(i) + 1;
      point.molID = point.atomID;
      if (const auto mol = topology[i].get<chemfiles::Property::DOUBLE>("mol")) {
        point.molID = static_cast<int>(std::lround(*mol));
      } else if (const auto residue = topology.residue_for_atom(i)) {
        if (const auto residueId = residue->id()) {
          point.molID = static_cast<int>(*residueId);
        }
      }
      point.x = positions[i][0];
      point.y = positions[i][1];
      point.z = positions[i][2];
      point.inSlice = true;
      yCloud.pts.push_back(point);
      rememberId(yCloud);
    }
    yCloud.nop = static_cast<int>(yCloud.pts.size());
    yCloud.currentFrame = request.frame;
    return yCloud;
  } catch (const chemfiles::SelectionError &error) {
    std::cerr << "chemfiles selection on " << filename << ": " << error.what()
              << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  } catch (const chemfiles::Error &) {
    chemfilesUnreadable = true;
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
}

bool sinp::forEachChemfilesFrame(
    const std::string &filename, int first, int last, int typeI,
    const std::string &selection, int part, int parts,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn) {
  if (extensionOf(filename) == "con") {
    return false;
  }
  const int steps = chemfilesSteps(filename);
  if (steps < 0) {
    return false;
  }
  const int begin = std::max(1, first);
  int end = last;
  if (end <= 0 || end > steps) {
    end = steps;
  }
  if (parts < 1) {
    parts = 1;
  }
  if (part < 0 || part >= parts) {
    return true;
  }
  for (int frame = begin; frame <= end; ++frame) {
    if (((frame - begin) % parts) != part) {
      continue;
    }
    auto cloud = openFrame(filename, frame, typeI, selection);
    fn(frame, cloud);
  }
  return true;
}

#endif

molSys::PointCloud<molSys::Point<double>, double>
sinp::openFrame(const std::string &filename, int frame, int typeI,
                const std::string &selection) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  if (extensionOf(filename) == "con") {
#ifdef SEAMS_HAS_READCON
    return readCon(filename, frame, cloud);
#else
    throw std::runtime_error(
        "readcon-core is not in this build; CON input needs it");
#endif
  }
#ifdef SEAMS_HAS_CHEMFILES
  ChemfilesRead request;
  request.frame = frame;
  request.selection = selection;
  if (!selection.empty()) {
    request.typeFilter = typeI > 0 ? typeI : -1;
  } else if (typeI < 0) {
    request.typeFilter = -1;
  } else if (typeI > 0) {
    request.typeFilter = typeI;
  } else {
    request.guessOxygen = true;
  }
  cloud = readChemfiles(filename, cloud, request);
  if (!chemfilesUnreadable) {
    return cloud;
  }
  {
    const auto ext = extensionOf(filename);
    if (ext == "xyz") {
      return readXYZ(filename);
    }
    if (ext == "lammpstrj" || ext == "dump") {
      if (typeI < 0) {
        return readLammpsTrj(filename, frame, cloud);
      }
      if (typeI > 0) {
        return readLammpsTrjO(filename, frame, cloud, typeI);
      }
      cloud = readLammpsTrjO(filename, frame, cloud, 2);
      if (cloud.nop == 0) {
        molSys::PointCloud<molSys::Point<double>, double> again;
        cloud = readLammpsTrjO(filename, frame, again, 1);
      }
      return cloud;
    }
    std::cerr << "chemfiles cannot read " << filename << "\n";
    return molSys::clearPointCloud(cloud);
  }
#else
  (void)selection;
  const auto ext = extensionOf(filename);
  if (ext == "xyz") {
    return readXYZ(filename);
  }
  if (typeI < 0) {
    return readLammpsTrj(filename, frame, cloud);
  }
  if (typeI > 0) {
    return readLammpsTrjO(filename, frame, cloud, typeI);
  }
  cloud = readLammpsTrjO(filename, frame, cloud, 2);
  if (cloud.nop > 0) {
    return cloud;
  }
  molSys::PointCloud<molSys::Point<double>, double> again;
  return readLammpsTrjO(filename, frame, again, 1);
#endif
}

#ifndef SEAMS_HAS_CHEMFILES
int sinp::chemfilesSteps(const std::string &) { return -1; }
#endif
