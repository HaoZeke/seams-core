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

// Per-frame neighbour list, primitive rings and published CHILL+ for one
// mW window, plus a readcon-db round trip of the .con trajectory.
// The cutoff and the ring depth are the ones the nucleation census uses.

#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <readcon-core.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using namespace readcon;

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr double kCutoffA = 3.5;
constexpr int kRingDepth = 7;
constexpr std::uint32_t kRngSeed = 20261118u;
constexpr int kMeltSteps = 10000;
constexpr double kStepFs = 10.0;
constexpr int kIceKinds = 9;
constexpr int kRingSizes = 8;
constexpr char kSymbol[] = "O";
constexpr double kMass = 18.015;

struct DumpFrame {
  long long step = 0;
  double lo[3] = {};
  double hi[3] = {};
  std::vector<int> id;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

struct Facts {
  int nop = 0;
  int rings[kRingSizes] = {};
  int ice[kIceKinds] = {};
  std::vector<int> pairs;
  std::vector<int> ring6;
  std::vector<int> iceId;
  std::vector<int> iceCode;
  double box[3] = {};
  double boxLow[3] = {};
};

std::string num(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

void require(bool ok, const std::string &msg) {
  if (!ok) {
    throw std::runtime_error(msg);
  }
}

std::vector<std::string> words(const std::string &line) {
  std::istringstream in(line);
  std::vector<std::string> out;
  std::string tok;
  while (in >> tok) {
    out.push_back(tok);
  }
  return out;
}

int column(const std::vector<std::string> &header, const std::string &name) {
  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    if (header[static_cast<std::size_t>(i)] == name) {
      return i;
    }
  }
  return -1;
}

std::vector<DumpFrame> readDump(const std::string &path) {
  std::ifstream in(path);
  require(static_cast<bool>(in), "cannot open dump " + path);
  std::vector<DumpFrame> frames;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) != 0) {
      continue;
    }
    DumpFrame frame;
    require(static_cast<bool>(std::getline(in, line)), "truncated timestep");
    frame.step = std::stoll(line);
    require(static_cast<bool>(std::getline(in, line)), "truncated atom header");
    require(line.rfind("ITEM: NUMBER OF ATOMS", 0) == 0, "missing atom count");
    require(static_cast<bool>(std::getline(in, line)), "truncated atom count");
    const int nop = std::stoi(line);
    require(nop > 0, "empty frame");
    require(static_cast<bool>(std::getline(in, line)), "truncated box header");
    require(line.rfind("ITEM: BOX BOUNDS", 0) == 0, "missing box");
    require(line.find("xy") == std::string::npos, "tilted box is not this deck");
    for (int axis = 0; axis < 3; ++axis) {
      require(static_cast<bool>(std::getline(in, line)), "truncated box");
      const auto tok = words(line);
      require(tok.size() >= 2, "box line");
      frame.lo[axis] = std::stod(tok[0]);
      frame.hi[axis] = std::stod(tok[1]);
    }
    require(static_cast<bool>(std::getline(in, line)), "truncated atom columns");
    require(line.rfind("ITEM: ATOMS", 0) == 0, "missing atoms");
    auto header = words(line);
    if (header.size() >= 2 && header[0] == "ITEM:" && header[1] == "ATOMS") {
      header.erase(header.begin(), header.begin() + 2);
    }
    const int idCol = column(header, "id");
    const int xCol = column(header, "x");
    const int yCol = column(header, "y");
    const int zCol = column(header, "z");
    require(idCol >= 0 && xCol >= 0 && yCol >= 0 && zCol >= 0, "dump columns");
    const int need = std::max(std::max(idCol, xCol), std::max(yCol, zCol));
    frame.id.reserve(static_cast<std::size_t>(nop));
    frame.x.reserve(static_cast<std::size_t>(nop));
    frame.y.reserve(static_cast<std::size_t>(nop));
    frame.z.reserve(static_cast<std::size_t>(nop));
    for (int i = 0; i < nop; ++i) {
      require(static_cast<bool>(std::getline(in, line)), "truncated atoms");
      const auto tok = words(line);
      require(static_cast<int>(tok.size()) > need, "short atom line");
      frame.id.push_back(std::stoi(tok[static_cast<std::size_t>(idCol)]));
      frame.x.push_back(std::stod(tok[static_cast<std::size_t>(xCol)]));
      frame.y.push_back(std::stod(tok[static_cast<std::size_t>(yCol)]));
      frame.z.push_back(std::stod(tok[static_cast<std::size_t>(zCol)]));
    }
    frames.push_back(std::move(frame));
  }
  require(!frames.empty(), "dump has no frames");
  return frames;
}

void writeCon(const std::string &path, const std::vector<DumpFrame> &frames) {
  std::filesystem::remove(path);
  // Decimal places, not the shortest scientific form: the con parser's
  // atom lines are read by a decimal kernel that rejects an exponent.
  RKRConFrameWriter *writer = create_writer_from_path_with_precision_c(path.c_str(), 17);
  require(writer != nullptr, "cannot create con writer");
  for (const DumpFrame &frame : frames) {
    const double cell[3] = {frame.hi[0] - frame.lo[0], frame.hi[1] - frame.lo[1],
                            frame.hi[2] - frame.lo[2]};
    const double angles[3] = {90.0, 90.0, 90.0};
    RKRConFrameBuilder *builder =
        rkr_frame_new(cell, angles, "mW", nullptr, nullptr, nullptr);
    require(builder != nullptr, "cannot build con frame");
    const auto nop = frame.id.size();
    for (std::size_t i = 0; i < nop; ++i) {
      const auto st = rkr_frame_add_atom(builder, kSymbol, frame.x[i], frame.y[i],
                                         frame.z[i], false,
                                         static_cast<std::uint64_t>(frame.id[i]), kMass);
      if (st != RKR_STATUS_SUCCESS) {
        free_rkr_frame_builder(builder);
        free_rkr_writer(writer);
        throw std::runtime_error("con atom rejected");
      }
    }
    RKRConFrame *built = rkr_frame_builder_build(builder);
    if (built == nullptr) {
      free_rkr_writer(writer);
      throw std::runtime_error("con frame rejected");
    }
    const RKRConFrame *one = built;
    const auto st = rkr_writer_extend(writer, &one, 1);
    free_rkr_frame(built);
    if (st != RKR_STATUS_SUCCESS) {
      free_rkr_writer(writer);
      throw std::runtime_error("con write failed");
    }
  }
  require(rkr_writer_flush(writer) == RKR_STATUS_SUCCESS, "con flush failed");
  free_rkr_writer(writer);
}

std::uint32_t frameSeed(int frame) {
  return kRngSeed + static_cast<std::uint32_t>(frame) * 0x9E3779B9u;
}

void writeXyzFrame(std::ostream &out, const Cloud &cloud) {
  const int n = cloud.nop;
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    order[static_cast<std::size_t>(i)] = i;
  }
  std::mt19937 rng(frameSeed(cloud.currentFrame));
  std::shuffle(order.begin(), order.end(), rng);
  out << n << "\n";
  out << "lattice";
  for (int k = 0; k < 3; ++k) {
    const double len = cloud.box.size() > static_cast<std::size_t>(k) ? cloud.box[k] : 0.0;
    const double lo =
        cloud.boxLow.size() > static_cast<std::size_t>(k) ? cloud.boxLow[k] : 0.0;
    out << " " << num(len);
    out << " " << num(lo);
  }
  out << "\n";
  for (const int i : order) {
    const auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    out << kSymbol << " " << num(pt.x) << " " << num(pt.y) << " " << num(pt.z) << "\n";
  }
}

void writeXyz(const std::string &conPath, const std::string &xyzPath, int nframes) {
  std::ofstream out(xyzPath);
  require(static_cast<bool>(out), "cannot write xyz");
  for (int frame = 1; frame <= nframes; ++frame) {
    Cloud cloud;
    cloud = sinp::readCon(conPath, frame, cloud);
    require(cloud.nop > 0, "con frame missing while writing xyz");
    writeXyzFrame(out, cloud);
  }
}

void applyLattice(Cloud &cloud, const std::string &comment) {
  const auto tok = words(comment);
  require(tok.size() == 7 && tok[0] == "lattice", "xyz comment is not a cell");
  cloud.box = {std::stod(tok[1]), std::stod(tok[3]), std::stod(tok[5])};
  cloud.boxLow = {std::stod(tok[2]), std::stod(tok[4]), std::stod(tok[6])};
}

Cloud coordinateOrder(const Cloud &cloud) {
  Cloud ordered = cloud;
  std::sort(ordered.pts.begin(), ordered.pts.end(),
            [](const molSys::Point<double> &a, const molSys::Point<double> &b) {
              if (a.x != b.x) {
                return a.x < b.x;
              }
              if (a.y != b.y) {
                return a.y < b.y;
              }
              return a.z < b.z;
            });
  ordered.idIndexMap.clear();
  for (int i = 0; i < ordered.nop; ++i) {
    ordered.idIndexMap.emplace(ordered.pts[static_cast<std::size_t>(i)].atomID, i);
  }
  return ordered;
}

Facts analyze(Cloud &cloud) {
  Facts facts;
  facts.nop = cloud.nop;
  for (int k = 0; k < 3; ++k) {
    facts.box[k] = cloud.box.size() > static_cast<std::size_t>(k) ? cloud.box[k] : 0.0;
    facts.boxLow[k] =
        cloud.boxLow.size() > static_cast<std::size_t>(k) ? cloud.boxLow[k] : 0.0;
  }
  if (cloud.nop <= 0) {
    return facts;
  }
  // The ring walk numbers vertices by index. Ordering by coordinates makes
  // the census a property of the geometry, then each vertex keeps its atom id.
  Cloud ordered = coordinateOrder(cloud);
  std::vector<int> offsets;
  std::vector<int> cols;
  require(nneigh::orthoFullCsr(ordered, nullptr, kCutoffA, offsets, cols, false),
          "cutoff neighbour list failed");
  const auto rings = primitive::ringNetwork(offsets, cols, kRingDepth);
  for (const auto &ring : rings) {
    const int sz = static_cast<int>(ring.size());
    if (sz >= 0 && sz < kRingSizes) {
      facts.rings[sz] += 1;
    }
    if (sz == 6) {
      for (const int a : ring) {
        if (a >= 0 && a < ordered.nop) {
          facts.ring6.push_back(ordered.pts[static_cast<std::size_t>(a)].atomID);
        }
      }
    }
  }
  std::sort(facts.ring6.begin(), facts.ring6.end());
  facts.ring6.erase(std::unique(facts.ring6.begin(), facts.ring6.end()), facts.ring6.end());

  offsets.clear();
  cols.clear();
  require(nneigh::orthoFullCsr(cloud, nullptr, kCutoffA, offsets, cols, false),
          "cutoff neighbour list failed");
  std::vector<std::pair<int, int>> pairs;
  pairs.reserve(cols.size() / 2);
  for (int i = 0; i < cloud.nop; ++i) {
    const int idI = cloud.pts[static_cast<std::size_t>(i)].atomID;
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      const int idJ = cloud.pts[static_cast<std::size_t>(j)].atomID;
      if (idI < idJ) {
        pairs.emplace_back(idI, idJ);
      }
    }
  }
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  facts.pairs.reserve(pairs.size() * 2);
  for (const auto &p : pairs) {
    facts.pairs.push_back(p.first);
    facts.pairs.push_back(p.second);
  }

  std::vector<int> indexToID(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    indexToID[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].atomID;
  }
  const auto idRows = nneigh::idRowsFromCsr(cloud, indexToID, offsets, cols);
  chill::applyChillPlus(cloud, idRows, chill::ChillPlusRule::Published, false);
  facts.iceId.reserve(static_cast<std::size_t>(cloud.nop));
  facts.iceCode.reserve(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    const int code = static_cast<int>(pt.iceType);
    if (code >= 0 && code < kIceKinds) {
      facts.ice[code] += 1;
    }
    facts.iceId.push_back(pt.atomID);
    facts.iceCode.push_back(code);
  }
  return facts;
}

void writeInts(std::ostream &out, const std::vector<int> &vals) {
  out << '[';
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if (i != 0) {
      out << ',';
    }
    out << vals[i];
  }
  out << ']';
}

void writeInts(std::ostream &out, const int *vals, int n) {
  out << '[';
  for (int i = 0; i < n; ++i) {
    if (i != 0) {
      out << ',';
    }
    out << vals[i];
  }
  out << ']';
}

void writeFacts(std::ostream &out, int frame, const Facts &facts) {
  out << "{\"frame\":" << frame << ",\"nop\":" << facts.nop << ",\"rings\":";
  writeInts(out, facts.rings, kRingSizes);
  out << ",\"ice\":";
  writeInts(out, facts.ice, kIceKinds);
  out << ",\"pairs\":";
  writeInts(out, facts.pairs);
  out << ",\"ring6\":";
  writeInts(out, facts.ring6);
  out << ",\"ice_id\":";
  writeInts(out, facts.iceId);
  out << ",\"ice_code\":";
  writeInts(out, facts.iceCode);
  out << ",\"box\":[" << num(facts.box[0]) << ',' << num(facts.box[1]) << ','
      << num(facts.box[2]) << ',' << num(facts.boxLow[0]) << ',' << num(facts.boxLow[1])
      << ',' << num(facts.boxLow[2]) << "]}\n";
}

void scanCon(const std::string &path, const std::string &outPath) {
  std::ofstream out(outPath);
  require(static_cast<bool>(out), "cannot write " + outPath);
  for (int frame = 1;; ++frame) {
    Cloud cloud;
    cloud = sinp::readCon(path, frame, cloud);
    if (cloud.nop <= 0) {
      require(frame > 1, "con file has no frames");
      break;
    }
    writeFacts(out, frame, analyze(cloud));
  }
}

void scanXyz(const std::string &path, const std::string &outPath) {
  std::ifstream in(path);
  require(static_cast<bool>(in), "cannot open xyz");
  std::ofstream out(outPath);
  require(static_cast<bool>(out), "cannot write " + outPath);
  const auto tmp = std::filesystem::path(path).parent_path() / "one-frame.xyz";
  int frame = 0;
  std::string countLine;
  while (std::getline(in, countLine)) {
    if (countLine.empty()) {
      continue;
    }
    const int nop = std::stoi(countLine);
    std::string comment;
    require(static_cast<bool>(std::getline(in, comment)), "xyz comment missing");
    std::ofstream one(tmp);
    require(static_cast<bool>(one), "cannot write a single xyz frame");
    one << nop << "\n" << comment << "\n";
    for (int i = 0; i < nop; ++i) {
      std::string atom;
      require(static_cast<bool>(std::getline(in, atom)), "xyz atom missing");
      one << atom << "\n";
    }
    one.close();
    Cloud cloud = sinp::readXYZ(tmp.string());
    require(cloud.nop == nop, "xyz reader kept a different atom count");
    applyLattice(cloud, comment);
    ++frame;
    writeFacts(out, frame, analyze(cloud));
  }
  std::filesystem::remove(tmp);
  require(frame > 0, "xyz file has no frames");
}

double maxAbsPos(const Cloud &a, const Cloud &b) {
  std::unordered_map<int, std::size_t> at;
  at.reserve(static_cast<std::size_t>(b.nop));
  for (int i = 0; i < b.nop; ++i) {
    at.emplace(b.pts[static_cast<std::size_t>(i)].atomID, static_cast<std::size_t>(i));
  }
  double worst = 0.0;
  for (int i = 0; i < a.nop; ++i) {
    const auto &pa = a.pts[static_cast<std::size_t>(i)];
    const auto found = at.find(pa.atomID);
    if (found == at.end()) {
      return -1.0;
    }
    const auto &pb = b.pts[found->second];
    worst = std::max(worst, std::abs(pa.x - pb.x));
    worst = std::max(worst, std::abs(pa.y - pb.y));
    worst = std::max(worst, std::abs(pa.z - pb.z));
  }
  return worst;
}

double maxAbsCell(const Cloud &a, const Cloud &b) {
  double worst = 0.0;
  const auto nBox = std::min(a.box.size(), b.box.size());
  if (a.box.size() != b.box.size() || a.boxLow.size() != b.boxLow.size()) {
    return -1.0;
  }
  for (std::size_t k = 0; k < nBox; ++k) {
    worst = std::max(worst, std::abs(a.box[k] - b.box[k]));
  }
  for (std::size_t k = 0; k < a.boxLow.size(); ++k) {
    worst = std::max(worst, std::abs(a.boxLow[k] - b.boxLow[k]));
  }
  return worst;
}

void roundTrip(const std::string &conPath, const std::string &corpus, int repeats) {
  require(repeats >= 8, "round trip wants at least 8 repeats");
  Cloud probe;
  probe = sinp::readCon(conPath, 1, probe);
  require(probe.nop > 0, "con file has no frames");
  const std::string formula = std::string(kSymbol) + ":" + std::to_string(probe.nop);
  const int cores = static_cast<int>(std::thread::hardware_concurrency());
  std::vector<double> seconds;
  seconds.reserve(static_cast<std::size_t>(repeats));
  double worstPos = 0.0;
  double worstCell = 0.0;
  int frames = 0;
  int exact = 1;
  for (int rep = 0; rep < repeats; ++rep) {
    const auto dir = std::filesystem::path(corpus) / ("rep-" + std::to_string(rep));
    std::filesystem::remove_all(dir);
    const auto t0 = std::chrono::steady_clock::now();
    const int n = sinp::ingestConCorpus(dir.string(), 1, conPath);
    if (n < 1) {
      exact = 0;
      frames = n;
      break;
    }
    frames = n;
    for (int frame = 1; frame <= n; ++frame) {
      Cloud fromFile;
      Cloud fromDb;
      fromFile = sinp::readCon(conPath, frame, fromFile);
      fromDb = sinp::readConCorpus(dir.string(), kSymbol, formula, frame, fromDb);
      if (fromFile.nop != fromDb.nop || fromFile.nop <= 0) {
        exact = 0;
        worstPos = -1.0;
        break;
      }
      bool ids = true;
      std::unordered_map<int, int> seen;
      for (int i = 0; i < fromFile.nop; ++i) {
        const int id = fromFile.pts[static_cast<std::size_t>(i)].atomID;
        if (!seen.emplace(id, 1).second) {
          ids = false;
        }
      }
      for (int i = 0; i < fromDb.nop; ++i) {
        const int id = fromDb.pts[static_cast<std::size_t>(i)].atomID;
        if (seen.find(id) == seen.end()) {
          ids = false;
        }
      }
      const double pos = maxAbsPos(fromFile, fromDb);
      const double cell = maxAbsCell(fromFile, fromDb);
      if (!ids || pos != 0.0 || cell != 0.0 || pos < 0.0 || cell < 0.0) {
        exact = 0;
      }
      if (pos > worstPos) {
        worstPos = pos;
      }
      if (cell > worstCell) {
        worstCell = cell;
      }
    }
    const auto t1 = std::chrono::steady_clock::now();
    seconds.push_back(std::chrono::duration<double>(t1 - t0).count());
    std::filesystem::remove_all(dir);
    if (exact == 0) {
      break;
    }
  }
  std::cout << "{\"frames\":" << frames << ",\"exact\":" << exact
            << ",\"n_cores\":" << cores << ",\"repeats\":" << static_cast<int>(seconds.size())
            << ",\"max_abs_position\":" << num(worstPos)
            << ",\"max_abs_cell\":" << num(worstCell) << ",\"seconds\":[";
  for (std::size_t i = 0; i < seconds.size(); ++i) {
    if (i != 0) {
      std::cout << ',';
    }
    std::cout << num(seconds[i]);
  }
  std::cout << "]}\n";
}

void emit(const std::string &dump, const std::string &conPath, const std::string &xyzPath) {
  const auto frames = readDump(dump);
  writeCon(conPath, frames);
  writeXyz(conPath, xyzPath, static_cast<int>(frames.size()));
  std::cout << "frames " << frames.size() << "\n";
  std::cout << "first_step " << frames.front().step << "\n";
  std::cout << "last_step " << frames.back().step << "\n";
  std::cout << "rng_seed " << kRngSeed << "\n";
  std::cout << "cutoff " << num(kCutoffA) << "\n";
  std::cout << "ring_depth " << kRingDepth << "\n";
  std::cout << "melt_steps " << kMeltSteps << "\n";
  std::cout << "step_fs " << num(kStepFs) << "\n";
  std::cout << "nop " << frames.front().id.size() << "\n";
}

}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2) {
      std::cerr << "usage: scan emit|scan-con|scan-xyz|roundtrip ...\n";
      return 2;
    }
    const std::string cmd = argv[1];
    if (cmd == "emit" && argc == 5) {
      emit(argv[2], argv[3], argv[4]);
      return 0;
    }
    if (cmd == "scan-con" && argc == 4) {
      scanCon(argv[2], argv[3]);
      return 0;
    }
    if (cmd == "scan-xyz" && argc == 4) {
      scanXyz(argv[2], argv[3]);
      return 0;
    }
    if (cmd == "roundtrip" && argc == 5) {
      roundTrip(argv[2], argv[3], std::stoi(argv[4]));
      return 0;
    }
    std::cerr << "usage: scan emit|scan-con|scan-xyz|roundtrip ...\n";
    return 2;
  } catch (const std::exception &err) {
    std::cerr << err.what() << "\n";
    return 1;
  }
}
