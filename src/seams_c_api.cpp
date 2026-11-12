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

#include <seams_c_api.h>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <diamond_cna.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <exception>
#include <string>
#include <vector>

int seams_chill_plus(const double *xyz, int n, const double *box, int *labels) {
  return seams_chill_plus_graph(xyz, n, box, 3.5, nullptr, labels);
}

namespace {

int chillLabels(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                double cutoff, const std::string &name, const char *rule,
                int *labels) {
  const auto nList =
      name == "cutoff"
          ? nneigh::neighListO(cutoff, cloud, 1)
          : nneigh::kNearestNeighbourList(cloud, 4, cutoff + 1.5, 1, true);
  chill::applyChillPlus(cloud, nList,
                        chill::chillPlusRule(rule == nullptr ? "" : rule),
                        false);
  for (int i = 0; i < cloud.nop; i++) {
    labels[i] = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
  }
  return 0;
}

molSys::PointCloud<molSys::Point<double>, double>
cloudFromXyz(const double *xyz, int n) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = xyz[3 * i];
    p.y = xyz[3 * i + 1];
    p.z = xyz[3 * i + 2];
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

} // namespace

int seams_chill_plus_graph(const double *xyz, int n, const double *box,
                           double cutoff, const char *graph, int *labels) {
  if (xyz == nullptr || box == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  if (!(box[3] > 0.0 && box[4] > 0.0 && box[5] > 0.0)) {
    return 1;
  }
  auto cloud = cloudFromXyz(xyz, n);
  cloud.box = {box[3], box[4], box[5]};
  cloud.boxLow = {box[0], box[1], box[2]};
  return chillLabels(cloud, cutoff, name, nullptr, labels);
}

int seams_chill_plus_rule(const double *xyz, int n, const double *box,
                          double cutoff, const char *graph, const char *rule,
                          int *labels) {
  if (xyz == nullptr || box == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  if (!(box[3] > 0.0 && box[4] > 0.0 && box[5] > 0.0)) {
    return 1;
  }
  try {
    auto cloud = cloudFromXyz(xyz, n);
    cloud.box = {box[3], box[4], box[5]};
    cloud.boxLow = {box[0], box[1], box[2]};
    return chillLabels(cloud, cutoff, name, rule, labels);
  } catch (const std::exception &) {
    return 1;
  }
}

int seams_chill_plus_graph_lammps(const double *xyz, int n,
                                  const double *box9, double cutoff,
                                  const char *graph, int *labels) {
  if (xyz == nullptr || box9 == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  auto cloud = cloudFromXyz(xyz, n);
  if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1],
                                    box9[2], box9[3], box9[4], box9[5],
                                    box9[6], box9[7], box9[8])) {
    return 1;
  }
  return chillLabels(cloud, cutoff, name, nullptr, labels);
}

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud cloudOf(const double *xyz, int n, const double *box9) {
  Cloud cloud;
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = xyz[3 * i];
    p.y = xyz[3 * i + 1];
    p.z = xyz[3 * i + 2];
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1],
                                    box9[2], box9[3], box9[4], box9[5], box9[6],
                                    box9[7], box9[8])) {
    throw std::invalid_argument("box edges are not positive");
  }
  return cloud;
}

bool badFrame(const double *xyz, int n, const double *box9, double cutoff) {
  return xyz == nullptr || box9 == nullptr || n <= 0 || !(cutoff > 0.0);
}

} // namespace

int seams_cutoff_csr(const double *xyz, int n, const double *box9, double cutoff,
                     int *offsets, int *cols, int cols_cap, int *n_cols) {
  if (badFrame(xyz, n, box9, cutoff) || offsets == nullptr || n_cols == nullptr) {
    return 1;
  }
  try {
    auto cloud = cloudOf(xyz, n, box9);
    std::vector<int> off;
    std::vector<int> col;
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, off, col, false)) {
      return 1;
    }
    const int count = static_cast<int>(col.size());
    *n_cols = count;
    if (static_cast<int>(off.size()) != n + 1) {
      return 1;
    }
    for (int i = 0; i < n + 1; i++) {
      offsets[i] = off[static_cast<std::size_t>(i)];
    }
    if (cols == nullptr) {
      return 0;
    }
    if (cols_cap < count) {
      return 2;
    }
    for (int i = 0; i < count; i++) {
      cols[i] = col[static_cast<std::size_t>(i)];
    }
    return 0;
  } catch (const std::exception &) {
    return 1;
  }
}

int seams_ql_cutoff(const double *xyz, int n, const double *box9, double cutoff,
                    int order_l, double *ql) {
  if (badFrame(xyz, n, box9, cutoff) || ql == nullptr) {
    return 1;
  }
  if (order_l != 3 && order_l != 4 && order_l != 6 && order_l != 8 &&
      order_l != 12) {
    return 1;
  }
  try {
    auto cloud = cloudOf(xyz, n, box9);
    const auto result = chill::steinhardtQlCutoff(cloud, cutoff, order_l, 1, false);
    if (static_cast<int>(result.ql.size()) != n) {
      return 1;
    }
    for (int i = 0; i < n; i++) {
      ql[i] = result.ql[static_cast<std::size_t>(i)];
    }
    return 0;
  } catch (const std::exception &) {
    return 1;
  }
}

int seams_diamond_cna(const double *xyz, int n, const double *box9, double cutoff,
                      int *labels) {
  if (badFrame(xyz, n, box9, cutoff) || labels == nullptr) {
    return 1;
  }
  try {
    auto cloud = cloudOf(xyz, n, box9);
    std::vector<int> off;
    std::vector<int> col;
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, off, col, false)) {
      return 1;
    }
    std::vector<int> kind;
    symm::diamondFromCsr(cloud, off, col, kind);
    for (int i = 0; i < n; i++) {
      labels[i] = kind[static_cast<std::size_t>(i)];
    }
    return 0;
  } catch (const std::exception &) {
    return 1;
  }
}

int seams_primitive_rings(const double *xyz, int n, const double *box9,
                          double cutoff, int max_depth, int *counts, int counts_len) {
  if (badFrame(xyz, n, box9, cutoff) || counts == nullptr || max_depth < 3 ||
      counts_len < max_depth + 1) {
    return 1;
  }
  try {
    auto cloud = cloudOf(xyz, n, box9);
    std::vector<int> off;
    std::vector<int> col;
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, off, col, false)) {
      return 1;
    }
    const auto rings = primitive::ringNetwork(off, col, max_depth);
    for (int i = 0; i < counts_len; i++) {
      counts[i] = 0;
    }
    for (const auto &ring : rings) {
      const int size = static_cast<int>(ring.size());
      if (size >= 0 && size < counts_len) {
        counts[size] += 1;
      }
    }
    return 0;
  } catch (const std::exception &) {
    return 1;
  }
}
