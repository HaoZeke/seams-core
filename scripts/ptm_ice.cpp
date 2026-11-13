// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Diamond-only PTM on an orthorhombic cell. The neighbour callback is the
// one the reference library asks for: the nearest points, self first, with
// the minimum image. Only cubic and hexagonal diamond are enabled, and a
// match is dropped when its RMSD exceeds the cutoff, which is how the
// OVITO modifier treats rmsd_cutoff.

#include "ptm_functions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <omp.h>

namespace {

struct Cloud {
  const double *xyz;
  int n;
  double box[3];
  int ncell[3];
  double cell[3];
  std::vector<int> head;
  std::vector<int> next;
};

double mic(double delta, double length) {
  return delta - length * std::nearbyint(delta / length);
}

int slot(const Cloud &cloud, int ix, int iy, int iz) {
  const int nx = cloud.ncell[0];
  const int ny = cloud.ncell[1];
  const int nz = cloud.ncell[2];
  ix = (ix % nx + nx) % nx;
  iy = (iy % ny + ny) % ny;
  iz = (iz % nz + nz) % nz;
  return (ix * ny + iy) * nz + iz;
}

void build_cells(Cloud *cloud) {
  const double reach = 4.0;
  for (int axis = 0; axis < 3; ++axis) {
    int count = static_cast<int>(cloud->box[axis] / reach);
    if (count < 1) {
      count = 1;
    }
    cloud->ncell[axis] = count;
    cloud->cell[axis] = cloud->box[axis] / static_cast<double>(count);
  }
  const int cells = cloud->ncell[0] * cloud->ncell[1] * cloud->ncell[2];
  cloud->head.assign(static_cast<std::size_t>(cells), -1);
  cloud->next.assign(static_cast<std::size_t>(cloud->n), -1);
  for (int i = 0; i < cloud->n; ++i) {
    const double *p = cloud->xyz + 3 * i;
    int ix = static_cast<int>(std::floor(p[0] / cloud->cell[0]));
    int iy = static_cast<int>(std::floor(p[1] / cloud->cell[1]));
    int iz = static_cast<int>(std::floor(p[2] / cloud->cell[2]));
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix >= cloud->ncell[0]) ix = cloud->ncell[0] - 1;
    if (iy >= cloud->ncell[1]) iy = cloud->ncell[1] - 1;
    if (iz >= cloud->ncell[2]) iz = cloud->ncell[2] - 1;
    const int bin = slot(*cloud, ix, iy, iz);
    cloud->next[static_cast<std::size_t>(i)] = cloud->head[static_cast<std::size_t>(bin)];
    cloud->head[static_cast<std::size_t>(bin)] = i;
  }
}

struct Near {
  double d2;
  double x, y, z;
  int index;
};

bool nearer(const Near &a, const Near &b) {
  if (a.d2 < b.d2) return true;
  if (b.d2 < a.d2) return false;
  if (a.z < b.z) return true;
  if (b.z < a.z) return false;
  if (a.y < b.y) return true;
  if (b.y < a.y) return false;
  return a.x < b.x;
}

int get_neighbours(void *vdata, size_t, size_t atom_index, int num, int *ordering,
                   size_t *nbr_indices, int32_t *numbers, double (*nbr_pos)[3]) {
  const Cloud *cloud = static_cast<const Cloud *>(vdata);
  if (num < 1 || atom_index >= static_cast<size_t>(cloud->n)) {
    return 0;
  }
  const double *origin = cloud->xyz + 3 * atom_index;
  const int ix = static_cast<int>(std::floor(origin[0] / cloud->cell[0]));
  const int iy = static_cast<int>(std::floor(origin[1] / cloud->cell[1]));
  const int iz = static_cast<int>(std::floor(origin[2] / cloud->cell[2]));
  // The 18th ice neighbour sits near 6 A. Two cells of the 4 A grid cover it,
  // and a wrapped cell is visited once.
  std::vector<Near> found;
  found.reserve(96);
  int seen[125];
  int nseen = 0;
  for (int dx = -2; dx <= 2; ++dx) {
    for (int dy = -2; dy <= 2; ++dy) {
      for (int dz = -2; dz <= 2; ++dz) {
        const int bin = slot(*cloud, ix + dx, iy + dy, iz + dz);
        bool again = false;
        for (int s = 0; s < nseen; ++s) {
          if (seen[s] == bin) {
            again = true;
            break;
          }
        }
        if (again) continue;
        seen[nseen++] = bin;
        for (int j = cloud->head[static_cast<std::size_t>(bin)]; j >= 0;
             j = cloud->next[static_cast<std::size_t>(j)]) {
          const double *q = cloud->xyz + 3 * j;
          const double ox = mic(q[0] - origin[0], cloud->box[0]);
          const double oy = mic(q[1] - origin[1], cloud->box[1]);
          const double oz = mic(q[2] - origin[2], cloud->box[2]);
          found.push_back(Near{ox * ox + oy * oy + oz * oz, ox, oy, oz, j});
        }
      }
    }
  }
  const int keep = std::min(num, static_cast<int>(found.size()));
  if (keep > 0 && keep < static_cast<int>(found.size())) {
    std::nth_element(found.begin(), found.begin() + (keep - 1), found.end(), nearer);
  }
  std::sort(found.begin(), found.begin() + keep, nearer);
  for (int k = 0; k < keep; ++k) {
    const Near &nbr = found[static_cast<std::size_t>(k)];
    nbr_indices[k] = static_cast<size_t>(nbr.index);
    numbers[k] = 1;
    ordering[k] = k;
    nbr_pos[k][0] = nbr.x;
    nbr_pos[k][1] = nbr.y;
    nbr_pos[k][2] = nbr.z;
  }
  return keep;
}

}  // namespace

extern "C" int ptm_ice_classify(const double *xyz, int n, const double box[3],
                                double rmsd_cut, int threads, int *types, double *rmsd_out) {
  if (xyz == nullptr || box == nullptr || types == nullptr || n < 1) {
    return 1;
  }
  if (ptm_initialize_global() != PTM_NO_ERROR) {
    return 2;
  }
  if (threads < 1) {
    threads = 1;
  }
  Cloud cloud;
  cloud.xyz = xyz;
  cloud.n = n;
  cloud.box[0] = box[0];
  cloud.box[1] = box[1];
  cloud.box[2] = box[2];
  build_cells(&cloud);
  const int32_t flags = PTM_CHECK_DCUB | PTM_CHECK_DHEX;
  omp_set_num_threads(threads);
#pragma omp parallel
  {
    ptm_local_handle_t handle = ptm_initialize_local();
#pragma omp for schedule(static)
    for (int i = 0; i < n; ++i) {
      int32_t kind = PTM_MATCH_NONE;
      int32_t alloy = PTM_ALLOY_NONE;
      double scale = 0;
      double rmsd = INFINITY;
      double quat[4] = {};
      ptm_index(handle, static_cast<size_t>(i), get_neighbours, &cloud, flags, false, &kind,
                &alloy, &scale, &rmsd, quat, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                nullptr, nullptr, nullptr);
      if (!(std::isfinite(rmsd) && rmsd <= rmsd_cut)) {
        kind = PTM_MATCH_NONE;
      }
      types[i] = static_cast<int>(kind);
      if (rmsd_out != nullptr) {
        rmsd_out[i] = rmsd;
      }
    }
    ptm_uninitialize_local(handle);
  }
  return 0;
}
