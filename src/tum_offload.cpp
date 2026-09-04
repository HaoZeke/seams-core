//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <tum_offload.hpp>

#include <franzblau.hpp>
#include <tum_device.hpp>

#include <algorithm>
#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

std::vector<std::vector<int>>
sixOf(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

void mapAtoms(const std::vector<std::vector<int>> &rings,
              const ring::CageAffiliation &aff, int nAtoms,
              std::vector<int> &atomHc, std::vector<int> &atomDdc) {
  atomHc.assign(static_cast<std::size_t>(nAtoms), 0);
  atomDdc.assign(static_cast<std::size_t>(nAtoms), 0);
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const bool hc = r < aff.hc.size() && aff.hc[r];
    const bool ddc = r < aff.ddc.size() && aff.ddc[r];
    for (const int a : rings[r]) {
      if (a < 0 || a >= nAtoms) {
        continue;
      }
      if (hc) {
        atomHc[static_cast<std::size_t>(a)] = 1;
      }
      if (ddc) {
        atomDdc[static_cast<std::size_t>(a)] = 1;
      }
    }
  }
}

#ifdef SEAMS_HAS_OFFLOAD
constexpr int kMaxPer = 16;

struct FlatGraph {
  std::vector<int> deg;
  std::vector<int> cols;
  int nAtoms = 0;
  int kMax = 0;
};

bool flattenIndexList(const std::vector<std::vector<int>> &nList,
                      FlatGraph &g) {
  g.nAtoms = static_cast<int>(nList.size());
  g.kMax = 0;
  if (g.nAtoms <= 0) {
    return false;
  }
  for (const auto &row : nList) {
    const int d = row.size() > 1 ? static_cast<int>(row.size()) - 1 : 0;
    if (d > g.kMax) {
      g.kMax = d;
    }
  }
  if (g.kMax < 2) {
    return false;
  }
  g.deg.assign(static_cast<std::size_t>(g.nAtoms), 0);
  g.cols.assign(static_cast<std::size_t>(g.nAtoms) *
                    static_cast<std::size_t>(g.kMax),
                -1);
  for (int i = 0; i < g.nAtoms; ++i) {
    const auto &row = nList[static_cast<std::size_t>(i)];
    int kept = 0;
    for (std::size_t j = 1; j < row.size(); ++j) {
      const int nb = row[j];
      if (nb < 0 || nb >= g.nAtoms || nb == i) {
        continue;
      }
      g.cols[static_cast<std::size_t>(i) * static_cast<std::size_t>(g.kMax) +
             static_cast<std::size_t>(kept)] = nb;
      ++kept;
    }
    g.deg[static_cast<std::size_t>(i)] = kept;
  }
  return true;
}

bool wantOffload() {
  const char *env = std::getenv("SEAMS_OFFLOAD");
  if (env != nullptr && env[0] == '0') {
    return false;
  }
  return omp_get_num_devices() > 0;
}

bool invertRings(const std::vector<int> &ringAtoms, int nRings, int nAtoms,
                 int maxPer, std::vector<int> &throughCount,
                 std::vector<int> &through) {
  throughCount.assign(static_cast<std::size_t>(nAtoms), 0);
  through.assign(static_cast<std::size_t>(nAtoms) *
                     static_cast<std::size_t>(maxPer),
                 -1);
  for (int r = 0; r < nRings; ++r) {
    for (int t = 0; t < 6; ++t) {
      const int a = ringAtoms[static_cast<std::size_t>(r) * 6 +
                              static_cast<std::size_t>(t)];
      if (a < 0 || a >= nAtoms) {
        continue;
      }
      const int slot = throughCount[static_cast<std::size_t>(a)];
      if (slot >= maxPer) {
        return false;
      }
      through[static_cast<std::size_t>(a) * static_cast<std::size_t>(maxPer) +
              static_cast<std::size_t>(slot)] = r;
      throughCount[static_cast<std::size_t>(a)] = slot + 1;
    }
  }
  for (int a = 0; a < nAtoms; ++a) {
    const int n = throughCount[static_cast<std::size_t>(a)];
    int *row = through.data() +
               static_cast<std::size_t>(a) * static_cast<std::size_t>(maxPer);
    for (int i = 1; i < n; ++i) {
      const int key = row[i];
      int p = i;
      while (p > 0 && row[p - 1] > key) {
        row[p] = row[p - 1];
        --p;
      }
      row[p] = key;
    }
  }
  return true;
}

bool sixRingsDevice(const FlatGraph &g, std::vector<int> &ringAtoms,
                    int &nRings, int &dropped) {
  const int n = g.nAtoms;
  const int kMax = g.kMax;
  const int maxPer = kMaxPer;
  std::vector<int> localCount(static_cast<std::size_t>(n), 0);
  std::vector<int> localRings(static_cast<std::size_t>(n) *
                                  static_cast<std::size_t>(maxPer) * 6,
                              -1);
  const int *deg = g.deg.data();
  const int *cols = g.cols.data();
  int *countP = localCount.data();
  int *ringsP = localRings.data();
  const int degN = n;
  const int colN = n * kMax;
  const int ringN = n * maxPer * 6;

#pragma omp target data map(to : deg[0 : degN], cols[0 : colN], n, kMax,      \
                                maxPer)                                       \
    map(from : countP[0 : degN], ringsP[0 : ringN])
  {
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; ++i) {
      countP[i] =
          tum::device::enumSixFrom(i, deg, cols, n, kMax, maxPer, ringsP);
    }
  }

  dropped = 0;
  nRings = 0;
  for (int i = 0; i < n; ++i) {
    if (localCount[static_cast<std::size_t>(i)] > maxPer) {
      ++dropped;
    } else {
      nRings += localCount[static_cast<std::size_t>(i)];
    }
  }
  if (dropped > 0) {
    return false;
  }
  ringAtoms.resize(static_cast<std::size_t>(nRings) * 6);
  int dest = 0;
  for (int i = 0; i < n; ++i) {
    const int got = localCount[static_cast<std::size_t>(i)];
    for (int s = 0; s < got; ++s) {
      const int src = (i * maxPer + s) * 6;
      for (int t = 0; t < 6; ++t) {
        ringAtoms[static_cast<std::size_t>(dest) * 6 +
                  static_cast<std::size_t>(t)] =
            localRings[static_cast<std::size_t>(src + t)];
      }
      ++dest;
    }
  }
  return true;
}

bool affiliateDevice(const FlatGraph &g, const std::vector<int> &ringAtoms,
                     int nRings, ring::CageAffiliation &out) {
  if (nRings <= 0) {
    out.hc.clear();
    out.ddc.clear();
    return true;
  }
  const int n = g.nAtoms;
  const int kMax = g.kMax;
  const int maxPer = kMaxPer;
  std::vector<int> throughCount;
  std::vector<int> through;
  if (!invertRings(ringAtoms, nRings, n, maxPer, throughCount, through)) {
    return false;
  }
  std::vector<int> hc(static_cast<std::size_t>(nRings), 0);
  std::vector<int> ddc(static_cast<std::size_t>(nRings), 0);
  const int *deg = g.deg.data();
  const int *cols = g.cols.data();
  const int *ringsP = ringAtoms.data();
  const int *tCount = throughCount.data();
  const int *tRow = through.data();
  int *hcP = hc.data();
  int *ddcP = ddc.data();
  const int degN = n;
  const int colN = n * kMax;
  const int ringN = nRings * 6;
  const int thruN = n * maxPer;

#pragma omp target data map(to : deg[0 : degN], cols[0 : colN],               \
                                ringsP[0 : ringN], tCount[0 : degN],          \
                                tRow[0 : thruN], n, kMax, maxPer, nRings)     \
    map(tofrom : hcP[0 : nRings], ddcP[0 : nRings])
  {
#pragma omp target teams distribute parallel for
    for (int i = 0; i < nRings; ++i) {
      hcP[i] = tum::device::hcAffiliated(i, ringsP, nRings, tCount, tRow, n,
                                         maxPer, deg, cols, kMax)
                   ? 1
                   : 0;
    }
#pragma omp target teams distribute parallel for
    for (int i = 0; i < nRings; ++i) {
      int peri[6];
      if (tum::device::equatorialPass(i, ringsP, nRings, tCount, tRow, n,
                                      maxPer, hcP, peri)) {
        ddcP[i] = 1;
        for (int t = 0; t < 6; ++t) {
          if (peri[t] >= 0 && peri[t] < nRings) {
            ddcP[peri[t]] = 1;
          }
        }
      }
    }
  }

  out.hc.assign(static_cast<std::size_t>(nRings), false);
  out.ddc.assign(static_cast<std::size_t>(nRings), false);
  for (int i = 0; i < nRings; ++i) {
    out.hc[static_cast<std::size_t>(i)] = hc[static_cast<std::size_t>(i)] != 0;
    out.ddc[static_cast<std::size_t>(i)] =
        ddc[static_cast<std::size_t>(i)] != 0;
  }
  return true;
}

bool allSix(const std::vector<std::vector<int>> &rings) {
  for (const auto &r : rings) {
    if (r.size() != 6) {
      return false;
    }
  }
  return true;
}

void packRings(const std::vector<std::vector<int>> &rings,
               std::vector<int> &flat) {
  flat.resize(rings.size() * 6);
  for (std::size_t i = 0; i < rings.size(); ++i) {
    for (int t = 0; t < 6; ++t) {
      flat[i * 6 + static_cast<std::size_t>(t)] = rings[i][static_cast<std::size_t>(t)];
    }
  }
}

void unpackRings(const std::vector<int> &flat, int nRings,
                 std::vector<std::vector<int>> &rings) {
  rings.assign(static_cast<std::size_t>(nRings), std::vector<int>(6));
  for (int i = 0; i < nRings; ++i) {
    for (int t = 0; t < 6; ++t) {
      rings[static_cast<std::size_t>(i)][static_cast<std::size_t>(t)] =
          flat[static_cast<std::size_t>(i) * 6 + static_cast<std::size_t>(t)];
    }
  }
}
#endif

} // namespace

std::vector<std::vector<int>>
primitive::sixRingNetwork(const std::vector<std::vector<int>> &nList) {
#ifdef SEAMS_HAS_OFFLOAD
  if (wantOffload()) {
    FlatGraph g;
    std::vector<int> flat;
    int nRings = 0;
    int dropped = 0;
    if (flattenIndexList(nList, g) && sixRingsDevice(g, flat, nRings, dropped) &&
        dropped == 0) {
      std::vector<std::vector<int>> rings;
      unpackRings(flat, nRings, rings);
      return rings;
    }
  }
#endif
  return sixOf(primitive::ringNetwork(nList, 6));
}

bool ring::cageAffiliationOffload(const std::vector<std::vector<int>> &rings,
                                  const std::vector<std::vector<int>> &nList,
                                  CageAffiliation &out) {
#ifdef SEAMS_HAS_OFFLOAD
  if (!wantOffload() || rings.empty() || !allSix(rings)) {
    return false;
  }
  FlatGraph g;
  if (!flattenIndexList(nList, g)) {
    return false;
  }
  std::vector<int> flat;
  packRings(rings, flat);
  return affiliateDevice(g, flat, static_cast<int>(rings.size()), out);
#else
  (void)rings;
  (void)nList;
  (void)out;
  return false;
#endif
}

ring::TumIceScore ring::tumIceScore(const std::vector<std::vector<int>> &nList) {
  TumIceScore out;
#ifdef SEAMS_HAS_OFFLOAD
  if (wantOffload()) {
    FlatGraph g;
    std::vector<int> flat;
    int nRings = 0;
    int dropped = 0;
    if (flattenIndexList(nList, g) &&
        sixRingsDevice(g, flat, nRings, dropped) && dropped == 0 &&
        affiliateDevice(g, flat, nRings, out.affiliation)) {
      unpackRings(flat, nRings, out.rings);
      out.ringsDropped = dropped;
      out.usedDevice = true;
      mapAtoms(out.rings, out.affiliation, g.nAtoms, out.atomHc, out.atomDdc);
      return out;
    }
  }
#endif
  out.rings = sixOf(primitive::ringNetwork(nList, 6));
  out.affiliation = cageAffiliation(out.rings, nList);
  out.usedDevice = false;
  mapAtoms(out.rings, out.affiliation, static_cast<int>(nList.size()),
           out.atomHc, out.atomDdc);
  return out;
}
