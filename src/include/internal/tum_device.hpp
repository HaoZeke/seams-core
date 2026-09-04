#ifndef SEAMS_TUM_DEVICE_H_
#define SEAMS_TUM_DEVICE_H_

/** @file tum_device.hpp
 *  @brief Device-safe TUM ice-score pieces: hop-bound primitive six-rings
 *  and claim-free HC/DDC affiliation.
 *
 *  No STL containers, no function-local statics. Host OpenMP and OpenMP
 *  target regions call the same functions. CHILL and q_lm are not here.
 *
 *  A six-cycle is primitive when the hop-bounded Franzblau test holds:
 *  every pair of non-adjacent members is at least as far through the
 *  graph as around the ring. On a 4-regular bond graph that is no
 *  chords (ring hops 2) and no common neighbour of opposite vertices
 *  (ring hops 3).
 */

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp declare target
#endif

namespace tum {
namespace device {

inline bool bonded(const int *deg, const int *cols, int nAtoms, int kMax,
                   int a, int b) {
  if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) {
    return false;
  }
  const int d = deg[a];
  const int row = a * kMax;
  for (int t = 0; t < d; ++t) {
    if (cols[row + t] == b) {
      return true;
    }
  }
  return false;
}

inline bool inSix(const int *r, int atom) {
  for (int t = 0; t < 6; ++t) {
    if (r[t] == atom) {
      return true;
    }
  }
  return false;
}

inline bool shareAtoms(const int *a, const int *b) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) {
      return true;
    }
  }
  return false;
}

inline int commonCount(const int *a, const int *b) {
  int n = 0;
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) {
      ++n;
    }
  }
  return n;
}

inline bool commonInThree(const int *a, const int *b, const int *c) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i]) && inSix(c, a[i])) {
      return true;
    }
  }
  return false;
}

inline bool shareNeigh(const int *deg, const int *cols, int nAtoms, int kMax,
                       int a, int b) {
  const int da = deg[a];
  const int ra = a * kMax;
  for (int t = 0; t < da; ++t) {
    if (bonded(deg, cols, nAtoms, kMax, b, cols[ra + t])) {
      return true;
    }
  }
  return false;
}

// Franzblau SP on a 6-cycle: no chords, opposite vertices at graph
// distance 3 (not bonded and no common neighbour).
inline bool primitiveSix(const int *r, const int *deg, const int *cols,
                         int nAtoms, int kMax) {
  const int d2[6][2] = {{0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {5, 1}};
  for (int t = 0; t < 6; ++t) {
    if (bonded(deg, cols, nAtoms, kMax, r[d2[t][0]], r[d2[t][1]])) {
      return false;
    }
  }
  const int opp[3][2] = {{0, 3}, {1, 4}, {2, 5}};
  for (int t = 0; t < 3; ++t) {
    const int u = r[opp[t][0]];
    const int v = r[opp[t][1]];
    if (bonded(deg, cols, nAtoms, kMax, u, v)) {
      return false;
    }
    if (shareNeigh(deg, cols, nAtoms, kMax, u, v)) {
      return false;
    }
  }
  return true;
}

inline bool basalNeighbours(const int *deg, const int *cols, int nAtoms,
                            int kMax, int n1, int n2, int atomOne,
                            int atomTwo) {
  const bool n1one = bonded(deg, cols, nAtoms, kMax, atomOne, n1);
  const bool n1two = bonded(deg, cols, nAtoms, kMax, atomTwo, n1);
  if (!n1one && !n1two) {
    return false;
  }
  if (n1one) {
    return bonded(deg, cols, nAtoms, kMax, atomTwo, n2);
  }
  return bonded(deg, cols, nAtoms, kMax, atomOne, n2);
}

inline bool notNeighboursOfRing(const int *deg, const int *cols, int nAtoms,
                                int kMax, const int *trip, const int *ring) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (bonded(deg, cols, nAtoms, kMax, ring[j], trip[i])) {
        return false;
      }
    }
  }
  return true;
}

inline bool basalConditions(const int *deg, const int *cols, int nAtoms,
                            int kMax, const int *b1, const int *b2) {
  int kIndex = -1;
  int compare1 = 0;
  int compare2 = 0;
  bool l1n = false;
  bool l2n = false;
  const int l1 = b1[0];
  const int l2 = b1[1];
  for (int k = 0; k < 6; ++k) {
    const int mk = b2[k];
    if (bonded(deg, cols, nAtoms, kMax, l1, mk)) {
      compare1 = b1[2];
      compare2 = b1[4];
      kIndex = k;
      l1n = true;
      break;
    }
    if (bonded(deg, cols, nAtoms, kMax, l2, mk)) {
      compare1 = b1[3];
      compare2 = b1[5];
      kIndex = k;
      l2n = true;
      break;
    }
  }
  if (!l1n && !l2n) {
    return false;
  }
  int evenT[3];
  int oddT[3];
  int ie = 0;
  int io = 0;
  for (int k = 0; k <= 5; ++k) {
    int ck = kIndex + k;
    if (ck >= 6) {
      ck -= 6;
    }
    if (k % 2 == 0) {
      evenT[ie++] = b2[ck];
    } else {
      oddT[io++] = b2[ck];
    }
  }
  if (!basalNeighbours(deg, cols, nAtoms, kMax, evenT[1], evenT[2], compare1,
                       compare2)) {
    return false;
  }
  return notNeighboursOfRing(deg, cols, nAtoms, kMax, oddT, b1);
}

inline bool basalPair(const int *deg, const int *cols, int nAtoms, int kMax,
                      const int *bi, const int *bj) {
  if (shareAtoms(bi, bj)) {
    return false;
  }
  return basalConditions(deg, cols, nAtoms, kMax, bi, bj);
}

inline int firstRingThrough(const int *A, int nA, const int *B, int nB,
                            const int *C, int nC, int skipA, int skipB) {
  int i = 0;
  int j = 0;
  int k = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i];
    const int y = B[j];
    const int z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB) {
        return x;
      }
      ++i;
      ++j;
      ++k;
      continue;
    }
    int lo = x;
    if (y < lo) {
      lo = y;
    }
    if (z < lo) {
      lo = z;
    }
    if (x == lo) {
      ++i;
    }
    if (y == lo) {
      ++j;
    }
    if (z == lo) {
      ++k;
    }
  }
  return -1;
}

inline int ringsThrough(const int *A, int nA, const int *B, int nB,
                        const int *C, int nC, int skipA, int skipB, int *out,
                        int cap) {
  int i = 0;
  int j = 0;
  int k = 0;
  int n = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i];
    const int y = B[j];
    const int z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB && n < cap) {
        out[n++] = x;
      }
      ++i;
      ++j;
      ++k;
      continue;
    }
    int lo = x;
    if (y < lo) {
      lo = y;
    }
    if (z < lo) {
      lo = z;
    }
    if (x == lo) {
      ++i;
    }
    if (y == lo) {
      ++j;
    }
    if (z == lo) {
      ++k;
    }
  }
  return n;
}

// Enumerate primitive six-rings whose lowest-indexed member is src.
// Writes at most maxPer rings into localRings[src * maxPer * 6].
// Returns the number written; if more exist, returns maxPer + 1.
inline int enumSixFrom(int src, const int *deg, const int *cols, int nAtoms,
                       int kMax, int maxPer, int *localRings) {
  int found = 0;
  const int di = deg[src];
  if (di < 2) {
    return 0;
  }
  const int irow = src * kMax;
  for (int ia = 0; ia < di; ++ia) {
    const int a = cols[irow + ia];
    if (a < 0 || a >= nAtoms) {
      continue;
    }
    for (int ib = ia + 1; ib < di; ++ib) {
      const int b = cols[irow + ib];
      if (b < 0 || b >= nAtoms) {
        continue;
      }
      const int da = deg[a];
      const int db = deg[b];
      const int arow = a * kMax;
      const int brow = b * kMax;
      for (int ix = 0; ix < da; ++ix) {
        const int x = cols[arow + ix];
        if (x == src || x == b || x < 0 || x >= nAtoms) {
          continue;
        }
        for (int iy = 0; iy < db; ++iy) {
          const int y = cols[brow + iy];
          if (y == src || y == a || y == x || y < 0 || y >= nAtoms) {
            continue;
          }
          const int dx = deg[x];
          const int xrow = x * kMax;
          for (int iz = 0; iz < dx; ++iz) {
            const int z = cols[xrow + iz];
            if (z == src || z == a || z == b || z == y || z < 0 ||
                z >= nAtoms) {
              continue;
            }
            if (!bonded(deg, cols, nAtoms, kMax, z, y)) {
              continue;
            }
            const int cyc[6] = {src, a, x, z, y, b};
            int mn = src;
            for (int t = 1; t < 6; ++t) {
              if (cyc[t] < mn) {
                mn = cyc[t];
              }
            }
            if (mn != src) {
              continue;
            }
            if (!primitiveSix(cyc, deg, cols, nAtoms, kMax)) {
              continue;
            }
            if (found >= maxPer) {
              return maxPer + 1;
            }
            int *dest = localRings + (src * maxPer + found) * 6;
            for (int t = 0; t < 6; ++t) {
              dest[t] = cyc[t];
            }
            ++found;
          }
        }
      }
    }
  }
  return found;
}

inline bool isPrismaticOf(int k, int i, int j, const int *ringAtoms,
                          int nRings, const int *throughCount,
                          const int *through, int nAtoms, int maxPer) {
  if (k < 0 || i < 0 || j < 0 || k >= nRings || i >= nRings || j >= nRings) {
    return false;
  }
  const int *bi = ringAtoms + i * 6;
  const int *bj = ringAtoms + j * 6;
  const int *bk = ringAtoms + k * 6;
  for (int q = 0; q < 6; ++q) {
    int trip[3];
    for (int m = 0; m < 3; ++m) {
      trip[m] = bi[(q + m) % 6];
    }
    if (trip[0] < 0 || trip[0] >= nAtoms || trip[1] < 0 || trip[1] >= nAtoms ||
        trip[2] < 0 || trip[2] >= nAtoms) {
      continue;
    }
    const int nA = throughCount[trip[0]];
    const int nB = throughCount[trip[1]];
    const int nC = throughCount[trip[2]];
    const int *A = through + trip[0] * maxPer;
    const int *B = through + trip[1] * maxPer;
    const int *C = through + trip[2] * maxPer;
    int cand[16];
    const int nc = ringsThrough(A, nA, B, nB, C, nC, i, j, cand, 16);
    for (int c = 0; c < nc; ++c) {
      if (cand[c] != k) {
        continue;
      }
      int rest[3];
      int nrst = 0;
      for (int u = 0; u < 6; ++u) {
        if (!inSix(trip, bk[u]) && nrst < 3) {
          rest[nrst++] = bk[u];
        }
      }
      if (nrst == 3 && commonCount(rest, bj) == 3) {
        return true;
      }
    }
  }
  return false;
}

inline bool hcAffiliated(int i, const int *ringAtoms, int nRings,
                         const int *throughCount, const int *through,
                         int nAtoms, int maxPer, const int *deg,
                         const int *cols, int kMax) {
  if (i < 0 || i >= nRings) {
    return false;
  }
  const int *bi = ringAtoms + i * 6;
  for (int slot = 0; slot < 2; ++slot) {
    const int anchor = bi[slot];
    if (anchor < 0 || anchor >= nAtoms) {
      continue;
    }
    const int da = deg[anchor];
    const int arow = anchor * kMax;
    for (int n = 0; n < da; ++n) {
      const int nb = cols[arow + n];
      if (nb < 0 || nb >= nAtoms) {
        continue;
      }
      const int nr = throughCount[nb];
      const int *row = through + nb * maxPer;
      for (int t = 0; t < nr; ++t) {
        const int j = row[t];
        if (j == i || j < 0 || j >= nRings) {
          continue;
        }
        const int *bj = ringAtoms + j * 6;
        if (basalPair(deg, cols, nAtoms, kMax, bi, bj) ||
            basalPair(deg, cols, nAtoms, kMax, bj, bi)) {
          return true;
        }
      }
    }
  }

  int adj[96];
  int nAdj = 0;
  for (int s = 0; s < 6; ++s) {
    const int atom = bi[s];
    if (atom < 0 || atom >= nAtoms) {
      continue;
    }
    const int nr = throughCount[atom];
    const int *row = through + atom * maxPer;
    for (int t = 0; t < nr; ++t) {
      const int r = row[t];
      if (r == i || r < 0 || r >= nRings) {
        continue;
      }
      bool seen = false;
      for (int p = 0; p < nAdj; ++p) {
        if (adj[p] == r) {
          seen = true;
          break;
        }
      }
      if (!seen && nAdj < 96) {
        adj[nAdj++] = r;
      }
    }
    const int da = deg[atom];
    const int arow = atom * kMax;
    for (int n = 0; n < da; ++n) {
      const int nb = cols[arow + n];
      if (nb < 0 || nb >= nAtoms) {
        continue;
      }
      const int nr2 = throughCount[nb];
      const int *row2 = through + nb * maxPer;
      for (int t = 0; t < nr2; ++t) {
        const int r = row2[t];
        if (r == i || r < 0 || r >= nRings) {
          continue;
        }
        bool seen = false;
        for (int p = 0; p < nAdj; ++p) {
          if (adj[p] == r) {
            seen = true;
            break;
          }
        }
        if (!seen && nAdj < 96) {
          adj[nAdj++] = r;
        }
      }
    }
  }

  for (int a = 0; a < nAdj; ++a) {
    const int p = adj[a];
    const int *bp = ringAtoms + p * 6;
    for (int slot = 0; slot < 2; ++slot) {
      const int anchor = bp[slot];
      if (anchor < 0 || anchor >= nAtoms) {
        continue;
      }
      const int da = deg[anchor];
      const int arow = anchor * kMax;
      for (int n = 0; n < da; ++n) {
        const int nb = cols[arow + n];
        if (nb < 0 || nb >= nAtoms) {
          continue;
        }
        const int nr = throughCount[nb];
        const int *row = through + nb * maxPer;
        for (int t = 0; t < nr; ++t) {
          const int q = row[t];
          if (q == i || q == p || q < 0 || q >= nRings) {
            continue;
          }
          bool inAdj = false;
          for (int u = 0; u < nAdj; ++u) {
            if (adj[u] == q) {
              inAdj = true;
              break;
            }
          }
          if (!inAdj) {
            continue;
          }
          const int *bq = ringAtoms + q * 6;
          if (!basalPair(deg, cols, nAtoms, kMax, bp, bq)) {
            continue;
          }
          if (isPrismaticOf(i, p, q, ringAtoms, nRings, throughCount, through,
                            nAtoms, maxPer)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

// Equatorial DDC test. On a pass, writes the six peripheral ring
// indices into peri[6] and returns true.
inline bool equatorialPass(int i, const int *ringAtoms, int nRings,
                           const int *throughCount, const int *through,
                           int nAtoms, int maxPer, const int *hc, int *peri) {
  if (i < 0 || i >= nRings) {
    return false;
  }
  if (hc[i]) {
    return false;
  }
  const int *bi = ringAtoms + i * 6;
  for (int m = 0; m < 6; ++m) {
    const int atom = bi[m];
    if (atom < 0 || atom >= nAtoms) {
      return false;
    }
    const int nr = throughCount[atom];
    int common = 0;
    const int *row = through + atom * maxPer;
    for (int t = 0; t < nr; ++t) {
      if (row[t] == i) {
        continue;
      }
      ++common;
    }
    if (common < 3) {
      return false;
    }
  }
  for (int k = 0; k < 6; ++k) {
    int trip[3];
    for (int t = 0; t < 3; ++t) {
      trip[t] = bi[(k + t) % 6];
    }
    if (trip[0] < 0 || trip[0] >= nAtoms || trip[1] < 0 || trip[1] >= nAtoms ||
        trip[2] < 0 || trip[2] >= nAtoms) {
      return false;
    }
    const int nA = throughCount[trip[0]];
    const int nB = throughCount[trip[1]];
    const int nC = throughCount[trip[2]];
    const int *A = through + trip[0] * maxPer;
    const int *B = through + trip[1] * maxPer;
    const int *C = through + trip[2] * maxPer;
    const int j = firstRingThrough(A, nA, B, nB, C, nC, i, -1);
    if (j < 0 || j >= nRings) {
      return false;
    }
    peri[k] = j;
  }
  const int *p0 = ringAtoms + peri[0] * 6;
  const int *p1 = ringAtoms + peri[1] * 6;
  const int *p2 = ringAtoms + peri[2] * 6;
  const int *p3 = ringAtoms + peri[3] * 6;
  const int *p4 = ringAtoms + peri[4] * 6;
  const int *p5 = ringAtoms + peri[5] * 6;
  if (!commonInThree(p0, p2, p4)) {
    return false;
  }
  if (!commonInThree(p1, p3, p5)) {
    return false;
  }
  if (commonCount(p0, p2) < 3 || commonCount(p1, p3) < 3 ||
      commonCount(p2, p4) < 3 || commonCount(p3, p5) < 3) {
    return false;
  }
  return true;
}

} // namespace device
} // namespace tum

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp end declare target
#endif

#endif
