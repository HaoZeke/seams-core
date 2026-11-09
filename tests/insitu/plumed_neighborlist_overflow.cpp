// The neighbour fault on a 64000-atom cell is in PLUMED, not in this action.
//
// plumed 2.10.2 src/tools/NeighborList.cpp, getIndexPair, line 138:
//
//   unsigned K = unsigned(std::floor((std::sqrt(double(8*ii+1))+1)/2));
//
// ii is unsigned, so 8*ii is computed in 32 bits and wraps once ii reaches
// 2^29. A single list of n atoms has n*(n-1)/2 pairs. That product passes
// 2^29 at 32769 atoms. The guard at lines 106-110 refuses only a list above
// 1250000000 pairs, which is about 50000 atoms, so every list from 32769 to
// 50000 atoms walks off the position array inside Pbc::distance.
//
// This program evaluates that expression. It does not link PLUMED.

#include <cmath>
#include <cstdint>
#include <cstdio>

int main() {
  const unsigned n = 40000;
  const unsigned nallpairs = n * (n - 1u) / 2u;
  const unsigned guard = 1250000000u;
  unsigned worst = 0;
  // The first pair is the one whose ii is nallpairs-1, the overflowing end.
  unsigned bad = 0;
  unsigned example = 0;
  // The same walk as NeighborList::initialize. Stepping keeps it short and
  // still meets an index that is not an atom.
  for (unsigned ipair = 0; ipair < nallpairs; ipair += 1000u) {
    const unsigned ii = nallpairs - 1u - ipair;
    const unsigned k = static_cast<unsigned>(
        std::floor((std::sqrt(static_cast<double>(8u * ii + 1u)) + 1.0) / 2.0));
    const unsigned jj = ii - k * (k - 1u) / 2u;
    const unsigned first = n - 1u - k;
    const unsigned second = n - 1u - jj;
    if (k > worst) {
      worst = k;
    }
    if (first >= n || second >= n) {
      ++bad;
      example = second >= n ? second : first;
    }
  }
  const unsigned ii = nallpairs - 1u;
  const unsigned wrapped = 8u * ii;
  std::printf("n %u pairs %u guard %u wrapped_8ii %u worst_K %u out_of_range %u example %u\n", n,
              nallpairs, guard, wrapped, worst, bad, example);
  // The list is under the guard, the unsigned product wrapped, and an index
  // falls outside the atom list.
  const bool overflow = (8ull * ii) != wrapped;
  return (nallpairs < guard && overflow && bad > 0) ? 0 : 1;
}
