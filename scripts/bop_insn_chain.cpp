// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// One public steinhardtQl call on a chain of N atoms. The instruction
// counter is applied by the caller, not by this program.

#include <bop.hpp>

#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

int main(int argc, char **argv) {
  const int nAtoms = argc > 1 ? std::atoi(argv[1]) : 0;
  if (nAtoms < 2) {
    return 2;
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {static_cast<double>(nAtoms) * 2.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.reserve(static_cast<std::size_t>(nAtoms));
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(nAtoms));
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.x = static_cast<double>(i) * 1.5;
    p.y = 0.0;
    p.z = 0.0;
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
    nList[static_cast<std::size_t>(i)].push_back(p.atomID);
    if (i + 1 < nAtoms) {
      nList[static_cast<std::size_t>(i)].push_back(i + 2);
    }
  }
  const auto q = chill::steinhardtQl(cloud, nList, 6);
  volatile double sink = q.ql.empty() ? 0.0 : q.ql[0];
  return sink == -1.0 ? 3 : 0;
}
