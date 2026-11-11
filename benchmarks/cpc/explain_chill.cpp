/* Copyright (c) d-SEAMS developers
**
** Per-molecule CHILL+ bond record for the liquid-agreement explanation.
**
** Prints one row: index shell nbonds stag ecl other label cij...
** Without SEAMS_TIP this is main's getCorrelPlus(..., 4) and
** getIceTypePlusNoPrint. With SEAMS_TIP, argv rule is chill+ or
** ovito-chill+ and applyChillPlus runs that rule.
*/

#include <bop.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#ifdef SEAMS_TIP
#include <chill_plus_scheme.hpp>
#endif

#include <cstdlib>
#include <iostream>
#include <string>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

int shellCount(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
               const std::vector<std::vector<int>> &nList, int iatom) {
  if (iatom < 0 || static_cast<std::size_t>(iatom) >= nList.size()) {
    return 0;
  }
  int shell = 0;
  const auto &row = nList[static_cast<std::size_t>(iatom)];
  for (std::size_t j = 1; j < row.size(); ++j) {
    if (cloud.idIndexMap.find(row[j]) != cloud.idIndexMap.end()) {
      ++shell;
    }
  }
  return shell;
}

const char *nameOf(molSys::atom_state_type kind) {
  switch (kind) {
  case molSys::atom_state_type::cubic:
    return "cubic";
  case molSys::atom_state_type::hexagonal:
    return "hexagonal";
  case molSys::atom_state_type::water:
    return "water";
  case molSys::atom_state_type::interfacial:
    return "interfacial";
  case molSys::atom_state_type::clathrate:
    return "clathrate";
  case molSys::atom_state_type::interClathrate:
    return "interClathrate";
  case molSys::atom_state_type::reCubic:
    return "reCubic";
  case molSys::atom_state_type::reHex:
    return "reHex";
  default:
    return "unclassified";
  }
}

} // namespace

int main(int argc, char **argv) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
  if (argc < 2) {
    std::cerr << "usage: explain_chill FILE [chill+|ovito-chill+]\n";
    return 2;
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(argv[1], 1, cloud, 1);
  if (cloud.nop <= 0) {
    std::cerr << "nop 0\n";
    return 2;
  }
  const auto nList = nneigh::neighListO(3.5, cloud, 1);
#ifdef SEAMS_TIP
  const std::string rule = argc > 2 ? argv[2] : "chill+";
  chill::applyChillPlus(cloud, nList, chill::chillPlusRule(rule), false);
#else
  chill::getCorrelPlus(cloud, nList, false, 4);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
#endif
  std::cout.setf(std::ios::fixed);
  std::cout.precision(8);
  std::cout << "index shell nbonds stag ecl other label cij\n";
  for (int i = 0; i < cloud.nop; i++) {
    const auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    int stag = 0;
    int ecl = 0;
    int other = 0;
    for (const auto &bond : pt.c_ij) {
      if (bond.classifier == molSys::bond_type::staggered) {
        ++stag;
      } else if (bond.classifier == molSys::bond_type::eclipsed) {
        ++ecl;
      } else {
        ++other;
      }
    }
    std::cout << i << " " << shellCount(cloud, nList, i) << " " << pt.c_ij.size()
              << " " << stag << " " << ecl << " " << other << " " << nameOf(pt.iceType);
    for (const auto &bond : pt.c_ij) {
      std::cout << " " << bond.c_value;
    }
    std::cout << "\n";
  }
  return 0;
}
