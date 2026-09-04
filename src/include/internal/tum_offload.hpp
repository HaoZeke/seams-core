#ifndef SEAMS_TUM_OFFLOAD_H_
#define SEAMS_TUM_OFFLOAD_H_

#include <cage_affiliation.hpp>

#include <vector>

/** @file tum_offload.hpp
 *  @brief OpenMP-target TUM ice score: hop-bound six-rings and HC/DDC
 *  affiliation. CHILL+ stays on the host. Steinhardt offload is
 *  unchanged.
 *
 *  SEAMS_OFFLOAD=0 forces the host Franzblau and cageAffiliation path.
 *  Unset, or any value other than 0, uses the device when a target
 *  exists. The cage counts match the host path on a perfect crystal.
 */

namespace primitive {

/** Primitive six-rings on an index neighbour list (leading self).
 *  Device hop-bound enum when SEAMS_OFFLOAD is not 0 and a target
 *  exists; otherwise Franzblau filtered to size 6. */
std::vector<std::vector<int>>
sixRingNetwork(const std::vector<std::vector<int>> &nList);

} // namespace primitive

namespace ring {

struct TumIceScore {
  std::vector<std::vector<int>> rings;
  CageAffiliation affiliation;
  std::vector<int> atomHc;
  std::vector<int> atomDdc;
  int ringsDropped = 0;
  bool usedDevice = false;
};

/** Six-rings plus HC/DDC affiliation. Device path when offload is
 *  live; otherwise host Franzblau and cageAffiliation. */
TumIceScore tumIceScore(const std::vector<std::vector<int>> &nList);

/** Device HC/DDC affiliation for six-rings. Returns false when the
 *  caller should use the host predicates (no device, SEAMS_OFFLOAD=0,
 *  a non-hexagon, or a capacity overflow). */
bool cageAffiliationOffload(const std::vector<std::vector<int>> &rings,
                            const std::vector<std::vector<int>> &nList,
                            CageAffiliation &out);

} // namespace ring

#endif
