//-----------------------------------------------------------------------------------
// d-SEAMS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __CAGE_H_
#define __CAGE_H_
#include <vector>

/** @file cage.hpp
 *  @brief File for cage types for topological network criteria.
 *
 */

/**
 *  @addtogroup cage
 *  @{
 */

/** @brief Functions for topological network criteria cage types.
 *         This namespace contains structs and enums for cage types
 *
 * Type definitions for atoms, rings and cages which are DDCs and HCs.
 *
 * ### Changelog ###
 *
 * - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 13, 2019
 * - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

// Namespace for cages
namespace cage {

// Type of a cage (a group of rings)
/** @enum cage::cageType
 * Qualifier for a cage.
 * According to the topological network criterion for DDCs and HCs
 *
 * @var cage::cageType HexC
 * The type for a hexagonal cage
 *
 * @var molSys::cageType DoubleDiaC
 * The type for a double-diamond cage
 */
enum cageType { HexC, DoubleDiaC };

// Type of ice for a particular atom. Dummy means that the atom is unclassified
// and is most probably water
/** @enum cage::iceType
 * Qualifier for an atom, based on whether it is part of cages,
 * according to the topological network criterion for DDCs and HCs
 *
 * @var cage::iceType dummy
 * Type for an atom which does not belong to any kind of cage
 *
 * @var cage::iceType hc
 * Type for an atom which belongs to an HC
 *
 * @var cage::iceType ddc
 * Type for an atom which belongs to a DDC
 *
 * @var cage::iceType mixed
 * Type for an atom which is part of a mixed ring, shared by both a DDC and an
 * HC
 *
 * @var cage::iceType pnc
 * Type for an atom which is part of a PNC
 *
 * @var cage::iceType mixed2
 * Type for an atom which is part of a mixed ring, shared by PNCs and DDCs/HCs
 *
 * @var cage::iceType clathS2
 * Type for an atom which is part of a clathrate S2 type cage 
 *
 */
enum iceType { dummy, hc, ddc, mixed, pnc, mixed2, clathS2 };

// Each DDC has one equatorial ring and 6 peripheral rings
// Each HC has two basal planes and 3 prismatic planes
/** @struct Cage
 * @brief This contains a cage, with the constituent rings
 *
 * Contains specifically the members:
 * - Cage classifier or qualifier, for each cage (can be a DDC or HC)
 * - Vector of rings in the cage
 */
struct Cage {
  cageType type;          //! type of the cage : can be DDC or HC
  std::vector<int> rings; //! coordinates
};

} // namespace cage

#endif // __CAGE_H_
