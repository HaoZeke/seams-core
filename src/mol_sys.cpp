//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <mol_sys.hpp>

/**
 * @details Function for clearing PointCloud if it is already
 *  filled. This should be called before every frame is read in.
 * @param[out] yCloud The cleared PointCloud
 */
molSys::PointCloud<molSys::Point<double>, double> molSys::clearPointCloud(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  //
  std::vector<molSys::Point<double>> tempPts;
  std::vector<double> tempBox;
  //
  std::vector<double> tempBox1;

  tempPts.swap(yCloud.pts);
  tempBox.swap(yCloud.box);
  tempBox1.swap(yCloud.boxLow);
  yCloud.idIndexMap.clear();
  yCloud.nop = 0;

  return yCloud;
}

/**
 * @details Function for creating an unordered map with the atomIDs in the
 *  pointCloud as the keys and the molecular IDs as the values
 */
std::unordered_map<int, int> molSys::createIDMolIDmap(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  std::unordered_map<int, int>
      idMolIDmap; // atom IDs as keys and mol IDs as values
  int iatomMolID; // molID of the current iatom
  int iatomID;    // atom ID of the current iatom

  // Loop through the atoms in yCloud
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    iatomID = yCloud.pts[iatom].atomID;   // atom ID
    iatomMolID = yCloud.pts[iatom].molID; // molecular ID
    // Update the unordered map
    idMolIDmap[iatomID] = iatomMolID;
  } // end of loop through every iatom in pointCloud

  return idMolIDmap;
}

/**
 * @details Function for creating an unordered map with the atomIDs in the
 *  pointCloud as the keys and the molecular IDs as the values. More than one atom
 * can have the same molecule ID. 
 */
std::unordered_multimap<int, int> molSys::createMolIDAtomIDMultiMap(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  std::unordered_multimap<int, int>
      molIDAtomIDmap; // atom IDs as keys and mol IDs as values
  int iatomMolID; // molID of the current iatom
  int iatomID;    // atom ID of the current iatom

  // Loop through the atoms in yCloud
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    iatomID = yCloud.pts[iatom].atomID;   // atom ID
    iatomMolID = yCloud.pts[iatom].molID; // molecular ID
    // Update the unordered multimap
    molIDAtomIDmap.emplace(iatomMolID,iatomID);
  } // end of loop through every iatom in pointCloud

  return molIDAtomIDmap;
}

/**
 * @details Function that returns a vector of vectors, which contains the
 *  hydrogen atoms for each molID in the oxygen atom pointCloud
 */
std::vector<std::vector<int>> molSys::hAtomMolList(
    molSys::PointCloud<molSys::Point<double>, double> &hCloud,
    molSys::PointCloud<molSys::Point<double>, double> &oCloud) {
  std::vector<std::vector<int>>
      hMolList; // the first column contains the molecular IDs, and the next
                // two elements in the row are the hydrogen bond atoms in the
                // molecule
  int iMolID;   // Current molecular ID
  int nHatoms;  // No. of h atoms found for a particular molID.

  std::unordered_map<int, std::vector<int>> hByMol;
  hByMol.reserve(static_cast<std::size_t>(std::max(hCloud.nop, 0)));
  for (int jatom = 0; jatom < hCloud.nop; jatom++) {
    hByMol[hCloud.pts[jatom].molID].push_back(jatom);
  }

  for (int iatom = 0; iatom < oCloud.nop; iatom++) {
    iMolID = oCloud.pts[iatom].molID;
    hMolList.push_back(std::vector<int>());
    hMolList[iatom].push_back(iMolID);
    nHatoms = 0;
    const auto found = hByMol.find(iMolID);
    if (found == hByMol.end()) {
      continue;
    }
    for (int jatom : found->second) {
      hMolList[iatom].push_back(jatom);
      nHatoms++;
      if (nHatoms == 2) {
        break;
      }
    }
  }

  return hMolList;
} // end of function

/**
 * @details Function for searching a vector of vectors for a particular
 * molecular ID, and
 * @returns the index found in molList
 * @returns -1 if not found
 */
int molSys::searchMolList(const std::vector<std::vector<int>> &molList,
                          int molIDtoFind) {
  int index = -1; // init invalid index

  for (int iatom = 0; iatom < molList.size(); iatom++) {
    // If the molecular ID is equal, return the index in the array
    if (molList[iatom][0] == molIDtoFind) {
      index = iatom;
      return index;
    } // end of check
  }   // end of looping through iatom

  return index;
}
