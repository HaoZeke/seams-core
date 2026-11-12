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

#ifndef SEAMS_INTERFACE_ORDER_H_
#define SEAMS_INTERFACE_ORDER_H_

#include <mol_sys.hpp>

#include <string>
#include <vector>

/** @file interface_order.hpp
 *  @brief The two published CHILL+ rules on an interface, and the ring
 *  labels of confined ice.
 *
 *  Thresholds are the 2015 values. Nothing here is fitted.
 *
 *  Nguyen and Molinero, J. Phys. Chem. B 119 (2015) 9369, equation 1
 *  and Table 1: q3 is the mean of Y_3m over the four nearest neighbours,
 *  a staggered bond is c(i,j) <= -0.8, an eclipsed bond is
 *  -0.35 <= c(i,j) <= 0.25, a crystal row has four neighbours inside
 *  the cutoff, and an interfacial row keeps the neighbour test in
 *  footnote a of that table.
 *
 *  OVITO 3.17.2 Chill+ (manual page chill_plus) and the mdapy kernel
 *  that was cross-checked against OVITO's ChillPlusModifier: q3 runs
 *  over every neighbour inside the cutoff, a staggered bond is
 *  c(i,j) < -0.8, an eclipsed bond is -0.35 < c(i,j) < 0.25, and an
 *  interfacial row is two staggered bonds, or three staggered and no
 *  eclipsed bond, with no neighbour test.
 */

namespace iface {

enum class Cause { agree = 0, neighbour = 1, q3Shell = 2, openEdge = 3, mixed = 4 };

enum class Kind { slab = 0, nucleus = 1 };

struct Front {
  double position = 0; //! 0.5 crossing, or the nucleus radius
  double width = 0;    //! 10-90 width along the same path
  bool defined = false;
};

struct Profile {
  double crystal = 0;      //! (cubic + hexagonal) / N
  double interfacial = 0;  //! interfacial / N
  double solid = 0;        //! crystal + interfacial, over N
  double thickness = 0;    //! ice span between the two 0.5 crossings, or 2 R
  std::vector<Front> fronts;
  bool defined = false;
};

struct Bin {
  double mid = 0;
  int n = 0;
  int disagree = 0;
  double crystalPub = 0;
  double solidPub = 0;
  double crystalOvi = 0;
  double solidOvi = 0;
};

struct Study {
  int n = 0;
  int nDisagree = 0;
  int causeCount[5] = {};
  int countPub[9] = {};
  int countOvi[9] = {};
  int transition[9][9] = {};
  int shellDisagree[8] = {}; //! 0..6 and 7 meaning 7 or more
  int shellAgree[8] = {};
  int ringDisagree[3] = {}; //! atoms in at least one ring of size 4, 5, 6
  int ringAgree[3] = {};
  int prismAtoms = 0;
  int prismBlocks = 0;
  int bilayerAtoms = 0;
  int polygonAtoms = 0;
  int thinSheets = 0;
  int diamondOther = 0;
  int diamondCubic = 0;
  int diamondHex = 0;
  int cnaOther = 0;
  int cnaCrystal = 0; //! FCC, HCP, icosahedral or BCC under adaptive CNA
  int ptmPass = 0;
  int ptmPassDisagree = 0;
  double q6Disagree = 0;
  double q6Agree = 0;
  double depthDisagree = 0; //! mean distance to the nearest front
  double depthAgree = 0;
  int depthDisagreeN = 0;
  int depthAgreeN = 0;
  Kind kind = Kind::nucleus;
  int axis = 2;
  Profile solidPub;    //! cubic + hexagonal + interfacial
  Profile crystalPub;  //! cubic + hexagonal
  Profile solidOvi;
  Profile crystalOvi;
  std::vector<Bin> bins;
};

//! Both rules, the ring census, and the interface or nucleus profile.
//! Cutoff and bond thresholds stay at the published values.
[[nodiscard]] Study studyCloud(
    molSys::PointCloud<molSys::Point<double>, double> cloud, double cutoff = 3.5,
    double binWidth = 1.0);

[[nodiscard]] Study studyDump(const std::string &path, int frame, double cutoff = 3.5,
                              double binWidth = 1.0);

[[nodiscard]] std::string formatStudy(const Study &study, const std::string &tag);

//! Ordinary least squares slope. Empty or a single point returns 0.
[[nodiscard]] double slope(const std::vector<double> &x, const std::vector<double> &y);

struct GrowthFrame {
  int frame = 0;
  double timePs = 0;
  Study study;
};

struct Growth {
  std::vector<GrowthFrame> frames;
  double rateSolidPub = 0; //! d(thickness)/dt in A/ps
  double rateSolidOvi = 0;
  double rateCrystalPub = 0;
  double rateCrystalOvi = 0;
};

//! Every frame of a dump. time of frame i is (i - 1) * dumpEveryPs.
[[nodiscard]] Growth studyTrajectory(const std::string &path, double dumpEveryPs,
                                     double cutoff = 3.5);

[[nodiscard]] std::string formatGrowth(const Growth &growth);

//! Ideal oxygen frameworks at a nearest-neighbour distance of `bond`.
//! The box is orthorhombic and large enough that a 3.5 A shell does not
//! see a lateral image.
[[nodiscard]] molSys::PointCloud<molSys::Point<double>, double>
nanotube(int ring, int stacks, double bond = 2.75);

[[nodiscard]] molSys::PointCloud<molSys::Point<double>, double>
honeycombSheet(int nx, int ny, int layers, double bond = 2.75);

} // namespace iface

#endif
