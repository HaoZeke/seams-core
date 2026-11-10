#include <bop.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: chill_labels chill|chillplus file frame type\n";
    return 2;
  }
  const std::string mode = argv[1];
  const std::string path = argv[2];
  const int frame = std::atoi(argv[3]);
  const int type = std::atoi(argv[4]);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(path, frame, cloud, type);
  if (cloud.nop <= 0) {
    std::cerr << "nop " << cloud.nop << "\n";
    return 2;
  }
  const auto nList = nneigh::neighListO(3.5, cloud, type);
  if (mode == "chill") {
    chill::getCorrel(cloud, nList, false);
    chill::getIceTypeNoPrint(cloud, nList, false);
  } else if (mode == "chillplus") {
    chill::getCorrelPlus(cloud, nList, false);
    chill::getIceTypePlusNoPrint(cloud, nList, false);
  } else {
    std::cerr << "unknown mode\n";
    return 2;
  }
  for (int i = 0; i < cloud.nop; ++i) {
    std::cout << static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType) << "\n";
  }
  return 0;
}
