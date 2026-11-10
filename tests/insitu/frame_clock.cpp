/* Time one dump frame the way the command line does, after the read.
 * Prints read_s and analysis_s at full precision. The analysis is the
 * neighbour list, the primitive rings, CHILL+ and the cage census.
 */

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <chrono>
#include <cstdio>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: frame_clock dump.lammpstrj\n");
    return 2;
  }
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  Cloud cloud;
  const auto t0 = std::chrono::steady_clock::now();
  cloud = sinp::readLammpsTrj(argv[1], 1, cloud);
  const auto t1 = std::chrono::steady_clock::now();
  // The same three passes the rings, chill-plus and cages commands make,
  // each with its own neighbour list, and none of the process or the read.
  if (cloud.nop > 0) {
    {
      auto nList = nneigh::neighListO(3.5, cloud, 1);
      auto idx = nneigh::neighbourListByIndex(cloud, nList);
      (void)primitive::ringNetwork(idx, 8);
    }
    {
      auto nList = nneigh::neighListO(3.5, cloud, 1);
      chill::getCorrelPlus(cloud, nList, false);
      chill::getIceTypePlusNoPrint(cloud, nList, false);
    }
    {
      auto nList = nneigh::neighListO(3.5, cloud, 1);
      auto idx = nneigh::neighbourListByIndex(cloud, nList);
      std::vector<std::vector<int>> six;
      for (const auto &ring : primitive::ringNetwork(idx, 6)) {
        if (ring.size() == 6) {
          six.push_back(ring);
        }
      }
      const auto faces = ring::contractibleRings(six, cloud);
      (void)ring::cageAffiliation(faces, idx);
    }
  }
  const auto t2 = std::chrono::steady_clock::now();
  const auto sec = [](auto a, auto b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::printf("read_s %.17g\nanalysis_s %.17g\nn %d\n", sec(t0, t1), sec(t1, t2),
              cloud.nop);
  return 0;
}
