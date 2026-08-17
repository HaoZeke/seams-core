/*
** Probe whether N frames of the TUM ice-score working set fit on the
** GPU, then time a cold launch and every warm repeat (mean and std).
**   bench_gpu_batch TRAJ [nFrames] [atomType] [repeats]
**   bench_gpu_batch synth:NATOMS [nFrames] [atomType] [repeats]
** synth: uses the same jittered mW-density lattice as bench_scaling.
*/

#include <gpu_batch.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr double kNumberDensity = 0.0332;

molSys::PointCloud<molSys::Point<double>, double> makeCloud(int nAtoms,
                                                            int typeI) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double boxLength = std::cbrt(nAtoms / kNumberDensity);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.pts.reserve(nAtoms);
  unsigned long long state = 88172645463325252ULL;
  auto jitter = [&state, spacing]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * 0.3 * spacing;
  };
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> point;
    point.type = typeI;
    point.atomID = i + 1;
    point.molID = i + 1;
    point.x = ((i % perSide) + 0.5) * spacing + jitter();
    point.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    point.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(point);
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int want = argc > 2 ? std::atoi(argv[2]) : 11;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;
  const int repeats = argc > 4 ? std::max(1, std::atoi(argv[4])) : 5;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const bool synth = traj.rfind("synth:", 0) == 0;
  if (synth) {
    const int nSynth = std::atoi(traj.c_str() + 6);
    if (nSynth <= 0) {
      std::fprintf(stderr, "synth: needs a positive atom count\n");
      return 1;
    }
    cloud = makeCloud(nSynth, typeI);
  } else {
    cloud = sinp::readLammpsTrjO(traj, 1, cloud, typeI);
  }
  if (cloud.nop == 0) {
    std::fprintf(stderr, "empty %s\n", traj.c_str());
    return 1;
  }
  const int nAtoms = cloud.nop;
  std::vector<double> xyz(static_cast<std::size_t>(want) *
                          static_cast<std::size_t>(nAtoms) * 3);
  const bool tiltDump = cloud.box.size() >= 6;
  std::vector<double> box;
  std::vector<double> boxLow;
  int nBox = 3;
  if (tiltDump) {
    nBox = static_cast<int>(cloud.box.size());
    box = cloud.box;
    boxLow = cloud.boxLow;
  } else {
    box.assign(static_cast<std::size_t>(want) * 3, 0.0);
  }
  int got = 0;
  for (int f = 1; f <= want; ++f) {
    if (!synth) {
      cloud = sinp::readLammpsTrjO(traj, f, cloud, typeI);
      if (cloud.nop != nAtoms) {
        break;
      }
    }
    const auto base = static_cast<std::size_t>(got) *
                      static_cast<std::size_t>(nAtoms) * 3;
    for (int i = 0; i < nAtoms; ++i) {
      xyz[base + static_cast<std::size_t>(i) * 3 + 0] = cloud.pts[i].x;
      xyz[base + static_cast<std::size_t>(i) * 3 + 1] = cloud.pts[i].y;
      xyz[base + static_cast<std::size_t>(i) * 3 + 2] = cloud.pts[i].z;
    }
    if (!tiltDump) {
      box[static_cast<std::size_t>(got) * 3 + 0] = cloud.box[0];
      box[static_cast<std::size_t>(got) * 3 + 1] = cloud.box[1];
      box[static_cast<std::size_t>(got) * 3 + 2] = cloud.box[2];
    }
    ++got;
  }

  const auto plan = gpu::planBatch(nAtoms, got);
  std::printf("device %s\n",
              plan.device.available ? plan.device.name.c_str() : "none");
  std::printf("reason %s\n", plan.reason.c_str());
  std::printf("sm %d.%d\n", plan.device.computeMajor,
              plan.device.computeMinor);
  std::printf("total_bytes %zu\n", plan.device.totalBytes);
  std::printf("free_bytes %zu\n", plan.device.freeBytes);
  std::printf("nAtoms %d\n", nAtoms);
  std::printf("requested_frames %d\n", got);
  std::printf("bytes_per_frame %zu\n",
              gpu::estimateFootprint(nAtoms, 1).totalBytes);
  std::printf("max_resident_frames %d\n",
              gpu::maxResidentFrames(plan.device, nAtoms));
  std::printf("plan_frames %d\n", plan.frames);
  std::printf("resident %s\n", plan.resident ? "yes" : "no");
  if (!plan.resident) {
    return 0;
  }
  const auto printRun = [](const char *tag, const gpu::BatchResult &r) {
    std::printf("%s_upload_ms %.3f\n%s_compute_ms %.3f\n%s_download_ms %.3f\n",
                tag, r.uploadMs, tag, r.computeMs, tag, r.downloadMs);
  };

  const auto cold = tiltDump
                         ? gpu::analyzeResident(xyz.data(), box.data(), nAtoms,
                                                got, 5.5, boxLow.data(), nBox)
                         : gpu::analyzeResident(xyz.data(), box.data(), nAtoms,
                                                got);
  if (!cold.error.empty()) {
    std::printf("error %s\n", cold.error.c_str());
    return 1;
  }
  printRun("cold", cold);

  gpu::BatchResult last = cold;
  double sum = 0.0;
  double sumsq = 0.0;
  double lo = std::numeric_limits<double>::max();
  double hi = 0.0;
  for (int i = 0; i < repeats; ++i) {
    const auto r =
        tiltDump ? gpu::analyzeResident(xyz.data(), box.data(), nAtoms, got,
                                        5.5, boxLow.data(), nBox)
                 : gpu::analyzeResident(xyz.data(), box.data(), nAtoms, got);
    if (!r.error.empty()) {
      std::printf("error %s\n", r.error.c_str());
      return 1;
    }
    last = r;
    const double t = r.computeMs;
    std::printf("warm_rep %d %.3f\n", i, t);
    sum += t;
    sumsq += t * t;
    if (t < lo) {
      lo = t;
    }
    if (t > hi) {
      hi = t;
    }
  }
  const double mean = sum / static_cast<double>(repeats);
  const double var = (repeats > 1)
                         ? (sumsq - sum * sum / static_cast<double>(repeats)) /
                               static_cast<double>(repeats - 1)
                         : 0.0;
  const double stdev = std::sqrt(var < 0.0 ? 0.0 : var);
  std::printf("warm_n %d\nwarm_mean_ms %.3f\nwarm_std_ms %.3f\n"
              "warm_min_ms %.3f\nwarm_max_ms %.3f\n",
              repeats, mean, stdev, lo, hi);
  printRun("warm", last);
  int nHc = 0;
  int nDdc = 0;
  int nSix = 0;
  for (int v : last.atomHc) {
    nHc += v;
  }
  for (int v : last.atomDdc) {
    nDdc += v;
  }
  for (int v : last.nRings) {
    nSix += v;
  }
  std::printf("ice_hc_atoms %d\nice_ddc_atoms %d\n", nHc, nDdc);
  std::printf("six_rings %d\nrings_dropped %d\n", nSix, last.ringsDropped);
  return 0;
}
