/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Benchmark neighbour-list stages only (cutoff and k-nearest) on one frame.
** Intended for head-to-head timing against external engines.
*/

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <neighbours.hpp>
#include <seams_input.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Rows = std::vector<std::vector<int>>;

constexpr int kBenchmarkRepeats = 6;
constexpr int kBenchmarkWarmupCalls = 2;
// Derived from input/traj/mW_cubic.lammpstrj by
// tests/derive_mw_gr_window_midpoint.py.
constexpr double kMwSecondShellWindowMidpointAngstrom = 5.360;
constexpr std::uint64_t kFnvOffsetBasis = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

struct Args {
  std::string traj;
  int frame{1};
  int typeI{1};
  int reps{kBenchmarkRepeats};
  int warmup{kBenchmarkWarmupCalls};
  double cutoff{3.5};
  int k{4};
  double kCand{kMwSecondShellWindowMidpointAngstrom};
  int nx{1};
  int ny{1};
  int nz{1};
  int threads{0};
  std::string dumpCutoffEdges;
};

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

double quantile(const std::vector<double> &sorted, double q) {
  if (sorted.empty()) {
    return 0.0;
  }
  if (sorted.size() == 1) {
    return sorted[0];
  }
  const double pos = q * static_cast<double>(sorted.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double frac = pos - static_cast<double>(i);
  if (i + 1 >= sorted.size()) {
    return sorted.back();
  }
  return sorted[i] * (1.0 - frac) + sorted[i + 1] * frac;
}

struct Summary {
  double median{0.0};
  double q1{0.0};
  double q3{0.0};
};

Summary summarize(std::vector<double> values) {
  if (values.empty()) {
    return {};
  }
  std::sort(values.begin(), values.end());
  Summary out;
  out.median = quantile(values, 0.5);
  out.q1 = quantile(values, 0.25);
  out.q3 = quantile(values, 0.75);
  return out;
}

long long directedEdges(const Rows &rows) {
  long long total = 0;
  for (const auto &row : rows) {
    if (row.size() > 1) {
      total += static_cast<long long>(row.size() - 1);
    }
  }
  return total;
}

std::uint64_t hashRows(const Rows &rows) {
  std::uint64_t h = kFnvOffsetBasis;
  for (const auto &row : rows) {
    h ^= static_cast<std::uint64_t>(row.size());
    h *= kFnvPrime;
    for (const int v : row) {
      h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(v));
      h *= kFnvPrime;
    }
  }
  return h;
}

std::string joinSamples(const std::vector<double> &values) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(6);
  for (std::size_t i = 0; i < values.size(); i++) {
    if (i > 0) {
      oss << ",";
    }
    oss << values[i];
  }
  return oss.str();
}

std::string joinNames(const std::vector<std::string> &values) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < values.size(); i++) {
    if (i > 0) {
      oss << ",";
    }
    oss << values[i];
  }
  return oss.str();
}

Cloud replicateCloud(const Cloud &base, int nx, int ny, int nz) {
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    return base;
  }

  Cloud out;
  out.currentFrame = base.currentFrame;
  out.nop = base.nop * nx * ny * nz;

  double H[3][3];
  double origin[3];
  nneigh::dumpBoundsToH(base.box, base.boxLow, H, origin);
  const double ax = H[0][0], ay = 0.0, az = 0.0;
  const double bx = H[1][0], by = H[1][1], bz = 0.0;
  const double cx = H[2][0], cy = H[2][1], cz = H[2][2];

  out.pts.reserve(static_cast<std::size_t>(out.nop));
  int nextID = 1;
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        const double sx = static_cast<double>(ix);
        const double sy = static_cast<double>(iy);
        const double sz = static_cast<double>(iz);
        const double tx = sx * ax + sy * bx + sz * cx;
        const double ty = sx * ay + sy * by + sz * cy;
        const double tz = sx * az + sy * bz + sz * cz;
        for (const auto &p : base.pts) {
          molSys::Point<double> q = p;
          q.atomID = nextID;
          q.molID = nextID;
          q.x += tx;
          q.y += ty;
          q.z += tz;
          out.idIndexMap[nextID] = nextID - 1;
          out.pts.push_back(q);
          nextID++;
        }
      }
    }
  }

  const double lx = ax * nx;
  const double ly = by * ny;
  const double lz = cz * nz;
  const double xy = bx * ny;
  const double xz = cx * nz;
  const double yz = cy * nz;

  if (base.box.size() >= 6) {
    const double xMinTilt = std::min({0.0, xy, xz, xy + xz});
    const double xMaxTilt = std::max({0.0, xy, xz, xy + xz});
    const double yMinTilt = std::min(0.0, yz);
    const double yMaxTilt = std::max(0.0, yz);

    const double xloBound = origin[0] + xMinTilt;
    const double yloBound = origin[1] + yMinTilt;
    const double zloBound = origin[2];

    out.box = {lx + xMaxTilt - xMinTilt, ly + yMaxTilt - yMinTilt, lz, xy, xz,
               yz};
    out.boxLow = {xloBound, yloBound, zloBound};
  } else {
    out.box = {lx, ly, lz};
    out.boxLow = {origin[0], origin[1], origin[2]};
  }

  return out;
}

bool dumpDirectedEdges(const std::string &path, const Rows &rows) {
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out.good()) {
    return false;
  }
  for (std::size_t i = 0; i < rows.size(); i++) {
    const auto &row = rows[i];
    if (row.size() <= 1) {
      continue;
    }
    for (std::size_t k = 1; k < row.size(); k++) {
      out << i << " " << row[k] << "\n";
    }
  }
  return out.good();
}

bool parseArgs(int argc, char **argv, Args &args) {
  for (int i = 1; i < argc; i++) {
    const std::string flag = argv[i];
    const auto need = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "missing value for %s\n", name);
        return nullptr;
      }
      i++;
      return argv[i];
    };
    if (flag == "--traj") {
      const char *v = need("--traj");
      if (!v) {
        return false;
      }
      args.traj = v;
    } else if (flag == "--frame") {
      const char *v = need("--frame");
      if (!v) {
        return false;
      }
      args.frame = std::atoi(v);
    } else if (flag == "--type") {
      const char *v = need("--type");
      if (!v) {
        return false;
      }
      args.typeI = std::atoi(v);
    } else if (flag == "--reps") {
      const char *v = need("--reps");
      if (!v) {
        return false;
      }
      args.reps = std::max(1, std::atoi(v));
    } else if (flag == "--warmup") {
      const char *v = need("--warmup");
      if (!v) {
        return false;
      }
      args.warmup = std::max(0, std::atoi(v));
    } else if (flag == "--cutoff") {
      const char *v = need("--cutoff");
      if (!v) {
        return false;
      }
      args.cutoff = std::atof(v);
    } else if (flag == "--k") {
      const char *v = need("--k");
      if (!v) {
        return false;
      }
      args.k = std::max(1, std::atoi(v));
    } else if (flag == "--k-cand") {
      const char *v = need("--k-cand");
      if (!v) {
        return false;
      }
      args.kCand = std::atof(v);
    } else if (flag == "--nx") {
      const char *v = need("--nx");
      if (!v) {
        return false;
      }
      args.nx = std::max(1, std::atoi(v));
    } else if (flag == "--ny") {
      const char *v = need("--ny");
      if (!v) {
        return false;
      }
      args.ny = std::max(1, std::atoi(v));
    } else if (flag == "--nz") {
      const char *v = need("--nz");
      if (!v) {
        return false;
      }
      args.nz = std::max(1, std::atoi(v));
    } else if (flag == "--threads") {
      const char *v = need("--threads");
      if (!v) {
        return false;
      }
      args.threads = std::max(0, std::atoi(v));
    } else if (flag == "--dump-cutoff-edges") {
      const char *v = need("--dump-cutoff-edges");
      if (!v) {
        return false;
      }
      args.dumpCutoffEdges = v;
    } else {
      std::fprintf(stderr, "unknown argument: %s\n", flag.c_str());
      return false;
    }
  }
  if (args.traj.empty()) {
    std::fprintf(stderr, "missing --traj\n");
    return false;
  }
  return true;
}

int effectiveThreads() {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

const char *metaTag(const char *name, const char *fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') {
    return fallback;
  }
  return raw;
}

} // namespace

int main(int argc, char **argv) {
  Args args;
  if (!parseArgs(argc, argv, args)) {
    std::fprintf(stderr,
                 "usage: bench_neigh_only --traj TRAJ [--frame 1] [--type 1] "
                 "[--reps 6] [--warmup 2] [--cutoff 3.5] [--k 4] "
                 "[--k-cand 5.360] [--nx 1] [--ny 1] [--nz 1] [--threads 0] "
                 "[--dump-cutoff-edges PATH]\n");
    return 2;
  }

#ifdef SEAMS_HAS_OPENMP
  if (args.threads > 0) {
    omp_set_num_threads(args.threads);
  }
#endif

  Cloud cloud;
  cloud = sinp::readLammpsTrjO(args.traj, args.frame, cloud, args.typeI);
  if (cloud.nop <= 0) {
    std::fprintf(stderr, "bench_neigh_only: empty frame\n");
    return 1;
  }
  if (args.nx > 1 || args.ny > 1 || args.nz > 1) {
    cloud = replicateCloud(cloud, args.nx, args.ny, args.nz);
  }

  // Cutoff stage: first-call timing and artefacts.
  const auto tCutFirst0 = Clock::now();
  Rows cutoffRows = nneigh::neighListO(args.cutoff, cloud, args.typeI);
  const auto tCutFirst1 = Clock::now();
  const double cutoffFirstMs = ms(tCutFirst0, tCutFirst1);
  const long long cutoffEdges = directedEdges(cutoffRows);
  const std::uint64_t cutoffHash = hashRows(cutoffRows);

  if (!args.dumpCutoffEdges.empty()) {
    const Rows byIndex = nneigh::neighbourListByIndex(cloud, cutoffRows);
    if (!dumpDirectedEdges(args.dumpCutoffEdges, byIndex)) {
      std::fprintf(stderr, "failed to write cutoff edge dump to %s\n",
                   args.dumpCutoffEdges.c_str());
      return 1;
    }
  }

  for (int i = 0; i < args.warmup; i++) {
    (void)nneigh::neighListO(args.cutoff, cloud, args.typeI);
  }

  std::vector<double> cutoffMs;
  cutoffMs.reserve(static_cast<std::size_t>(args.reps));
  for (int i = 0; i < args.reps; i++) {
    const auto t0 = Clock::now();
    (void)nneigh::neighListO(args.cutoff, cloud, args.typeI);
    const auto t1 = Clock::now();
    cutoffMs.push_back(ms(t0, t1));
  }

  // k-nearest stage: first-call timing, backend path and artefacts.
  const auto tKnnFirst0 = Clock::now();
  Rows knnRows =
      nneigh::kNearestNeighbourList(cloud, args.k, args.kCand, args.typeI, true);
  const auto tKnnFirst1 = Clock::now();
  const double knnFirstMs = ms(tKnnFirst0, tKnnFirst1);
  const long long knnEdges = directedEdges(knnRows);
  const std::uint64_t knnHash = hashRows(knnRows);
  const std::string firstPathName = "k-nearest";

  for (int i = 0; i < args.warmup; i++) {
    (void)nneigh::kNearestNeighbourList(cloud, args.k, args.kCand, args.typeI,
                                        true);
  }

  std::vector<double> knnMs;
  knnMs.reserve(static_cast<std::size_t>(args.reps));
  std::vector<std::string> knnPathSamples;
  knnPathSamples.reserve(static_cast<std::size_t>(args.reps));
  std::map<std::string, int> knnPathCounts;
  for (int i = 0; i < args.reps; i++) {
    const auto t0 = Clock::now();
    (void)nneigh::kNearestNeighbourList(cloud, args.k, args.kCand, args.typeI,
                                        true);
    const auto t1 = Clock::now();
    knnMs.push_back(ms(t0, t1));
    const std::string pathName = "k-nearest";
    knnPathSamples.push_back(pathName);
    knnPathCounts[pathName]++;
  }

  const Summary cutoffSummary = summarize(cutoffMs);
  const Summary knnSummary = summarize(knnMs);
  const int threadsUsed = effectiveThreads();
  const char *gitSha = metaTag("SEAMS_GIT_SHA", "unknown");
  const char *dirty = metaTag("SEAMS_GIT_DIRTY", "unknown");

  std::printf(
      "meta traj %s frame %d nop %d reps %d warmup %d cutoff %.6f k %d "
      "kCand %.6f rep %dx%dx%d threads_req %d threads_used %d git_sha %s "
      "dirty %s\n",
      args.traj.c_str(), args.frame, cloud.nop, args.reps, args.warmup,
      args.cutoff, args.k, args.kCand, args.nx, args.ny, args.nz, args.threads,
      threadsUsed, gitSha, dirty);
  std::printf("cutoff_first_ms %.6f\n", cutoffFirstMs);
  std::printf("cutoff_samples_ms %s\n", joinSamples(cutoffMs).c_str());
  std::printf("cutoff_ms median %.6f q1 %.6f q3 %.6f iqr %.6f\n",
              cutoffSummary.median, cutoffSummary.q1, cutoffSummary.q3,
              cutoffSummary.q3 - cutoffSummary.q1);
  std::printf("cutoff_directed_edges %lld\n", cutoffEdges);
  std::printf("cutoff_hash %llu\n",
              static_cast<unsigned long long>(cutoffHash));

  std::printf("knn_first_ms %.6f\n", knnFirstMs);
  std::printf("knn_samples_ms %s\n", joinSamples(knnMs).c_str());
  std::printf("knn_ms median %.6f q1 %.6f q3 %.6f iqr %.6f\n",
              knnSummary.median, knnSummary.q1, knnSummary.q3,
              knnSummary.q3 - knnSummary.q1);
  std::printf("knn_directed_edges %lld\n", knnEdges);
  std::printf("knn_hash %llu\n", static_cast<unsigned long long>(knnHash));
  std::printf("knn_path_first %s\n", firstPathName.c_str());
  std::printf("knn_path_samples %s\n", joinNames(knnPathSamples).c_str());
  std::printf("knn_path_counts");
  for (const auto &[name, count] : knnPathCounts) {
    std::printf(" %s=%d", name.c_str(), count);
  }
  std::printf("\n");

  return 0;
}
