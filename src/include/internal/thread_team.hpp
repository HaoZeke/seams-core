#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <vector>

// One width for OpenMP and Rayon. OpenMP already follows
// OMP_NUM_THREADS. Linkcell's pool does not, so a process whose
// hardware reports one cpu keeps the neighbour walk on one thread
// even when OpenMP was given several. lc_set_threads is read by the
// next linkcell search.
//
// OMP_PROC_BIND makes the OpenMP runtime pin the calling thread to its
// first place when the library loads, and threads created afterwards, a
// Rayon pool among them, inherit that one cpu. widenAffinity gives the
// caller back the cpus of every place, which the runtime took from the
// launch mask (taskset, numactl, a cgroup cpuset) and OMP_PLACES. It adds
// no cpu outside those places. A caller the runtime did not bind, and a
// one-thread width, keep their mask. The pool is built once and its
// workers sleep between searches. The OpenMP team is not paused: a pause
// makes the next region pay a wake-up that moves by tens of milliseconds.

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

extern "C" void lc_set_threads(int n);
extern "C" int lc_num_threads(void);

namespace seams {

inline void widenAffinity() {
#if defined(__linux__) && defined(SEAMS_HAS_OPENMP)
  static std::mutex lock;
  static bool done = false;
  if (omp_in_parallel() || omp_get_max_threads() < 2 || omp_get_place_num() < 0) {
    return;
  }
  const std::lock_guard<std::mutex> hold(lock);
  if (done) {
    return;
  }
  done = true;
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) != 0) {
    return;
  }
  const int places = omp_get_num_places();
  for (int place = 0; place < places; ++place) {
    const int count = omp_get_place_num_procs(place);
    if (count <= 0) {
      continue;
    }
    std::vector<int> ids(static_cast<std::size_t>(count));
    omp_get_place_proc_ids(place, ids.data());
    for (const int cpu : ids) {
      if (cpu >= 0 && cpu < CPU_SETSIZE) {
        CPU_SET(cpu, &mask);
      }
    }
  }
  sched_setaffinity(0, sizeof(mask), &mask);
#endif
}

inline void releaseTeam() {
#ifdef SEAMS_HAS_OPENMP
  omp_pause_resource_all(omp_pause_soft);
#endif
}

// OpenMP work that does not hand the cores to Rayon. Pausing the team
// here parks every worker and the next parallel region pays a wake-up
// that moves by several milliseconds from call to call. Widen the mask
// and freeze the team size, and leave the workers asleep only under
// the passive wait policy.
inline void prepareCompute() {
  // Linkcell's pool sleeps on a condvar between searches. Forcing OpenMP
  // passive as well made the next region wake on a spread of tens of
  // milliseconds, which is the slow 8-thread sample. Leave OpenMP on the
  // active wait the benchmark already uses.
  widenAffinity();
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
#endif
}

inline void adoptThreads() {
  // Record the Rayon width and leave the OpenMP team running. Pausing
  // it, or forcing a passive wait, made the next region late by tens of
  // milliseconds. The pool sleeps on its own condvar between searches.
  prepareCompute();
#ifdef SEAMS_HAS_OPENMP
  if (omp_in_parallel()) {
#ifdef SEAMS_HAS_LINKCELL
    lc_set_threads(1);
#endif
    return;
  }
  int n = omp_get_max_threads();
#else
  int n = 1;
#endif
  if (n < 1) {
    n = 1;
  }
#ifdef SEAMS_HAS_LINKCELL
  lc_set_threads(n);
#else
  (void)n;
#endif
}

// Thread count for a call that is too small for the fork. Below 2048
// atoms an eight-thread region is slower than the serial walk. Below
// 8192 atoms four threads beat eight: the extra workers do not have a
// full share of the walk.
inline int teamFor(int n) {
#ifdef SEAMS_HAS_OPENMP
  const int team = std::max(1, omp_get_max_threads());
  if (n < 2048) {
    return 1;
  }
  if (n < 8192 && team > 4) {
    return 4;
  }
  return team;
#else
  (void)n;
  return 1;
#endif
}

// Sets the OpenMP width from teamFor for the rest of the scope, then
// puts the caller's width back. An OpenMP-only kernel uses this instead
// of adoptThreads, which would also start the Rayon pool.
struct TeamWidth {
  int saved = 1;
  bool on = false;
#ifdef __linux__
  cpu_set_t oldMask;
#endif
  explicit TeamWidth(int n) {
#ifdef SEAMS_HAS_OPENMP
    if (!omp_in_parallel()) {
      saved = std::max(1, omp_get_max_threads());
      const int want = teamFor(n);
      if (want != saved) {
        omp_set_num_threads(want);
        on = true;
      }
#ifdef __linux__
      // A four-thread team on an eight-core mask still migrates across
      // all eight. Bind it to the first cores the process is allowed,
      // before the workers are created, so they inherit that mask.
      if (want > 0 && want < saved &&
          sched_getaffinity(0, sizeof(oldMask), &oldMask) == 0) {
        int cpus[256];
        int ncpu = 0;
        for (int cpu = 0; cpu < CPU_SETSIZE && ncpu < want && ncpu < 256; ++cpu) {
          if (CPU_ISSET(cpu, &oldMask)) {
            cpus[ncpu++] = cpu;
          }
        }
        if (ncpu > 0) {
          // Once per process. Repeating the bind on a call of a few
          // hundred microseconds was the whole gap between four threads
          // requested directly and eight narrowed to four. OpenMP's own
          // bind would put the workers back on the whole allowed set,
          // so it is turned off before the team is created.
          // Thread 0 is the caller and stays on its own mask: a pinned
          // caller hands its one cpu to every OpenMP worker and Rayon
          // thread it creates later, and the next eight-thread region
          // then runs five threads on that cpu.
          static int bound = 0;
          if (bound != want) {
            setenv("OMP_PROC_BIND", "false", 1);
#pragma omp parallel num_threads(want)
            {
              const int tid = omp_get_thread_num();
              if (tid != 0) {
                cpu_set_t mask;
                CPU_ZERO(&mask);
                CPU_SET(cpus[tid < ncpu ? tid : 0], &mask);
                sched_setaffinity(0, sizeof(mask), &mask);
              }
            }
            bound = want;
          }
        }
      }
#endif
    }
#else
    (void)n;
#endif
  }
  ~TeamWidth() {
#ifdef SEAMS_HAS_OPENMP
    if (on) {
      omp_set_num_threads(saved);
    }
#endif
  }
  TeamWidth(const TeamWidth &) = delete;
  TeamWidth &operator=(const TeamWidth &) = delete;
};

// Team width for one library call on a machine other processes share.
// When their threads hold the cores, a full-width team waits at each fork
// that follows serial work for workers the scheduler has parked, a slice
// at a time, and a call of a few milliseconds grows tenfold. A narrower
// team leaves cores free and does not wait. One small region at entry,
// after the caller's own gap, takes at most a few hundred microseconds on
// an idle machine, where sleeping workers wake, and milliseconds on a
// contended one. Two slow entries in a row halve the team, down to a
// quarter of the requested width; at most once a second a narrowed team
// tries twice its width again. Results do not depend on the width.
// SEAMS_FAIR_WIDTH=0 keeps the requested width.
struct FairState {
  std::mutex lock;
  int base = 0;
  int granted = 0;
  int slow = 0;
  double tried = 0.0;
};

inline FairState &fairState() {
  static FairState state;
  return state;
}

inline bool fairEnabled() {
  static const bool on = [] {
    const char *flag = std::getenv("SEAMS_FAIR_WIDTH");
    return flag == nullptr || std::strcmp(flag, "0") != 0;
  }();
  return on;
}

class FairWidth {
public:
  explicit FairWidth(bool adopt = false) {
#ifdef SEAMS_HAS_OPENMP
    if (depth()++ == 0 && !omp_in_parallel() && fairEnabled()) {
      choose();
    }
#endif
    if (adopt) {
      adoptThreads();
    }
  }
  ~FairWidth() {
#ifdef SEAMS_HAS_OPENMP
    if (set_) {
      omp_set_num_threads(asked_);
    }
    --depth();
#endif
  }
  FairWidth(const FairWidth &) = delete;
  FairWidth &operator=(const FairWidth &) = delete;

private:
  static constexpr double kSlowFork = 1.0e-3;
  static constexpr double kRetry = 1.0;
  int asked_ = 1;
  bool set_ = false;

  static int &depth() {
    thread_local int level = 0;
    return level;
  }

#ifdef SEAMS_HAS_OPENMP
  static double forkSeconds(int width) {
    const double start = omp_get_wtime();
#pragma omp parallel num_threads(width)
    {
      [[maybe_unused]] volatile int touch = omp_get_thread_num();
    }
    return omp_get_wtime() - start;
  }

  void choose() {
    asked_ = std::max(1, omp_get_max_threads());
    if (asked_ < 2) {
      return;
    }
    FairState &state = fairState();
    const std::lock_guard<std::mutex> hold(state.lock);
    if (state.base != asked_) {
      state.base = asked_;
      state.granted = asked_;
      state.slow = 0;
      state.tried = 0.0;
    }
    const int floor = std::max(1, state.base / 4);
    const double now = omp_get_wtime();
    int width = state.granted;
    if (width < state.base && now - state.tried >= kRetry) {
      width = std::min(state.base, 2 * width);
      state.tried = now;
    }
    if (forkSeconds(width) <= kSlowFork) {
      state.slow = 0;
      state.granted = width;
    } else if (width == state.granted && ++state.slow >= 2 && width > floor) {
      state.granted = std::max(floor, width / 2);
      state.slow = 0;
      state.tried = now;
    }
    if (state.granted != asked_) {
      omp_set_num_threads(state.granted);
      set_ = true;
    }
  }
#endif
};

} // namespace seams
