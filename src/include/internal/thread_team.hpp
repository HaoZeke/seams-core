#pragma once

#include <cstdlib>

// One width for OpenMP and Rayon. OpenMP already follows
// OMP_NUM_THREADS. Linkcell's pool does not, so a process whose
// hardware reports one cpu keeps the neighbour walk on one thread
// even when OpenMP was given several. lc_set_threads is read by the
// next linkcell search.
//
// The process mask is whatever the user or the scheduler inherited
// to us. A taskset or a cpuset stays in force: the pool is counted
// from that mask and is never placed on a cpu outside it. A live
// OpenMP team is asked to stop spinning so the search can use the
// cpus the mask already names.

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

extern "C" void lc_set_threads(int n);
extern "C" int lc_num_threads(void);

namespace seams {

// Cpus in the mask inherited from the parent. Zero when the kernel
// has no affinity call. The mask is only read.
inline int affinityCpus() {
#ifdef __linux__
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) != 0) {
    return 0;
  }
  int n = 0;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, &set)) {
      ++n;
    }
  }
  return n;
#else
  return 0;
#endif
}

// OpenMP's requested width, cut down to the inherited mask so a
// restricted taskset is not oversubscribed. A proc-bind policy pins
// this thread to one place, so this thread's mask is a single cpu
// while the place count is still the mask the process started with.
inline int poolWidth() {
#ifdef SEAMS_HAS_OPENMP
  int n = omp_get_max_threads();
#else
  int n = 1;
#endif
  if (n < 1) {
    n = 1;
  }
  int limit = affinityCpus();
#ifdef SEAMS_HAS_OPENMP
  const int places = omp_get_num_procs();
  const char *bind = std::getenv("OMP_PROC_BIND");
  const bool procBind = bind != nullptr && bind[0] != '\0' && bind[0] != 'f' &&
                        bind[0] != 'F' && bind[0] != '0';
  if (places > 0 && places < limit) {
    limit = places;
  } else if (procBind && places > limit) {
    limit = places;
  }
#endif
  if (limit > 0 && n > limit) {
    n = limit;
  }
  return n;
}

inline void releaseTeam() {
#ifdef SEAMS_HAS_OPENMP
  omp_pause_resource_all(omp_pause_soft);
#endif
}

// OpenMP work that does not hand the cores to Rayon. Freeze the team
// size inside the inherited mask. Do not force a passive wait: a team
// that sleeps between calls wakes up on a spread of several
// milliseconds.
inline void prepareCompute() {
  // Do not force a passive wait. A team that sleeps between calls wakes
  // up on a spread of several milliseconds; the benchmark keeps the
  // team with OMP_WAIT_POLICY=active.
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  if (!omp_in_parallel()) {
    omp_set_num_threads(poolWidth());
  }
#endif
}

inline void adoptThreads() {
  // Idle OpenMP threads spin by default and take cores from Rayon.
  setenv("OMP_WAIT_POLICY", "PASSIVE", 0);
  releaseTeam();
#ifdef SEAMS_HAS_OPENMP
  // A domain is searched by the thread that owns it. Rayon stays at one
  // so the OpenMP team is the only pool.
  if (omp_in_parallel()) {
#ifdef SEAMS_HAS_LINKCELL
    lc_set_threads(1);
#endif
    return;
  }
#endif
  const int n = poolWidth();
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(n);
#endif
#ifdef SEAMS_HAS_LINKCELL
  lc_set_threads(n);
#else
  (void)n;
#endif
}

// Thread count for a call that is too small for the fork. Measured on
// this host: below 2048 atoms an eight-thread region is slower than the
// serial walk, for both the cutoff list and the CHILL+ cell walk.
inline int teamFor(int n) {
  const int team = poolWidth();
  if (n < 2048) {
    return 1;
  }
  return team;
}

} // namespace seams
