/* Map CLOCK_REALTIME onto an anchor plus elapsed CLOCK_MONOTONIC.
   Meson reads CLOCK_REALTIME when it records a test. A step of that
   clock makes the record disagree with a monotonic interval taken
   around the same process. CLOCK_MONOTONIC keeps advancing while the
   process waits or is not scheduled. The mapping hides those realtime
   steps. It is not a reading of civil time.
   The timing script sets LD_PRELOAD for the pixi process. Meson and
   every child it starts inherit that setting. A result that does not
   read civil time is unchanged.
   The anchor is taken once. pthread_once covers two threads that
   enter the first call together. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

static int (*real_clock_gettime)(clockid_t, struct timespec *);
static struct timespec rt0, mono0;
static pthread_once_t once = PTHREAD_ONCE_INIT;

static void init_clock_once(void) {
  real_clock_gettime = dlsym(RTLD_NEXT, "clock_gettime");
  real_clock_gettime(CLOCK_REALTIME, &rt0);
  real_clock_gettime(CLOCK_MONOTONIC, &mono0);
}

static void init_clock(void) { pthread_once(&once, init_clock_once); }

static void fold_nsec(long *sec, long *nsec) {
  if (*nsec >= 1000000000L) {
    *sec += 1;
    *nsec -= 1000000000L;
  } else if (*nsec < 0) {
    *sec -= 1;
    *nsec += 1000000000L;
  }
}

#ifdef SEAMS_CLOCK_TEST
void seams_clock_fold_test(long *sec, long *nsec) { fold_nsec(sec, nsec); }
#endif

int clock_gettime(clockid_t id, struct timespec *tp) {
  init_clock();
#if defined(CLOCK_REALTIME_COARSE)
  if (id == CLOCK_REALTIME || id == CLOCK_REALTIME_COARSE) {
#else
  /* Apple's clockid_t has CLOCK_REALTIME and no coarse realtime clock. */
  if (id == CLOCK_REALTIME) {
#endif
    struct timespec mono;
    long sec;
    long nsec;
    real_clock_gettime(CLOCK_MONOTONIC, &mono);
    sec = rt0.tv_sec + (mono.tv_sec - mono0.tv_sec);
    nsec = rt0.tv_nsec + (mono.tv_nsec - mono0.tv_nsec);
    fold_nsec(&sec, &nsec);
    tp->tv_sec = sec;
    tp->tv_nsec = nsec;
    return 0;
  }
  return real_clock_gettime(id, tp);
}

int gettimeofday(struct timeval *tv, struct timezone *tz) {
  struct timespec ts;
  (void)tz;
  clock_gettime(CLOCK_REALTIME, &ts);
  tv->tv_sec = ts.tv_sec;
  tv->tv_usec = ts.tv_nsec / 1000;
  return 0;
}
