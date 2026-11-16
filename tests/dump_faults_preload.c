/* Copyright (c) d-SEAMS developers */
/* Test-side faults for pread and fstat. The library keeps calling libc.
 * ELF loaders take the overrides from LD_PRELOAD by symbol name. dyld
 * resolves two-level names, so on macOS the same functions are registered
 * in __DATA,__interpose and dump_fault_wrap.sh inserts the library with
 * DYLD_INSERT_LIBRARIES. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __APPLE__
#include <stdint.h>
#define SEAMS_INTERPOSE(replacement, replacee)                                 \
  __attribute__((used)) static const struct {                                 \
    const void *with;                                                         \
    const void *what;                                                         \
  } interpose_##replacee __attribute__((section("__DATA,__interpose"))) = {   \
      (const void *)(uintptr_t)&replacement, (const void *)(uintptr_t)&replacee}
#endif

static int pread_faults_left(void) {
  const char *n = getenv("DUMP_PREAD_EINTR");
  if (n == NULL || n[0] == '\0') {
    return 0;
  }
  return atoi(n);
}

#ifdef __APPLE__
/* Interposed calls from this image reach libc itself. */
static ssize_t fault_pread(int fd, void *buf, size_t count, off_t offset) {
  static int left = -1;
  if (left < 0) {
    left = pread_faults_left();
  }
  if (getenv("DUMP_PREAD_EIO") != NULL) {
    errno = EIO;
    return -1;
  }
  if (left > 0) {
    left -= 1;
    errno = EINTR;
    return -1;
  }
  return pread(fd, buf, count, offset);
}
SEAMS_INTERPOSE(fault_pread, pread);

static int fstat_seen = 0;

/* dyld interposes libSystem's internal calls too: stdio stats a stream
 * before it buffers it. glibc's own calls bypass LD_PRELOAD, so only
 * calls from outside the system libraries count here as well. */
static int fault_fstat(int fd, struct stat *buf) {
  Dl_info caller;
  const int system_call =
      dladdr(__builtin_return_address(0), &caller) != 0 && caller.dli_fname != NULL &&
      strncmp(caller.dli_fname, "/usr/lib/", 9) == 0;
  if (!system_call) {
    fstat_seen += 1;
    /* The parse records one successful stat, then the change check fails. */
    if (getenv("DUMP_FSTAT_EIO") != NULL && fstat_seen > 1) {
      errno = EIO;
      return -1;
    }
  }
  return fstat(fd, buf);
}
SEAMS_INTERPOSE(fault_fstat, fstat);
#else
ssize_t pread(int fd, void *buf, size_t count, off_t offset) {
  static ssize_t (*real_pread)(int, void *, size_t, off_t) = NULL;
  static int left = -1;
  if (real_pread == NULL) {
    real_pread = dlsym(RTLD_NEXT, "pread");
  }
  if (left < 0) {
    left = pread_faults_left();
  }
  if (getenv("DUMP_PREAD_EIO") != NULL) {
    errno = EIO;
    return -1;
  }
  if (left > 0) {
    left -= 1;
    errno = EINTR;
    return -1;
  }
  return real_pread(fd, buf, count, offset);
}

static int fstat_seen = 0;

static int fstat_fault(void) {
  fstat_seen += 1;
  /* The parse records one successful stat, then the change check fails. */
  if (getenv("DUMP_FSTAT_EIO") != NULL && fstat_seen > 1) {
    errno = EIO;
    return -1;
  }
  return 0;
}

int fstat(int fd, struct stat *buf) {
  static int (*real_fstat)(int, struct stat *) = NULL;
  if (real_fstat == NULL) {
    real_fstat = dlsym(RTLD_NEXT, "fstat");
  }
  if (fstat_fault() != 0) {
    return -1;
  }
  return real_fstat(fd, buf);
}

/* glibc with _FILE_OFFSET_BITS=64 calls __fxstat64, not fstat. */
int __fxstat64(int ver, int fd, struct stat64 *buf) {
  static int (*real_fxstat)(int, int, struct stat64 *) = NULL;
  if (real_fxstat == NULL) {
    real_fxstat = dlsym(RTLD_NEXT, "__fxstat64");
  }
  if (fstat_fault() != 0) {
    return -1;
  }
  return real_fxstat(ver, fd, buf);
}
#endif
