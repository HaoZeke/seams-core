// Test-only allocator. The library does not know it exists.
// Arm it, call one library function, then read the counters.

#include <dlfcn.h>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <new>

namespace {

using NewFn = void *(*)(std::size_t);
using DeleteFn = void (*)(void *);
using SizedDeleteFn = void (*)(void *, std::size_t);

NewFn real_new() {
  static NewFn fn = reinterpret_cast<NewFn>(dlsym(RTLD_NEXT, "_Znwm"));
  return fn;
}
NewFn real_new_array() {
  static NewFn fn = reinterpret_cast<NewFn>(dlsym(RTLD_NEXT, "_Znam"));
  return fn;
}
DeleteFn real_delete() {
  static DeleteFn fn = reinterpret_cast<DeleteFn>(dlsym(RTLD_NEXT, "_ZdlPv"));
  return fn;
}
DeleteFn real_delete_array() {
  static DeleteFn fn = reinterpret_cast<DeleteFn>(dlsym(RTLD_NEXT, "_ZdaPv"));
  return fn;
}
SizedDeleteFn real_sized_delete() {
  static SizedDeleteFn fn = reinterpret_cast<SizedDeleteFn>(dlsym(RTLD_NEXT, "_ZdlPvm"));
  return fn;
}
SizedDeleteFn real_sized_delete_array() {
  static SizedDeleteFn fn = reinterpret_cast<SizedDeleteFn>(dlsym(RTLD_NEXT, "_ZdaPvm"));
  return fn;
}

std::atomic<long> g_live{0};
std::atomic<long> g_peak{0};
std::atomic<long> g_serial{0};
std::atomic<int> g_skip{0};
std::atomic<int> g_kind{0};
std::atomic<long> g_throw_at{-1};
std::atomic<int> g_env_ready{0};

void read_env() {
  if (g_env_ready.exchange(1) != 0) {
    return;
  }
  const char *kind = getenv("SEAMS_TEST_KIND");
  const char *at = getenv("SEAMS_TEST_AT");
  if (kind != nullptr && std::strcmp(kind, "bad_alloc") == 0) {
    g_kind.store(1);
  } else if (kind != nullptr && std::strcmp(kind, "int") == 0) {
    g_kind.store(2);
  }
  if (at != nullptr) {
    g_throw_at.store(std::strtol(at, nullptr, 10));
  }
}

void note_alloc() {
  const long now = g_live.fetch_add(1) + 1;
  long seen = g_peak.load();
  while (now > seen && !g_peak.compare_exchange_weak(seen, now)) {
  }
}

void note_free() { g_live.fetch_sub(1); }

void *allocate(std::size_t n, NewFn real) {
  read_env();
  const long serial = g_serial.fetch_add(1) + 1;
  const int kind = g_kind.load();
  const long at = g_throw_at.load();
  if (kind != 0 && at >= 0 && serial == at) {
    g_kind.store(0);
    if (kind == 2) {
      throw 1;
    }
    throw std::bad_alloc();
  }
  if (kind != 0 && at < 0) {
    const int left = g_skip.load();
    if (left <= 0) {
      g_kind.store(0);
      if (kind == 2) {
        throw 1;
      }
      throw std::bad_alloc();
    }
    g_skip.store(left - 1);
  }
  if (std::getenv("SEAMS_TRACE_SIZE") != nullptr && n >= 64 && n <= 128) {
    std::fprintf(stderr, "size-serial %ld %zu\n", serial, n);
  }
  void *pointer = real(n);
  note_alloc();
  return pointer;
}

} // namespace

extern "C" {

void seams_test_arm(int skip, int kind) {
  g_skip.store(skip);
  g_kind.store(kind);
}

long seams_test_live() { return g_live.load(); }

long seams_test_serial() { return g_serial.load(); }

long seams_test_peak() { return g_peak.load(); }

void seams_test_reset_peak() { g_peak.store(g_live.load()); }

}

void *operator new(std::size_t n) { return allocate(n, real_new()); }

void *operator new[](std::size_t n) { return allocate(n, real_new_array()); }

void operator delete(void *pointer) noexcept {
  if (pointer != nullptr) {
    note_free();
  }
  real_delete()(pointer);
}

void operator delete[](void *pointer) noexcept {
  if (pointer != nullptr) {
    note_free();
  }
  real_delete_array()(pointer);
}

void operator delete(void *pointer, std::size_t n) noexcept {
  if (pointer != nullptr) {
    note_free();
  }
  real_sized_delete()(pointer, n);
}

void operator delete[](void *pointer, std::size_t n) noexcept {
  if (pointer != nullptr) {
    note_free();
  }
  real_sized_delete_array()(pointer, n);
}

namespace {
struct FaultReporter {
  ~FaultReporter() {
    const char *kind = std::getenv("SEAMS_TEST_KIND");
    if (kind != nullptr && std::strcmp(kind, "count") == 0) {
      std::fprintf(stderr, "fault-serial %ld\n", static_cast<long>(g_serial.load()));
    }
  }
};
FaultReporter fault_reporter;
}
