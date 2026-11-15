/* GCC's InitTlsSize calls dlsym(RTLD_NEXT, "_dl_get_tls_static_info")
   and jumps to the result. RTLD_NEXT does not see ld-linux once libasan
   is the last object, which is what happens when a sanitized library is
   dlopened into R. Preload this file ahead of that dlopen. It is not
   itself sanitized. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <string.h>

void *dlsym(void *handle, const char *symbol) {
  static void *(*real_dlsym)(void *, const char *) = 0;
  if (real_dlsym == 0) {
    real_dlsym = (void *(*)(void *, const char *))dlvsym(
        RTLD_NEXT, "dlsym", "GLIBC_2.2.5");
  }
  if (strcmp(symbol, "_dl_get_tls_static_info") == 0) {
    return dlvsym(RTLD_DEFAULT, symbol, "GLIBC_PRIVATE");
  }
  return real_dlsym(handle, symbol);
}
